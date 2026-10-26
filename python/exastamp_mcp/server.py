"""
exaStamp MCP server (stdio): tools that let an LLM agent write exaStamp input
files that run the first time.

    exastamp-mcp [--workspace DIR] [--doc-dir doc_exaStamp/docs] [--max-run-steps N]
    exastamp-mcp --index-examples      # build every example once, cache the status

The workspace (default: the current directory) is where inputs are written
and run; paths outside it are refused. --max-run-steps caps every run_input
to N time steps (evaluations: an agent's run only checks that the input starts).
"""

import argparse
import functools
import os
import sys
import threading
from typing import Any, Literal, Optional

import yaml as _yaml
from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

from onika_mcp.catalog import Catalog
from onika_mcp.jobs import JobManager
from onika_mcp.smoke import check_build

from . import PACKAGE
from .corpus import Corpus
from .data import POTENTIAL_FORMATS, potential_files
from .knowledge import SCHEMAS, doc_url
from .render import RenderError, render_input
from .spec import SimulationSpec
from .validate import validate_input as _validate

__all__ = ["build_server", "main", "WORKFLOW"]

WORKFLOW = """\
exaStamp input files (.msp) are YAML. Follow these steps to get an input that runs the first time:

1. list_capabilities: what this build supports (potentials, integrators, lattices, outputs).
2. If the request fits the structured description (lattice or extended XYZ system, EAM alloy or
   Lennard-Jones potential, NVE / Langevin / Berendsen), call build_input with a SimulationSpec:
   the server derives the domain, the lattice types, the grid flavor and the block order, and
   validates the result. Prefer this path.
3. Otherwise start from a working example: search_examples (prefer results with builds: true),
   get_example, and describe_block for every block you change. For an operator's meaning,
   parameters and units, or a feature no example shows: search_docs, then get_doc on the page.
4. validate_input before handing over an input; fix every error it reports (each issue has a
   path, a message and often a suggestion) and validate again. Warnings deserve a look.
5. write_input saves the file in the workspace (refused while validation fails);
   run_input / job_status run it locally.

Rules:
- Quantities carry units, written without spaces inside the unit: 3.3 ang, 300 K, 1.0e-3 ps,
  0.1 ps^-1, 50. m/s, 180.95 Da, 0 e-, 0.0104 eV. Valid symbols: m mm um nm ang kg g Da h s us ns
  ps fs C e- K mol particle cd rad degree J kJ eV cal kcal. 'A', 'Angstrom', 'amu' are not units.
- Never guess a potential file name: list_potential_files gives the files of this install, and
  species must be named like the elements in an EAM file.
- Parameter names are checked: an unknown one is an error, not ignored. Use describe_block.
- setup_system for a lattice is: domain, init_rcb_grid, lattice, init_temperature_new.
  domain bounds must equal grid_dims x cell_size exactly along each axis.
- A lattice needs one type per atom of the unit cell: SC 1, BCC 2, FCC/HCP/2BCT 4, diamond 8.
"""


def _resolve(workspace, path):
    full = os.path.realpath(os.path.join(workspace, path))
    if not (full == workspace or full.startswith(workspace + os.sep)):
        raise ValueError(f"'{path}' is outside the workspace {workspace}")
    return full


class _State:
    """Lazy catalog and corpus (the first catalog load starts exaStamp once)."""

    def __init__(self, workspace, doc_dir, max_run_steps=None):
        self.workspace = os.path.realpath(workspace)
        self.doc_dir = doc_dir
        self.max_run_steps = max_run_steps
        self.catalog = Catalog(PACKAGE)
        self.jobs = JobManager(PACKAGE)
        self._corpus = None
        self._lock = threading.Lock()

    @property
    def corpus(self):
        with self._lock:
            if self._corpus is None:
                c = Corpus(self.catalog.paths["XSP_DATA_DIR"], self.doc_dir)
                c.load_build_status(os.path.join(self.catalog.cache_dir, f"{self.catalog.cache_key}-examples.json"))
                self._corpus = c
            return self._corpus


def build_server(workspace=".", doc_dir=None, max_run_steps=None):
    state = _State(workspace, doc_dir or os.environ.get("EXASTAMP_DOC_DIR"), max_run_steps)
    server = MCPServer("exastamp", title="exaStamp input builder", instructions=WORKFLOW,
                       description="Write, check and run exaStamp molecular dynamics input files")
    ro = ToolAnnotations(read_only_hint=True, open_world_hint=False)

    def tool(annotations):
        """server.tool with structured results, and expected failures (bad name, missing
        file...) reported to the model with their message instead of a generic error."""
        def register(fn):
            @functools.wraps(fn)
            def wrapper(*args, **kwargs):
                try:
                    return fn(*args, **kwargs)
                except (ValueError, FileNotFoundError, KeyError, RuntimeError) as e:
                    message = e.args[0] if isinstance(e, KeyError) and e.args else str(e)
                    raise ToolError(message) from e
            return server.tool(annotations=annotations, structured_output=True)(wrapper)
        return register

    @tool(ro)
    def get_workflow() -> str:
        """The recommended sequence of tool calls and the input-file rules. Read it first."""
        return WORKFLOW

    @tool(ro)
    def list_capabilities() -> dict[str, Any]:
        """What this exaStamp build supports: potentials (with whether build_input handles them and
        how many examples use them), numerical schemes, lattice structures, snapshot writers, grid flavors."""
        cat, corp = state.catalog, state.corpus
        used = corp.tag_values("potential")
        potentials = {}
        for name in sorted(set(SCHEMAS["compute_force"]["operators"]) | set(POTENTIAL_FORMATS.values()) | set(used)):
            if cat.get(name) is None and name not in cat.defaults:
                continue
            potentials[name] = {"examples": used.get(name, 0), "build_input": name in ("eam_alloy_force", "lj_compute_force"),
                                "summary": SCHEMAS.get(name, {}).get("summary", "")}
        return {
            "potentials": dict(sorted(potentials.items(), key=lambda kv: -kv[1]["examples"])),
            "numerical_schemes": SCHEMAS["numerical_scheme"]["values"],
            "build_input_ensembles": ["nve", "langevin", "berendsen"],
            "lattice_structures": SCHEMAS["lattice"]["structure_types_count"],
            "snapshot_writers": SCHEMAS["write_snapshot"]["values"],
            "grid_flavors": SCHEMAS["grid_flavor"]["values"],
            "documentation": bool(corp.doc_dir),
        }

    @tool(ro)
    def describe_block(name: str) -> dict[str, Any]:
        """Everything known about an input block or operator (e.g. 'domain', 'lattice', 'global',
        'langevin_thermostat', 'eam_alloy_force'): parameters with type, unit dimension, default and
        allowed values; required parameters; a verified snippet; the documentation page."""
        schema = SCHEMAS.get(name)
        op = state.catalog.get(name)
        if schema is None and op is None:
            close = state.catalog.closest(name) + [k for k in SCHEMAS if name in k][:3]
            raise ValueError(f"unknown block or operator '{name}'" + (f"; did you mean {', '.join(close)}?" if close else ""))
        out = {"name": name}
        if schema:
            out.update({k: v for k, v in schema.items() if k != "doc"})
            if schema.get("doc"):
                out["doc"] = doc_url(schema["doc"])
        if op:
            if op.get("documentation"):
                out["operator_documentation"] = op["documentation"][:3000]
            slots = [s for s in op.get("slots", []) if s["direction"] != "out"]
            out["operator_parameters"] = [
                {k: v for k, v in s.items() if k in ("name", "type", "required", "documentation") and v not in ("", None)}
                for s in slots if not s["type"].startswith(("Grid", "onika::", "exanb::Grid", "ompi_"))]
        if name in state.catalog.defaults:
            out["default_definition"] = state.catalog.defaults[name]
        return out

    @tool(ro)
    def search_operators(query: str, limit: int = 20) -> list[dict[str, Any]]:
        """Registered operators whose name or documentation contains every word of query."""
        return state.catalog.search(query, limit)

    @tool(ro)
    def search_examples(query: str = "", tags: Optional[list[str]] = None, limit: int = 8) -> list[dict[str, Any]]:
        """Shipped example inputs matching the query words (e.g. 'eam Ta langevin', 'water rigid').
        tags filters on exact tag values (operator, scheme or structure names, e.g. ['verlet_nhnpt']).
        Results with builds: false do not work with this build; prefer builds: true."""
        return state.corpus.search_examples(query, tags, limit)

    @tool(ro)
    def get_example(path: str) -> dict[str, Any]:
        """Full text of an example (path as returned by search_examples) and the files next to it."""
        return state.corpus.get_example(path)

    @tool(ro)
    def search_docs(query: str, limit: int = 5) -> list[dict[str, Any]]:
        """Sections of the exaStamp documentation matching the query words, each with the list of
        sections of its page. Use get_doc for the whole page when a hit is truncated or when the
        syntax, parameters or examples are in another section."""
        if not state.corpus.doc_dir:
            raise ValueError("documentation not installed: start the server with --doc-dir or EXASTAMP_DOC_DIR")
        return state.corpus.search_docs(query, limit)

    @tool(ro)
    def get_doc(page: str) -> dict[str, Any]:
        """Full markdown of a documentation page (page as returned by search_docs): operator
        syntax, parameter tables with defaults and units, usage examples."""
        return state.corpus.get_doc(page)

    @tool(ro)
    def list_potential_files(operator: Optional[str] = None, element: Optional[str] = None) -> list[dict[str, Any]]:
        """Potential files of this install, with elements and cutoff for EAM alloy files.
        operator filters by force operator (e.g. 'eam_alloy_force'), element by chemical symbol.
        Files in the workspace can be used too."""
        found = potential_files(state.catalog.data_dirs, operator, element)
        return [{k: v for k, v in f.items() if k != "directory"} for f in found]

    @tool(ro)
    def build_input(spec: SimulationSpec, validate: Literal["none", "static", "build", "smoke"] = "smoke") -> dict[str, Any]:
        """Turn a structured simulation description into a complete, commented input file and
        validate it (default: static checks and a 2-step smoke run, a few seconds). Returns the
        YAML, notes on derived choices (cell size, masses, grid flavor...) and the validation report."""
        try:
            text, notes = render_input(spec, state.catalog.data_dirs, state.catalog.find_data_file)
        except RenderError as e:
            return {"ok": False, "error": str(e)}
        out = {"ok": True, "yaml": text, "notes": notes}
        if validate != "none":
            report = _validate(text, state.catalog, state.workspace, level=validate)
            out["ok"] = report["ok"]
            out["validation"] = report
        return out

    @tool(ro)
    def validate_input(yaml: str, level: Literal["static", "build", "smoke"] = "smoke", steps: int = 2) -> dict[str, Any]:
        """Check an input file. static: parameters, units, files, consistency rules (instant).
        build: exaStamp builds the operator graph in strict mode. smoke: runs `steps` time steps with
        outputs disabled and checks the result (finite energies, no temperature explosion, atom count).
        Each layer runs only if the previous one passed. Relative file names resolve in the workspace."""
        return _validate(yaml, state.catalog, state.workspace, level=level, steps=steps)

    @tool(ToolAnnotations(destructive_hint=False, idempotent_hint=True, open_world_hint=False))
    def write_input(path: str, yaml: str, force: bool = False) -> dict[str, Any]:
        """Save an input file in the workspace after validating it (smoke level); refused while
        errors remain, unless force is true."""
        full = _resolve(state.workspace, path)
        report = _validate(yaml, state.catalog, os.path.dirname(full), level="smoke")
        if not report["ok"] and not force:
            return {"written": False, "reason": "validation failed", "validation": report}
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, "w") as f:
            f.write(yaml)
        return {"written": True, "path": full, "validation": report}

    @tool(ToolAnnotations(destructive_hint=False, open_world_hint=False))
    def run_input(path: str, nprocs: int = 1, omp_threads: Optional[int] = None,
                  max_iteration: Optional[int] = None) -> dict[str, Any]:
        """Run an input file of the workspace in the background (output in <path>.log next to it).
        max_iteration overrides the number of steps. Poll with job_status."""
        full = _resolve(state.workspace, path)
        if not os.path.isfile(full):
            raise FileNotFoundError(f"no file '{path}' in the workspace")
        note = None
        if state.max_run_steps:
            with open(full) as f:
                g = (_yaml.safe_load(f) or {}).get("global")
            start = int((g.get("start_iteration") if isinstance(g, dict) else 0) or 0)
            capped = start + state.max_run_steps
            if not max_iteration or max_iteration > capped:
                max_iteration = capped
                note = f"this server caps runs to {state.max_run_steps} step(s): max_iteration set to {capped}"
        overrides = {"global": {"max_iteration": max_iteration}} if max_iteration else None
        job = state.jobs.start(full, nprocs=nprocs, overrides=overrides, omp_threads=omp_threads)
        status = state.jobs.status(job, tail=5)
        return dict(status, note=note) if note else status

    @tool(ro)
    def job_status(job_id: str, tail: int = 40) -> dict[str, Any]:
        """State (running, finished, failed), exit code and the end of the log of a run."""
        return state.jobs.status(job_id, tail)

    @tool(ToolAnnotations(destructive_hint=True, open_world_hint=False))
    def stop_job(job_id: str) -> dict[str, Any]:
        """Stop a running job."""
        return state.jobs.stop(job_id)

    @tool(ro)
    def index_examples() -> dict[str, Any]:
        """Build every shipped example once (about a minute, cached per install) so that
        search_examples can tell which examples work with this build."""
        built, failed = state.corpus.index_build_status(
            lambda text, d: check_build(PACKAGE, text, d, timeout=120))
        return {"build": built, "fail": failed}

    @server.resource("exastamp://workflow", mime_type="text/plain")
    def workflow_resource() -> str:
        """Workflow and rules for writing exaStamp inputs."""
        return WORKFLOW

    @server.resource("exastamp://block/{name}", mime_type="application/json")
    def block_resource(name: str) -> dict[str, Any]:
        """Curated schema of an input block."""
        if name not in SCHEMAS:
            raise ValueError(f"no schema for '{name}'")
        return SCHEMAS[name]

    @server.prompt()
    def build_exastamp_input(request: str) -> str:
        """Guided workflow to turn a simulation request into a validated exaStamp input."""
        return (f"Write an exaStamp input file for this request:\n\n{request}\n\n{WORKFLOW}\n"
                "Finish with write_input, then give the file path, the validation summary and the notes "
                "on derived choices.")

    server.state = state
    return server


def main(argv=None):
    parser = argparse.ArgumentParser(prog="exastamp-mcp", description=__doc__.split("\n\n")[0])
    parser.add_argument("--workspace", default=os.environ.get("EXASTAMP_MCP_WORKSPACE", "."))
    parser.add_argument("--doc-dir", default=os.environ.get("EXASTAMP_DOC_DIR"),
                        help="doc_exaStamp/docs directory (markdown), for search_docs")
    parser.add_argument("--index-examples", action="store_true",
                        help="build every example once, cache which ones work, and exit")
    parser.add_argument("--max-run-steps", type=int, default=os.environ.get("EXASTAMP_MCP_MAX_RUN_STEPS"),
                        help="cap every run_input to this number of time steps (evaluations)")
    args = parser.parse_args(argv)
    server = build_server(args.workspace, args.doc_dir, args.max_run_steps)
    if args.index_examples:
        built, failed = server.state.corpus.index_build_status(
            lambda text, d: check_build(PACKAGE, text, d, timeout=120))
        print(f"{built} examples build, {failed} do not", file=sys.stderr)
        return 0
    server.run("stdio")
    return 0
