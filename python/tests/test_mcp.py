"""exaStamp MCP server: checker rules, renderer, validation layers and tools."""

import asyncio
import os

import pytest

mcp = pytest.importorskip("mcp", reason="the MCP server needs the 'mcp' package")

from onika_mcp.catalog import Catalog
from onika_mcp.units import check_quantity
from onika_mcp.worker import run_worker
from exastamp_mcp.checks import check_input, has_errors
from exastamp_mcp.knowledge import SCHEMAS
from exastamp_mcp.render import RenderError, render_input
from exastamp_mcp.validate import validate_input

# doc_exaStamp Starter Pack, 8_example.md
EXAMPLE = """\
species:
  - Ta: { mass: 180.95 Da, z: 73, charge: 0 e- }
setup_system:
  - domain:
      cell_size: 3.3 ang
      grid_dims: [ 20, 20, 20 ]
      bounds: [ [ 0 ang, 0 ang, 0 ang ], [ 66 ang, 66 ang, 66 ang ] ]
      periodic: [ true, true, true ]
      expandable: false
  - init_rcb_grid
  - lattice:
      structure: BCC
      types: [ Ta, Ta ]
      size: [ 3.3 ang, 3.3 ang, 3.3 ang ]
  - init_temperature_new:
      T: 600 K
eam_alloy_force:
  rcut: 5.30 ang
  parameters:
    file: "Ta1_Ravelo_2013.eam.alloy"
compute_force:
  - eam_alloy_force
numerical_scheme: verlet_lnvt
langevin_thermostat:
  T: 300 K
  gamma: 0.1 ps^-1
global:
  max_iteration: 1000
  dt: 1.0e-3 ps
  rcut_inc: 1.0 ang
  simulation_thermostate_screen_frequency: 10
  simulation_thermostate_file_frequency: 10
  simulation_restart_frequency: 1000
  simulation_snapshot_frequency: 100
  simulation_analysis_frequency: 0
write_snapshot: write_snapshot_xyz
"""


@pytest.fixture(scope="module")
def catalog():
    return Catalog("pyexastamp")


def _errors(text, catalog):
    return [(i["path"], i["message"]) for i in check_input(text, catalog) if i["level"] == "error"]


# --------------------------------------------------------------------------- units

@pytest.mark.parametrize("text,dim,ok", [
    ("3.3 ang", "length", True), ("5 A", "length", False), ("1.0e-3 pss", "time", False),
    ("1.0e-3 ps", "length", False), ("0.1 ps^-1", "inverse_time", True), ("50. m/s", "velocity", True),
    ("1.0e-3 fs.", "time", False), ("12.011 g/mol", "mass", True), ("0.2 kcal/mol", "energy", True), (300, "temperature", True),
])
def test_quantities(text, dim, ok):
    assert (check_quantity(text, dim) is None) == ok


# --------------------------------------------------------------------------- static checks

def test_docs_example_is_clean(catalog):
    assert _errors(EXAMPLE, catalog) == []


def test_schema_snippets_are_clean(catalog):
    for name, schema in SCHEMAS.items():
        if "snippet" in schema:
            assert _errors(schema["snippet"], catalog) == [], name


@pytest.mark.parametrize("old,new,path,message", [
    ("T: 300 K", "Tstrat: 300 K", "langevin_thermostat.Tstrat", "unknown parameter 'Tstrat'"),
    ("dt: 1.0e-3 ps", "dt: 1.0e-3 pss", "global.dt", "unknown unit 'pss'"),
    ("size: [ 3.3 ang, 3.3 ang, 3.3 ang ]", "size: [ 3.3 A, 3.3 ang, 3.3 ang ]", "setup_system[2].lattice.size", "unknown unit 'A'"),
    ("grid_dims: [ 20, 20, 20 ]", "grid_dims: [ 20, 20, 21 ]", "setup_system[0].domain.bounds", "along z"),
    ("types: [ Ta, Ta ]", "types: [ Ta ]", "setup_system[2].lattice.types", "BCC needs 2 types"),
    ("  - init_rcb_grid\n", "", "setup_system[1].lattice", "init_rcb_grid must be called"),
    ("Ta1_Ravelo_2013.eam.alloy", "Ta1_Ravelo.eam.alloy", "eam_alloy_force.parameters.file", "not found"),
    ("numerical_scheme: verlet_lnvt", "numerical_scheme: verlet_lnvtt", "numerical_scheme", "is not one of"),
    ("- Ta: { mass", "- Tantalum: { mass", "setup_system[2].lattice.types", "not a declared species"),
])
def test_static_rule(catalog, old, new, path, message):
    bad = EXAMPLE.replace(old, new)
    assert bad != EXAMPLE
    errors = _errors(bad, catalog)
    assert any(p == path and message in m for p, m in errors), errors


def test_includes_are_merged(catalog, tmp_path):
    (tmp_path / "base.yaml").write_text(EXAMPLE.split("eam_alloy_force:")[0])
    rest = "includes:\n  - base.yaml\neam_alloy_force:" + EXAMPLE.split("eam_alloy_force:")[1]
    assert [i for i in check_input(rest, catalog, str(tmp_path)) if i["level"] != "info"] == []


def test_parent_include_builds(catalog, tmp_path):
    # physics_regression layout: <case>/run.msp includes ../common.yaml
    (tmp_path / "base.yaml").write_text(EXAMPLE.split("eam_alloy_force:")[0])
    (tmp_path / "case").mkdir()
    rest = "includes:\n  - ../base.yaml\neam_alloy_force:" + EXAMPLE.split("eam_alloy_force:")[1]
    report = validate_input(rest, catalog, str(tmp_path / "case"), level="build")
    assert report["ok"], report["issues"]


# --------------------------------------------------------------------------- renderer + validation

SPECS = {
    "ta_bcc_langevin": ({"species": [{"name": "Ta"}],
                         "system": {"kind": "lattice", "structure": "BCC", "lattice_parameter": "3.3026 ang", "repeats": [10, 10, 10]},
                         "potential": {"kind": "eam_alloy", "file": "Ta1_Ravelo_2013.eam.alloy"},
                         "ensemble": {"kind": "langevin", "temperature": "300 K"},
                         "outputs": {"snapshot": "xyz", "snapshot_every": 100}}, 2000),
    "cu_fcc_berendsen": ({"species": [{"name": "Cu"}],
                          "system": {"kind": "lattice", "structure": "FCC", "lattice_parameter": "3.615 ang", "repeats": [8, 8, 8]},
                          "potential": {"kind": "eam_alloy", "file": "Cu.eam.alloy"},
                          "ensemble": {"kind": "berendsen", "temperature_start": "100 K", "temperature_stop": "500 K"}}, 2048),
    "ar_lj_nve": ({"species": [{"name": "Ar"}],
                   "system": {"kind": "lattice", "structure": "FCC", "lattice_parameter": "5.26 ang", "repeats": [6, 6, 6]},
                   "potential": {"kind": "lj", "epsilon": "0.0104 eV", "sigma": "3.4 ang", "rcut": "8.5 ang"},
                   "initial_temperature": "80 K"}, 864),
    "alcu_alloy": ({"species": [{"name": "Al"}, {"name": "Cu"}],
                    "system": {"kind": "lattice", "structure": "FCC", "lattice_parameter": "3.8 ang", "repeats": [7, 7, 7],
                               "types": ["Al", "Cu", "Cu", "Al"]},
                    "potential": {"kind": "eam_alloy", "file": "AlCu.eam.alloy"},
                    "ensemble": {"kind": "langevin", "temperature": "600 K"}}, 1372),
    "sn_2bct_xform": ({"species": [{"name": "Sn"}],
                       "system": {"kind": "lattice", "structure": "2BCT", "lattice_parameter": ["5.83 ang", "5.83 ang", "3.18 ang"],
                                  "repeats": [8, 8, 14]},
                       "potential": {"kind": "eam_alloy", "file": "Sn.eam.alloy"}, "initial_temperature": "300 K"}, 3584),
}


@pytest.mark.parametrize("name", sorted(SPECS))
def test_rendered_input_runs(catalog, name, tmp_path):
    spec, atoms = SPECS[name]
    text, notes = render_input(spec, catalog.data_dirs, catalog.find_data_file)
    report = validate_input(text, catalog, str(tmp_path))
    assert report["ok"], report["issues"]
    assert report["thermo_after_smoke"]["particle_count"] == atoms


@pytest.mark.parametrize("change,message", [
    ({"species": [{"name": "Tantalum", "mass": "180.95 Da"}]}, "not in Ta1_Ravelo_2013.eam.alloy"),
    ({"potential": {"kind": "eam_alloy", "file": "Ta_missing.eam.alloy"}}, "not found"),
    ({"system": {"kind": "lattice", "structure": "BCC", "lattice_parameter": "3.3 ang", "repeats": [4, 4, 4], "types": ["Ta"]}}, "needs 2 types"),
    ({"ensemble": {"kind": "langevin"}}, "needs temperature"),
])
def test_render_errors(catalog, change, message):
    spec = dict(SPECS["ta_bcc_langevin"][0], **change)
    with pytest.raises(RenderError, match=message):
        render_input(spec, catalog.data_dirs, catalog.find_data_file)


def test_smoke_detects_duplicated_atoms(catalog, tmp_path):
    # 20 x 3.3026 rounds just below 66.052: without a lattice shift a layer is generated twice
    text = (EXAMPLE.replace("3.3 ang, 3.3 ang, 3.3 ang ]", "3.3026 ang, 3.3026 ang, 3.3026 ang ]")
                   .replace("cell_size: 3.3 ang", "cell_size: 3.3026 ang")
                   .replace("66 ang, 66 ang, 66 ang", "66.052 ang, 66.052 ang, 66.052 ang"))
    static = check_input(text, catalog)
    assert any("duplicate a layer" in i["message"] for i in static)
    report = validate_input(text, catalog, str(tmp_path))
    assert not report["ok"]
    assert any("atoms were generated" in i["message"] or "temperature" in i["message"] for i in report["issues"])


def test_build_layer_reports_runtime_error(catalog, tmp_path):
    # passes the static checks (the operator exists) but fails when the graph is built
    text = EXAMPLE + "chunk_neighbors:\n  chunk_size: 2\n"
    report = validate_input(text, catalog, str(tmp_path), level="build")
    assert report["stage_reached"] in ("static", "build") and not report["ok"]


def test_worker_timeout_and_bad_package():
    assert "timeout" in run_worker("pyexastamp", "catalog", timeout=0.01)["error"]
    r = run_worker("no_such_package_xyz", "env", timeout=60)
    assert not r["ok"] and "no_such_package_xyz" in r["error"]


# --------------------------------------------------------------------------- MCP tools

def test_server_tools(tmp_path):
    from mcp import Client
    from exastamp_mcp.server import build_server

    async def run():
        async with Client(build_server(str(tmp_path))) as c:
            names = {t.name for t in (await c.list_tools()).tools}
            assert {"build_input", "validate_input", "describe_block", "search_examples", "write_input", "run_input"} <= names
            r = await c.call_tool("describe_block", {"name": "langevin_thermostatt"})
            assert r.is_error and "did you mean langevin_thermostat" in r.content[0].text
            r = await c.call_tool("build_input", {"spec": SPECS["ar_lj_nve"][0]})
            assert r.structured_content["ok"], r.structured_content
            r = await c.call_tool("write_input", {"path": "ar/ar.msp", "yaml": r.structured_content["yaml"]})
            assert r.structured_content["written"]
            assert os.path.isfile(tmp_path / "ar" / "ar.msp")
            r = await c.call_tool("write_input", {"path": "../outside.msp", "yaml": "species: []"})
            assert r.is_error and "outside the workspace" in r.content[0].text

    asyncio.run(run())
