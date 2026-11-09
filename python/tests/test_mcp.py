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
    ("numerical_scheme: verlet_lnvt", "numerical_scheme: verlet_nhnpt\ngrid_flavor: grid_flavor_multimat", "grid_flavor", "no virial"),
    ("- Ta: { mass", "- Tantalum: { mass", "setup_system[2].lattice.types", "not a declared species"),
    ("numerical_scheme: verlet_lnvt", "numerical_scheme: verlet_nhnpt\ninit_nose_hoover: { algo: NPT, Tstart: 300 K }",
     "includes", "Nose-Hoover configuration"),
    ("numerical_scheme: verlet_lnvt", "includes: [ config_nose_hoover.msp ]\nnumerical_scheme: verlet_nhnpt\ninit_nose_hoover: { Tstart: 300 K }",
     "init_nose_hoover.algo", "algo: NPT"),
    ("numerical_scheme: verlet_lnvt", "includes: [ config_nose_hoover.msp ]\nnumerical_scheme: verlet_nhnpt\n"
     "init_nose_hoover: { algo: NPT, Tstart: 300 K, Pstart: 1 GPa }", "init_nose_hoover.Pstart", "plain number"),
])
def test_static_rule(catalog, old, new, path, message):
    bad = EXAMPLE.replace(old, new)
    assert bad != EXAMPLE
    errors = _errors(bad, catalog)
    assert any(p == path and message in m for p, m in errors), errors


def test_doc_sections(tmp_path):
    from exastamp_mcp.corpus import Corpus
    (tmp_path / "sw.md").write_text("# **Stillinger-Weber**\n\n## **Description**\n\nThree-body silicon.\n\n"
                                    "## Parameters\n\nfile: parameter file\n")
    hits = Corpus(str(tmp_path), str(tmp_path)).search_docs("stillinger silicon")
    # the page title stays with its first section instead of being an empty hit
    assert hits[0]["section"] == "Stillinger-Weber > Description"
    assert hits[0]["page_sections"] == ["Stillinger-Weber > Description", "Parameters"]
    assert "get_doc('sw.md')" in hits[0]["more"]


def test_includes_are_merged(catalog, tmp_path):
    (tmp_path / "base.yaml").write_text(EXAMPLE.split("eam_alloy_force:")[0])
    rest = "includes:\n  - base.yaml\neam_alloy_force:" + EXAMPLE.split("eam_alloy_force:")[1]
    assert [i for i in check_input(rest, catalog, str(tmp_path)) if i["level"] != "info"] == []


def test_batches_without_false_notes(catalog, tmp_path):
    # '+batch' appends, batches defined by an include or used by another batch, species from a generator batch
    (tmp_path / "base.yaml").write_text("my_included_batch:\n  - nop\n")
    species = EXAMPLE.split("setup_system:")[0].replace("species:\n", "")
    text = ("includes:\n  - base.yaml\n" + "+init_epilog:\n  - nop\n"
            + "my_generator:\n  - species:\n      species:\n" + species.replace("  - ", "        - ")
            + "setup_system:" + EXAMPLE.split("setup_system:")[1].replace("  - init_rcb_grid\n", "  - init_rcb_grid\n  - my_generator\n"))
    issues = [(i["path"], i["message"]) for i in check_input(text, catalog, str(tmp_path)) if i["level"] != "error"]
    assert issues == []
    # an unused batch is still noted
    assert any(i["path"] == "my_unused" for i in check_input(text + "my_unused:\n  - nop\n", catalog, str(tmp_path)))


def test_observed_parameters(catalog):
    from exastamp_mcp.corpus import Corpus
    corpus = Corpus(catalog.paths["XSP_DATA_DIR"])
    names = set(catalog.names())
    ljwolf = corpus.describe_observed("ljwolf_compute_force", names)
    assert "parameters.cw.alpha" in ljwolf["parameters"]
    multimat = corpus.describe_observed("compute_force_pair_multimat", names)["parameters"]
    assert "energy" in multimat["potentials[potential=lj].parameters.epsilon"]

    observed = corpus.observed_parameters(names)
    def soft(text):
        return [(i["path"], i["message"]) for i in check_input(text, catalog, None, observed) if "shipped example" in i["message"]]
    wolf = "ljwolf_compute_force:\n  rcut: 12.5 ang\n  parameters:\n    lj: { epsilon: 0.01 kcal/mol, sigma: 3.8 ang }\n" \
           "    WOLF: { alpha: 0.2 1/ang, rc: 12.5 ang }\n"
    assert soft(wolf.replace("WOLF", "cw")) == []
    assert [p for p, _ in soft(wolf.replace("WOLF", "coul_wolf"))] == ["ljwolf_compute_force.parameters.coul_wolf"]
    pair = "compute_force_pair_multimat:\n  potentials:\n" \
           "    - { type_a: Cu, type_b: Cu, potential: lj, rcut: 6 ang, parameters: { epsilon: EPS, sigma: 2.3 ang } }\n"
    assert soft(pair.replace("EPS", "0.4 eV")) == [] and soft(pair.replace("EPS", "0.2 kcal/mol")) == []
    assert [p for p, _ in soft(pair.replace("EPS", "0.4 ang"))] == ["compute_force_pair_multimat.potentials[0].parameters.epsilon"]
    # maps keyed by data (species names) are not checked key by key
    types = "setup_system:\n  - particle_types:\n      particle_type_map: { Xx: 0, Yy: 1 }\n"
    assert soft(types) == []


@pytest.mark.parametrize("text,path,message", [
    ("species:\n  - W: { mass: 183.84 Da, z: 74, charge: 0 e- }\ncompute_force:\n  - snap_force_fp64:\n"
     "      parameters: { nt: 1, param: Ta06A.snapparam, coef: Ta06A.snapcoeff }\n",
     "snap_force_fp64.parameters.coef", "species 'W' is not in Ta06A.snapcoeff"),
    ("species:\n  - Ta: { mass: 180.95 Da, z: 73, charge: 0 e- }\ncompute_force:\n  - snap_force_fp64:\n"
     "      parameters: { nt: 1, param: Ta06A.snapparam, coef: Ta6A.snapcoeff }\n",
     "snap_force_fp64.parameters.coef", "not found"),
    ("species:\n  - Ge: { mass: 72.63 Da, z: 32, charge: 0 e- }\ninit_parameters:\n  - species\n"
     "  - sw_init:\n      parameters: { file: Si.sw }\n", "sw_init.parameters.file", "species 'Ge' is not in Si.sw"),
    ("species:\n  - Cu: { mass: 63.546 Da, z: 29, charge: 0 e- }\ninit_parameters:\n  - species\n"
     "  - pace_init:\n      parameters: { coef: CH_ACE_Willman_JCP2024.yace }\n", "pace_init.parameters.coef", "species 'Cu' is not in"),
])
def test_potential_files(catalog, text, path, message):
    issues = [(i["path"], i["message"]) for i in check_input(text, catalog) if i["level"] == "error"]
    assert any(p == path and message in m for p, m in issues), issues


def test_missing_operators(tmp_path):
    import json
    from exastamp_mcp.corpus import Corpus
    (tmp_path / "regression_new").mkdir()
    example = tmp_path / "regression_new" / "rsa.msp"
    example.write_text("setup_system:\n  - init_rsa\n")
    status = {"regression_new/rsa.msp": {"mtime": int(example.stat().st_mtime), "builds": False,
              "error": "YAML node for batch simulation.setup_system.init_rsa is not a list but does not contain a 'body' entry"}}
    (tmp_path / "status.json").write_text(json.dumps(status))
    corpus = Corpus(str(tmp_path))
    corpus.load_build_status(str(tmp_path / "status.json"))
    assert corpus.missing_operators() == {"init_rsa": 1}


def test_potential_file_info(catalog):
    from exastamp_mcp.data import potential_file_info
    info = {os.path.basename(f): potential_file_info(f) for f in map(catalog.find_data_file,
            ["WBe_Wood_PRB2019.snapcoeff", "GaN.sw", "Ta_param.pod", "pot.almtp", "Cu-PBE-core-rep.ace", "Ta06A.snapparam"])}
    assert info["WBe_Wood_PRB2019.snapcoeff"]["elements"] == ["W", "Be"]
    assert info["GaN.sw"]["elements"] == ["Ga", "N"]
    assert info["Ta_param.pod"]["elements"] == ["Ta"] and info["pot.almtp"]["species_count"] == 1
    assert info["Cu-PBE-core-rep.ace"]["elements"] == ["Cu"] and info["Ta06A.snapparam"]["rcutfac"] > 4


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
            # defaults exported, wiring slots hidden, parameters seen in the examples
            r = (await c.call_tool("describe_block", {"name": "gaussian_noise_v"})).structured_content
            params = {p["name"]: p for p in r["operator_parameters"]}
            assert params["sigma"]["default"] == 1 and params["ghost"]["default"] is False
            assert "grid" not in params and "domain" not in params and "region" in params
            assert "sigma" in r["observed_parameters"]["parameters"]
            r = await c.call_tool("build_input", {"spec": SPECS["ar_lj_nve"][0]})
            assert r.structured_content["ok"], r.structured_content
            r = await c.call_tool("write_input", {"path": "ar/ar.msp", "yaml": r.structured_content["yaml"]})
            assert r.structured_content["written"]
            assert os.path.isfile(tmp_path / "ar" / "ar.msp")
            r = await c.call_tool("write_input", {"path": "../outside.msp", "yaml": "species: []"})
            assert r.is_error and "outside the workspace" in r.content[0].text

    asyncio.run(run())


def test_run_input_step_cap(tmp_path):
    from mcp import Client
    from exastamp_mcp.server import build_server

    async def run():
        async with Client(build_server(str(tmp_path), max_run_steps=1)) as c:
            r = await c.call_tool("build_input", {"spec": SPECS["ar_lj_nve"][0]})
            await c.call_tool("write_input", {"path": "ar.msp", "yaml": r.structured_content["yaml"]})
            r = await c.call_tool("run_input", {"path": "ar.msp"})
            assert "max_iteration set to 1" in r.structured_content["note"]
            job = r.structured_content["job_id"]
            for _ in range(600):
                status = (await c.call_tool("job_status", {"job_id": job})).structured_content
                if status["state"] != "running":
                    break
                await asyncio.sleep(0.5)
            assert status["state"] == "finished", status

    asyncio.run(run())
