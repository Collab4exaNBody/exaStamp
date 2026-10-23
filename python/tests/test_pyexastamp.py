"""pyexastamp: run an EAM sample and read thermodynamics, species and particles."""

import numpy as np
import pytest

import pyexastamp as xs

EAM_SAMPLE = "regression_new/potentials/eam/eam_vniitf/single_specy.msp"


def _start(tmp_path, monkeypatch, **overrides):
    monkeypatch.chdir(tmp_path)  # plugins.db and outputs are written in the working directory
    base = {"global": {"max_iteration": 5}}
    base.update(overrides)
    ctx = xs.init(xs.command_line(xs.find_example(EAM_SAMPLE), program="test", overrides=base))
    assert ctx.error_code < 0
    return ctx


@pytest.fixture(scope="module")
def eam(tmp_path_factory):
    """One short run shared by the tests of this module."""
    with pytest.MonkeyPatch.context() as mp:
        ctx = _start(tmp_path_factory.mktemp("eam"), mp,
                     grid_flavor="grid_flavor_full",  # field set with virial (Mat3d) and cmol (uint64 x4)
                     simulation_epilog=[{"extract_particle_fields": {}}])
        xs.run(ctx)
        yield ctx
        xs.end(ctx)


def test_thermo(eam):
    t = xs.read_thermo(eam)
    assert t["particle_count"] > 0
    for k in ("temperature", "pressure", "total_energy", "kinetic_energy", "potential_energy", "volume"):
        assert np.isfinite(t[k]), k
    assert t["temperature"] > 0
    assert t["potential_energy"] < 0  # cohesive EAM crystal
    assert t["stress"].shape == (3, 3)
    assert t["internal"]["virial"].shape == (3, 3)
    assert t["units"]["temperature"] == "K"


def test_species(eam):
    species = xs.read_species(eam)
    assert [s["name"] for s in species] == ["Sn"]
    assert species[0]["mass"] == pytest.approx(118.71, rel=1e-6)
    assert species[0]["z"] == 50


def test_particles(eam):
    p = xs.particles(eam)
    n = xs.read_thermo(eam)["particle_count"]
    assert p["id"].shape == (n,) and len(np.unique(p["id"])) == n
    assert p["rx"].dtype == np.float64 and p["rx"].shape == (n,)
    assert p["virial"].shape == (n, 3, 3)
    assert p["cmol"].shape == (n, 4) and p["cmol"].dtype == np.uint64
    assert np.all(np.isfinite(p["virial"]))


def test_bad_input_raises(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    bad = tmp_path / "bad.msp"
    bad.write_text("global: { dt: [ unclosed\n")
    with pytest.raises(RuntimeError):
        xs.init(["test", str(bad)])


def test_examples():
    assert EAM_SAMPLE in xs.examples(contains="eam_vniitf")
    assert len(xs.examples("samples")) > 0


RIGID_SAMPLE = "regression_new/rigid_molecule/O2/o2.msp"


@pytest.fixture(scope="module")
def o2(tmp_path_factory):
    """Rigid O2 molecules: grid_flavor_rigidmol with orient/angmom/couple fields."""
    with pytest.MonkeyPatch.context() as mp:
        mp.chdir(tmp_path_factory.mktemp("o2"))
        ctx = xs.init(xs.command_line(xs.find_example(RIGID_SAMPLE), program="test", overrides={
            "global": {"max_iteration": 3},
            "simulation_epilog": [{"extract_particle_fields": {"fields": ["orient", "angmom", "couple"]}}],
        }))
        assert ctx.error_code < 0
        xs.run(ctx)
        yield ctx
        xs.end(ctx)


def test_rigid_molecule_fields(o2):
    p = xs.particles(o2)
    n = len(p["id"])
    assert n > 0
    assert p["orient"].shape == (n, 4) and p["orient"].dtype == np.float64
    assert p["angmom"].shape == (n, 3)
    assert p["couple"].shape == (n, 3)
    # orientations are unit quaternions (w, x, y, z)
    assert np.allclose(np.linalg.norm(p["orient"], axis=1), 1.0, atol=1e-6)


def test_rigid_molecule_species(o2):
    species = xs.read_species(o2)
    rigid = [s for s in species if s.get("rigid_atoms")]
    assert rigid, species
    assert all(a["position"].shape == (3,) for a in rigid[0]["rigid_atoms"])
    assert np.isfinite(xs.read_thermo(o2)["temperature"])
