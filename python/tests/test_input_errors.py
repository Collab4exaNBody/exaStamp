"""Input mistakes must give an explicit error, not be ignored or produce NaN."""

import pytest

import pyexastamp as xs

EAM_SAMPLE = "regression_new/potentials/eam/eam_vniitf/single_specy.msp"
STRICT = {"debug": {"strict_slots": True}}


def _run(tmp_path, monkeypatch, overrides, config=None):
    monkeypatch.chdir(tmp_path)
    base = {"global": {"max_iteration": 1}}
    base.update(overrides)
    ctx = xs.init(xs.command_line(xs.find_example(EAM_SAMPLE), program="test", overrides=base, config=config))
    try:
        xs.run(ctx)
    finally:
        xs.end(ctx)


def test_default_config_is_strict_clean(tmp_path, monkeypatch):
    _run(tmp_path, monkeypatch, {}, STRICT)


def test_unknown_parameter_strict(tmp_path, monkeypatch):
    with pytest.raises(xs.OnikaError, match="has no parameter 'Tstrat'.*Valid parameters are:.* Tstart"):
        _run(tmp_path, monkeypatch, {"numerical_scheme": "verlet_lnvt", "langevin_thermostat": {"Tstrat": "300 K"}}, STRICT)


def test_unknown_unit(tmp_path, monkeypatch):
    with pytest.raises(xs.OnikaError, match="unknown unit 'pss'.*Valid units are:.* ps"):
        _run(tmp_path, monkeypatch, {"global": {"max_iteration": 1, "dt": "1.0e-3 pss"}})


def test_lattice_unknown_type(tmp_path, monkeypatch):
    setup = [
        {"particle_types": {"particle_type_map": {"Sn": 0}}},
        {"domain": {"cell_size": "6.0 ang", "grid_dims": [10, 10, 10],
                    "bounds": [["0 ang", "0 ang", "0 ang"], ["60 ang", "60 ang", "60 ang"]],
                    "periodic": [True, True, True]}},
        "init_rcb_grid",
        {"lattice": {"structure": "BCC", "types": ["Sn", "Xx"], "size": ["3 ang", "3 ang", "3 ang"]}},
    ]
    with pytest.raises(xs.OnikaError, match="type 'Xx' is not a declared particle type"):
        _run(tmp_path, monkeypatch, {"setup_system": setup})
