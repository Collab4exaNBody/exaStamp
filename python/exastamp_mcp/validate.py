"""
Layered validation of an exaStamp input, cheapest layer first:

1. ``static``: checks.check_input (milliseconds, no exaStamp process);
2. ``build``: the simulation graph is built in strict mode in a worker;
3. ``smoke``: a few time steps are run with outputs disabled, then the
   thermodynamic state is checked (finite values, no temperature explosion,
   expected number of atoms for a lattice).

A layer runs only when the previous one found no error. Every issue is
``{"stage", "level", "path", "message", "suggestion"?, "doc"?}``.
"""

import math
import re

import yaml

from onika_mcp.smoke import check_build, smoke_run
from onika_mcp.units import UnitError, parse_quantity

from . import PACKAGE
from .checks import check_input, has_errors
from .knowledge import SCHEMAS

__all__ = ["validate_input", "LEVELS"]

LEVELS = ("static", "build", "smoke")


def _kelvin(q):
    try:
        v, _, _ = parse_quantity(q)
        return v
    except (UnitError, TypeError):
        return None


def _length(q):
    try:
        v, _, has_units = parse_quantity(q)
        return v * 1e10 if has_units else v
    except (UnitError, TypeError):
        return None


def _setup(doc):
    for item in doc.get("setup_system") or []:
        if isinstance(item, dict) and len(item) == 1:
            name, params = next(iter(item.items()))
            yield name, params if isinstance(params, dict) else {}
        elif isinstance(item, str):
            yield item, {}


def _reference_temperature(doc):
    temps = []
    for name, params in _setup(doc):
        if name == "init_temperature_new" and "T" in params:
            temps.append(_kelvin(params["T"]))
    for block in ("langevin_thermostat", "berendsen_thermostat"):
        p = doc.get(block)
        if isinstance(p, dict):
            temps += [_kelvin(p.get(k)) for k in ("T", "Tstart", "Tstop") if k in p]
    temps = [t for t in temps if t is not None]
    return max(temps) if temps else None


def _expected_atoms(doc):
    """Atoms a lattice-only setup must generate, or None when it cannot be known."""
    domain, count = None, None
    for name, params in _setup(doc):
        if name == "domain":
            domain = params
        elif name == "lattice":
            if count is not None or any(k in params for k in ("region", "void_mode", "void_porosity", "void_radius")):
                return None
            n = SCHEMAS["lattice"]["structure_types_count"].get(params.get("structure"))
            size = params.get("size")
            if not (n and domain and domain.get("bounds") and isinstance(size, list) and len(size) == 3):
                return None
            if "xform" in domain or not all(domain.get("periodic") or [False]):
                return None
            try:
                box = [_length(h) - _length(l) for l, h in zip(*domain["bounds"])]
                cells = [b / _length(s) for b, s in zip(box, size)]
            except (TypeError, ValueError, ZeroDivisionError):
                return None
            if any(abs(c - round(c)) > 1e-6 for c in cells):
                return None
            count = n * math.prod(round(c) for c in cells)
        elif name.startswith(("read_", "replicate")) or "molecule" in name:
            return None
    return count


def _locate(message, doc):
    """Top-level block or setup_system item an error message points at."""
    names = re.findall(r"'([A-Za-z_][\w\-]*)'", message) + re.findall(r"operator ([A-Za-z_]\w*)", message)
    for i, (name, _) in enumerate(_setup(doc)):
        if name in names:
            return f"setup_system[{i}].{name}"
    for key in doc:
        if key in names:
            return key
    return ""


def _run_issues(stage, result, doc):
    issues = []
    if not result.get("ok"):
        msg = result.get("error") or "failed"
        issue = {"stage": stage, "level": "error", "path": _locate(msg, doc), "message": msg}
        if "timeout" in msg:
            issue["suggestion"] = "the system may be too large for a smoke run; check with a smaller box first"
        elif "crashed" in msg:
            issue["suggestion"] = "see log_tail; a crash usually comes from an inconsistent setup (empty domain, overlapping atoms)"
        issues.append(issue)
    for w in result.get("warnings", []):
        issues.append({"stage": stage, "level": "warning", "path": _locate(w, doc), "message": w})
    return issues


def _sanity(thermo, doc):
    issues = []
    if not thermo:
        return [{"stage": "smoke", "level": "warning", "path": "",
                 "message": "no thermodynamic state was computed during the smoke run"}]
    for k, v in thermo.items():
        if isinstance(v, float) and not math.isfinite(v):
            issues.append({"stage": "smoke", "level": "error", "path": "",
                           "message": f"{k} is {v} after the smoke run: the simulation is numerically unstable",
                           "suggestion": "check units, the time step and overlapping atoms (lattice parameter, duplicated layers)"})
    t, t_ref = thermo.get("temperature"), _reference_temperature(doc)
    if t is not None and math.isfinite(t):
        limit = max(10 * (t_ref or 300.0), (t_ref or 300.0) + 5000.0)
        if t > limit:
            issues.append({"stage": "smoke", "level": "error", "path": "",
                           "message": f"temperature {t:.3g} K after a few steps (expected about {t_ref or 0:g} K): atoms overlap or the time step is too large",
                           "suggestion": "check the lattice parameter against the potential, duplicated atoms at the box faces (add a lattice shift), and dt"})
    expected, got = _expected_atoms(doc), thermo.get("particle_count")
    if expected and got is not None and int(got) != expected:
        issues.append({"stage": "smoke", "level": "error", "path": "setup_system",
                       "message": f"{int(got)} atoms were generated, {expected} expected from the box and the lattice",
                       "suggestion": "atoms on the upper box faces are duplicated or lost: add a small lattice shift, e.g. 0.001 x the lattice parameter"})
    return issues


def validate_input(text, catalog, workdir=None, level="smoke", steps=2, timeout=300, observed=None):
    if level not in LEVELS:
        raise ValueError(f"level must be one of {LEVELS}")
    issues = [dict(i, stage="static") for i in check_input(text, catalog, workdir, observed)]
    report = {"ok": False, "stage_reached": "static", "issues": issues}
    if has_errors(issues) or level == "static":
        report["ok"] = not has_errors(issues)
        return report
    doc = yaml.safe_load(text) or {}

    if level == "build":
        result = check_build(PACKAGE, text, workdir, timeout=timeout)
        report["stage_reached"] = "build"
    else:
        start = (doc.get("global") or {}).get("start_iteration", 0) if isinstance(doc.get("global"), dict) else 0
        overrides = {"global": {"max_iteration": int(start or 0) + steps, "simulation_thermostate_file_frequency": -1},
                     "write_snapshot": "nop", "write_restart": "nop"}
        result = smoke_run(PACKAGE, text, workdir, overrides=overrides, timeout=timeout, thermo_reader="read_thermo")
        report["stage_reached"] = "smoke"
    issues += _run_issues(report["stage_reached"], result, doc)
    if result.get("ok") and level == "smoke":
        issues += _sanity(result.get("thermo"), doc)
        report["thermo_after_smoke"] = result.get("thermo")
    report["elapsed_s"] = result.get("elapsed")
    if not result.get("ok"):
        report["log_tail"] = (result.get("log") or "")[-2000:]
    report["ok"] = not has_errors(issues)
    return report
