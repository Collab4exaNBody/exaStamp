"""
Static checks of exaStamp input files (no simulation is run).

Catches, in milliseconds, the mistakes exaStamp only reports at run time or
not at all: unknown parameters, wrong or misspelled units, inconsistent domain
sizes, lattice types that do not match the structure or the species, missing
potential files, thermostat set-up errors.

``check_input(text, catalog, workdir)`` returns a list of issues
``{"level": "error"|"warning"|"info", "path", "message", "suggestion", "doc"}``.
"""

import difflib
import math
import os

import yaml

from onika_mcp.units import DIMENSIONS, check_quantity, parse_quantity, UnitError

from .corpus import is_map_list, is_open, item_key, operator_uses, quantity_dimension
from .data import FILE_PARAMETERS, POTENTIAL_FORMATS, potential_file_info
from .knowledge import SCHEMAS, doc_url

__all__ = ["check_input", "has_errors"]


def has_errors(issues):
    return any(i["level"] == "error" for i in issues)


class _Report:
    def __init__(self, observed=None):
        self.issues = []
        self.observed = observed or {}   # operator -> parameters seen in the examples (corpus)

    def add(self, level, path, message, suggestion=None, block=None):
        issue = {"level": level, "path": path, "message": message}
        if suggestion:
            issue["suggestion"] = suggestion
        doc = SCHEMAS.get(block, {}).get("doc") if block else None
        if doc:
            issue["doc"] = doc_url(doc)
        self.issues.append(issue)


def _closest(word, candidates):
    match = difflib.get_close_matches(str(word), [str(c) for c in candidates], n=1, cutoff=0.6)
    return match[0] if match else None


def _op_name_and_params(item):
    """A batch item is 'name' or {name: params}."""
    if isinstance(item, str):
        return item, None
    if isinstance(item, dict) and len(item) == 1:
        name, params = next(iter(item.items()))
        return name, params
    return None, None


# ---------------------------------------------------------------------------
# values
# ---------------------------------------------------------------------------

def _check_value(rep, path, spec, value, block):
    t = spec.get("type")
    dim = spec.get("dimension")
    if t == "quantity":
        if isinstance(value, bool) or not isinstance(value, (int, float, str)):
            rep.add("error", path, f"expected a quantity like '1.0 ang', got {value!r}", block=block)
            return
        err = check_quantity(value, dim)
        if err:
            rep.add("error", path, err, block=block)
        elif isinstance(value, (int, float)) and dim and dim != "none" and value != 0:
            rep.add("warning", path, f"{value} has no unit: it is read in internal units",
                    suggestion=f"write the unit explicitly, e.g. '{value} {_unit_hint(dim)}'", block=block)
    elif t in ("vec3", "aabb") and dim:
        flat = []
        def _walk(v):
            if isinstance(v, list):
                for x in v:
                    _walk(x)
            else:
                flat.append(v)
        _walk(value)
        expected = 3 if t == "vec3" else 6
        if len(flat) != expected:
            rep.add("error", path, f"expected {expected} values, got {len(flat)}", block=block)
        for v in flat:
            err = check_quantity(v, dim)
            if err:
                rep.add("error", path, err, block=block)
                break
    elif t == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            rep.add("error", path, f"expected an integer, got {value!r}", block=block)
    elif t == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)) and not _is_number(value):
            unit = f" in {spec['doc']}" if spec.get("doc") else ""
            rep.add("error", path, f"expected a plain number{unit} (no unit), got {value!r}", block=block)
    elif t == "bool":
        if not isinstance(value, bool):
            rep.add("error", path, f"expected true or false, got {value!r}", block=block)
    elif t == "ijk":
        if not (isinstance(value, list) and len(value) == 3 and all(isinstance(v, int) and not isinstance(v, bool) for v in value)):
            rep.add("error", path, f"expected 3 integers, got {value!r}", block=block)
    allowed = spec.get("values")
    if allowed and not spec.get("open_values") and isinstance(value, str) and value not in allowed:
        rep.add("error", path, f"'{value}' is not one of {allowed}", suggestion=_closest(value, allowed), block=block)
    sub = spec.get("keys")
    if sub and isinstance(value, dict):
        _check_keys(rep, path, sub, value, block, closed=True)


def _is_number(text):
    """yaml-cpp reads '1.0e9' as a number; PyYAML (YAML 1.1) leaves it a string."""
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


def _unit_hint(dim):
    return {"length": "ang", "time": "ps", "mass": "Da", "temperature": "K", "energy": "eV",
            "charge": "e-", "inverse_time": "ps^-1", "velocity": "m/s"}.get(dim, "")


# keys of an operator node read by onika itself, not by the operator
_STRUCTURAL_KEYS = {"rebind", "condition", "body", "name", "profiling", "parallel", "loop", "unlocked"}


def _check_keys(rep, path, keys_spec, params, block, closed, catalog_slots=None, definition=False):
    if not isinstance(params, dict):
        rep.add("error", path, f"expected a map of parameters, got {type(params).__name__}", block=block)
        return
    for k, v in params.items():
        p = f"{path}.{k}"
        if k in _STRUCTURAL_KEYS:
            continue
        if k in keys_spec:
            _check_value(rep, p, keys_spec[k], v, block)
        elif catalog_slots is not None and k in catalog_slots:
            pass
        elif closed or catalog_slots is not None:
            known = set(keys_spec) | set(catalog_slots or ())
            rep.add("error", p, f"unknown parameter '{k}'", suggestion=_closest(k, known)
                    or f"valid parameters: {', '.join(sorted(known))}", block=block)
    for k, spec in keys_spec.items():
        if spec.get("required") and k not in params:
            if definition:   # top-level definition: the value may come where it is used
                rep.add("warning", f"{path}.{k}", f"'{k}' is required: give it here or where '{block}' is used", block=block)
            else:
                rep.add("error", f"{path}.{k}", f"missing required parameter '{k}'", block=block)


def _check_operator(rep, path, name, params, catalog, definition=False):
    """Parameters of one operator, from its schema and the catalog's slot names."""
    schema = SCHEMAS.get(name, {})
    op = catalog.get(name)
    slots = catalog.slot_names(name) if op else None
    if op is None and name not in SCHEMAS:
        return
    if params is None:
        params = {}
    if not isinstance(params, dict):
        # short form: a scalar or a list given for the operator's main slot
        # (e.g. particle_regions: [...], langevin_thermostat: 300 K)
        return
    if name in rep.observed:
        _check_observed(rep, path, params, rep.observed[name], name)
    if schema.get("open_keys"):
        return
    _check_keys(rep, path, schema.get("keys", {}), params, name, closed=bool(schema) and slots is None,
                catalog_slots=slots, definition=definition)


def _check_observed(rep, path, value, node, block, depth=0):
    """Soft checks against the parameters seen in the shipped examples: nested keys
    never used there (top-level keys are checked against the operator slots), and
    quantities whose dimension differs from every example."""
    if isinstance(value, dict):
        if is_open(node):
            return
        for k, v in value.items():
            child = node["keys"].get(str(k))
            if child is not None:
                _check_observed(rep, f"{path}.{k}", v, child, block, depth + 1)
            elif depth > 0 and node["n"] >= 2 and node["keys"] and k not in _STRUCTURAL_KEYS:
                seen = sorted(node["keys"], key=lambda c: -node["keys"][c]["n"])
                rep.add("warning", f"{path}.{k}", f"'{k}' is never used here in the shipped examples of {block} "
                        f"({node['n']} uses): misspelled or misplaced? The key may be ignored",
                        suggestion=_closest(k, seen) or f"keys used here: {', '.join(seen[:12])}", block=block)
    elif is_map_list(value):
        for i, item in enumerate(value):
            child = node["keys"].get(item_key(item))
            if child is not None:
                _check_observed(rep, f"{path}[{i}]", item, child, block, depth + 1)
    elif node["dims"] and len(node["dims"]) == 1:
        expected = next(iter(node["dims"]))
        for q in value if isinstance(value, list) else [value]:
            dim = quantity_dimension(q)
            # named dimensions only; check_quantity also accepts molar forms (kcal/mol for an energy)
            if dim and expected in DIMENSIONS and check_quantity(q, expected):
                rep.add("warning", path, f"'{q}' is a {dim}; every shipped example gives {'an' if expected[0] in 'aeiou' else 'a'} {expected} here", block=block)
                break


# ---------------------------------------------------------------------------
# blocks
# ---------------------------------------------------------------------------

def _check_batch(rep, path, items, catalog, known_names):
    if not isinstance(items, list):
        rep.add("error", path, "expected a list of operators", block="setup_system")
        return
    for i, item in enumerate(items):
        name, params = _op_name_and_params(item)
        p = f"{path}[{i}]"
        if name is None:
            rep.add("error", p, "each item must be an operator name or a one-key map 'name: {parameters}'")
            continue
        if isinstance(params, dict) and "body" in params and catalog.get(name) is None:
            # inline batch definition: name: { rebind: ..., condition: ..., body: [...] }
            _check_batch(rep, f"{p}.{name}.body", params["body"], catalog, known_names | {name})
            continue
        if catalog.get(name) is None and name not in known_names:
            rep.add("error", f"{p}.{name}", f"unknown operator '{name}'",
                    suggestion=_closest(name, list(catalog.names()) + list(known_names)))
            continue
        _check_operator(rep, f"{p}.{name}", name, params, catalog)


def _batch_items(doc):
    """(batch, operator, params) for every item of every batch (top-level list) of doc."""
    for key, value in doc.items():
        if isinstance(value, list):
            for item in value:
                name, params = _op_name_and_params(item)
                if name:
                    yield key, name, params


def _species_names(doc):
    """Species of the species block, or of a species operator in any batch (e.g. a molecule generator)."""
    lists = [doc.get("species") or []]
    lists += [p.get("species") or [] for b, n, p in _batch_items(doc) if n == "species" and isinstance(p, dict)]
    names = []
    for items in lists:
        for item in items if isinstance(items, list) else []:
            if isinstance(item, dict) and len(item) == 1:
                names.append(next(iter(item)))
    return names


def _setup_ops(doc):
    ops = []
    for item in doc.get("setup_system") or []:
        name, params = _op_name_and_params(item)
        if name:
            ops.append((name, params if isinstance(params, dict) else {}))
    return ops


def _length(q):
    try:
        v, _, has_units = parse_quantity(q)
    except UnitError:
        return None
    return v * 1e10 if has_units else float(q)   # in angstrom


def _check_domain(rep, path, dom):
    if not isinstance(dom, dict):
        return
    cell, dims, bounds = dom.get("cell_size"), dom.get("grid_dims"), dom.get("bounds")
    if cell is None or dims is None or bounds is None:
        return
    c = _length(cell)
    try:
        lo = [_length(x) for x in bounds[0]]
        hi = [_length(x) for x in bounds[1]]
    except (TypeError, IndexError):
        return
    if c is None or None in lo or None in hi or not (isinstance(dims, list) and len(dims) == 3):
        return
    for axis, d, a, b in zip("xyz", dims, lo, hi):
        if not isinstance(d, int) or d <= 0:
            continue
        if abs((b - a) - d * c) > 1e-6 * max(1.0, abs(b - a)):
            rep.add("error", f"{path}.bounds", f"along {axis}: bounds size {b - a:g} ang != grid_dims {d} x cell_size {c:g} ang = {d * c:g} ang",
                    suggestion=f"cell_size {(b - a) / d:.10g} ang (keeps the bounds), or bounds max {a + d * c:g} ang", block="domain")
    if dom.get("expandable") and dom.get("periodic") == [True, True, True]:
        rep.add("warning", f"{path}.expandable", "a fully periodic domain cannot be expandable (exaStamp disables it)", block="domain")


def _check_lattice(rep, path, params, species, ops_before, domain):
    structure = params.get("structure")
    types = params.get("types")
    counts = SCHEMAS["lattice"]["structure_types_count"]
    if isinstance(structure, str) and structure in counts and isinstance(types, list) and len(types) != counts[structure]:
        rep.add("error", f"{path}.types", f"{structure} needs {counts[structure]} types (one per atom of the unit cell), got {len(types)}",
                suggestion=f"types: {[types[0]] * counts[structure] if types else '[...]'}", block="lattice")
    if isinstance(types, list) and species:
        for t in types:
            if t not in species:
                rep.add("error", f"{path}.types", f"type '{t}' is not a declared species ({', '.join(species)})",
                        suggestion=_closest(t, species), block="lattice")
    if "init_rcb_grid" not in ops_before:
        rep.add("error", path, "init_rcb_grid must be called between domain and lattice in setup_system", block="lattice")
    size = params.get("size")
    if isinstance(domain, dict) and isinstance(size, list) and len(size) == 3 and domain.get("bounds") and domain.get("periodic"):
        try:
            box = [_length(h) - _length(l) for l, h in zip(domain["bounds"][0], domain["bounds"][1])]
            cell = [_length(s) for s in size]
        except (TypeError, IndexError):
            return
        xform = domain.get("xform")
        if xform is not None:
            # grid space -> physical space: only diagonal transforms can be checked
            try:
                if any(float(xform[i][j]) != 0.0 for i in range(3) for j in range(3) if i != j):
                    return
                box = [b * float(xform[i][i]) for i, b in enumerate(box)]
            except (TypeError, ValueError, IndexError):
                return
        shift = params.get("shift")
        for axis, b, s, per in zip("xyz", box, cell, domain["periodic"]):
            if per and s and abs(b / s - round(b / s)) <= 1e-6 and round(b / s) * s < b and not shift:
                rep.add("warning", f"{path}.size", f"along periodic {axis}, {round(b / s)} x {s:g} ang rounds just below the box edge ({b:g} ang): the generator may duplicate a layer of atoms",
                        suggestion=f"add shift: [ {1e-3 * cell[0]:.4g} ang, {1e-3 * cell[1]:.4g} ang, {1e-3 * cell[2]:.4g} ang ] to the lattice", block="lattice")
            if per and s and abs(b / s - round(b / s)) > 1e-6:
                rep.add("warning", f"{path}.size", f"along periodic {axis}, the box ({b:g} ang) is not a multiple of the lattice parameter ({s:g} ang): the crystal will have a defect at the boundary",
                        suggestion=f"box size {round(b / s) * s:g} ang", block="lattice")


def _check_potential_files(rep, doc, species, catalog, workdir):
    """Potential files named in the parameters of the operators reading them: the file
    exists, and the species are elements it defines."""
    uses = []
    operator_uses(doc, set(FILE_PARAMETERS), uses)
    for op, params in uses:
        files = params.get("parameters") if isinstance(params.get("parameters"), dict) else {}
        for key in FILE_PARAMETERS[op]:
            fname = files.get(key)
            if not isinstance(fname, str):
                continue
            path = f"{op}.parameters.{key}"
            found = os.path.join(workdir, fname) if workdir and os.path.isfile(os.path.join(workdir, fname)) else None
            found = found or catalog.find_data_file(fname)
            if not found:
                ext = next((e for e in POTENTIAL_FORMATS if fname.endswith(e)), os.path.splitext(fname)[1])
                available = sorted({os.path.basename(f) for d in catalog.data_dirs if os.path.isdir(d)
                                    for f in os.listdir(d) if ext and f.endswith(ext)})
                close = difflib.get_close_matches(fname, available, n=3, cutoff=0.5)
                rep.add("error", path, f"potential file '{fname}' not found in the workspace nor in the data directories",
                        suggestion=(f"one of: {', '.join(close)}" if close else f"available: {', '.join(available) or 'none'}")
                        + " (list_potential_files gives elements and cutoffs)", block=op)
                continue
            try:
                info = potential_file_info(found)
            except (OSError, StopIteration, ValueError, IndexError, KeyError):
                continue
            elements = info.get("elements")
            for s in species if elements else []:
                if s not in elements:
                    rep.add("error", path, f"species '{s}' is not in {os.path.basename(found)} (elements: {', '.join(elements)})",
                            suggestion="name the species like the file's elements", block=op)
            count = info.get("species_count")
            if count and species and len(species) != count:
                rep.add("warning", path, f"{os.path.basename(found)} is fitted for {count} species, the input declares {len(species)}",
                        block=op)


def _check_thermostat(rep, doc):
    scheme = doc.get("numerical_scheme")
    req = SCHEMAS["numerical_scheme"].get("requires", {})
    if isinstance(scheme, str) and scheme in req:
        block = req[scheme]
        params = doc.get(block)
        if not isinstance(params, dict):
            rep.add("error", block, f"numerical_scheme {scheme} needs a '{block}' block with the target temperature",
                    suggestion=SCHEMAS[block].get("snippet"), block=block)
            return
        groups = SCHEMAS[block]["one_of"]
        given = [g for g in groups if all(k in params for k in g)]
        if len(given) != 1:
            rep.add("error", block, f"give exactly one of: {' | '.join('+'.join(g) for g in groups)}", block=block)
    for block in ("langevin_thermostat", "berendsen_thermostat"):
        if block in doc and scheme not in [s for s, b in req.items() if b == block]:
            want = [s for s, b in req.items() if b == block][0]
            rep.add("warning", block, f"'{block}' is defined but numerical_scheme is {scheme or 'verlet_nve (default)'}: it has no effect",
                    suggestion=f"numerical_scheme: {want}", block="numerical_scheme")


def _check_nose_hoover(rep, doc):
    scheme = doc.get("numerical_scheme")
    if scheme not in ("verlet_nhnvt", "verlet_nhnpt"):
        return
    if not any(n == "init_nose_hoover" for _, n, _ in _batch_items(doc)):
        rep.add("error", "includes", f"numerical_scheme {scheme} needs the Nose-Hoover configuration file",
                suggestion="includes:\n  - config_nose_hoover.msp", block="init_nose_hoover")
    params = doc.get("init_nose_hoover")
    if not isinstance(params, dict):
        rep.add("error", "init_nose_hoover", f"numerical_scheme {scheme} needs an 'init_nose_hoover' block (Tstart, ...)",
                suggestion=SCHEMAS["init_nose_hoover"]["snippet"], block="init_nose_hoover")
        return
    algo = "NPT" if scheme == "verlet_nhnpt" else "NVT"
    if params.get("algo", "NVT") != algo:
        rep.add("error", "init_nose_hoover.algo", f"numerical_scheme {scheme} runs with algo: {algo}", block="init_nose_hoover")


# grid flavors without the virial field, and the closest flavor that has it
_NO_VIRIAL = {"grid_flavor_multimat": "grid_flavor_multimat_mechanics", "grid_flavor_minimal": "grid_flavor_multimat_mechanics",
              "grid_flavor_atomic": "grid_flavor_multimat_mechanics", "grid_flavor_atomic_charge": "grid_flavor_atomic_full",
              "grid_flavor_rigidmol": None}
_NPT_SCHEMES = ("verlet_nhnpt", "numerical_scheme_parrinellorahman")


def _check_pressure(rep, doc):
    """A barostat on a grid without virial sees the kinetic pressure only, silently."""
    scheme, flavor = doc.get("numerical_scheme"), doc.get("grid_flavor")
    if scheme in _NPT_SCHEMES and isinstance(flavor, str) and flavor in _NO_VIRIAL:
        better = _NO_VIRIAL[flavor]
        rep.add("error", "grid_flavor", f"{flavor} has no virial field: the pressure is kinetic only and "
                f"{scheme} would drive the wrong pressure",
                suggestion=f"grid_flavor: {better}" if better else "constant pressure is not available for rigid molecules",
                block="grid_flavor")


def _check_global(rep, params, catalog):
    if not isinstance(params, dict):
        return
    keys = SCHEMAS["global"]["keys"]
    known = set(keys) | set((catalog.defaults.get("global") or {}).keys())
    for k, v in params.items():
        if k in keys:
            _check_value(rep, f"global.{k}", keys[k], v, "global")
        elif k not in known:
            rep.add("warning", f"global.{k}", f"'{k}' is not a known global parameter: it defines a new global variable",
                    suggestion=_closest(k, known), block="global")


# ---------------------------------------------------------------------------
# entry point
# ---------------------------------------------------------------------------

def _merge(a, b):
    """onika's merge_nodes: maps merge recursively, '+key' appends to list 'key', anything else is replaced."""
    if not isinstance(a, dict) or not isinstance(b, dict):
        return a if b is None else b
    out = dict(a)
    for k, v in b.items():
        if isinstance(k, str) and k.startswith("+") and isinstance(out.get(k[1:]), list) and isinstance(v, list):
            out[k[1:]] = out[k[1:]] + v
        else:
            out[k] = _merge(out[k], v) if k in out else v
    return out


def _with_includes(rep, doc, base_dir, catalog, seen=None):
    """``doc`` merged on top of the files it includes (recursively)."""
    seen = seen if seen is not None else set()
    includes = doc.get("includes")
    if not includes:
        return doc
    if not isinstance(includes, list):
        rep.add("error", "includes", "includes must be a list of file names")
        return doc
    config_dir = catalog.data["env"].get("ONIKA_CONFIG_PATH", "")
    merged = {}
    for name in includes:
        path = next((p for p in (os.path.join(base_dir or ".", str(name)), os.path.join(config_dir, str(name)))
                     if os.path.isfile(p)), None)
        if path is None:
            rep.add("error", "includes", f"included file '{name}' not found next to the input nor in {config_dir}")
            continue
        path = os.path.realpath(path)
        if path in seen:
            continue
        seen.add(path)
        try:
            with open(path) as f:
                sub = yaml.safe_load(f) or {}
        except (OSError, yaml.YAMLError) as e:
            rep.add("error", "includes", f"cannot read included file '{name}': {e}")
            continue
        if isinstance(sub, dict):
            merged = _merge(merged, _with_includes(rep, sub, os.path.dirname(path), catalog, seen))
    rest = {k: v for k, v in doc.items() if k != "includes"}
    return _merge(merged, rest)


def check_input(text, catalog, workdir=None, observed=None):
    """Static checks. ``observed``: operator -> parameters seen in the shipped
    examples (``Corpus.observed_parameters``), for the soft checks."""
    rep = _Report(observed)
    try:
        doc = yaml.safe_load(text)
    except yaml.YAMLError as e:
        # yaml-cpp (used by exaStamp) tolerates tabs where PyYAML does not
        try:
            doc = yaml.safe_load(text.expandtabs(2)) if "\t" in text else None
        except yaml.YAMLError:
            doc = None
        if doc is None:
            rep.add("error", "", f"YAML syntax error: {e}")
            return rep.issues
        rep.add("warning", "", "the file contains tab characters: use spaces for YAML indentation")
    if not isinstance(doc, dict):
        rep.add("error", "", "an input file is a YAML map of top-level blocks")
        return rep.issues
    own_keys = set(doc)
    own_species = _species_names(doc)   # an included file may declare every element (config_species.msp)
    doc = _with_includes(rep, doc, workdir, catalog)

    defaults = catalog.defaults
    known_names = set(defaults) | set(doc)   # names defined by the configuration or the input itself
    referenced = {n for _, n, _ in _batch_items(doc)}

    # top-level blocks
    for key, value in doc.items():
        if key in ("includes", "configuration", "simulation"):
            continue
        if key == "global":
            _check_global(rep, value, catalog)
        elif key == "species":
            if not isinstance(value, list) or not value:
                rep.add("error", "species", "species must be a non-empty list of 'Name: {mass, z, charge}'", block="species")
                continue
            for i, item in enumerate(value):
                if not (isinstance(item, dict) and len(item) == 1):
                    rep.add("error", f"species[{i}]", "each species is a one-key map 'Name: {mass: ..., z: ..., charge: ...}'", block="species")
                    continue
                name, props = next(iter(item.items()))
                if isinstance(props, dict) and "rigid_molecule" not in props:
                    _check_keys(rep, f"species[{i}].{name}", SCHEMAS["species"]["item_keys"], props, "species", closed=True)
        elif key == "setup_system":
            _check_batch(rep, "setup_system", value, catalog, known_names)
        elif key in ("grid_flavor", "numerical_scheme", "write_snapshot", "write_restart"):
            allowed = SCHEMAS[key]["values"]
            if not isinstance(value, str):
                continue   # redefined as a batch of operators
            if not (value in allowed or value in doc or value in defaults):
                rep.add("error", key, f"'{value}' is not one of {allowed}", suggestion=_closest(value, allowed), block=key)
        elif key == "compute_force":
            refs = value if isinstance(value, list) else [value]
            for r in refs:
                name, _ = _op_name_and_params(r)
                if name and name not in doc and catalog.get(name) is None and name not in defaults:
                    rep.add("error", f"compute_force.{name}", f"unknown force operator '{name}': misspelled, or its plugin is not part of this build",
                            suggestion=_closest(name, SCHEMAS["compute_force"]["operators"]), block="compute_force")
        elif key == "domain":
            _check_operator(rep, "domain", "domain", value, catalog)
            _check_domain(rep, "domain", value)
        elif catalog.get(key) is not None:
            # parameters of a registered operator (e.g. eam_alloy_force, langevin_thermostat)
            _check_operator(rep, key, key, value, catalog, definition=SCHEMAS.get(key, {}).get("kind") != "potential")
        elif key in defaults or key in SCHEMAS:
            if key in SCHEMAS and isinstance(value, dict):
                _check_operator(rep, key, key, value, catalog, definition=SCHEMAS.get(key, {}).get("kind") != "potential")
        elif key.startswith("+") and key[1:] in known_names:
            continue   # appends to an existing batch
        elif key not in own_keys or key in referenced:
            continue   # a batch defined by an included file, or used by another batch
        elif isinstance(value, (list, dict)):
            rep.add("info", key, f"'{key}' is not a known operator or configuration block: it defines a new batch, used only if referenced",
                    suggestion=_closest(key, list(catalog.names()) + list(defaults)))
        else:
            rep.add("warning", key, f"unknown top-level key '{key}'", suggestion=_closest(key, list(catalog.names()) + list(defaults)))

    # cross-block rules
    species = _species_names(doc)
    if not species and not any(n in ("read_dump_atoms", "read_extended_xyz_file") for n, _ in _setup_ops(doc)):
        rep.add("warning", "species", "no species block: species must come from particle_types or a restart file",
                suggestion=SCHEMAS["species"]["snippet"], block="species")
    ops = _setup_ops(doc)
    domain = None
    for i, (name, params) in enumerate(ops):
        if name == "domain":
            domain = params
            _check_domain(rep, f"setup_system[{i}].domain", params)
        if name == "lattice":
            _check_lattice(rep, f"setup_system[{i}].lattice", params, species, [n for n, _ in ops[:i]], domain or doc.get("domain"))
    _check_potential_files(rep, doc, own_species, catalog, workdir)
    if "compute_force" not in doc:
        rep.add("warning", "compute_force", "no interatomic potential: compute_force is 'nop' by default",
                suggestion=SCHEMAS["eam_alloy_force"]["snippet"], block="compute_force")
    _check_thermostat(rep, doc)
    _check_nose_hoover(rep, doc)
    _check_pressure(rep, doc)
    return rep.issues
