"""
``SimulationSpec`` -> commented exaStamp input file.

Everything that must agree is derived here, so it is right by construction:
* domain: bounds = repeats x lattice parameter, cell_size dividing every box
  edge, grid_dims = bounds / cell_size;
* lattice types: one per atom of the unit cell, each a declared species;
* setup order: domain, init_rcb_grid, lattice, init_temperature_new;
* grid flavor: multimat when several species;
* thermostat block matching the numerical scheme;
* EAM cutoff, species masses and atomic numbers from the potential file or the
  periodic table.
"""

import math

import yaml

from onika_mcp.units import UnitError, parse_quantity

from .data import element, read_eam_alloy_header
from .knowledge import SCHEMAS, doc_url
from .spec import SimulationSpec

__all__ = ["render_input", "RenderError"]

_SCHEMES = {"nve": "verlet_nve", "langevin": "verlet_lnvt", "berendsen": "verlet_bnvt"}
_SNAPSHOTS = {"xyz": "write_snapshot_xyz", "paraview": "write_snapshot_paraview",
              "lammps_data": "write_snapshot_lammps_data"}


class RenderError(ValueError):
    """The description is inconsistent; the message says what to change."""


def _si(q, what):
    try:
        v, _, has_units = parse_quantity(q)
    except UnitError as e:
        raise RenderError(f"{what}: {e}")
    return v, has_units


def _ang(q, what):
    v, has_units = _si(q, what)
    return v * 1e10 if has_units else v


def _fmt(x):
    """12 significant digits: exaStamp compares grid_dims x cell_size to the bounds within 1e-15."""
    x = float(x)
    return f"{x:.1f}" if x.is_integer() else f"{x:.12g}"


def _q(x, unit):
    return f"{_fmt(x)} {unit}"


def _cell_size(box, target):
    """Largest common divisor c of the box edges with c close to ``target``:
    every edge must be an integer number of cells (exaStamp checks exactly)."""
    best = None
    for k in range(1, 100000):
        c = box[0] / k
        if c < target / 2.5:
            break
        dims = [b / c for b in box]
        if all(abs(d - round(d)) < 1e-9 and round(d) >= 1 for d in dims):
            score = abs(math.log(c / target))
            if best is None or score < best[0]:
                best = (score, c, [round(d) for d in dims])
    return best


def _rcut(spec, file_info):
    p = spec.potential
    if p.kind == "eam_alloy":
        if p.rcut:
            r = _ang(p.rcut, "potential.rcut")
            if file_info and r > file_info["cutoff_ang"] + 1e-9:
                raise RenderError(f"potential.rcut {r:g} ang is larger than the file cutoff {file_info['cutoff_ang']:g} ang")
            return r
        if file_info:
            return file_info["cutoff_ang"]
        raise RenderError("potential.rcut is needed when the potential file cannot be read")
    return _ang(p.rcut, "potential.rcut")


class _Dumper(yaml.SafeDumper):
    """Block style for maps, flow style for lists of scalars ([ 20, 20, 20 ])."""


def _represent_list(dumper, data):
    flow = all(not isinstance(v, (list, dict)) for v in data)
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=flow)


_Dumper.add_representer(list, _represent_list)


class _Doc:
    """Input file assembled block by block, each with its comment."""

    def __init__(self):
        self.parts = []

    def block(self, comment, data, doc=None):
        text = yaml.dump(data, Dumper=_Dumper, sort_keys=False, default_flow_style=False, width=100, allow_unicode=True)
        head = f"# --- {comment}" + (f" ({doc_url(doc)})" if doc else "") + "\n"
        self.parts.append(head + text)

    def raw(self, comment, text):
        self.parts.append(f"# --- {comment}\n{text.rstrip()}\n")

    def text(self, title=None):
        return (f"# {title}\n\n" if title else "") + "\n".join(self.parts)


def _doc(block):
    return SCHEMAS.get(block, {}).get("doc")


def render_input(spec, data_dirs, find_data_file):
    """Return (yaml_text, notes) for ``spec`` (a SimulationSpec or a dict).

    ``data_dirs``/``find_data_file`` come from the catalog (install data paths).
    ``notes`` explains the derived choices (cell size, flavor, masses...).
    """
    if not isinstance(spec, SimulationSpec):
        spec = SimulationSpec.model_validate(spec)
    notes = []
    d = _Doc()

    # --- potential file first: it gives elements, masses and the cutoff
    file_info = None
    if spec.potential.kind == "eam_alloy":
        path = find_data_file(spec.potential.file)
        if path is None:
            raise RenderError(f"potential file '{spec.potential.file}' not found in the data directories (use list_potential_files)")
        file_info = read_eam_alloy_header(path)
        for s in spec.species:
            if s.name not in file_info["elements"]:
                raise RenderError(f"species '{s.name}' is not in {spec.potential.file} (elements: {', '.join(file_info['elements'])}); name species like the file's elements")
    rcut = _rcut(spec, file_info)
    rcut_inc = _ang(spec.rcut_inc, "rcut_inc")

    # --- species
    species = []
    for s in spec.species:
        props = {}
        el = element(s.name, data_dirs)
        if s.mass:
            props["mass"] = s.mass
        elif el:
            props["mass"] = _q(el["mass"], "Da")
            notes.append(f"mass of {s.name}: {el['mass']} Da (periodic table)")
        else:
            raise RenderError(f"species '{s.name}': give its mass (not a chemical symbol)")
        z = s.z if s.z is not None else (el["z"] if el else None)
        if z is not None:
            props["z"] = z
        props["charge"] = s.charge
        species.append({s.name: props})
    names = [s.name for s in spec.species]
    d.block("Particle species", {"species": species}, _doc("species"))

    # --- system
    setup = []
    sysd = spec.system
    if sysd.kind == "lattice":
        a = sysd.lattice_parameter if isinstance(sysd.lattice_parameter, list) else [sysd.lattice_parameter] * 3
        if len(a) != 3:
            raise RenderError("system.lattice_parameter: give one value or three [a, b, c]")
        a = [_ang(x, "system.lattice_parameter") for x in a]
        if sysd.structure in ("SC", "BCC", "FCC", "c-DIA") and len(set(round(x, 9) for x in a)) > 1:
            notes.append(f"{sysd.structure} with unequal edges {a}: the cell is strained")
        if any(n < 1 for n in sysd.repeats):
            raise RenderError("system.repeats must be positive")
        box = [n * x for n, x in zip(sysd.repeats, a)]
        target = rcut + rcut_inc
        found = _cell_size(box, target)
        xform = None
        if found is not None:
            _, cell, dims = found
        else:
            # incommensurate edges (e.g. 2BCT, HCP): cubic cells in grid space,
            # stretched to the box by a diagonal xform
            cell = max(box) / max(1, round(max(box) / target))
            dims = [max(1, round(b / cell)) for b in box]
            cell = float(f"{cell:.10g}")
            xform = [[0.0] * 3 for _ in range(3)]
            for i in range(3):
                xform[i][i] = float(f"{box[i] / (dims[i] * cell):.15g}")
            notes.append(f"box edges {[round(b, 6) for b in box]} ang have no common cell size: grid cells are stretched by xform diagonal {[xform[i][i] for i in range(3)]}")
        # grid_dims x cell_size must equal the bounds exactly: round the cell
        # size, then write the bounds as that product (same float as exaStamp's)
        cell = float(f"{cell:.10g}")
        grid_box = [n * cell for n in dims]
        count = SCHEMAS["lattice"]["structure_types_count"][sysd.structure]
        types = sysd.types or [names[0]] * count
        if len(types) != count:
            raise RenderError(f"{sysd.structure} needs {count} types (one per atom of the unit cell), got {len(types)}")
        for t in types:
            if t not in names:
                raise RenderError(f"system.types: '{t}' is not a declared species ({', '.join(names)})")
        notes.append(f"box {' x '.join(f'{b:g}' for b in box)} ang, cell_size {cell:g} ang ({' x '.join(map(str, dims))} cells), rcut + rcut_inc = {target:g} ang")
        domain = {"cell_size": _q(cell, "ang"), "grid_dims": dims,
                  "bounds": [[_q(0, "ang")] * 3, [_q(b, "ang") for b in grid_box]]}
        if xform:
            domain["xform"] = xform
        domain.update(periodic=sysd.periodic, expandable=False)
        setup.append({"domain": domain})
        setup.append("init_rcb_grid")
        # a small shift keeps the atoms off the upper box faces: without it,
        # n * a rounding just below the box edge duplicates a layer of atoms
        setup.append({"lattice": {"structure": sysd.structure, "types": types, "size": [_q(x, "ang") for x in a],
                                  "shift": [_q(round(1e-3 * x, 6), "ang") for x in a]}})
        t0 = spec.initial_temperature or spec.ensemble.temperature or spec.ensemble.temperature_start
        if t0:
            setup.append({"init_temperature_new": {"T": t0}})
            if not spec.initial_temperature:
                notes.append(f"initial velocities at {t0} (ensemble temperature)")
        n_atoms = count * math.prod(sysd.repeats)
        notes.append(f"{n_atoms} atoms")
    else:
        cell = max(rcut + rcut_inc, 2.0)
        setup.append({"domain": {"cell_size": _q(round(cell, 3), "ang"), "periodic": sysd.periodic, "expandable": False}})
        setup.append({"read_extended_xyz_file": {"filename": sysd.filename, "bounds_mode": "FILE"}})
        if spec.initial_temperature:
            setup.append({"init_temperature_new": {"T": spec.initial_temperature}})
    d.block("Domain and initial system", {"setup_system": setup}, _doc("setup_system"))

    # --- grid flavor
    flavor = spec.grid_flavor or ("grid_flavor_multimat" if len(names) > 1 else None)
    if flavor:
        d.block("Per-particle fields", {"grid_flavor": flavor}, _doc("grid_flavor"))
        if not spec.grid_flavor:
            notes.append("grid_flavor_multimat because there are several species")

    # --- potential
    p = spec.potential
    if p.kind == "eam_alloy":
        force = {"eam_alloy_force": {"rcut": _q(rcut, "ang"), "parameters": {"file": p.file}}}
        name = "eam_alloy_force"
    else:
        if len(names) > 1:
            raise RenderError("the lj potential handles a single species; use extra with lj_multi_force for mixtures")
        force = {"lj_compute_force": {"rcut": _q(rcut, "ang"), "parameters": {"epsilon": p.epsilon, "sigma": p.sigma}}}
        name = "lj_compute_force"
    force["compute_force"] = [name]
    d.block("Interatomic potential", force, _doc(name))

    # --- integrator and thermostat
    e = spec.ensemble
    scheme = {"numerical_scheme": _SCHEMES[e.kind]}
    if e.kind in ("langevin", "berendsen"):
        if e.temperature and (e.temperature_start or e.temperature_stop):
            raise RenderError("ensemble: give temperature, or temperature_start and temperature_stop, not both")
        if e.temperature:
            target = {"T": e.temperature}
        elif e.temperature_start and e.temperature_stop:
            target = {"Tstart": e.temperature_start, "Tstop": e.temperature_stop}
        else:
            raise RenderError(f"ensemble {e.kind} needs temperature (or temperature_start and temperature_stop)")
        if e.kind == "langevin":
            scheme["langevin_thermostat"] = dict(target, gamma=e.gamma)
        else:
            scheme["berendsen_thermostat"] = dict(target, tau=e.tau)
    elif e.temperature or e.temperature_start:
        notes.append("ensemble nve: the temperature is only used for the initial velocities")
    d.block("Time integration", scheme, _doc("numerical_scheme"))

    # --- control
    o = spec.outputs
    g = {"max_iteration": spec.steps, "dt": spec.timestep, "rcut_inc": spec.rcut_inc,
         "simulation_thermostate_screen_frequency": o.thermo_every,
         "simulation_thermostate_file_frequency": o.thermo_every,
         "simulation_snapshot_frequency": o.snapshot_every if o.snapshot != "none" else -1,
         "simulation_restart_frequency": o.restart_every if o.restart_every else -1,
         "simulation_analysis_frequency": 0}
    d.block("Global control", {"global": g}, _doc("global"))

    outputs = {}
    if o.snapshot != "none":
        outputs["write_snapshot"] = _SNAPSHOTS[o.snapshot]
    if o.restart_every:
        outputs["write_restart"] = "write_restart_atoms"
    if outputs:
        d.block("Outputs", outputs, _doc("write_snapshot"))

    if spec.extra:
        try:
            yaml.safe_load(spec.extra)
        except yaml.YAMLError as err:
            raise RenderError(f"extra is not valid YAML: {err}")
        d.raw("Extra blocks (from the description's extra field)", spec.extra)

    return d.text(spec.title), notes
