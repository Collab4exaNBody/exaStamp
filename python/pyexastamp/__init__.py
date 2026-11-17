"""
pyexastamp -- Python interface to the exaStamp molecular dynamics code.

Importing this package sets up the onika runtime for exaStamp (onika, exaNBody
and exaStamp plugins, exaStamp's ``main-config.msp``, data and potential
directories), registers extractors for exaStamp types and re-exports the
pyonika / pyexanbody API.

Typical usage::

    import pyexastamp as xs

    # run an input file end to end, overriding some operator definitions
    xs.run_file("my_input.msp", overrides={"global": {"simulation_end_iteration": 100}})

    # or keep the simulation alive to read results before end()
    ctx = xs.init(xs.command_line(xs.find_example("single_specy.msp")))
    xs.run(ctx)
    thermo = xs.read_thermo(ctx)      # temperature (K), pressure (Pa), energies (eV/atom)...
    species = xs.read_species(ctx)    # [{"name": "Sn", "mass": 118.71, ...}]
    atoms = xs.particles(ctx)         # needs an extract_particle_fields node in the graph
    xs.end(ctx)
"""

import os

import pyexanbody
from pyexanbody import (  # noqa: F401  (re-exported API)
    init,
    run,
    run_node,
    end,
    make_operator,
    available_operators,
    get_operator_defaults,
    set_operator_defaults,
    build_simulation_graph,
    slot_as_array,
    setup_env,
    particles,
    read_cell_values,
    ApplicationContext,
    OperatorNode,
    OperatorSlotBase,
    OnikaError,
)
from pyonika import command_line as _command_line
from ._paths import ONIKA_INSTALL_DIR, XNB_INSTALL_DIR, APP_INSTALL_DIR as XSP_INSTALL_DIR

# exaStamp's data/config (main-config.msp), plugins and data directories
setup_env([ONIKA_INSTALL_DIR, XNB_INSTALL_DIR, XSP_INSTALL_DIR])

# Registers extractors for ThermodynamicState, ParticleSpecies, MoleculeSpeciesVector.
from . import _exastamp_data  # noqa: E402,F401

#: Root of the installed exaStamp data tree (config, samples, tutorial, regression_new, potentials).
XSP_DATA_DIR = os.path.join(XSP_INSTALL_DIR, "data")

#: Data subdirectories holding example input files.
EXAMPLE_CATEGORIES = ("samples", "tutorial", "regression_new")

__all__ = [
    "init", "run", "run_node", "end", "make_operator", "available_operators",
    "get_operator_defaults", "set_operator_defaults", "build_simulation_graph", "slot_as_array",
    "setup_env", "ApplicationContext", "OperatorNode", "OperatorSlotBase", "OnikaError",
    "XSP_INSTALL_DIR", "XSP_DATA_DIR", "EXAMPLE_CATEGORIES",
    "command_line", "run_file", "examples", "find_example",
    "read_thermo", "read_species", "read_molecules", "particles", "read_cell_values",
]


# ---------------------------------------------------------------------------
# Running input files
# ---------------------------------------------------------------------------

def command_line(msp_file, *args, program="pyexastamp", overrides=None, config=None):
    """Build the argv list for ``init()``.

    ``overrides`` are merged on top of the input file's operator definitions,
    e.g. ``{"global": {"dt": "1.0e-3 ps", "simulation_end_iteration": 100}}``;
    ``config`` sets ``configuration:`` options, e.g. ``{"omp_num_threads": 4}``.
    """
    return _command_line(msp_file, *args, program=program, overrides=overrides, config=config)


def run_file(msp_file, overrides=None, program="pyexastamp", **config):
    """Run an input file end to end (init, run, end) and return 0 on success.

    ``overrides``: operator definitions merged on top of the file.
    ``**config``: ``configuration:`` options (e.g. ``omp_num_threads=4``).
    """
    return pyexanbody.pyonika.run_file(os.path.abspath(msp_file), program=program,
                                       overrides=overrides, config=config)


# ---------------------------------------------------------------------------
# Installed examples
# ---------------------------------------------------------------------------

def examples(category=None, contains=None):
    """List installed example input files, as paths relative to ``XSP_DATA_DIR``.

    ``category`` restricts to one of ``EXAMPLE_CATEGORIES``; ``contains`` keeps
    paths containing that substring (e.g. ``"eam"``).
    """
    categories = [category] if category else EXAMPLE_CATEGORIES
    found = []
    for cat in categories:
        root = os.path.join(XSP_DATA_DIR, cat)
        for dirpath, _, filenames in os.walk(root):
            for f in filenames:
                if f.endswith(".msp"):
                    rel = os.path.relpath(os.path.join(dirpath, f), XSP_DATA_DIR)
                    if contains is None or contains in rel:
                        found.append(rel)
    return sorted(found)


def find_example(name):
    """Absolute path of an installed example, given a path relative to
    ``XSP_DATA_DIR`` or a unique file name (e.g. ``"single_specy.msp"``)."""
    path = os.path.join(XSP_DATA_DIR, name)
    if os.path.isfile(path):
        return path
    matches = [e for e in examples() if os.path.basename(e) == os.path.basename(name) and e.endswith(name)]
    if len(matches) == 1:
        return os.path.join(XSP_DATA_DIR, matches[0])
    if not matches:
        raise FileNotFoundError(f"no example named '{name}' under {XSP_DATA_DIR}")
    raise FileNotFoundError(f"'{name}' is ambiguous, candidates: {matches}")


# ---------------------------------------------------------------------------
# Reading results (call between run() and end())
# ---------------------------------------------------------------------------

def _graph(graph):
    return graph.graph if isinstance(graph, ApplicationContext) else graph


def _read_slot(graph, slot, node=None):
    """Value of slot ``slot`` in the last node, in graph order, that holds one.

    Several operator instances may own the slot (e.g. the thermodynamic state
    computed in first_iteration and in the compute loop); the last one in graph
    order is the most recently updated.  ``node`` selects an operator by name or
    pathname instead.
    """
    result = None

    def _collect(op):
        nonlocal result
        if node is not None and op.pathname() != node and op.name() != node:
            return
        value = op.slot_as_array(slot)
        if value is not None:
            result = value

    _graph(graph).apply_graph(_collect)
    return result


def read_thermo(graph, node=None):
    """Last thermodynamic state computed by ``simulation_thermodynamic_state``.

    Returns a dict with the values printed in the log: ``temperature`` (K),
    ``temperature_xyz``, ``pressure`` (Pa), ``stress`` (3x3, Pa),
    ``total_energy`` / ``kinetic_energy`` / ``potential_energy`` (eV/atom),
    ``volume``, ``density`` (g/cm^3), ``particle_count``; ``units`` names the
    units and ``internal`` holds raw internal-unit totals (``virial``,
    ``ke_tensor``, energies, mass...).  ``None`` if not computed yet.
    """
    return _read_slot(graph, "thermodynamic_state", node)


def read_species(graph, node=None):
    """Particle species: list of dicts with ``type``, ``name``, ``mass`` (Da),
    ``charge`` (e), ``z`` and rigid molecule data when relevant."""
    return _read_slot(graph, "species", node)


def read_molecules(graph, node=None):
    """Flexible molecule definitions: ``{"molecules": [...], "bridge_molecules": [...]}``,
    each with ``name``, ``atom_types``, ``bonds``, ``angles``, ``torsions``, ``impropers``."""
    return _read_slot(graph, "molecules", node)
