# pyexastamp

Python interface to exaStamp. It runs exaStamp input files (`.msp`) from Python and reads the results back as Python and numpy objects:
- thermodynamic state;
- particle species and molecules;
- per-particle fields, including vector and tensor fields.

It is built on `pyexanbody` and `pyonika` and re-exports their API.

## Build and install

Configure exaStamp with `-DEXASTAMP_BUILD_PYTHON=ON` against an exaNBody and onika built with their Python packages (`-DEXANB_BUILD_PYTHON=ON`, `-DONIKA_BUILD_PYTHON=ON`). All three must use the same Python interpreter.

The install tree then contains a pip-installable package:

```bash
pip install <onika-prefix>/python <exanbody-prefix>/python <exastamp-prefix>/python
```

Install paths are baked into the packages, so no `setup-env.sh` is needed. If you do source one, the `ONIKA_*` variables it sets take precedence. `builds/build_install_*.sh` do all of this into the `pyexastamp` virtual environment.

## Quick start

```python
import pyexastamp as xs

# run end to end, overriding operator definitions of the input file
xs.run_file("my_input.msp", overrides={"global": {"max_iteration": 100}}, omp_num_threads=4)

# keep the run alive to read results (read before end(): it frees the graph)
msp = xs.find_example("regression_new/potentials/eam/eam_vniitf/single_specy.msp")
ctx = xs.init(xs.command_line(msp, overrides={
    "simulation_epilog": [{"extract_particle_fields": {"fields": ["rx", "ry", "rz", "vx", "vy", "vz"]}}],
}))
xs.run(ctx)
thermo = xs.read_thermo(ctx)
atoms = xs.particles(ctx)
xs.end(ctx)
```

## API

| Function | Returns |
|---|---|
| `command_line(msp, overrides=None, config=None)` | argv for `init()`. `overrides` are merged on top of the file's operators (via `configuration.set`); `config` sets `configuration:` options. |
| `run_file(msp, overrides=None, **config)` | init, run, end; returns 0 on success |
| `examples(category=None, contains=None)` | installed `.msp` files under `data/{samples,tutorial,regression_new}` |
| `find_example(name)` | absolute path from a relative path or a unique file name |
| `read_thermo(ctx_or_graph)` | dict, described below |
| `read_species(ctx_or_graph)` | list of `{type, name, mass (Da), charge (e), z, ...}` |
| `read_molecules(ctx_or_graph)` | `{molecules, bridge_molecules}`, each with `atom_types`, `bonds`, `angles`, `torsions`, `impropers` |
| `particles(ctx_or_graph, node=None)` | `{field: ndarray}` from an `extract_particle_fields` operator |

**`read_thermo` result.** The values match what the thermodynamic log prints:
- `temperature` (K) and `temperature_xyz`;
- `pressure` (Pa) and `stress` (3×3, Pa);
- `total_energy`, `kinetic_energy`, `potential_energy` (eV/atom);
- `volume`, `density` (g/cm³) and `particle_count`.

`units` names the unit of each value, and `internal` holds raw internal-unit totals (including `virial` and `ke_tensor`).

**`particles` shapes.** Every field keeps its native dtype, and all fields share the same particle order:

| Field type | Shape |
|---|---|
| scalar | `(N,)` |
| `Vec3d` (`angmom`, `couple`) | `(N,3)` |
| `Mat3d` (`virial`) | `(N,3,3)` |
| `Quaternion` (`orient`) | `(N,4)`, ordered w, x, y, z |
| `cmol` | `(N,4)` uint64 |

`id` is always included. Which fields exist depends on the `grid_flavor`. An empty `fields` list copies every field.

The full pyonika API is also available: `init`, `run`, `run_node`, `end`, `make_operator`, `available_operators`, `get_operator_defaults`, `set_operator_defaults`, `build_simulation_graph`, `slot_as_array`, `OnikaError`.

## Notes

- **Errors:** errors in the input (missing file, malformed YAML, unknown operator or field, and so on) raise `pyonika.OnikaError` (a `RuntimeError`), and the Python session survives. This holds for single-process runs; with several MPI ranks onika still aborts.
- **Log:** `pyonika.capture_output(stream)` redirects the C++ log written during a block to a Python stream (Jupyter, loggers).
- **MPI:** run scripts with `mpirun -n N python script.py`. Thermodynamic values are global. Particle fields are local to each process, unless `extract_particle_fields` gets `gather: true`, which collects them on rank 0. An uncaught Python exception on one process aborts all of them. `tests/mpi/` holds the multi-process checks (ctest `pyexastamp_mpi*`).
- **Several runs:** you can run several simulations in one process, e.g. `examples/temperature_sweep.py`.
- **Working directory:** onika writes `plugins.db` and the simulation outputs there.

## Examples and tests

- `examples/run_eam_sample.py` runs a sample and prints its thermodynamics, species and particle arrays.
- `examples/temperature_sweep.py` runs several simulations in one process.
- `tests/` holds the pytest suite, also registered in ctest as `pyexastamp_tests`.

## MCP server for LLM agents

`exastamp_mcp` is an [MCP](https://modelcontextprotocol.io) server that helps an LLM agent (any model) write input files that run the first time. Install it with the `mcp` extra:

```bash
pip install "<prefix>/python[mcp]"
exastamp-mcp --workspace my_runs --doc-dir path/to/doc_exaStamp/docs   # stdio server
exastamp-mcp --index-examples    # once per install: records which shipped examples build
```

Register it in an MCP client, e.g. `.mcp.json` for Claude Code:

```json
{ "mcpServers": { "exastamp": {
    "command": "<venv>/bin/exastamp-mcp",
    "args": ["--workspace", "."],
    "env": { "EXASTAMP_DOC_DIR": "<path>/doc_exaStamp/docs" } } } }
```

What it does:

- **Grounding:** `list_capabilities`, `describe_block` (curated schemas of the blocks users write, plus the operator's real parameters), `search_examples` / `get_example`, `search_docs`, `list_potential_files`.
- **Structured description:** `build_input` turns a `SimulationSpec` (species, lattice or XYZ file, EAM alloy or Lennard-Jones, NVE / Langevin / Berendsen, outputs) into a commented input. Domain arithmetic, lattice types, grid flavor and block order are derived, so they are right by construction.
- **Layered validation:** `validate_input` runs static checks (unknown parameters with suggestions, units and dimensions, domain and lattice consistency, potential files, thermostat set-up), then a strict graph build, then a 2-step smoke run in a child process, checking the result for NaN, temperature explosions and the expected atom count. Issues are `{stage, level, path, message, suggestion, doc}`.
- **Workspace:** `write_input` (refused while validation fails), `run_input` / `job_status` / `stop_job` for local background runs.

The server never loads exaStamp itself. Catalog, build and runs go through child processes (`onika_mcp.worker`), so an abort cannot take it down. `mcp_eval/` holds 20 requests with expected results: `python mcp_eval/run_eval.py --reference` checks they can all be met, `--claude-code` measures Claude Code on them with your Claude account (no API billing), and `--api-model <model>` does the same through the Anthropic API (billed per token, needs `pip install "anthropic[mcp]"`). `run_eval.grade()` can score workspaces produced by another framework.
