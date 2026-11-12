# polymlp_export

Converts a polynomial machine-learning potential of A. Seko (`polymlp.lammps` or `polymlp.yaml`, from pypolymlp or the
Polynomial MLP Repository) into the table file read by exaStamp's `polymlp_init`. The invariant features are
enumerated by the C++ core of [pypolymlp](https://github.com/sekocha/pypolymlp) (BSD-3-Clause), which this tool links;
exaStamp itself does not depend on it.

```bash
pip install pypolymlp
make -j4 POLYMLP_SRC=$(python -c "import pypolymlp, os; print(os.path.dirname(pypolymlp.__file__))")/cxx/src
./polymlp_export polymlp.lammps polymlp.exastamp.yaml
```

On macOS, pip may need an OpenMP-capable compiler to build pypolymlp (e.g. `CC=$(brew --prefix llvm)/bin/clang CXX=$(brew --prefix llvm)/bin/clang++`);
the exporter itself needs none. Tested with pypolymlp 0.22.1. Gaussian radial functions only; one file per potential
(hybrid potentials: not supported yet).
