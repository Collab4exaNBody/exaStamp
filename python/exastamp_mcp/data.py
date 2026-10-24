"""
Data files shipped with exaStamp: periodic table and potential files.

Read straight from the install's data directories (found through the catalog),
so the server never imports pyexastamp.
"""

import csv
import os

__all__ = ["element", "read_eam_alloy_header", "potential_files"]

# potential formats, by file name ending -> force operator using them
POTENTIAL_FORMATS = {
    ".eam.alloy": "eam_alloy_force",
    ".yace": "pace_force",
    ".ace": "pace_force",
    ".snapcoeff": "snap_force",
    ".snapparam": "snap_force",
    ".pod": "pod_force",
    ".rebo": "rebo_pairwise",
}

_PERIODIC_TABLE = None


def _periodic_table(data_dirs):
    global _PERIODIC_TABLE
    if _PERIODIC_TABLE is None:
        _PERIODIC_TABLE = {}
        for d in data_dirs:
            p = os.path.join(d, "periodic_table.csv")
            if os.path.isfile(p):
                with open(p, newline="") as f:
                    for row in csv.DictReader(f):
                        _PERIODIC_TABLE[row["Symbol"].strip()] = {
                            "z": int(row["AtomicNumber"]), "mass": float(row["AtomicMass"]),
                            "name": row["Element"].strip()}
                break
    return _PERIODIC_TABLE


def element(symbol, data_dirs):
    """{'z', 'mass' (Da), 'name'} of a chemical element, or None."""
    return _periodic_table(data_dirs).get(symbol)


def read_eam_alloy_header(path):
    """Elements, masses, atomic numbers and cutoff of a LAMMPS eam/alloy (setfl) file."""
    with open(path, errors="replace") as f:
        lines = [next(f) for _ in range(6)]
    elements = lines[3].split()[1:]
    cutoff = float(lines[4].split()[4])
    z, mass = (float(x) for x in lines[5].split()[:2])
    info = {"elements": elements, "cutoff_ang": cutoff, "first_element": {"z": int(z), "mass": mass}}
    if len(lines[5].split()) >= 4:
        info["first_element"]["lattice"] = lines[5].split()[3]
        info["first_element"]["lattice_parameter_ang"] = float(lines[5].split()[2])
    return info


def potential_files(data_dirs, operator=None, element_name=None):
    """Potential files found in the data directories, with what can be read from them."""
    found, seen = [], set()
    for d in data_dirs:
        if not os.path.isdir(d) or os.path.abspath(d) in seen:
            continue
        seen.add(os.path.abspath(d))
        for f in sorted(os.listdir(d)):
            op = next((o for ext, o in POTENTIAL_FORMATS.items() if f.endswith(ext)), None)
            if op is None or (operator and op != operator):
                continue
            entry = {"file": f, "operator": op, "directory": d}
            if f.endswith(".eam.alloy"):
                try:
                    entry.update(read_eam_alloy_header(os.path.join(d, f)))
                except (OSError, StopIteration, ValueError, IndexError):
                    entry["error"] = "unreadable header"
            if element_name and "elements" in entry and element_name not in entry["elements"]:
                continue
            if element_name and "elements" not in entry and element_name not in f:
                continue
            found.append(entry)
    return found
