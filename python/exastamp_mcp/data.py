"""
Data files shipped with exaStamp: periodic table and potential files.

Read straight from the install's data directories (found through the catalog),
so the server never imports pyexastamp.
"""

import csv
import os

__all__ = ["element", "read_eam_alloy_header", "potential_file_info", "potential_files",
           "POTENTIAL_FORMATS", "FILE_PARAMETERS"]

# potential formats, by file name ending -> operator reading them
POTENTIAL_FORMATS = {
    ".eam.alloy": "eam_alloy_force",
    ".yace": "pace_init",
    ".ace": "pace_init",
    ".snapcoeff": "snap_force_fp64",
    ".snapparam": "snap_force_fp64",
    ".pod": "pod_init",
    ".almtp": "mtp_init",
    ".mtp": "mtp_init",
    ".sw": "sw_init",
    ".dpnative.yaml": "deepmd_init",
    ".rebo": "rebo_pairwise",
}

# operator -> keys of its 'parameters' map that name a potential file
FILE_PARAMETERS = {
    "eam_alloy_force": ("file",),
    "snap_force_fp64": ("param", "coef"),
    "snap_force_fp32": ("param", "coef"),
    "pace_init": ("coef",),
    "pod_init": ("pod_file", "coeff_file"),
    "mtp_init": ("mtp_file",),
    "sw_init": ("file",),
    "deepmd_init": ("file",),
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


def _lines(path, limit=None):
    """Non-empty, non-comment lines of a text file."""
    out = []
    with open(path, errors="replace") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                out.append(line)
                if limit and len(out) >= limit:
                    break
    return out


def potential_file_info(path):
    """What can be read from a potential file without the potential code: 'elements'
    (species names it defines), 'species_count', 'cutoff_ang'. Empty if unknown format."""
    name = os.path.basename(path)
    if name.endswith(".eam.alloy"):
        h = read_eam_alloy_header(path)
        return {"elements": h["elements"], "cutoff_ang": h["cutoff_ang"]}
    if name.endswith(".snapcoeff"):
        lines = _lines(path)
        nelem, ncoeff = (int(x) for x in lines[0].split()[:2])
        return {"elements": [lines[1 + i * (ncoeff + 1)].split()[0] for i in range(nelem)]}
    if name.endswith(".snapparam"):
        values = dict(l.split()[:2] for l in _lines(path) if len(l.split()) >= 2)
        return {"rcutfac": float(values["rcutfac"])} if "rcutfac" in values else {}
    if name.endswith((".ace", ".yace")):
        for line in _lines(path, 20):
            if line.startswith("elements:"):
                return {"elements": line.split(":", 1)[1].strip(" []").replace(",", " ").split()}
        return {}
    if name.endswith(".pod"):
        for line in _lines(path, 20):
            if line.split()[0] == "species":
                return {"elements": line.split()[1:]}
        return {}
    if name.endswith((".almtp", ".mtp")):
        for line in _lines(path, 20):
            if line.replace(" ", "").startswith("species_count="):
                return {"species_count": int(line.split("=")[1])}
        return {}
    if name.endswith(".dpnative.yaml"):
        for line in _lines(path, 20):
            if line.startswith("type_map:"):
                return {"elements": line.split(":", 1)[1].strip(" []").replace(",", " ").split()}
        return {}
    if name.endswith(".sw"):
        elements = []
        for line in _lines(path):
            words = line.split()
            if len(words) >= 3 and all(w[0].isalpha() for w in words[:3]):
                elements += [w for w in dict.fromkeys(words[:3]) if w not in elements]
        return {"elements": elements}
    return {}


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
            try:
                entry.update(read_eam_alloy_header(os.path.join(d, f)) if f.endswith(".eam.alloy")
                             else potential_file_info(os.path.join(d, f)))
            except (OSError, StopIteration, ValueError, IndexError, KeyError):
                entry["error"] = "unreadable header"
            if element_name and "elements" in entry and element_name not in entry["elements"]:
                continue
            if element_name and "elements" not in entry and element_name not in f:
                continue
            found.append(entry)
    return found
