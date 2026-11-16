#!/usr/bin/env python3
#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements. See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership. The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
"""Small random tagged-Cauchy correction (ye3t tagged slice v2) and composite manifest for the exaStamp ye3t tests.

ye3t_tagged_test_model.py <backbone.yace> <out prefix> [seed]

Writes <prefix>.tagged.json (schema ye3t_tagged_cauchy_slice_v2) and <prefix>.ye3t.json (schema
ye3t_tagged_cauchy_composite_v1: the .yace backbone plus the tagged correction), for the elements and the cutoff of
the backbone (largest rcut of its bonds, dcut and lambda of its first bond). Channels per neighbour element: l = 0 (radials 0, 1), l = 1
(radial 0), l = 2 (radial 1), real orthonormal tesseral harmonics. Terms are rotation invariant by construction:
scalar densities and moments, dot products sum_a X_a Y_a of real components of the same l (densities, and moments of
one edge), weighted by the falling factorial of the edge count (tag_count 2, p = 0, 1, 2). Random weights, per
element: not a fitted model, for regression tests only. Only the fields read by exaStamp ye3t_init are written (no
compiler artefacts or hashes: the files are not meant for ye3t-lammps).
"""

import json
import os
import sys
import re
import numpy as np


def real_form(l):
    """Real orthonormal tesseral harmonics of degree l from the complex ones (Y_{l,-m} = (-1)^m conj Y_lm).

    Real coordinates, in order: cos(m) for m = l .. 1, m = 0, sin(m) for m = 1 .. l. Returns (magnetic order,
    real coordinate order, real_to_complex matrix U with complex[row] = sum_a U[row][a] real[a])."""
    mo = list(range(-l, l + 1))
    coords = [("cosine", m) for m in range(l, 0, -1)] + [("m_zero", 0)] + [("sine", m) for m in range(1, l + 1)]
    R = np.zeros((2 * l + 1, 2 * l + 1), dtype=complex)   # real = R complex
    for a, (kind, m) in enumerate(coords):
        if kind == "m_zero":
            R[a, mo.index(0)] = 1.0
        elif kind == "cosine":   # sqrt(2) (-1)^m Re Y_m = ((-1)^m Y_m + Y_-m) / sqrt(2)
            R[a, mo.index(m)] = (-1) ** m / np.sqrt(2)
            R[a, mo.index(-m)] = 1 / np.sqrt(2)
        else:                    # sqrt(2) (-1)^m Im Y_m = i (Y_-m - (-1)^m Y_m) / sqrt(2)
            R[a, mo.index(-m)] = 1j / np.sqrt(2)
            R[a, mo.index(m)] = -1j * (-1) ** m / np.sqrt(2)
    U = np.linalg.inv(R)
    return mo, [{"kind": k, "m": m} for k, m in coords], [[{"binary64": [float(U[r, a].real), float(U[r, a].imag)]} for a in range(2 * l + 1)] for r in range(2 * l + 1)]


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    backbone, prefix = sys.argv[1], sys.argv[2]
    rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 2026)
    text = open(backbone).read()   # the .yace keys of the bonds are lists (no plain YAML loader): fields by pattern
    m = re.search(r"^elements:\s*\[([^\]]*)\]", text, re.M)   # flow or block list
    elements = re.findall(r"[A-Z][a-z]?", m.group(1)) if m else re.findall(r"^- ([A-Z][a-z]?)\s*$", re.search(r"^elements:\s*\n((?:- .*\n)+)", text, re.M).group(1), re.M)
    rc = max(float(x) for x in re.findall(r"^\s+rcut:\s*(\S+)", text, re.M))   # largest cutoff of the bonds
    dcut = float(re.search(r"^\s+dcut:\s*(\S+)", text, re.M).group(1))
    lam = float(re.search(r"radparameters:\s*\n\s*-\s*(\S+)", text).group(1))

    forms, records = {}, {}
    for l in (0, 1, 2):
        mo, coords, U = real_form(l)
        fid = f"real_form_test_l{l}"
        forms[l] = fid
        records[fid] = {"angular_l": l, "magnetic_order": mo, "real_coordinate_order": coords, "real_to_complex_matrix": U}

    channels, chan = [], {}   # chan[(nu, l, n)] = index
    for nu in elements:
        for l, n in ((0, 0), (0, 1), (1, 0), (2, 1)):
            chan[(nu, l, n)] = len(channels)
            channels.append({"channel_index": len(channels), "l": l, "neighbor_species": nu, "radial_channel": n, "real_form_id": forms[l]})

    dens, dkey = [], {}
    def d(c, a):
        if (c, a) not in dkey:
            dkey[(c, a)] = len(dens); dens.append([c, a])
        return dkey[(c, a)]
    moms, mkey = [], {}
    def m(*factors):
        k = tuple(sorted(factors))
        if k not in mkey:
            mkey[k] = len(moms); moms.append([list(f) for f in k])
        return mkey[k]

    terms, nf = [], 0
    def feature(rows, p):   # rows: list of (coefficient, density indices, moment indices), one invariant feature
        nonlocal nf
        for c, di, mi in rows:
            terms.append({"feature_index": nf, "coefficient": c, "density_factor_indices": di, "moment_indices": mi, "p": p})
        nf += 1

    e0, e1 = elements[0], elements[-1]
    feature([(1.0, [], [m()])], 2)                                                         # z (edge count)
    feature([(1.0, [d(chan[(e0, 0, 0)], 0), d(chan[(e1, 0, 0)], 0)], [])], 0)              # A_00 A_00'
    feature([(1.0, [d(chan[(e0, 0, 1)], 0)], [])], 1)                                      # A_01
    for l, n in ((1, 0), (2, 1)):
        for x, w in ((e0, e1), (e1, e1)):
            feature([(1.0, [d(chan[(x, l, n)], a), d(chan[(w, l, n)], a)], []) for a in range(2 * l + 1)], 1)   # A . A'
            feature([(1.0, [], [m((chan[(x, l, n)], a), (chan[(w, l, n)], a))]) for a in range(2 * l + 1)], 2)  # sum_j phi . phi'
    feature([(1.0, [], [m((chan[(e0, 0, 1)], 0)), m((chan[(e1, 0, 0)], 0))])], 1)          # M M'
    feature([(1.0, [d(chan[(e0, 0, 0)], 0)], [m((chan[(e0, 2, 1)], a), (chan[(e0, 2, 1)], a))]) for a in range(5)], 0)
    feature([(1.0, [], [m((chan[(e0, 1, 0)], a), (chan[(e0, 1, 0)], a), (chan[(e1, 0, 1)], 0))]) for a in range(3)], 2)

    # weights: features of order ~ 1 eV in a bcc environment (z ~ 14, F(z, 0) ~ 200)
    scale = np.array([2e-4, 2e-2, 2e-2, 1e-1, 2e-3, 1e-1, 2e-3, 1e-1, 2e-3, 1e-1, 2e-3, 2e-2, 2e-2, 2e-3])[:nf]
    beta = {e: [float(x) for x in scale * rng.uniform(-1.0, 1.0, nf)] for e in elements}
    tagged = {
        "schema": "ye3t_tagged_cauchy_slice_v2",
        "synthetic_fixture_note": "random test model of exaStamp scripts/python/ye3t_tagged_test_model.py, not fitted",
        "species_order": elements, "cutoff": rc, "tag_count": 2,
        "pooling": "ordered_sum", "context_policy": "inclusive", "combination_matrix": None, "descriptor_selection": None,
        "offsets": {e: 0.0 for e in elements}, "beta": beta,
        "radial_definition": {"kind": "pace_cheb_exp_cos", "parameters": {"rc": rc, "cutoff_width": dcut, "lmbda": lam, "radial_count": 2}},
        "real_form_records": records, "channel_real_forms": channels,
        "real_moment_program": {"feature_count": nf, "tag_count": 2, "real_density_keys": dens, "real_moment_keys": moms, "terms": terms},
    }
    tpath = prefix + ".tagged.json"
    json.dump(tagged, open(tpath, "w"), indent=1)
    comp = {
        "schema": "ye3t_tagged_cauchy_composite_v1", "species_order": elements,
        "ordinary_component": {"path": os.path.relpath(backbone, os.path.dirname(os.path.abspath(prefix)))},
        "tagged_component": {"path": os.path.basename(tpath)},
    }
    json.dump(comp, open(prefix + ".ye3t.json", "w"), indent=1)
    print(f"{prefix}.ye3t.json: {len(elements)} elements, rc {rc}, {len(channels)} channels, {len(dens)} densities, {len(moms)} moments, "
          f"{nf} features, {len(terms)} terms")


if __name__ == "__main__":
    main()
