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
"""Small random linear ACE model (.yace, ye3t / PACE format) for the exaStamp ye3t regression tests.

ye3t_test_model.py out.yace [seed]

Two elements (Ta, W), ChebExpCos radial basis, rotation-invariant functions of rank 1 to 4 built from Wigner 3j
couplings with random weights (not fitted: for regression tests only). The model loads in LAMMPS pair_style ye3t and
pace, which can serve as references.
"""

import sys
from math import factorial, sqrt
import numpy as np


def w3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol (Racah formula)."""
    if m1 + m2 + m3 != 0 or j3 < abs(j1 - j2) or j3 > j1 + j2 or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    f = factorial
    tri = f(j1 + j2 - j3) * f(j1 - j2 + j3) * f(-j1 + j2 + j3) / f(j1 + j2 + j3 + 1)
    pre = sqrt(tri * f(j1 + m1) * f(j1 - m1) * f(j2 + m2) * f(j2 - m2) * f(j3 + m3) * f(j3 - m3))
    s = 0.0
    for k in range(0, j1 + j2 + j3 + 1):
        d = [k, j1 + j2 - j3 - k, j1 - m1 - k, j2 + m2 - k, j3 - j2 + m1 + k, j3 - j1 - m2 + k]
        if min(d) < 0:
            continue
        s += (-1) ** k / np.prod([f(x) for x in d])
    return (-1) ** (j1 - j2 - m3) * pre * s


def rows_rank2(l):
    return [((m, -m), (-1) ** m / sqrt(2 * l + 1)) for m in range(-l, l + 1)]


def rows_rank3(l1, l2, l3):
    out = []
    for m1 in range(-l1, l1 + 1):
        for m2 in range(-l2, l2 + 1):
            m3 = -m1 - m2
            c = w3j(l1, l2, l3, m1, m2, m3)
            if abs(m3) <= l3 and abs(c) > 1e-14:
                out.append(((m1, m2, m3), c))
    return out


def rows_rank4(l1, l2, l3, l4, L):
    out = []
    for m1 in range(-l1, l1 + 1):
        for m2 in range(-l2, l2 + 1):
            for m3 in range(-l3, l3 + 1):
                m4 = -m1 - m2 - m3
                if abs(m4) > l4:
                    continue
                M = m3 + m4
                c = (-1) ** M * w3j(l1, l2, L, m1, m2, -M) * w3j(l3, l4, L, m3, m4, M)
                if abs(c) > 1e-14:
                    out.append(((m1, m2, m3, m4), c))
    return out


def main():
    out = sys.argv[1]
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 2026)
    el = ["Ta", "W"]
    nrad, lmax, nbase = 2, 3, 3
    L = []
    L.append("elements:\n- Ta\n- W\nE0:\n- -5.0\n- -6.0\nembeddings:\n")
    for k in range(2):
        L.append(f"  {k}:\n    ndensity: 1\n    FS_parameters:\n    - 1.0\n    - 1.0\n    npoti: FinnisSinclair\n"
                 f"    rho_core_cutoff: 100000.0\n    drho_core_cutoff: 250.0\n")
    L.append("bonds:\n")
    for i in range(2):
        for j in range(2):
            rc = 4.6 + 0.1 * (i + j)
            L.append(f"  ? - {i}\n    - {j}\n  : nradmax: {nrad}\n    lmax: {lmax}\n    nradbasemax: {nbase}\n"
                     f"    radbasename: ChebExpCos\n    radparameters:\n    - {0.5 + 0.1 * (i + j)}\n")
            rc3 = []
            for n in range(nrad):
                rc3.append([])
                for l in range(lmax + 1):
                    c = rng.normal(0, 0.5, nbase)
                    c[n] += 1.0
                    rc3[-1].append("[" + ", ".join(f"{x:.12f}" for x in c) + "]")
            L.append("    radcoefficients: [" + ", ".join("[" + ", ".join(r) + "]" for r in rc3) + "]\n")
            L.append(f"    prehc: 0.0\n    lambdahc: 0.5\n    rcut: {rc}\n    dcut: 0.01\n    rcut_in: 0.0\n    dcut_in: 0.0\n"
                     f"    inner_cutoff_type: distance\n")
    L.append("functions:\n")

    def func(mu0, mus, ns, ls, rows, w):
        rank = len(mus)
        s = (f"  - mu0: {mu0}\n    rank: {rank}\n    ndensity: 1\n    num_ms_combs: {len(rows)}\n"
             "    mus:\n" + "".join(f"    - {x}\n" for x in mus) + "    ns:\n" + "".join(f"    - {x}\n" for x in ns)
             + "    ls:\n" + "".join(f"    - {x}\n" for x in ls) + "    ms_combs:\n"
             + "".join(f"    - {m}\n" for r in rows for m in r[0]) + "    ctildes:\n" + "".join(f"    - {w * r[1]:.14e}\n" for r in rows))
        return s

    for mu in range(2):
        L.append(f"  {mu}:\n")
        fl = []
        for nu in range(2):
            for k in range(1, nbase + 1):
                fl.append(func(mu, [nu], [k], [0], [((0,), 1.0)], rng.normal(0, 1.0)))
        for nu in range(2):
            for n in range(1, nrad + 1):
                for l in range(0, lmax + 1):
                    fl.append(func(mu, [nu, nu], [n, n], [l, l], rows_rank2(l), rng.normal(0, 0.3)))
        fl.append(func(mu, [0, 1], [1, 1], [1, 1], rows_rank2(1), rng.normal(0, 0.3)))
        for ls in [(1, 1, 2), (2, 2, 2), (1, 2, 3), (0, 1, 1)]:
            fl.append(func(mu, [mu, mu, mu], [1, 1, 2], list(ls), rows_rank3(*ls), rng.normal(0, 0.2)))
        fl.append(func(mu, [0, 1, 1], [1, 2, 1], [2, 2, 2], rows_rank3(2, 2, 2), rng.normal(0, 0.2)))
        for ls, Lc in [((1, 1, 1, 1), 0), ((1, 1, 1, 1), 2), ((2, 2, 1, 1), 2)]:
            fl.append(func(mu, [mu] * 4, [1, 1, 1, 2], list(ls), rows_rank4(*ls, Lc), rng.normal(0, 0.1)))
        L.append("".join(fl))
    L.append("deltaSplineBins: 0.001\n")
    open(out, "w").write("".join(L))


if __name__ == "__main__":
    main()
