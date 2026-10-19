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
#!/usr/bin/python3
"""Converts a flexible molecule system in Stamp FAtomes format (+ its DONNEES init file) into exaStamp input files :

 - a YAML force field file : species, and a molecule_force_field block (intramolecular potentials, intramolecular
   pair weights, pair potentials, reaction field and long range correction flags), to include in an input file
   that reads the system with read_xyz_molecules followed by molecule_force_field ;
 - optionally, the system itself in read_xyz_molecules format (cell, positions, bonds, charges).

The conversions are those of read_fatomes_mol, so that both ways of reading the system give the same forces and
energies. Atom velocities are not converted (read_xyz_molecules starts from zero velocities).

usage:
  fatomes_to_yaml.py FAtomes [--donnees DONNEES] [--xyz out.xyz] [--molecule NAME] [--strip-prefix PREFIX] [-o out.cfg]

  --donnees       DONNEES init file : reaction field (ReactionField, ReactionField_epsilon, ReactionField_rc,
                  ReactionField_selfpairs) and long range correction (CorrectionVdwLongueDistance)
  --xyz           also write the system in read_xyz_molecules format
  --molecule      molecule type name (default MOL, as read_fatomes_mol)
  --strip-prefix  remove this prefix from FAtomes atom names to make species names (e.g. TATB_ : TATB_O -> O),
                  force field type names (nomFF) are unchanged

example:
  fatomes_to_yaml.py FAtomes_tatb_gas --donnees DONNEES_tatb_gas --molecule TATB --strip-prefix TATB_ \\
                     --xyz tatb.xyz -o tatb_ff.cfg
"""
import argparse
import math
import sys

# unit strings rewritten as read_fatomes_mol does before parsing them
UNIT_REPLACE = [("J.m6", "J*m^6"), ("ang2", "ang^2"), ("ang-1", "ang^-1"), ("deg", "degree"), ("degreere", "degree"),
                ("m6", "m^6"), ("m-1", "m^-1"), ("metre", "m"), ("kcal/mol.ang^6", "kcal/mol*ang^6"),
                ("kcal/mol.ang6", "kcal/mol*ang^6")]
LENGTH_TO_ANG = {"ang": 1.0, "A": 1.0, "nm": 10.0, "m": 1.0e10}


def read_unit(u):
    for a, b in UNIT_REPLACE:
        p = u.find(a)
        while p != -1:
            u = u[:p] + b + u[p + len(a):]
            p = u.find(a, p + len(b))
    return u.strip()


def q(v, u):
    """YAML quantity"""
    return f"{v!r} {u}" if isinstance(v, float) else f"{v} {u}"


def lines_of(path):
    """non empty lines, '*' starts a comment"""
    out = []
    for l in open(path, errors="replace"):
        l = l.split("*")[0].rstrip()
        if l.strip():
            out.append(l)
    return out


def fail(msg):
    sys.exit("fatomes_to_yaml: " + msg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fatomes")
    ap.add_argument("--donnees")
    ap.add_argument("--xyz")
    ap.add_argument("--molecule", default="MOL")
    ap.add_argument("--strip-prefix", default="")
    ap.add_argument("-o", "--output", default="-")
    args = ap.parse_args()

    L = lines_of(args.fatomes)
    words = [l.split() for l in L]

    def section(key):
        for i, w in enumerate(words):
            if w[0] == key:
                return i
        return None

    # ---------------- species ----------------
    species = []  # FAtomes name, force field name, mass (string quantity), charge (string quantity)
    i = 0
    while i < len(words):
        if words[i][0] == "nom" and len(words[i]) > 1:
            d = {}
            j = i + 1
            while j < len(words) and words[j][0] in ("nomXYZ", "nomFF", "type", "masse", "charge"):
                d[words[j][0]] = words[j][1:]
                j += 1
            if d.get("type", ["Atome"])[0] != "Atome":
                fail(f"atom type '{' '.join(d['type'])}' of {words[i][1]} is not supported")
            for k in ("nomFF", "masse", "charge"):
                if k not in d:
                    fail(f"missing '{k}' for atom {words[i][1]}")
            species.append((words[i][1], d["nomFF"][0], q(float(d["masse"][0]), read_unit(" ".join(d["masse"][1:]))),
                            q(float(d["charge"][0]), read_unit(" ".join(d["charge"][1:])))))
            i = j
        else:
            i += 1
    if not species:
        fail("no atom type found")

    def sp_name(n):
        if args.strip_prefix and not n.startswith(args.strip_prefix):
            fail(f"atom name {n} does not start with {args.strip_prefix}")
        return n[len(args.strip_prefix):]

    names = [sp_name(s[0]) for s in species]
    if len(set(names)) != len(names):
        fail("species names are not unique : " + " ".join(names))
    fat_index = {s[0]: k for k, s in enumerate(species)}

    def pair_type(t):
        """pair potential type, given as an atom type index or an atom name"""
        return names[int(t)] if t[0].isdigit() else sp_name(t)

    # ---------------- pair potentials ----------------
    pairs = []  # [type_a, type_b, {"rcut": (value, unit) or None, "exp6"/"lj": dict of quantities, "rf": ...}]
    rcut_max_ang = 0.0
    for w in words:
        if w[0] != "Potentiel":
            continue
        ta, tb, kind = pair_type(w[1]), pair_type(w[2]), w[3]
        p = {}
        if kind == "LJ":
            v = w[4:13]
            prm = {}
            prop_sigma = False
            for k in range(3):
                name, val, unit = v[3 * k], float(v[3 * k + 1]), v[3 * k + 2]
                if unit == "-":
                    prop_sigma, unit = True, "ang"
                prm[name] = (val, read_unit(unit))
            sigma, eps, rc = prm["sigma"], prm["epsilon"], prm["rc"]
            sigma_ang = sigma[0] * LENGTH_TO_ANG.get(sigma[1], float("nan"))
            if prop_sigma:  # rc given in units of sigma
                if math.isnan(sigma_ang):
                    fail(f"unsupported sigma unit {sigma[1]}")
                rc = (rc[0] * sigma_ang, "ang")
            ext = w[13:]
            if ext and ext[0] == "CONV_EXP6v1":  # Lennard-Jones parameters converted to Exp6
                alpha, D = float(ext[2]), (float(ext[4]), read_unit(ext[5]))
                if math.isnan(sigma_ang):
                    fail(f"unsupported sigma unit {sigma[1]}")
                rmin = sigma_ang * 2.0 ** (1.0 / 6.0)
                A = 6.0 * eps[0] * math.exp(alpha) / (alpha - 6.0)
                B = alpha / rmin
                C = eps[0] * alpha / (alpha - 6.0) * rmin ** 6
                p["exp6"] = {"A": (A, eps[1]), "B": (B, "1/ang"), "C": (C, f"{eps[1]}*ang^6"), "D": D}
            elif ext:
                fail(f"unexpected LJ modifier '{ext[0]}'")
            else:
                p["lj"] = {"epsilon": eps, "sigma": sigma}
            p["rcut"] = rc
        elif kind == "Exp6v1":
            v = w[4:19]
            prm = {v[3 * k]: (float(v[3 * k + 1]), read_unit(v[3 * k + 2])) for k in range(5)}
            p["exp6"] = {"A": prm["a"], "B": prm["b"], "C": prm["c"], "D": prm["d"]}
            p["rcut"] = prm["rc"]
        else:
            fail(f"pair potential {kind} is not supported")
        rc_ang = p["rcut"][0] * LENGTH_TO_ANG.get(p["rcut"][1], float("nan"))
        rcut_max_ang = max(rcut_max_ang, rc_ang)
        pairs.append([ta, tb, p])

    # ---------------- cross pairs from a mixing rule (Regle_melange), as read_fatomes_mol ----------------
    rule = None
    for w in words:
        if w[0] == "Regle_melange" and len(w) > 1:
            rule = w[1]
    defined = {(a, b) for a, b, _ in pairs} | {(b, a) for a, b, _ in pairs}
    self_pair = {a: p for a, b, p in pairs if a == b}

    def same_unit(x, y, what):
        if x[1] != y[1]:
            fail(f"mixing rule : {what} given in different units ({x[1]}, {y[1]})")
        return x[1]

    for ia in range(len(names) - 1):
        for ib in range(ia + 1, len(names)):
            a, b = names[ia], names[ib]
            if (a, b) in defined:
                continue
            if a not in self_pair or b not in self_pair:
                fail(f"cross pair potential {a} / {b} is not defined and cannot be built by a mixing rule")
            pa, pb = self_pair[a], self_pair[b]
            ru = same_unit(pa["rcut"], pb["rcut"], "rcut")
            ra, rb = pa["rcut"][0], pb["rcut"][0]
            if "lj" in pa and "lj" in pb:
                r = rule or "Lorentz-Berthelot"
                ea, eb, sa, sb = pa["lj"]["epsilon"], pb["lj"]["epsilon"], pa["lj"]["sigma"], pb["lj"]["sigma"]
                eu, su = same_unit(ea, eb, "epsilon"), same_unit(sa, sb, "sigma")
                if r == "geometrique":
                    e, sg, rc = math.sqrt(ea[0] * eb[0]), math.sqrt(sa[0] * sb[0]), math.sqrt(ra * rb)
                elif r in ("Lorentz-Berthelot", "LB"):
                    e, sg, rc = math.sqrt(ea[0] * eb[0]), 0.5 * (sa[0] + sb[0]), 0.5 * (ra + rb)
                else:
                    fail(f"unexpected mixing rule {r} for LJ cross pair {a} / {b}")
                pairs.append([a, b, {"lj": {"epsilon": (e, eu), "sigma": (sg, su)}, "rcut": (rc, ru)}])
            elif "exp6" in pa and "exp6" in pb:
                r = rule or "Waldman-Hagler"
                xa, xb = pa["exp6"], pb["exp6"]
                u = {k: same_unit(xa[k], xb[k], k) for k in "ABCD"}
                Aa, Ab, Ba, Bb = xa["A"][0], xb["A"][0], xa["B"][0], xb["B"][0]
                C, D = math.sqrt(xa["C"][0] * xb["C"][0]), math.sqrt(xa["D"][0] * xb["D"][0])
                rc = math.sqrt(ra * rb)
                if r == "geometrique":
                    A, B = math.sqrt(Aa * Ab), math.sqrt(Ba * Bb)
                elif r == "Waldman-Hagler":
                    B = (2.0 / (Ba ** -6 + Bb ** -6)) ** (1.0 / 6.0)
                    A = math.sqrt(Aa * Ab) * B ** 6 / Ba ** 3 / Bb ** 3
                else:
                    fail(f"unexpected mixing rule {r} for Exp-6 cross pair {a} / {b}")
                pairs.append([a, b, {"exp6": {"A": (A, u["A"]), "B": (B, u["B"]), "C": (C, u["C"]), "D": (D, u["D"])}, "rcut": (rc, ru)}])
            else:
                fail(f"atom types {a} and {b} do not have the same pair potential, the cross term cannot be defined")

    # ---------------- reaction field and long range correction (DONNEES) ----------------
    rf_self_pairs, long_range_correction = False, False
    if args.donnees:
        init = {}
        for w in (l.split() for l in lines_of(args.donnees)):
            if len(w) >= 2:
                try:
                    init[w[0]] = float(w[1])
                except ValueError:
                    pass
        long_range_correction = init.get("CorrectionVdwLongueDistance", 0) == 1
        if init.get("ReactionField", 0) == 1:
            rf_self_pairs = init.get("ReactionField_selfpairs", 0) == 1
            rf_eps = init.get("ReactionField_epsilon", 0.0)
            rf_rc = init.get("ReactionField_rc", 0.0) * 1.0e10  # implicit meter unit
            if rf_rc < 0.0:
                rf_rc = rcut_max_ang
            existing = {(a, b) for a, b, _ in pairs} | {(b, a) for a, b, _ in pairs}
            for pr in pairs:
                pr[2]["rf"] = (rf_eps, rf_rc)
            for a in range(len(names)):
                for b in range(a, len(names)):
                    if (names[a], names[b]) not in existing:
                        pairs.append([names[a], names[b], {"rf": (rf_eps, rf_rc)}])

    # ---------------- intramolecular potentials ----------------
    bonds, bends, torsions, impropers = [], [], [], []
    k0 = section("ChampDeForces")
    stop = section("ContribDispRepIntra")
    ff_lines = words[k0 + 1:stop] if k0 is not None else []
    il = 0
    while il < len(ff_lines):
        w = ff_lines[il]
        il += 1
        kind = w[0]
        if kind == "torsion_gaff":
            # up to 3 consecutive lines of the same torsion, sum of k/m (1 + cos(n phi - phi0)), folded into the OPLS form
            # plus a constant energy offset, as read_fatomes_mol does (k1, k2, k3 with the sign given by n and phi0)
            terms = [w]
            while len(terms) < 3 and il < len(ff_lines) and ff_lines[il][0] == "torsion_gaff" and ff_lines[il][1:5] == w[1:5]:
                terms.append(ff_lines[il])
                il += 1
            kopls, offset, eunit = [0.0, 0.0, 0.0], 0.0, read_unit(w[6])
            for tw in terms:
                if read_unit(tw[6]) != eunit:
                    fail(f"torsion_gaff {' '.join(w[1:5])} : mixed energy units")
                k, m, n = float(tw[5]), float(tw[7]), int(abs(float(tw[9])))
                au = read_unit(tw[12])
                phi0 = float(tw[11]) * (math.pi / 180.0 if au == "degree" else 1.0)
                if au not in ("degree", "rad"):
                    fail(f"unsupported angle unit {au}")
                if n == 4:
                    fail("n = 4 for force field type torsion_gaff is not available")
                idx = n - 1
                if idx > -1:
                    kopls[idx] = 2.0 * k / m
                    if n in (1, 3) and phi0 > 0.0:
                        kopls[idx] *= -1.0
                    elif n == 2 and phi0 == 0.0:
                        kopls[idx] *= -1.0
                    offset += kopls[idx] * (-1) ** idx * math.cos(phi0 / 2.0 + (idx + 1) * math.pi / 2.0)
            torsions.append(f"{{ types: [{', '.join(w[1:5])}], potential: opls_torsion, parameters: {{ k1: {q(kopls[0], eunit)}, k2: {q(kopls[1], eunit)}, k3: {q(kopls[2], eunit)}, e0: {q(offset, eunit)} }} }}")
        elif kind in ("bond_opls", "bond_gaff", "bond_0.5_harm", "bond_harm", "bond_charmm"):
            t, r0, k = w[1:3], (w[3], read_unit(w[4])), (float(w[5]), read_unit(w[6]))
            if kind in ("bond_0.5_harm", "bond_charmm"):
                k = (k[0] / 2.0, k[1])
            bonds.append(f"{{ types: [{', '.join(t)}], potential: opls_bond, parameters: {{ r0: {q(*r0)}, k: {q(*k)} }} }}")
        elif kind in ("angle_opls", "angle_harm", "angle_gaff", "angle_0.5_harm", "angle_charmm"):
            t, t0, k = w[1:4], (w[4], read_unit(w[5])), (float(w[6]), read_unit(w[7]))
            if kind in ("angle_0.5_harm", "angle_charmm"):
                k = (k[0] / 2.0, k[1])
            bends.append(f"{{ types: [{', '.join(t)}], potential: opls_bend, parameters: {{ theta0: {q(*t0)}, k: {q(*k)} }} }}")
        elif kind in ("torsion_opls", "torsion_0.5_3cos-", "torsion_compass", "torsion_3cos-", "torsion_0.5_ncos"):
            t = w[1:5]
            if kind == "torsion_0.5_ncos":
                n = (int(w[5]), int(w[9]), int(w[13]))
                p = [(w[7], read_unit(w[8])), (w[11], read_unit(w[12])), (w[15], read_unit(w[16]))]
                if n == (1, 2, 3):
                    pot = "0.5compass_torsion"
                elif n[0] == 3:
                    pot, p = "0.5compass_torsion", [("0.0", p[0][1]), ("0.0", p[0][1]), p[0]]
                else:
                    fail(f"torsion_0.5_ncos with n = {n} has no input file form")
            else:
                p = [(w[5], read_unit(w[6])), (w[7], read_unit(w[8])), (w[9], read_unit(w[10]))]
                pot = {"torsion_opls": "opls_torsion", "torsion_0.5_3cos-": "0.5compass_torsion",
                       "torsion_compass": "compass_torsion", "torsion_3cos-": "compass_torsion"}[kind]
            torsions.append(f"{{ types: [{', '.join(t)}], potential: {pot}, parameters: {{ k1: {q(*p[0])}, k2: {q(*p[1])}, k3: {q(*p[2])} }} }}")
        elif kind == "impropre_0.5_harm":
            t, c0, k = w[1:5], (w[5], read_unit(w[6])), (w[7], read_unit(w[8]))
            impropers.append(f"{{ types: [{', '.join(t)}], potential: harm_improper, parameters: {{ chi0: {q(*c0)}, k: {q(*k)} }} }}")
        elif kind == "impropre_opls":
            t, p = w[1:5], [(w[5], read_unit(w[6])), (w[7], read_unit(w[8])), (w[9], read_unit(w[10]))]
            impropers.append(f"{{ types: [{', '.join(t)}], potential: opls_improper, parameters: {{ k1: {q(*p[0])}, k2: {q(*p[1])}, k3: {q(*p[2])} }} }}")
        elif kind == "impropre_gaff":
            t, u = w[1:5], read_unit(w[6])
            impropers.append(f"{{ types: [{', '.join(t)}], potential: opls_improper, parameters: {{ k1: 0.0 {u}, k2: {q(2.0 * float(w[5]), u)}, k3: 0.0 {u} }} }}")
        elif kind.startswith(("bond_", "angle_", "torsion_", "impropre_")):
            fail(f"intramolecular potential {kind} has no input file form")

    # ---------------- intramolecular pair weights ----------------
    weights = []
    if stop is not None:
        n = int(words[stop + 1][0])
        for w in words[stop + 2:stop + 2 + n]:
            a, b = sp_name(w[0]), sp_name(w[1])
            weights.append(f"{a}_{b}: {{ bond: {w[2]}, bend: {w[3]}, torsion: {w[4]}, bond_rf: {w[5]}, bend_rf: {w[6]}, torsion_rf: {w[7]} }}")

    # ---------------- force field file ----------------
    out = [f"# Flexible molecule force field converted from {args.fatomes} by scripts/python/fatomes_to_yaml.py.",
           "# Use it with read_xyz_molecules followed by molecule_force_field in setup_system.",
           "species:"]
    for (fname, ff, m, c), n in zip(species, names):
        ffkey = f", ff: {ff}" if ff != n else ""
        out.append(f"  - {n}: {{ mass: {m}, charge: {c}, molecule: {args.molecule}{ffkey} }}")
    out += ["", "molecule_force_field:", f"  rf_self_pairs: {str(rf_self_pairs).lower()}",
            f"  long_range_correction: {str(long_range_correction).lower()}"]
    for key, items in (("potentials_for_bonds", bonds), ("potentials_for_angles", bends),
                       ("potentials_for_torsions", torsions), ("potentials_for_impropers", impropers)):
        if items:
            out.append(f"  {key}:")
            out += ["    - " + it for it in items]
    if weights:
        out.append("  mol_pair_weights:")
        out += ["    " + it for it in weights]
    if pairs:
        out.append("  potentials_for_pairs:")
        for a, b, p in pairs:
            fields = []
            if "rcut" in p:
                fields.append(f"rcut: {q(*p['rcut'])}")
            for kind in ("exp6", "lj"):
                if kind in p:
                    fields.append(f"{kind}: {{ " + ", ".join(f"{k}: {q(*v)}" for k, v in p[kind].items()) + " }")
            if "rf" in p:
                fields.append(f"rf: {{ epsilon: {p['rf'][0]!r}, rc: {p['rf'][1]!r} ang }}")
            out.append(f"    - {{ type_a: {a}, type_b: {b}, potential: {{ {', '.join(fields)} }} }}")
    text = "\n".join(out) + "\n"
    if args.output == "-":
        sys.stdout.write(text)
    else:
        open(args.output, "w").write(text)

    # ---------------- system in read_xyz_molecules format ----------------
    if args.xyz:
        k = section("PositionDesAtomesCart")
        if k is None:
            fail("--xyz needs a PositionDesAtomesCart section")
        unit = words[k][1] if len(words[k]) > 1 else "ang"
        f = LENGTH_TO_ANG.get(unit)
        if f is None:
            fail(f"unsupported position unit {unit}")
        na = int(words[k + 1][0])
        atoms = words[k + 2:k + 2 + na]
        # cell : columns are the cell vectors, as read_fatomes_mol builds it
        i = section("maille_long")
        lx, ly, lz = (float(x) for x in words[i][1:4])
        i = section("maille_angle")
        right = all(float(x) == 90.0 for x in words[i][1:4])
        al, bt, gm = (float(x) * math.pi / 180.0 for x in words[i][1:4])
        if right:
            H = [[lx, 0.0, 0.0], [0.0, ly, 0.0], [0.0, 0.0, lz]]
        else:
            a = [1.0, 0.0, 0.0]
            b = [math.cos(gm), math.sin(gm), 0.0]
            cy = (math.cos(al) - math.cos(gm) * math.cos(bt)) / math.sin(gm)
            c = [math.cos(bt), cy, math.sqrt(1.0 - math.cos(bt) ** 2 - cy ** 2)]
            H = [[a[r] * lx, b[r] * ly, c[r] * lz] for r in range(3)]
        # bonds (0-based atom ids, as read_fatomes_mol)
        con = {}
        z = section("Zmatrice")
        if z is not None:
            for w in words[z + 2:z + 2 + int(words[z + 1][0])]:
                c = [int(x) for x in w[1:5]]
                con[int(w[0])] = c + [-1] * (4 - len(c))
        # per atom charges
        charge = {}
        z = section("ModificationChargeDesAtomes")
        if z is not None:
            for w in words[z + 2:z + 2 + int(words[z + 1][0])]:
                charge[int(w[0])] = float(w[1])
        xyz = [str(na), " ".join(repr(H[r][s]) for r in range(3) for s in range(3))]
        for idx, w in enumerate(atoms):
            if w[0] not in fat_index:
                fail(f"unknown atom name {w[0]} in positions")
            sp = species[fat_index[w[0]]]
            qa = charge.get(idx, float(sp[3].split()[0]))
            c = con.get(idx, [-1, -1, -1, -1])
            xyz.append(f"{args.molecule} {sp_name(w[0])} 0 {idx} " + " ".join(repr(float(x) * f) for x in w[1:4])
                       + " " + " ".join(str(x) for x in c) + f" -1 {qa!r}")
        open(args.xyz, "w").write("\n".join(xyz) + "\n")


if __name__ == "__main__":
    main()
