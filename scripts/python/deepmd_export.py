#!/usr/bin/env python3
"""
Export a DeePMD-kit model to the native exaStamp format read by deepmd_init.

    python deepmd_export.py model.pth model.dpnative.yaml [--check frame.npz [--reference model.pth]]

Runs in a Python environment where deepmd-kit is installed (any backend that
can load the model: .pth, .pt2, .dp, .yaml...).  The model is read through
deepmd's backend-neutral serialization, so the output does not depend on the
training backend.

Supported (v1): standard energy model, descriptors se_e2_a (type_one_side true
or false), se_e2_r (se_r) and se_e3 (se_t), 'ener' fitting with one network per
type or a single network (mixed_types), any activation listed in ACTIVATIONS.  Anything else (other
descriptors, fparam/aparam, excluded types, vacuum references, compression)
is refused rather than exported with different physics.

--check evaluates the exported model with plain numpy, following exactly the
algorithm of the exaStamp kernel, and compares energies to deepmd's own
DeepPot on one periodic configuration (npz with coord (N,3), box (3,3), atype (N,)
as model type indices, in Angstrom).
"""

import argparse
import sys

import numpy as np
import yaml

ACTIVATIONS = ("tanh", "relu", "gelu", "gelu_tf", "softplus", "sigmoid", "silu", "none", "linear")


def fail(msg):
    sys.exit(f"deepmd_export: {msg}")


def _layers(net):
    out = []
    for layer in net["layers"]:
        var = layer["@variables"]
        act = layer.get("activation_function") or "none"
        if act not in ACTIVATIONS:
            fail(f"activation '{act}' is not supported")
        w = np.asarray(var["w"], dtype=np.float64)
        b = var.get("b")
        idt = var.get("idt")
        out.append({
            "nin": int(w.shape[0]), "nout": int(w.shape[1]), "activation": act,
            "resnet": bool(layer.get("resnet", False)),
            "w": w.ravel().tolist(),                       # row-major, nin x nout
            "b": np.zeros(w.shape[1]).tolist() if b is None else np.asarray(b, dtype=np.float64).ravel().tolist(),
            **({"idt": np.asarray(idt, dtype=np.float64).ravel().tolist()} if idt is not None else {}),
        })
    return {"layers": out}


def _collection(coll, ntypes, ndim):
    """Networks of a NetworkCollection, flattened by type index (ndim 0, 1 or 2)."""
    if coll["ndim"] != ndim:
        fail(f"unexpected network collection with ndim {coll['ndim']} (expected {ndim})")
    nets = coll["networks"]
    expect = ntypes ** ndim
    if len(nets) != expect or any(n is None for n in nets):
        fail(f"network collection has {len(nets)} networks, {expect} expected")
    if ndim == 2:   # deepmd flattens key (a, b) as a + b * ntypes, the native file as a * ntypes + b
        nets = [nets[a + b * ntypes] for a in range(ntypes) for b in range(ntypes)]
    return [_layers(n) for n in nets]


# deepmd descriptor type -> native name
DESCRIPTORS = {"se_e2_a": "se_e2_a", "se_r": "se_e2_r", "se_e2_r": "se_e2_r", "se_e3": "se_e3", "se_t": "se_e3",
               "dpa1": "se_atten", "se_atten": "se_atten", "se_atten_v2": "se_atten",
               "se_e3_tebd": "se_t_tebd", "se_t_tebd": "se_t_tebd"}


def _mixed_stats(blk, name):
    """per center type statistics of a mixed neighbor list (constant over the slots), repeated per neighbor type"""
    a = np.asarray(getattr(blk, name), dtype=np.float64)
    if not np.allclose(a, a[:, :1, :], rtol=0, atol=1e-12 * max(1.0, np.abs(a).max())):
        fail(f"{name} varies over the neighbor slots: not supported")
    return a[:, 0, :]


def _fold_concat(first, rest, tebd, pairs):
    """embedding networks of a concat type-embedding input [x, tebd_a, tebd_b, ...]: the type part folded into the
    first layer bias, one network per entry of pairs (lists of type embedding rows)"""
    w = np.asarray(first["w"]).reshape(first["nin"], first["nout"])
    if first["resnet"] and first["nout"] in (first["nin"], 2 * first["nin"]):
        fail("embedding with a residual first layer is not supported")
    nets = []
    for t in pairs:
        t = np.concatenate(t)
        if 1 + len(t) != first["nin"]:
            fail("unexpected embedding input width")
        f = dict(first, nin=1, resnet=False, w=w[:1].ravel().tolist(), b=(np.asarray(first["b"]) + t @ w[1:]).tolist())
        nets.append({"layers": [f] + rest})
    return nets


def _se_t_tebd(d, nt, blk=None, tebd_table=None, concat_output_tebd=None):
    """se_e3_tebd / se_t_tebd: three-body angles over the pairs of a mixed neighbor list, type embeddings of both
    neighbors (concat: folded into a per type pair first layer bias; strip: type pair table, times sw_j sw_k if smooth)"""
    if blk is None:
        from deepmd.dpmodel.descriptor.base_descriptor import BaseDescriptor
        obj = BaseDescriptor.deserialize(d)
        blk, tebd_table, concat_output_tebd = obj.se_ttebd, obj.type_embedding.call(), obj.concat_output_tebd
    if blk.exclude_types:
        fail("se_t_tebd with excluded types is not supported")
    mode = blk.tebd_input_mode
    if mode not in ("concat", "strip"):
        fail(f"tebd_input_mode '{mode}' is not supported")
    tebd = np.asarray(tebd_table, dtype=np.float64)[:nt]
    davg, dstd = _mixed_stats(blk, "mean"), _mixed_stats(blk, "stddev")
    if np.any(davg[:, 1:] != 0.0):
        fail("se_t_tebd with non-zero angular averages is not supported (empty slots would depend on the neighbor type)")
    emb = _layers(blk.embeddings.serialize()["networks"][0])
    if mode == "concat":
        embedding = _fold_concat(emb["layers"][0], emb["layers"][1:], tebd, [(tebd[a], tebd[b]) for a in range(nt) for b in range(nt)])
        strip = None
    else:
        embedding = [emb]
        pairs = np.array([np.concatenate([tebd[a], tebd[b]]) for a in range(nt) for b in range(nt)])
        strip = np.asarray(blk.cal_g_strip(pairs, 0), dtype=np.float64)          # [tj * nt + tk][ng]
    return {
        "nnei": int(blk.nnei), "sel": [int(blk.nnei)], "embedding": embedding,
        "davg": np.repeat(davg[:, None, :], nt, axis=1), "dstd": np.repeat(dstd[:, None, :], nt, axis=1),
        "extra": {"tebd": tebd.tolist(), "tebd_input_mode": mode, "concat_output_tebd": bool(concat_output_tebd),
                  "smooth": bool(blk.smooth), **({"strip_table": strip.tolist()} if strip is not None else {})},
    }


def _atten(d, nt, blk=None, tebd_table=None, concat_output_tebd=None):
    """DPA1 / se_atten(_v2): mixed neighbor list, type embedding tables, gated self-attention layers.
    The type parts of the embedding input are tabulated with deepmd's own networks: concat mode folds them into a
    per type pair first-layer bias, strip mode precomputes the type pair factor of the embedding."""
    if blk is None:
        from deepmd.dpmodel.descriptor.base_descriptor import BaseDescriptor
        obj = BaseDescriptor.deserialize(d)
        blk, tebd_table, concat_output_tebd = obj.se_atten, obj.type_embedding.call(), obj.concat_output_tebd
    if getattr(blk, "lmax", 1) != 1:
        fail("DPA1 with lmax > 1 (moment basis) is not supported")
    if blk.exclude_types or getattr(blk, "geo_compress", False) or getattr(blk, "tebd_compress", False):
        fail("DPA1 with excluded types or compression is not supported")
    mode = blk.tebd_input_mode
    if mode not in ("concat", "strip"):
        fail(f"tebd_input_mode '{mode}' is not supported")
    tebd = np.asarray(tebd_table, dtype=np.float64)[:nt]      # real types (the last row is the padding type)
    tdim = tebd.shape[1]
    nnei = int(blk.nnei)
    davg = np.asarray(blk.mean, dtype=np.float64)
    dstd = np.asarray(blk.stddev, dtype=np.float64)
    for name, a in (("davg", davg), ("dstd", dstd)):
        if not np.allclose(a, a[:, :1, :], rtol=0, atol=1e-12 * max(1.0, np.abs(a).max())):
            fail(f"DPA1 {name} varies over the neighbor slots: not supported")
    davg, dstd = davg[:, 0, :], dstd[:, 0, :]
    # empty neighbor slots take the type of extended atom 0 in deepmd: only type independent if they are zero rows
    # (zero averages) or their type factor vanishes (strip mode, smooth)
    if np.any(davg != 0.0) and not (mode == "strip" and blk.smooth):
        fail("DPA1 with non-zero averages (davg) is only supported in strip mode with smooth_type_embedding (train with set_davg_zero)")
    one_side = bool(blk.type_one_side)
    emb = _layers(blk.embeddings.serialize()["networks"][0])
    if mode == "concat":
        first = emb["layers"][0]
        w = np.asarray(first["w"]).reshape(first["nin"], first["nout"])
        tin = tdim if one_side else 2 * tdim
        if first["nin"] != 1 + tin:
            fail("unexpected DPA1 embedding input width")
        if first["resnet"] and first["nout"] in (first["nin"], 2 * first["nin"]):
            fail("DPA1 embedding with a residual first layer is not supported")
        nets = []
        for a in range(nt):
            for b in range(nt):
                if one_side and a > 0:
                    break
                t = tebd[b] if one_side else np.concatenate([tebd[b], tebd[a]])   # input [s, tebd_j, tebd_i]
                f = dict(first, nin=1, resnet=False, w=w[:1].ravel().tolist(), b=(np.asarray(first["b"]) + t @ w[1:]).tolist())
                nets.append({"layers": [f] + emb["layers"][1:]})
        embedding, strip = nets, None
    else:
        embedding = [emb]
        if one_side:
            strip = np.asarray(blk.cal_g_strip(tebd, 0), dtype=np.float64)                     # [neighbor type][ng]
        else:
            pairs = np.array([np.concatenate([tebd[b], tebd[a]]) for a in range(nt) for b in range(nt)])
            strip = np.asarray(blk.cal_g_strip(pairs, 0), dtype=np.float64)                    # [center * nt + neighbor][ng]
    layers = []
    for lay in blk.dpa1_attention.attention_layers:
        at, ln = lay.attention_layer, lay.attn_layer_norm
        if at.num_heads != 1:
            fail("DPA1 attention with several heads is not supported")
        layers.append({"in_w": np.asarray(at.in_proj.w, np.float64).ravel().tolist(), "in_b": np.asarray(at.in_proj.b, np.float64).tolist(),
                       "out_w": np.asarray(at.out_proj.w, np.float64).ravel().tolist(), "out_b": np.asarray(at.out_proj.b, np.float64).tolist(),
                       "ln_w": np.asarray(ln.w, np.float64).tolist(), "ln_b": np.asarray(ln.b, np.float64).tolist(), "ln_eps": float(ln.eps)})
    att = blk.dpa1_attention
    a0 = att.attention_layers[0].attention_layer if att.attention_layers else None
    return {
        "nnei": nnei, "sel": [nnei], "one_side": one_side, "embedding": embedding,
        "davg": np.repeat(davg[:, None, :], nt, axis=1), "dstd": np.repeat(dstd[:, None, :], nt, axis=1),
        "extra": {
            "tebd": tebd.tolist(), "tebd_input_mode": mode, "concat_output_tebd": bool(concat_output_tebd),
            **({"strip_table": strip.tolist()} if strip is not None else {}),
            "attention": {"layers": layers, "hidden": int(att.hidden_dim), "scaling": float(a0.scaling) if a0 else 1.0,
                          "normalize": bool(att.normalize), "smooth": bool(blk.smooth), "dotr": bool(att.dotr), "shift": 20.0},
        },
    }


def _per_type_stats(name, arr, sel, fill):
    """davg / dstd: (ntypes, nnei, 4 or 1) -> (ntypes, ntypes, 4), constant over the slots of a neighbor type.
    A radial-only array (se_e2_r) is completed with ``fill``."""
    arr = np.asarray(arr, dtype=np.float64)
    if arr.shape[-1] == 1:
        arr = np.concatenate([arr, np.full(arr.shape[:-1] + (3,), fill)], axis=-1)
    out = np.zeros((arr.shape[0], len(sel), 4))
    start = 0
    for t, s in enumerate(sel):
        block = arr[:, start:start + s, :]
        if not np.allclose(block, block[:, :1, :], rtol=0, atol=1e-12 * max(1.0, np.abs(block).max())):
            fail(f"{name} is not constant over the neighbor slots of type {t}: not supported")
        out[:, t, :] = block[:, 0, :]
        start += s
    return out


def block_descriptor(kind, blk, tebd_table, nt):
    """native fields of a DPA2 sub-block: kind 'se_atten' (repinit) or 'se_t_tebd' (three-body), no type embedding output"""
    A = _atten(None, nt, blk, tebd_table, False) if kind == "se_atten" else _se_t_tebd(None, nt, blk, tebd_table, False)
    ng = A["embedding"][0]["layers"][-1]["nout"]
    axis = int(blk.axis_neuron) if kind == "se_atten" else 1
    return {"descriptor": kind, "rcut": float(blk.rcut), "rcut_smth": float(blk.rcut_smth), "protection": float(blk.env_protection),
            "sel": A["sel"], "axis_neuron": axis, "type_one_side": A.get("one_side", False),
            "davg": A["davg"].tolist(), "dstd": A["dstd"].tolist(), "embedding": A["embedding"], **A["extra"]}, ng * axis


def _descriptor(d, nt):
    """the native fields of one (non-hybrid) descriptor, and its output width"""
    descriptor = DESCRIPTORS.get(d.get("type"))
    if descriptor is None:
        fail(f"descriptor '{d.get('type')}' is not supported yet ({', '.join(sorted(set(DESCRIPTORS) | {'hybrid'}))})")
    if d.get("exclude_types") or d.get("spin") or d.get("compress"):
        fail(f"{descriptor} with excluded types, spin or compression is not supported")
    if d.get("env_mat", {}).get("use_exp_switch") or d.get("use_exp_switch"):
        fail("exponential switch is not supported")
    extra = {}
    if descriptor == "se_atten":
        A = _atten(d, nt)
        sel, one_side, embeddings, extra = A["sel"], A["one_side"], A["embedding"], A["extra"]
    elif descriptor == "se_t_tebd":
        A = _se_t_tebd(d, nt)
        sel, one_side, embeddings, extra = A["sel"], False, A["embedding"], A["extra"]
    else:
        sel = [int(s) for s in d["sel"]]
        one_side = descriptor == "se_e2_r" or ( descriptor == "se_e2_a" and bool(d.get("type_one_side", True)) )
        embeddings = _collection(d["embeddings"], nt, 1 if one_side else 2)
    if descriptor in ("se_atten", "se_t_tebd"):
        davg, dstd = A["davg"], A["dstd"]
    else:
        davg = _per_type_stats("davg", d["@variables"]["davg"], sel, 0.0)
        dstd = _per_type_stats("dstd", d["@variables"]["dstd"], sel, 1.0)
    if descriptor == "se_e3" and np.any(davg[:, :, 1:] != 0.0):
        fail("se_e3 with non-zero angular averages (davg[..., 1:]) is not supported")
    ng = embeddings[0]["layers"][-1]["nout"]
    axis = int(d["axis_neuron"]) if descriptor in ("se_e2_a", "se_atten") else 1    # se_t_tebd: ng values
    dim = ng * axis + (len(extra["tebd"][0]) if extra.get("concat_output_tebd") else 0)
    return {
        "descriptor": descriptor,
        "rcut": float(d["rcut"]),
        "rcut_smth": float(d["rcut_smth"]),
        "protection": float(d.get("env_protection", 0.0)),
        "sel": sel,
        "axis_neuron": axis,
        "type_one_side": one_side,
        "davg": davg.tolist(),          # [center type][neighbor type][4]
        "dstd": dstd.tolist(),
        "embedding": embeddings,        # by neighbor type (one side), else type pair a*ntypes+b (se_e2_a: center, neighbor)
        **extra,                        # se_atten / se_t_tebd: tebd, tebd_input_mode, concat_output_tebd, strip_table, attention
    }, dim


def export(model_file):
    from deepmd.backend.backend import Backend
    data = Backend.detect_backend_by_model(model_file)().serialize_hook(model_file)
    m = data["model"]
    if m.get("type", "standard") != "standard" or m.get("@class") != "Model":
        fail(f"model type '{m.get('type')}' is not supported (standard energy models only)")
    if m.get("atom_exclude_types") or m.get("pair_exclude_types"):
        fail("excluded types are not supported")
    d, f = m["descriptor"], m["fitting"]
    if f.get("type") != "ener" or f.get("dim_out", 1) != 1:
        fail(f"fitting '{f.get('type')}' is not supported (ener only)")
    if f.get("numb_fparam", 0) or f.get("numb_aparam", 0) or f.get("exclude_types") or f.get("spin"):
        fail("fitting with fparam, aparam, excluded types or spin is not supported")
    if f.get("vacuum_ref") or any(f.get("remove_vaccum_contribution") or []) or f.get("atom_ener"):
        fail("fitting with vacuum references or atom_ener is not supported")

    type_map = list(m["type_map"])
    nt = len(type_map)
    if d.get("type") == "dpa2":
        import os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from deepmd_dpa2 import export_dpa2
        desc, dim = export_dpa2(d, nt)
    elif d.get("type") == "dpa3":
        import os
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from deepmd_dpa3 import export_dpa3
        desc, dim = export_dpa3(d, nt)
    elif d.get("type") == "hybrid":
        # sub-descriptors on the hybrid neighbor list, each cut to its sel: neighbors beyond a sub-descriptor's rcut have
        # a zero switch, as empty slots (exception: attention without smoothing counts them as neighbors)
        subs = [_descriptor(x, nt) for x in d["list"]]
        rmax = max(x[0]["rcut"] for x in subs)
        for x, _ in subs:
            if x["descriptor"] == "se_atten" and not x["attention"]["smooth"] and x["rcut"] < rmax:
                fail("hybrid with a non-smooth DPA1 of a smaller rcut is not supported")
            if x["descriptor"] == "hybrid":
                fail("nested hybrid descriptors are not supported")
        desc = {"descriptor": "hybrid", "rcut": rmax, "list": [x for x, _ in subs]}
        dim = sum(n for _, n in subs)
    else:
        desc, dim = _descriptor(d, nt)
    mixed = bool(f.get("mixed_types", False))
    fittings = _collection(f["nets"], nt, 0 if mixed else 1)

    bias = np.asarray(f["@variables"]["bias_atom_e"], dtype=np.float64).reshape(nt, -1)[:, 0]
    out_bias = m.get("@variables", {}).get("out_bias")
    if out_bias is not None:
        bias = bias + np.asarray(out_bias, dtype=np.float64).reshape(-1, nt)[0]
    if fittings[0]["layers"][0]["nin"] != dim:
        fail(f"fitting input {fittings[0]['layers'][0]['nin']} != {dim}")

    return {
        "format": "exastamp-deepmd",
        "version": 1,
        "type_map": type_map,
        **desc,
        "fitting_mixed_types": mixed,
        "fitting": fittings,            # by center type, or a single network
        "energy_bias": bias.tolist(),   # eV, per type
        "units": {"energy": "eV", "length": "ang"},
    }


# --------------------------------------------------------------------------- numpy evaluation

def _act(name, x):
    if name == "tanh":
        return np.tanh(x)
    if name == "relu":
        return np.maximum(x, 0)
    if name in ("gelu", "gelu_tf"):
        return 0.5 * x * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x ** 3)))
    if name == "softplus":
        return np.log1p(np.exp(x))
    if name == "sigmoid":
        return 1 / (1 + np.exp(-x))
    if name == "silu":
        return x / (1 + np.exp(-x))
    return x


def _mlp(net, x):
    for L in net["layers"]:
        w = np.asarray(L["w"]).reshape(L["nin"], L["nout"])
        y = _act(L["activation"], x @ w + np.asarray(L["b"]))
        if "idt" in L:
            y = y * np.asarray(L["idt"])
        if L["resnet"] and L["nout"] == L["nin"]:
            y = y + x
        elif L["resnet"] and L["nout"] == 2 * L["nin"]:
            y = y + np.concatenate([x, x], axis=-1)
        x = y
    return x


def _switch(r, rs, rc):
    u = np.clip((r - rs) / (rc - rs), 0, 1)
    return u ** 3 * (-6 * u ** 2 + 15 * u - 10) + 1


def _atten_energy(model, ti, d, r, atype, p):
    """Per-atom energy of se_atten (DPA1): mixed list of the nnei nearest neighbors, empty slots handled explicitly."""
    A, nt, rc, nnei = model["attention"], len(model["davg"]), model["rcut"], model["sel"][0]
    davg, dstd = np.asarray(model["davg"])[ti, 0], np.asarray(model["dstd"])[ti, 0]
    tebd = np.asarray(model["tebd"])
    idx = [j for j in np.argsort(r) if r[j] < rc and r[j] > 0][:nnei]
    n, npad = len(idx), nnei - len(idx)
    rr = r[idx][:, None]
    sw = _switch(rr, model["rcut_smth"], rc)
    R = (np.concatenate([1 / (rr + p), d[idx] / (rr + p) ** 2], axis=1) * sw - davg) / dstd
    Rpad = (np.zeros(4) - davg) / dstd
    tj = atype[idx]
    one_side = model["type_one_side"]
    if model["tebd_input_mode"] == "concat":
        G = np.array([_mlp(model["embedding"][t if one_side else ti * nt + t], R[q, :1]) for q, t in enumerate(tj)])
        Gpad = None
    else:
        T = np.asarray(model["strip_table"])
        Tq = T[tj] if one_side else T[ti * nt + tj]
        gs = _mlp(model["embedding"][0], R[:, :1])
        G = gs * (Tq * sw if A["smooth"] else Tq) + gs
        Gpad = _mlp(model["embedding"][0], Rpad[None, :1])[0]     # smooth: the type factor vanishes (sw = 0)
    v = R[:, 1:]
    u = v / np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-12)
    sw1 = sw[:, 0]
    X = G
    for L in A["layers"]:
        ng = X.shape[1]
        hd = A["hidden"]
        P = X @ np.asarray(L["in_w"]).reshape(ng, 3 * hd) + np.asarray(L["in_b"])
        q, k, vv = P[:, :hd], P[:, hd:2 * hd], P[:, 2 * hd:]
        if A["normalize"]:
            q = q / np.linalg.norm(q, axis=1, keepdims=True)
            k = k / np.linalg.norm(k, axis=1, keepdims=True)
            vv = vv / np.linalg.norm(vv, axis=1, keepdims=True)
        lg = (q * A["scaling"]) @ k.T
        if A["smooth"]:
            S = sw1[:, None] * sw1[None, :]
            lg = (lg + A["shift"]) * S - A["shift"]
            mx = lg.max(axis=1, keepdims=True)
            e = np.exp(lg - mx)
            Z = e.sum(axis=1, keepdims=True) + npad * np.exp(-A["shift"] - mx)    # empty slots: logit -shift
            W = e / Z * S
        else:
            mx = lg.max(axis=1, keepdims=True)
            e = np.exp(lg - mx)
            W = e / e.sum(axis=1, keepdims=True)
        if A["dotr"]:
            W = W * (u @ u.T)
        out = (W @ vv) @ np.asarray(L["out_w"]).reshape(hd, ng) + np.asarray(L["out_b"])
        Y = X + out
        mu, var = Y.mean(axis=1, keepdims=True), Y.var(axis=1, keepdims=True)
        X = (Y - mu) / np.sqrt(var + L["ln_eps"]) * np.asarray(L["ln_w"]) + np.asarray(L["ln_b"])
        if Gpad is not None:    # empty slots: no attention (masked rows), residual + output bias + layer norm
            Yp = Gpad + np.asarray(L["out_b"])
            Gpad = (Yp - Yp.mean()) / np.sqrt(Yp.var() + L["ln_eps"]) * np.asarray(L["ln_w"]) + np.asarray(L["ln_b"])
    gr = X.T @ R
    if Gpad is not None:
        gr = gr + npad * np.outer(Gpad, Rpad)
    gr /= nnei
    m2 = model["axis_neuron"]
    D = (gr @ gr[:m2].T).ravel()
    if model["concat_output_tebd"]:
        D = np.concatenate([D, tebd[ti]])
    return D


def _t_tebd_energy(model, ti, d, r, atype, p):
    nt, rc, nnei = len(model["davg"]), model["rcut"], model["sel"][0]
    davg, dstd = np.asarray(model["davg"])[ti, 0], np.asarray(model["dstd"])[ti, 0]
    idx = [j for j in np.argsort(r) if r[j] < rc and r[j] > 0][:nnei]
    rr = r[idx][:, None]
    sw = _switch(rr, model["rcut_smth"], rc)
    R = (np.concatenate([1 / (rr + p), d[idx] / (rr + p) ** 2], axis=1) * sw - davg) / dstd
    tj = atype[idx]
    th = R[:, 1:] @ R[:, 1:].T
    n = len(idx)
    ng = model["embedding"][0]["layers"][-1]["nout"]
    D = np.zeros(ng)
    for j in range(n):
        for k in range(n):
            if model["tebd_input_mode"] == "concat":
                G = _mlp(model["embedding"][tj[j] * nt + tj[k]], th[j, k:k + 1])
            else:
                g = _mlp(model["embedding"][0], th[j, k:k + 1])
                T = np.asarray(model["strip_table"])[tj[j] * nt + tj[k]]
                G = g * (T * sw[j, 0] * sw[k, 0] if model["smooth"] else T) + g
            D += th[j, k] * G
    D /= nnei * nnei
    if model["concat_output_tebd"]:
        D = np.concatenate([D, np.asarray(model["tebd"])[ti]])
    return D


def native_descriptor(model, ti, d, r, atype):
    """descriptor of a central atom of type ti (d: minimum image vectors to all atoms, r: their lengths)"""
    if model["descriptor"] == "hybrid":
        return np.concatenate([native_descriptor(sub, ti, d, r, atype) for sub in model["list"]])
    nt, sel, rc = len(model["davg"]), model["sel"], model["rcut"]
    p = model["protection"]
    if model["descriptor"] == "se_atten":
        return _atten_energy(model, ti, d, r, atype, p)
    if model["descriptor"] == "se_t_tebd":
        return _t_tebd_energy(model, ti, d, r, atype, p)
    davg, dstd = np.asarray(model["davg"]), np.asarray(model["dstd"])
    ng = model["embedding"][0]["layers"][-1]["nout"]
    m2, nnei = model["axis_neuron"], sum(sel)
    rows = []   # normalized rows of the kept neighbors and empty slots, by neighbor type
    for t in range(nt):
        idx = [j for j in np.argsort(r) if r[j] > 0 and atype[j] == t and r[j] < rc][:sel[t]]
        R = np.zeros((len(idx), 4))
        if idx:
            rr = r[idx][:, None]
            sw = _switch(rr, model["rcut_smth"], rc)
            R = np.concatenate([1 / (rr + p), d[idx] / (rr + p) ** 2], axis=1) * sw
        npad = sel[t] - len(idx)
        rows.append((np.concatenate([R, np.zeros((npad, 4))]) - davg[ti, t]) / dstd[ti, t])
    desc = model["descriptor"]
    if desc == "se_e2_a":
        gr = np.zeros((ng, 4))
        for t in range(nt):
            emb = model["embedding"][t if model["type_one_side"] else ti * nt + t]
            gr += _mlp(emb, rows[t][:, :1]).T @ rows[t]
        gr /= nnei
        return (gr @ gr[:m2].T).ravel()
    if desc == "se_e2_r":
        return sum(_mlp(model["embedding"][t], rows[t][:, :1]).sum(0) for t in range(nt)) / (5.0 * nnei)
    D = np.zeros(ng)    # se_e3
    for a in range(nt):
        for b in range(a, nt):
            env = rows[a][:, 1:] @ rows[b][:, 1:].T
            G = _mlp(model["embedding"][a * nt + b], env.reshape(-1, 1))
            D += (env.reshape(-1, 1) * G).sum(0) / (sel[a] * sel[b])
    return D


def native_energies(model, coord, box, atype):
    """Per-atom energies (eV) with the exaStamp kernel algorithm, orthorhombic periodic box."""
    L = np.diag(box)
    energies = np.zeros(len(coord))
    if model["descriptor"] in ("dpa3", "dpa2"):
        if model["descriptor"] == "dpa3":
            from deepmd_dpa3 import dpa3_descriptors
            D = dpa3_descriptors(model, coord, box, atype)
        else:
            from deepmd_dpa2 import dpa2_descriptors
            D = dpa2_descriptors(model, coord, box, atype)
        for i in range(len(coord)):
            fit = model["fitting"][0 if model["fitting_mixed_types"] else atype[i]]
            energies[i] = _mlp(fit, D[i][None, :])[0, 0] + model["energy_bias"][atype[i]]
        return energies
    for i in range(len(coord)):
        ti = atype[i]
        d = coord - coord[i]
        d -= L * np.round(d / L)
        r = np.linalg.norm(d, axis=1)
        D = native_descriptor(model, ti, d, r, atype)
        fit = model["fitting"][0 if model["fitting_mixed_types"] else ti]
        energies[i] = _mlp(fit, D[None, :])[0, 0] + model["energy_bias"][ti]
    return energies


def check(model, model_file, frame):
    from deepmd.infer import DeepPot
    z = np.load(frame)
    coord, box, atype = z["coord"], z["box"], z["atype"]
    dp = DeepPot(model_file)
    e, f, v, ae, av = dp.eval(coord[None].reshape(1, -1), box.reshape(1, 9), atype, atomic=True)
    mine = native_energies(model, coord, box, atype)
    err = np.abs(mine - ae.ravel()).max() / np.abs(ae).max()
    print(f"energy: deepmd {e.ravel()[0]:.12f} eV, native numpy {mine.sum():.12f} eV; max per-atom relative diff {err:.2e}")
    if err > 1e-10:
        fail("the exported model does not reproduce deepmd's energies")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("output")
    ap.add_argument("--check", metavar="FRAME.npz", help="compare the exported model with deepmd on this configuration")
    ap.add_argument("--reference", metavar="MODEL", help="model file DeepPot loads for --check (default: the exported one)")
    args = ap.parse_args(argv)
    model = export(args.model)
    with open(args.output, "w") as fh:
        fh.write(f"# exaStamp native DeePMD model, exported from {args.model} by deepmd_export.py\n")
        yaml.safe_dump(model, fh, default_flow_style=None, sort_keys=False, width=1 << 20)
    print(f"wrote {args.output}: {model['descriptor']}, types {model['type_map']}, rcut {model['rcut']} ang")
    if args.check:
        check(model, args.reference or args.model, args.check)


if __name__ == "__main__":
    main()
