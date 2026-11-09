"""DPA3 (repflows) for the exaStamp native DeePMD kernels: export of a deepmd DPA3 descriptor to the native format, and a
numpy evaluation of the native algorithm (whole periodic system, message passing between atoms) for --check.

Model (per layer, static neighbor selection; the dynamic selection only changes the normalizations):
  node   n'  = n + r0 act(Wns n) + r1 act(Wsym [grrg(e sw, h2), grrg(n_j sw, h2)]) + r2 sum_j act(Wne [n, n_j, e]) sw / nnei
  edge   e'  = e + r0 act(Wes [n, n_j, e]) + r1 act(W2 sum_k a_sw_j a_sw_k act(W1 [a, nc, ec_ik, ec_ij]) / sqrt(a_sel))
  angle  a'  = a + r0 act(Was [a, nc, ec_ik, ec_ij])
with grrg(g, h) = (H^T G)[:, :axis]^T (H^T G) / 3, H^T G = sum_j h2_j x g_j / sqrt(nnei); nc, ec the compressed node / edge
features (first channels, or a linear map). Only the neighbors' node features n_j cross atoms.
"""
import numpy as np


def fail(msg):
    raise SystemExit(f"deepmd_export (dpa3): {msg}")


SUPPORTED_ACT = ("silu", "tanh", "relu", "gelu", "gelu_tf", "sigmoid", "softplus", "linear", "none")


def _act_name(a):
    a = a.lower()
    if a.startswith("silut") or a.startswith("custom_silu"):
        return "silut", float(a.split(":")[-1]) if ":" in a else 3.0
    if a not in SUPPORTED_ACT:
        fail(f"activation '{a}' is not supported")
    return a, 0.0


def act(name, T, x):
    if name == "silut":
        sg = lambda z: 1 / (1 + np.exp(-z))
        slope = sg(T) + T * sg(T) * (1 - sg(T))
        const = T * sg(T)
        return np.where(x < T, x * sg(x), np.tanh(slope * (x - T)) + const)
    if name == "silu":
        return x / (1 + np.exp(-x))
    if name == "tanh":
        return np.tanh(x)
    if name == "relu":
        return np.maximum(x, 0)
    if name in ("gelu", "gelu_tf"):
        return 0.5 * x * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x ** 3)))
    if name == "sigmoid":
        return 1 / (1 + np.exp(-x))
    if name == "softplus":
        return np.log1p(np.exp(x))
    return x


def _lin(layer):
    if layer is None:
        return None
    w = np.asarray(layer.w, dtype=np.float64)
    b = None if layer.b is None else np.asarray(layer.b, dtype=np.float64)
    return {"w": w.ravel().tolist(), "nin": int(w.shape[0]), "nout": int(w.shape[1]), "b": None if b is None else b.tolist()}


def export_dpa3(d, nt):
    from deepmd.dpmodel.descriptor.base_descriptor import BaseDescriptor
    obj = BaseDescriptor.deserialize(d)
    blk = obj.repflows
    if obj.add_chg_spin_ebd:
        fail("charge / spin embeddings are not supported")
    if blk.exclude_types:
        fail("excluded types are not supported")
    if blk.update_style != "res_residual":
        fail(f"update_style '{blk.update_style}' is not supported (res_residual)")
    if getattr(blk, "sequential_update", False):
        fail("sequential_update is not supported")
    if blk.n_multi_edge_message != 1:
        fail("n_multi_edge_message > 1 is not supported")
    if not blk.update_angle:
        fail("update_angle false is not supported")
    if blk.use_dynamic_sel and not blk.smooth_edge_update:
        fail("dynamic selection requires smooth_edge_update")
    an, aT = _act_name(obj.activation_function)
    davg = np.asarray(blk.mean, dtype=np.float64)
    dstd = np.asarray(blk.stddev, dtype=np.float64)
    for name, a in (("davg", davg), ("dstd", dstd)):
        if not np.allclose(a, a[:, :1, :], rtol=0, atol=1e-12 * max(1.0, np.abs(a).max())):
            fail(f"{name} varies over the neighbor slots: not supported")
    tebd = np.asarray(obj.type_embedding.call(), dtype=np.float64)[:nt]
    layers = []
    for L in blk.layers:
        layers.append({
            "node_self": _lin(L.node_self_mlp), "node_sym": _lin(L.node_sym_linear), "node_edge": _lin(L.node_edge_linear),
            "edge_self": _lin(L.edge_self_linear), "edge_angle1": _lin(L.edge_angle_linear1), "edge_angle2": _lin(L.edge_angle_linear2),
            "angle_self": _lin(L.angle_self_linear),
            "compress_n": _lin(L.a_compress_n_linear), "compress_e": _lin(L.a_compress_e_linear),
            "n_residual": [np.asarray(r, np.float64).ravel().tolist() for r in L.n_residual],
            "e_residual": [np.asarray(r, np.float64).ravel().tolist() for r in L.e_residual],
            "a_residual": [np.asarray(r, np.float64).ravel().tolist() for r in L.a_residual],
        })
    L0 = blk.layers[0]
    out = {
        "descriptor": "dpa3", "rcut": float(blk.e_rcut), "rcut_smth": float(blk.e_rcut_smth), "protection": float(blk.env_protection),
        "sel": [int(blk.nnei)], "a_rcut": float(blk.a_rcut), "a_rcut_smth": float(blk.a_rcut_smth), "a_sel": int(blk.a_sel),
        "n_dim": int(blk.n_dim), "e_dim": int(blk.e_dim), "a_dim": int(blk.a_dim), "axis_neuron": int(blk.axis_neuron),
        "activation": an, "act_threshold": aT,
        "use_exp_switch": bool(blk.use_exp_switch), "edge_init_use_dist": bool(blk.edge_init_use_dist),
        "smooth_edge_update": bool(blk.smooth_edge_update), "use_dynamic_sel": bool(blk.use_dynamic_sel),
        "sel_reduce_factor": float(blk.sel_reduce_factor),
        "a_compress_rate": int(blk.a_compress_rate), "a_compress_use_split": bool(blk.a_compress_use_split),
        "n_a_compress_dim": int(L0.n_a_compress_dim), "e_a_compress_dim": int(L0.e_a_compress_dim),
        "concat_output_tebd": bool(obj.concat_output_tebd),
        "tebd": tebd.tolist(), "davg": davg[:, 0, :].tolist(), "dstd": dstd[:, 0, :].tolist(),
        "edge_embd": _lin(blk.edge_embd), "angle_embd": _lin(blk.angle_embd), "layers": layers,
    }
    dim = out["n_dim"] + (tebd.shape[1] if out["concat_output_tebd"] else 0)
    return out, dim


# --------------------------------------------------------------------------- numpy evaluation

def _apply(lin, x):
    w = np.asarray(lin["w"]).reshape(lin["nin"], lin["nout"])
    y = x @ w
    return y + np.asarray(lin["b"]) if lin["b"] is not None else y


def _poly_sw(r, rs, rc):
    u = np.clip((r - rs) / (rc - rs), 0, 1)
    return u ** 3 * (-6 * u ** 2 + 15 * u - 10) + 1


def _exp_sw(r, rs, rc):
    rb = np.clip(r, 0.0, rc)
    return np.exp(-np.exp(20.0 / rs * (rb - rs)))


def dpa3_descriptors(m, coord, box, atype):
    """descriptors of all atoms of an orthorhombic periodic box"""
    N = len(coord)
    Lb = np.diag(box)
    A = lambda x: act(m["activation"], m["act_threshold"], x)
    sw = _exp_sw if m["use_exp_switch"] else _poly_sw
    rc, rs, ac, ars = m["rcut"], m["rcut_smth"], m["a_rcut"], m["a_rcut_smth"]
    nnei, asel, axis, p = m["sel"][0], m["a_sel"], m["axis_neuron"], m["protection"]
    davg, dstd = np.asarray(m["davg"]), np.asarray(m["dstd"])
    if m["use_dynamic_sel"]:
        f = m["sel_reduce_factor"]
        s_sym, s_msg, s_ang = (nnei / f) ** -0.5, 1.0 / (nnei / f), (asel / f) ** -0.5
    else:
        s_sym, s_msg, s_ang = nnei ** -0.5, 1.0 / nnei, asel ** -0.5
    # neighbors (nearest nnei within rcut), geometry, initial features
    nb, D, R, SW, H2, E, AN, ASW, ANG = [], [], [], [], [], [], [], [], []
    for i in range(N):
        d = coord - coord[i]
        d -= Lb * np.round(d / Lb)
        r = np.linalg.norm(d, axis=1)
        idx = [j for j in np.argsort(r, kind="stable") if 0 < r[j] < rc][:nnei]
        nb.append(np.array(idx, dtype=int))
        dd, rr = d[idx], r[idx]
        w = sw(rr, rs, rc)
        q = rr + p
        ti = atype[i]
        Rn = (np.concatenate([(w / q)[:, None], dd * (w / q ** 2)[:, None]], axis=1) - davg[ti]) / dstd[ti]
        SW.append(w); H2.append(Rn[:, 1:])
        ein = rr[:, None] if m["edge_init_use_dist"] else Rn[:, :1]
        e0 = _apply(m["edge_embd"], ein)
        E.append(e0 if m["edge_init_use_dist"] else A(e0))
        na = [k for k in range(min(asel, len(idx))) if rr[k] < ac]
        AN.append(np.array(na, dtype=int))
        ASW.append(sw(rr[na], ars, ac))
        nv = dd[na] / (rr[na][:, None] + 1e-6)
        cos = (nv @ nv.T) * (1 - 1e-6)
        ANG.append(_apply(m["angle_embd"], (cos / np.sqrt(np.pi))[..., None]))
    tebd = np.asarray(m["tebd"])
    node = A(tebd[atype])
    def grrg(hg):
        return (hg[:, :axis].T @ hg / 3.0).ravel()
    for L in m["layers"]:
        new_node = np.zeros_like(node)
        new_E, new_ANG = [], []
        rn, re, ra = [np.asarray(x) for x in L["n_residual"]], [np.asarray(x) for x in L["e_residual"]], [np.asarray(x) for x in L["a_residual"]]
        for i in range(N):
            n, e, a, j = node[i], E[i], ANG[i], nb[i]
            nj, w, h2 = node[j], SW[i], H2[i]
            u_self = A(_apply(L["node_self"], n))
            hg_e = h2.T @ (e * w[:, None]) * s_sym
            hg_n = h2.T @ (nj * w[:, None]) * s_sym
            u_sym = A(_apply(L["node_sym"], np.concatenate([grrg(hg_e), grrg(hg_n)])))
            info = np.concatenate([np.repeat(n[None], len(j), 0), nj, e], axis=1)
            u_msg = (A(_apply(L["node_edge"], info)) * w[:, None]).sum(0) * s_msg
            new_node[i] = n + rn[0] * u_self + rn[1] * u_sym + rn[2] * u_msg
            e_self = A(_apply(L["edge_self"], info))
            na = AN[i]
            if m["a_compress_rate"] == 0:
                nc, ec = n, e[na]
            elif m["a_compress_use_split"]:
                nc, ec = n[:m["n_a_compress_dim"]], e[na, :m["e_a_compress_dim"]]
            else:
                nc, ec = _apply(L["compress_n"], n), _apply(L["compress_e"], e[na])
            k = len(na)
            ainfo = np.concatenate([a, np.broadcast_to(nc, (k, k, len(nc))), np.broadcast_to(ec[None, :, :], (k, k, ec.shape[1])),
                                    np.broadcast_to(ec[:, None, :], (k, k, ec.shape[1]))], axis=2)   # [j][k]: (a, nc, ec_k, ec_j)
            upd = A(_apply(L["edge_angle1"], ainfo))
            aw = ASW[i]
            red = (aw[:, None, None] * aw[None, :, None] * upd).sum(1) * s_ang
            full = np.zeros((len(j), m["e_dim"]))
            if not m["smooth_edge_update"]:
                full[:] = e
            full[na] = red
            e_ang = A(_apply(L["edge_angle2"], full))
            new_E.append(e + re[0] * e_self + re[1] * e_ang)
            new_ANG.append(a + ra[0] * A(_apply(L["angle_self"], ainfo)))
        node, E, ANG = new_node, new_E, new_ANG
    if m["concat_output_tebd"]:
        node = np.concatenate([node, tebd[atype]], axis=1)
    return node
