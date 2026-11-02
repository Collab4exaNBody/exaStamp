"""DPA2 for the exaStamp native DeePMD kernels: export of a deepmd DPA2 descriptor to the native format, and a numpy
evaluation of the native algorithm (whole periodic system, message passing between atoms) for --check.

  D0  = [ repinit(se_atten without attention), three-body(se_t_tebd) ]          per atom, its own neighbor lists
  g1  = act( T D0 (+ Tt tebd) ),  g2 = act( Wg2 s ), h2 (normalized env rows) on the repformer neighbors
  layer:  g1' = g1 + sum_m r_m u_m,  u = [ act(Wself g1), conv, act(W1 [grrg(g2), drrd(g1_j)]) ]  (deepmd's order)
          g2' = g2 + r0 act(W2 g2) + r1 LN( Wh attn(g2, h2) )                   (not in the last layer)
  descriptor = g1 (+ tebd)
"""
import numpy as np

from deepmd_dpa3 import act, _act_name, _lin, _apply, _poly_sw


def fail(msg):
    raise SystemExit(f"deepmd_export (dpa2): {msg}")


def export_dpa2(d, nt):
    from deepmd.dpmodel.descriptor.base_descriptor import BaseDescriptor
    import deepmd_export as de
    obj = BaseDescriptor.deserialize(d)
    if obj.exclude_types:
        fail("excluded types are not supported")
    tebd_table = np.asarray(obj.type_embedding.call(), dtype=np.float64)
    rep, dim_rep = de.block_descriptor("se_atten", obj.repinit, tebd_table, nt)
    three, dim_3b = (None, 0)
    if obj.use_three_body:
        three, dim_3b = de.block_descriptor("se_t_tebd", obj.repinit_three_body, tebd_table, nt)
    R = obj.repformers
    for flag in ("update_g1_has_attn", "update_g2_has_g1g1", "update_h2", "direct_dist"):
        if getattr(R, flag, False):
            fail(f"repformer option {flag} is not supported")
    if R.update_style != "res_residual":
        fail(f"update_style '{R.update_style}' is not supported (res_residual)")
    if not R.smooth:
        fail("non-smooth repformers are not supported")
    davg = np.asarray(R.mean, dtype=np.float64)
    dstd = np.asarray(R.stddev, dtype=np.float64)
    for name, a in (("davg", davg), ("dstd", dstd)):
        if not np.allclose(a, a[:, :1, :], rtol=0, atol=1e-12 * max(1.0, np.abs(a).max())):
            fail(f"repformer {name} varies over the neighbor slots: not supported")
    if np.any(davg[:, 0, 1:] != 0.0):
        fail("repformer with non-zero angular averages is not supported (empty slots would enter the attention gate)")
    an, aT = _act_name(R.activation_function)
    layers = []
    for L in R.layers:
        layers.append({
            "update_chnnl_2": bool(L.update_chnnl_2), "update_g2_has_attn": bool(L.update_g2_has_attn),
            "linear1": _lin(L.linear1), "linear2": _lin(L.linear2), "g1_self_mlp": _lin(L.g1_self_mlp), "proj_g1g2": _lin(L.proj_g1g2),
            "mapqk": _lin(L.attn2g_map.mapqk) if L.attn2g_map is not None else None,
            "mapv": _lin(L.attn2_mh_apply.mapv) if L.attn2_mh_apply is not None else None,
            "head_map": _lin(L.attn2_mh_apply.head_map) if L.attn2_mh_apply is not None else None,
            "attn2_lm": {"w": np.asarray(L.attn2_lm.w, np.float64).tolist(), "b": np.asarray(L.attn2_lm.b, np.float64).tolist(),
                         "eps": float(L.attn2_lm.eps)} if L.attn2_lm is not None else None,
            "g1_residual": [np.asarray(r, np.float64).ravel().tolist() for r in L.g1_residual],
            "g2_residual": [np.asarray(r, np.float64).ravel().tolist() for r in L.g2_residual],
        })
    L0 = R.layers[0]
    shape = obj.g1_shape_tranform
    out = {
        "descriptor": "dpa2", "rcut": max(rep["rcut"], three["rcut"] if three else 0.0, float(R.rcut)),
        "repinit": rep, "three_body": three, "repinit_dim": dim_rep + dim_3b,
        "g1_shape": _lin(shape) if hasattr(shape, "w") else None,
        "tebd_transform": _lin(obj.tebd_transform) if obj.add_tebd_to_repinit_out else None,
        "tebd": tebd_table[:nt].tolist(), "concat_output_tebd": bool(obj.concat_output_tebd),
        "repformer": {
            "rcut": float(R.rcut), "rcut_smth": float(R.rcut_smth), "sel": [int(R.nnei)], "protection": float(R.env_protection),
            "g1_dim": int(R.g1_dim), "g2_dim": int(R.g2_dim), "axis_neuron": int(R.axis_neuron),
            "activation": an, "act_threshold": aT, "use_sqrt_nnei": bool(R.use_sqrt_nnei),
            "g1_out_conv": bool(L0.g1_out_conv), "g1_out_mlp": bool(L0.g1_out_mlp),
            "update_g1_has_conv": bool(L0.update_g1_has_conv), "update_g1_has_grrg": bool(L0.update_g1_has_grrg),
            "update_g1_has_drrd": bool(L0.update_g1_has_drrd),
            "attn2_hidden": int(L0.attn2_hidden), "attn2_nhead": int(L0.attn2_nhead), "attn2_has_gate": bool(L0.attn2_has_gate),
            "attnw_shift": 20.0,
            "davg": davg[:, 0, :].tolist(), "dstd": dstd[:, 0, :].tolist(), "g2_embd": _lin(R.g2_embd), "layers": layers,
        },
    }
    dim = int(R.g1_dim) + (tebd_table.shape[1] if obj.concat_output_tebd else 0)
    return out, dim


# --------------------------------------------------------------------------- numpy evaluation

def _ln(x, w, b, eps):
    mu, var = x.mean(-1, keepdims=True), x.var(-1, keepdims=True)
    return (x - mu) / np.sqrt(var + eps) * np.asarray(w) + np.asarray(b)


def dpa2_descriptors(m, coord, box, atype):
    import deepmd_export as de
    N = len(coord)
    Lb = np.diag(box)
    R = m["repformer"]
    A = lambda x: act(R["activation"], R["act_threshold"], x)
    nnei, axis, g1d, g2d = R["sel"][0], R["axis_neuron"], R["g1_dim"], R["g2_dim"]
    s_sym = nnei ** -0.5 if R["use_sqrt_nnei"] else 1.0 / nnei
    tebd = np.asarray(m["tebd"])
    # repinit descriptors, initial g1
    g1 = np.zeros((N, g1d))
    nb, SW, H2, G2 = [], [], [], []
    davg, dstd = np.asarray(R["davg"]), np.asarray(R["dstd"])
    for i in range(N):
        d = coord - coord[i]
        d -= Lb * np.round(d / Lb)
        r = np.linalg.norm(d, axis=1)
        ti = atype[i]
        D0 = de._atten_energy(m["repinit"], ti, d, r, atype, m["repinit"]["protection"])
        if m["three_body"] is not None:
            D0 = np.concatenate([D0, de._t_tebd_energy(m["three_body"], ti, d, r, atype, m["three_body"]["protection"])])
        x = _apply(m["g1_shape"], D0) if m["g1_shape"] is not None else D0
        if m["tebd_transform"] is not None:
            x = x + _apply(m["tebd_transform"], tebd[ti])
        g1[i] = A(x)
        idx = [j for j in np.argsort(r, kind="stable") if 0 < r[j] < R["rcut"]][:nnei]
        nb.append(np.array(idx, dtype=int))
        rr, dd = r[idx], d[idx]
        w = _poly_sw(rr, R["rcut_smth"], R["rcut"])
        q = rr + R["protection"]
        Rn = (np.concatenate([(w / q)[:, None], dd * (w / q ** 2)[:, None]], axis=1) - davg[ti]) / dstd[ti]
        SW.append(w); H2.append(Rn[:, 1:])
        G2.append(A(_apply(R["g2_embd"], Rn[:, :1])))
    def grrg(hg):
        return (hg[:, :axis].T @ hg / 3.0).ravel()
    for L in R["layers"]:
        new_g1 = np.zeros_like(g1)
        new_G2 = []
        for i in range(N):
            j, w, h2, g2, gi = nb[i], SW[i], H2[i], G2[i], g1[i]
            n, npad = len(j), nnei - len(j)
            gg1 = g1[j]
            upd, mlp = [], ([] if R["g1_out_mlp"] else [gi])
            if R["g1_out_mlp"]:
                upd.append(A(_apply(L["g1_self_mlp"], gi)))
            if R["update_g1_has_conv"]:
                if R["g1_out_conv"]:
                    upd.append((_apply(L["proj_g1g2"], g2) * gg1 * w[:, None]).sum(0) / nnei)
                else:
                    mlp.append((g2 * _apply(L["proj_g1g2"], gg1) * w[:, None]).sum(0) / nnei)
            if R["update_g1_has_grrg"]:
                mlp.append(grrg(h2.T @ (g2 * w[:, None]) * s_sym))
            if R["update_g1_has_drrd"]:
                mlp.append(grrg(h2.T @ (gg1 * w[:, None]) * s_sym))
            upd.append(A(_apply(L["linear1"], np.concatenate(mlp))))
            new_g1[i] = gi + sum(np.asarray(r) * u for r, u in zip(L["g1_residual"], upd))
            if L["update_chnnl_2"]:
                upd2 = [A(_apply(L["linear2"], g2))]
                if L["update_g2_has_attn"]:
                    nh, hid = R["attn2_nhead"], R["attn2_hidden"]
                    qk = _apply(L["mapqk"], g2).reshape(n, hid, 2 * nh)
                    v = _apply(L["mapv"], g2).reshape(n, g2d, nh)
                    hh = h2 @ h2.T
                    S = w[:, None] * w[None, :]
                    out = np.zeros((n, g2d, nh))
                    for h in range(nh):
                        lg = qk[:, :, h] @ qk[:, :, nh + h].T / np.sqrt(hid)
                        if R["attn2_has_gate"]:
                            lg = lg * hh
                        lg = (lg + 20.0) * S - 20.0
                        mx = lg.max(1, keepdims=True)
                        e = np.exp(lg - mx)
                        Z = e.sum(1, keepdims=True) + npad * np.exp(-20.0 - mx)
                        Aw = e / Z * S * hh / np.sqrt(3.0)
                        out[:, :, h] = Aw @ v[:, :, h]
                    lm = L["attn2_lm"]
                    upd2.append(_ln(_apply(L["head_map"], out.reshape(n, g2d * nh)), lm["w"], lm["b"], lm["eps"]))
                new_G2.append(g2 + sum(np.asarray(r) * u for r, u in zip(L["g2_residual"], upd2)))
            else:
                new_G2.append(g2)
        g1, G2 = new_g1, new_G2
    if m["concat_output_tebd"]:
        g1 = np.concatenate([g1, tebd[atype]], axis=1)
    return g1
