"""Plain numpy model of the DPA4 (SeZM) descriptor, the specification of the exaStamp kernels.

Covers one architecture, the one of deepmd's dpa4 presets: one attention head, one focus, Lebedev SO(3) grid
FFN (GridBranch) and message-node grid product, MLP read-out, env seed FiLM, geometric initial embedding.
The model is read from a deserialized deepmd descriptor object (DescrptDPA4).
Check: python deepmd_dpa4.py   (compares with deepmd's own numpy descriptor on a random water-like frame)
"""
import numpy as np


def sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def silu(x):
    return x * sig(x)


def lin(L, x):
    """NativeLayer without residual: act( x W + b )."""
    y = x @ L.w
    if L.b is not None:
        y = y + L.b
    if L.activation_function == "silu":
        return silu(y)
    assert L.activation_function in ("none", "linear"), L.activation_function
    return y


def norm_vec(x, eps):
    return x / np.sqrt(np.sum(x * x, axis=-1, keepdims=True) + eps * eps)


def smooth_step(x):
    x = np.clip(x, 0.0, 1.0)
    e = float(np.finfo(x.dtype).eps)
    l = np.exp(-1.0 / np.maximum(x, e))
    r = np.exp(-1.0 / np.maximum(1.0 - x, e))
    return np.where(x <= 0.0, 0.0, np.where(x >= 1.0, 1.0, l / (l + r)))


def edge_quaternion(v, r, eps):
    """Rotation sending the edge direction to +z, two charts blended away from the poles."""
    rr = np.sqrt(r * r + eps * eps)[..., None]
    u = v / rr
    x, y, z = u[..., 0], u[..., 1], u[..., 2]
    qp = norm_vec(np.stack([1.0 + z, y, -x, np.zeros_like(x)], -1), eps)
    qn = norm_vec(np.stack([-x, np.zeros_like(x), 1.0 - z, y], -1), eps)
    w = smooth_step(0.5 * (z + 1.0))[..., None]
    qp = np.where(np.sum(qn * qp, -1, keepdims=True) < 0.0, -qp, qp)
    return norm_vec((1.0 - w) * qn + w * qp, eps)


def rotation_matrix(q):
    w, x, y, z = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
                     2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
                     2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)], -1).reshape(q.shape[:-1] + (3, 3))


def wigner(q, lmax, kern, eps):
    """Block diagonal real Wigner D (E, D, D), lmax <= 3."""
    assert lmax <= 3
    C_l2 = kern.C_l2
    q = norm_vec(q, eps)
    E = q.shape[0]
    Dd = (lmax + 1) ** 2
    D = np.zeros((E, Dd, Dd))
    D[:, 0, 0] = 1.0
    if lmax >= 1:
        perm = [1, 2, 0]
        s = np.array([-1.0, -1.0, 1.0])
        D[:, 1:4, 1:4] = rotation_matrix(q)[:, perm][:, :, perm] * np.outer(s, s)
    if lmax >= 2:
        q4 = np.einsum("ei,ej,ek,el->eijkl", q, q, q, q).reshape(E, 256)
        D[:, 4:9, 4:9] = (q4 @ C_l2.reshape(25, 256).T).reshape(E, 5, 5)
    if lmax >= 3:
        mono = np.prod(q[:, None, :] ** np.asarray(kern.exp_l3)[None], -1)
        D[:, 9:16, 9:16] = (mono @ np.asarray(kern.C_l3).T).reshape(E, 7, 7)
    return D


def envelope(r, env):
    u = np.clip((env.rcut - r) / env.rcut, 0.0, 1.0)
    x = 1.0 - u
    c = env._series_coefficients
    s = np.full(x.shape, c[-1])
    for a in reversed(c[:-1]):
        s = a + x * s
    return u ** 4 * s


def rms_scalar(N, x):
    """ScalarRMSNorm, x (.., C)."""
    inv = 1.0 / np.sqrt(np.mean(x * x, -1, keepdims=True) + N.eps)
    return x * inv * N.adam_scale[0]


def so3_lin(L, x):
    """SO3Linear: per degree weights, x (N, D, Cin) -> (N, D, Cout)."""
    W = L.weight.reshape(L.lmax + 1, L.in_channels, L.out_channels)[L.expand_index]
    return np.einsum("ndi,dio->ndo", x, W)


def eq_rms(N, x):
    """EquivariantRMSNorm, x (N, D, C)."""
    x0 = x[:, :1] - np.mean(x[:, :1], -1, keepdims=True)
    bw = N.balance_weight
    mv = np.sum(x0 * x0, axis=(1, 2)) * bw[0] + np.sum((x[:, 1:] * x[:, 1:]) * bw[1:][None, :, None], axis=(1, 2))
    inv = (1.0 / np.sqrt(mv + N.eps))[:, None, None]
    sc = N.adam_scale[N.expand_index][:, 0][None]     # (1, D, C)
    out = np.concatenate([x0, x[:, 1:]], 1) * inv * sc
    out[:, 0] += N.bias.reshape(-1)
    return out


def so2_linear(L, x, f=0):
    """SO2Linear of focus f, x (E, Dm, C): m = 0 block, then m = 1 (neg, pos) with u, v weights."""
    E, Dm, C = x.shape
    nl0 = L.lmax + 1
    out = np.zeros((E, Dm, L.out_channels))
    W0 = L.weight_m0.reshape(nl0 * L.in_channels, -1, nl0 * L.out_channels)[:, f]
    out[:, :nl0] = (x[:, :nl0].reshape(E, -1) @ W0).reshape(E, nl0, L.out_channels)
    off = nl0
    for m, w in enumerate(L.weight_m):
        nl = L.lmax - m
        ib, ob = nl * L.in_channels, nl * L.out_channels
        w = w.reshape(ib, -1, 2 * ob)[:, f]
        wu, wv = w[:, :ob], w[:, ob:]
        neg = x[:, off:off + nl].reshape(E, ib)
        pos = x[:, off + nl:off + 2 * nl].reshape(E, ib)
        out[:, off:off + nl] = (neg @ wu - pos @ wv).reshape(E, nl, L.out_channels)
        out[:, off + nl:off + 2 * nl] = (neg @ wv + pos @ wu).reshape(E, nl, L.out_channels)
        off += 2 * nl
    return out


def gated(A, x, f=0):
    """GatedActivation of focus f, x (E, Dm, C): silu on the scalar row, sigmoid gates per l on the others."""
    y = np.empty_like(x)
    y[:, 0] = silu(x[:, 0])
    g = sig(x[:, 0] @ A.gate_linear.weight.reshape(A.channels, -1, A.lmax * A.channels)[:, f]).reshape(x.shape[0], A.lmax, A.channels)
    y[:, 1:] = x[:, 1:] * g[:, A.expand_index]
    return y


def mixer(M, xl, rad):
    """DynamicRadialDegreeMixer, rank 1 degree_channel: xl (E, Dm, C), rad (E, Dm, C)."""
    E, Dm, C = xl.shape
    k = rad[:, :M.lmax + 1].reshape(E, M.input_dim) @ M.weight
    K = np.zeros((E, Dm * Dm))
    K[:, M.kernel_dense_index] = k[:, M.kernel_compact_index]
    return (K.reshape(E, Dm, Dm) @ xl) * M.channel_basis.reshape(1, 1, C)


def frames(x, K):
    return x.reshape(x.shape[0], x.shape[1], K, x.shape[2] // K)


def to_grid(G, x):
    """x (N, D, K C) -> (N, Ng, C), coefficients flattened (D, K)."""
    N, D, Cw = x.shape
    C = Cw // G.n_frames
    return G.projector.to_grid_mat @ x.reshape(N, D * G.n_frames, C)


def from_grid(G, g):
    N, Ng, C = g.shape
    D = G.projector.coeff_dim // G.n_frames
    return (G.projector.from_grid_mat @ g).reshape(N, D, G.n_frames * C)


def scalar_path(G, coeff, pair, f=0):
    """scalar_pair (N, 2C): gate on all coefficients, SwiGLU added to the scalar of frame zero (gate weights of focus f)."""
    nc = pair.shape[-1] // 2
    so = pair[:, :nc] * sig(pair[:, :nc]) * pair[:, nc:]
    gate = sig(pair @ G.scalar_gate.weight.reshape(2 * nc, -1, nc)[:, f])
    N, D, Cw = coeff.shape
    v = coeff.reshape(N, D, G.n_frames, nc) * gate[:, None, None, :]
    v[:, 0, G.frame_zero_index] += so
    return v.reshape(N, D, Cw)


def message_grid(G, out, x, f=0):
    """Cross GridNet with GridProduct: message-node product, residual scale (query out, context x), one focus stream."""
    W = G.frame_expand.weight[G.frame_expand.degree_index]
    Wc = G.frame_contract.weight[G.frame_contract.degree_index]
    pair = np.concatenate([out[:, 0], x[:, 0]], -1)
    left = np.einsum("ndi,dio->ndo", out, W)
    right = np.einsum("ndi,dio->ndo", x, W)
    c = from_grid(G, to_grid(G, left) * to_grid(G, right))
    c = scalar_path(G, c, pair, f)
    return np.einsum("ndi,dio->ndo", c, Wc) * G.residual_scale[f].reshape(1, 1, -1)


def ffn(F, x):
    """EquivariantFFN with a self GridNet and GridBranch (one branch), x (N, D, C)."""
    h = so3_lin(F.so3_linear_1, x)
    G = F.act
    Op = G.grid_op
    assert Op.n_branches == 1
    Ce = G.expanded_channels
    left, right = h[..., :Ce], h[..., Ce:]
    K = G.n_frames
    pair = np.concatenate([left[:, 0, :Ce // K], right[:, 0, :Ce // K]], -1)      # frame zero first (asserted)
    assert G.frame_zero_index == 0
    N, D, _ = left.shape
    l = (frames(left, K) @ Op.left_proj.weight).reshape(N, D, Ce)
    r = (frames(right, K) @ Op.right_proj.weight).reshape(N, D, Ce)
    c = from_grid(G, to_grid(G, l) * to_grid(G, r))
    c = (frames(c, K) @ Op.out_proj.weight).reshape(N, D, Ce)
    c = scalar_path(G, c, pair)
    return so3_lin(F.so3_linear_2, c)


def readout(F, x):
    """Scalar read-out: output_ffn.call_scalar, GridMLP closed form scalar product, x (N, D, C) -> (N, C)."""
    h = so3_lin(F.so3_linear_1, x)
    G = F.act
    Op = G.grid_op
    Ce = G.expanded_channels
    K = G.n_frames
    assert G.frame_zero_index == 0
    left, right = h[..., :Ce], h[..., Ce:]
    pair = np.concatenate([left[:, 0, :Ce // K], right[:, 0, :Ce // K]], -1)
    fused = np.concatenate([frames(left, K), frames(right, K)], -1)        # (N, D, K, 2C)
    l = fused @ Op.left_proj.weight
    r = fused @ Op.right_proj.weight
    w = G._scalar_product_weight                                           # (D, K)
    s = np.sum(l * w[None, :, :, None] * r, axis=(1, 2))                   # (N, Chid)
    c = s @ Op.out_proj.weight
    nc = pair.shape[-1] // 2
    so = pair[:, :nc] * sig(pair[:, :nc]) * pair[:, nc:]
    gate = sig(pair @ G.scalar_gate.weight.reshape(2 * nc, nc))
    z = c * gate + so                                                      # frame zero of the l = 0 row
    W0 = F.so3_linear_2.weight.reshape(F.so3_linear_2.lmax + 1, F.so3_linear_2.in_channels, F.so3_linear_2.out_channels)[0]
    return z @ W0[:nc]                                                     # frames in order, frame zero first


def segment_softmax(logit, env, dst, nn, z_raw, eps):
    """Softmax per destination, weights env^2, one null logit (softplus(z) + eps)."""
    pos = env > 0.0
    lw = 2.0 * np.log(np.where(pos, env, 1.0))
    eff = np.where(pos, logit + lw, -np.inf)
    # softplus_t of deepmd: log(1 + exp(z))
    null = np.log(np.logaddexp(0.0, float(np.ravel(z_raw)[0])) + eps)
    mx = np.full(nn, null)
    np.maximum.at(mx, dst, eff)
    ex = np.exp(eff - mx[dst])
    den = np.zeros(nn)
    np.add.at(den, dst, ex)
    den += np.exp(null - mx)
    return ex / den[dst]


def focus_alpha(c, src):
    """Cross-focus competition, src (E, F, Cf) -> (E, F)."""
    if c.focus_norm:
        src = src / np.sqrt(np.mean(src * src, -1, keepdims=True) + c.focus_compete_norm.eps) * c.focus_compete_norm.adam_scale[None]
    logit = np.sum(src * c.adamw_focus_compete_w.T[None], 2) / c.focus_softmax_tau
    a = np.exp(logit - np.max(logit, 1, keepdims=True))
    a = a / np.sum(a, 1, keepdims=True)
    return a * (1.0 - c.focus_label_smoothing) + c.focus_label_smoothing / c.n_focus


def so2_conv(c, x, q, rad_block, env, inv_sqrt_deg, src, dst, D):
    """SO2Convolution: x (N, Dd, C) node features, rad_block (E, l+1, C). The hidden width C_wide packs F focus streams."""
    N, Dd, _ = x.shape
    E = src.shape[0]
    F, Cf, Cw = c.n_focus, c.so2_focus_dim, c.hidden_channels
    xp = so3_lin(c.pre_focus_mix, x)                                       # (N, Dd, Cw)
    if c.radial_hidden_proj is not None:
        rad_block = rad_block @ c.radial_hidden_proj.weight
    # edge message
    Dm = D[:, c.coeff_index_m]                                             # (E, Dm, Dd)
    xl = Dm @ xp[src]
    rad = rad_block[:, c.degree_index_m]
    xl = mixer(c.radial_degree_mixer, xl, rad)
    alpha_f = focus_alpha(c, xl[:, 0].reshape(E, F, Cf)) if c.focus_compete and F > 1 else None
    msg_l = np.empty_like(xl)
    for f in range(F):
        xf = xl[:, :, f * Cf:(f + 1) * Cf]
        for k in range(len(c.so2_linears)):
            y = so2_linear(c.so2_linears[k], xf, f)
            if k < len(c.so2_linears) - 1:
                y = gated(c.non_linearities[k], y, f)
            xf = xf + y
        msg_l[:, :, f * Cf:(f + 1) * Cf] = xf if alpha_f is None else xf * alpha_f[:, f, None, None]
    msg = np.transpose(Dm, (0, 2, 1)) @ msg_l * c.rotate_inv_rescale_full[None, :, None]
    # attention and message-node product, per focus
    out_all = np.zeros((N, Dd, Cw))
    for f in range(F):
        sl = slice(f * Cf, (f + 1) * Cf)
        x0 = xp[:, 0, sl]
        l0 = rms_focus(c.attn_qk_norm, x0, f)
        Q = l0 @ c.attn_q_proj.weight.reshape(Cf, F, Cf)[:, f]
        Kk = l0 @ c.attn_k_proj.weight.reshape(Cf, F, Cf)[:, f]
        logit = np.sum(Q[dst] * Kk[src], -1) * Cf ** -0.5 + (rad[:, 0, sl] @ c.adamw_attn_logit_w.reshape(Cf, F)[:, f])
        alpha = segment_softmax(logit, env, dst, N, c.adamw_attn_z_bias_raw.reshape(F)[f], c.eps)
        gate = sig(rms_focus(c.attn_output_gate_norm, x0, f) @ c.adamw_attn_gate_w.reshape(Cf, F)[:, f])
        out = np.zeros((N, Dd, Cf))
        np.add.at(out, dst, msg[:, :, sl] * alpha[:, None, None])
        out = out * gate[:, None, None]
        out_all[:, :, sl] = out + message_grid(c.message_node_grid_product, out, xp[:, :, sl], f)
    return so3_lin(c.post_focus_mix, out_all)


def rms_focus(N, x, f):
    """ScalarRMSNorm of focus f, x (.., Cf)."""
    return x / np.sqrt(np.mean(x * x, -1, keepdims=True) + N.eps) * N.adam_scale[f]


def dpa4_descriptor(m, coord_ext, atype_ext, nlist, mapping):
    """Descriptor (nloc, C) of one frame; ghosts are images (mapping), as a single domain."""
    nloc = nlist.shape[0]
    C = m.channels
    eps = m.eps
    # edges: j (source, image of a local atom) -> i (destination)
    ii, kk = np.nonzero(nlist >= 0)
    jj = nlist[ii, kk]
    v = coord_ext[jj] - coord_ext[ii]
    dst = ii
    src = mapping[jj]
    E = len(dst)
    r = np.sqrt(np.sum(v * v, -1) + eps * eps)
    atype = atype_ext[:nloc]
    te = m.type_embedding.adam_type_embedding[atype]
    env = envelope(r, m.edge_envelope)
    # radial basis
    f = m.radial_basis.adam_freqs[0]
    arg = r[:, None] * f
    rbf = f * np.sin(arg) / arg * envelope(r, m.radial_basis.envelope)[:, None]
    # Wigner
    q = edge_quaternion(v, r, eps)
    D = wigner(q, m.lmax, m.wigner_calc.small_order_kernels, m.wigner_calc.eps)
    # degree
    deg = np.zeros(nloc)
    np.add.at(deg, dst, env * env)
    isd = 1.0 / np.sqrt(deg + m.deg_norm_floor)
    # radial features
    net = m.radial_embedding.net
    h = lin(net[0], rbf)
    if m.radial_norm:
        h = rms_scalar_plain(net[1], h)
    rf = (silu(h) @ net[-1].w).reshape(E, m.node_init_lmax + 1, C) * env[:, None, None]
    # env seed
    S = m.env_seed_embedding
    r2 = np.sum(v * v, -1)
    inv_r = 1.0 / np.sqrt(r2 + S.eps * S.eps)
    s = env * inv_r
    rt = np.concatenate([s[:, None], s[:, None] * (v * inv_r[:, None])], -1)
    tt = S.env_type_embed.adam_type_embedding
    g = lin(S.g_layer2, lin(S.g_layer1, np.concatenate([lin(S.rbf_proj_layer2, lin(S.rbf_proj_layer1, rbf)), tt[atype[src]], tt[atype[dst]]], -1)))
    agg = np.zeros((nloc, S.coord_dim, S.embed_dim))
    np.add.at(agg, dst, rt[:, :, None] * g[:, None, :])
    agg = agg * isd[:, None, None]
    Dd = np.transpose(agg, (0, 2, 1)) @ agg[:, :, :S.axis_dim]
    film = lin(S.output_proj, Dd.reshape(nloc, S.embed_dim * S.axis_dim))
    sh = rms_scalar(m.film_scale_norm, film[:, :C]) if m.film_norm else film[:, :C]
    sf = rms_scalar(m.film_shift_norm, film[:, C:]) if m.film_norm else film[:, C:]
    scale = 1.0 + np.exp(m.film_scale_strength_log[0]) * np.tanh(sh)
    shift = np.exp(m.film_shift_strength_log[0]) * np.tanh(sf)
    x = np.zeros((nloc, m.node_init_dim, C))
    x[:, 0] = te * scale + shift
    # geometric initial embedding
    G = m.gie
    zonal = D[:, G.zonal_m0_col_index_for_row, G.non_scalar_row_index]     # Dt[row, col] = D[col, row]
    msg = zonal[:, :, None] * rf[:, 1:][:, G.radial_slot_index_for_row]
    gie = np.zeros((nloc, m.node_init_dim - 1, C))
    np.add.at(gie, dst, msg)
    x[:, 1:] += gie * isd[:, None, None]
    # radial features fused with the edge types
    rf = rf + (te[src] + te[dst])[:, None, :]
    # blocks
    for k, b in enumerate(m.blocks):
        so = so2_conv(b.so2_conv, x, None, rf[:, :m.rad_sizes_per_block[k]], env, isd, src, dst, D)
        x = x + eq_rms(b.post_so2_norm, so)
        for i in range(b.ffn_blocks):
            x = x + ffn(b.ffns[i], eq_rms(b.pre_ffn_norms[i], x))
    return x[:, 0] + readout(m.output_ffn, x)


def rms_scalar_plain(N, x):
    inv = 1.0 / np.sqrt(np.mean(x * x, -1, keepdims=True) + N.eps)
    return x * inv * N.adam_scale


def fail(msg):
    raise SystemExit(f"deepmd_export (dpa4): {msg}")


def _f(a):
    return np.asarray(a, dtype=np.float64).ravel().tolist()


def _layer(L, acts):
    """Weights of a NativeLayer without bias, residual or time step; returns (weights, is_silu)."""
    if L.b is not None or L.idt is not None or L.resnet:
        fail("layer with bias, idt or residual is not supported")
    if L.activation_function not in acts:
        fail(f"activation '{L.activation_function}' is not supported here")
    return _f(L.w), L.activation_function == "silu"


def _eq_norm(N, lmax):
    return {"scale": _f(N.adam_scale), "bias": _f(N.bias), "bw": _f(N.balance_weight), "eps": float(N.eps)}


def _check(D):
    if D.exclude_types:
        fail("excluded types are not supported")
    if D.lmax not in (1, 2, 3):
        fail(f"lmax {D.lmax} is not supported (1, 2 or 3)")
    if D.n_atten_head != 1 or D.focus_dim != 0:
        fail("n_atten_head != 1 or focus_dim are not supported")
    if any(b.so2_conv.focus_compete_bias is not None for b in D.blocks):
        fail("focus competition bias is not supported")
    if D.charge_spin_embedding is not None or D.spin_embedding is not None or D.bridging_switch is not None or D.bridging_clamp is not None:
        fail("charge, spin embeddings and zone bridging are not supported")
    if not D.use_env_seed or not D.use_gie or D.use_full_attn_res or D.use_block_attn_res:
        fail("this architecture is not supported (env seed, GIE required, no attention residuals)")
    if D.so3_readout != "mlp" or D.readout_layers != 1 or D.readout_pre_layers:
        fail("read-out other than one mlp layer is not supported")
    if D.node_init_lmax != D.lmax or any(l != D.lmax for l in D.l_schedule) or any(m != 1 for m in D.m_schedule):
        fail("degree schedules are not supported")
    if D.radial_basis.basis_family != "bessel":
        fail("radial basis other than bessel is not supported")
    for b in D.blocks:
        c = b.so2_conv
        if b.ffn_blocks != 1 or type(b.pre_so2_norm).__name__ != "Identity" or type(b.post_ffn_norms[0]).__name__ != "Identity" or b.layer_scale:
            fail("block norms other than the presets' are not supported")
        if c.edge_cartesian or c.node_cartesian_tp is not None or c.node_wise_grid_product is not None or c.use_so2_attn_res or c.layer_scale \
                or c.attn_v_proj is not None or c.attn_o_proj is not None or c.attn_focus_mix is not None:
            fail("SO(2) convolution options other than the presets' are not supported")
        if c.message_node_grid_product is None or type(c.message_node_grid_product.grid_op).__name__ != "GridProduct":
            fail("a message-node grid product is required")
        M = c.radial_degree_mixer
        if M is not None and (M.mode != "degree_channel" or M.rank != 1):
            fail("radial mixer other than degree_channel rank 1 is not supported")
        for L in c.so2_linears:
            if L.bias0 is not None or len(L.weight_m) != 1:
                fail("SO(2) linear with bias or mmax != 1 is not supported")
        if not (c.non_linearities[-1].__class__.__name__ == "Identity" and all(n.activation_function == "silu" for n in c.non_linearities[:-1])):
            fail("SO(2) activations other than the presets' are not supported")
        G = b.ffns[0].act
        if type(G).__name__ != "SO3GridNet" or type(G.grid_op).__name__ != "GridBranch" or G.grid_op.n_branches != 1 or G.frame_zero_index != 0:
            fail("FFN grid other than one-branch SO(3) GridBranch is not supported")
    G = D.output_ffn.act
    if type(G).__name__ != "SO3GridNet" or type(G.grid_op).__name__ != "GridMLP" or G.grid_op.mode != "self" or G.frame_zero_index != 0:
        fail("read-out grid other than SO(3) GridMLP is not supported")


def export_dpa4(m):
    """Native YAML of a serialized DPA4 model (dict): descriptor, GLU fitting, bias."""
    import itertools
    from deepmd.dpmodel.model.base_model import BaseModel
    M = BaseModel.deserialize(m)
    am = M.atomic_model
    D, F = am.descriptor, am.fitting_net
    _check(D)
    if getattr(F, "numb_fparam", 0) or getattr(F, "numb_aparam", 0) or not F.mixed_types or len(F.neuron) != 1:
        fail("fitting with fparam, aparam, per-type networks or several layers is not supported")
    net = F.nets[()]
    lmax = D.lmax
    # degree 4 monomials of the quaternion
    Cl2 = np.zeros((25, 35))
    mono = list(itertools.combinations_with_replacement(range(4), 4))
    if lmax >= 2:
        full = np.asarray(D.wigner_calc.small_order_kernels.C_l2).reshape(25, 4, 4, 4, 4)
        for k, t in enumerate(mono):
            for pm in set(itertools.permutations(t)):
                Cl2[:, k] += full[(slice(None),) + pm]
    c0 = D.blocks[0].so2_conv
    Dd = (lmax + 1) ** 2
    deg = [int(np.floor(np.sqrt(i))) for i in range(Dd)]
    base = {0: 0, 1: 1, 2: 10, 3: 35}

    def packed(row, col):
        l = deg[row]
        return base[l] + (row - l * l) * (2 * l + 1) + (col - l * l)
    coeff = [int(x) for x in c0.coeff_index_m]
    rot = [(r, d, packed(coeff[r], d)) for r in range(len(coeff)) for d in range(deg[coeff[r]] ** 2, (deg[coeff[r]] + 1) ** 2)]
    if [int(x) for x in c0.degree_index_m[:lmax + 1]] != list(range(lmax + 1)):
        fail("unexpected coefficient order")
    G0 = c0.message_node_grid_product
    T, Fm = np.asarray(G0.projector.to_grid_mat), np.asarray(G0.projector.from_grid_mat)
    blocks = []
    for b in D.blocks:
        c = b.so2_conv
        G = c.message_node_grid_product
        Fg = b.ffns[0].act
        for g in (G, Fg):
            if not (np.array_equal(g.projector.to_grid_mat, T) and np.array_equal(g.projector.from_grid_mat, Fm)):
                fail("grids differ between the nets")
        Mx = c.radial_degree_mixer
        blocks.append({
            "pre_mix": _f(c.pre_focus_mix.weight), "post_mix": _f(c.post_focus_mix.weight),
            "lin0": [_f(L.weight_m0) for L in c.so2_linears], "linm": [_f(L.weight_m[0]) for L in c.so2_linears],
            "gate": [_f(n.gate_linear.weight) for n in c.non_linearities[:-1]],
            "mixw": _f(Mx.weight) if Mx is not None else [], "mixb": _f(Mx.channel_basis) if Mx is not None else [],
            "qk_scale": _f(c.attn_qk_norm.adam_scale), "qk_eps": float(c.attn_qk_norm.eps),
            "gate_scale": _f(c.attn_output_gate_norm.adam_scale), "gate_eps": float(c.attn_output_gate_norm.eps),
            "q": _f(c.attn_q_proj.weight), "k": _f(c.attn_k_proj.weight), "logit_w": _f(c.adamw_attn_logit_w), "gate_w": _f(c.adamw_attn_gate_w),
            "log_null": [float(np.log(np.logaddexp(0.0, float(z)) + c.eps)) for z in np.ravel(c.adamw_attn_z_bias_raw)],
            **({"rproj": _f(c.radial_hidden_proj.weight)} if c.radial_hidden_proj is not None else {}),
            **({"comp_w": _f(c.adamw_focus_compete_w)} if c.focus_compete and c.n_focus > 1 else {}),
            **({"comp_scale": _f(c.focus_compete_norm.adam_scale), "comp_eps": float(c.focus_compete_norm.eps)} if c.focus_compete and c.n_focus > 1 and c.focus_norm else {}),
            "gexp": _f(G.frame_expand.weight), "gcon": _f(G.frame_contract.weight), "ggate": _f(G.scalar_gate.weight), "gres": _f(G.residual_scale),
            "post_norm": _eq_norm(b.post_so2_norm, lmax), "ffn_norm": _eq_norm(b.pre_ffn_norms[0], lmax),
            "flin1": _f(b.ffns[0].so3_linear_1.weight), "flin2": _f(b.ffns[0].so3_linear_2.weight),
            "fl": _f(Fg.grid_op.left_proj.weight), "fr": _f(Fg.grid_op.right_proj.weight), "fo": _f(Fg.grid_op.out_proj.weight), "fgate": _f(Fg.scalar_gate.weight),
        })
    Mx = c0.radial_degree_mixer
    S = D.env_seed_embedding
    net_r = D.radial_embedding.net
    if len(net_r) != (4 if D.radial_norm else 3) or (D.radial_norm and type(net_r[1]).__name__ != "RMSNorm"):
        fail("unexpected radial network")
    Ro = D.output_ffn
    RG = Ro.act
    Ch = RG.channels
    W2 = np.asarray(Ro.so3_linear_2.weight).reshape(lmax + 1, RG.expanded_channels, D.channels)[0][:Ch]
    r1, a1 = _layer(S.rbf_proj_layer1, ("silu", "none", "linear"))
    r2, a2 = _layer(S.rbf_proj_layer2, ("silu", "none", "linear"))
    g1, a3 = _layer(S.g_layer1, ("silu", "none", "linear"))
    g2, a4 = _layer(S.g_layer2, ("silu", "none", "linear"))
    so, a5 = _layer(S.output_proj, ("silu", "none", "linear"))
    desc = {
        "n_focus": int(c0.n_focus), "focus_tau": float(c0.focus_softmax_tau), "focus_ls": float(c0.focus_label_smoothing),
        "descriptor": "dpa4", "rcut": float(D.rcut), "sel": [int(D.nnei)], "channels": int(D.channels), "lmax": lmax,
        "eps": float(D.eps), "wigner_eps": float(D.wigner_calc.eps), "deg_floor": float(D.deg_norm_floor if D.version >= 1.1 else D.eps),
        "env_series": list(map(float, D.edge_envelope._series_coefficients)), "rbf_series": list(map(float, D.radial_basis.envelope._series_coefficients)),
        "freq": _f(D.radial_basis.adam_freqs),
        "te": _f(D.type_embedding.adam_type_embedding), "tt": _f(S.env_type_embed.adam_type_embedding), "tdim": int(S.type_dim),
        "rf1": _f(net_r[0].w), "rf2": _f(net_r[-1].w), "hrf": int(np.asarray(net_r[0].w).shape[1]),
        **({"rf_scale": _f(net_r[1].adam_scale), "rf_eps": float(net_r[1].eps)} if D.radial_norm else {}),
        "seed_eps": float(S.eps), "emb": int(S.embed_dim), "axis": int(S.axis_dim), "s_r": int(S.rbf_out_dim), "s_gh": int(S.hidden_dim),
        "s_r1": r1, "s_r2": r2, "s_g1": g1, "s_g2": g2, "s_out": so, "s_act": [a1, a2, a3, a4, a5],
        "film_scale": float(np.exp(np.ravel(D.film_scale_strength_log)[0])), "film_shift": float(np.exp(np.ravel(D.film_shift_strength_log)[0])),
        **({"film_s": _f(D.film_scale_norm.adam_scale), "film_s_eps": float(D.film_scale_norm.eps), "film_f": _f(D.film_shift_norm.adam_scale), "film_f_eps": float(D.film_shift_norm.eps)} if D.film_norm else {}),
        "C_l2": _f(Cl2), "mono": [int(x) for t in mono for x in t],
        **({"C_l3": _f(D.wigner_calc.small_order_kernels.C_l3), "exp_l3": [int(x) for x in np.ravel(D.wigner_calc.small_order_kernels.exp_l3)]} if lmax >= 3 else {}),
        "to_grid": _f(T), "from_grid": _f(Fm), "ng": int(T.shape[0]), "frames": int(G0.n_frames), "rescale": _f(c0.rotate_inv_rescale_full), "sweight": _f(RG._scalar_product_weight),
        "coeff": coeff, "deg_m": [int(x) for x in c0.degree_index_m], "gate_group": [0] + [int(x) for x in c0.non_linearities[0].expand_index],
        "mix_compact": [int(x) for x in Mx.kernel_compact_index] if Mx is not None else [],
        "mix_r": [int(x) // len(coeff) for x in Mx.kernel_dense_index] if Mx is not None else [],
        "mix_s": [int(x) % len(coeff) for x in Mx.kernel_dense_index] if Mx is not None else [],
        "nk": int(Mx.degree_kernel_size) if Mx is not None else 0,
        "z_offset": [packed(int(c), int(r)) for c, r in zip(D.gie.zonal_m0_col_index_for_row, D.gie.non_scalar_row_index)],
        "z_slot": [int(x) + 1 for x in D.gie.radial_slot_index_for_row],
        "rot_r": [t[0] for t in rot], "rot_d": [t[1] for t in rot], "rot_o": [t[2] for t in rot],
        "blocks": blocks,
        "ffn_channels": int(RG.channels), "rhid": int(RG.grid_op.hidden_channels),
        "rlin1": _f(Ro.so3_linear_1.weight), "rl": _f(RG.grid_op.left_proj.weight), "rr": _f(RG.grid_op.right_proj.weight), "ro": _f(RG.grid_op.out_proj.weight), "rgate": _f(RG.scalar_gate.weight), "rlin2": _f(W2),
        "fit_h": int(net.neuron[0]), "fw": _f(net.hidden_layers[0].w), "fb": _f(net.hidden_layers[0].b), "fout": _f(net.output_layer.w),
        "fout_b": float(np.ravel(net.output_layer.b)[0]) if net.output_layer.b is not None else 0.0,
    }
    bias = np.ravel(np.asarray(am.out_bias, dtype=np.float64)[0]) + np.ravel(np.asarray(F.bias_atom_e, dtype=np.float64))
    return desc, bias


def fit_net(net, d):
    """GLU fitting network: y = x W + b, split value / gate, value * silu(gate), linear output."""
    x = d
    for nh, L in zip(net.neuron, net.hidden_layers):
        y = x @ L.w + L.b
        x = y[..., :nh] * silu(y[..., nh:])
    y = x @ net.output_layer.w
    if net.output_layer.b is not None:
        y = y + net.output_layer.b
    return y[..., 0]


def zero_arrays(m):
    seen = set()
    def walk(o, name):
        if id(o) in seen or not hasattr(o, "__dict__") or isinstance(o, type):
            return
        seen.add(id(o))
        for k, v in vars(o).items():
            if isinstance(v, np.ndarray) and v.dtype.kind == "f" and v.size > 0 and np.all(v == 0.0):
                yield name + "." + k
            elif isinstance(v, (list, tuple)):
                for i, x in enumerate(v):
                    yield from walk(x, f"{name}.{k}[{i}]")
            elif isinstance(v, dict):
                for kk, x in v.items():
                    yield from walk(x, f"{name}.{k}[{kk}]")
            elif hasattr(v, "__dict__"):
                yield from walk(v, name + "." + k)
    return list(walk(m, "m"))


def randomize(m, seed=7):
    """Fresh DPA4 models zero-initialize their output layers: fill zero weights and residual scales with random
    values, so that every path of the network reaches the output (test models only)."""
    rng = np.random.default_rng(seed)
    seen, nrand = set(), [0]

    def walk(o):
        if id(o) in seen or not hasattr(o, "__dict__") or isinstance(o, type):
            return
        seen.add(id(o))
        for k, v in vars(o).items():
            if isinstance(v, np.ndarray) and v.dtype.kind == "f" and v.size > 0:
                if k == "residual_scale" or ( np.all(v == 0.0) and ( k in ("weight", "w", "weight_m0", "bias") or k.startswith("adam") ) ):
                    setattr(o, k, rng.normal(0.0, 0.2, v.shape))
                    nrand[0] += 1
            elif isinstance(v, (list, tuple)):
                for i, x in enumerate(v):
                    if isinstance(x, np.ndarray) and x.dtype.kind == "f" and np.all(x == 0.0):
                        v[i] = rng.normal(0.0, 0.2, x.shape)
                        nrand[0] += 1
                    else:
                        walk(x)
            elif isinstance(v, dict):
                for x in v.values():
                    walk(x)
            elif hasattr(v, "__dict__") and not callable(getattr(v, "__func__", None)):
                walk(v)
    walk(m)
    return nrand[0]


if __name__ == "__main__":
    from deepmd.dpmodel.descriptor.dpa4 import DescrptDPA4
    from deepmd.dpmodel.utils.nlist import extend_input_and_build_neighbor_list
    import sys
    args = dict(rcut=6.0, n_radial=16, use_env_seed=True, mmax=1, radial_so2_mode="degree_channel", focus_dim=0, n_atten_head=1,
                message_node_so3=True, ffn_neurons=0, ffn_so3_grid=True, grid_mlp=False, grid_branch=[0, 0, 1], ffn_blocks=1,
                so3_readout="mlp", precision="float64", channels=32, lmax=2, n_blocks=2, mixing_layers=3, radial_so2_rank=1, n_focus=1, seed=1)
    m = DescrptDPA4(ntypes=2, sel=120, type_map=["O", "H"], **args)
    randomize(m)
    z = np.load(sys.argv[1]) if len(sys.argv) > 1 else None
    if z is not None:
        coord, box, atype = z["coord"].reshape(1, -1, 3), z["box"].reshape(1, 9), z["atype"].reshape(1, -1)
    else:
        rng = np.random.default_rng(3)
        n = 48
        coord = rng.uniform(0.0, 10.0, (1, n, 3))
        box = (np.eye(3) * 10.0).reshape(1, 9)
        atype = (np.arange(n) % 3 > 0).astype(int).reshape(1, -1)
    ce, ae, map_, nl = None, None, None, None
    ce, ae, map_, nl = extend_input_and_build_neighbor_list(coord, atype, m.get_rcut(), m.get_sel(), mixed_types=True, box=box)
    ref = m.call(ce, ae, nl, mapping=map_)[0][0]
    got = dpa4_descriptor(m, ce[0].reshape(-1, 3), ae[0], nl[0], map_[0])
    print("descriptor shape", got.shape, "max|ref|", np.abs(ref).max(), "max|diff|", np.abs(got - ref).max())
    t0 = atype[0] == 0
    print("spread between atoms of one type", np.abs(ref[t0] - ref[t0].mean(0)).max(), "(0 would mean an empty test)")
    print("zero weight arrays left:", sum(1 for x in zero_arrays(m)))
