"""Random-weight DeePMD models of several descriptors, with non-trivial statistics, for the native-kernel tests.

    python deepmd_test_models.py se_e2_r se_e3     # writes se_e2_r.dp, se_e3.dp (convert with dp convert-backend)
"""
import sys, numpy as np
from deepmd.dpmodel.model.model import get_model
from deepmd.dpmodel.utils.serialization import save_dp_model

DESCRIPTORS = {
    "se_e2_r": {"type": "se_e2_r", "sel": [46, 92], "rcut": 6.0, "rcut_smth": 0.5, "neuron": [5, 10, 20]},
    "se_e3": {"type": "se_e3", "sel": [20, 40], "rcut": 4.5, "rcut_smth": 0.5, "neuron": [4, 8, 16]},
    "se_e2_a_tt": {"type": "se_e2_a", "sel": [46, 92], "rcut": 6.0, "rcut_smth": 0.5, "neuron": [6, 12, 24], "axis_neuron": 4, "type_one_side": False},
    "se_e2_a": {"type": "se_e2_a", "sel": [46, 92], "rcut": 6.0, "rcut_smth": 0.5, "neuron": [6, 12, 24], "axis_neuron": 4},
    # DPA1 / se_atten: mixed neighbor list, gated self-attention between the neighbors of an atom
    "dpa1": {"type": "dpa1", "sel": 120, "rcut": 6.0, "rcut_smth": 0.5, "neuron": [8, 16, 32], "axis_neuron": 4, "tebd_dim": 8,
             "attn": 16, "attn_layer": 2, "attn_dotr": True, "type_one_side": False, "tebd_input_mode": "concat"},
    "dpa1_strip": {"type": "dpa1", "sel": 120, "rcut": 6.0, "rcut_smth": 0.5, "neuron": [8, 16, 32], "axis_neuron": 4, "tebd_dim": 6,
                   "attn": 24, "attn_layer": 1, "attn_dotr": False, "normalize": False, "type_one_side": True, "tebd_input_mode": "strip",
                   "concat_output_tebd": False},
    "dpa1_nosmooth": {"type": "dpa1", "sel": 120, "rcut": 6.0, "rcut_smth": 0.5, "neuron": [8, 16, 32], "axis_neuron": 4, "tebd_dim": 8,
                      "attn": 16, "attn_layer": 2, "attn_dotr": True, "smooth_type_embedding": False, "temperature": 0.7},
    "dpa1_sel40": {"type": "dpa1", "sel": 40, "rcut": 6.0, "rcut_smth": 0.5, "neuron": [8, 16, 32], "axis_neuron": 4, "tebd_dim": 8,
                   "attn": 16, "attn_layer": 2, "attn_dotr": True, "type_one_side": True},
    "se_t_tebd": {"type": "se_e3_tebd", "sel": 40, "rcut": 4.5, "rcut_smth": 0.5, "neuron": [4, 8, 16], "tebd_dim": 6},
    "se_t_tebd_strip": {"type": "se_e3_tebd", "sel": 40, "rcut": 4.5, "rcut_smth": 0.5, "neuron": [4, 8, 16], "tebd_dim": 6,
                        "tebd_input_mode": "strip", "concat_output_tebd": False},
    "hybrid": {"type": "hybrid", "list": [
        {"type": "se_e2_a", "sel": [20, 40], "rcut": 4.0, "rcut_smth": 0.5, "neuron": [4, 8, 12], "axis_neuron": 4},
        {"type": "dpa1", "sel": 100, "rcut": 6.0, "rcut_smth": 0.5, "neuron": [6, 12, 16], "axis_neuron": 4, "tebd_dim": 4,
         "attn": 8, "attn_layer": 1, "attn_dotr": True},
        {"type": "se_e3_tebd", "sel": 30, "rcut": 3.5, "rcut_smth": 0.5, "neuron": [4, 8], "tebd_dim": 4, "concat_output_tebd": False}]},
    # DPA3 (repflows): the options of deepmd's examples/water/dpa3 (small widths), float64
    "dpa3": {"type": "dpa3", "activation_function": "silut:10.0", "precision": "float64", "concat_output_tebd": False,
             "repflow": {"n_dim": 16, "e_dim": 8, "a_dim": 8, "nlayers": 3, "e_rcut": 6.0, "e_rcut_smth": 5.3, "e_sel": 120,
                         "a_rcut": 4.0, "a_rcut_smth": 3.5, "a_sel": 30, "axis_neuron": 4, "fix_stat_std": 0.3,
                         "a_compress_rate": 1, "a_compress_e_rate": 2, "a_compress_use_split": True, "update_angle": True,
                         "smooth_edge_update": True, "edge_init_use_dist": True, "use_exp_switch": True,
                         "update_style": "res_residual", "update_residual": 0.1, "update_residual_init": "const"}},
    # DPA2: repinit (se_atten, no attention) + three-body (se_t_tebd), repformer layers (message passing)
    "dpa2": {"type": "dpa2", "precision": "float64", "add_tebd_to_repinit_out": False,
             "repinit": {"tebd_dim": 8, "rcut": 6.0, "rcut_smth": 0.5, "nsel": 120, "neuron": [8, 16, 24], "axis_neuron": 4,
                         "activation_function": "tanh", "three_body_sel": 48, "three_body_rcut": 4.0, "three_body_rcut_smth": 3.5,
                         "use_three_body": True, "three_body_neuron": [4, 8]},
             "repformer": {"rcut": 4.0, "rcut_smth": 3.5, "nsel": 48, "nlayers": 3, "g1_dim": 16, "g2_dim": 8, "attn2_hidden": 8,
                           "attn2_nhead": 2, "attn1_hidden": 16, "attn1_nhead": 2, "axis_neuron": 4, "update_h2": False,
                           "update_g1_has_conv": True, "update_g1_has_grrg": True, "update_g1_has_drrd": True, "update_g1_has_attn": False,
                           "update_g2_has_g1g1": False, "update_g2_has_attn": True, "update_style": "res_residual", "update_residual": 0.01,
                           "update_residual_init": "norm", "attn2_has_gate": True, "use_sqrt_nnei": True, "g1_out_conv": True, "g1_out_mlp": True}},
    # deepmd's water example (input_torch_medium.json), for timings
    "dpa2_full": {"type": "dpa2", "precision": "float64", "add_tebd_to_repinit_out": False,
                  "repinit": {"tebd_dim": 8, "rcut": 6.0, "rcut_smth": 0.5, "nsel": 120, "neuron": [25, 50, 100], "axis_neuron": 12,
                              "activation_function": "tanh", "three_body_sel": 48, "three_body_rcut": 4.0, "three_body_rcut_smth": 3.5,
                              "use_three_body": True},
                  "repformer": {"rcut": 4.0, "rcut_smth": 3.5, "nsel": 48, "nlayers": 6, "g1_dim": 128, "g2_dim": 32, "attn2_hidden": 32,
                                "attn2_nhead": 4, "attn1_hidden": 128, "attn1_nhead": 4, "axis_neuron": 4, "update_h2": False,
                                "update_g1_has_conv": True, "update_g1_has_grrg": True, "update_g1_has_drrd": True, "update_g1_has_attn": False,
                                "update_g2_has_g1g1": False, "update_g2_has_attn": True, "update_style": "res_residual", "update_residual": 0.01,
                                "update_residual_init": "norm", "attn2_has_gate": True, "use_sqrt_nnei": True, "g1_out_conv": True, "g1_out_mlp": True}},
    "dpa2_b": {"type": "dpa2", "precision": "float64", "add_tebd_to_repinit_out": True, "concat_output_tebd": True,
               "repinit": {"tebd_dim": 6, "rcut": 6.0, "rcut_smth": 0.5, "nsel": 120, "neuron": [8, 16], "axis_neuron": 3,
                           "activation_function": "tanh", "use_three_body": False, "tebd_input_mode": "strip"},
               "repformer": {"rcut": 4.0, "rcut_smth": 3.5, "nsel": 48, "nlayers": 2, "g1_dim": 12, "g2_dim": 6, "attn2_hidden": 6,
                             "attn2_nhead": 2, "axis_neuron": 3, "update_h2": False, "activation_function": "silu",
                             "update_g1_has_conv": True, "update_g1_has_grrg": True, "update_g1_has_drrd": True, "update_g1_has_attn": False,
                             "update_g2_has_g1g1": False, "update_g2_has_attn": False, "update_style": "res_residual", "update_residual": 0.01,
                             "update_residual_init": "const", "attn2_has_gate": True, "use_sqrt_nnei": True, "g1_out_conv": False, "g1_out_mlp": False}},
    # DPA2 variant: no gate, 1/nnei normalization, no conv / drrd, three-body, non-zero repformer radial averages
    "dpa2_c": {"type": "dpa2", "precision": "float64", "add_tebd_to_repinit_out": True, "concat_output_tebd": False,
               "repinit": {"tebd_dim": 5, "rcut": 5.0, "rcut_smth": 0.5, "nsel": 100, "neuron": [6, 12], "axis_neuron": 2,
                           "activation_function": "tanh", "three_body_sel": 40, "three_body_rcut": 3.5, "three_body_rcut_smth": 3.0,
                           "use_three_body": True, "three_body_neuron": [4, 6]},
               "repformer": {"rcut": 4.5, "rcut_smth": 3.0, "nsel": 60, "nlayers": 3, "g1_dim": 10, "g2_dim": 6, "attn2_hidden": 5,
                             "attn2_nhead": 3, "axis_neuron": 2, "update_h2": False, "activation_function": "gelu",
                             "update_g1_has_conv": False, "update_g1_has_grrg": True, "update_g1_has_drrd": False, "update_g1_has_attn": False,
                             "update_g2_has_g1g1": False, "update_g2_has_attn": True, "update_style": "res_residual", "update_residual": 0.01,
                             "update_residual_init": "norm", "attn2_has_gate": False, "use_sqrt_nnei": False, "g1_out_conv": True, "g1_out_mlp": True}},
    "dpa3_full": {"type": "dpa3", "activation_function": "silut:10.0", "precision": "float64", "concat_output_tebd": False,
                  "repflow": {"n_dim": 128, "e_dim": 64, "a_dim": 32, "nlayers": 6, "e_rcut": 6.0, "e_rcut_smth": 5.3, "e_sel": 120,
                              "a_rcut": 4.0, "a_rcut_smth": 3.5, "a_sel": 30, "axis_neuron": 4, "fix_stat_std": 0.3,
                              "a_compress_rate": 1, "a_compress_e_rate": 2, "a_compress_use_split": True, "update_angle": True,
                              "smooth_edge_update": True, "edge_init_use_dist": True, "use_exp_switch": True,
                              "update_style": "res_residual", "update_residual": 0.1, "update_residual_init": "const"}},
    "dpa3_b": {"type": "dpa3", "activation_function": "silu", "precision": "float64", "concat_output_tebd": True,
               "repflow": {"n_dim": 12, "e_dim": 8, "a_dim": 6, "nlayers": 2, "e_rcut": 6.0, "e_rcut_smth": 0.5, "e_sel": 120,
                           "a_rcut": 4.0, "a_rcut_smth": 0.5, "a_sel": 30, "axis_neuron": 3, "a_compress_rate": 0, "update_angle": True,
                           "smooth_edge_update": False, "edge_init_use_dist": False, "use_exp_switch": False,
                           "update_style": "res_residual", "update_residual": 0.1, "update_residual_init": "const"}},
    "dpa3_c": {"type": "dpa3", "activation_function": "silut:3.0", "precision": "float64", "concat_output_tebd": False,
               "repflow": {"n_dim": 12, "e_dim": 8, "a_dim": 8, "nlayers": 2, "e_rcut": 6.0, "e_rcut_smth": 5.3, "e_sel": 1200,
                           "a_rcut": 4.0, "a_rcut_smth": 3.5, "a_sel": 300, "axis_neuron": 4, "a_compress_rate": 2, "a_compress_e_rate": 1,
                           "a_compress_use_split": False, "update_angle": True, "smooth_edge_update": True, "edge_init_use_dist": True,
                           "use_exp_switch": True, "use_dynamic_sel": True, "sel_reduce_factor": 10.0,
                           "update_style": "res_residual", "update_residual": 0.1, "update_residual_init": "const"}},
    "se_atten_v2": {"type": "se_atten_v2", "sel": 120, "rcut": 6.0, "rcut_smth": 0.5, "neuron": [8, 16, 32], "axis_neuron": 4, "tebd_dim": 8,
                    "attn": 16, "attn_layer": 2, "attn_dotr": True},
}
# descriptors with a single mixed neighbor list: zero averages (DPA1 is trained with set_davg_zero), random radial
# averages only where padded neighbor slots stay type independent (strip mode, smooth)
MIXED = {"se_t_tebd": True, "se_t_tebd_strip": True, "dpa1_sel40": False, "dpa1": False, "dpa1_strip": True, "dpa1_nosmooth": False, "se_atten_v2": True}

def _dpa4_base(**kw):
    d = dict(type="dpa4", sel=120, rcut=6.0, n_radial=16, use_env_seed=True, mmax=1, radial_so2_mode="degree_channel", focus_dim=0, n_atten_head=1,
             message_node_so3=True, ffn_neurons=0, ffn_so3_grid=True, grid_mlp=False, grid_branch=[0, 0, 1], ffn_blocks=1, so3_readout="mlp",
             precision="float64", channels=32, lmax=2, n_blocks=2, mixing_layers=3, radial_so2_rank=1, n_focus=1)
    d.update(kw)
    return d


DPA4 = {"dpa4_mini": _dpa4_base(),                                                     # the grade 'mini' of deepmd's presets
        "dpa4_nano": _dpa4_base(lmax=1, radial_so2_mode="none"),                       # the grade 'nano'
        "dpa4_b": _dpa4_base(channels=16, n_blocks=3, mixing_layers=2, rcut=5.0),      # other sizes
        "dpa4_s": _dpa4_base(channels=8),                                              # small: physics regressions
        "dpa4_s1": _dpa4_base(channels=8, lmax=1, radial_so2_mode="none")}


def build_dpa4(name, seed=1):
    """DPA4: fresh models zero-initialize their output layers: random weights (reference with `dp convert-backend x.dp x.pt2`)."""
    import os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from deepmd_dpa4 import randomize
    d = dict(DPA4[name], seed=seed)
    params = {"type_map": ["O", "H"], "descriptor": d, "fitting_net": {"type": "dpa4_ener", "neuron": [0], "precision": "float64", "seed": seed}}
    m = get_model(params)
    randomize(m.atomic_model.descriptor, seed + 6)
    randomize(m.atomic_model.fitting_net, seed + 7)
    data = m.serialize()
    data["@variables"]["out_bias"] = np.array([[[-90.0], [-180.0]]])
    return {"model": data, "model_def_script": params, "@variables": {}, "backend": "DPModel", "software": "deepmd-kit"}


def build(name, seed=1):
    if name in DPA4:
        return build_dpa4(name, seed)
    rng = np.random.default_rng(seed)
    d = dict(DESCRIPTORS[name])
    if d["type"] in ("dpa3", "dpa2"):
        d["seed"] = seed
    elif d["type"] == "hybrid":
        d["list"] = [dict(x, seed=seed + k) for k, x in enumerate(d["list"])]
    elif d["type"] not in ("dpa3", "dpa2"):
        d["seed"] = seed
    params = {"type_map": ["O", "H"], "descriptor": d,
              "fitting_net": {"neuron": [16, 16, 16], "resnet_dt": True, "seed": seed}}
    m = get_model(params)
    data = m.serialize()
    if name.startswith("dpa2"):
        dv = data["descriptor"]
        for blk in ("repinit_variable", "repformers_variable", "repinit_three_body_variable"):
            if blk not in dv or dv[blk] is None:
                continue
            v = dv[blk]["@variables"]
            a, sd = np.asarray(v["davg"]), np.asarray(v["dstd"])
            a[...] = 0.0
            sd[..., 0] = rng.uniform(0.2, 0.4); sd[..., 1:] = rng.uniform(0.2, 0.4)
            if name == "dpa2_c" and blk == "repformers_variable":
                a[..., 0] = 0.15
            v["davg"], v["dstd"] = a, sd
        for lay in dv["repformers_variable"]["repformer_layers"]:
            for k, x in lay["@variables"].items():
                if "residual" in k:
                    lay["@variables"][k] = [rng.uniform(0.05, 0.3, size=np.asarray(r).shape) for r in x]
        data["fitting"]["@variables"]["bias_atom_e"] = np.array([[-90.0], [-180.0]])
        return {"model": data, "model_def_script": params, "@variables": {}, "backend": "DPModel", "software": "deepmd-kit"}
    if name.startswith("dpa3"):
        v = data["descriptor"]["repflow_variable"]["@variables"]
        a, sd = np.asarray(v["davg"]), np.asarray(v["dstd"])
        a[...] = 0.0
        sd[..., 0] = rng.uniform(0.2, 0.4); sd[..., 1:] = rng.uniform(0.2, 0.4)
        v["davg"], v["dstd"] = a, sd
        # random residual weights (deepmd initializes them to update_residual)
        for lay in data["descriptor"]["repflow_variable"]["repflow_layers"]:
            for k, x in lay["@variables"].items():
                if "residual" in k:
                    lay["@variables"][k] = [rng.uniform(0.05, 0.3, size=np.asarray(r).shape) for r in x] if isinstance(x, list) else rng.uniform(0.05, 0.3, size=np.asarray(x).shape)
        data["fitting"]["@variables"]["bias_atom_e"] = np.array([[-90.0], [-180.0]])
        return {"model": data, "model_def_script": params, "@variables": {}, "backend": "DPModel", "software": "deepmd-kit"}
    if name == "hybrid":
        for k, sub in enumerate(data["descriptor"]["list"]):
            v = sub["@variables"]; a, sd = np.asarray(v["davg"]), np.asarray(v["dstd"])
            a[...] = 0.0
            sd[..., 0] = rng.uniform(0.05, 0.3); sd[..., 1:] = rng.uniform(0.05, 0.3)
            if sub["type"] == "se_e2_a": a[..., 0] = 0.1
            v["davg"], v["dstd"] = a, sd
        data["fitting"]["@variables"]["bias_atom_e"] = np.array([[-90.0], [-180.0]])
        return {"model": data, "model_def_script": params, "@variables": {}, "backend": "DPModel", "software": "deepmd-kit"}
    var = data["descriptor"]["@variables"]
    avg, std = np.asarray(var["davg"]), np.asarray(var["dstd"])
    nt, nnei, nc = avg.shape
    if name in MIXED:
        for i in range(nt):
            avg[i, :, :] = 0.0
            if MIXED[name]: avg[i, :, 0] = rng.uniform(0.0, 0.2)
            std[i, :, 0] = rng.uniform(0.05, 0.3)
            std[i, :, 1:] = rng.uniform(0.05, 0.3)
        var["davg"], var["dstd"] = avg, std
        data["fitting"]["@variables"]["bias_atom_e"] = np.array([[-90.0], [-180.0]])
        return {"model": data, "model_def_script": params, "@variables": {}, "backend": "DPModel", "software": "deepmd-kit"}
    # per (center type, neighbor type) constant statistics, angular averages zero (as deepmd's)
    sel = d["sel"]; start = 0
    for t, s in enumerate(sel):
        for i in range(nt):
            avg[i, start:start + s, 0] = rng.uniform(0.0, 0.2)
            std[i, start:start + s, 0] = rng.uniform(0.05, 0.3)
            if nc == 4: std[i, start:start + s, 1:] = rng.uniform(0.05, 0.3)
        start += s
    var["davg"], var["dstd"] = avg, std
    data["fitting"]["@variables"]["bias_atom_e"] = np.array([[-90.0], [-180.0]])
    return {"model": data, "model_def_script": params, "@variables": {}, "backend": "DPModel", "software": "deepmd-kit"}

for name in sys.argv[1:]:
    save_dp_model(f"{name}.dp", build(name))
    print("wrote", f"{name}.dp")
