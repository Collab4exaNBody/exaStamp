#!/usr/bin/env python3
"""Run the same EAM sample at several initial temperatures in one Python process
and collect the final temperature and pressure of each run."""

import pyexastamp as xs

msp = xs.find_example("regression_new/potentials/eam/eam_vniitf/single_specy.msp")

results = []
for sigma_v in (25.0, 50.0, 100.0):  # width of the initial velocity noise, m/s
    # setup_system is a list: redefine it with the new noise amplitude
    ctx = xs.init(xs.command_line(msp, overrides={"global": {"max_iteration": 20}}))
    setup = xs.get_operator_defaults()["setup_system"]
    for item in setup:
        if isinstance(item, dict) and "gaussian_noise_v" in item:
            item["gaussian_noise_v"]["sigma"] = f"{sigma_v} m/s"
    xs.set_operator_defaults({"setup_system": setup})
    graph = xs.build_simulation_graph(ctx)
    xs.run_node(ctx, graph)
    t = xs.read_thermo(graph)
    results.append((sigma_v, t["temperature"], t["pressure"]))
    xs.end(ctx)

if ctx.mpi_rank != 0:   # thermodynamic values are global: print once
    raise SystemExit(0)

print(f"\n{'sigma_v (m/s)':>14} {'T (K)':>10} {'P (Pa)':>12}")
for sigma_v, temp, press in results:
    print(f"{sigma_v:14.1f} {temp:10.2f} {press:12.4g}")
