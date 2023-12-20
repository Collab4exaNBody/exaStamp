ccc_mprun -psklb -n1 -c48 /usr/bin/env OMP_NUM_THREADS=48 \
    ./xstampv2 stampv3_to_grid_vtk.msp \
        --set-global-default_cell_size "40.0 ang" --set-user_cell_size nop \
        --set-global-ghost_dist "40.0 ang" --set-user_ghost_dist nop \
        --set-global-grid_subdiv 3 --set-user_grid_subdiv nop \
        --set-global-splat_size "6.0 ang" --set-user_splat_size nop \
        --set-global-file "microjet.mpio" --set-user_source_file nop \
        --debug-config

