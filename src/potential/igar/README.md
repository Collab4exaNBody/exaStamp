# IGAR — Force from Grid Cell Values

## Overview

The `igar_force` operator computes per-particle potential energy and forces from a scalar energy field stored on the simulation grid (`grid_cell_values`). For each particle, it locates the particle's subcell, reads the 3×3×3 neighborhood of **subcell** energy values centered on it, then applies **triquadratic Lagrange interpolation** to estimate the energy and its gradient at the particle's exact position. The interpolation operates at full subcell resolution for any value of `grid_subdiv`.

---

## Algorithm

### 1. Input data

- `grid_cell_values`: holds a scalar field `"ep"` — one energy value per cell (or `subdiv³` subcell values per cell when `grid_subdiv > 1`).
- `grid`: the particle grid, providing cell dimensions, ghost layers, and per-cell particle arrays (`rx/ry/rz`, `fx/fy/fz`, `ep`).

---

### 2. Per-particle subcell stencil construction

For each particle `p` in cell `C`:

1. Compute `rco = r - cell_origin` (position relative to cell corner).
2. Call `localize_subcell(rco, ...)` to find which subcell `(si, sj, sk)` contains the particle.
3. Build a 3×3×3 stencil of **subcell** energy values centered on `(si, sj, sk)`:

```
stencil[ck+1][cj+1][ci+1]  for  ci,cj,ck ∈ {-1, 0, +1}
```

Each entry is a single subcell energy value. `subcell_neighbor` handles cross-cell-boundary indexing transparently when the neighboring subcell falls in an adjacent cell. Out-of-domain entries are set to 0.

> The stencil is per-particle because different particles in the same cell may occupy different subcells.

---

### 3. Triquadratic Lagrange interpolation

#### Normalized coordinates

For a particle at position **r**, the normalized offset from the **subcell center** `(si, sj, sk)` is:

```
u = rco.x / subcell_size - (si + 0.5)
v = rco.y / subcell_size - (sj + 0.5)
w = rco.z / subcell_size - (sk + 0.5)
```

with `u, v, w ∈ [-0.5, 0.5]` and `subcell_size = cell_size / subdiv`.

The three neighboring subcell centers in these units are located at `-1`, `0`, `+1` — exactly the Lagrange interpolation nodes.

#### Lagrange basis functions (per axis, nodes at −1, 0, +1)

| Basis | Value | Derivative |
|-------|-------|------------|
| `L₀(t)` | `t(t−1)/2` | `t − 1/2` |
| `L₁(t)` | `1 − t²` | `−2t` |
| `L₂(t)` | `t(t+1)/2` | `t + 1/2` |

#### Tensor-product interpolation

Energy at the particle position:

```
E(u,v,w) = Σᵢⱼₖ  Lᵢ(u) · Lⱼ(v) · Lₖ(w) · stencil[k][j][i]
```

Gradient (in normalized coordinates):

```
∂E/∂u = Σᵢⱼₖ  L'ᵢ(u) · Lⱼ(v) · Lₖ(w) · stencil[k][j][i]
∂E/∂v = Σᵢⱼₖ  Lᵢ(u) · L'ⱼ(v) · Lₖ(w) · stencil[k][j][i]
∂E/∂w = Σᵢⱼₖ  Lᵢ(u) · Lⱼ(v) · L'ₖ(w) · stencil[k][j][i]
```

Chain rule converts to physical coordinates:

```
∂E/∂x = (1 / subcell_size) · ∂E/∂u
```

---

### 4. Force and energy accumulation

```
ep[p] += E(u,v,w)
fx[p] += −∂E/∂x
fy[p] += −∂E/∂y
fz[p] += −∂E/∂z
```

Forces are **added** (not assigned), consistent with multi-potential accumulation in exaStamp.

---

## Properties

| Property | Value |
|---|---|
| Interpolation order | Quadratic (degree 2 in each direction) |
| Continuity | C¹ — energy and forces are continuous across cell boundaries |
| Stencil | 3×3×3 = 27 cell-averaged energy values |
| Exact at nodes | Yes — interpolant reproduces exact cell-center values |
| Subcell support | Full — interpolates at native subcell resolution for any `subdiv` |
| Complexity | O(27) stencil reads + O(27) FMAs per particle |

---

## Coordinate system

```
cell_origin    = grid.cell_position(cell_loc)        # corner of cell, not center
subcell_size   = cell_size / subdiv
(si, sj, sk)   = center_subcell_loc                  # particle's subcell index in [0, subdiv-1]^3

rco = r - cell_origin                                # ∈ [0, cell_size]^3

u = rco.x / subcell_size - (si + 0.5)               # ∈ [-0.5, 0.5]
                                                     #    ^
                                                     # (si+0.5) shifts origin to subcell center
```

Stencil index mapping:

```
stencil[ck+1][cj+1][ci+1]
         ^                   → axis z (w coordinate)
              ^              → axis y (v coordinate)
                   ^         → axis x (u coordinate)
stencil[1][1][1]             → particle's own subcell (ci=cj=ck=0)
```

---

## Files

| File | Role |
|---|---|
| `igar.cu` | Operator node, slot declarations, `execute()` entry point |
| `include/force_from_grid.h` | `get_particle_force_from_grid()` — stencil build + interpolation |

---

## Known limitations

- Boundary particles whose subcell stencil extends outside the domain have missing entries set to 0, introducing an artificial gradient toward the boundary. Ghost layers normally cover all non-ghost subcells; verify ghost layer depth ≥ 1.
- Energy field name is hardcoded to `"ep"`.
