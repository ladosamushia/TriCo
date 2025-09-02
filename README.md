# TriCo

**TriCo** is a Julia package for fast triangle counting in 3D point sets.  
It is designed for cosmology applications (e.g. galaxy clustering) but can be applied
to any large point cloud where three-point statistics are needed.

Features:

- Counts triangles with constraints on side lengths `r_ij` and line-of-sight cosines `μ_ij`.
- Bins triangles into histograms of `(r12, r23, μ12, μ13)` (4D joint histogram).
- Supports both **non-periodic** and **periodic** boxes (with different lengths in x, y, z).
- Periodic mode treats **z as the line of sight** (`μ² = Δz² / r²` with minimum-image distances).
- Parallelized with `Threads.@threads` and designed for HPC workloads.

---

## Installation

Clone the repository and develop locally:

```bash
git clone https://github.com/yourusername/TriCo.git
cd TriCo
julia --project=.
```

In the Julia REPL:

```julia
using Pkg
Pkg.instantiate()   # install dependencies if needed
```

TriCo is not in the Julia registry; use it via `--project` or `Pkg.develop`.

---

## Usage

### Non-periodic cube

```julia
using TriCo

# Generate some points
N = 10_000
X, Y, Z = randn(N), randn(N), randn(N)

# Count triangles
H = count_triangles!(X, Y, Z; rmin=5.0, rmax=60.0, Nr=55,
                     μmax=0.9, Nμ=2, cellsize=60.0)

@show size(H.h)     # (Nr, Nr, Nμ, Nμ)
@show sum(H.h)      # total number of triangles
```

### Periodic cube (different lengths allowed)

```julia
using TriCo

N = 10_000
Lx, Ly, Lz = 2000.0, 2000.0, 2000.0
X, Y, Z = Lx .* rand(N), Ly .* rand(N), Lz .* rand(N)

H = count_triangles_periodic!(X, Y, Z;
                              Lx=Lx, Ly=Ly, Lz=Lz,
                              rmin=5.0, rmax=60.0, Nr=55,
                              μmax=0.9, Nμ=2, cellsize=60.0)

@show size(H.h)
@show sum(H.h)
```

---

## Example scripts

We provide two ready-to-run scripts in `scripts/`:

- **`scripts/random_cube.jl`**  
  Non-periodic cube with random points.

- **`scripts/random_cube_periodic.jl`**  
  Periodic cube (with optional different `Lx, Ly, Lz`).

Run them with multiple threads, e.g.:

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/random_cube.jl
JULIA_NUM_THREADS=8 julia --project=. scripts/random_cube_periodic.jl
```

Parameters (number of points, box size, bins, etc.) can be overridden via environment variables, e.g.:

```bash
TRICO_N=20000 TRICO_RMAX=80 TRICO_NR=60 JULIA_NUM_THREADS=auto julia --project=. scripts/random_cube_periodic.jl
```

---

## Development notes

- TriCo is structured as a Julia package:
  ```
  src/
    TriCo.jl
    geometry.jl
    binning.jl
    grid.jl
    triangles.jl
  scripts/
    random_cube.jl
    random_cube_periodic.jl
  test/
    runtests.jl
  ```

- Only edit `Project.toml` when adding/removing dependencies.
- `Manifest.toml` is machine-generated — no need to edit it manually.
- For smooth development, install [Revise.jl](https://timholy.github.io/Revise.jl/stable/)
  so code changes are picked up without restarting Julia.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

