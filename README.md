# TriCo

**TriCo** is a Julia package for fast triangle counting in 3D point sets.  
It is designed for cosmology applications (e.g. galaxy clustering) but can be applied
to any large point cloud where three-point statistics are needed.

---

## Features

- Counts triangles with constraints on side lengths `r_ij` and line-of-sight cosines `μ_ij`.
- Bins triangles into a **4D histogram** of `(r12, r23, μ12, μ13)`.
- Supports both:
  - **Non-periodic boxes** with true line-of-sight (pair midpoint definition).
  - **Periodic boxes** with arbitrary lengths in `x,y,z` and **z as LOS** (`μ² = Δz² / r²`, min-image convention).
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

TriCo is not (yet) in the Julia registry; use it via `--project` or `Pkg.develop`.

---

## Usage from Julia

### Non-periodic example

```julia
using TriCo

N = 10_000
X, Y, Z = randn(N), randn(N), randn(N)

H = count_triangles_grid!(X, Y, Z;
                          rmin=5.0, rmax=60.0, Nr=55,
                          μmax=0.9, Nμ=2, cellsize=60.0)

@show size(H.h)   # (Nr, Nr, Nμ, Nμ)
@show sum(H.h)    # total number of triangles
```

### Periodic example

```julia
using TriCo

N = 10_000
Lx, Ly, Lz = 2000.0, 2000.0, 2000.0
X, Y, Z = Lx .* rand(N), Ly .* rand(N), Lz .* rand(N)

H = count_triangles_periodic_grid!(X, Y, Z;
                                   Lx=Lx, Ly=Ly, Lz=Lz,
                                   rmin=5.0, rmax=60.0, Nr=55,
                                   μmax=0.9, Nμ=2, cellsize=60.0)

@show size(H.h)
@show sum(H.h)
```

---

## Example scripts

Scripts in `scripts/` provide ready-to-run examples:

- **`scripts/random_cube.jl`**  
  Non-periodic cube with random points.

- **`scripts/random_cube_periodic.jl`**  
  Periodic cube (with arbitrary box sizes).

- **`scripts/fits_to_hist.jl`**  
  Build histograms directly from FITS files and save to `.npz`.

Run them with multiple threads:

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/random_cube.jl
JULIA_NUM_THREADS=8 julia --project=. scripts/random_cube_periodic.jl
```

Parameters (e.g., number of points, box size, bins) can be overridden via environment variables or command-line options.

---

## Command-line usage: `fits_to_hist.jl`

This script builds histograms directly from FITS files.

**Non-periodic example:**
```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_hist.jl   --fits galaxies.fits --xcol X --ycol Y --zcol Z   --rmin 5 --rmax 60 --Nr 55 --mumax 0.9 --Nmu 2   --out triangles.npz
```

**Periodic example (z is LOS):**
```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_hist.jl   --fits galaxies.fits --xcol X --ycol Y --zcol Z   --periodic --Lx 2000 --Ly 2000 --Lz 2000   --rmin 5 --rmax 60 --Nr 55 --mumax 0.9 --Nmu 2   --out triangles_box.npz
```

### Output format

The `.npz` file contains:

- `hist` — the 4D histogram array of shape `(Nr, Nr, Nμ, Nμ)`

No metadata or bin edges are saved in this version.  
You can compute bin edges separately in Julia/Python if needed.

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
    io_save.jl
  scripts/
    random_cube.jl
    random_cube_periodic.jl
    fits_to_hist.jl
  test/
    runtests.jl
  ```

- Only edit `Project.toml` when adding/removing dependencies.
- `Manifest.toml` is auto-generated.
- For smoother development, use [Revise.jl](https://timholy.github.io/Revise.jl/stable/) to pick up changes automatically.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

