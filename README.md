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
- Supports **mixed catalogs**: triangles with vertices drawn from different point sets  
  (e.g. `ABC`, `AAB`, or `ABB`). Triangles are deduplicated when catalogs repeat  
  (e.g. only `i<j` counted if two vertices are from the same catalog).

---

## Installation

Clone the repository and develop locally:

```bash
git clone https://github.com/ladosamushia/TriCo.git
cd TriCo
julia --project=.
```

In the Julia REPL:

```julia
using Pkg
Pkg.instantiate()
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

### Mixed-catalog examples

You can also count triangles using multiple catalogs. Wrap each point set with `TriCat(X,Y,Z)`:

```julia
using TriCo

# Three independent catalogs
A = TriCo.TriCat(randn(1000), randn(1000), randn(1000))
B = TriCo.TriCat(randn(1000), randn(1000), randn(1000))
C = TriCo.TriCat(randn(1000), randn(1000), randn(1000))

# One from each (ABC)
Habc = count_triangles_mixed!(A,B,C; rmin=5.0, rmax=60.0, Nr=55, μmax=0.9, Nμ=2)

# Two from A, one from B (AAB)
Haab = count_triangles_mixed!(A,A,B; rmin=5.0, rmax=60.0, Nr=55, μmax=0.9, Nμ=2)

# One from A, two from B (ABB)
Habb = count_triangles_mixed!(A,B,B; rmin=5.0, rmax=60.0, Nr=55, μmax=0.9, Nμ=2)

@show sum(Habc.h), sum(Haab.h), sum(Habb.h)
```

The mixed functions also have periodic versions:

```julia
HabcP = count_triangles_mixed_periodic!(A,B,C;
                                        Lx=2000, Ly=2000, Lz=2000,
                                        rmin=5.0, rmax=60.0, Nr=55,
                                        μmax=0.9, Nμ=2)
```

*Implementation details:*  
Triangles are deduplicated by enforcing index order only within repeated catalogs  
(e.g. if `A===B`, then only `i<j` is accepted). After passing cuts, triangle sides  
are sorted by length before binning, so the histogram is invariant to vertex order.

---

## Example scripts

Scripts in `scripts/` provide ready-to-run examples:

- **`scripts/random_cube.jl`**  
  Generate a cube of random points and count triangles.  
  Controlled via environment variables (see below).

- **`scripts/fits_to_hist.jl`**  
  Build histograms directly from FITS files and save to `.npz`.

Run them with multiple threads:

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/random_cube.jl
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_hist.jl --fits galaxies.fits
```

---

## `random_cube.jl` usage

This script uses environment variables for configuration.

**Example (all options set):**

```bash
TRICO_N=100000 TRICO_L=1500 TRICO_RMIN=2.0 TRICO_RMAX=50.0 TRICO_NR=40 TRICO_MUMAX=0.8 TRICO_NMU=4 TRICO_CELL=50.0 TRICO_SEED=42 TRICO_PERIODIC=1 JULIA_NUM_THREADS=8 julia --project=. scripts/random_cube.jl
```

**Defaults (if not set):**

- `TRICO_N=200000` — number of points
- `TRICO_L=2000.0` — box length (same for x,y,z)
- `TRICO_RMIN=5.0`, `TRICO_RMAX=60.0` — radial bin range
- `TRICO_NR=55` — number of r bins
- `TRICO_MUMAX=0.9` — maximum μ
- `TRICO_NMU=2` — number of μ bins
- `TRICO_CELL=rmax` — neighbor grid cell size
- `TRICO_SEED=12345` — RNG seed
- `TRICO_PERIODIC=0` — 0 = non-periodic, 1 = periodic

---

## `fits_to_hist.jl` usage

This script builds histograms directly from FITS files.
It supports both **single-catalog** and **mixed-catalog** counting.

### Single catalog (AAA)

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_hist.jl   --fits galaxies.fits   --xcol X --ycol Y --zcol Z --hdu 2   --rmin 5 --rmax 60 --Nr 55   --mumax 0.9 --Nmu 2   --cellsize 60   --out triangles.npz
```

### Mixed catalogs (AAB / ABB / ABC)

Two catalogs (default = AAB if `--pattern` omitted):

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_hist.jl   --fitsA catA.fits --fitsB catB.fits   --pattern AAB   --rmin 5 --rmax 60 --Nr 55   --mumax 0.9 --Nmu 2   --cellsize 60   --out triangles_AAB.npz
```

Three catalogs (ABC):

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_hist.jl   --fitsA catA.fits --fitsB catB.fits --fitsC catC.fits   --pattern ABC   --rmin 5 --rmax 60 --Nr 55   --mumax 0.9 --Nmu 2   --cellsize 60   --out triangles_ABC.npz
```

Periodic mode (works in both single and mixed):

```bash
--periodic --Lx 2000 --Ly 2000 --Lz 2000
```

---

### Options

**Catalog inputs**
- `--fits PATH` — single-catalog mode
- `--fitsA PATH` — mixed mode, catalog A
- `--fitsB PATH` — mixed mode, catalog B (optional, required for ABB/ABC)
- `--fitsC PATH` — mixed mode, catalog C (optional, required for ABC)

**Mixed composition**
- `--pattern AAB|ABB|ABC` — composition of triangle vertices
  Defaults:
  • 2 catalogs → `AAB`
  • 3 catalogs → `ABC`

**FITS columns (apply globally unless per-catalog given)**
- `--xcol=X`, `--ycol=Y`, `--zcol=Z` — column names (defaults: `X,Y,Z`)
- `--hdu=2` — HDU index (default: 2)
- `--xcolA`, `--ycolB`, `--zcolC`, `--hduB`, etc. — per-catalog overrides

**Selection & binning**
- `--rmin=5.0`, `--rmax=60.0` — radial bin range
- `--Nr=55` — number of r bins
- `--mumax=0.9` — maximum μ (0 < μ ≤ 1)
- `--Nmu=2` — number of μ bins

**Neighbor grid**
- `--cellsize=rmax` — grid cell size (must be ≥ rmax, default = rmax)

**Periodic box (z = LOS)**
- `--periodic` — enable periodic mode (default: off)
- `--Lx=2000`, `--Ly=2000`, `--Lz=2000` — box sizes (required if `--periodic`)

**Output**
- `--out FILE.npz` — save histogram to NPZ (default: none, must end with `.npz`)

## `fits_to_pairs.jl` usage

This script builds **pair histograms** directly from FITS files.
It supports both **single-catalog** (auto-pairs) and **cross-catalog** counting.

### Single catalog

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_pairs.jl \
  --fits galaxies.fits \
  --xcol X --ycol Y --zcol Z --hdu 2 \
  --rmin 5 --rmax 60 --Nr 55 \
  --mumax 0.9 --Nmu 2 \
  --cellsize 60 \
  --out pairs.npz
```

### Cross catalogs (A vs B)

```bash
JULIA_NUM_THREADS=8 julia --project=. scripts/fits_to_pairs.jl \
  --fitsA galaxies_A.fits --fitsB galaxies_B.fits \
  --xcolA X --ycolA Y --zcolA Z --hduA 2 \
  --xcolB X --ycolB Y --zcolB Z --hduB 2 \
  --rmin 5 --rmax 60 --Nr 55 \
  --mumax 0.9 --Nmu 2 \
  --cellsize 60 \
  --out pairs_AB.npz
```

### Periodic mode (works in both single and cross)

```bash
--periodic --Lx 2000 --Ly 2000 --Lz 2000
```

---

### Options

**Catalog inputs**
- `--fits PATH` — single-catalog mode
- `--fitsA PATH` — cross mode, catalog A
- `--fitsB PATH` — cross mode, catalog B

**FITS columns (apply globally unless per-catalog given)**
- `--xcol=X`, `--ycol=Y`, `--zcol=Z` — column names (defaults: `X,Y,Z`)
- `--hdu=2` — HDU index (default: 2)
- `--xcolA`, `--ycolB`, `--zcolA`, `--hduA`, `--hduB`, etc. — per-catalog overrides

**Selection & binning**
- `--rmin=5.0`, `--rmax=60.0` — radial bin range
- `--Nr=55` — number of r bins
- `--mumax=0.9` — maximum μ (0 < μ ≤ 1)
- `--Nmu=2` — number of μ bins

**Neighbor grid**
- `--cellsize=rmax` — grid cell size (must be ≥ rmax, default = rmax)

**Periodic box (z = LOS)**
- `--periodic` — enable periodic mode (default: off)
- `--Lx=2000`, `--Ly=2000`, `--Lz=2000` — box sizes (required if `--periodic`)

**Output**
- `--out FILE.npz` — save histogram to NPZ (default: none, must end with `.npz`)

---

### Defaults summary

- Catalog: `--fits` required for single; `--fitsA` and `--fitsB` required for cross
- Columns: `X,Y,Z`
- HDU: `2`
- rmin: `5.0`
- rmax: `60.0`
- Nr: `55`
- μmax: `0.9`
- Nμ: `2`
- cellsize: `rmax`
- periodic: `false`
- Lx,Ly,Lz: `2000.0` (only used if periodic)
- out: none (no file saved unless `--out` specified)

---

### Output format

The `.npz` file contains:

- `hist` — the 2D histogram array of shape `(Nr, Nμ)`

No metadata or bin edges are saved in this version.
You can compute bin edges separately in Julia/Python if needed.

---

### Defaults summary

- Catalog: `--fits` required for single mode; `--fitsA/B/C` required for mixed
- Columns: `X,Y,Z`
- HDU: `2`
- rmin: `5.0`
- rmax: `60.0`
- Nr: `55`
- μmax: `0.9`
- Nμ: `2`
- cellsize: `rmax`
- periodic: `false`
- Lx,Ly,Lz: `2000.0` (only used if periodic)
- pattern: depends on number of catalogs (AAB if 2, ABC if 3)
- out: none (no file saved unless `--out` specified)

---

### Output format

The `.npz` file contains:

- `hist` — the 4D histogram array of shape `(Nr, Nr, Nμ, Nμ)`

No metadata or bin edges are saved in this version.  
You can compute bin edges separately in Julia/Python if needed.

---

## Tests

Unit tests are in `test/`. They include:

- **Brute-force comparison tests** that generate small random point sets, compute triangle histograms by exhaustive search, and verify **bin-by-bin equality** against the accelerated grid-based implementation (`count_triangles_grid!`, `count_triangles_periodic_grid!`).

Run tests with:

```bash
JULIA_NUM_THREADS=1 julia --project=. test/runtests.jl
```

---

## Development notes

Current repo structure:

```
src/
  TriCo.jl
  triangles.jl
  triangles_mixed.jl   # ← new
  io_save.jl
scripts/
  random_cube.jl
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

