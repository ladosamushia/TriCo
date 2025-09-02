# Simple uniform grid (cell-linked list) with CSR storage for indices.

struct Grid
    xmin::Float64; ymin::Float64; zmin::Float64
    invh::Float64
    nx::Int; ny::Int; nz::Int
    offs::Vector{Int}    # length ncell+1 (CSR)
    idxs::Vector{Int}    # point indices sorted into cells
end

function build_grid(X::AbstractVector{<:Real},
                    Y::AbstractVector{<:Real},
                    Z::AbstractVector{<:Real},
                    h::Real)
    N = length(X)
    @assert length(Y) == N == length(Z)
    xmin = float(minimum(X)); xmax = float(maximum(X))
    ymin = float(minimum(Y)); ymax = float(maximum(Y))
    zmin = float(minimum(Z)); zmax = float(maximum(Z))
    invh = 1 / float(h + eps())  # tiny eps to keep right edge closed

    nx = max(1, Int(floor((xmax - xmin)*invh)) + 1)
    ny = max(1, Int(floor((ymax - ymin)*invh)) + 1)
    nz = max(1, Int(floor((zmax - zmin)*invh)) + 1)
    ncell = nx*ny*nz

    cellkey = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        ix = Int(floor((float(X[i]) - xmin)*invh)); ix = ifelse(ix>=nx, nx-1, max(ix,0))
        iy = Int(floor((float(Y[i]) - ymin)*invh)); iy = ifelse(iy>=ny, ny-1, max(iy,0))
        iz = Int(floor((float(Z[i]) - zmin)*invh)); iz = ifelse(iz>=nz, nz-1, max(iz,0))
        cellkey[i] = ix + nx*(iy + ny*iz)    # 0-based
    end

    counts = zeros(Int, ncell)
    @inbounds for k in cellkey
        counts[k+1] += 1
    end

    offs = Vector{Int}(undef, ncell+1)
    offs[1] = 1
    @inbounds for c in 1:ncell
        offs[c+1] = offs[c] + counts[c]
    end

    writepos = copy(offs)
    idxs = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        cid = cellkey[i] + 1
        p = writepos[cid]
        idxs[p] = i
        writepos[cid] = p + 1
    end

    Grid(xmin, ymin, zmin, invh, nx, ny, nz, offs, idxs)
end

@inline function cell_of(g::Grid, x::Float64, y::Float64, z::Float64)
    ix = Int(floor((x - g.xmin)*g.invh)); ix = ifelse(ix>=g.nx, g.nx-1, ifelse(ix<0,0,ix))
    iy = Int(floor((y - g.ymin)*g.invh)); iy = ifelse(iy>=g.ny, g.ny-1, ifelse(iy<0,0,iy))
    iz = Int(floor((z - g.zmin)*g.invh)); iz = ifelse(iz>=g.nz, g.nz-1, ifelse(iz<0,0,iz))
    return ix,iy,iz
end

@inline cell_id(g::Grid, ix,iy,iz) = ix + g.nx*(iy + g.ny*iz) + 1  # 1-based

struct GridPeriodic
    Lx::Float64; Ly::Float64; Lz::Float64
    invhx::Float64; invhy::Float64; invhz::Float64
    nx::Int; ny::Int; nz::Int
    offs::Vector{Int}    # CSR offsets, length ncell+1
    idxs::Vector{Int}    # point indices sorted by cell
end

function build_grid_periodic(X::AbstractVector{<:Real},
                             Y::AbstractVector{<:Real},
                             Z::AbstractVector{<:Real},
                             Lx::Real, Ly::Real, Lz::Real,
                             h::Real)
    N = length(X); @assert length(Y)==N==length(Z)
    Lx = float(Lx); Ly = float(Ly); Lz = float(Lz)
    invhx = 1/float(h + eps()); invhy = invhx; invhz = invhx
    nx = max(1, Int(ceil(Lx*invhx))); ny = max(1, Int(ceil(Ly*invhy))); nz = max(1, Int(ceil(Lz*invhz)))
    ncell = nx*ny*nz

    # map point -> cell (0-based, modulo box)
    cellkey = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        ix = Int(floor(float(X[i])*invhx)) % nx
        iy = Int(floor(float(Y[i])*invhy)) % ny
        iz = Int(floor(float(Z[i])*invhz)) % nz
        cellkey[i] = ix + nx*(iy + ny*iz)
    end

    counts = zeros(Int, ncell)
    @inbounds for k in cellkey
        counts[k+1] += 1
    end

    offs = Vector{Int}(undef, ncell+1)
    offs[1] = 1
    @inbounds for c in 1:ncell
        offs[c+1] = offs[c] + counts[c]
    end

    writepos = copy(offs)
    idxs = Vector{Int}(undef, N)
    @inbounds for i in 1:N
        cid = cellkey[i] + 1
        p = writepos[cid]
        idxs[p] = i
        writepos[cid] = p + 1
    end

    GridPeriodic(Lx, Ly, Lz, invhx, invhy, invhz, nx, ny, nz, offs, idxs)
end

@inline function cell_of(g::GridPeriodic, x::Float64, y::Float64, z::Float64)
    cx = Int(floor(x*g.invhx)) % g.nx
    cy = Int(floor(y*g.invhy)) % g.ny
    cz = Int(floor(z*g.invhz)) % g.nz
    return cx, cy, cz
end

@inline cell_id(g::GridPeriodic, cx,cy,cz) = (cx + g.nx*(cy + g.ny*cz)) + 1

# periodic wrap of neighbor cell index
@inline wrap(i::Int, n::Int) = (i % n + n) % n

