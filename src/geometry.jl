# --- basic squared distance ---
@inline dist2(x1::Real,y1::Real,z1::Real, x2::Real,y2::Real,z2::Real) =
    (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1)

# --- true line-of-sight μ^2 for pair (ri,rj) ---
# μ² = [(s ⋅ (ri + rj))²] / (|s|² · |ri + rj|²), where s = rj - ri
@inline function mu2_true_los(xi::Real,yi::Real,zi::Real,
                              xj::Real,yj::Real,zj::Real,
                              r2_ij::Real)
    sx = xj - xi; sy = yj - yi; sz = zj - zi         # separation
    lx = xi + xj; ly = yi + yj; lz = zi + zj         # LOS ≈ midpoint vector * 2
    num = (sx*lx + sy*ly + sz*lz)^2
    den = r2_ij * (lx*lx + ly*ly + lz*lz)
    return num / den
end

# -------- Periodic minimum-image helpers --------
@inline minimage(Δ::Float64, L::Float64) = Δ - L*round(Δ / L)

@inline function dist2_periodic(
    x1::Float64,y1::Float64,z1::Float64,
    x2::Float64,y2::Float64,z2::Float64,
    Lx::Float64,Ly::Float64,Lz::Float64
)
    dx = minimage(x2 - x1, Lx)
    dy = minimage(y2 - y1, Ly)
    dz = minimage(z2 - z1, Lz)
    return dx*dx + dy*dy + dz*dz, dx, dy, dz
end

# μ^2 for z-LOS under periodic minimum image: μ² = Δz² / r²
@inline mu2_zlos_from_deltas(dx::Float64, dy::Float64, dz::Float64, r2::Float64) =
    (dz*dz) / r2

