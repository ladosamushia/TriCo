"""
    triplet_geometry_periodic(x1,y1,z1, x2,y2,z2, x3,y3,z3; L=(2000.0, 2000.0,
2000.0))

Return `(r12_sq, r23_sq, r31_sq, μ12_sq, μ23_sq, μ31_sq)` where:
- `rij_sq` = squared distance between points i and j.
- `μij_sq` = (Δz_ij)^2 / rij_sq = squared cosine of the angle of vector i→j w.r.t. z-axis.
  If `rij_sq == 0`, the corresponding `μij_sq` is `NaN`.

Anisotropic minimum-image convention:
- Pass `L = (Lx, Ly, Lz)`; each axis uses its own period. Default is 2000.0 for
  all three axes (AbacusSummit size).
"""
function triplet_geometry_periodic(x1, y1, z1, x2, y2, z2, x3, y3, z3;
L::NTuple{3,<:Real}=(2000.0,2000.0,2000.0))
    Lx, Ly, Lz = L
    @inline wrap_axis(d, L) = (d - L*round(d/L))
    @inline r2_mu2(dx, dy, dz) = begin
        r2 = dx*dx + dy*dy + dz*dz
        mu2 = (dz*dz) / r2
        r2, mu2
    end

    # 1 → 2
    dx12 = wrap_axis(x2 - x1, Lx)
    dy12 = wrap_axis(y2 - y1, Ly)
    dz12 = wrap_axis(z2 - z1, Lz)
    r12_sq, μ12_sq = r2_mu2(dx12, dy12, dz12)

    # 2 → 3
    dx23 = wrap_axis(x3 - x2, Lx)
    dy23 = wrap_axis(y3 - y2, Ly)
    dz23 = wrap_axis(z3 - z2, Lz)
    r23_sq, μ23_sq = r2_mu2(dx23, dy23, dz23)

    # 3 → 1
    dx31 = wrap_axis(x1 - x3, Lx)
    dy31 = wrap_axis(y1 - y3, Ly)
    dz31 = wrap_axis(z1 - z3, Lz)
    r31_sq, μ31_sq = r2_mu2(dx31, dy31, dz31)

    return r12_sq, r23_sq, r31_sq, μ12_sq, μ23_sq, μ31_sq
end

"""
    triplet_geometry_los(x1,y1,z1, x2,y2,z2, x3,y3,z3)

Return `(r12_sq, r23_sq, r31_sq, μ12_sq, μ23_sq, μ31_sq)` for a non-periodic volume, where:
- `rij_sq` = squared distance between points i and j.
- `μij_sq` = squared cosine between the separation vector `s_ij = r_j - r_i`
  and the **true line of sight** from the origin to the pair midpoint.
  Concretely, `μij_sq = ((s_ij ⋅ (r_i + r_j))^2) / (|s_ij|^2 * |r_i + r_j|^2)`.
"""
function triplet_geometry_true_los(x1, y1, z1, x2, y2, z2, x3, y3, z3)
    # positions
    # r1 = (x1,y1,z1); r2 = (x2,y2,z2); r3 = (x3,y3,z3)

    @inline function r2_mu2_true(dx, dy, dz, mx, my, mz)
        r2 = dx*dx + dy*dy + dz*dz          # |s_ij|^2
        m2 = mx*mx + my*my + mz*mz          # |r_i + r_j|^2  (midpoint direction, scale cancels)
        if r2 == 0 || m2 == 0
            return r2, NaN
        end
        dot = dx*mx + dy*my + dz*mz
        μ2 = (dot*dot) / (r2*m2)
        return r2, μ2
    end

    # 1 → 2
    dx12 = x2 - x1; dy12 = y2 - y1; dz12 = z2 - z1
    mx12 = x1 + x2; my12 = y1 + y2; mz12 = z1 + z2   # midpoint LoS ~ (r1 + r2)
    r12_sq, μ12_sq = r2_mu2_true(dx12, dy12, dz12, mx12, my12, mz12)

    # 2 → 3
    dx23 = x3 - x2; dy23 = y3 - y2; dz23 = z3 - z2
    mx23 = x2 + x3; my23 = y2 + y3; mz23 = z2 + z3
    r23_sq, μ23_sq = r2_mu2_true(dx23, dy23, dz23, mx23, my23, mz23)

    # 3 → 1
    dx31 = x1 - x3; dy31 = y1 - y3; dz31 = z1 - z3
    mx31 = x3 + x1; my31 = y3 + y1; mz31 = z3 + z1
    r31_sq, μ31_sq = r2_mu2_true(dx31, dy31, dz31, mx31, my31, mz31)

    return r12_sq, r23_sq, r31_sq, μ12_sq, μ23_sq, μ31_sq
end

