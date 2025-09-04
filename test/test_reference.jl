module TestReference
using TriCo
export reference_hist_nonperiodic, reference_hist_periodic, sumhist

sumhist(H) = sum(H.h)

@inline function cswap!(r_a, r_b, m_a, m_b, id_a, id_b)
    if (r_a > r_b) || (r_a == r_b && (m_a > m_b || (m_a == m_b && id_a > id_b)))
        r_a, r_b = r_b, r_a
        m_a, m_b = m_b, m_a
        id_a, id_b = id_b, id_a
    end
    return r_a, r_b, m_a, m_b, id_a, id_b
end

function reference_hist_nonperiodic(X, Y, Z; rmin, rmax, Nr, μmax, Nμ)
    N = length(X)
    H = zeros(Int, Nr, Nr, Nμ, Nμ)
    invΔr = Nr / (rmax - rmin)
    invΔμ = Nμ / μmax
    rmin2, rmax2 = rmin^2, rmax^2
    μmax2 = μmax^2

    @inbounds for i in 1:N-2, j in i+1:N-1, k in j+1:N
        x1=Float64(X[i]); y1=Float64(Y[i]); z1=Float64(Z[i])
        x2=Float64(X[j]); y2=Float64(Y[j]); z2=Float64(Z[j])
        x3=Float64(X[k]); y3=Float64(Y[k]); z3=Float64(Z[k])

        r12_sq = TriCo.dist2(x1,y1,z1, x2,y2,z2)
        r13_sq = TriCo.dist2(x1,y1,z1, x3,y3,z3)
        r23_sq = TriCo.dist2(x2,y2,z2, x3,y3,z3)
        (r12_sq <= rmin2 || r12_sq >= rmax2) && continue
        (r13_sq <= rmin2 || r13_sq >= rmax2) && continue
        (r23_sq <= rmin2 || r23_sq >= rmax2) && continue

        μ12_sq = TriCo.mu2_true_los(x1,y1,z1, x2,y2,z2, r12_sq)
        μ13_sq = TriCo.mu2_true_los(x1,y1,z1, x3,y3,z3, r13_sq)
        μ23_sq = TriCo.mu2_true_los(x2,y2,z2, x3,y3,z3, r23_sq)
        (μ12_sq >= μmax2 || μ13_sq >= μmax2 || μ23_sq >= μmax2) && continue

        r12b = Int(floor((sqrt(r12_sq) - rmin) * invΔr)) + 1
        r13b = Int(floor((sqrt(r13_sq) - rmin) * invΔr)) + 1
        r23b = Int(floor((sqrt(r23_sq) - rmin) * invΔr)) + 1
        μ12b = Int(floor(sqrt(μ12_sq) * invΔμ)) + 1
        μ13b = Int(floor(sqrt(μ13_sq) * invΔμ)) + 1
        μ23b = Int(floor(sqrt(μ23_sq) * invΔμ)) + 1
        (1 <= r12b <= Nr && 1 <= r13b <= Nr && 1 <= r23b <= Nr &&
         1 <= μ12b <= Nμ && 1 <= μ13b <= Nμ && 1 <= μ23b <= Nμ) || continue

        # sort by r-bin (short→long), carry μ with it; tie-break by μ then pair-id
        r1, r2, r3 = r12b, r23b, r13b
        m1, m2, m3 = μ12b, μ23b, μ13b
        id1, id2, id3 = 1, 2, 3

        r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)
        r2, r3, m2, m3, id2, id3 = cswap!(r2, r3, m2, m3, id2, id3)
        r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)

        H[r1, r2, m1, m2] += 1
    end
    H
end

function reference_hist_periodic(X, Y, Z; Lx, Ly, Lz, rmin, rmax, Nr, μmax, Nμ)
    N = length(X)
    H = zeros(Int, Nr, Nr, Nμ, Nμ)
    invΔr = Nr / (rmax - rmin)
    invΔμ = Nμ / μmax
    rmin2, rmax2 = rmin^2, rmax^2
    μmax2 = μmax^2

    @inbounds for i in 1:N-2, j in i+1:N-1, k in j+1:N
        x1=Float64(X[i]); y1=Float64(Y[i]); z1=Float64(Z[i])
        x2=Float64(X[j]); y2=Float64(Y[j]); z2=Float64(Z[j])
        x3=Float64(X[k]); y3=Float64(Y[k]); z3=Float64(Z[k])

        r12_sq, dx12, dy12, dz12 = TriCo.dist2_periodic(x1,y1,z1, x2,y2,z2, Lx,Ly,Lz)
        r13_sq, dx13, dy13, dz13 = TriCo.dist2_periodic(x1,y1,z1, x3,y3,z3, Lx,Ly,Lz)
        r23_sq, dx23, dy23, dz23 = TriCo.dist2_periodic(x2,y2,z2, x3,y3,z3, Lx,Ly,Lz)
        (r12_sq <= rmin2 || r12_sq >= rmax2) && continue
        (r13_sq <= rmin2 || r13_sq >= rmax2) && continue
        (r23_sq <= rmin2 || r23_sq >= rmax2) && continue

        μ12_sq = TriCo.mu2_zlos_from_deltas(dx12,dy12,dz12, r12_sq)
        μ13_sq = TriCo.mu2_zlos_from_deltas(dx13,dy13,dz13, r13_sq)
        μ23_sq = TriCo.mu2_zlos_from_deltas(dx23,dy23,dz23, r23_sq)
        (μ12_sq >= μmax2 || μ13_sq >= μmax2 || μ23_sq >= μmax2) && continue

        r12b = Int(floor((sqrt(r12_sq) - rmin) * invΔr)) + 1
        r13b = Int(floor((sqrt(r13_sq) - rmin) * invΔr)) + 1
        r23b = Int(floor((sqrt(r23_sq) - rmin) * invΔr)) + 1
        μ12b = Int(floor(sqrt(μ12_sq) * invΔμ)) + 1
        μ13b = Int(floor(sqrt(μ13_sq) * invΔμ)) + 1
        μ23b = Int(floor(sqrt(μ23_sq) * invΔμ)) + 1
        (1 <= r12b <= Nr && 1 <= r13b <= Nr && 1 <= r23b <= Nr &&
         1 <= μ12b <= Nμ && 1 <= μ13b <= Nμ && 1 <= μ23b <= Nμ) || continue

        r1, r2, r3 = r12b, r23b, r13b
        m1, m2, m3 = μ12b, μ23b, μ13b
        id1, id2, id3 = 1, 2, 3

        r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)
        r2, r3, m2, m3, id2, id3 = cswap!(r2, r3, m2, m3, id2, id3)
        r1, r2, m1, m2, id1, id2 = cswap!(r1, r2, m1, m2, id1, id2)

        H[r1, r2, m1, m2] += 1
    end
    H
end

end # module

