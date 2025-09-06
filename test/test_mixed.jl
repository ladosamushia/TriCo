using Test
using TriCo

@testset "Mixed triangles basics" begin
    Xa,Ya,Za = randn(20), randn(20), randn(20)
    Xb,Yb,Zb = randn(20), randn(20), randn(20)
    Xc,Yc,Zc = randn(20), randn(20), randn(20)

    A = TriCo.TriCat(Xa,Ya,Za)
    B = TriCo.TriCat(Xb,Yb,Zb)
    C = TriCo.TriCat(Xc,Yc,Zc)

    H1 = count_triangles_mixed!(A,B,C; rmin=0.5, rmax=3.0, Nr=5, μmax=0.9, Nμ=2)
    @test sum(H1.h) > 0

    # AAA should reduce to the existing single-catalog API
    H2 = count_triangles_grid!(Xa,Ya,Za; rmin=0.5, rmax=3.0, Nr=5, μmax=0.9, Nμ=2)
    H3 = count_triangles_mixed!(A,A,A; rmin=0.5, rmax=3.0, Nr=5, μmax=0.9, Nμ=2)
    @test H2.h == H3.h
end

