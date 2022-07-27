
@testset "unitcircle" begin
    x,w = unitcircle(100)
    @test sum(w) ≈ 2*pi rtol=1e-14
end
@testset "starfish" begin
    x,w = starfish(300)
    x2,w2 = starfish(400)
    @test sum(w) ≈ sum(w2) rtol=1e-14   # perim convergence spectral
end
