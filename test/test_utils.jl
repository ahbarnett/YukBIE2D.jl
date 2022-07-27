
@info "Testing utils..."

@testset "unitcircle" begin
    x,w = unitcircle(100)
    @test sum(w) ≈ 2*pi rtol=1e-14
end
@testset "starfish" begin
    x,w = starfish(300)
    x2,w2 = starfish(400)
    @test sum(w) ≈ sum(w2) rtol=1e-14   # perim convergence, is spectral
end
@testset "clencurt" begin
    x,w = clencurt(20)
    f(x) = sin(1+3*x)
    fp(x) = 3*cos(1+3*x)
    Iex = f(1.0)-f(-1.0)
    I = sum(w.*fp.(x))
    @debug "clencurt err" I-Iex
    @test I ≈ Iex rtol=1e-14
end
