@info "Testing utils..."

using LinearAlgebra

@testset "unitcircle" begin
    x,w,nx,kap = unitcircle(100)
    @test sum(w) ≈ 2*pi rtol=1e-14
end
@testset "starfish" begin
    x,w,nx,kap = starfish(300)
    @test sum(w.*nx) ≈ 0  atol=1e-13    # integral of nx vector
    x2,w2,nx2,kap2 = starfish(400)
    @test sum(w) ≈ sum(w2) rtol=1e-14   # perim convergence, is spectral
end

@testset "spectral" begin
    f(x) = sin(1+3*x)        # pick analytic func
    fp(x) = 3*cos(1+3*x)
    @testset "clencurt" begin
        x,w = clencurt(20)        
        Iex = f(1.0)-f(-1.0)
        I = sum(w.*fp.(x))
        @debug "clencurt err" I-Iex
        @test I ≈ Iex rtol=1e-14
    end
    @testset "chebydiffmat" begin
        D,x = chebydiffmat(20)
        @test norm(D*f.(x)-fp.(x)) < 1e-13      # actually 1e-14 rel to size 3
    end
end
