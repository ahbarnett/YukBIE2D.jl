@info "Testing utils..."

using LinearAlgebra

@testset "unitcircle" begin
    x,w,nx,curv = unitcircle(100)
    @test sum(w) ≈ 2*pi rtol=1e-14
    @test norm(sum(w'.*nx,dims=2)) ≈ 0  atol=1e-14    # integral of nx vector
    @test sum(w.*curv) ≈ -2*pi  rtol=1e-13    # integral of kappa
end
@testset "starfish" begin
    x,w,nx,curv = starfish(300)
    @test norm(sum(w'.*nx,dims=2)) ≈ 0  atol=1e-14    # integral of nx vector
    @test sum(w.*curv) ≈ -2*pi  rtol=1e-10    # integral of kappa
    _,w2,_,_ = starfish(350)
    @test sum(w) ≈ sum(w2) rtol=1e-13   # perim convergence, is spectral
end

@testset "spectral" begin
    f(x) = sin(1+3*x)        # pick analytic func
    fp(x) = 3*cos(1+3*x)
    @testset "clencurt" begin
        Iex = f(1.0)-f(-1.0)
        x,w = clencurt(20)        
        I = sum(w.*fp.(x))
        @debug "clencurt err (even)" I-Iex
        @test I ≈ Iex rtol=1e-14
        x,w = clencurt(21)        
        I = sum(w.*fp.(x))
        @debug "clencurt err (odd)" I-Iex
        @test I ≈ Iex rtol=1e-14
    end
    @testset "chebydiffmat" begin
        D,x = chebydiffmat(20)
        @test norm(D*f.(x)-fp.(x)) < 1e-13      # actually 1e-14 rel to size 3
    end
end
