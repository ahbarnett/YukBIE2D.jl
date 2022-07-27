using Bessels
using LinearAlgebra

ka = 0.8       # kappa = phi = Yukawa inverse length
r = 1.5; th0 = 1.2
tx = [r*cos(th0),r*sin(th0)]   # targ loc
N = 80
sx,sw = unitcircle(N)

@testset "YukSLP and YukSLPmat far eval unit dens on unit circle" begin
    dens = ones(N)
    u = YukSLP(tx,sx,sw,dens,ka)
    # solve m=0 Fourier mode case by jump matching a.I0(ka.r) r<1, b.K0(ka.r) r>1:
    bcoeff = 1/(ka*(besselk(1,ka)+besselk(0,ka)*besseli(1,ka)/besseli(0,ka)))
    uex = bcoeff*besselk(0,r*ka)
    #println("check SLP far from unit dens, unit circle...")
    @test uex ≈ u[1] rtol=1e-15
    A = YukSLPmat(tx,sx,sw,ka)
    @test norm(A*dens-u)/norm(u) < 1e-15
end
