using Bessels
using LinearAlgebra

@info "testing layerpot (YukSLP and related)..."

ka = 0.8       # kappa = phi = Yukawa inverse length
r = 1.5; th0 = 1.2
tx = [r*cos(th0),r*sin(th0)]   # targ loc
N = 80
sx,sw = unitcircle(N)

# Notes about test framework:
# * is the testset name supposed to be long description or a keyword?
# * each @testset has its own local scope (not good for keeping setups)
# * verbose=true only prints nested testset names, *not* individual tests.
@testset "YukSLP" verbose=true begin
    dens = ones(N)
    u = YukSLP(tx,sx,sw,dens,ka)
    # solve m=0 Fourier mode case by jump matching a.I0(ka.r) r<1, b.K0(ka.r) r>1:
    bcoeff = 1/(ka*(besselk(1,ka)+besselk(0,ka)*besseli(1,ka)/besseli(0,ka)))
    uex = bcoeff*besselk(0,r*ka)
    abserr = uex-u[1]
    @test abserr<1e-15             # quadrature error
    # try Julia's logger as good way for feedback during tests?
    @info "SLP vs analytic far from unit dens on unit circle:" abserr
    A = YukSLPmat(tx,sx,sw,ka)
    @test norm(A*dens-u)/norm(u) < 1e-15   # is pure rounding error
end
