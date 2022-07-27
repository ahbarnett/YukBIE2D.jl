# test convergence of Yukawa SLP self-evaluation schemes, eg crude O(h).
# plots not a unit-test.
# Barnett 7/25/22

using YukBIE2D
using Bessels
using Printf
using Gnuplot

ka = 0.8       # kappa = phi = Yukawa inverse length
# solve m=0 Fourier mode case by jump matching a.I0(ka.r) r<1, b.K0(ka.r) r>1:
bcoeff = 1/(ka*(besselk(1,ka)+besselk(0,ka)*besseli(1,ka)/besseli(0,ka)))
uex = bcoeff*besselk(0,ka)  # exact on-surf SLP for dens=1

println("self-mat conv test wrt N, to known dens=1 on unit circle...")
Ns = 50:50:300
errs = zeros(size(Ns))
for j in eachindex(Ns)
    N = Ns[j]
    sx,sw = unitcircle(N)
    A = YukSLPmat_selfcrude(sx,sw,ka)
    dens = ones(N)
    u = A*dens
    errs[j] = abs(u[1]-uex)
    @printf "N=%d:\terr=%.3g\n" N errs[j]
end
@gp Ns errs "w lp t 'on-surf err'"
@gp :- "set logscale x" "set logscale y" xlabel="N"
@gp :- Ns 0.2*Ns.^(-1) "w l dashtype 2 t 'O(h)'"  # verify 1st-order
