# test convergence of Yukawa SLP self-evaluation schemes, eg crude O(h).
# plots not a unit-test.
# Barnett 7/25/22

using YukBIE2D
using Bessels
using Printf
using Gnuplot
using ColorSchemes
using LinearAlgebra

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

println("self-conv test wrt N, to unknown on closed starfish...")
Ns = 100:100:500
for j in eachindex(Ns)
    N = Ns[j]
    sx,sw = starfish(N,5,0.3)
    A = YukSLPmat_selfcrude(sx,sw,ka)
    dens = 0.6 .+ sx[1,:]       # dens = x-coord
    u = A*dens
    @printf "N=%d:\tperim=%.12g\tusurf=%.6g\n" N sum(sw) u[1]  # the fixed node!
end

# add 2D potential image plot corresp to on-surf convergence...
N=200
sx,sw = starfish(N)
dens = 0.6 .+ sx[1,:]       # dens = x-coord
ng = 300           # plot grid size (N=200 ng=30 takes 0.04s to eval :)
g = range(-2,2,ng)       # grid in each dim
o = ones(size(g))           # also row-vec
tx = [kron(o,g)';kron(g,o)']  # fill grid of targs (ok to fill, sim size to u)
# also consider LazyGrids.jl to save RAM
u = YukSLP(tx,sx,sw,dens,ka)
u = reshape(u,(ng,ng))
@gp :pot g g u "w image notit" "set size square" palette(:jet1) xlab="x" ylab="y"
@gp :pot :- sx[1,:] sx[2,:] "w lp pt 7 ps 0.3 lc '#000000'"
# to kill the window...
#Gnuplot.quit(:pot)
# to kill all windows...
#for i in session_names(); Gnuplot.quit(i); end

# convergence of solving 1st-kind lin sys for Yukawa Dirichlet BVP u=f
xfar = [1.5;1.0]    # target
Ns = 100:100:500
for j in eachindex(Ns)
    N = Ns[j]
    sx,sw = starfish(N,5,0.3)
    f = 0.5 .+ sx[2,:]      # RHS func
    A = YukSLPmat_selfcrude(sx,sw,ka)
    dens = A\f
    ufar = YukSLP(xfar,sx,sw,dens,ka)
    # exploit that 1st (=0th) node fixed in space -> meaningful convergence...
    @printf "N=%d:\tcond=%.6g\tdens[1]=%.6g\tufar=%.6g\n" N cond(A) dens[1] ufar[1]
    if j==3
        u = YukSLP(tx,sx,sw,dens,ka)
        u = reshape(u,(ng,ng))
        @gp :ufar g g u "w image notit" "set size square" palette(:jet1) xlab="x" ylab="y"
        @gp :ufar :- sx[1,:] sx[2,:] "w lp pt 7 ps 0.3 lc '#000000'"
    end
end
# we see cond# grow like N

