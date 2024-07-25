# simple solve of modified Helmholtz PDE with Dirichlet BC on an arc, using
# SLP on the arc. Barnett 7/25/24
using YukBIE2D
using LinearAlgebra
using Printf
using Gnuplot
using ColorSchemes  # for gnuplot heatmaps

verb = 0
ka = 0.01      # aka phi, inverse decay length
# arc param on [-1,1]
a = 0.5; b=0.1          # angular half-width, angular offset
seg(t) = [-1.0+cos(a*t+b),sin(a*t+b)]   # fit inside unit circle
segp(t) = [-a*sin(a*t+b),a*cos(a*t+b)]

function arcquad(seg,segp,N)   # discretize arc param by seg and deriv segp
    t,w = clencurt(N-1)       # N nodes
    sx = reduce(hcat, seg.(t))  # 2*N
    sxp = reduce(hcat, segp.(t))  # 2*N
    speed = vec(sqrt.(sum(sxp.^2,dims=1)))
    sw = w.*speed
    sx,sw,t
end

function testsegp(seg,segp,N)      # check funcs segp and seg
    t,w = clencurt(N-1)       # N nodes
    sx = reduce(hcat, seg.(t))  # 2*N
    sxp = reduce(hcat, segp.(t))  # 2*N
    D,_ = chebydiffmat(N-1)   # N nodes
    sxpa = sx*D'          # deriv of rows
    @printf "deriv max err on coords %.3g\n" norm(vec(sxp-sxpa),Inf)
end
if verb>0 testsegp(100); end

f(t) = 1.0    # RHS Dirichlet data on arc

xtest = [-0.4,0.3]   # test pt
ng = 300           # plot grid size (N=200 ng=30 takes 0.04s to eval :)
g = range(-1,1,ng)       # grid in each dim
o = ones(size(g))           # also row-vec
tx = [kron(o,g)';kron(g,o)']  # fill grid of targs (ok to fill, sim size to u)

doplainbvp = (verb>0)
if doplainbvp
Ns=20:20:100
for (k,N) in enumerate(Ns)    # ...... conv
    sx,sw,t = arcquad(seg,segp,N)
    A = YukSLPmat_selfcrude(sx,sw,ka)
    rhs = f.(t)
    dens = A\rhs
    utest,_ = YukSLPeval(xtest,[],sx,sw,dens,ka,grad=false)
    @printf "N=%d\tutest=%.12g\n" N utest[1]
    if N==Ns[end]
        u,_ = YukSLPeval(tx,[],sx,sw,dens,ka,grad=false)   # grid eval
        u = reshape(u,(ng,ng))
        @gp g g u "w image notit" palette(:jet1) xlab="x" ylab="y"
        @gp :- sx[1,:] sx[2,:] "w p" "set size ratio -1"
        @gp :- [xtest[1]] [xtest[2]] "w p pt 7"
        Gnuplot.save("pics/openarc_bvp_u.png",term="pngcairo")
        @gp dens "w lp t 'dens vs j'"
        Gnuplot.save("pics/openarc_bvp_dens.png",term="pngcairo")
    end
end                            # .......
end

# ========== add in Neumann (zero-flux) outer bdry (circle)
Ns=20:20:100               # arc nodes (-1)
for (k,N) in enumerate(Ns)    # ...... conv
    sx,sw,t = arcquad(seg,segp,N)
    Nb = 3*N                   # since bdry longer than seg
    bx,bw,bnx,bcurv = unitcircle(Nb)     # discretize bdry
    # setup 2x2 lin sys:       A [dens1; dens2] = [f; 0]
    # where dens1 is on arc, dens2 on outer bdry
    S11 = YukSLPmat_selfcrude(sx,sw,ka)     # arc to self
    S12,_ = YukSLPmats(sx,[],bx,bw,ka,grad=false)    # arc from bdry
    _,DT21 = YukSLPmats(bx,bnx,sx,sw,ka)
    DT22 = YukSLPdermat_selfcrude(bx,bnx,bcurv,bw,ka)
    A = [S11 S12; DT21 0.5*I+DT22]
    #println(svdvals(0.5*I + DT22))   # bdry only, int Neu BVP -> nullity 1
    rhs = [f.(t); zeros(Nb)]     # stack the RHS vectors
    dens = A\rhs              # solve
    dens1,dens2 = dens[1:N],dens[N+1:end]    # extract density vecs
    #show(varinfo())
    utest,_ = YukSLPeval(xtest,[],sx,sw,dens1,ka,grad=false) .+ YukSLPeval(xtest,[],bx,bw,dens2,ka,grad=false)
    @printf "N=%d\tutest=%.12g\n" N utest[1]
    if N==Ns[end]
        u,_ = YukSLPeval(tx,[],sx,sw,dens1,ka,grad=false) .+ YukSLPeval(tx,[],bx,bw,dens2,ka,grad=false)
        u = reshape(u,(ng,ng))
        @gp g g u "w image notit" palette(:jet1) xlab="x" ylab="y"
        @gp :- sx[1,:] sx[2,:] "w p ps 0.3 pt 7 lc '#000000'" "set size ratio -1"
        @gp :- bx[1,:] bx[2,:] "w p ps 0.3 pt 7 lc '#000000'"
        @gp :- [xtest[1]] [xtest[2]] "w p"
        Gnuplot.save("pics/openarc_noflux_u.png",term="pngcairo")
        @gp dens "w lp t 'dens vs j'"
        Gnuplot.save("pics/openarc_noflux_dens.png",term="pngcairo")
    end
end                            # .......


