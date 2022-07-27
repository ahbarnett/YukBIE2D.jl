using Bessels
using Base.Threads

"""
    YukSLP
     
evaluate Yukawa (aka modified Helmholtz) SLP at targs tx, given source
density dens sampled at nodes sx with quadrature weights sw.
Nodes and weights must be good for line integrals on the curve.

Inputs:
    `tx` : 2*N targ coords
    `sx` : 2*N src node coords
    `sw` : N src node (speed) weights
    `dens` : N vec of densities
    `kappa` : modified Helmholtz param, ie, PDE is (Delta - kappa^2)u = 0
"""
function YukSLP(tx,sx,sw,dens,kappa)
    T = eltype(dens)
    if kappa<=0
        return throw(DomainError(kappa, "`kappa` must be positive"))
    end
    Nt=size(tx,2)          # num targs
    Ns=size(sx,2)
    @assert length(sw)==Ns
    pot = zeros(T,Nt)
    prefac = 1/(2*pi)
    @threads for i in eachindex(pot)       # @tturbo LoopVec segfaults!
        for j in eachindex(dens)
            r2 = (tx[1,i]-sx[1,j])^2 + (tx[2,i]-sx[2,j])^2
            pot[i] += prefac * dens[j] * sw[j] * besselk0(kappa*sqrt(r2))
        end
    end
    pot
end

"""
    YukSLPmat
     
matrix mapping density values to potentials, for Yukawa (aka modified
Helmholtz) SLP at targs tx, given source nodes sx with quadrature weights sw.
Nodes and weights must be good for line integrals on the curve.

Inputs:
    `tx` : 2*N targ coords
    `sx` : 2*N src node coords
    `sw` : N src node (speed) weights
    `kappa` : modified Helmholtz param, ie, PDE is (Delta - kappa^2)u = 0
"""
function YukSLPmat(tx,sx,sw,kappa)
    T = eltype(tx)
    if kappa<=0
        return throw(DomainError(kappa, "`kappa` must be positive"))
    end
    Nt=size(tx,2)
    Ns=size(sx,2)
    @assert length(sw)==Ns
    A = zeros(T,Nt,Ns)
    prefac = 1/(2*pi)
    @threads for j in eachindex(sw)
        for i=1:Nt
            r2 = (tx[1,i]-sx[1,j])^2 + (tx[2,i]-sx[2,j])^2
            if r2==0.0
                A[i,j] = Inf
            else
                A[i,j] = prefac * sw[j] * besselk0(kappa*sqrt(r2))
            end
        end
    end
    A
end
# maybe have dens as optional arg, returns mat if absent?

"""
    YukSLPmat_selfcrude
     
Crude (1st-order) version of square special sel-evaluation quadrature matrix
mapping density values to potentials at nodes on a curve, for Yukawa
(aka modified Helmholtz) S operator.
Nodes sx and weights sw must be good for line integrals on the curve.
This works equally well for closed curve or open arc.
Needs h (node spacing) << 1/kappa (Yukawa lengthscale).
Returns N*N matrix.

Inputs:
    `sx` : 2*N src node coords
    `sw` : N src node (speed) weights
    `kappa` : modified Helmholtz param, ie, PDE is (Delta - kappa^2)u = 0
"""
function YukSLPmat_selfcrude(sx,sw,kappa)
    # first fill self-eval mat w/ plain quadr rule (diag = Inf)
    A = YukSLPmat(sx,sx,sw,kappa)
    # override just diag entries
    for j in eachindex(sw)
        # analytic formula for const-dens, const-speed, straight-segment
        A[j,j] = -sw[j]/(2*pi) * (log(kappa*sw[j]/2)-1)
    end
    A
end
