using Bessels
using Base.Threads

"""
    u, un = YukSLPeval(tx,tnx,sx,sw,dens,kappa; grad=true)
     
evaluate potential and its target-normal derivative for Yukawa (aka modified
H elmholtz) SLP at targs tx, given source density dens sampled at nodes sx 
with quadrature weights sw. Nodes and weights must be good for line integrals
on the curve.

Inputs:
    `tx` : 2*N targ coords
    `tnx` : 2*N targ normals (only used if grad true)
    `sx` : 2*N src node coords
    `sw` : N src node (speed) weights
    `dens` : N vec of densities
    `kappa` : modified Helmholtz param, ie, PDE is (Delta - kappa^2)u = 0
    `grad` : optional kwarg, whether to eval gradients.
Outputs:
    `u` : N vec of potentials at targs
    `un` : if grad=true, N vec of target normal derivs du/dn, else empty
"""
function YukSLPeval(tx,tnx,sx,sw,dens::Vector{T},kappa; grad=true) where T
    if kappa<=0
        return throw(DomainError(kappa, "`kappa` must be positive"))
    end
    Nt=size(tx,2)          # num targs
    Ns=size(sx,2)
    @assert length(sw)==Ns
    @assert length(dens)==Ns
    u = zeros(T,Nt)
    un::Vector{Float64} = zeros(T,0)    # un "captured var" => explicit type! 
    if grad
        un = zeros(T,Nt)    # the reallocation led to un being "captured" :(
        @assert size(tnx)==size(tx)
    end
    prefac = 1/(2*pi)
    @threads for i in eachindex(u)
        for j in eachindex(dens)
            d1 = tx[1,i]-sx[1,j]; d2 = tx[2,i]-sx[2,j] 
            r2 = d1^2 + d2^2
            r = sqrt(r2)
            pdw = prefac * dens[j] * sw[j]
            u[i] += pdw * besselk0(kappa*r)
            if grad                  # not much perf hit for conditional
                costh = (d1*tnx[1,i] + d2*tnx[2,i])/r
                un[i] += pdw * costh * kappa * besselk1(kappa*r)
            end
        end
    end
    u, un
end
# check K_0' = -K_1 :
# h=1e-5; x=0.7; (besselk0(x+h/2)-besselk0(x-h/2))/h + besselk1(x)

"""
    pot = YukSLP(tx,sx,sw,dens,kappa)
     
evaluate Yukawa (aka modified Helmholtz) SLP at targs tx, given source
density dens sampled at nodes sx with quadrature weights sw.
Nodes and weights must be good for line integrals on the curve.

Note: *** Legacy interface, to remove. ***

Inputs:
    `tx` : 2*N targ coords
    `sx` : 2*N src node coords
    `sw` : N src node (speed) weights
    `dens` : N vec of densities
    `kappa` : modified Helmholtz param, ie, PDE is (Delta - kappa^2)u = 0
Outputs:
    `pot` : N vec of potentials at targs
"""
YukSLP(tx,sx,sw,dens,kappa) = YukSLPeval(tx,[],sx,sw,dens,kappa,grad=false)[1]


"""
    A, An = YukSLPmats(tx,tnx,sx,sw,kappa; grad=true)
     
fill matrix mapping density values to potentials, and optionally to normal
gradients, for Yukawa (aka modified Helmholtz) SLP at targs tx
(optionally with target normals tnx), given source nodes sx with quadrature
weights sw. Nodes and weights must be good for line integrals on the curve.

Inputs:
    `tx` : 2*M targ coords
    `tnx` : 2*M targ normal coords (only used if grad=true)
    `sx` : 2*N src node coords
    `sw` : N src node (speed) weights
    `kappa` : modified Helmholtz param, ie, PDE is (Delta - kappa^2)u = 0
    `grad` : optional kwarg, whether to eval gradients.
Outputs:
    `A` : M*N potential-evaluation matrix
    `An` : M*N target-normal gradient evaluation matrix, if grad true, else []
"""
function YukSLPmats(tx,tnx,sx,sw,kappa; grad=true)
    T = eltype(tx)
    if kappa<=0
        return throw(DomainError(kappa, "`kappa` must be positive"))
    end
    Nt=size(tx,2)
    Ns=size(sx,2)
    @assert length(sw)==Ns
    A = Matrix{T}(undef,Nt,Ns)
    An::Matrix{Float64} = zeros(T,0,0)   # "captured var" => explicit type!
    if grad
        @assert size(tnx)==size(tx)
        An = Matrix{T}(undef,Nt,Ns)
    end
    prefac = 1/(2*pi)
    @threads for j in eachindex(sw)
        for i=1:Nt
            d1 = tx[1,i]-sx[1,j]; d2 = tx[2,i]-sx[2,j] 
            r2 = d1^2 + d2^2
            if r2==0.0
                A[i,j] = Inf
                if grad An[i,j] = Inf; end
            else
                r = sqrt(r2)
                A[i,j] = prefac*sw[j] * besselk0(kappa*r)
                if grad                  # not much perf hit for conditional
                    costh = (d1*tnx[1,i] + d2*tnx[2,i])/r
                    An[i,j] = prefac*sw[j] * costh * kappa*besselk1(kappa*r)
                end
            end
        end
    end
    A, An
end
# maybe have dens as optional arg, returns mat if absent?

"""
    YukSLPmat
     
matrix mapping density values to potentials, for Yukawa (aka modified
Helmholtz) SLP at targs tx, given source nodes sx with quadrature weights sw.
Nodes and weights must be good for line integrals on the curve.

Note: *** Legacy interface, to remove. ***

Inputs:
    `tx` : 2*N targ coords
    `sx` : 2*N src node coords
    `sw` : N src node (speed) weights
    `kappa` : modified Helmholtz param, ie, PDE is (Delta - kappa^2)u = 0
"""
YukSLPmat(tx,sx,sw,kappa) = YukSLPmats(tx,[],sx,sw,kappa,grad=false)[1]



"""
    YukSLPmat_selfcrude
     
Crude (1st-order) version of square special self-evaluation quadrature matrix
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
