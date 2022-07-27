# random Alex utils

using LinearAlgebra    # needed for I

function di(x)
    """Print condensed form of a variable (REPL style) to screen"""
    io=IOBuffer();
    show(IOContext(io, :limit => true, :displaysize => (24,80)), "text/plain",x);
    println(String(take!(io)))  # yuk
end


"""
interpmat1d(t,s) interpolation matrix from source nodes s to target nodes t

Only stable for small node numbers (<100).
"""
function interpmat1d(t,s)



end


"""
clencurt(n::Integer) -> nodes, weights
for (n+1)-point Clenshaw-Curtis quadrature on [-1,1].
Taken from Toby Driscoll's RNC book. n must be even.
"""
function clencurt(n)
    @assert iseven(n) "Value of `n` must be an even integer."
    # set x as Chebyshev extreme nodes...
    θ = [ i*π/n for i in 0:n ]
    x = -cos.(θ)
    # Compute the C-C weights c...
    c = similar(θ)
    c[[1,n+1]] .= 1/(n^2-1)        # end weights
    s = sum( cos.(2k*θ[2:n])/(4k^2-1) for k in 1:n/2-1 )
    v = @. 1 - 2s - cos(n*θ[2:n])/(n^2-1)
    c[2:n] = 2v/n
    x, c
end


"""
    D,x = chebydiffmat(N) returns dense (N+1)*(N+1) differentiation matrix 
    `D` and N+1 Chebychev nodes `x`, for the standard 1D interval [-1,1].
    The matrix multiplies a vector of function values
    at these nodes to give an approximation to the vector of derivative values.
    Note, nodes are in descending order starting at 1.0 and ending at -1.0.
    Adapted from L N Trefethen's cheb.m from "Spectral Methods in MATLAB" book.
"""
function chebydiffmat(N)
    if N==0
        x=1; D=0
    else
        x = cos.(π*(0:N)/N)
        c = [2; ones(N-1); 2] .* (-1).^(0:N)
        X = x*ones(N+1)'              # duplicates nodes in each col
        dX = X-X'                     # matrix of pairwise node differences
        D = (c*(1.0./c)') ./ (dX + I)
        D = D - diagm(vec(sum(D,dims=2)))
    end
    return D,x
end



"""
unitcircle(N::Integer)
    
Periodic trap rule quadrature for unit circle. Returns 2*N node coords,
N weights.
The first node is at [1;0]
"""
function unitcircle(N::Integer)
    th = (0:N-1)/N*2*pi
    x = [cos.(th)'; sin.(th)']
    w = (2*pi/N)*ones(N)
    x,w
end

"""
starfish(N::Integer=100,freq::Integer=5,ampl=0.3,rot=1.0) -> x,w
    
Periodic trap rule quadrature for smooth starfish. Returns 2*N node coords,
N weights. The first node is at theta=0.

*** to doc, output struct instead, etc.
"""
function starfish(N::Integer=100,freq::Integer=5,ampl=0.3,rot=1.0)
    th = (0:N-1)'/N*2*pi                      # col vec
    r = 1.0 .+ ampl*cos.(freq*th .- rot)   # set of radii. Note .operators
    rp = -ampl*freq*sin.(freq*th .- rot)
    c = cos.(th); s = sin.(th)
    x = [r.*c; r.*s]
    xp = [rp.*c.-r.*s; rp.*s.+r.*c]  # velocity
    sp = sqrt.(sum(xp.^2,dims=1))        # speed at nodes
    w = (2*pi/N)*sp                   # PTR
    x,w
end
