# benchmarking K_0 implementations available in Julia. Barnett 7/25/22

using BenchmarkTools
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1

using SpecialFunctions
using BesselK
using Bessels

println("single-thread real-arg K_0 Bessel tests...")
n = Int(1e4); x = 10*rand(Float64,n)         # small arrays stay in cache?
println("warm up besselk"); y = besselk.(0,x);
println("warm up adbesselk"); yad = adbesselk.(0,x);   # NaN for x>8.5 !
println("warm up besselk0"); y0 = besselk0.(x);
@printf "l2 rel diff ad vs SF: %.3g\n" norm(y-yad)/norm(y)
@printf "l2 rel diff heltonmc vs SF: %.3g\n" norm(y-y0)/norm(y)
#x[.!isfinite.(yad)] # check args giving NaNs

tobj = @benchmark besselk.(0,x)
@printf "SF besselk    K_0 time: %.3g ns\n" median(tobj.times)/n
tobj = @benchmark adbesselk.(0,x)
@printf "adbesselk     K_0 time: %.3g ns\n" median(tobj.times)/n
tobj = @benchmark besselk0.(x)
@printf "heltonmc besselk0 time: %.3g ns\n" median(tobj.times)/n

# results: SF ~ 280ns, ad ~ 200ns, heltonmc ~ 15ns

# singlethread SF.jl: 0.005 Gpair/s ~ 200ns per K0 eval!
# SpecialFunctions.jl calls (lib)openspecfun which wraps AMOS (Donald E Amos)
# Fortran codes, namely zbesk (which is complex args, not needed here!)
# This destroys the advantage of native Julia compilation (SIMD, etc)

# And BesselK.jl is unusable due to NaNs and no faster than SF.

# cf MATLAB R2021b: n=1e7; x=10*rand(n,1); tic; y=besselk(0,x); t=toc; t/n*1e9
# gives: 34 ns (16 threads, 8 cores), about same per core
# Octave v6.1.1: 300 ns (single-threaded)
