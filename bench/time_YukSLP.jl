# basic direct summation potential eval timer

using YukBIE2D
using Printf
using LinearAlgebra

ka = 1.6
N = 100
sx,sw = unitcircle(N)
Nt = Int(1e6)   # big enough to time single call (checked gives same as btime)
tx=rand(2,Nt)
dens = rand(N)
t = @elapsed u = YukSLP(tx,sx,sw,dens,ka)
@printf "YukSLP time  %.3g s\t%.3g Gpair/s\n" t N*Nt/t/1e9
t = @elapsed A = YukSLPmat(tx,sx,sw,ka)
@printf "YukSLPmat time %.3g s\t%.3g Gpair/s\n" t N*Nt/t/1e9
@printf "mat vs eval err %.3g\n" norm(A*dens-u)
