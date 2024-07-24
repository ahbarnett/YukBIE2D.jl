# basic direct summation potential eval timer

using YukBIE2D
using Printf
using LinearAlgebra

ka = 1.6
N = 100
sx,sw = unitcircle(N)
Nt = Int(1e6)   # big enough to time single call (checked gives same as btime)
tx=rand(2,Nt)
tnx = rand(2,Nt)
dens = rand(N)
t = @elapsed u = YukSLP(tx,sx,sw,dens,ka)
@printf "YukSLP     \t%.3g s\t\t%.3g Gpair/s\n" t N*Nt/t/1e9
t = @elapsed _,un = YukSLPeval(tx,tnx,sx,sw,dens,ka)
@printf "YukSLPeval \t%.3g s\t\t%.3g Gpair/s\n" t N*Nt/t/1e9
# debug perf: @code_warntype YukSLPeval(tx,tnx,sx,sw,dens,ka);
t = @elapsed A = YukSLPmat(tx,sx,sw,ka)
@printf "YukSLPmat  \t%.3g s\t\t%.3g Gpair/s\n" t N*Nt/t/1e9
t = @elapsed _,An = YukSLPmats(tx,tnx,sx,sw,ka)
@printf "YukSLPmats \t%.3g s\t\t%.3g Gpair/s\n" t N*Nt/t/1e9
@printf "mat vs pot err  %.3g\n" norm(A*dens-u)
@printf "mat vs grad err %.3g\n" norm(An*dens-un)
