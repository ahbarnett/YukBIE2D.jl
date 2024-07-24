# getting D^T Nystrom matrix right. 7/24/24
using YukBIE2D
using Gnuplot

ka = 0.4
N = 100
x,w,nx,curv = unitcircle(N)
DT = YukSLPdermat_selfcrude(x,nx,curv,w,ka)
@gp DT "w image notit" "set size square" palette(:jet1) xlab="i" ylab="j"
Gnuplot.save("pics/DT.png", term="pngcairo")   # verifies diag value by eye
