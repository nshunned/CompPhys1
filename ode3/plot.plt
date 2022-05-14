# macro plt.plt

plot "terminal.dat" u 1:2 w lp lt 7
set xlabel "mass (kg)"
set ylabel "terminal velocity (m/s)"
set title "Terminal Velocity vs. Mass"

load "format.plt"

set term png
set output "vt_vs_m.png"
replot
unset term