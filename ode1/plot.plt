# macro plt.plt

# -------------------------------------------------------
# landing x vs. dt

plot "landing.dat" u 1:2 w lp lt 7 t "x_f"
set xlabel "dt"
set ylabel "x_f"
set title "Calculated Landing Range vs. Time Step"

load "format.plt"

set term png
set output "x_f_vs_dt.png"
replot
unset term

# -------------------------------------------------------
# landing y vs. dt

plot "landing.dat" u 1:3 w lp lt 7 t "y_f"
set xlabel "dt"
set ylabel "y_f"
set title "Calculated Landing Height vs. Time Step"

load "format.plt"

set term png
set output "y_f_vs_dt.png"
replot
unset term

# -------------------------------------------------------
# landing x error vs. dt

plot "landing.dat" u 1:4 w lp lt 7 t "x_{f,err}"
set xlabel "dt"
set ylabel "x_{f,err}"
set title "Error in Calculated Range vs. Time Step"

load "format.plt"

set term png
set output "x_f_err_vs_dt.png"
replot
unset term

# -------------------------------------------------------
# energies vs. t

plot "euler.dat" u 1:6 w l lw 2.5 lc "red" t "K(t)", \
     "euler.dat" u 1:7 w l lw 2.5 lc "green" t "U(t)", \
     "euler.dat" u 1:8 w l lw 2.5 lc "blue" t "E(t)"
set xlabel "t"
set ylabel "energy"
set title "Calculated Energies vs. Time"

load "format.plt"

set term png
set output "e_vs_t.png"
replot
unset term