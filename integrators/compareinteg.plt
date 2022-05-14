set xlabel "# Trials"
set ylabel "Relative error (estimate) per trial"
set logscale x
set logscale y
set format y "%1.2e"


set title "Error analysis per trial for integral of exp(-x)"
set label "points: delta trial, lines: true err" at 2,.2

if ( system("[ ! -e gaus2.dat  ]; echo $?")) { \
plot [][1e-17:] "trap.dat" lc 1 title "Trap rule", \
"trap.dat" using 1:3 w line lc 1 ti "", \
"simp.dat" lc 2 title "Simpson", \
"simp.dat" using 1:3 w line lc 2 ti "", \
"gaus.dat" lc 4 pt 7 title "Gauss (from text)", \
"gaus.dat" using 1:3 w line lc 4 ti "", \
"gaus2.dat" lc 6 title "Gauss (improved)", \
"gaus2.dat" using 1:3 w line lc 6 ti "" }
else { \
plot "trap.dat" title "trap rule",\
"simp.dat" title "Simpson",\
"gaus.dat" title "Gauss" }


# make a pdf file
set term pdf
set output "compareinteg.pdf"
replot

# make a png file
set term png
set output "compareinteg.png"
replot


# reset gnuplot state to continue using interactively
set output
set term qt
