set terminal pngcairo enhanced
set output 'Res_plot.png'

set title "Inf-Norm of Residuals in Inner Loops"
set xlabel "Iteration"
set ylabel "Residual Inf Norm"

set xtics add (2 2)
set xrange [2:*]
set format y "10^{%L}"
set format x "10^{%L}"
set logscale xy
set yrange[0.00000002: 0.00000005]

set ytics add ('2x10^{-8}' 2e-8, '3x10^{-8}' 3e-8, '4x10^{-8}' 4e-8, '5x10^{-8}' 5e-8)
set xtics ("2" 2, "10^{1}" 10, "10^{2}" 100, "10^{3}" 1000, "10^{4}" 10000)

set style line 1 linetype 1 linewidth 2 linecolor 'blue' pointtype 7 pointsize 1.5 dashtype (10,5)
set style line 2 linetype 1 linewidth 2 linecolor 'red'

plot 'Norms_Re100.dat' every ::1 using 0:3 with lines linestyle 1 title 'p (Re = 100)', \
     'Norms_Re1000.dat' every ::1 using 0:3 with lines linestyle 2 title 'p (Re = 1000)'