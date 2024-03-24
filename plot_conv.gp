set terminal pngcairo enhanced
set output 'conv_plot_100.png'

set title "Inf-Norm of Greatest Absolute Change Inner and Outer"
set xlabel "Iteration"
set ylabel "|∆|_{∞}"

set logscale xy
set format y "10^{%L}"
set format x "10^{%L}"
set xtics add (2 2)

set xtics ("2" 2, "10^{1}" 10, "10^{2}" 100, "10^{3}" 1000)
set xrange [2:*]

set style line 1 linetype 1 linewidth 2 linecolor 'blue' pointtype 7 pointsize 1.5 dashtype (10,5)
set style line 2 linetype 1 linewidth 2 linecolor 'red' pointtype 7 pointsize 1.5 dashtype (10,5)
set style line 3 linetype 1 linewidth 2 linecolor 'green' pointtype 7 pointsize 1.5 dashtype (10,5)
set style line 4 linetype 1 linewidth 2 linecolor 'purple'
set style line 5 linetype 1 linewidth 2 linecolor 'orange'
set style line 6 linetype 1 linewidth 2 linecolor 'cyan'

plot 'Norms_Re100.dat' every ::1 using 0:4 with lines linestyle 1 title 'vx(Out)', \
     'Norms_Re100.dat' every ::1 using 0:5 with lines linestyle 2 title 'vy(Out)', \
     'Norms_Re100.dat' every ::1 using 0:6 with lines linestyle 3 title 'p(Out)', \
     'Norms_Re100.dat' every ::1 using 0:7 with lines linestyle 4 title 'vx(In)', \
     'Norms_Re100.dat' every ::1 using 0:8 with lines linestyle 5 title 'vy(In)', \
     'Norms_Re100.dat' every ::1 using 0:9 with lines linestyle 6 title 'p(In)'