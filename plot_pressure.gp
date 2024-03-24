set terminal pngcairo enhanced

set output 'Sols_Re100.png'

set pm3d map

set palette defined ( \
    -14 'blue', \
    0 'white', \
    20 'red' \
)
set xtics rotate by -45

set cbrange [-14:20]
set title "Pressure for Re = 100, 200 by 200 cells"
set size square

splot 'pressure_Re100.dat' using 1:2:5 with image notitle
