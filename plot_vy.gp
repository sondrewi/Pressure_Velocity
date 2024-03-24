set terminal pngcairo enhanced
set output 'Sols_Re100.png'

set pm3d map

set palette defined ( \
    -0.5 'dark-blue', \
    -0.3333 'blue', \
    -0.1666 'light-blue', \
    0 'white', \
    0.1 'light-red', \
    0.2 'red', \
    0.3 'dark-red' \
)

set cbrange [-0.5:0.3]
set title "Y-velocity at Re = 100, 200 by 200 cells"
set size square
set xtics rotate by -45

splot 'vy_Re2.dat' using 1:2:4 with image notitle