
set terminal pngcairo enhanced
set output 'Sols_Re100.png'

set pm3d map

set palette defined ( \
    -0.2 'dark-blue', \
    -0.133333 'blue', \
    -0.066666 'light-blue', \
    0 'white', \
    0.333333 'light-red', \
    0.666666 'red', \
    1 'dark-red' \
)

set cbrange [-0.2:1]
set title "X-velocity at Re = 100, 200 by 200 cells
set size square
set xtics rotate by -45


splot 'vx_Re100.dat' using 1:2:3 with image notitle