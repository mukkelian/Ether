set nokey
set terminal png font "Times-Roman,10"
set output 'fig1.png'
set multiplot                      
set size squar 0.6,0.52
set nokey
set format y "%g"
set xtics 0,20,200
set mxtics 10
set origin 0.0, 0.0
set grid
set xlabel "Temp. (K)"
set ylabel "Magnetization"
plot 'magnetization.dat' u 1:2 w lp lt 6 notitle

set origin 0.5, 0.0
set xlabel "Temp. (K)"
set ylabel "Susceptibility"
plot 'magnetization.dat' u 1:3 w lp lt 6 notitle

set origin 0.0,0.5
set xlabel "Temp. (K)"
set ylabel "Energy"
plot 'energy.dat' u 1:2 w lp lt 6 notitle

set origin 0.5,0.5
set xlabel "Temp. (K)"
set ylabel "Specific heat"
plot 'energy.dat' u 1:3 w lp lt 6 notitle
unset multiplot
