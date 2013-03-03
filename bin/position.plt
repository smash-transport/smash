set xlabel 'x1'
set ylabel 'x2'
set zlabel 'x3'
splot 'data/position.dat' using 2:3:4 title '3d IC positions' with points lc rgb '#8b1a0e'
pause -1
