set xlabel "t"
set ylabel "scattering rate"
set yrange [7:9]
# emulate error bands
set bars small 
# prettier default colors
set style line 1 lt rgb "#A00000" lw 2 pt 1
set style line 2 lt rgb "#00A000" lw 2 pt 2
set style line 3 lt rgb "#5060D0" lw 2 pt 0
set style line 4 lt rgb "#F25900" lw 4 pt 3
set style line 5 lt rgb "#ffabab" lw 1 pt 0
set style line 6 lt rgb "#aaffaa" lw 1 pt 0
set style line 7 lt rgb "#99a5ff" lw 1 pt 0
set ytics 0.25
set terminal pdfcairo font "Gill Sans,12" size 15cm,12cm
set output "sigma_volume_06.pdf"
set title 'pions T = 0.6 [GeV], sigma = 10 [mbarn], eps = 0.001, 50 runs'
plot "s_mash_a_10_T_6_s_10_t_1_avg.dat" using 1:2:3 with yerrorb ls 5 notitle, \
     "s_mash_a_10_T_6_s_10_t_1_avg.dat" title "V = 10 [fm^3]"  with lines ls 1, \
     8.43548 title "analytic estimate" with line ls 4
#     "/tmp/sj_a_10_eps_001_T_3_avg.dat" using 1:2:3 with yerrorb ls 6 notitle, \
#     "/tmp/sj_a_20_eps_001_T_3_avg.dat" using 1:2:3 with yerrorb ls 7 notitle, \
#     "/tmp/sj_a_10_eps_001_T_3_avg.dat" title "V = 10 [fm^3]"  with lines ls 2, \
#     "/tmp/sj_a_20_eps_001_T_3_avg.dat" title "V = 20 [fm^3]"  with lines ls 3, \
