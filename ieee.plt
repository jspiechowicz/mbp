set terminal pngcairo enhanced
set key outside center top horizontal samplen 2
set autoscale fix
set xzeroaxis
set ylabel '{/Symbol \341}v{/Symbol \361}' rotate by 360

set output 'fig3a.png'
plot [0:0.2] 'fig3a_lmd2.dat' using (sqrt($1*2)):2 title '{/Symbol l} = 2' with lines linetype 1 linewidth 3,\
'fig3a_lmd20.dat' using (sqrt($1*20)):2 title '{/Symbol l} = 20' with lines linetype 2 linewidth 3,\
'fig3a_lmd200.dat' using (sqrt($1*200)):2 title '{/Symbol l} = 200' with lines linetype 3 linewidth 3,\
'fig3a_lmd2000.dat' using (sqrt($1*2000)):2 title '{/Symbol l} = 2000' with lines linetype 4 linewidth 3
unset output

set logscale x
set format x '10^{%L}'
set output 'fig3b.png'
set xlabel 'D_G'
plot 'fig3b_lmd2.dat' using 1:2 title '{/Symbol l} = 2' with lines linetype 1 linewidth 3,\
'fig3b_lmd20.dat' using 1:2 title '{/Symbol l} = 20' with lines linetype 2 linewidth 3,\
'fig3b_lmd200.dat' using 1:2 title '{/Symbol l} = 200' with lines linetype 3 linewidth 3,\
'fig3b_lmd2000.dat' using 1:2 title '{/Symbol l} = 2000' with lines linetype 4 linewidth 3
unset output

set output 'fig5.13a.png'
unset logscale x
set format x '%.2f'
set xlabel '{/Symbol \341}{/Symbol h}{/Symbol \361}'
plot [0:0.95] 'fig5.13a_lmd20.dat' using (sqrt($1*20)):2 title '{/Symbol l} = 20' with lines linetype 1 linewidth 3,\
'fig5.13a_lmd200.dat' using (sqrt($1*200)):2 title '{/Symbol l} = 200' with lines linetype 2 linewidth 3,\
'fig5.13a_lmd2000.dat' using (sqrt($1*2000)):2 title '{/Symbol l} = 2000' with lines linetype 3 linewidth 3
unset output

set output 'fig5.13b.png'
set xlabel '{/Symbol \341}{/Symbol h}{/Symbol \361}'
plot [0:0.95] 'fig5.13a_lmd2000_Dg1e-06.dat' using (sqrt($1*2000)):2 title 'D_G = 10^{-6}' with lines linetype 1 linewidth 3,\
'fig5.13a_lmd2000.dat' using (sqrt($1*2000)):2 title 'D_G = 3*10^{-4}' with lines linetype 2 linewidth 3,\
'fig5.13a_lmd2000_Dg0.01.dat' using (sqrt($1*2000)):2 title 'D_G = 10^{-2}' with lines linetype 3 linewidth 3
unset output
exit gnuplot
