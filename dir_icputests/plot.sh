#!/bin/bash

ref=-0.01549966

file=$1

for i in 10 100 1000 10000 100000 1000000; 
do
  grep " $i " $file > buf_paths$i.dat
done

rm test*.png

gnuplot<<EOF
se t pngcairo enhanced
set key outside center top horizontal samplen 2
set autoscale fix
set xzeroaxis
se output 'test01.png'
se xlabel 'paths'
se ylabel 'time [s]'
se log x 2
se log y
p [][1:] 'buf_paths10.dat' u 1:5 t 'paths 10' w lp,\
'buf_paths100.dat' u 1:5 t 'paths 100' w lp,\
'buf_paths1000.dat' u 1:5 t 'paths 1000' w lp,\
'buf_paths10000.dat' u 1:5 t 'paths 10000' w lp,\
'buf_paths100000.dat' u 1:5 t 'paths 100000' w lp,\
'buf_paths1000000.dat' u 1:5 t 'paths 10000000' w lp
EOF

gnuplot<<EOF
se t pngcairo enhanced
set key outside center top horizontal samplen 2
set autoscale fix
set xzeroaxis
se output 'test02.png'
se xlabel 'paths'
se ylabel '|(v_r - <<v>>)/v_r|'
se log x 2
se log y
p\
'buf_paths100.dat' u 1:(abs(\$3/$ref)) t 'paths 100' w lp,\
'buf_paths1000.dat' u 1:(abs(\$3/$ref)) t 'paths 1000' w lp,\
'buf_paths10000.dat' u 1:(abs(\$3/$ref)) t 'paths 10000' w lp,\
'buf_paths100000.dat' u 1:(abs(\$3/$ref)) t 'paths 100000' w lp,\
'buf_paths1000000.dat' u 1:(abs(\$3/$ref)) t 'paths 1000000' w lp,\
1 not w l lt 0
EOF

convert -append test0?.png test.png

rm buf*.dat

#'buf_paths10.dat' u 1:5 t 'paths 10' w lp,\
