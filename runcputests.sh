#!/bin/bash

TDIR=CPUTESTS
test -d $TDIR || mkdir -p $TDIR

force=0.1
eta=$force
for paths in 4 8 16 32 64 128 256 512 1024 2048 4096 8162
do
  ./cputests --force=$force --paths=$paths --periods=1000000 > ${TDIR}/gauss_paths${paths}.dat
  for lam in 4 16 64 512
  do
    Dp=$(echo ${eta}*${eta}/${lam}|bc -l)
    ./cputests --Dp=$Dp --paths=$paths --periods=1000000 --lambda=$lam > ${TDIR}/poisson_lambda${lam}_paths${paths}.dat
  done
done
