#!/bin/bash

prog=icputests
test -e $prog || make $prog

TDIR=dir_$prog
test -d $TDIR || mkdir -p $TDIR

force=0.1
eta=$force
paths=8192
periods=1000000
GSL_RNG_SEED=$RANDOM ./$prog --force=$force --paths=$paths --periods=$periods > ${TDIR}/gauss_paths${paths}.dat
for lam in 4 16 64 512
do
  Dp=$(echo ${eta}*${eta}/${lam}|bc -l)
  GSL_RNG_SEED=$RANDOM ./$prog --Dp=$Dp --paths=$paths --periods=$periods --lambda=$lam > ${TDIR}/poisson_lambda${lam}_paths${paths}.dat
done
