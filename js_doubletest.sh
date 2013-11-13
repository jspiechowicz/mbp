#!/bin/bash

a=10.1882478074  #= 1.62151*2*pi
w=3.76991118431  # = 0.6*2*pi
gamma=0.81681408993 # = 0.13*2*pi
Dg=0.0001
Dp=0.0000064
lambda=1000

prog=idoublecputests
test -e $prog || make $prog

TDIR=dir_$prog
test -d $TDIR || mkdir -p $TDIR

force=0
paths=8192
periods=1000000

GSL_RNG_SEED=$RANDOM ./$prog --amp=$a --omega=$w --force=$force --gam=$gamma --Dg=$Dg --Dp=$Dp --lambda=$lambda --paths=$paths --periods=$periods --algorithm=euler --spp=500 > ${TDIR}/gauss_paths${paths}.dat


#Usage: ./idoublecputests <params> 
#
#Model params:
#    -a, --amp=FLOAT         set the AC driving amplitude 'amp' to FLOAT
#    -b, --omega=FLOAT       set the AC driving frequency '\omega' to FLOAT
#    -c, --force=FLOAT       set the external bias 'force' to FLOAT
#    -d, --gam=FLOAT         set the viscosity '\gamma' to FLOAT
#    -e, --Dg=FLOAT          set the Gaussian noise intensity 'Dg' to FLOAT
#    -f, --Dp=FLOAT          set the Poissonian noise intensity 'Dp' to FLOAT
#    -g, --lambda=FLOAT      set the Poissonian kicks frequency '\lambda' to FLOAT
#
#    -h, --comp=INT          choose between biased and unbiased Poissonian noise. INT can be one of:
#                            0: biased (default); 1: unbiased
#Simulation params:
#    -k, --paths=LONG        set the number of paths to LONG
#    -l, --periods=LONG      set the number of periods to LONG
#    -m, --trans=LONG        set the number of periods which stands for transients to LONG
#    -n, --spp=INT           specify how many integration steps should be calculated
#                            for a single period of the driving force
#
#    -o, --algorithm=STRING  sets the algorithm. STRING can be one of:
#                            predcorr: simplified weak order 2.0 adapted predictor-corrector
#                            euler: simplified weak order 1.0 regular euler-maruyama
#    -p, --help		prints this help and exits
#
