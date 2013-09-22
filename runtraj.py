#!/usr/bin/python
import commands, os
import numpy

#Model
amp = 4.2
omega = 4.9
force = 0.1
gam = 0.9
Dg = 0.001
Dp = 0
lmd = 0
comp = 0

#Simulation
dev = 3
block = 64
paths = 2**12 #4096
periods = 500
spp = 100
algorithm = 'predcorr'
trans = 0.1

#Output
mode = 'trajectory'
DIRNAME='./tests/traj/'
os.system('mkdir -p %s' % DIRNAME)

os.system('rm -v %s*.dat %s*.png' % (DIRNAME, DIRNAME))

out = 'traj'
_cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%d --spp=%d --algorithm=%s --mode=%s >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, mode, out)
output = open('%s.dat' % out, 'w')
print >>output, '#%s' % _cmd
output.close()
print _cmd
cmd = commands.getoutput(_cmd)
os.system('mv -v %s.dat %s.png %s' % (out, out, DIRNAME))
