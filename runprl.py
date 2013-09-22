#!/usr/bin/python
import commands, os
import numpy

#Model
amp = 4.2
omega = 4.9
force = 0.0
gam = 0.9
Dg = 0.001
Dp = 0
lmd = 0
comp = 0

#Simulation
dev = 3
block = 64
paths = 4096
periods = 2000
spp = 100
algorithm = 'predcorr'
trans = 0.1

#Output
mode = 'moments'
points = 100
beginx = 0
endx = 0.2
domain = '1d'
domainx = 'f'
logx = 0
DIRNAME='./tests/prl/'
os.system('mkdir -p %s' % DIRNAME)

os.system('rm -v %s*.dat %s*.png' % (DIRNAME, DIRNAME))

out = 'prl'
_cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --domain=%s --domainx=%s --logx=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, domain, domainx, logx, out)
output = open('%s.dat' % out, 'w')
print >>output, '#%s' % _cmd
output.close()
print _cmd
cmd = commands.getoutput(_cmd)
os.system('gnuplot prl.plt')
os.system('mv -v %s.dat %s.png %s' % (out, out, DIRNAME))
