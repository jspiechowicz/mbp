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
paths = 8192
periods = 4000
spp = 800
algorithm = 'euler'
trans = 0.1

#Output
mode = 'moments'
points = 100
beginxx = 0
endxx = 0.2
domain = '1d'
domainx = 'f'
logx = 0
DIRNAME='./tests/ieee/'
os.system('mkdir -p %s' % DIRNAME)

os.system('rm -v %s*.dat %s*.png' % (DIRNAME, DIRNAME))

#fig 3a

beginx = 0
c = 0.2
domainx = 'p'
for lmd in [0.2, 2, 20, 200, 2000]:
    endx = c**2/lmd
    out = 'fig3a_lmd%s' % lmd
    _cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --domain=%s --domainx=%s --logx=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, domain, domainx, logx, out)
    output = open('%s.dat' % out, 'w')
    print >>output, '#%s' % _cmd
    output.close()
    print _cmd
    cmd = commands.getoutput(_cmd)

#fig 3b

c = 0.087
domainx = 'D'
beginx = -6
endx = 0
logx = 1
for lmd in [2, 20, 200, 2000]:
    if lmd == 2:
        Dp = 8.0e-4
    elif lmd == 20:
        Dp = 3.8e-4
    elif lmd == 200:
        Dp = 3.6e-5
    elif lmd == 2000:
        Dp = 3.8e-6
    out = 'fig3b_lmd%s' % lmd
    _cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --domain=%s --domainx=%s --logx=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, domain, domainx, logx, out)
    output = open('%s.dat' % out, 'w')
    print >>output, '#%s' % _cmd
    output.close()
    print _cmd
    cmd = commands.getoutput(_cmd)

#fig 5.13a

amp = 4.2
omega = 4.9
force = 0
gam = 0.9
Dg = 0.0003
Dp = 0
lmd = 0
beginx = 0
c = 0.95
domainx = 'p'
logx = 0
spp = 2000
for lmd in [20, 200, 2000]:
    endx = c**2/lmd
    out = 'fig5.13a_lmd%s' % lmd
    _cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --domain=%s --domainx=%s --logx=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, domain, domainx, logx, out)
    output = open('%s.dat' % out, 'w')
    print >>output, '#%s' % _cmd
    output.close()
    print _cmd
    cmd = commands.getoutput(_cmd)

#fig 5.13b

amp = 4.2
omega = 4.9
force = 0
gam = 0.9
Dg = 0.0003
Dp = 0
lmd = 2000
beginx = 0
c = 0.95
domainx = 'p'
logx = 0
spp = 2000
for Dg in [1.0e-6, 1.0e-4, 1.0e-2]:
    endx = c**2/lmd
    out = 'fig5.13a_lmd%s_Dg%s' % (lmd, Dg)
    _cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --domain=%s --domainx=%s --logx=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, domain, domainx, logx, out)
    output = open('%s.dat' % out, 'w')
    print >>output, '#%s' % _cmd
    output.close()
    print _cmd
    cmd = commands.getoutput(_cmd)

os.system('gnuplot ieee.plt')
os.system('mv -vf fig*.dat fig*.png %s' % DIRNAME)
