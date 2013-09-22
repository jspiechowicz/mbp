#!/usr/bin/python
import commands, os
import numpy

pipi = 2.0*numpy.pi

#Model
amp = 0
omega = 0.6*pipi
force = 0.05*pipi
gam = 0
Dg = 0
Dp = 0
lmd = 0
comp = 0

#Simulation
dev = 3
block = 64
paths = 4096
periods = 2000
spp = 100
algorithm= 'predcorr'
trans = 0.1

#Output
mode = 'moments'
points = 100
beginx = 0.1*pipi
endx = 1.0*pipi
beginy = 0.0
endy = 6.0*pipi
domain = '2d'
domainx = 'g'
domainy = 'a'
logx = 0
logy = 0
DIRNAME='./tests/mj2d/'
os.system('mkdir -p %s' % DIRNAME)

os.system('rm -v %s*.dat %s*.png' % (DIRNAME, DIRNAME))

#fig 5.1a

out = 'fig5.1a'
_cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --beginy=%s --endy=%s --domain=%s --domainx=%s --domainy=%s --logx=%d --logy=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, beginy, endy, domain, domainx, domainy, logx, logy, out)
output = open('%s.dat' % out, 'w')
print >>output, '#%s' % _cmd
output.close()
print _cmd
cmd = commands.getoutput(_cmd)

#fig 5.1b

gam = 0.19*pipi
beginx = 0.1*pipi
endx = 2.0*pipi
domainx = 'w'

out = 'fig5.1b'
_cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --beginy=%s --endy=%s --domain=%s --domainx=%s --domainy=%s --logx=%d --logy=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, beginy, endy, domain, domainx, domainy, logx, logy, out)
output = open('%s.dat' % out, 'w')
print >>output, '#%s' % _cmd
output.close()
print _cmd
cmd = commands.getoutput(_cmd)

#fig 5.2b

omega = 0.6*pipi
force = 0.001*pipi
beginx = 1.4*pipi
endx = 2.0*pipi
beginy = 0.05*pipi
endy = 0.25*pipi
domainx = 'a'
domainy = 'g'

out = 'fig5.2b'
_cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --beginy=%s --endy=%s --domain=%s --domainx=%s --domainy=%s --logx=%d --logy=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, beginy, endy, domain, domainx, domainy, logx, logy, out)
output = open('%s.dat' % out, 'w')
print >>output, '#%s' % _cmd
output.close()
print _cmd
cmd = commands.getoutput(_cmd)

#fig 5.2c

force = 0.01*pipi

out = 'fig5.2c'
_cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --beginy=%s --endy=%s --domain=%s --domainx=%s --domainy=%s --logx=%d --logy=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, beginy, endy, domain, domainx, domainy, logx, logy, out)
output = open('%s.dat' % out, 'w')
print >>output, '#%s' % _cmd
output.close()
print _cmd
cmd = commands.getoutput(_cmd)

#fig 5.2d

force = 0.06*pipi

out = 'fig5.2d'
_cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --beginy=%s --endy=%s --domain=%s --domainx=%s --domainy=%s --logx=%d --logy=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, beginy, endy, domain, domainx, domainy, logx, logy, out)
output = open('%s.dat' % out, 'w')
print >>output, '#%s' % _cmd
output.close()
print _cmd
cmd = commands.getoutput(_cmd)

#fig 5.6b-d

force = 0.01*pipi

for Dg in [0.0001, 0.0005, 0.001]:
    out = 'fig5.6_Dg%s' % Dg
    _cmd = './prog --dev=%d --amp=%s --omega=%s --force=%s --gam=%s --Dg=%s --Dp=%s --lambda=%s --comp=%d --block=%d --paths=%d --periods=%s --spp=%d --algorithm=%s --trans=%s --mode=%s --points=%d --beginx=%s --endx=%s --beginy=%s --endy=%s --domain=%s --domainx=%s --domainy=%s --logx=%d --logy=%d >> %s.dat' % (dev, amp, omega, force, gam, Dg, Dp, lmd, comp, block, paths, periods, spp, algorithm, trans, mode, points, beginx, endx, beginy, endy, domain, domainx, domainy, logx, logy, out)
    output = open('%s.dat' % out, 'w')
    print >>output, '#%s' % _cmd
    output.close()
    print _cmd
    cmd = commands.getoutput(_cmd)

os.system('gnuplot mj2d.plt')
os.system('mv -v *.dat *.png %s' % DIRNAME)
