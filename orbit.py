
import sys
sys.path.insert(0, './SOURCE')
sys.path.insert(0, './WEIGHT')
sys.path.insert(0, './INOUT')
sys.path.insert(0, './C_RESULT')

import numpy as np
import trace
#import trace2
import gl
import matplotlib.pylab as plt
from include import *

plt.clf()

R0 = 1.15
Z0 = -0.49
pitch0 = -0.9
E0 = 60

nseg   = 1000
nstep  = 2000000

# note that for trapped particle,
# tstep will be automatically 0.1 of this value
# that is, tstep = 0.1*tstep
tstep = 3e-10

switch_full_orbit=1
cal_prob = 1

import timeit
t0 = timeit.default_timer()
out = trace.main(
        R0,Z0,E0,pitch0,
        nseg,nstep,tstep,
        switch_full_orbit,
        cal_prob,
        R,Z,RR,ZZ,
        br,bz,bt,b,
        bc,rmaxis,sibry,simag,rbbbs,zbbbs,
        nr,nz,psi,
        f_br,f_bz,f_bt,f_b,f_psi,
        equ_curve,equ_grad
        )

t1 = timeit.default_timer()
print(t1-t0)

obr = out['obr']
obz = out['obz']
rho = out['rho']

if cal_prob == 1 and out['ob']>3:
   obxyz = out['obxyz']
   steps = out['steps']
   pout = trace.probability(R0,Z0,steps,obxyz)
   p = pout['p']
   s = pout['s']

   sol = out['sol']
   tsol = out['tsol']













