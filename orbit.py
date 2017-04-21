
import sys
sys.path.insert(0, '/home/duxiaodi/GCOM_v3/GCOM/SOURCE')
sys.path.insert(0, '/home/duxiaodi/GCOM_v3/GCOM/RESONANCE')

import numpy as np
import creatobj
import geqdsk
import bgrid
import ugrid
import trace

#                   (R,    Z, phi, pitch, E)
inpu = creatobj.inpu(1.90, 0, 0.0,  0.8,  60)
inpu.nseg = 1000
inpu.nstep = 6000000
inpu.tstep = 1e-10
inpu.charge = 1.6*(1e-19)
inpu.mass = 2*1.67*(1e-27)
# note that for trapped particle,
# tstep will be automatically 0.1 of this value
# that is, tstep = 0.1*tstep
inpu.switch_full_orbit=1
inpu.switch_calc_prob = 1
#inpu.fn = '/home/duxiaodi/thomek/g165042.03705'
inpu.fn='/home/duxiaodi/GCOM/GCOM_v2/IN/g165037.03705'

# calculate equilibrum
g = geqdsk.read(inpu.fn)
g = bgrid.main(g)
g = ugrid.equ(g)

outpu = creatobj.outpu('NEED MORE BENCHMARK')
import timeit
t0 = timeit.default_timer()
outpu = trace.main(g,inpu,outpu)
t1 = timeit.default_timer()
print('consumed time [s] :', t1-t0 )



#obr = out['obr']
#obz = out['obz']
#
#if cal_prob == 1 and out['ob']>3:
#   obxyz = out['obxyz']
#   steps = out['steps']
#   pout = trace.probability(R0,Z0,steps,obxyz)
#   p = pout['p']
#   s = pout['s']
#
#   sol = out['sol']
#   tsol = out['tsol']













