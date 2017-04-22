import sys
import numpy as np
from joblib import Parallel, delayed
sys.path.insert(0, '/home/duxiaodi/GCOM/SOURCE')
sys.path.insert(0, '/home/duxiaodi/GCOM/RESONANCE')
import creatobj
import geqdsk
import bgrid
import ugrid
import trace

"""
-------------------
DEVELOPER MESSAGE |
-------------------
+ 20170421 SOVLE THE DATA SCATTER FOR FISHBONE RESONANCE
- Scan all R,Z,PITCH is dumb one. But It should not be different from R, pitch only
But it turns out that the it indeed different.

- Scan R, Z, Pitch give scatters data. This might be due to the starting point is close to banana tip, where particle stays there for lots of time. If the steps is not short enough, this will gives significant error bar!i

- Instead, if starting point from midplane, when the time fraction spent there in mid-plane is fairly short, therefore the uncertanties due to the steps is fairly small, indeed negligible.
Therefore make sure that we should always scan from the mid=plane

CALCULATION HISTORY
20170412 (for Kathereen)
#165037 fishbone, # 165865 AE

"""

# dummy, usually do not need to change
inpu = creatobj.inpu(0, 0, 0, 0, 0)
inpu.shot = 165037
inpu.cores = 12
inpu.nseg = 1000
inpu.nstep = 200000
inpu.tstep = 3e-10
inpu.charge = 1.6*(1e-19)
inpu.mass = 2*1.67*(1e-27)
inpu.switch_full_orbit = 1
inpu.switch_calc_freq = 1
inpu.gfile='/home/duxiaodi/gfile/g'+np.str(inpu.shot)+'.03705'

# generate scan step
fn_nubeam = '/home/duxiaodi/thomek/dist_60keV_1M_'+np.str(inpu.shot)+'.dat'
tmp = np.loadtxt(fn_nubeam)
inpu.mc = 200000
inpu.R_array = tmp[0:inpu.mc,0]/1e2
inpu.Z_array = tmp[0:inpu.mc,1]/1e2
inpu.pitch_array = tmp[0:inpu.mc,2]
inpu.E_array = tmp[0:inpu.mc,3]/1e3
inpu.phi0 = 0

# prepare output
outpu = creatobj.outpu('Kathreen_fi_'+np.str(inpu.shot)+'_'+np.str(inpu.E0)+'keV')
outpu.path = '/home/duxiaodi/GCOM_v3/GCOM/RESONANCE/output/'

# calculate equilibrum
g = geqdsk.read(inpu.gfile)
g = bgrid.main(g)
g = ugrid.equ(g)

# real calculation begins from here
output = np.zeros((inpu.mc,11))
def PARA_SCAN(g,inpu,outpu,i):

    inpu.E0 = inpu.E_array[i]
    inpu.R0 = inpu.R_array[i]
    inpu.Z0 = inpu.Z_array[i]
    inpu.pitch0 = inpu.pitch_array[i]

    outpu = trace.main(g,inpu,outpu)

    output[0]  = inpu.R0
    output[1]  = inpu.Z0
    output[2]  = inpu.pitch0
    output[3]  = inpu.E0
    output[4]  = outpu.ob
    output[5]  = outpu.pphi
    output[6]  = outpu.mu_E
    output[7]  = outpu.f_phi
    output[8]  = outpu.f_theta
    output[9]  = np.min(outpu.rho)
    output[10] = np.max(outpu.rho)

    return output


if __name__ == '__main__':
    output = Parallel(n_jobs=inpu.cores,verbose=5)(delayed(PARA_SCAN)(g,inpu,outpu,i) \
             for i in range(0,inpu.mc))

    outpu.fn = outpu.path+outpu.comment+'.npz'
    np.savez(outpu.fn,output=output)












