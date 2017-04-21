import sys
import numpy as np
from joblib import Parallel, delayed
sys.path.insert(0, '/home/duxiaodi/GCOM_v3/GCOM/SOURCE')
sys.path.insert(0, '/home/duxiaodi/GCOM_v3/GCOM/RESONANCE')
import creatobj
import geqdsk
import bgrid
import ugrid
import trace


# dummy, usually do not need to change
inpu = creatobj.inpu(0, 0, 0, 0, 0)
inpu.gfile='/home/duxiaodi/GCOM/GCOM_v2/IN/g165037.03705'
inpu.cores = 12
inpu.nseg = 1000
inpu.nstep = 200000
inpu.tstep = 3e-10
inpu.charge = 1.6*(1e-19)
inpu.mass = 2*1.67*(1e-27)
inpu.switch_full_orbit = 1
inpu.switch_calc_freq = 1

# real input for scanning
inpu.R_array = np.arange(1.2,2.25,0.01)
inpu.Z_array = np.arange(0,0.1,0.1)
inpu.pitch_array = np.arange(-1,1.05,0.01)
inpu.phi0 = 0
inpu.E0 = 60

# prepare output
# generate scan step
outpu = creatobj.outpu('Kathreen_165037_60keV')
outpu.path = '/home/duxiaodi/GCOM_v3/GCOM/RESONANCE/output/'

# calculate equilibrum
g = geqdsk.read(inpu.gfile)
g = bgrid.main(g)
g = ugrid.equ(g)

# real calculation begins from here
output = np.zeros((len(inpu.Z_array),len(inpu.pitch_array),11))
def PARA_SCAN(g,inpu,outpu,i):

    for j in range(0,len(inpu.Z_array)):
        for k in range(0,len(inpu.pitch_array)):

            inpu.R0 = inpu.R_array[i]
            inpu.Z0 = inpu.Z_array[j]
            inpu.pitch0 = inpu.pitch_array[k]

            outpu = trace.main(g,inpu,outpu)

            output[j,k,0]  = inpu.R0
            output[j,k,1]  = inpu.Z0
            output[j,k,2]  = inpu.pitch0
            output[j,k,3]  = inpu.E0
            output[j,k,4]  = outpu.ob
            output[j,k,5]  = outpu.pphi
            output[j,k,6]  = outpu.mu_E
            output[j,k,7]  = outpu.f_phi
            output[j,k,8]  = outpu.f_theta
            output[j,k,9]  = np.min(outpu.rho)
            output[j,k,10] = np.max(outpu.rho)

    return output


if __name__ == '__main__':
    output = Parallel(n_jobs=inpu.cores,verbose=5)(delayed(PARA_SCAN)(g,inpu,outpu,i) \
             for i in range(0,len(inpu.R_array)))

    outpu.fn = outpu.path+outpu.comment+'.npz'
    np.savez(outpu.fn,output=output)












