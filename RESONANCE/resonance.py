
import sys
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/SOURCE')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/RESONANCE')

import numpy as np
import trace
import gl
import timeit
from include import *
from joblib import Parallel, delayed
from progressbar import ProgressBar


R0 = np.arange(1.05,2.35,0.05)
Z0 = np.arange(-1.2,1.2,0.1)
pitch0 = np.arange(-1,1.05,0.05)

E0 = 60
#comment = '165042'
#comment = '165037'
comment = '165037_benchmark_final_old'
#comment = '159243_p'

nseg   = 1000
nstep  = 200000

# note that for trapped particle,
# tstep will be automatically 0.1 of this value
# that is, tstep = 0.1*tstep
tstep = 3e-10

switch_full_orbit=1
cal_prob = 1

output = np.zeros((len(Z0),len(pitch0),11))
def PARA_RESONANCE(i):
    for j in range(0,len(Z0)):
        for k in range(0,len(pitch0)):

            out = trace.main(
                    R0[i],Z0[j],E0,pitch0[k],
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

            obr = out['obr']
            obz = out['obz']
            rho = out['rho']

            print(out['ob'])
            output[j,k,0]  = R0[i]
            output[j,k,1]  = Z0[j]
            output[j,k,2]  = pitch0[k]
            output[j,k,3]  = E0
            output[j,k,4]  = out['ob']
            output[j,k,5]  = out['pphi']
            output[j,k,6]  = out['mu_E']
            output[j,k,7]  = out['f_phi']
            output[j,k,8]  = out['f_theta']
            output[j,k,9]  = np.min(rho)
            output[j,k,10] = np.max(rho)

    return output

if __name__ == '__main__':
    output = Parallel(n_jobs=10,verbose=5)(delayed(PARA_RESONANCE)(i) \
             for i in range(0,len(R0)) \
             )
np.savez('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/Kathreen_'+np.str(E0)+'_'+np.str(comment)+'.npz',
      output = output)












