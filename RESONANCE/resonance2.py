
import sys
sys.path.insert(0, './SOURCE')
sys.path.insert(0, './WEIGHT')
sys.path.insert(0, './INOUT')
sys.path.insert(0, './C_RESULT')
sys.path.insert(0, '../')

import numpy as np
import trace
import gl
import timeit
from include import *
from joblib import Parallel, delayed
from progressbar import ProgressBar

R0 = np.arange(1.0,2.3,0.05)
Z0 = np.arange(-0.8,0.8,0.1)
E0 = np.arange(30,90,1)

nseg   = 1000
nstep  = 200000

# note that for trapped particle,
# tstep will be automatically 0.1 of this value
# that is, tstep = 0.1*tstep
tstep = 3e-10

switch_full_orbit=1
cal_prob = 1

output = np.zeros((len(Z0),len(E0),2,11))
def PARA_RESONANCE(i):
    for j in range(0,len(Z0)):
        for k in range(0,len(E0)):
            for l in range(0,2):

                pitch0 = (-1)**l*np.sqrt(1-30./E0[k])

                out = trace.main(
                        R0[i],Z0[j],E0[k],pitch0,
                        nseg,nstep,tstep,
                        switch_full_orbit,
                        cal_prob,
                        R,Z,RR,ZZ,
                        br,bz,bt,b,
                        bc,rmaxis,sibry,rbbbs,zbbbs,
                        nr,nz,psi,
                        f_br,f_bz,f_bt,f_b,
                        equ_curve,equ_grad
                        )

                obr = out['obr']
                obz = out['obz']
                rho = trace.RZ_to_RHO(obr,obz,f_psi,simag,sibry)

                output[j,k,l,0]  = R0[i]
                output[j,k,l,1]  = Z0[j]
                output[j,k,l,2]  = pitch0
                output[j,k,l,3]  = E0[k]
                output[j,k,l,4]  = out['ob']
                output[j,k,l,5]  = out['pphi']
                output[j,k,l,6]  = out['mu_E']
                output[j,k,l,7]  = out['f_phi']
                output[j,k,l,8]  = out['f_theta']
                output[j,k,l,9]  = np.min(rho)
                output[j,k,l,10] = np.max(rho)

    return output

if __name__ == '__main__':
    output = Parallel(n_jobs=16)(delayed(PARA_RESONANCE)(i) \
             for i in range(0,len(R0)) \
             )

np.savez('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/PARA_E_PPHI.npz',
      output = output)












