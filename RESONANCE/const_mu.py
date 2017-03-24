
import sys
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/SOURCE')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/WEIGHT')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/INOUT')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/C_RESULT')

import numpy as np
import trace
import gl
import timeit
from include import *
from joblib import Parallel, delayed
from progressbar import ProgressBar


cores = 16
R0 = np.arange(1.0,2.3,0.05)
Z0 = np.arange(-0.8,0.8,0.1)
E0 = np.arange(20,81,1)
mu_ini = 4

# hard coated in loop ...
#pitch0 = np.sqrt(1-mu_ini*f_b(Z0,R0)/E0)

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

                if l == 0:
                   pitch0 = np.sqrt(1-mu_ini*f_b(Z0[j],R0[i])/E0[k])
                if l == 1:
                   pitch0 = -1.0*np.sqrt(1-mu_ini*f_b(Z0[j],R0[i])/E0[k])


                if pitch0 >= -1 and pitch0 <= 1:
                   # debug   print(R0[i],Z0[j],E0[k],f_b(Z0[j],R0[i]),pitch0)
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
    output = Parallel(n_jobs=cores)(delayed(PARA_RESONANCE)(i) \
             for i in range(0,len(R0)) \
             )

np.savez('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/const_mu_'+np.str(mu_ini)+'.npz',
      output = output)












