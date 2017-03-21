
import sys
sys.path.insert(0, './SOURCE')
sys.path.insert(0, './WEIGHT')
sys.path.insert(0, './INOUT')
sys.path.insert(0, './C_RESULT')
sys.path.insert(0, '../')

import numpy as np
import trace
import gl
import matplotlib.pylab as plt
import timeit
from include import *
import netCDF4
from joblib import Parallel, delayed
from progressbar import ProgressBar


plt.clf()
"""local lib"""

fn_nubeam = './IN/159243H06_fi_9.cdf'
fn = netCDF4.Dataset(fn_nubeam)
r2d = fn.variables['R2D'][:]/100. # [m]
z2d = fn.variables['Z2D'][:]/100. # [m]
energy = fn.variables['E_D_NBI'][:]/1000. # [keV]
pitch = fn.variables['A_D_NBI'][:]
fi = fn.variables['F_D_NBI'][:]

nseg   = 1000
nstep  = 400000
tstep = 3e-10
switch_full_orbit = 0
cal_prob = 0


output = np.zeros((len(pitch),len(energy),10))
def PARA_FI(i):
    for j in range(0,len(pitch)):
        for k in range(0,len(energy)):
            R0 = r2d[i]
            Z0 = z2d[i]
            pitch0 = pitch[j]
            E0 = energy[k]

            out = trace.main(
                  R0,Z0,E0,pitch0,
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

#           obtype = out['ob']
#
#            if obtype>-1:

            output[j,k,0]  = R0
            output[j,k,1]  = Z0
            output[j,k,2]  = pitch0
            output[j,k,3]  = E0
            output[j,k,4]  = out['ob']
            output[j,k,5]  = out['pphi']
            output[j,k,6]  = out['mu_E']
            output[j,k,7]  = fi[i,j,k]
            output[j,k,8]  = out['Rmax']
            output[j,k,9]  = out['Pmax']

    return output


if __name__ == '__main__':
   #st = timeit.default_timer()
   output = Parallel(n_jobs=16)(delayed(PARA_FI)(i) \
        for i in range(0,len(r2d)) \
        )
   #et = timeit.default_timer()
   #print(et-st)

np.savez('./RESONANCE/OUT/FI_NUBEAM.npz',
        output = output,)






