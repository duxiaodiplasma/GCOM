
import sys
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/SOURCE')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/RESONANCE')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/WEIGHT')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/INOUT')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/C_RESULT')
sys.path.insert(0, '/home/duxiaodi/GCOM/GCOM_v2/')

import numpy as np
import trace
import gl
import timeit
from include import *
import netCDF4
import mcparticles_from_df
from joblib import Parallel, delayed


"""local lib"""
#fn_nubeam = '/home/duxiaodi/GCOM/GCOM_v2/IN/history/159243H06_fi_9.cdf'
#fn_nubeam = '/home/duxiaodi/thomek/165037Z46_fi_1.cdf'
#mcnum = 10000
#R_array,Z_array,E_array,P_array =  \
#mcparticles_from_df.main(fn_nubeam,mcnum,Eini)

Eini = 80
mcnum=20000
fn_nubeam = '/home/duxiaodi/thomek/165037/dist_80kev_1M_165037.dat'
tmp = np.loadtxt(fn_nubeam)
R_array = tmp[0:mcnum,0]/1e2
Z_array = tmp[0:mcnum,1]/1e2
P_array = tmp[0:mcnum,2]
E_array = tmp[0:mcnum,3]/1e3

nseg   = 1000
nstep  = 400000
tstep = 3e-10
switch_full_orbit = 0
cal_prob = 0


def PARA_FI(R_array,Z_array,E_array,P_array,i):
    output = np.zeros(10)
    R0 = R_array[i]
    Z0 = Z_array[i]
    E0 = E_array[i]
    pitch0 = P_array[i]

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

    output[0]  = R0
    output[1]  = Z0
    output[2]  = pitch0
    output[3]  = E0
    output[4]  = out['ob']
    output[5]  = out['pphi']
    output[6]  = out['mu_E']
    output[7]  = 0
    output[8]  = out['Rmax']
    output[9]  = out['Pmax']

    return output


if __name__ == '__main__':
   output = Parallel(n_jobs=12,verbose=5)(delayed(PARA_FI)(R_array,Z_array,E_array,P_array,i) \
        for i in range(0,mcnum) \
        )

   np.savez('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/FI_'+np.str(Eini)+'_165037a.npz',
           output = np.asarray(output))






