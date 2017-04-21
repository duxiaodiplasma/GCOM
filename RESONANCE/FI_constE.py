
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
import matplotlib.pylab as plt
import timeit
from include import *
import netCDF4
import mcparticles_from_df
from joblib import Parallel, delayed


plt.close('all')
"""local lib"""
fn_nubeam = '/home/duxiaodi/GCOM/GCOM_v2/IN/history/159243H06_fi_9.cdf'
#fn_nubeam = '/home/duxiaodi/thomek/165037Z46_birth.cdf1'
mcnum = 5000
Eini = 80
R_array,Z_array,E_array,P_array =  \
mcparticles_from_df.main(fn_nubeam,mcnum,Eini)

nseg   = 1000
nstep  = 400000
tstep = 3e-10
switch_full_orbit = 0
cal_prob = 0


output = np.zeros((mcnum,10))
#def PARA_FI(R_array,Z_array,E_array,P_array,i):
for i in range(0,mcnum):
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
          bc,rmaxis,sibry,simag,rbbbs,zbbbs,
          nr,nz,psi,
          f_br,f_bz,f_bt,f_b,f_psi,
          equ_curve,equ_grad
          )

#   obtype = out['ob']
#    if obtype>-1:

    output[i,0]  = R0
    output[i,1]  = Z0
    output[i,2]  = pitch0
    output[i,3]  = E0
    output[i,4]  = out['ob']
    output[i,5]  = out['pphi']
    output[i,6]  = out['mu_E']
    output[i,7]  = 0
    output[i,8]  = out['Rmax']
    output[i,9]  = out['Pmax']

    #return output


#if __name__ == '__main__':
#   #st = timeit.default_timer()
#   output_all = Parallel(n_jobs=1)(delayed(PARA_FI)(R_array,Z_array,E_array,P_array,i) \
#        for i in range(0,mcnum) \
#        )
#   #et = timeit.default_timer()
#   #print(et-st)

np.savez('./RESONANCE/OUT/FI_'+np.str(Eini)+'.npz',
        output = output)






