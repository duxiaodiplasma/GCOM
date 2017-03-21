import sys
sys.path.insert(0, './SOURCE')
sys.path.insert(0, './WEIGHT')
sys.path.insert(0, './INOUT')
sys.path.insert(0, './C_RESULT')


import wght_ini
import numpy as np
import trace
import gl
from include import *
from joblib import Parallel, delayed

def p_ob(ch):

    fn_fida = './INOUT/159243.inpa801_npa_weights.h5'
    fn_grid = './INOUT/gridf.h5'
    ini = wght_ini.main(fn_fida,fn_grid)

    gx = ini['x']
    gy = ini['y']
    gz = ini['z']
    ngx = ini['nx']
    ngy = ini['ny']
    ngz = ini['nz']

    E = ini['E']
    phit = ini['phit']
    pitch = ini['pitch']
    nch = phit.shape[0]

    nseg   = 1000
    nstep  = 1000000
    tstep = 1e-8
    switch_full_orbit = 0
    cal_prob = 0

    # RESULTS in these guys

#    for ch in range(0,nch):
    for i in range(0,ngz):
        for j in range(0,ngy):
            for k in range(0,ngx):
#
#                    if phit[ch,i,j,k]>0:
#                       # loop Energy
                for l in range(0,len(E)):

                    R0 = np.sqrt(gx[i,j,k]**2+gy[i,j,k]**2)
                    Z0 = gz[i,j,k]
                    E0 = E[l]
                    print(R0,Z0,E0)
                    print(ch,i,j,k)

    out = trace.main(
          R0,Z0,E0,0.5,
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


    return {'done'}

if __name__ == '__main__':
   result = Parallel(n_jobs=3)(delayed(p_ob)(ch) \
        for ch in range(0,1))
#        for i in range(0,50)
#        for j in range(0,50)
#        for k in range(0,50)
#        for l in range(0,1))



