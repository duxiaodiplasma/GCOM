import sys
sys.path.insert(0, './SOURCE')
sys.path.insert(0, './WEIGHT')
sys.path.insert(0, './INOUT')
sys.path.insert(0, './C_RESULT')

import numpy as np
import netCDF4
from progressbar import ProgressBar

"""local lib"""
import trace
import gl
from include import *

fn_nubeam = './IN/159243H06_fi_9.cdf'
fn = netCDF4.Dataset(fn_nubeam)
r2d = fn.variables['R2D'][:]/100. # [m]
z2d = fn.variables['Z2D'][:]/100. # [m]
energy = fn.variables['E_D_NBI'][:]/1000. # [keV]
pitch = fn.variables['A_D_NBI'][:]
fi = fn.variables['F_D_NBI'][:]

nseg   = 1000
nstep  = 400000
tstep = 3e-9
switch_full_orbit = 0
cal_prob = 0

output = np.zeros((5,len(r2d),len(pitch),len(energy)))

pbar = ProgressBar()
for i in pbar(range(0,len(r2d))):
    for j in range(0,len(pitch)):
        for k in range(0,len(energy)): # energy 12keV to 50keV

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

            obtype = out['ob']

            if obtype>3:
               output[0,i,j,k] = obtype
               output[1,i,j,k] = out['Rmax']
               output[2,i,j,k] = out['Pmax']
               output[3,i,j,k] = E0
               output[4,i,j,k] = fi[i,j,k]
               #print(obtype,out['Rmax'],out['Pmax'],E0,fi[i,j,k])

np.savez('./OUT/FI_RMAX2.npz',
          output = output,)






