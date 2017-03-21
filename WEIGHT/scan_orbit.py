import sys
sys.path.insert(0, './SOURCE')
sys.path.insert(0, './WEIGHT')
sys.path.insert(0, './INOUT')
sys.path.insert(0, './C_RESULT')

import numpy as np
from progressbar import ProgressBar

"""local lib"""
import trace
import gl
from include import *
import wght_ini

"""THE CODE CALCULATE THE PROBABILITY OF ORBIT AT THE GRID POINT"""

fn_fida = './IN/159243.inpa801_npa_weights.h5'
fn_grid = './IN/beamgrid.h5'
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
tstep = 5e-9
switch_full_orbit = 1
cal_prob = 1


# RESULTS in these guys
#all_obs = np.zeros((len(E),ngx,ngy,ngz,2,nseg))
#all_pitch = np.zeros((len(E),ngx,ngy,ngz,nseg))
#p_all_obs = np.zeros((nch,len(E),ngx,ngy,ngz))


#with open("wght_"+str(nch)+".dat",'w') as myfile:
with open("./OUT/ch1_ssnpa.wght",'w') as myfile:
   pbar = ProgressBar()
   for ch in range(0,1):
       print('Channel: ** ('+str(ch)+') **')
       for i in pbar(range(0,ngz)):
           for j in range(0,ngy):
               for k in range(0,ngx):
                   for l in range(10,len(E)-10):

                       if phit[ch,i,j,k]>0:

                          R0 = np.sqrt(gx[i,j,k]**2+gy[i,j,k]**2)
                          Z0 = gz[i,j,k]
                          E0 = E[l]
                          pitch0 = pitch[ch,i,j,k]

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
                          # if it is a confined orbit
                          if obtype >3:
                             obr = out['obr']
                             obz = out['obz']
                             obpitch = out['pitch']
                             obxyz = out['obxyz']
                             steps = out['steps']
                             pout = trace.probability(R0,Z0,steps,obxyz)
                             p = pout['p']
                             s = pout['s']

                             lines = []
                             lines.append(
                                 '{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(
                                 ch,l,k,j,i,R0,Z0,E0,pitch0,p[0],obtype[0],
                                 np.max(obr),obpitch[np.argmax(obr)],obz[np.argmax(obr)])
                                 )
                             myfile.writelines(lines)

                             print(ch,l,k,j,i,R0,Z0,E0,pitch0,p[0],obtype[0],
                                   np.max(obr),obpitch[np.argmax(obr)],obz[np.argmax(obr)])






