"""
;+
; NAME:
;     bfield
; PURPOSE:
      calcualte magnetic field strength from G-file
; EXPLANATION:
;
; CALLING SEQUENCE:
	geqdsk.py   (read gfile)
	util_tokamak.py	(tools for reading gfile)
#;
; OPTIONAL INPUT:

; OPTIONAL KEYWORD INPUT:
;
; OPTIONAL KEYWORD OUTPUTS:

; PROCEDURE:

; EXAMPLE:

; RESTRICTIONS:
;
; PROCEDURES USED:
;
; FIRST USABLE VERSION
	10-06-2016 by X.D.DU
; REVISION HISTORY:
"""

import gl
import geqdsk
import numpy as np
from scipy import interpolate
import sys

def main(lr,lz,phi,fn):

   # read gfile
   res = geqdsk.read(fn)

   nr = res['nw']
   nz = res['nh']
   r = res['r']
   z = res['z']
   psirz = res['psirz']
   rc = res['rcentr']
   bc = res['bcentr']
   rbbbs = res['rbbbs']
   zbbbs = res['zbbbs']
   sibry = res['sibry']
   simag = res['simag']

   # first deravitive of the psi in RZ direction
   dpsi_dr = np.zeros([nr,nz])
   dpsi_dz = np.zeros([nr,nz])

   for i in range(1, nr-1):
       for j in range(1, nz-1):
           dpsi_dr[i,j] = (psirz[i+1,j]-psirz[i-1,j])/(r[i+1,j]-r[i-1,j])
           dpsi_dz[i,j] = (psirz[i,j+1]-psirz[i,j-1])/(z[i,j+1]-z[i,j-1])

   # 2D interpolation for dpsirz
   rr = np.zeros(nr)
   zz = np.zeros(nz)
   for i in range(0, nr):
       rr[i] = r[i, 0]
   for i in range(0, nz):
       zz[i] = z[0, i]

   fr = interpolate.interp2d(zz, rr, dpsi_dr, kind='cubic')
   fz = interpolate.interp2d(zz, rr, dpsi_dz, kind='cubic')

   # local magnetic field strength
   l_dpsir = fr(lz,lr)
   l_dpsiz = fz(lz,lr)

   # local psi
   #fpsi = interpolate.interp2d(zz, rr, (psirz-simag)/(sibry-simag), kind='cubic')
   fpsi = interpolate.interp2d(zz, rr, psirz, kind='cubic')
   psi = fpsi(lz,lr)

   #br = -1.0*-(1/lr)*l_dpsiz
   #bz = -1.0* (1/lr)*l_dpsir # -1.0 for unify the defination
   # right hand?
   br = -1.0*-(1/lr)*l_dpsiz
   bz = -1.0* (1/lr)*l_dpsir # -1.0 for unify the defination
			     # of the plus current direction in DIII-D,$
			     # which is CCW direction for plus

   # -1 for right-handed coordinate
   bt = np.asarray(rc*bc/lr)
   # print(br,bz,bt)

   # rz coordinate to xyz coordinate
   bx = -1.0*bt*np.cos(phi)-br*np.sin(phi)
   by = -1.0*bt*np.sin(phi)+br*np.cos(phi)
   bz = bz

   b = np.sqrt(bx**2+by**2+bz**2)
   # print(bx,by,bz)

   #return np.concatenate([bx, by, bz])
   return {
           'bx' :bx,
           'by' :by,
           'bz' :bz,
           'br' :br,
           'bt' :bt,
           'b'  :b,
           'psi':psi,
           'fr' :fr,
           'fz' :fz,
           'fpsi':fpsi,
           'rcbc':rc*bc,
           'rbbbs':rbbbs,
           'zbbbs':zbbbs,
           'simag':simag,
           'sibry':sibry
           }

