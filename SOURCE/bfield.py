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

import numpy as np

def main(lr,lz,phi,f_br,f_bz,f_bt,f_psi):

   # right hand?
   br = f_br(lz,lr)
   bz = f_bz(lz,lr) # -1.0 for unify the defination
			     # of the plus current direction in DIII-D,$
			     # which is CCW direction for plus

   # -1 for right-handed coordinate
   bt = f_bt(lz,lr)
   # print(br,bz,bt)
   psi = f_psi(lz,lr)

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
           }

