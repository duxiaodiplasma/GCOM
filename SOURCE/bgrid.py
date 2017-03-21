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

def main(fn):

   # read gfile
   res = geqdsk.read(fn)

   nr = res['nw']
   nz = res['nh']
   r = res['r']
   z = res['z']
   psirz = res['psirz']
   rc = res['rcentr']
   bc = res['bcentr']

   # ---------
   rr = np.zeros(nr)
   zz = np.zeros(nz)
   for i in range(0, nr):
       rr[i] = r[i, 0]
   for i in range(0, nz):
       zz[i] = z[0, i]

   f_psirz = interpolate.interp2d(zz, rr, psirz, kind='cubic')
#   #
   new_nr = 204
   new_nz = 204
   new_rr = np.linspace(np.min(rr),np.max(rr),new_nr)
   new_zz = np.linspace(np.min(zz),np.max(zz),new_nz)
   new_psirz = np.zeros((new_nr,new_nz))
   new_psirz = f_psirz(new_zz,new_rr)
   new_r = np.tile(new_rr,new_nr).reshape((new_nr,new_nz)).transpose()
   new_z = np.tile(new_zz,new_nz).reshape((new_nr,new_nz))

   nr = np.copy(new_nr)
   nz = np.copy(new_nz)
   rr = np.copy(new_rr)
   zz = np.copy(new_zz)
   r  = np.copy(new_r)
   z  = np.copy(new_z)
   psirz = np.copy(new_psirz)

   # ---------


   # first deravitive of the psi in RZ direction
   dpsi_dr = np.zeros((nr,nz))
   dpsi_dz = np.zeros((nr,nz))

   br = np.zeros((nr,nz))
   bz = np.zeros((nr,nz))
   bt = np.zeros((nr,nz))
   b  = np.zeros((nr,nz))

   dbz_dr = np.zeros((nr,nz))
   dbz_dz = np.zeros((nr,nz))
   dbr_dr = np.zeros((nr,nz))
   dbr_dz = np.zeros((nr,nz))
   dbt_dr = np.zeros((nr,nz))
   dbt_dz = np.zeros((nr,nz))
   db_dr  = np.zeros((nr,nz))
   db_dz  = np.zeros((nr,nz))

   for i in range(1, nr-1):
       for j in range(1, nz-1):
           dpsi_dr[i,j] = (psirz[i+1,j]-psirz[i-1,j])/(r[i+1,j]-r[i-1,j])
           dpsi_dz[i,j] = (psirz[i,j+1]-psirz[i,j-1])/(z[i,j+1]-z[i,j-1])

   for i in range(1, nr-1):
       for j in range(1, nz-1):
           br[i,j]      = -1.0*-(1/r[i,j])*dpsi_dz[i,j]
           bz[i,j]      = -1.0* (1/r[i,j])*dpsi_dr[i,j]
           bt[i,j]      = rc*bc/r[i,j]
           b[i,j]       = np.sqrt(br[i,j]**2+bz[i,j]**2+bt[i,j]**2)

   for i in range(2, nr-2):
       for j in range(2, nz-2):
           dbr_dr[i,j]  = (br[i+1,j  ]-br[i-1,j  ])/(r[i+1,j  ]-r[i-1,j  ])
           dbt_dr[i,j]  = (bt[i+1,j  ]-bt[i-1,j  ])/(r[i+1,j  ]-r[i-1,j  ])
           dbz_dr[i,j]  = (bz[i+1,j  ]-bz[i-1,j  ])/(r[i+1,j  ]-r[i-1,j  ])

           dbr_dz[i,j]  = (br[i  ,j+1]-br[i  ,j-1])/(z[i  ,j+1]-z[i  ,j-1])
           dbt_dz[i,j]  = (bt[i  ,j+1]-bt[i  ,j-1])/(z[i  ,j+1]-z[i  ,j-1])
           dbz_dz[i,j]  = (bz[i  ,j+1]-bz[i  ,j-1])/(z[i  ,j+1]-z[i  ,j-1])

           db_dr[i,j]  = ( b[i+1,j  ]- b[i-1,j  ])/(r[i+1,j  ] -r[i-1,j  ])
           db_dz[i,j]  = ( b[i  ,j+1]- b[i  ,j-1])/(z[i  ,j+1] -z[i  ,j-1])

   # reshape the mesh

   r     = np.copy(r[2:nr-2,2:nz-2])
   z     = np.copy(z[2:nr-2,2:nz-2])
   rr    = np.copy(rr[2:nr-2])
   zz    = np.copy(zz[2:nz-2])
   psirz = np.copy(psirz[2:nr-2,2:nz-2])


   br = np.copy(br[2:nr-2,2:nz-2])
   bz = np.copy(bz[2:nr-2,2:nz-2])
   bt = np.copy(bt[2:nr-2,2:nz-2])
   b  = np.copy( b[2:nr-2,2:nz-2])

   dbz_dr = np.copy(dbz_dr[2:nr-2,2:nz-2])
   dbz_dz = np.copy(dbz_dz[2:nr-2,2:nz-2])

   dbr_dr = np.copy(dbr_dr[2:nr-2,2:nz-2])
   dbr_dz = np.copy(dbr_dz[2:nr-2,2:nz-2])

   dbt_dr = np.copy(dbt_dr[2:nr-2,2:nz-2])
   dbt_dz = np.copy(dbt_dz[2:nr-2,2:nz-2])

   db_dr  = np.copy( db_dr[2:nr-2,2:nz-2])
   db_dz  = np.copy( db_dz[2:nr-2,2:nz-2])

   # new dimension
   nr,nz = np.shape(db_dr)

   # 2D interpolation for everything

   f_br = interpolate.interp2d(zz, rr, br, kind='cubic')
   f_bz = interpolate.interp2d(zz, rr, bz, kind='cubic')
   f_bt = interpolate.interp2d(zz, rr, bt, kind='cubic')
   f_b  = interpolate.interp2d(zz, rr, b,  kind='cubic')

   f_dbz_dr = interpolate.interp2d(zz, rr, dbz_dr, kind='cubic')
   f_dbz_dz = interpolate.interp2d(zz, rr, dbz_dz, kind='cubic')

   f_dbr_dr = interpolate.interp2d(zz, rr, dbr_dr, kind='cubic')
   f_dbr_dz = interpolate.interp2d(zz, rr, dbr_dz, kind='cubic')

   f_dbt_dr = interpolate.interp2d(zz, rr, dbt_dr, kind='cubic')
   f_dbt_dz = interpolate.interp2d(zz, rr, dbt_dz, kind='cubic')

   f_db_dr =  interpolate.interp2d(zz, rr, db_dr, kind='cubic')
   f_db_dz =  interpolate.interp2d(zz, rr, db_dz, kind='cubic')


   return {
           'br' : br,
           'bz' : bz,
           'bt' : bt,
           'b'  : b,

           'f_br' : f_br,
           'f_bz' : f_bz,
           'f_bt' : f_bt,
           'f_b'  : f_b,

           'dbz_dr' : dbz_dr,
           'dbz_dz' : dbz_dz,

           'dbr_dr' : dbr_dr,
           'dbr_dz' : dbr_dz,

           'dbt_dr' : dbt_dr,
           'dbt_dz' : dbt_dz,

           'db_dr'  : db_dr,
           'db_dz'  : db_dz,

           'f_dbz_dr':f_dbz_dr,
           'f_dbz_dz':f_dbz_dz,

           'f_dbr_dr':f_dbr_dr,
           'f_dbr_dz':f_dbr_dz,

           'f_dbt_dr':f_dbt_dr,
           'f_dbt_dz':f_dbt_dz,

           'f_db_dr' :f_db_dr,
           'f_db_dz' :f_db_dz,

           'rr':rr,
           'zz':zz,
           'r' :r,
           'z' :z,

           'psirz':psirz,

           'nr':nr,
           'nz':nz
           }

