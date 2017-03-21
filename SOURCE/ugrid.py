import bgrid
import numpy as np
import gl
from scipy import interpolate

def equ(fn):

    out = bgrid.main(fn)

    nr = out['nr']
    nz = out['nz']

    br = out['br']
    bz = out['bz']
    bt = out['bt']
    b  = out['b']

    r  = out['r']
    z  = out['z']
    rr = out['rr']
    zz = out['zz']

    dbz_dr = out['dbz_dr']
    dbz_dz = out['dbz_dz']

    dbr_dr = out['dbr_dr']
    dbr_dz = out['dbr_dz']

    dbt_dr = out['dbt_dr']
    dbt_dz = out['dbt_dz']

    db_dr = out['db_dr']
    db_dz = out['db_dz']

    equ_curve = np.zeros((nr,nz,3))
    equ_grad  = np.zeros((nr,nz,3))


    # R
    equ_curve[:,:,0] = 1/(b**4)* \
         (bt*(br*dbz_dr+bz*dbz_dz) - bz*(br*dbt_dr+bz*dbt_dz +bt*br/r))
    # phi
    equ_curve[:,:,1] = 1/(b**4)* \
         (bz*(br*dbr_dr+bz*dbr_dz) - br*(br*dbz_dr+bz*dbz_dz))
    # Z
    equ_curve[:,:,2] = 1/(b**4)* \
         (br*(br*dbt_dr+bz*dbt_dz) - bt*(br*dbr_dr+bz*dbr_dz-bt**2/r))

    # R
    equ_grad[:,:,0] = 1/(b**3)*(bt*db_dz)
    # phi
    equ_grad[:,:,1] = 1/(b**3)*(bz*db_dr - br*db_dz)
    # Z
    equ_grad[:,:,2] = 1/(b**3)*(-1.0*bt*db_dr)


    return{
          'equ_curve' : equ_curve,
          'equ_grad'  : equ_grad,

          'rr':rr,
          'zz':zz,
          'r' : r,
          'z' : z
          }

def main(nr,nz,rr,zz,Epara,Eperp,equ_curve,equ_grad):
        #nr,nz,RR,ZZ,Epara,Eperp,equ_curve,equ_grad

    u       = np.zeros((nr,nz,3))
    u_curve = np.zeros((nr,nz,3))
    u_grad  = np.zeros((nr,nz,3))

    u_curve[:,:,0] = 2*Epara*equ_curve[:,:,0]/gl.q  # R
    u_curve[:,:,1] = 2*Epara*equ_curve[:,:,1]/gl.q  # phi
    u_curve[:,:,2] = 2*Epara*equ_curve[:,:,2]/gl.q  # Z

    u_grad[:,:,0]  = Eperp*equ_grad[:,:,0]/gl.q # R
    u_grad[:,:,1]  = Eperp*equ_grad[:,:,1]/gl.q # phi
    u_grad[:,:,2]  = Eperp*equ_grad[:,:,2]/gl.q # Z

    u[:,:,0] = u_curve[:,:,0] + u_grad[:,:,0] # uR
    u[:,:,1] = u_curve[:,:,1] + u_grad[:,:,1] # uphi
    u[:,:,2] = u_curve[:,:,2] + u_grad[:,:,2] # uZ


    # R
    f_u1 = interpolate.interp2d(zz, rr, u[:,:,0], kind='cubic') # R
    # phi
    f_u2 = interpolate.interp2d(zz, rr, u[:,:,1], kind='cubic') # phi
    # Z
    f_u3 = interpolate.interp2d(zz, rr, u[:,:,2], kind='cubic') # Z

    return{
        'u_curve' : u_curve,
        'u_grad'  : u_grad,
        'u'       : u,

        'f_u1'    : f_u1,
        'f_u2'    : f_u2,
        'f_u3'    : f_u3

    #    'f_u_curve0' : f_u_curve0,
    #    'f_u_curve1' : f_u_curve1,
    #    'f_u_curve2' : f_u_curve2,

    #    'f_u_grad0' : f_u_grad0,
    #    'f_u_grad1' : f_u_grad1,
    #    'f_u_grad2' : f_u_grad2
        }













