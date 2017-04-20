
import gl
import bfield
import numpy as np
import geqdsk
import ugrid
import com
from scipy.integrate import odeint


def main(
        R0,Z0,E0,pitch0,
        nseg,nstep,tstep,
        switch_full_orbit,
        calc_prob,
        R,Z,RR,ZZ,
        br,bz,bt,b,
        bc,rmaxis,sibry,simag,rbbbs,zbbbs,
        nr,nz,psi,
        f_br,f_bz,f_bt,f_b,f_psi,
        equ_curve,equ_grad
        ):
    """  trace the orbit  """
    # ---- calculation ----
    ini = com.f_ini(R0,Z0,E0,pitch0)
    # all in SI unit
    pphi0 = ini['pphi0'] # ini. canonical toroidal momentum
    E0    = ini['E0']    # ini. energy [SI unit]
    # mu0 in unit of keV
    mu0   = ini['mu0']   # ini. magnetic momentum

    murz = com.f_mu(R,E0,pphi0,b,bt,psi)
    vpararz = com.f_vpara(R,b,bt,pphi0,psi)
    Epara = 1/2.*gl.m*vpararz**2
    Eperp = E0 - Epara
    ob_contour \
    = com.f_mucontour(RR,ZZ,murz,vpararz,rbbbs,zbbbs,mu0,rmaxis,R0,Z0,nseg,f_psi,sibry,simag)

    # normalized pphi
    pphi = (pphi0/gl.q-simag)/(sibry-simag)
    # normalized mu
    mu_E = mu0*np.abs(bc)/(E0/1e3/gl.q)
    ob   = ob_contour['orbit_class']
    fvpara = ob_contour['fvpara']

    obr = ob_contour['obr']
    obz = ob_contour['obz']
    rho = ob_contour['rho']

    ob_vpara = ob_contour['ob_vpara']
    pitch = ob_vpara/(np.sqrt(2*E0/gl.m))

    # these guys will have ...
    vperp = np.zeros((nstep,3)) # perp. velocity
    vpara = np.zeros((nstep,3)) # para. velocity
    vgc   = np.zeros((nstep,3)) # guiding center motion

    brzt  = np.zeros((nstep,3))     # B field (vec)
    unit_brzt = np.zeros((nstep,3)) # unit B field (vec)


    obxyz = np.zeros((np.int(nstep),3))# orbit for pts
    # this guy will have real steps
    steps = 0

    sol = 0
    tsol = 0
    f_phi = 0
    f_theta = 0

    # if the particle is not lost and want the full orbit
    if switch_full_orbit==1 and ob>3:

       # drift velocity perpendicular to B
       u_grid = ugrid.main(nr,nz,RR,ZZ,Epara,Eperp,equ_curve,equ_grad)
       f_u1 = u_grid['f_u1']
       f_u2 = u_grid['f_u2']
       f_u3 = u_grid['f_u3']

       # INITIAL CONDITION FOR INTEGRATOR
       ini_y0 = [R0,0,Z0]

       # time slices for integrator
       tsol = np.arange(0,tstep*nstep,tstep)

       # ODE INTEGRATOR
       sol = odeint(ode_equ, ini_y0, tsol, args=(f_u1,f_u2,f_u3,f_br,f_bt,f_bz,fvpara))

       # transfer RPHIZ coordinate to XYZ coordinate
       obxyz[:,0] = sol[:,0]*np.cos(sol[:,1])
       obxyz[:,1] = sol[:,0]*np.sin(sol[:,1])
       obxyz[:,2] = sol[:,2]


       # USED IN WEIGHT FUNCTION CALCULATION
       # calculate probability for the orbit
       if calc_prob == 1:
          # get index where the orbit accomplish one full poloidal projection
          index = full_proj_pol(ob,obr,obz,ob_vpara,sol)

          # if the orbit has not finished one full poloidal projection
          # trace again with nstep+200000
          if index == -1:
             nstep = nstep+500000
             obxyz = np.zeros((np.int(nstep),3))# orbit for pts
             tsol = np.arange(0,tstep*nstep,tstep)
             sol = odeint(ode_equ, ini_y0, tsol, args=(f_u1,f_u2,f_u3,f_br,f_bt,f_bz,fvpara))
             obxyz[:,0] = sol[:,0]*np.cos(sol[:,1])
             obxyz[:,1] = sol[:,0]*np.sin(sol[:,1])
             obxyz[:,2] = sol[:,2]
             index = full_proj_pol(ob,obr,obz,ob_vpara,sol)

          # total index for one poloidal projection
          steps = index

          # toroidal precession frequency
          # in unit of kHz
          f_phi = sol[steps,1]/tsol[steps]/1e3/(2.*np.pi)

          # transit freqency/ bounce frequency
          # in unit of kHz
          f_theta = 1./tsol[steps]/1e3

          # reshape the orbit array
          # to save the computing time...
          obxyz = np.copy(obxyz[0:index+1,:])


    # OUTPUT BELOW
    if switch_full_orbit == 0:

       return {
           'pphi' :  pphi,
           'mu_E' :  mu_E,
           'ob'   :  ob,
           'obr'  :  obr,
           'obz'  :  obz,
           'pitch': pitch,
           'Rmax' : np.max(obr),
           'Pmax' : pitch[np.argmax(obr)],
           'Zmax' : obz[np.argmax(obr)],
           'rho'  : rho
           }

    elif switch_full_orbit ==1:

       return {
           'obxyz': obxyz,
           'obr'  : obr,
           'obz'  : obz,
           'pitch': pitch,
           'Rmax' : np.max(obr),
           'Pmax' : pitch[np.argmax(obr)],
           'Zmax' : obz[np.argmax(obr)],
           'rho'  : rho,

           'pphi' :  pphi,
           'mu_E' :  mu_E,
           'ob'   :  ob,
           'steps':  steps,

           # OUTPUT OF FREQUENCY
           'sol'  : sol,
           'tsol' : tsol,
           'f_phi':f_phi,
           'f_theta':f_theta,
           }

# PREPARE ODE EQUATION for ODE INTEGRATOR
def ode_equ(y,t,f_u1,f_u2,f_u3,f_br,f_bt,f_bz,fvpara):
    """ODE EQUATION having dydt = f(t,y)"""

    vperp = np.zeros(3)
    brzt = np.zeros(3)

    R_i0 = y[0]
    Z_i0 = y[2]

    # perpendicular drift velocity in RPHIZ coordinate
    vperp[0] = f_u1(Z_i0,R_i0) # R
    vperp[1] = f_u2(Z_i0,R_i0) # phi
    vperp[2] = f_u3(Z_i0,R_i0) # Z

    # magnetic fieed vecotr in RPHIZ coordinate
    brzt[0]  = f_br(Z_i0,R_i0) # br
    brzt[1]  = f_bt(Z_i0,R_i0) # bphi
    brzt[2]  = f_bz(Z_i0,R_i0) # bz
    unit_brzt = brzt/np.linalg.norm(brzt)

    # parallel velocity in RPHIZ coordinate
    vpara = -1.0*fvpara(Z_i0,R_i0)*unit_brzt

    # guiding center motion
    # drift motion + parallel velocity
    v = vpara+vperp
    # PHI
    v[1] = v[1]/R_i0

    return v


def full_proj_pol(ob,obr,obz,ob_vpara,sol):

    # find dummy center
    if ob == 4:
       mask1 = np.where(ob_vpara>0)
       mask2 = np.where(ob_vpara<0)

       cent = (np.max(obr[mask1])+(np.max(obr[mask2])))/2.
    elif ob > 4 :
       cent = (np.min(obr)+np.max(obr))/2.

    # calculate the angle
    angle = np.angle(sol[:,0]-cent+1j*sol[:,2],deg=True)

    # find next period (jump from 180 to -180 degree)
    mask = np.where(np.abs(np.diff(angle))>180)[0]

    if np.size(mask) > 1:
       # In the next period, find the minimum difference value from index
       index = mask[0]+1+np.argmin(np.abs(angle[mask[0]+1:mask[1]]-angle[0]))
    else:
       print('not finish one full poloidal projection')
       index = -1

    return index



def probability(R0,Z0,steps,obxyz):

    if steps!= 0:
       R = np.sqrt(obxyz[:,0]**2+obxyz[:,1]**2)
       Z = np.copy(obxyz[:,2])
       n = np.copy(steps)

       out = com.ob_stream(R,Z)

       # s: distance
       s = out['s_flow']

       # total grid number
       # dummy value
       ngrid = np.sum(s)

       # fraction in total length
       sf = s/np.sum(s)

       # fractions in total grids
       sf = sf*ngrid

       # probability per grid
       p = 1./steps/sf

       p = np.copy(p[0:len(p)-1])
       s = np.copy(s[0:len(p)-1])
       sf = np.copy(sf[0:len(p)-1])

    if steps==0:
       # tracing failed! 'bad' calculation!
       # put dummy values for output
       p = np.array([0])
       s = np.array([0])
       sf = np.array([0])

    return{
        'p'  :  p,
        's'  :  s,
        'sf' :  sf,
            }

def RZ_to_RHO(obr,obz,f_psi,simag,sibry):
    psi = np.zeros(len(obr))
    rho = np.copy(psi)
    for i in range(0,len(obz)):
        psi[i] = f_psi(obz[i],obr[i])
        rho[i] = np.sqrt((psi[i]-simag)/(sibry-simag))

    return rho






















