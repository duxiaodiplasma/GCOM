import gl
import bfield
import matplotlib
matplotlib.use( "Agg" )
import matplotlib.pylab as plt
import numpy as np
import geqdsk
import ugrid
import trace
from progressbar import ProgressBar
from scipy import interpolate
from scipy.signal import savgol_filter

def f_pphi(R,bphi,b,v_para,psi):
    """ calculate canonical toroial momentum """
    pphi \
    = gl.m*R*(bphi/b)*v_para+gl.q*psi
    # not sure plus or minus
    #= gl.m*R*(bphi/b)*v_para+gl.q*psi
    return pphi


def f_mu(R,E,pphi,btotal,bphi,psi):
    """ claculate magnetic momentum """
    mu = E/btotal - (btotal/(2*gl.m))*((pphi-gl.q*psi)/(R*bphi))**2
    mu = mu/1e3/gl.q
    return mu


def f_ini(R0,Z0,E0,pitch0):
    """ initialize the inputs """
    # SI unit
    E0 = E0*gl.q*1e3

    # velocity
    v = np.sqrt(2*E0/gl.m)
    v_para = v*pitch0
    v_perp = np.sqrt(v**2-v_para**2)

    # magnetic field in the equilibirum
    b = bfield.main(R0,Z0,0,gl.fn)
    bphi = b['bt']
    btotal = b['b']
    psi   = b['psi']

    # fixed equilibrium
    f_br  = b['fr']
    f_bz  = b['fz']
    f_psi = b['fpsi']
    rcbc  = b['rcbc']
    rbbbs = b['rbbbs']
    zbbbs = b['zbbbs']
    simag = b['simag']
    sibry = b['sibry']

    # initial toroidal canonical momentum
    pphi0 = f_pphi(R0,bphi,btotal,v_para,psi)

    # intial mangetic momentum
    mu0 = E0/btotal*(1-pitch0**2)/1e3/gl.q

    return {
       	    'E0'    :E0,
            'v'     :v,
            'v_para':v_para,
            'v_perp':v_perp,
            'pphi0' :pphi0,
            'mu0'    :mu0,
            'f_br'  :f_br,
            'f_bz'  :f_bz,
            'f_psi' :f_psi,
            'rcbc'  :rcbc,
            'rbbbs' :rbbbs,
            'zbbbs' :zbbbs
            }

def f_vpara(R,btotal,bphi,pphi,psi):
    """ calculate vparallel velocity from pphi  """
    vpara = btotal/(gl.m*R*bphi)*(pphi-gl.q*psi)
    return vpara

def f_mucontour(r,z,murz,vpararz,rbbbs,zbbbs,mu0,rmaxis,R0,Z0,nseg,f_psi,sibry,simag):
    """  get the specific mu contour """

    # Do not plot the contour in X window
#    plt.ioff()

    # this guy will have orbit type ...
    orbit_class = np.zeros(1)

    # prepare to clarify the orbit type
    # interpolate v_para
    fvpara = interpolate.interp2d(z,r,vpararz, kind='cubic')

    # which mu levels to get the orbit
    levels = mu0

    # get the contour at level mu0
    cs = plt.contour(r,z,np.transpose(murz),levels=levels)

    # Get the path of the contour line
    p = cs.collections[0].get_paths()#[0]

    # dummy array for output
    tmpr = np.zeros(1)
    tmpz = np.zeros(1)
    obr = np.zeros(1)
    obz = np.zeros(1)
    ob_vpara = np.zeros(1)

    # dummy array
    allobr = np.zeros(1)
    allobz = np.zeros(1)
    minall = np.zeros(len(p))

    if len(p) >0:

       for subline in np.arange(len(p)):
           allobr = p[subline].vertices[:,0]
           allobz = p[subline].vertices[:,1]
           minall[subline] = np.min(np.sqrt((allobr-R0)**2+(allobz-Z0)**2))

       # find the index of minimum
       min_index = np.argmin(minall)

       # this gives the real poloidal projection
       obr = p[min_index].vertices[:,0]
       obz = p[min_index].vertices[:,1]

       # will have vpara in contour path
       ob_vpara_tmp = np.zeros(len(obr))
       for i in np.arange(len(obr)):
           ob_vpara_tmp[i] = fvpara(obz[i],obr[i])

       # give the orbit class
       orbit_class,rho \
       = orbit_check(ob_vpara_tmp,obr,obz,rbbbs,zbbbs,rmaxis,f_psi,sibry,simag)

       # interpolatation for fine resolution
       tck, u = interpolate.splprep([obr,obz],s=0)
       obr,obz= interpolate.splev(np.linspace(0,1,nseg),tck)

       # new vapra for fine obr, obz
       ob_vpara = np.zeros(len(obr))
       for i in np.arange(len(obr)):
           ob_vpara[i] = fvpara(obz[i],obr[i])
# debug        print(1/2.*gl.m*ob_vpara[i]**2/gl.q/1e3)
#
    return {
            'obr'         : obr,
            'obz'         : obz,
            'orbit_class' : orbit_class,
            'ob_vpara'    : ob_vpara,
            'fvpara'      : fvpara,
            'levels'      : levels,
            'rho'         : rho
           }


def orbit_check(ob_vpara,obr,obz,rbbbs,zbbbs,rmaxis,f_psi,sibry,simag):
    """ identify the orbit type from vpara """
    # Lost  -  Trapped             1
    # Lost  -  counter-Ip Passing  2
    # lost  -  co-Ip Passing       3
    # Conf. -  Trapped             4
    # Conf. -  counter-Ip Passing  5
    # Conf. -  co-Ip Passing       6

#    outboundary =  \
#    [obr>np.max(rbbbs), obr<np.min(rbbbs), \
#     obz>np.max(zbbbs), obz<np.min(zbbbs)]
    rho = trace.RZ_to_RHO(obr,obz,f_psi,simag,sibry)
    if np.max(np.abs(rho)) > 1:
       outboundary = True
    else:
       outboundary = False

    if outboundary: \
       # hit the wall, lost particles
       if np.any(np.sign(ob_vpara) ==1) and np.any(np.sign(ob_vpara) ==-1):
          # lost, trapped particles
          orb = 1
       elif np.any(np.sign(ob_vpara)==1):
          # lost, counter-Ip particles
          orb = 2
       elif np.any(np.sign(ob_vpara)==-1):
          # lost, co-Ip passing particles
          orb = 3
    else:
       if np.any(np.sign(ob_vpara) ==1) and np.any(np.sign(ob_vpara) ==-1):
          # Confined trapped particle, bounced back in v_para
          orb = 4
       elif np.any(np.sign(ob_vpara)==1):
          # Confined counter-Ip particle, always one direction in v_para
          orb = 5
       elif np.any(np.sign(ob_vpara)==-1):
          # Confined co-Ip passing particle, always one direction in v_para
          orb = 6

    return orb,rho

def ob_stream(obr,obz):

    """ calculate the unit vector along the orbit
        and the distance between every two neightbor points """

    # these guys will have ...
    n = len(obr) # dimention of the orbit
    flow = np.zeros((n,3))  # direction vector between two neighbour points
    unit_flow = np.zeros((n,3)) # unit vector between two neighbour points
    s_flow = np.zeros(n) # distance between every two neighbouring points

    for i in range(0,n):
        if i<n-1:
            flow[i,0] = obr[i+1]-obr[i]
            flow[i,2] = obz[i+1]-obz[i]
        else:
            flow[i,0] = obr[0]-obr[i]
            flow[i,2] = obz[0]-obz[i]

    for i in range(0,n):
        s_flow[i] = np.linalg.norm(flow[i,:])
        unit_flow[i,:] = flow[i,:]/s_flow[i]

    tmp_flow = np.zeros((n+1,3))
    tmp_flow[:,0] = np.append(unit_flow[:,0],unit_flow[0,0])
    tmp_flow[:,2] = np.append(unit_flow[:,2],unit_flow[0,2])

    return{
        'flow'  : flow,
        's_flow': s_flow, # intermitent distance between points in orbit
        'unit_flow':unit_flow,
        }

def py_ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2' """
    cosang = np.dot(v1, v2)
    sinang = np.linalg.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

