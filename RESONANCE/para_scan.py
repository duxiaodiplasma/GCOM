import sys
import numpy as np
from joblib import Parallel, delayed
sys.path.insert(0, '/home/duxiaodi/GCOM_v3/GCOM/SOURCE')
import creatobj
import geqdsk
import bgrid
import ugrid
import trace

def constE2(g,inpu,Rstep=0.01,Pstep=0.01,E0=60):
    R_array = np.arange(np.min(g.rbbbs),np.max(g.rbbbs),Rstep)
    pitch_array = np.arange(-1,1.05,Pstep)
    E_array = E0
    Z_array = 0

    rv,zv,ev,pv = np.meshgrid(R_array,Z_array,E_array,pitch_array)
    inpu.R_array = np.reshape(rv,rv.size)
    inpu.Z_array = np.reshape(zv,zv.size)
    inpu.E_array = np.reshape(ev,ev.size)
    inpu.pitch_array = np.reshape(pv,pv.size)
    inpu.scan_steps = len(inpu.R_array)
    inpu.E0 = E0
    return inpu

def constE(g,inpu,Rstep=0.05,Zstep=0.1,Pstep=0.05,E0=60):
    R_array = np.arange(np.min(g.rbbbs),np.max(g.rbbbs),Rstep)
    Z_array = np.arange(np.min(g.zbbbs),np.max(g.zbbbs),Zstep)
    E_array = E0
    pitch_array = np.arange(-1,1.05,Pstep)

    rv,zv,ev,pv = np.meshgrid(R_array,Z_array,E_array,pitch_array)
    inpu.R_array = np.reshape(rv,rv.size)
    inpu.Z_array = np.reshape(zv,zv.size)
    inpu.E_array = np.reshape(ev,ev.size)
    inpu.pitch_array = np.reshape(pv,pv.size)
    inpu.scan_steps = len(inpu.R_array)
    inpu.E0 = E0
    return inpu

def main(g,inpu,outpu,i):

    inpu.R0 = inpu.R_array[i]
    inpu.Z0 = inpu.Z_array[i]
    inpu.pitch0 = inpu.pitch_array[i]
    inpu.E0 = inpu.E_array[i]
    inpu.phi0 = 0

    out = trace.main(g,inpu,outpu)

    output = np.zeros(11)
    output[0]  = inpu.R0
    output[1]  = inpu.Z0
    output[2]  = inpu.pitch0
    output[3]  = inpu.E0
    output[4]  = outpu.ob[0]
    output[5]  = outpu.pphi[0]
    output[6]  = outpu.mu_E[0]
    output[7]  = outpu.f_phi[0]
    output[8]  = outpu.f_theta[0]
    output[9]  = np.min(outpu.rho[0])
    output[10] = np.max(outpu.rho[0])

    return output











