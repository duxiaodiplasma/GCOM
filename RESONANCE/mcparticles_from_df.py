import netCDF4
import numpy as np

def read(fn_nubeam):
    fn = netCDF4.Dataset(fn_nubeam)
    r2d = fn.variables['R2D'][:]/100. # [m]
    z2d = fn.variables['Z2D'][:]/100. # [m]
    energy = fn.variables['E_D_NBI'][:]/1000. # [keV]
    pitch = fn.variables['A_D_NBI'][:]
    fi = fn.variables['F_D_NBI'][:]
    return r2d,z2d,energy,pitch,fi

def probability_constE(fi,pitch,energy,constE):

    iE = np.argmin(np.abs(energy-constE))
    prz = np.sum(fi,axis=1)[:,iE]
    prz = prz/np.sum(prz)

    dimrz = np.shape(fi)[0]
    dimP = np.shape(fi)[1]
    dimE = np.shape(fi)[2]
    #energy probability
    ppitch = np.zeros((dimrz,dimP))
    for i in np.arange(dimrz):
        ppitch[i,:] = fi[i,:,iE]/np.sum(fi[i,:,iE])

    return prz,ppitch

def generator_constE(fi,r2d,z2d,energy,pitch,prz,ppitch,E0):
    dimrz = np.shape(fi)[0]
    dimP = np.shape(fi)[1]
#    dimE = np.shape(fi)[2]

    iRZ = np.random.choice(np.arange(dimrz),p=prz)
    R0 = r2d[iRZ]
    Z0 = z2d[iRZ]
    ipitch = np.random.choice(np.arange(dimP),p=ppitch[iRZ,:])
    P0 = pitch[ipitch]

    return R0,Z0,E0,P0


def main(fn,num,constE):
    r2d,z2d,energy,pitch,fi = read(fn)
    prz,ppitch = probability_constE(fi,pitch,energy,constE)

    R0 = np.zeros(num)
    Z0 = np.zeros(num)
    E0 = np.zeros(num)
    P0 = np.zeros(num)
    for i in range(0,num):
        R0[i],Z0[i],E0[i],P0[i] = generator_constE(fi,r2d,z2d,energy,pitch,prz,ppitch,constE)

    return R0,Z0,E0,P0


