import h5py as h5
import numpy as np


def main(fn_fida,fn_grid):

    """THIS SUBROUTINE READ IN FIDASIM OUPUT AND CALCULATE ORBIT WEIGHT"""

    """
    READ IN
    INCLUDE BEAMGRID,CHORD,PITCH
    INCLUDE PROBABILITY FOR CX ATTENUATION
    """

    fgrid = h5.File(fn_grid)

    nx = fgrid['nx'].value
    ny = fgrid['ny'].value
    nz = fgrid['nz'].value

    # cm -> m
    x = fgrid['x_grid'].value/100.
    y = fgrid['y_grid'].value/100.
    z = fgrid['z_grid'].value/100.

    fp = h5.File(fn_fida)

    E = fp['energy'].value

    #In [141]: fp['phit'].attrs.values()
    #Out[141]: ['Probability of hitting the detector given an isotropic source: phit(chan,z,y,x)']
    phit = fp['phit'].value

    #In [142]: fp['cx'].attrs.values()
    #Out[142]: ['Charge-exchange rate: cx(energy,x,y,z,chan)', 's^(-1)']
    # beam charge-exchange rate ; you want to increase this probability to have more neutrals generated.
    # ??wait a second... is there factor 1e14?
    cx   = fp['cx'].value

    #In [10]: fp['attenuation'].attrs.values()
    #Out[10]: ['Attenuation factor i.e. survival probability: attenuation(energy,x,y,z,chan)']
    #Survival probability for neutrals to reach detectors
    att  = fp['attenuation'].value

    #In [46]: fp['pitchxd'].attrs.values()
    #Out[46]: ['pitch for EP to enter detector from each grid']
    pitch = fp['pitchxd'].value

    return{
          # BEAM GRID
          'nx' : nx,
          'ny' : ny,
          'nz' : nz,
          'x'  :  x,
          'y'  :  y,
          'z'  :  z,

          # FIDASIM OUTPUT
          'E'   :   E,
          'phit': phit,
          'cx'  :  cx,
          'att' : att,
          'pitch' : pitch
            }


