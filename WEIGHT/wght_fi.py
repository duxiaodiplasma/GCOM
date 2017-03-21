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
import matplotlib.pylab as plt

def rebin(x,y,z,binx,biny,xnew,ynew):

    # bin irregular datas to regular grid binx, biny
    mask_x = np.digitize(x,binx)-1
    mask_y = np.digitize(y,biny)-1

    # Z value correpsonding to binXY
    binz = np.zeros((len(binx),len(biny)))
    for i in range(0,len(x)):
        binz[mask_x[i],mask_y[i]] \
        = z[i] + binz[mask_x[i],mask_y[i]]

    # interpolate data for any desired output
    f = interpolate.interp2d(binx,biny,binz.transpose())

    znew = f(xnew,ynew)

    return znew


"""THE CODE CALCULATE THE PROBABILITY OF ORBIT AT THE GRID POINT"""

fn_fida = './INOUT/159243.inpa801_npa_weights.h5'
fn_grid = './INOUT/gridf.h5'
out_fidasim = wght_ini.main(fn_fida,fn_grid)

cx  = out_fidasim['cx']
att = out_fidasim['att']
gz  = out_fidasim['z']
ngx = out_fidasim['nx']
ngy = out_fidasim['ny']
ngz = out_fidasim['nz']
En  = out_fidasim['E']
pitch = out_fidasim['pitch']
nEn = En.shape[0]

fn_ob = "./INOUT/OUT/ch1_ssnpa.wght"
out_ob = np.loadtxt(fn_ob)

E = out_ob[:,7]  # Energy
P = out_ob[:,12] # pitch
W = out_ob[:,9]  #
R = out_ob[:,11]

# energy 21keV to 79keV
bin_E = np.copy(En[10:40])

# pitch -1 to 1
bin_P = np.linspace(-1,1,50)

bin_R = np.linspace(1.0,2.5,100)


outfi_Rmax = np.load('./INOUT/PARA_NU.npz')

E_new = np.linspace(np.min(bin_E),np.max(bin_E),100)
P_new = np.linspace(-1,1,100)
R_new = np.linspace(np.min(bin_R),np.max(bin_R),100)

W_EP = rebin(E,P,W,bin_E,bin_P,E_new,P_new)
W_RP = rebin(R,P,W,bin_R,bin_P,R_new,P_new)
W_ER = rebin(E,R,W,bin_E,bin_R,E_new,R_new)

plt.clf()
f,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,5))
ax1.contourf(E_new,P_new,W_EP,100,cmap=plt.cm.plasma)
ax2.contourf(R_new,P_new,W_RP,100,cmap=plt.cm.plasma)
ax3.contourf(E_new,R_new,W_ER,100,cmap=plt.cm.plasma)

ax1.set_xlabel('$E [keV]$')
ax1.set_ylabel('$(v_{\parallel}/v)_{mid} $')
ax2.set_xlabel('$R_{max} [m]$ at midplane')
ax2.set_ylabel('$(v_{\parallel}/v)_{mid}$')
ax3.set_xlabel('$E [keV]$')
ax3.set_ylabel('$R_{max} [m]$ at midplane')


