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
import netCDF4
import bindata


### -- read in orbit weight for NPA --
fn_fida = './IN/159243.inpa801_npa_weights.h5'
fn_grid = './IN/beamgrid.h5'
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

fn_ob = "./OUT/ch1_ssnpa.wght"
out_ob = np.loadtxt(fn_ob)

E = out_ob[:,7]  # Energy
P = out_ob[:,12] # pitch
W = out_ob[:,9]
# normalized Weight ??
#W = np.copy(W/np.sum(W))
R = out_ob[:,11]




# -- read in distribution in orbit coordinate --
fnfi = np.load('./OUT/PARA_FI_RMAX.npz')
out_fi = fnfi['output']
ob_fi = np.reshape(out_fi[:,:,:,0],220*50*75)
R_fi  = np.reshape(out_fi[:,:,:,1],220*50*75)
P_fi  = np.reshape(out_fi[:,:,:,2],220*50*75)
E_fi  = np.reshape(out_fi[:,:,:,3],220*50*75)
W_fi  = np.reshape(out_fi[:,:,:,4],220*50*75)

mask = np.where(ob_fi>0)
ob_fi = np.copy(ob_fi[mask])
R_fi = np.copy(R_fi[mask])
P_fi = np.copy(P_fi[mask])
E_fi = np.copy(E_fi[mask])
W_fi = np.copy(W_fi[mask])

mask_E = np.where((E_fi<80) & (E_fi>60))
mask_R = np.where((R_fi<2.0) &(R_fi>1.8))
mask_P = np.where((P_fi<1) & (P_fi>0.5))

# -- gives the grids for bin3d --
bin_P = np.arange(-1,1.05,0.05)
bin_R = np.arange(1.6,2.45,0.05)
bin_E = np.arange(0,90,5)


bin_W0 =bindata.TwoD(E,P,W,bin_E,bin_P,ppbin=False,binval='median')

bin_W = bindata.TwoD(R_fi[mask_E],P_fi[mask_E],W_fi[mask_E],bin_R,bin_P,ppbin=False,binval='median')
bin_W2 = bindata.TwoD(E_fi[mask_R],P_fi[mask_R],W_fi[mask_R],bin_E,bin_P,ppbin=False,binval='median')
bin_W3 = bindata.TwoD(E_fi[mask_P],R_fi[mask_P],W_fi[mask_P],bin_E,bin_R,ppbin=False,binval='median')
# -- bin everything to above grids using bin3d

f,axarr = plt.subplots(2,3,figsize=(15,10))
axarr[0,0].contourf(bin_E,bin_P,bin_W0,100,cmap=plt.cm.gnuplot2)
#axarr[0,0].scatter(E,P,c=W,cmap=plt.cm.gnuplot2)
#axarr[0,1].contourf(bin_E,bin_R,np.sum(bin_ob,axis=1).transpose(),100,cmap=plt.cm.gnuplot2)
#axarr[0,2].contourf(bin_P,bin_R,np.sum(bin_ob,axis=0).transpose(),100,cmap=plt.cm.gnuplot2)

axarr[1,0].contourf(bin_E,bin_R,bin_W3,100,cmap=plt.cm.gnuplot2)
axarr[1,1].contourf(bin_E,bin_P,bin_W2,100,cmap=plt.cm.gnuplot2)
axarr[1,2].contourf(bin_R,bin_P,bin_W,100,cmap=plt.cm.gnuplot2)

axarr[0,0].set_xlabel('$E [keV]$')
axarr[0,0].set_ylabel('$(v_{\parallel}/v)_{mid} $')
axarr[0,1].set_xlabel('$E [keV]$')
axarr[0,1].set_ylabel('$R_{max} [m]$')
axarr[0,2].set_xlabel('$(v_{\parallel}/v)_{mid}$')
axarr[0,2].set_ylabel('$R_{max} [m]$')

#axarr[1,0].set_xlabel('$E [keV]$')
#axarr[1,0].set_ylabel('$(v_{\parallel}/v)_{mid} $')
#axarr[1,1].set_xlabel('$E [keV]$')
#axarr[1,1].set_ylabel('$R_{max} [m]$')
#axarr[1,2].set_xlabel('$(v_{\parallel}/v)_{mid}$')
#axarr[1,2].set_ylabel('$R_{max} [m]$')
