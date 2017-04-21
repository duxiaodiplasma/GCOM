import sys
sys.path.insert(0, '/home/duxiaodi/GCOM_v3/GCOM/RESONANCE')
import matplotlib.pylab as plt
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
from scipy.stats import binned_statistic_2d
import creatobj



def read(fn):
    out = np.load(fn)
    output = out['output']
    r = creatobj.outpu(fn)
    r.R0 = np.reshape(output[:,:,:,0],np.size(output)/len(output[0,0,0,:]))
    r.Z0 = np.reshape(output[:,:,:,1],np.size(output)/len(output[0,0,0,:]))
    r.ob = np.reshape(output[:,:,:,4],np.size(output)/len(output[0,0,0,:]))
    r.pphi = np.reshape(output[:,:,:,5],np.size(output)/len(output[0,0,0,:]))
    r.mu_E = np.reshape(output[:,:,:,6],np.size(output)/len(output[0,0,0,:]))
    r.fphi = np.reshape(output[:,:,:,7],np.size(output)/len(output[0,0,0,:]))
    r.ftheta = np.reshape(output[:,:,:,8],np.size(output)/len(output[0,0,0,:]))
    r.rhomin = np.reshape(output[:,:,:,9],np.size(output)/len(output[0,0,0,:]))
    r.rhomax = np.reshape(output[:,:,:,10],np.size(output)/len(output[0,0,0,:]))
    return r

def fi(fn):
    out = np.load(fn)
    output = out['output']
    pphi_f = output[:,5]
    mu_E_f = output[:,6]

    binfi =  binned_statistic_2d(pphi_f,mu_E_f,None,'count',bins = [30,30])
    x = binfi[1][0:30]
    y = binfi[2][0:30]
    z = binfi[0]

    from scipy import interpolate
    f = interpolate.interp2d(x, y, z, kind='cubic')
    xx = np.linspace(np.min(x),np.max(x),500)
    yy = np.linspace(np.min(y),np.max(y),500)
    zz = f(xx,yy)
    return xx,yy,zz

def resonance(r,fishbone=0,AE=0,fmode=6,tolerance=0.2,m=1,n=0,rhorange=[0,1]):
    if fishbone ==1:
       phgh = ((fmode+tolerance) - n*r.fphi )
       plow = ((fmode-tolerance) - n*r.fphi )
       resonance_order = [0]

    if AE ==1:
       phgh = ((fmode+tolerance) - n*r.fphi )/r.ftheta-m
       plow = ((fmode-tolerance) - n*r.fphi )/r.ftheta-m
       resonance_order = np.arange(-5,5)

    a = np.array([],dtype='int')
    for i in resonance_order:
         mask = np.where(
                    ((r.ob > 3) & (r.mu_E>0))
                  & (((i>plow) & (i<phgh)))
                  & (((r.rhomin<rhorange[1]) & (r.rhomin>rhorange[0]))
                  | ((r.rhomax<rhorange[1]) & (r.rhomax>rhorange[0])))
                  )
         a = np.concatenate((a,mask[0]))
    return a

def plot(r,mask=None,fi=None,loss=None):

    pphi = -r.pphi
    mu_E = r.mu_E

    plt.clf()
    matplotlib.interactive('t')
    ax = plt.subplot(111)

    if loss != None:
       m1 = np.where(r.ob==1)
       f1 = ax.scatter(pphi[m1],mu_E[m1],c='purple',alpha=0.1,marker='x')

       m2 = np.where(r.ob==2)
       f2 = ax.scatter(pphi[m2],mu_E[m2],c='black',alpha=0.1,marker='x')

       m3 = np.where(r.ob==3)
       f3 = ax.scatter(pphi[m3],mu_E[m3],c ='orange',alpha=0.1,marker='x')

    m5 = np.where(r.ob==5)
    f5 = ax.plot(pphi[m5],mu_E[m5],
        alpha=0.05,color='red')

    m4 = np.where(r.ob==4)
    f4 = ax.plot(pphi[m4],mu_E[m4],
           alpha=0.4,color='green')

    m6 = np.where(r.ob==6)
    f6 = ax.plot(pphi[m6],mu_E[m6],alpha=0.2,color='blue')

    sc=ax.scatter(pphi[mask],mu_E[mask],color='black',marker='o',s=40)

    if fi != None:
       xx,yy,zz = fi(fi)
       ax.contour(xx,yy,zz)

    ax.set_xlim([-1.5,1])
    ax.set_ylim([0.,1.4])
    ax.set_xlabel('$P_{\phi}/\Psi_{wall}$')
    ax.set_ylabel('$\mu B/E$')

    return


fn = '/home/duxiaodi/GCOM_v3/GCOM/RESONANCE/output/Kathreen_165037_60keV.npz'
r = read(fn)
mask = resonance(r,fishbone=1,AE=0,fmode=6,tolerance=0.5,m=0,n=1,rhorange=[0,0.6])
plot(r,mask)
