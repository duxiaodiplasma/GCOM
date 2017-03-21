import sys

sys.path.insert(0, './WEIGHT')

import matplotlib.pylab as plt
import numpy as np
import matplotlib
import matplotlib.mlab as mlab
#import bindata


fn = np.load('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/PARA_RESONANCE_30_152932.npz')
output = fn['output']

R0 = np.reshape(output[:,:,:,0],np.size(output)/len(output[0,0,0,:]))
Z0 = np.reshape(output[:,:,:,1],np.size(output)/len(output[0,0,0,:]))
ob = np.reshape(output[:,:,:,4],np.size(output)/len(output[0,0,0,:]))
pphi = np.reshape(output[:,:,:,5],np.size(output)/len(output[0,0,0,:]))
mu_E = np.reshape(output[:,:,:,6],np.size(output)/len(output[0,0,0,:]))
fphi = np.reshape(output[:,:,:,7],np.size(output)/len(output[0,0,0,:]))
ftheta = np.reshape(output[:,:,:,8],np.size(output)/len(output[0,0,0,:]))
rhomin = np.reshape(output[:,:,:,9],np.size(output)/len(output[0,0,0,:]))
rhomax = np.reshape(output[:,:,:,10],np.size(output)/len(output[0,0,0,:]))


def bin(x,y,z,nx,ny):
    xmin = np.min(x)
    ymin = np.min(y)
    xmax = np.max(x)
    ymax = np.max(y)

    xi = np.linspace(xmin,xmax,nx)
    yi = np.linspace(ymin,ymax,ny)

    xi,yi = np.meshgrid(xi,yi)

    x_new, xi_new = normalize(x), normalize(xi)
    y_new, yi_new = normalize(y), normalize(yi)

    zi = mlab.griddata(x_new, y_new, z, xi_new, yi_new,interp='linear')

    return {'xx':xi,
            'yy':yi,
            'zz':zi}

# Normalize coordinate system
def normalize(data):
   data = data.astype(np.float)
   xmin = np.min(data)
   xmax = np.max(data)
   return (data - xmin) / (xmax - xmin)

def int_p(p,intp,msrho):
    mask = np.where(((ob > 3) & (mu_E>0))
                  & ((p>intp*0.99) & (p<1.01*intp))
                  & (((rhomin<msrho[1]) & (rhomin>msrho[0]))
                  | ((rhomax<msrho[1]) & (rhomax>msrho[0])))
                    )
    return mask



msmin,msmax = 0.3,0.8
fmode = 100
nmode = -4
p = (fmode - nmode*fphi )/ftheta

width = np.concatenate(([rhomax-msmin],[rhomax-rhomin],[np.repeat(msmax-msmin,len(rhomin))],[msmax-rhomin]),axis=0)
growth = (msmax-msmin)/np.max(width,axis=0)

plt.clf()
plt.title('$f_{mode}=$'+np.str(fmode)+'$kHz$'+'   $n=$'+np.str(nmode)+
          '  E=30keV')

for i in range(-20,20):
    mask = int_p(p,i,[msmin,msmax])
    #sc = plt.scatter(pphi[mask],mu_E[mask],c=growth[mask],vmin=0.4,vmax=1.0,cmap=plt.cm.rainbow)
    #print(p[mask])
    sc = plt.scatter(pphi[mask],mu_E[mask],c=p[mask]-15,vmin=-5,vmax=5,cmap=plt.cm.rainbow)
cb = plt.colorbar(sc,orientation='horizontal')

#nu = np.load('RESONANCE/OUT/FI_NUBEAM.npz')
#nubeam = nu['output']
#fi_pphi = -1.0*np.reshape(nubeam[:,:,:,5],220*50*75)
#fi_muE  = np.reshape(nubeam[:,:,:,6],220*50*75)
#fi      = np.reshape(nubeam[:,:,:,7],220*50*75)
#fi_ob   = np.reshape(nubeam[:,:,:,4],220*50*75)
#mask_fi = np.where(fi>5e6)
#bin_fi_x = np.linspace(-1.5,0.6,20)
#bin_fi_y = np.linspace(0,1.4,20)
#binfi = bindata.TwoD(fi_pphi,fi_muE,fi,bin_fi_x,bin_fi_y,ppbin=False,binval='median')
#plt.contour(bin_fi_x,bin_fi_y,binfi,20,alpha=0.3,cmap=plt.cm.plasma)
#

matplotlib.interactive('t')
ax = plt.subplot(111)

color = [str(item/255.) for item in ob]
m1 = np.where(ob==1)
f1 = ax.scatter(pphi[m1],mu_E[m1],c='purple',alpha=0.1,marker='x')

# L,co passing
m2 = np.where(ob==2)
f2 = ax.scatter(pphi[m2],mu_E[m2],c='black',alpha=0.1,marker='x')

m3 = np.where(ob==3)
f3 = ax.scatter(pphi[m3],mu_E[m3],c ='orange',alpha=0.1,marker='x')

m4 = np.where(ob==4)
f4 = ax.scatter(pphi[m4],mu_E[m4],
    marker='s',alpha=0.1,s=20,facecolors='none',edgecolors='r')

# C,co passing
m5 = np.where(ob==5)
f5 = ax.scatter(pphi[m5],mu_E[m5],
    alpha=0.1,marker='s',s=20,facecolors='none',edgecolors='b')

m6 = np.where(ob==6)
f6 = ax.scatter(pphi[m6],mu_E[m6],alpha=0.1,marker='s',facecolors='none',edgecolors='g')

plt.legend((f1,f2,f3,f4,f5,f6),
            ('L,T','L,CP','L,CTP','C,T','C,CP','C,CTP'),
            scatterpoints=1,ncol=3,fontsize=12,loc=4)

plt.ylabel('$\mu B/E$')
plt.xlabel('$P_{\phi}/\Psi_{wall}$')

