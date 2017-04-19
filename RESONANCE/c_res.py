import sys

sys.path.insert(0, './WEIGHT')

import matplotlib.pylab as plt
import numpy as np
import matplotlib
import matplotlib.mlab as mlab


#fn = np.load('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/PARA_RESONANCE_60_152932.npz')
Eini =60
#comment = '165037'
#comment = '159243_p'
#comment = '165042'
comment = '165037'
comment = '165037_won'
fn = np.load('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/Kathreen_'
              +np.str(Eini)+'_'+np.str(comment)+'.npz')
output = fn['output']
#fn2 = np.load('/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/FI_'+np.str(Eini)+'_'+comment+'.npz')
#output2 = fn2['output']

R0 = np.reshape(output[:,:,:,0],np.size(output)/len(output[0,0,0,:]))
Z0 = np.reshape(output[:,:,:,1],np.size(output)/len(output[0,0,0,:]))
ob = np.reshape(output[:,:,:,4],np.size(output)/len(output[0,0,0,:]))
pphi = np.reshape(output[:,:,:,5],np.size(output)/len(output[0,0,0,:]))
mu_E = np.reshape(output[:,:,:,6],np.size(output)/len(output[0,0,0,:]))
fphi = np.reshape(output[:,:,:,7],np.size(output)/len(output[0,0,0,:]))
ftheta = np.reshape(output[:,:,:,8],np.size(output)/len(output[0,0,0,:]))
rhomin = np.reshape(output[:,:,:,9],np.size(output)/len(output[0,0,0,:]))
rhomax = np.reshape(output[:,:,:,10],np.size(output)/len(output[0,0,0,:]))


# Normalize coordinate system
def normalize(data):
   data = data.astype(np.float)
   xmin = np.min(data)
   xmax = np.max(data)
   return (data - xmin) / (xmax - xmin)

def int_p(p,intp,msrho):
    mask = np.where(((ob > 3) & (mu_E>0))
                  & ((p>(intp-0.6)) & (p<(intp+0.6)))
                  & (((rhomin<msrho[1]) & (rhomin>msrho[0]))
                  | ((rhomax<msrho[1]) & (rhomax>msrho[0])))
                    )
    print(intp)
    return mask[0]



msmin,msmax = 0.,1.0
fmode = 6
nmode = 1
mmode=1
# fishbone
p = (fmode - nmode*fphi )
# AE
#p = (fmode - nmode*fphi )/ftheta - mmode

width = np.concatenate(([rhomax-msmin],[rhomax-rhomin],[np.repeat(msmax-msmin,len(rhomin))],[msmax-rhomin]),axis=0)
growth = (msmax-msmin)/np.max(width,axis=0)

plt.clf()
plt.title('$f_{mode}=$'+np.str(fmode)+'$kHz$'+'   $n=$'+np.str(nmode)+
          '  E='+np.str(Eini)+'keV'+'   '+comment)

a = np.array([],dtype='int')
#intp = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
intp = [0]
for i in intp:
    mask = int_p(p,i,[msmin,msmax])
    a = np.concatenate((a,mask))


mask = np.copy(a)
#sc2 = plt.scatter(pphi_f[0:100000],mu_E_f[0:100000],marker='+',alpha=1,color='black')
sc=plt.scatter(pphi[mask],mu_E[mask],c=p[mask],vmin=0,vmax=5,cmap=plt.cm.rainbow,marker='o')
#sc=plt.scatter(pphi[mask],mu_E[mask],color='red')
#cb = plt.colorbar(sc,orientation='horizontal')

plt.xlim(-2,9)
plt.ylim(-0.2,1.4)

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

#plt.legend((f1,f2,f3,f4,f5,f6),
#            ('L,T','L,CP','L,CTP','C,T','C,CP','C,CTP'),
#            scatterpoints=1,ncol=3,fontsize=12,loc=4)

plt.ylabel('$\mu B/E$')
plt.xlabel('$P_{\phi}/\Psi_{wall}$')

