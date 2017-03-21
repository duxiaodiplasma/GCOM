import numpy as np
import gl
import prepare
import hickle as hkl
from scipy import interpolate

# prepare everything necessay
fn1  = prepare.save(gl.fn)

out = hkl.load(fn1)
equ   = out[0]
bmesh = out[1]
umesh = out[2]

bc    = equ['bcentr']     # center bt
rmaxis= equ['rmaxis']
simag = equ['simag']
sibry = equ['sibry']
rbbbs = equ['rbbbs']
zbbbs = equ['zbbbs']

nr  = bmesh['nr']
nz  = bmesh['nz']
psi = bmesh['psirz']
R   = bmesh['r']
Z   = bmesh['z']
RR  = bmesh['rr']
ZZ  = bmesh['zz']
br  = bmesh['br']
bz  = bmesh['bz']
bt  = bmesh['bt']
b   = bmesh['b']

equ_curve = umesh['equ_curve']
equ_grad  = umesh['equ_grad']

f_br = interpolate.interp2d(ZZ, RR, br, kind='cubic')
f_bz = interpolate.interp2d(ZZ, RR, bz, kind='cubic')
f_bt = interpolate.interp2d(ZZ, RR, bt, kind='cubic')
f_b  = interpolate.interp2d(ZZ, RR, b,  kind='cubic')
f_psi  = interpolate.interp2d(ZZ, RR, psi,  kind='cubic')

print('INCLUDE EVERYTING ...')

