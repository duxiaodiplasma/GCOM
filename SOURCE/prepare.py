import numpy as np
import bgrid
import ugrid
import geqdsk
import os.path
import cPickle as pickle
import hickle as hkl

def save(fn):
    outfn1 = fn+'.orbit'+'.hkl'

    if os.path.isfile(outfn1) and os.path.isfile(outfn1):
       print('USE EXISTING FILES ...')
    else:
       equ   = geqdsk.read(fn)
       bmesh = bgrid.main(fn)
       umesh = ugrid.equ(fn)

       dics = [equ,bmesh,umesh]
       hkl.dump(dics,outfn1)

       print('** Good to Go! **')

    return outfn1#,outfn2






