"""
Read and write 'G' formatted equilibria. This is an R-Z free boundary
format.

Format of G-EQDSK file is specified here:
  https://fusion.gat.com/THEORY/efit/g_eqdsk.html


HISTORY:
  09-xx-2016 made by S.Ohdachi

  10-06-2016 modified by X.D. Du
  (1) psirz is read correctly, do not need transpose.
  (2) unify the names of variables to be the same with EFIT

"""

import numpy as np
import creatobj
from util_tokamak import file_numbers, writef

def read(f):
    g = creatobj.gfile(f)
    """ Reads a G-EQDSK file

    Parameters
    ----------

    f = Input file. Can either be a file-like object,
        or a string. If a string, then treated as a file name
        and opened.

    Returns
    -------

    """

    if isinstance(f, basestring):
        # If the input is a string, treat as file name
        with open(f) as fh: # Ensure file is closed
            return read(fh) # Call again with file object

    # Read the first line, which should contain the mesh sizes
    desc = f.readline()
    if not desc:
        raise IOError("Cannot read from input file")

    s = desc.split() # Split by whitespace
    if len(s) < 3:
        raise IOError("First line must contain at least 3 numbers")

    idum = int(s[-3])
    nw = int(s[-2])
    nh = int(s[-1])

    # Use a generator to read numbers
    token = file_numbers(f)

    rdim   = float(token.next())
    zdim   = float(token.next())
    rcentr = float(token.next())
    rleft  = float(token.next())
    zmid   = float(token.next())

    rmaxis = float(token.next())
    zmaxis = float(token.next())
    simag  = float(token.next())
    sibry  = float(token.next())
    bcentr = float(token.next())

    current= float(token.next())
    simag  = float(token.next())
    xdum   = float(token.next())
    rmaxis = float(token.next())
    xdum   = float(token.next())

    zmaxis = float(token.next())
    xdum   = float(token.next())
    sibry  = float(token.next())
    xdum   = float(token.next())
    xdum   = float(token.next())

    # Read arrays
    def read_array(n, name="Unknown"):
        data = np.zeros([n])
        try:
            for i in np.arange(n):
                data[i] = float(token.next())
        except:
            raise IOError("Failed reading array '"+name+"' of size ", n)
        return data

    # read 2d array
    def read_2d(nw, nh, name="Unknown"):
        data = np.zeros([nw, nh])
        for j in np.arange(nh):
            for i in np.arange(nw):
                data[i,j] = float(token.next())
        return data

    fpol   = read_array(nw, "fpol")
    pres   = read_array(nw, "pres")
    ffprim = read_array(nw, "ffprim")
    pprime = read_array(nw, "pprime")
    psirz  = read_2d(nw, nh, "psirz")
    qpsi   = read_array(nw, "qpsi")

    # Read boundary and limiters, if present
    nbbbs  = int(token.next())
    limitr = int(token.next())

    if nbbbs > 0:
        rbbbs = np.zeros([nbbbs])
        zbbbs = np.zeros([nbbbs])
        for i in range(nbbbs):
            rbbbs[i] = float(token.next())
            zbbbs[i] = float(token.next())
    else:
        rbbbs = [0]
        zbbbs = [0]

    if limitr > 0:
        rlim = np.zeros([limitr])
        zlim = np.zeros([limitr])
        for i in range(limitr):
            rlim[i] = float(token.next())
            zlim[i] = float(token.next())
    else:
        rlim = [0]
        zlim = [0]

    # Construct R-Z mesh
    r = np.zeros([nw, nh])
    z = r.copy()
    for i in range(nw):
        r[i,:] = rleft + rdim*i/float(nw-1)
    for j in range(nh):
        z[:,j] = (zmid-0.5*zdim) + zdim*j/float(nh-1)

    # Create dictionary of values to return
    g.nw = nw
    g.nh = nh
    g.r = r
    g.z = z
    g.rdim = rdim
    g.zdim = zdim
    g.rcentr = rcentr
    g.bcentr = bcentr
    g.rleft = rleft
    g.zmid = zmid
    g.rmaxis = rmaxis
    g.zmaxis = zmaxis
    g.simag = simag
    g.sibry = sibry
    g.current = current
    g.psirz = psirz
    g.fpol = fpol
    g.pres = pres
    g.qpsi = qpsi
    g.nbbbs = nbbbs
    g.rbbbs = rbbbs
    g.zbbbs = zbbbs
    g.limitr = limitr
    g.rlim = rlim
    g.zlim = zlim

    return g
