import numpy as np
class gfile(object):
    """
    Attributes:
    """
    def __init__(self,fn):
        self.fn = fn

class inpu(object):
    """
    Attributes:
    """
    def __init__(self,R0,Z0,phi0,pitch0,E0):
        self.R0 = R0
        self.Z0 = Z0
        self.phi0 = phi0
        self.pitch0 = pitch0
        self.E0 = E0

class outpu(object):
    """
    Attributes:
    """
    def __init__(self,comment):
        self.comment = comment


