import edges_swig
from numpy import zeros_like


def edgeNms2d(E0, O, r=1, s=5, m=1.01):
    E = zeros_like(E0).astype(float)
    edges_swig._edgeNms2d(E0, O, E, r, s, m)
    return E

def edgeNms3d(E0, dX, dY, dZ, r=1, s=5, m=1.01):
    E = zeros_like(E0).astype(float)
    edges_swig._edgeNms3d(E0, dX, dY, dZ, E, r, s, m)
    return E

def interp3(*args):
    return edges_swig._interp3(*args)