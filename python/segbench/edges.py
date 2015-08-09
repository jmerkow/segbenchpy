import edges_swig
from numpy import zeros_like


def edgeNms2d(E0, O, r=1, s=5, m=1.01):
    E = zeros_like(E0).astype(float)
    edges_swig._edgeNms2d(E0, O, E, r, s, m)
    return E