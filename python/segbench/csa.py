import csa_swig
from numpy import array

def csaAssign(n, g):
    return array (csa_swig._csaAssign(n, g))