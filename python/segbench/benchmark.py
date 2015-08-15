import benchmark_swig
from numpy import zeros, zeros_like

def correspondPixels(bmap1, bmap2, maxDist=0.005):
    match1 = zeros_like(bmap1)
    match2 = zeros_like(bmap2)
    cost,oc = benchmark_swig._correspondPixels(bmap1, bmap2, match1, match2, maxDist)
    return cost,oc,match1,match2
def correspondVoxels(bmap1, bmap2, maxDist=0.05,degree=6):
    match1 = zeros_like(bmap1)
    match2 = zeros_like(bmap2)
    cost,oc = benchmark_swig._correspondVoxels(bmap1, bmap2, match1, match2, maxDist, degree)
    return cost,oc,match1,match2

def matchEdgeMaps2D(bmap1, bmap2, maxDist=0.005, outlierCost=100, degree=6):
    m1 = zeros_like(bmap1)
    m2 = zeros_like(bmap2)
    cost = benchmark_swig._matchEdgeMaps2D(bmap1, bmap2, m1, m2, 
        maxDist, outlierCost, degree)
    return cost,m1,m2

def matchEdgeMaps3D(bmap1, bmap2, maxDist=.005, outlierCost=100, degree=6):
    m1 = zeros_like(bmap1)
    m2 = zeros_like(bmap2)
    cost = benchmark_swig._matchEdgeMaps3D(bmap1, bmap2, m1, m2, 
        maxDist, outlierCost, degree)
    return cost,m1,m2