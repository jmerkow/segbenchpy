from __future__ import print_function
import benchmark_swig
import numpy as np
from numpy import zeros, zeros_like
import sys, time, os
from util import normal_img

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
###############################################################
def eval_edgemap(E,GT,id=None,res_dir=None,K=99):
    thrs= np.linspace(1/K,1,K)
    ret = list()
    if id is None:
        id = 'results'
    if res_dir is not None:
        fname = os.path.join(res_dir,str(id)+'.txt')
    else:
        fname = None
    print("[",id,"]",end='',sep='')
    for t,th in enumerate(thrs):
            Et = (normalimg(E)>th).astype(float)
            _,_,matchE,matchG = correspondVoxels((Et>0).astype(float),
                                                   (GT>0).astype(float),.0075, 8)
            matchE,matchG = matchE>0,matchG>0
            d = dict(
                th=th,                           #th
                cntR=np.sum(matchG.astype(int)), #cntR
                sumR=np.sum(GT.astype(int)),     #sumR
                cntP=np.sum(matchE.astype(int)), #cntP
                sumP=np.sum(Et.astype(int))      #sumP
                )
            print(".",end='')        
            ret.append(d)
    print("!",end='')
    if fname is not None:
        print(fname)
        keys = ['th','cntR','sumR','cntP','sumP']
        ret2 = [ [r[k] for k in keys] for r in ret]
        header = "\t\t{:10s}{:10s}{:10s}{:10s}{:10s}".format(*keys)
        np.savetxt(fname,ret2,fmt="%10g %10g %10g %10g %10g",
                   header=header)
    return ret