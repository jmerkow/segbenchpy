import edges_swig
import numpy as np
from numpy import zeros_like
from util import normal_img, conv_tri_3D, maxeps
import os

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

###############################################################
def edge_direction_2D(E,r=5):
    if(r<=1):
        p=12/r/(r+2)-2;f=np.array([1,p,1])/(2+p); r=1;
    else:
        f = np.array(hstack((range(1,r),r+1,range(r,1,-1))),dtype=float)/(r+1)**2
    F = np.zeros((len(f),len(f)))
    F[:,len(f)/2] = f
    F = signal.convolve2d(F,F.T,mode='full')
    En = signal.convolve2d(E,F,mode='same')
    Ox,Oy = np.gradient(En)
    Oxx,_ = np.gradient(Ox);
    Oxy,Oyy = np.gradient(Oy);
    return np.arctan(Oyy*np.sign(-Oxy)/maxeps(Oxx,1e-5))

def edge_direction_3D(E,r=1):
    D = conv_tri_3D(E,r)
    Ox,Oy,Oz = np.gradient(D)
    M = np.sqrt(Ox*Ox+Oy*Oy+Oz*Oz)
    Ox,Oy,Oz = (Ox,Oy,Oz)/maxeps(M,1e-5)
    theta = np.arctan(Oy/maxeps(Oz,1e-5))
    phi = np.arctan(Oz/maxeps(M,1e-5))
    Px= np.sin(theta)*np.cos(phi)
    Py= np.sin(theta)*np.sin(phi)
    Pz= np.cos(phi)
    return Px,Py,Pz,D
def edge_auto_nms_3D(E,sigma1=1,sigma2=5,r=1,m=1.001,border_s=0,
    debug_save=False, size_check=True):
    E = normal_img(E)
    if sigma1>0:
        Eb = conv_tri_3D(E,sigma1)
    else:
        Eb = E
    vmin,vmax = E.min(),E.max()
    dX, dY, dZ, _ = edge_direction_3D(Eb,r=sigma2)
    if size_check: print(E.shape)
    if debug_save:
        savearr(dX,"dX.txt")
        savearr(dY,"dY.txt")
        savearr(dZ,"dZ.txt")
        savearr(E,"E.txt")
    Enms = edgeNms3d(Eb, dY,dX, dZ, r=r, s=border_s, m=m)
    return normal_img(Enms)

def savearr(arr,fn):
    shfn = os.path.splitext(fn)[0]+'-dim.txt'
    dims = arr.shape
    arr = arr.ravel()
    np.savetxt(shfn,dims,fmt='%d')
    np.savetxt(fn,arr,fmt='%.10f')