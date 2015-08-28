from __future__ import print_function
from scipy.ndimage.filters import convolve1d as convn
import numpy as np

def maxeps(arr,eps=np.finfo(float).eps):
    if hasattr(arr, '__iter__'):
        arr=np.array(arr)
        sh=arr.shape
        return np.array([max(eps,a) for a in arr.ravel()]).reshape(sh)                 
    else:
        return max(eps,arr)

def normal_img(E):
    if (E.max()-E.min())>0:
        return (E-E.min())/(E.max()-E.min()+1e-10)
    else:
        return E

def triangle_filt(r):
    if(r<=1):
        p=12.0/r/(r+2)-2;f=np.array([1,p,1])/(2+p); r=1;
    else:
        f = np.array(np.hstack((np.linspace(1,r,r),r+1,np.linspace(r,1,r))),
                     dtype=float)/(r-1)**2
    return f
def conv_tri_3D(img,r=5):
    f = triangle_filt(r)
    return convn(convn(convn(img,f,axis=0),f,axis=1),f,axis=2)