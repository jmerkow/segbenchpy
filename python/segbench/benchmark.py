from __future__ import print_function
import benchmark_swig
import numpy as np
from numpy import zeros, zeros_like
import sys, time, os
from util import normal_img, maxeps
import fnmatch, re
import matplotlib.pyplot as plt
from tabulate import tabulate
npa=np.array
pjoin = os.path.join

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
    thrs= np.linspace(1.0/K,1,K)
    ret = list()
    if id is None:
        id = 'results'
    if res_dir is not None:
        fname = os.path.join(res_dir,str(id)+'.txt')
    else:
        fname = None
    print("[",id,"]",end='',sep='')
    num_dots = 30
    dot_on = np.ceil((1.0+K)/num_dots)
    for t,th in enumerate(thrs):
            Et = (normal_img(E)>th).astype(float)
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
            if t%dot_on==0: print(".",end='')        
            ret.append(d)
    print("!",end='')
    if fname is not None:
        keys = ['th','cntR','sumR','cntP','sumP']
        ret2 = [ [r[k] for k in keys] for r in ret]
        header = "\t\t{:10s}{:10s}{:10s}{:10s}{:10s}".format(*keys)
        np.savetxt(fname,ret2,fmt="%10g %10g %10g %10g %10g",
                   header=header)
        print("!")
    return ret

def computeRPF(cntR,sumR,cntP,sumP):
    R = np.divide(cntR.astype(float),sumR+1e-5)
    P = np.divide(cntP.astype(float),sumP+1e-5)
    d = npa(maxeps(P+R))
    F =(2*R*P)/d
    return R,P,F

def findBestRPF(T,R,P):
    A = np.linspace(0,1,100)
    B = np.add(1,-A)
    l = len(R)
    def interpt(G,A,B):
        return np.array([G[j]*A + G[j-1]*B for j in range(1,len(G))])
    Tj,Rj,Pj =  interpt(T,A,B), interpt(R,A,B), interpt(P,A,B)
    Fj = (2*Rj*Pj)/(Rj+Pj+1e-5)
    k = np.unravel_index(np.argmax(Fj),Fj.shape)
    bstR,bstP,bstF,bstT = Rj[k],Pj[k],Fj[k],Tj[k]
    return bstR,bstP,bstF,bstT
def collect_result_files(res_dir):
    allresfiles = sorted(os.listdir(res_dir))
    resfiles = dict()
    reg = '([\d]*)\.txt'
    reobj = re.compile(reg)
    for f in allresfiles:
        m = reobj.match(f)
        if m:
            resfiles[int(m.groups()[0])] = f
    resfiles=list(sorted(resfiles.itervalues()))
    return resfiles,len(resfiles)
def evaluate_edge_dir(res_dir):
    resfiles,num_res_files = collect_result_files(res_dir)
    res = npa([np.loadtxt(pjoin(res_dir,r)) for r in resfiles]).transpose(2,1,0)
    Ts = npa([res[0,:,r] for r in range(num_res_files)]).T
    allCnts = npa([res[1:,:,r] for r in range(num_res_files)]).transpose(1,2,0)

    Rs,Ps,Fs = npa([computeRPF(*allCnts[...,r]) for r in range(num_res_files)]).transpose(1,2,0)
    Ks = np.argmax(Fs,axis=0)
    bstRs,bstPs,bstFs,bstTs = npa([findBestRPF(Ts[:,i],Rs[:,i],Ps[:,i]) 
                                       for i in range(num_res_files)]).T

    oisCnts = np.sum([allCnts[:,Ks[i],i] for i in range(num_res_files)],axis=0)
    oisR,oisP,oisF  = computeRPF(*(oisCnts))
    T = Ts[:,0]
    cntSum = np.sum(allCnts,axis=2)
    R,P,F  = computeRPF(*(cntSum))
    odsR,odsP,odsF,odsT = findBestRPF(T,R,P)
    AP=np.interp(np.arange(0,1,.01),P[P>0],R[P>0],left=0,right=0); AP = np.sum(AP)/100

    num_out = dict(odsR=odsR,odsP=odsP,odsF=odsF,odsT=odsT,
                  oisR=oisR,oisP=oisP,oisF=oisF,AP=AP)
    
    graph_out = dict(R=R,P=P,F=F,T=T,
                     Rs=Rs,Ps=Ps,Fs=Fs,Ts=Ts,
                     bstRs=bstRs, bstPs=bstPs, bstTs=bstTs, bstFs=bstFs)
    return num_out,graph_out

def plot_results(res_dirs,row_labels,headers = ['AP','odsF','oisF','nimg'],
                show_per_image=False,fig_size=(5,10),title=None):

        fig,(axRP,axF) = plt.subplots(1,2,sharey='col')
        fig.set_figwidth(fig_size[1]);fig.set_figheight(fig_size[0])
        cmap = plt.get_cmap('gnuplot2')
        layer_colors = [cmap(i) for i in np.linspace(.1, .8, len(res_dirs.keys()))]
        num_outc=dict()
        tab = list()
        for c,layer in zip(layer_colors,sorted(res_dirs.keys())):
            try:
                num_outc[layer],graph_out=evaluate_edge_dir(res_dirs[layer])
                axRP.plot(graph_out['R'][:-1],graph_out['P'][:-1],'-',color=c,label=row_labels[layer])
                axRP.plot(num_outc[layer]['odsR'],num_outc[layer]['odsP'],'o',color=c)
                axF.plot(graph_out['T'],graph_out['F'],'-',color=c)
                axF.plot(num_outc[layer]['odsT'],num_outc[layer]['odsF'],'o',color=c)
                num_outc[layer]['nimg'] = len(graph_out['bstTs'])

            except:
                num_outc[layer] = dict()
                for h in headers:
                    num_outc[layer][h]= "N/A"
            row = [row_labels[layer]]+[num_outc[layer][h] for h in headers]
            tab.append(row)
        axRP.legend(ncol=len(layer_colors),loc='upper center', bbox_to_anchor=(1.1, -0.15))
        axRP.set_xlim(0,1);axRP.set_ylim(0,1)
        axRP.set_ylabel('Recall');axRP.set_xlabel('Precision')
        axF.set_xlim(0,1);axF.set_ylim(0,1)
        axF.set_ylabel('F-Score');axF.set_xlabel('Threshold')
        if title is not None:
            plt.suptitle(title)
        plt.show()                                     
        print(tabulate(tab,headers=headers,tablefmt="fancy_grid",
                       floatfmt="0.4f",numalign="decimal",stralign='right'))
        if show_per_image:
            fig2,(axRPs,axFs) = plt.subplots(1,2)
            fig2.set_figwidth(10);fig2.set_figheight(5)
            cmap = plt.get_cmap('jet')
            colors = [cmap(i) for i in np.linspace(.1, .8, len(graph_out['bstTs']))]
            for i,c in enumerate(colors):
                axRPs.plot(graph_out['Rs'][:-1,i],graph_out['Ps'][:-1,i],'-',color=c,label=str(i))
                axRPs.plot(graph_out['bstRs'][i],graph_out['bstPs'][i],'o',color=c)
                axFs.plot(graph_out['Ts'][:,i],graph_out['Fs'][:,i],'-',color=c)
                axFs.plot(graph_out['bstTs'][i],graph_out['bstFs'][i],'o',color=c)
            axRPs.legend(ncol=5,loc='upper center', bbox_to_anchor=(1, -0.05))
        print("")