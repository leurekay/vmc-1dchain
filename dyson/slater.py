# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 17:59:42 2016
Ratio2是错误的，会导致负符号
Ratio:real-space skill
Ratio1:real-space brute
Ratio3:k-space brute
Ratio4:k-space Dyson
@author: aa
"""

from __future__ import division
import numpy as np
from hamiltonian import H
import scipy
from math import *
from configuration import conf_initial

def Ratio(conf,before,after,L,u):  #before-th site to the after-th site
    I=np.identity(L)
    #h,dia,u=H(L)
    M=u[:,np.arange(0,L,1)]
    occupy_index=[]
    for i in range(np.size(conf)):
        if i==before:
            beta=len(occupy_index)    #beta :correspoding position in D matrix or the particle index
        if conf[0][i]!=0:
            occupy_index.append(i)
    D=M[occupy_index]+0.00000000001*I
    #D=M[occupy_index]
    #print "##########",np.linalg.det(D)
    W=np.dot(M,np.linalg.inv(D))
    #print W[after][beta]
    
    W=(scipy.linalg.solve(D.T,M.T)).T
    #print W[after][beta]
    return [W[after][beta],occupy_index]

def Ratio2(conf,before,after,L,u):  #before-th site to the after-th site
    #h,dia,u=H(L)
    M=u[:,np.arange(0,L,1)]
    occupy_index=[]
    for i in range(np.size(conf)):       
        if conf[0][i]!=0:
            occupy_index.append(i)
    D=M[occupy_index]
    
    conf[0][after]=conf[0][before]
    conf[0][before]=0
    occupy_index=[]
    for i in range(np.size(conf)):
        if conf[0][i]!=0:
            occupy_index.append(i)
    D1=M[occupy_index]
    ratio=np.linalg.det(D1)/np.linalg.det(D)
    #print ratio
    return [ratio,0,0,0,0]

def Ratio1(conf,before,after,L,u):
    M=u[:,np.arange(0,L,1)]
    occupy_index=[]
    for i in range(np.size(conf)):
        if i==before:
            beta=len(occupy_index)    #beta :correspoding position in D matrix or the particle index
        if conf[0][i]!=0:
            occupy_index.append(i)
    D=M[occupy_index]
    det_old=np.linalg.det(D)
    D[beta,:]=M[after,:]
    det_new=np.linalg.det(D)
    #print det_old,det_new
    return det_new/det_old

def Det(conf,L,u):
    #h,dia,u=H(L)
    M=u[:,np.arange(0,L,1)]
    occupy_index=[]
    for i in range(np.size(conf)):       
        if conf[0][i]!=0:
            occupy_index.append(i)
    D=M[occupy_index]
    return np.linalg.det(D)
    
def Ratio3(conf,before,after,L,N_up):
    N_down=L-N_up
    k_range=np.linspace(-pi,pi,L+1)
    
    dispersion=-np.cos(k_range)
    idx=dispersion.argsort()
    k_up=idx[0:N_up]
    k_down=idx[0:N_down]
    k_up.sort()
    k_down.sort()
    
    occupy_index=[]
    for i in range(np.size(conf)):
        if i==before:
            beta=len(occupy_index)    #beta :correspoding position in D matrix or the particle index
        if conf[0][i]!=0:
            occupy_index.append(i)    
    slater_up=np.zeros((N_up,N_up),"complex")
    slater_down=np.zeros((N_down,N_down),"complex")
    if before<L:
        #indict hopping in spin up zone
        for x in range(N_up):
            for k in range(N_up):
                slater_up[k][x]=scipy.exp(1j*k_range[k_up[k]]*occupy_index[x])
        det_old=np.linalg.det(slater_up)
        for k in range(N_up):
            slater_up[k][beta]=scipy.exp(1j*k_range[k_up[k]]*after)        
        det_new=np.linalg.det(slater_up)      
        return det_new/det_old
    else:
        after=after-L
        beta=beta-N_up
        for x in range(N_down):
            for k in range(N_down):
                slater_down[k][x]=scipy.exp(1j*k_range[k_down[k]]*(occupy_index[x+N_up]-L))
        det_old=np.linalg.det(slater_down)
        for k in range(N_down):
            slater_down[k][beta]=scipy.exp(1j*k_range[k_down[k]]*after)        
        det_new=np.linalg.det(slater_down)
        return det_new/det_old     


def Ratio4(conf,before,after,L,N_up):
    N_down=L-N_up
    k_range=np.linspace(-pi,pi,L+1)
    dispersion=-np.cos(k_range)
    idx=dispersion.argsort()
    k_up=idx[0:N_up]
    k_down=idx[0:N_down]
    k_up.sort()
    k_down.sort()
    occupy_index=[]
    for i in range(np.size(conf)):
        if i==before:
            beta=len(occupy_index)    #beta :correspoding position in D matrix or the particle index
        if conf[0][i]!=0:
            occupy_index.append(i)    
    slater_up=np.zeros((N_up,N_up),"complex")
    slater_down=np.zeros((N_down,N_down),"complex")
    
    if before<L:
        #indict hopping in spin up zone
        for x in range(N_up):
            for k in range(N_up):
                slater_up[k][x]=scipy.exp(1j*k_range[k_up[k]]*occupy_index[x])
        Sigma=np.zeros((N_up,),"complex")
        for k in range(N_up):
            Sigma[k]=slater_up[k][beta]-scipy.exp(1j*k_range[k_up[k]]*after)
        X=scipy.linalg.solve(slater_up,Sigma)
        return 1-X[beta]
    else:
        after=after-L
        beta=beta-N_up
        for x in range(N_down):
            for k in range(N_down):
                slater_down[k][x]=scipy.exp(1j*k_range[k_down[k]]*(occupy_index[x+N_up]-L))
        Sigma=np.zeros((N_down,),"complex")
        for k in range(N_down):
            Sigma[k]=slater_down[k][beta]-scipy.exp(1j*k_range[k_up[k]]*after)
        X=scipy.linalg.solve(slater_down,Sigma)
        return 1-X[beta]


    
if __name__=='__main__' :
    L=30
    N_up=15
    h,dia,u=H(L)
    conf=conf_initial(L,N_up)
    
    #conf=np.array([[0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0]])
    print Ratio(conf,22,23,L,u)[0]
    print Ratio1(conf,22,23,L,u)
    #print Ratio2(conf,15,8,L,u)
    print Ratio3(conf,22,23,L,N_up)
    print Ratio4(conf,22,23,L,N_up)