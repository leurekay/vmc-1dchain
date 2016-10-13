# -*- coding: utf-8 -*-
"""
Created on Sun Aug 28 17:59:42 2016
Ratio2是错误的，会导致负符号
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
    
def Ratio3(conf,before,after,L,u):
    k_range=np.linspace(-pi,pi,L)
    #k_range=a[2:6]
    #k_range=np.append(k_range,k_range)
    slater=np.zeros((L,L),"complex")
    occupy_index=[]
    for i in range(np.size(conf)):
        if i==before:
            beta=len(occupy_index)    #beta :correspoding position in D matrix or the particle index
        if conf[0][i]!=0:
            occupy_index.append(i)
    for x in range(L):
        for k in range(L):
            slater[k][x]=scipy.exp(1j*k_range[k]*occupy_index[x])
    det_old=np.linalg.det(slater)
    for k in range(L):
        slater[k][beta]=scipy.exp(1j*k_range[k]*after)
    det_new=np.linalg.det(slater)
    return det_new/det_old

    
if __name__=='__main__' :
    L=8
    h,dia,u=H(L)
    conf=conf_initial(L,4)
    
    conf=np.array([[0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0]])
    print Ratio(conf,2,1,L,u)
    print Ratio1(conf,2,1,L,u)
    #print Ratio2(conf,15,8,L,u)
    print Ratio3(conf,2,1,L,u)