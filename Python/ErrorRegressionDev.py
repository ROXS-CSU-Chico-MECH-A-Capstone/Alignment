# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:56:54 2023

@author: westj
"""
from robodk.robolink import *       # import the robolink library (bridge with RoboDK)
from robodk.robomath import * 
RDK = Robolink()   
import math
import numpy as np
import matplotlib.pyplot as plt 
import array
import pandas as pd
import requests
import re
from scipy.signal import find_peaks, peak_prominences
from scipy.stats import linregress
from scipy.ndimage import uniform_filter1d
import json
import time
from art import *



Ex=Ey=Ez=0

def BuildRowScan(zi,zf,points):
    PRZ=np.linspace(zi,zf,points)
    fl=1000#focal length in mm
    Rlist=[]
    X=[]
    Y=[]
    Zl=[]
    D=[]
    for i in range(0,len(PRZ)):
        Z=-1*PRZ[i] #call coord for particular scan
        theta=np.arccos(Z/(2*fl))
        R=Z/2*np.tan(theta) #distance from emitter detector axis
        coord=(-Z/2+Ex,0+Ey,-R+Ez)#rowland circle scan position offset from LEDT
        
        d=np.sqrt(coord[0]**2+coord[2]**2)
        X.append(coord[0])
        Y.append(coord[1])
        Zl.append(coord[2])
        D.append(d)
        Rlist.append(coord)
        
    Plot3DS(X,Y,Zl)
    print(D)
    return Rlist

def Plot3DS(X,Y,I,prominence,width,height):

    ps=find_peaks(I,prominence=[prominence],width=width,height=height)[0][0]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(X, Y, I)
    ax.plot3D(X[ps], Y[ps], I[ps],'x')
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show() 
    
def Plot2DS(Y,I,title):
    #fig1 = plt.figure()
    #plt.style.use('seaborn-notebook')
    #theta=np.array(np.linspace(0,revs*2*np.pi,points))
    peaks, _ = find_peaks(I)

    # height=.1,distance=3,threshold=0.01,
    prominences = peak_prominences(I, peaks)[0]
    contour_heights = np.array(I)[peaks] - prominences
    plt.figure(np.random.rand())
    plt.plot(I)
    plt.plot(peaks, np.array(I)[peaks], "x")
    plt.vlines(x=peaks, ymin=contour_heights, ymax=np.array(I)[peaks])
    
    plt.title(title)
    plt.xlabel('Translation')
    plt.ylabel('Intensity')
    plt.show()
    
    
    
    
    
    LEDT= RDK.Item('LEDT',ITEM_TYPE_TARGET)
    LEDT.setPose(Pose(0,0,0,0,90,0))
    tx,ty,tz,ta,tb,tc=Pose_2_KUKA(TPP.Pose())
    LEDT.setPose(LEDT.Pose()*rotx(np.deg2rad(-ta)))
    
    
    #%%
# E=  [[array([100., 200., 300.])],[(-0.4013895581238693, 0.22907231479107754, 7.25),(0.029686906117525898, -1.253417329681019, 5.25),(-0.3073011508995377, -1.0086050071245312, 0.75)], []]
  
  
  
#   LRR=linregress(X,E) #run regression
  
  
#   for i in range(1,len(E[0][0])):
      
#       xl=E[0][i-1]
#       xh=E[0][i]
  
#   A[1][1]
  
# if zl <= number <= zh:
    #%%
import numpy as np
A = 2*np.arange(10)

idx = (A>zl)*(A<zh)
w=np.where(idx)
print(w)

#%%
#Error=[[],[],[]]
E=[[100.0, 200.0, 300.0],
 [(-0.4013895581238693, 0.22907231479107754, 7.25),
  (0.029686906117525898, -1.253417329681019, 5.25),
  (-0.3073011508995377, -1.0086050071245312, 0.75)],
 []]

def ErrorReg(E,R_index):
    ERL=[E[0],[]]
    
    for i in range (1,len(E[0])):
        LRL=[]
        print('i is'+str(i))
        print(E[0])
        zlow=E[0][i-1]
        zhigh=E[0][i] #gantry position
        X=[zlow,zhigh]
        clow=E[R_index][i-1]
        chigh=E[R_index][i]
        
        for j in range(0,3):
            print(j)
            print(clow,chigh)
            el=clow[j] #coordinate error
            eh=chigh[j]
            Er=[el,eh]
            print(X,E)
            LR=linregress(X,Er)
            LRL.append(LR)
        ERL[1].append(LRL)
            
    return ERL

            
    def BuildRowScan(self,zi,zf,points):
        PRZ=np.linspace(zi,zf,points)
        fl=1000#focal length in mm
        Rlist=[]
        X=[]
        Y=[]
        Zl=[]
        D=[]
        for i in range(0,len(PRZ)):
            Z=-1*PRZ[i] #call coord for particular scan
            theta=np.arccos(Z/(2*fl))
            R=Z/2*np.tan(theta) #distance from emitter detector axis
            coord=(-Z/2+self.Ex,0+self.Ey,-R+self.Ez)#rowland circle scan position offset from LEDT
            
            d=np.sqrt(coord[0]**2+coord[2]**2)
            X.append(coord[0])
            Y.append(coord[1])
            Zl.append(coord[2])
            D.append(d)
            Rlist.append(coord)
            
        self.Plot3DS(X,Y,Zl)
        print(D)
        return Rlist
    
    
def Plot3DS(X,Y,I):

    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(X, Y, I)
    
    
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    
    # ax.set_xlim3d(-100, -500)
    # ax.set_ylim3d(-100, -500)
    # ax.set_zlim3d(-100, -500)
    
    
    ax.set_aspect('equal')
    plt.show() 
    
def BuildRowScanADV(ERL,zi,zf,points):
        PRZ=np.linspace(zi,zf,points)
        fl=1000#focal length in mm
        Rlist=[]
        X=[]
        Y=[]
        Zl=[]
        D=[]
        
        for i in range(0,len(PRZ)):
            Z=-1*PRZ[i] #call coord for particular scan
            print(Z)
            
            for j in range(0,len(ERL)):
                print(j)
                if -ERL[0][j] >= Z >= -ERL[0][j+1]:
                    print('wahoo')
                    
                    
                    mx=ERL[1][j][0][0] #1 regressions j index for which position 0 x regression 0 slope
                    bx=ERL[1][j][0][1] #1 regressions j index for which position 0 x regression 1 b offset
                    
                    my=ERL[1][j][1][0] #1 regressions j index for which position 0 y regression 0 slope
                    by=ERL[1][j][1][1] #1 regressions j index for which position 0 y regression 1 b offset

                    mz=ERL[1][j][2][0] #1 regressions j index for which position 0 z regression 0 slope
                    bz=ERL[1][j][2][1] #1 regressions j index for which position 0 z regression 1 b offset

                    
                    
                    Ex=mx*Z+bx
                    Ey=my*Z+by
                    Ez=mz*Z+bz
                    print('Errors'+str((Ex,Ey,Ez)))
                    break
                    
                    

            theta=np.arccos(Z/(2*fl))
            R=Z/2*np.tan(theta) #distance from emitter detector axis
            
            coord=(-Z/2+Ex,0+Ey,-R+Ez)#rowland circle scan position offset from LEDT
            
            d=np.sqrt(coord[0]**2+coord[2]**2)
            X.append(coord[0])
            Y.append(coord[1])
            Zl.append(coord[2])
            D.append(d)
            Rlist.append(coord)
            
        Plot3DS(X,Y,Zl)
        print(D)
        return Rlist       
            
             
ERL=ErrorReg(E,1)
print(ERL[1])
BuildRowScanADV(ERL,100,300,200)


    
    
    
    