# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 12:41:48 2023

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
from scipy.signal import find_peaks, peak_prominences, savgol_filter
from scipy.stats import linregress
import json
import time
import threading




class Robot:
    def __init__(self,rdf,index) -> None:
        self.rdf=rdf
        self.index=index
        # self.Name=self.rdf['RID'][self.index]
        # self.IP=self.rdf['IP'][self.index]
        # self.Tool=self.rdf['Tool'][self.index]
        # self.TPP=self.rdf['TPP'][self.index]
        # self.APP=self.rdf['APP'][self.index]
        # self.Ex=self.rdf['Ex'][self.index]
        # self.Ey=self.rdf['Ey'][self.index]
        # self.Ez=self.rdf['Ez'][self.index]
        
        
    def RScan(self,lock,Rlist,t):
        global state
        i=0  
        while i< len(Rlist):
            time.sleep(0.1)
            #with lock:
            if state[t]==0:
                #safeMJ(Rlist[i])
                print('thread'+str(t),i)
                state[t]=1
                i+=1
                time.sleep(0.1)
                if state==[1,1,1]:
                    state=[0,0,0]
                
            elif state[t]==1:
                time.sleep(0.1)
#%%
IPs = ['192.168.0.100','192.168.0.101']
PORT = [30002,30003]
R_NAMES = ['1 Mecademic Meca500 R3','2 Mecademic Meca500 R3']
RTOOLS = ['1 CrystalToolSample','2 CrystalToolSample']

obj=[]

for nrobot in range(0,len(R_NAMES)):
    obj.append("robot"+str(nrobot+1))
    
robots={ 'OBJ':obj,'Names':R_NAMES,'Tools':RTOOLS,'IP': IPs,'Port': PORT}
rdf = pd.DataFrame(robots)

RID=[]
TPP=[]
APP=[]
Tool=[]
for nrobot in range(0,len(R_NAMES)):
    RID.append(RDK.Item(R_NAMES[nrobot]))
    Tool.append(RDK.Item(str(nrobot+1)+' CrystalToolSample'))
    TPP.append(RDK.Item('TPP'+str(nrobot+1)))
    APP.append(RDK.Item('APP'+str(nrobot+1)))

    
rdf["RID"]=RID
rdf["Tool"]=Tool
rdf["TPP"]=TPP
rdf["APP"]=APP

Ex=[]
Ey=[]
Ez=[]

for n in range(0,len(R_NAMES)):
    APP=rdf['APP'][n]
    ex,ey,ez,ea,eb,ec=Pose_2_KUKA(APP.Pose())
    Ex.append(ex)
    Ey.append(ey)
    Ez.append(ez)
    
rdf['Ex']=Ex
rdf['Ey']=Ey
rdf['Ez']=Ez







#%%

lock = threading.Lock()
state=[0,0,0]

R=Re(1,1)
l=np.linspace(1,10,10)
# Create two threads
tg = threading.Thread(target=R.RScan,args=(lock,l,0))
t1 = threading.Thread(target=R.RScan,args=(lock,l,1))
t2 = threading.Thread(target=R.RScan,args=(lock,l,2))

# Start the threads
tg.start()
t1.start()
t2.start()
print('Threads started')

# Wait for the threads to finish
t1.join()
t2.join()
tg.join()

print("Done!")