# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 11:28:24 2023

@author: westj
"""
from robodk.robolink import *       # import the robolink library (bridge with RoboDK)
from robodk.robomath import * 
RDK = Robolink()
import pandas as pd



A=range(0,2)
B=[RDK.Item('TPP1'),RDK.Item('TPP2')]
C=range(4,6)
d={ 'A':A,'B':B,'C':C}
tdf = pd.DataFrame(d)


class Test:
    def __init__(self,rdf,index) -> None:
        self.rdf=rdf
        self.index=index
        self.Name=self.rdf['A'][self.index]
        self.IP=self.rdf['B'][self.index]
        self.Tool=self.rdf['C'][self.index]



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





class Robot:
    def __init__(self,rdf,index) -> None:
        self.rdf=rdf
        self.index=index
        self.Name=self.rdf['RID'][self.index]
        self.IP=self.rdf['IP'][self.index]
        self.Tool=self.rdf['Tool'][self.index]
        self.TPP=self.rdf['TPP'][self.index]
        self.APP=self.rdf['APP'][self.index]
        self.Ex=self.rdf['Ex'][self.index]
        self.Ey=self.rdf['Ey'][self.index]
        self.Ez=self.rdf['Ez'][self.index]
        
    def P(self,x):
        
        print(x)
        
    def hope(x):
        Robot.P(x)
        
        
r=Robot(rdf,1)
r.P(3)

print(rdf)
test=Test(tdf,0)



