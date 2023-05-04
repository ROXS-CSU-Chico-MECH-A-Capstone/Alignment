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
from scipy.signal import find_peaks, peak_prominences
from scipy.stats import linregress
from scipy.ndimage import uniform_filter1d
import json
import time

class Robot:
    def __init__(self,Name,IP,Tool,TPP,APP,LEDT,Home,ESP):
        self.Name=RDK.Item(Name)
        self.IP=IP
        self.Tool=RDK.Item(Tool)
        self.TPP=RDK.Item(TPP)
        self.APP=RDK.Item(APP)
        self.LEDT=RDK.Item(LEDT)
        self.Home=RDK.Item(Home)
        Error=Pose_2_KUKA(self.APP.Pose())
        self.Ex=Error[0]
        self.Ey=Error[1]
        self.Ez=Error[2]
        self.ESP=ESP
        self.JPos=self.Name.Joints()
        
    def SafeLM(self,target,trans_coord):
        Final_Pose=target.Pose()*transl(trans_coord)
        
        if self.Name.MoveL_Test(self.JPos,Final_Pose) == 0:         # linear move to the approach position
            
            self.Name.MoveL(Final_Pose)
            self.JPos=self.Name.Joints()
        else :
            print("collision avoided at "+ str(trans_coord))
            
    def SafeJM(self,target,trans_coord):
        Final_Pose=target.Pose()*transl(trans_coord)
        self.Name.setSpeed(-1,30)
        if self.Name.MoveJ_Test(self.JPos,Final_Pose) == 0:         # linear move to the approach position
        
            self.Name.MoveJ(Final_Pose)
            self.JPos=self.Name.Joints()
        else :
            print("collision avoided")
            
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
    
    def ErrorReg(self,E,R_index):
        ERL=[E[0],[]]
        
        for i in range (1,len(E[0])):
            LRL=[]
            #print('i is'+str(i))
            #print(E[0])
            zlow=-E[0][i-1]
            zhigh=-E[0][i] #gantry position
            X=[zlow,zhigh]
            clow=E[R_index][i-1]
            chigh=E[R_index][i]
            
            for j in range(0,3):
                #print(j)
                #print(clow,chigh)
                el=clow[j] #coordinate error
                eh=chigh[j]
                Er=[el,eh]
                #print(X,E)
                LR=linregress(X,Er)
                LRL.append(LR)
            ERL[1].append(LRL)
        return ERL
    
    def BuildRowScanADV(self,ERL,zi,zf,points):
            PRZ=np.linspace(zi,zf,points)
            fl=1000#focal length in mm
            Rlist=[]
            X=[]
            Y=[]
            Zl=[]
            D=[]
            
            for i in range(0,len(PRZ)):
                Z=-1*PRZ[i] #call coord for particular scan
                #print(Z)
                
                for j in range(0,len(ERL)):
                    #print(j)
                    if -ERL[0][j] >= Z >= -ERL[0][j+1]:
                        print('used eq set'+str(j))
                        
                        
                        mx=ERL[1][j][0][0] #1 regressions j index for which position 0 x regression 0 slope
                        bx=ERL[1][j][0][1] #1 regressions j index for which position 0 x regression 1 b offset
                        
                        my=ERL[1][j][1][0] #1 regressions j index for which position 0 y regression 0 slope
                        by=ERL[1][j][1][1] #1 regressions j index for which position 0 y regression 1 b offset
    
                        mz=ERL[1][j][2][0] #1 regressions j index for which position 0 z regression 0 slope
                        bz=ERL[1][j][2][1] #1 regressions j index for which position 0 z regression 1 b offset
    
                        
                        
                        Ex=mx*Z+bx
                        Ey=my*Z+by
                        Ez=mz*Z+bz
                        #print('Errors'+str((Ex,Ey,Ez)))
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
                
            self.Plot3DS(X,Y,Zl)
            print(D)
            return Rlist       
    
    def Plot3DS(self,X,Y,I):
        ps=I.index(max(I))
        #ps=find_peaks(I,prominence=[prominence],width=width,height=height)[0][0]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D(X, Y, I)
        ax.plot3D(X[ps], Y[ps], I[ps],'x')
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        plt.show() 
        
    def Plot2DS(self,Y,I,title):
        #fig1 = plt.figure()
        #plt.style.use('seaborn-notebook')
        #theta=np.array(np.linspace(0,revs*2*np.pi,points))
        peaks, _ = find_peaks(I)

        # height=.1,distance=3,threshold=0.01,
        prominences = peak_prominences(I, peaks)[0]
        contour_heights = np.array(I)[peaks] - prominences
        plt.plot(I)
        plt.plot(peaks, np.array(I)[peaks], "x")
        plt.vlines(x=peaks, ymin=contour_heights, ymax=np.array(I)[peaks])
        
        plt.title(title)
        plt.xlabel('Translation')
        plt.ylabel('Intensity')
        plt.show()
    
    
    def Connect_and_Home(self):
        print(self.Name)
        
        self.Name.ConnectSafe(self.IP)
        
        
        self.Name.setPoseTool(self.Tool) # Update the TCP
        self.Name.setSpeed(-1,25)  # Set linear speed in mm/s
        self.Name.setSpeed(25)  # Set linear speed in mm/s
        #print(self.Name.Joints())
        #print(self.Home.Pose())
        self.SafeJM(self.Home,(0,0,0))
        self.Name.setJoints([0,0,0,0,0,0])
        self.Name.WaitFinished()
        
    def SetTargets(self,PR_pos):
        # Get a reference to the target object
        PR = RDK.Item('Photo Resistor',ITEM_TYPE_TARGET)
        CT= RDK.Item('Crystal Target',ITEM_TYPE_TARGET)
        TPP1= RDK.Item('TPP1',ITEM_TYPE_TARGET)
        TPP2= RDK.Item('TPP2',ITEM_TYPE_TARGET)
        RB1 = RDK.Item('1 Mecademic Meca500 R3 Base',ITEM_TYPE_FRAME)
        RB2 = RDK.Item('2 Mecademic Meca500 R3 Base',ITEM_TYPE_FRAME)
        LED= RDK.Item('LED',ITEM_TYPE_FRAME)
        BI=RDK.Item('Bing',ITEM_TYPE_TARGET)
        LEDT1= RDK.Item('LEDT1',ITEM_TYPE_TARGET)
        LEDT2=RDK.Item('LEDT2')
        
        
        
        
        
        
        TPP1.setParent(LED)
        TPP2.setParent(LED)
        
        # Define the position of the photoresistor
        zPR=-PR_pos
        #set pose of Hardware
        PR.setPose(Pose(0,0,zPR,0,0,0)) #photo Resistor
        CT.setPose(Pose(0,0,zPR/2,0,0,0)) #crystal Target point
        
        
        #Define offsets of robot 1
        x1=-1145+16#-20
        y1=-6*25.4-9#+18-2.2556
        z1=-407.6+0.16
        a1=0
        b1=0
        c1=0.0
        
        #Define offsets of robot 2
        x2=-1145+12.4-6.5-12.975
        y2=6*25.4+10-.42-3-11.383
        z2=-409.6+10+3.8-8.006
        a2=0
        b2=0
        c2=0.0
        
        x1=-1145+16-22.75#-20
        y1=-6*25.4-9-2.067#+18-2.2556
        z1=-407.6+0.16+2.5524
        a1=0
        b1=0
        c1=0.0
        
        #Define offsets of robot 2
        x2=-1145+12.4-6.5-12.975-12.75
        y2=6*25.4+10-.42-3-11.383+0.644251028441243+1.2456
        z2=-409.6+10+3.8-8.006+2.896558268590142
        a2=0
        b2=0
        c2=0.0
        
        FD=1000 #focal length
        YLT1=0*25.4+ y1#crystal distance from base
        YLT2=-0*25.4+ y2 #crystal distance from base
        
        bt1=np.rad2deg(np.arcsin(YLT1/FD)) #b tool angle offset for robot 1
        bt2=np.rad2deg(np.arcsin(YLT2/FD)) #b tool angle offset for robot 1
        
        
        LEDT1.setPose(Pose(0,0,0,0,90,0))
        LEDT2.setPose(Pose(0,0,0,0,90,0))
        
        LEDT1.setPose(LEDT1.Pose()*rotx(np.deg2rad(-bt1)))
        LEDT2.setPose(LEDT2.Pose()*rotx(np.deg2rad(-bt2)))
        
        XLT1=-np.sqrt(FD**2-YLT1**2) #x tool offset for robot 1
        XLT2=-np.sqrt(FD**2-YLT2**2)  #x tool offset for robot 2
        
        TPP1.setPose(Pose(XLT1,YLT1,zPR/2,90,90-bt1,-90))
        TPP2.setPose(Pose(XLT2,YLT2,zPR/2,90,90-bt2,-90))
        
        RB1.setPose(Pose(x1,y1,z1,a1,b1,c1))#set postion of robot 1 base
        RB2.setPose(Pose(x2,y2,z2,a2,b2,c2))#set postion of robot 2 base
        
        BI.setParent(TPP1)
        BI.setPose(Pose(-zPR/2,0,1000,0,0,0))
        
class Webserver:
#webserver returns json string that looks like {"speed":25,"goalpos":140,"zCurrent":140.09935,"PRInt":4,"ledStatus":false,"zero":1}
    def __init__(self,IP,extenstion)->None:
        self.IP=IP #x offset
        self.ex=extenstion #y offset
        
        
    def patch(self, obj, value):
        url = "http://" + self.IP + "/" +self.ex #address of server
        obj=str(obj)
        
        payload = {obj: value} #rewrite the obj to value
        headers = {"Content-Type": "application/json"} #formatting the message as a json string don't change this

        response = requests.patch(url, data=json.dumps(payload), headers=headers) #patch in or change the value

        if response.status_code == 200: #if 204 which is http good patch
            print(obj + "Updated successfully")
        else:
            print(obj + "Error: ", response.status_code)
        return 
            
    def get(self, obj):
        url = "http://" + self.IP + "/" +self.ex #address of server
        obj=str(obj)
        
        response = requests.get(url)
        if response.status_code == 200:
            value=response.json()[obj]
            print(obj, value)
        else:
           print("Error: ", response.status_code)
        
        return value
    
    def Zero(self):
        self.patch("speed","100")
        self.patch("zero",1) 
        
    def EJog(self,pos):
        self.patch("speed","100")
        self.patch("goalposLED",pos)

class MyThreads:
    def __init__(self,threads):
        self.state=[0]*threads
        self.threads=[]
    
    def RScan(self,Robot,Rlist,t):
        Robot.Name.link=Robolink()
        i=0  
        while i< len(Rlist):
            time.sleep(0.001)
            if self.state[t]==0:
                
                Robot.SafeJM(Robot.LEDT,Rlist[i])
                #print('thread'+str(t),i)
                #print('list reads: '+str(Rlist[i]))
                self.state[t]=1
                i+=1
                #print(self.state)
                if self.state==[1]*len(self.state):
                    self.state=[0]*len(self.state)
                
            # elif self.state[t]==1:
            #     print('Race Prevented')
                
    def GScan(self,ESP,Rlist,t):
        i=0  
        while i< len(Rlist):
            time.sleep(0.001)
            if self.state[t]==0:
                ESP.patch("goalposLED",Rlist[i])
                ESP.get("zCurrent")
                time.sleep(2)
                
                #print('thread'+str(t),i)
                #print('list reads: '+str(Rlist[i]))
                self.state[t]=1
                i+=1
                #print(self.state)
                if self.state==[1]*len(self.state):
                    self.state=[0]*len(self.state)
                
            # elif self.state[t]==1:
            #     print('Race Prevented')        
                          
                
    def SyncThreads(self,fun,args):
        for t in range(0,len(args)):  
        
            thread = threading.Thread(target=fun[t],args=args[t])
            thread.daemon = True
            thread.start()
            self.threads.append(thread)
        
        print('Threads started')
        
        for t in self.threads:
            t.join()
        
        print("Done!")

#%%

ESP=Webserver('192.168.0.99','values')  
R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','LEDT1','Home1',ESP)

R1.BuildRowScan(100,400,20)

#%%
########################
#peak detection doesnt find max which is introducing error
#a peak is being detected but not the right peak
#maybe peak detection but find the max height peak
#need consider changing to max indexing
#will have trouble with plateaus but should be better than whats happening now
#

ESP=Webserver('192.168.0.99','values')  
  
R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','LEDT1','Home1',ESP)
R2=Robot('2 Mecademic Meca500 R3','192.168.0.101','2 CrystalToolSample','TPP2','APP2','LEDT2','Home2',ESP)


ESP.patch("ledStatus",False)
R1.Connect_and_Home()
R2.Connect_and_Home()

zi=100
zf=500
points=10
Rlist1=R1.BuildRowScan(zi,zf,points)
Rlist2=R2.BuildRowScan(zi,zf,points)
Glist=np.linspace(zi,zf,points)
#%%
T=MyThreads(3)
fun=[T.RScan,T.RScan,T.GScan]
args=[[R1,Rlist1,0],[R2,Rlist2,1],[ESP,Glist,2]] 

ESP.patch("ledStatus",True)

T.SyncThreads(fun,args)

ESP.patch("ledStatus",False)
#R1.Connect_and_Home()
#R2.Connect_and_Home()
#%%
ESP=Webserver('192.168.0.99','values')  
  
R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','LEDT1','Home1',ESP)
R2=Robot('2 Mecademic Meca500 R3','192.168.0.101','2 CrystalToolSample','TPP2','APP2','LEDT2','Home2',ESP)


ESP.patch("ledStatus",False)
R1.Connect_and_Home()
R2.Connect_and_Home()

zi=150
zf=300
points=10
Rlist1=R1.BuildRowScan(zi,zf,points)
Rlist2=R2.BuildRowScan(zi,zf,points)
Glist=np.linspace(zi,zf,points)
#%%
T=MyThreads(2)
fun=[T.RScan,T.RScan]
args=[[R1,Rlist1,0],[R2,Rlist2,1]] 

ESP.patch("ledStatus",True)

T.SyncThreads(fun,args)

ESP.patch("ledStatus",False)
R1.Connect_and_Home()
R2.Connect_and_Home()

#%%
ESP=Webserver('192.168.0.99','values')  
  
R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','LEDT1','Home1',ESP)
R2=Robot('2 Mecademic Meca500 R3','192.168.0.101','2 CrystalToolSample','TPP2','APP2','LEDT2','Home2',ESP)


ESP.patch("ledStatus",False)
R1.Connect_and_Home()
R2.Connect_and_Home()

zi=100
zf=300
points=20

Error=[[100., 250., 400.],
  [(-0.43206826952948507, 0.13363120837748896, 4.0),
   (0.029686906117525898, -1.253417329681019, 3.0),
   (-0.615752255473664, -2.6427825445537203, 3.0)],
  [(1.2983878132051547, -3.0491394495704838, 3.0),
   (0.8371217487513968, -2.151956559525346, 2.0),
   (0.0, 0.0, 3.0)]]

ERL1=R1.ErrorReg(Error,1)
ERL2=R1.ErrorReg(Error,2)
#print(ERL[1])

Rlist1=R1.BuildRowScanADV(ERL1,zi,zf,points)
Rlist2=R1.BuildRowScanADV(ERL2,zi,zf,points)


#Rlist1=R1.BuildRowScan(zi,zf,points)
#Rlist2=R2.BuildRowScan(zi,zf,points)
Glist=np.linspace(zi,zf,points)

#%%
T=MyThreads(3)
fun=[T.RScan,T.RScan,T.GScan]
args=[[R1,Rlist1,0],[R2,Rlist2,1],[ESP,Glist,2]] 
ESP.patch("ledStatus",True)

T.SyncThreads(fun,args)

ESP.patch("ledStatus",False)
R1.Connect_and_Home()
R2.Connect_and_Home()
