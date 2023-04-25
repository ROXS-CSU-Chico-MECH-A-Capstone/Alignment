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
    def __init__(self,Name,IP,Tool,TPP,APP,Home,ESP):
        self.Name=RDK.Item(Name)
        self.IP=IP
        self.Tool=RDK.Item(Tool)
        self.TPP=RDK.Item(TPP)
        self.APP=RDK.Item(APP)
        self.Home=RDK.Item(Home)
        self.LEDT=RDK.Item(LEDT)
        # Error=Pose_2_KUKA(self.APP.Pose())
        # self.Ex=Error[0]
        # self.Ey=Error[1]
        # self.Ez=Error[2]
        self.Error=Pose_2_KUKA(self.APP.Pose()) #set error offsets
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
        
        if self.Name.MoveJ_Test(self.JPos,Final_Pose) == 0:         # linear move to the approach position
        
            self.Name.MoveJ(Final_Pose)
            self.JPos=self.Name.Joints()
        else :
            print("collision avoided")
            
    def Check_Alignment_Status(data,plateau_size,height,xl,yl,zl):
        if len(find_peaks(data['Intensity'],plateau_size=plateau_size,height=height)[0])>0: #check to see if we are at the peak
            PI=find_peaks(data['Intensity'],plateau_size=plateau_size,height=height)[0]
            xl=data['X'][PI]
            yl=data['Y'][PI]
            coord=(xl,yl,zl)
            self.SafeJM(self.TPP,coord)
            print('SSAP Complete!')
            
    def Spiral_scan(self,spiral,coord,prominence,width,height):
        x,y=spiral
        s_center=coord
        Iv=[]
        X=[]
        Y=[]
        Z=[]
        xl,yl,zl=coord
        for i in range(0,len(x)): #loop to get coord and corresponding intensity
            if len(find_peaks(Iv,prominence=[prominence],width=width,height=height)[0])<=0: #check if we found a peak
                if i==len(x):
                    print('No peaks found in spiral: Spiral centered')
                    coord=(xl,yl,zl)
                    
                    self.SafeJM(self.TPP,coord)
                    
                    raise Exception('Alignment Acheived')
                    break
                
                coord=(x[i],y[i],zl)
                X.append(x[i])
                Y.append(y[i])
                Z.append(0)
                
                
                self.SafeJM(self.TPP,coord)
                
                
                self.Name.WaitFinished()
                
                Ic=ESP.get("PRInt")
                Iv.append(np.array(Ic))
            else:
                
                
                self.Plot2DS(Y,Iv)
                self.Plot3DS(X,Y,Iv,prominence,width,height)
                
                break
        
        
        #pandas dictionary
        titled_columns={'X': X,'Y': Y,'Z': Z,
                            'Intensity': Iv}
        data = pd.DataFrame(titled_columns)
        print('Spiral ' +str(L)+ ' Completed')
    
        peaks= find_peaks(data['Intensity'],height=height,prominence=[prominence],width=width) 
        IPeaks=peaks[0][0]
    
        xL=data['X'][IPeaks]
        yL=data['Y'][IPeaks]
        #yL=CL.getC(xL,xl,yl) #start y
        
        print(data['X'][IPeaks],data['Y'][IPeaks])
        #Xl=[xL]  
        #Yl=[yL]
    
        # PRS=data['Intensity'][IPeaks]
        # Il=[PRS]
        
        coord=(xL,yL,zl)
        self.SafeLM(self.TPP,coord)
        
        print('at max of spiral')
        print('Position is', coord)
        
        return coord,s_center,data
    
    def Plot3DS(self,X,Y,I,prominence,width,height):

        ps=find_peaks(I,prominence=[prominence],width=width,height=height)[0][0]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D(X, Y, I)
        ax.plot3D(X[ps], Y[ps], I[ps],'x')
        
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        plt.show() 
        
    def Plot2DS(self,Y,I):
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
        plt.show()
    
    
    
    
    def Linear_scan(self,data,step,coord,s_center,prominence,width,height):
        xl,yl,zl=coord
        print('1',str(coord))
        peaks= find_peaks(data['Intensity'],height=height,prominence=[prominence],width=width) 
        IPeaks=peaks[0][0]
        CL=Linear(data,2,prominence,width,height) #current linear line definition
        
        xL=data['X'][IPeaks]
        
        print('precurser',str((xL,xl,yl)))
        yL=CL.getC(xL,s_center) #start y
        print(yL)
        
        Xl=[xL]  
        Yl=[yL]
        PRS=data['Intensity'][IPeaks]
        Il=[PRS]
        
        
        if  data['X'][IPeaks] > xl:  #find direction to move
            step=step
        else:
            step=-step
        
        print('Linear move: STARTED')
        
        #Linear search
        while len(find_peaks(Il,prominence=[prominence],width=width,height=height)[0])<=0: #linear move
            xL=xL+step
            yL=CL.getC(xL,s_center)
            PRLT=[]
            for j in range(3):
                #PRL=PR.get()
                PRL=ESP.get("PRInt")
                PRLT.append(PRL)
                
            LI=np.average(PRLT)
            print(PRLT)
            print(LI)
            #PRS=TestInt.get(xL,yL)
            Xl.append(xL)
            Yl.append(yL)
            Il.append(LI)
            coord=(xL,yL,zl)
            print('1',str(coord))
            
            self.SafeLM(self.TPP,coord)
    
            if len(find_peaks(-np.array(Il),plateau_size=4,prominence=30)[0])<=0:
                print('Linear moved past alignment')
                break
                
        print('Linear move '+str(L)+ ': COMPLETED')
        
        LLdata={'X':Xl,'Y':Yl,'I':Il}
        Ldata=pd.DataFrame(LLdata)
        
        Lmax_I=Il.index(max(Ldata['I']))
        xl=Xl[Lmax_I]
        yl=Yl[Lmax_I]
        print('2',str(coord))
        coord=(xl,yl,zl)
        #plot goes here
        print('3',str(coord))
        return coord, Ldata
        
    def PushPull_scan(self,data,step,coord,samples,prominence,width,height):
        xl,yl,zl=coord
        IPP=[]
        ZLPP=[]
        print('Push pull test: Start')
        for i in range(3):
            coord=(xl,yl,zl)
            ZLPP.append(zl)
            #robot.MoveJ(target.Pose()*transl(coord))
            
            self.SafeJM(self.TPP,coord)
            
            PRPPT=[]
            for j in range(samples):
                #PRPP=PR.get()
                PRPP=ESP.get("PRInt")
                PRPPT.append(PRPP)
            PPI=np.average(PRPPT)
            IPP.append(PPI)
    
            zl=zl+(step)
        print('Push pull test: Complete')
        plt.plot(ZLPP,IPP)
        plt.show()
    
        if IPP.index(max(IPP))>0:
            step=step
        else:
            step=-step
    
        coord=(xl,yl,zl)
        #robot.MoveJ(target.Pose()*transl(coord))
        
        self.SafeJM(self.TPP,coord) 
        
        ZLPP=[zl]
        IPP=[ESP.get("PRInt")]

        print('Push pull translation:Start')
        while len(find_peaks(IPP,prominence=prominence,height=height)[0])<=0:  #Push pull
            if PPI > .80*max(IPP) and IPP.count(ESP.get("PRInt"))<10:
                coord=(xl,yl,zl)
                ZLPP.append(zl)
                #robot.MoveL(target.Pose()*transl(coord))
                self.SafeLM(self.TPP,coord)
                
                PRPPT=[]
                for j in range(samples):
                    #PRPP=PR.get()
                    PRPP=ESP.get("PRInt")
                    PRPPT.append(PRPP)
                PPI=np.average(PRPPT)
                IPP.append(PPI)
                zl=zl+(step)
            else:
                break
        print('Push pull translation:Complete')
        print(1)
        zl=ZLPP[IPP.index(max(IPP))]
        print(2)
     
        plt.plot(ZLPP,IPP)
        plt.show()

        coord=(xl,yl,zl)

        self.SafeJM(self.TPP,coord)
        print('Done')
        return coord
    
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


    def RScan(self,lock,Rlist,t):
        global state
        i=0  
        while i< len(Rlist):
            time.sleep(0.1)
            if state[t]==0:
                self.safeMJ(self.LEDT,Rlist[i])
                print('thread'+str(t),i)
                state[t]=1
                i+=1
                #time.sleep(0.1)
                print(state)
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

loops=2
ESP=Webserver('192.168.0.99','values')
R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','Home1',ESP)
R2=Robot('2 Mecademic Meca500 R3','192.168.0.101','2 CrystalToolSample','TPP2','APP2','Home2',ESP)

ESP.patch('ledStatus',False)
if ESP.get('zero')==0:#If gantry hasn't been zeroed zero it
    ESP.Zero()
ESP.EJog(150) #jog gantry to position for alignment
ESP.get('zCurrent') #wait until the gantry move has stopped and can process a get before continuing



Name=RDK.Item('1 Mecademic Meca500 R3')
Name.ConnectSafe('192.168.0.100')
R1.Connect_and_Home()
R2.Connect_and_Home()

R1.SafeJM(RDK.Item('TPP1'),(0,0,0))


Robots=[R1,R2]


for R in Robots:





#%%
def Rscan(fun1,fun2,fun3,t1,t2,t3)
    lock = threading.Lock()
    state=[0,0,0]
    
    R=Robot(1,1)
    l=np.linspace(1,10,10)
    # Create two threads
    
    
    for t in range(0,threads):    
    tg = threading.Thread(target=fun1,args=(lock,l,t1))
    t1 = threading.Thread(target=fun2,args=(lock,l,t2))
    t2 = threading.Thread(target=fun3,args=(lock,l,t3))
    
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