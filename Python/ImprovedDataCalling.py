# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 23:05:45 2023

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

class Spiral:
    def __init__(self,x0:float,y0:float,points:float) -> None:
        self.x0=x0 #x offset
        self.y0=y0 #y offset
        self.points=points
        
    def build(self, radius:float, spacing:float):
        revs=radius/spacing            #calc number of revolutions
        k=spacing/(2*np.pi)                       
        theta=np.array(np.linspace(0,revs*2*np.pi,self.points)) #how many points and angle coord
        r=k * theta                       #radius coord
        x=(r*np.cos(theta)+self.x0)       #convert to x cartesian
        y=(r*np.sin(theta)+self.y0)       #convert to y cartesian
        return x, y 
    
class Linear: 
    def __init__(self,data,points):
        self.data=data #data
        self.points=points #number of points to fit
        
    def get(self,x): 
        peaks= find_peaks(self.data['Intensity'], prominence=[prominence],width=width,height=height)
        IPeaks=peaks[0]
        H=peaks[1]['peak_heights']
        A=[]
        for i in range(len(IPeaks)):
            corr=[IPeaks[i],H[i]]
            A.append(corr)
        SA=sorted(A, key=lambda x: x[1], reverse=True) #sort by peak height largest to smallest
        if len(SA) <= 1:
             
            E=np.sqrt((Ix-xl)**2+(Iy-yl)**2)
            print('***Linear Translation********************************************************')
            print(Ldata)
            print('_________SSAP Report:_________')
            print('Calculated peak is at',xl,yl)
            print('True peak is at',Ix,Iy)
            print('Error is:',E)
            print('Linear step size is:',step)
            raise Exception("Spiral did not find enough peaks")
        XPP=[]
        YPP=[]
        for i in range(self.points):
            #ICP.append(SA[i][0])
            XPP.append(self.data['X'][SA[i][0]])
            YPP.append(self.data['Y'][SA[i][0]])
        LRR=linregress(XPP,YPP)

        return LRR[0]*x+LRR[1]
    def getC(self,x,xl,yl):
        peaks= find_peaks(self.data['Intensity'], prominence=[prominence],width=width,height=height) #find peaks
        IPeaks=Iv.index(max(Iv))
        XPP=[xl]  #use origin of previous spiral as a point
        YPP=[yl]
        XPP.append(self.data['X'][peaks[0][0]])
        YPP.append(self.data['Y'][peaks[0][0]])
        LRR=linregress(XPP,YPP) #run regression
        return LRR[0]*x+LRR[1]


    def getSA(self):
        peaks= find_peaks(self.data['Intensity'],prominence=[prominence],width=width,height=height)
        IPeaks=peaks[0]
        H=peaks[1]['peak_heights']
        A=[]
        for i in range(len(IPeaks)):
            corr=[IPeaks[i],H[i]]
            A.append(corr)
        SA=sorted(A, key=lambda g: g[1], reverse=True) #sort by peak height largest to smallest
        return SA,IPeaks

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

        if response.status_code == 204: #if 204 which is http good patch
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


class Intensity: #intensity function is a gaussian curve 
    def __init__(self,x0:float,y0:float, mu_x, mu_y, sigma_x, sigma_y):
        self.x0=-x0 #x offset
        self.y0=-y0 #y offset
        self.mu_x=mu_x
        self.mu_y=mu_y
        self.sigma_x=sigma_x #x squish/stretch
        self.sigma_y=sigma_y #y sqush/stretch
        
        
    def get(self,x, y): #gaussian curve definition
        I=np.exp(-0.5 * (((x +self.x0- self.mu_x) / self.sigma_x)**2 + ((y+self.y0 - self.mu_y) / self.sigma_y)**2))
        return I


# Robot Initialization

TPP= RDK.Item('TPP',ITEM_TYPE_TARGET)
LED= RDK.Item('LED',ITEM_TYPE_FRAME)
LEDT= RDK.Item('LEDT',ITEM_TYPE_TARGET)
LEDT2=RDK.Item('LEDT2')

RCS= RDK.Item('RCS',ITEM_TYPE_TARGET)
PR= RDK.Item('Photo Resistor',ITEM_TYPE_TARGET)
APP= RDK.Item('APP',ITEM_TYPE_TARGET)
#RCS.setParent(LEDT)

#Webserver Initilization
IP="192.168.0.99"
ext="values"
ESP=Webserver(IP,ext)

#Robot Home

# Set up the IP address and port number of the robot
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

def SafeLM(joints,Final_Pose):
    if self.Name.MoveL_Test(joints,Final_Pose) == 0:         # linear move to the approach position
        self.Name.MoveL(Final_Pose)
    else :
        print("collision avoided")
        
def SafeJM(joints,Final_Pose):
    if self.Name.MoveJ_Test(joints,Final_Pose) == 0:         # linear move to the approach position
        self.Name.MoveJ(Final_Pose)
    else :
        print("collision avoided")
        
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
        
    def SafeLM(self,joints,Final_Pose):
        if self.Name.MoveL_Test(joints,Final_Pose) == 0:         # linear move to the approach position
            self.Name.MoveL(Final_Pose)
        else :
            print("collision avoided")
            
    def SafeJM(self,joints,Final_Pose):
        if self.Name.MoveJ_Test(joints,Final_Pose) == 0:         # linear move to the approach position
            self.Name.MoveJ(Final_Pose)
        else :
            print("collision avoided")
            
    def spiralscan(self,spiral,xl,yl,zl):
        x,y=spiral
        for i in range(0,len(x)): #loop to get coord and corresponding intensity
            if len(find_peaks(Iv,prominence=[prominence],width=width,height=height)[0])<=0: #check if we found a peak
                if i==len(x):
                    print('No peaks found in spiral: Spiral centered')
                    coord=(xl,yl,zl)
                    
                    SafeJM(self.Name.Joints(),self.TPP.Pose()*transl(coord))
                    
                    raise Exception('Alignment Acheived')
                    break
                
                coord=(x[i],y[i],zl)
                X.append(x[i])
                Y.append(y[i])
                Z.append(0)
                #Ic=PR.get() #get intensity for PR for online run

    
                #self.Name.MoveL(self.TPP.Pose()*transl(coord))
                

                SafeLM(self.Name.Joints(),self.TPP.Pose()*transl(coord))
                print(Final)
                self.Name.WaitFinished()
                
                Ic=ESP.get("PRInt")
                #Ic=TestInt.get(x[i],y[i]) #test intensity get
                Iv.append(np.array(Ic))
            else:
                ps=find_peaks(Iv,prominence=[prominence],width=width,height=height)[0][0]
                fig1 = plt.figure()
                #plt.style.use('seaborn-notebook')
                theta=np.array(np.linspace(0,revs*2*np.pi,points))
                peaks, _ = find_peaks(Iv)
    
                # height=.1,distance=3,threshold=0.01,
                prominences = peak_prominences(Iv, peaks)[0]
                contour_heights = np.array(Iv)[peaks] - prominences
                plt.plot(Iv)
                plt.plot(peaks, np.array(Iv)[peaks], "x")
                plt.vlines(x=peaks, ymin=contour_heights, ymax=np.array(Iv)[peaks])
                plt.show()
    
                fig = plt.figure()
    
                ax = fig.add_subplot(111, projection='3d')
                ax.plot3D(X, Y, Iv)
                ax.plot3D(X[ps], Y[ps], Iv[ps],'x')
                
                ax.set_xlabel('X axis')
                ax.set_ylabel('Y axis')
                ax.set_zlabel('Z axis')
                plt.show() 
                break
        
        
        #pandas dictionary
        titled_columns={'X': X,'Y': Y,'Z': Z,
                            'Intensity': Iv}
        data = pd.DataFrame(titled_columns)
        print('Spiral ' +str(L)+ ' Completed')
    
        if len(find_peaks(data['Intensity'],plateau_size=plateau_size,height=200)[0])>0: #check to see if we are at the peak
            PI=find_peaks(data['Intensity'],plateau_size=plateau_size,height=height)[0]
            xl=data['X'][PI]
            yl=data['Y'][PI]
            coord=(xl,yl,zl)
            print(coord)
            #self.Name.MoveJ(self.TPP.Pose()*transl(coord))
            
            SafeJM(self.Name.Joints(),self.TPP.Pose()*transl(coord))
            
            print('SSAP Complete!')
    
       
    
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
        SafeLM(self.Name.Joints(),self.TPP.Pose()*transl(coord))
        
        print('at max of spiral')
        print('Position is', coord)
        
        return data,coord
    
    def linearscan(self,data,step,xl,yl,zl):
        peaks= find_peaks(data['Intensity'],height=height,prominence=[prominence],width=width) 
        IPeaks=peaks[0][0]
        CL=Linear(data,2) #current linear line definition
        
        xL=data['X'][IPeaks]
        yL=CL.getC(xL,xl,yl) #start y
        
        
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
            yL=CL.getC(xL,xl,yl)
            PRLT=[]
            for j in range(IAS):
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
            
            SafeLM(self.Name.Joints(),self.TPP.Pose()*transl(coord))
    
            if len(find_peaks(-np.array(Il),plateau_size=4,prominence=30)[0])<=0:
                print('Linear moved past alignment')
                break
                
        print('Linear move '+str(L)+ ': COMPLETED')
        
        LLdata={'X':Xl,'Y':Yl,'I':Il}
        Ldata=pd.DataFrame(LLdata)
        
        Lmax_I=Il.index(max(Ldata['I']))
        xl=Xl[Lmax_I]
        yl=Yl[Lmax_I]
    
        #plot goes here
        
    def PPscan(self,data,step,xl,yl,zl):
        IPP=[]
        ZLPP=[]
        print('Push pull test: Start')
        for i in range(3):
            coord=(xl,yl,zl)
            ZLPP.append(zl)
            #robot.MoveJ(target.Pose()*transl(coord))
            
            SafeJM(self.Name.Joints(),self.TPP.Pose()*transl(coord))
            
            PRPPT=[]
            for j in range(IAS):
                #PRPP=PR.get()
                PRPP=ESP.get("PRInt")
                PRPPT.append(PRPP)
            PPI=np.average(PRPPT)
            IPP.append(PPI)
    
            zl=zl+(zstep)
        print('Push pull test: Complete')
        plt.plot(ZLPP,IPP)
        plt.show()
    
        if IPP.index(max(IPP))>0:
            zstep=zstep
        else:
            zstep=-zstep
    
        coord=(xl,yl,zl)
        #robot.MoveJ(target.Pose()*transl(coord))
        
        SafeJM(self.Name.Joints(),self.TPP.Pose()*transl(coord)) 
        
        ZLPP=[zl]
        IPP=[ESP.get("PRInt")]

        print('Push pull translation:Start')
        while len(find_peaks(IPP,prominence=5,height=height)[0])<=0:  #Push pull
            if PPI > .80*max(IPP) and IPP.count(ESP.get("PRInt"))<10:
                coord=(xl,yl,zl)
                ZLPP.append(zl)
                #robot.MoveL(target.Pose()*transl(coord))
                SafeLM(self.Name.Joints(),self.TPP.Pose()*transl(coord))
                
                PRPPT=[]
                for j in range(IAS):
                    #PRPP=PR.get()
                    PRPP=ESP.get("PRInt")
                    PRPPT.append(PRPP)
                PPI=np.average(PRPPT)
                IPP.append(PPI)
                zl=zl+(zstep)
            else:
                break
        print('Push pull translation:Complete')
        print(1)
        zl=ZLPP[IPP.index(max(IPP))]
        print(2)
     
        plt.plot(ZLPP,IPP)
        plt.show()
        print(3)
    print(4)
    coord=(xl,yl,zl)
    #robot.MoveL(target.Pose()*transl(coord))
    rj=robot.Joints()
    Final = target.Pose()*transl(coord)
    SafeJM(rj,Final)  
    
    def Connect_and_Home_All(rdf):
        for nrobot in range(0,len(rdf)):
            print(rdf['Names'][nrobot])
            
            robot=self.Name
            IP=rdf['IP'][nrobot]
            robot.ConnectSafe(IP)
            
            Home='Home'+str(nrobot+1)
            Home=RDK.Item(Home)
            
            #robot.setConnectionParams(IP, Port, 'socket')
            tool=RDK.Item(rdf['Tools'][nrobot])
            robot.setPoseTool(tool) # Update the TCP
            robot.setSpeed(-1,25)  # Set linear speed in mm/s
            robot.setSpeed(25)  # Set linear speed in mm/s
            
            rj=robot.Joints()
            Final = Home.Pose()
            SafeJM(rj,Final)
        
        
        
        
        
def my3dplot(xl,yl,XL,YL,data,Ldata):
        peaks= find_peaks(data['Intensity'],height=height,prominence=[prominence],width=width) 
        IPeaks=peaks[0][0]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D(data['X'], data['Y'], data['Intensity'])
        ax.plot3D(data['X'][IPeaks], data['Y'][IPeaks], data['Intensity'][IPeaks],"x")
        ax.plot3D(xl, yl, max(Ldata['I']),"x")
        #ax.plot3D(Xl, Yl, Il,'o')
    
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        plt.show()
    






#SSAP
#connect and home
for R in (range(0,len(rdf))):
    print('R'+str(R+1))
    
    robot=Robot(rdf,R)
    robot.Connect_and_Home(rdf)







loops=2
for R in (range(0,len(rdf))):
    print('R'+str(R+1))
    
    robot=Robot(rdf,R)
    

    
    for L in (range(1,loops+1)):
        print('L'+str(L))
        robot.spiralscan(spiral, xl, yl, zl)
        robot.linearscan(data, step, xl, yl, zl)
        robot.ppscan()
        
        