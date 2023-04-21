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
from scipy.signal import find_peaks, peak_prominences
from scipy.stats import linregress
import json
import time

class Spiral:
    def __init__(self,radius:float, spacing:float,points:float,x0:float,y0:float) -> None:
        self.x0=x0 #x offset
        self.y0=y0 #y offset
        self.points=points
        self.radius=radius
        self.spacing=spacing
        
    def build(self):
        revs= self.radius/ self.spacing            #calc number of revolutions
        k= self.spacing/(2*np.pi)                       
        theta=np.array(np.linspace(0,revs*2*np.pi,self.points)) #how many points and angle coord
        r=k * theta                       #radius coord
        x=(r*np.cos(theta)+self.x0)       #convert to x cartesian
        y=(r*np.sin(theta)+self.y0)       #convert to y cartesian
        return x, y 
    
class Linear: 
    def __init__(self,data,points,prominence,width,height):
        self.data=data #data
        self.points=points #number of points to fit
        self.prominence=prominence
        self.width=width
        self.height=height
        
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
    
    
    def getC(self,x,s_center):
        xl,yl,zl=s_center
        peaks= find_peaks(self.data['Intensity'], prominence=[self.prominence],width=self.width,height=self.height) #find peaks
        #IPeaks=Iv.index(max(Iv))
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
    
    def Zero(self):
        self.patch("speed","100")
        self.patch("zero",1) 
        
    def EJog(self,pos):
        self.patch("speed","100")
        self.patch("goalposLED",pos)

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
  
class Robot:
    def __init__(self,Name,IP,Tool,TPP,APP,Home,ESP):
        self.Name=RDK.Item(Name)
        self.IP=IP
        self.Tool=RDK.Item(Tool)
        self.TPP=RDK.Item(TPP)
        self.APP=RDK.Item(APP)
        self.Home=RDK.Item(Home)
        # Error=Pose_2_KUKA(self.APP.Pose())
        # self.Ex=Error[0]
        # self.Ey=Error[1]
        # self.Ez=Error[2]
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
    print(R)

    coord=(0,0,0)

    s_prom=10
    s_width=2
    s_height=90
    
    ESP.patch('ledStatus',True)
    for L in (range(1,loops+1)):
        cur_spiral=Spiral(10,1/L,200,0,0)
        spiral=cur_spiral.build()
        
        Lstep=0.5/L**2
        Pstep=1/L**2
        
        print('Loop'+str(L)+'On Robot'+R)
        
        
        coord,s_center,Sdata=R.Spiral_scan(spiral,coord,s_prom,s_width,s_height)
        coord,Ldata=R.Linear_scan(Sdata, Lstep, coord,s_center,s_prom,s_width,s_height)
        coord=R.PushPull_scan(Ldata, Pstep, coord,3,s_prom,s_width,s_height)
        
    ESP.patch('ledStatus',False)    
    R.Connect_and_Home()
        
        