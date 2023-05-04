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
from scipy.ndimage import uniform_filter1d
import json
import time
from art import *

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
    
    
    def getC(self,d,s_center):
        xl,yl,zl=s_center
        peaks= find_peaks(self.data['Intensity'], prominence=[self.prominence],width=self.width,height=self.height) #find peaks
        
        i_peaks,_=find_peaks(self.data['Intensity']) #index of all peak positions
        IPeaks=i_peaks[np.argmax(self.data['Intensity'][i_peaks])] #find max
        
        
        #IPeaks=Iv.index(max(Iv))
        XPP=[xl]  #use origin of previous spiral as a point
        YPP=[yl]
        XPP.append(self.data['X'][IPeaks])
        YPP.append(self.data['Y'][IPeaks])
        LRR=linregress(XPP,YPP) #run regression
        if d >= 0:
            T=1
        else:
            T=-1
        
        m=LRR[0]
        b=LRR[1]
        
        #print(m,b,d)
        x=(-(2*m*b)+T*(np.sqrt((4*m**2*b**2)-4*(1+m**2)*(b**2-d**2))))/(2*(1+m**2))
        
        
        
        return x,LRR[0]*x+LRR[1]
    
    def getx(self,d,H,m,b):
        xl,yl,zl=s_center
        
        
        x=(-(2*m*b)+H*(np.sqrt((4*m**2*b**2)-4*(1+m**2)*(b**2-d**2))))/(2*(1+m**2))
        
        
        
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
        self.Name.setSpeed(-1,50)
        if self.Name.MoveL_Test(self.JPos,Final_Pose) == 0:         # linear move to the approach position
            
            self.Name.MoveL(Final_Pose)
            self.JPos=self.Name.Joints()
        else :
            print("collision avoided at "+ str(trans_coord))
            
    def SafeJM(self,target,trans_coord):
        Final_Pose=target.Pose()*transl(trans_coord)
        self.Name.setSpeed(-1,50)
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
        tprint("SPIRAL: \n-STARTED-",font='larry3d')
        print('Spiral: Started')
        print(coord)
        
        self.PlotSpiral(spiral,'Current Spiral')
        
        x,y=spiral
        s_center=coord
        Iv=[0]
        X=[0]
        Y=[0]
        Z=[0]
        xl,yl,zl=coord
        for i in range(0,len(x)): #loop to get coord and corresponding intensity
            if len(find_peaks(Iv,prominence=[prominence],width=width,height=height)[0])<=0: #check if we found a peak
                if self.OvershootCheck(Iv,3,100)==1 and ESP.get("PRInt")<.80*max(Iv):
                    print('Spiral moved past alignment')
                    coord=(xl,yl,zl)
                    
                    self.SafeJM(self.TPP,coord)
                    break   
            
                coord=(x[i],y[i],zl)
                print(coord)
                X.append(x[i])
                Y.append(y[i])
                Z.append(0)
                
                self.SafeJM(self.TPP,coord)
                self.Name.WaitFinished()
                
                Ic=ESP.get("PRInt")
                #Iv.append(np.array(Ic))
                Iv.append(Ic)
                #self.Plot2DS(Y,Iv,'Spiral Scan')
            else:

                self.Plot2DS(Y,Iv,'Spiral Scan')
                self.Plot3DS(X,Y,Iv,prominence,width,height)
                break
        
        
        #pandas dictionary
        titled_columns={'X': X,'Y': Y,'Z': Z,
                            'Intensity': Iv}
        data = pd.DataFrame(titled_columns)
        
        #peaks= find_peaks(data['Intensity'],height=height,prominence=[prominence],width=width) 
        #IPeaks=peaks[0][0]
    
        i_peaks,_=find_peaks(data['Intensity'])
        IPeaks=i_peaks[np.argmax(data['Intensity'][i_peaks])]  
    
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
        print('Spiral Complete')
        tprint("SPIRAL: \n-COMPLETED-",font='larry3d')
        print(coord)
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
        
    def Plot2DS(self,Y,I,title):
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
    
    
    
    
    def Linear_scan(self,data,step,coord,s_center,prominence,width,height):
        print(coord)
        
        
        
        tprint("LINEAR: \n-STARTED-",font='larry3d')
        global Il,PRS
        xl,yl,zl=coord
        sxl,syl,szl=s_center
        
        i_peaks,_=find_peaks(data['Intensity'])
        IPeaks=i_peaks[np.argmax(data['Intensity'][i_peaks])]
         
        CL=Linear(data,2,prominence,width,height) #current linear line definition
        
        xL=data['X'][IPeaks]
        yL=data['Y'][IPeaks]
      
        Xl=[xL]  
        Yl=[yL]
        PRS=data['Intensity'][IPeaks]
        Il=[PRS]
        
        d=np.sqrt((xl-sxl)**2+(yl-syl)**2)
        D=[d]
        if  data['X'][IPeaks] > xl:  #find direction to move
            step=-step
            d=-d
        else:
            step=step
            d=d
        
        print('Linear move: STARTED')
        
        #Linear search
        while len(find_peaks(Il,prominence=prominence,width=width,height=height)[0]) == 0: #linear move
            if coord == s_center:
                print('Spiral on Center')
                break
        
            xL,yL=CL.getC(d,s_center)
            coord=(xL,yL,zl)
            print(coord)
            self.SafeLM(self.TPP,coord)
            
            # if input('before step (y/n)') != 'y':
            #     ESP.patch('ledStatus',False)
            #     break
            
            d +=step
        
            xL,yL=CL.getC(d,s_center)
            coord=(xL,yL,zl)
            print(coord)
            
            
            self.SafeLM(self.TPP,coord)
            
            PRLT=[]
            for j in range(3):
                #PRL=PR.get()
                PRL=ESP.get("PRInt")
                PRLT.append(PRL)
                
            LI=np.average(PRLT)
            
            #PRS=TestInt.get(xL,yL)
            Xl.append(xL)
            Yl.append(yL)
            Il.append(LI)
            D.append(d)
            
            if self.OvershootCheck(Il,3,4)==1:
                print('Linear moved past alignment')
                break
            else:
                print('Linear not past alignment')
                
            
            # if input('step taken (y/n)') != 'y':
            #     ESP.patch('ledStatus',False)
            #     break
            
            print('length'+str(len(find_peaks(Il,prominence=[prominence],width=width,height=height)[0])))
            # if len(find_peaks(-np.array(Il),plateau_size=4,prominence=30)[0])<=0:
            #     print('Linear moved past alignment')
            #     break
                
        #print('Linear move '+str(L)+ ': COMPLETED')
        
        LLdata={'X':Xl,'Y':Yl,'D':D,'I':Il}
        Ldata=pd.DataFrame(LLdata)
        
        self.Plot2DS(Ldata['D'],Ldata['I'],'Linear Scan')
        
        Lmax_I=Il.index(max(Ldata['I']))
        xl=Xl[Lmax_I]
        yl=Yl[Lmax_I]
        
        coord=(xl,yl,zl)
        #plot goes here
        
        
        tprint("LINEAR: \n-COMPLETED-",font='larry3d')
        print(coord)
        return coord, Ldata
 
    def PushPull_scan(self,data,step,coord,samples,prominence,width,height):
        tprint("PUSH PULL: \n-STARTED-",font='larry3d')
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
        self.Plot2DS(ZLPP,IPP,'Push Pull Test')
    
        if IPP.index(max(IPP))>0:
            step=step
        else:
            step=-step
    
        coord=(xl,yl,zl)
        #robot.MoveJ(target.Pose()*transl(coord))
        
        self.SafeJM(self.TPP,coord) 
        
        ZLPP=[zl]
        #IPP=[ESP.get("PRInt")]
        IPP=[0]

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
                
                if self.OvershootCheck(IPP,3,5)==1:
                    print('PushPull moved past alignment')
                    break
            else:
                break
        print('Push pull translation:Complete')
        print(1)
        zl=ZLPP[IPP.index(max(IPP))]
        print(2)
        
        self.Plot2DS(ZLPP,IPP,'Push Pull Scan')
        
        #plt.plot(ZLPP,IPP)
        #plt.show()

        coord=(xl,yl,zl)

        self.SafeJM(self.TPP,coord)
        print('Done')
        print(coord)
        tprint("PUSH PULL: \n-COMPLETED-",font='larry3d')
        return coord
    
    def OvershootCheck(self,data,window,threshold):
    #overshoot check evaluates if the system moved past alighnment by counting instances of 0 slope dI/dx
       AIdx=uniform_filter1d(np.diff(data),window) 
       AIdx=AIdx.tolist()
       AIdx=np.round(AIdx).tolist()
       # print('averaged intensity slope',str(AIdx))
       # plt.plot(AIdx)
       # plt.show()
       print('Intensity slope zeros:',str(AIdx.count(0)))
       if AIdx.count(0)<threshold:
           state=0
       else:
           state=1    
           return state
        
    def PlotSpiral(self,spiral,title):
        x,y=spiral
        plt.plot(x,y)
        plt.axis('equal')
        plt.title(title)
        plt.xlabel('X Translation mm')
        plt.ylabel('Y Translation mm')
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
        
    def SSAP_ALL(self,Robots,PR_pos,loops,tests):
        self.SetTargets(PR_pos)
        global Sdata, Ldata
        Time=[]
        for t in (range(1,tests+1)):
            
            ESP=Webserver('192.168.0.99','values')

            ESP.patch('ledStatus',False)
            if ESP.get('zero')==0:#If gantry hasn't been zeroed zero it
                ESP.Zero()
            ESP.EJog(PR_pos) #jog gantry to position for alignment
            ESP.get('zCurrent') #wait until the gantry move has stopped and can process a get before continuing
            
            
            for R in Robots:
                R.Connect_and_Home()

            
            Robots=[R1,R2]
            
            for R in Robots:
                print(R)
            
                coord=(0,0,0)
            
                s_prom=10
                s_width=1
                s_height=50
                
                ESP.patch('ledStatus',True)
                tic()
                for L in (range(1,loops+1)):
                    
                    cur_spiral=Spiral(10,1/L,200*L,coord[0],coord[1])
                    spiral=cur_spiral.build()
                    
                    Lstep=0.5/L**2
                    Pstep=1/L**2
                    
                    print('Loop'+str(L))
                    
                    print('start'+str(coord))
                    coord,s_center,Sdata=R.Spiral_scan(spiral,coord,s_prom,s_width,s_height)
                    coord,Ldata=R.Linear_scan(Sdata, Lstep, coord,s_center,s_prom,s_width,s_height)
                    coord=R.PushPull_scan(Ldata, Pstep, coord,3,s_prom,s_width,s_height)
                    print('end'+str(coord))
                time=toc()
                Time.append(time)
                print('Time'+str(time))
                if input('Do You Want To Continue? (y/n)') != 'y':
                    ESP.patch('ledStatus',False)
                    break
                ESP.patch('ledStatus',False)    
                R.Connect_and_Home()
                
            if input('Do You Want To Continue? (y/n)') != 'y':
                break
            return Time
    
    def SSAP(self,PR_pos,loops,tests):
        self.SetTargets(PR_pos)
        
        if ESP.get('zero')==0:#If gantry hasn't been zeroed zero it
            ESP.Zero()
        ESP.EJog(PR_pos) #jog gantry to position for alignment
        ESP.get('zCurrent') #wait until the gantry move has stopped and can process a get before continuing
        
        Time=[]
        coord=(0,0,0)
        for t in (range(1,tests+1)):
            ESP.patch('ledStatus',False)
            if ESP.get('zero')==0:#If gantry hasn't been zeroed zero it
                ESP.Zero()
            ESP.EJog(PR_pos) #jog gantry to position for alignment
            ESP.get('zCurrent') #wait until the gantry move has stopped and can process a get before continuing
            
            
            self.Connect_and_Home()

            s_prom=5
            s_width=1
            s_height=50
            
            ESP.patch('ledStatus',True)
            tic()
            for L in (range(1,loops+1)):
                
                cur_spiral=Spiral(10,1/L,200*L,coord[0],coord[1])
                spiral=cur_spiral.build()
                
                Lstep=0.5/L**2
                Pstep=1/L**2
                
                print('Loop'+str(L))
                
                print('start'+str(coord))
                coord,s_center,Sdata=self.Spiral_scan(spiral,coord,s_prom,s_width,s_height)
                coord,Ldata=self.Linear_scan(Sdata, Lstep, coord,s_center,s_prom,s_width,s_height)
                coord=self.PushPull_scan(Ldata, Pstep, coord,3,s_prom,s_width,s_height)
                print('end'+str(coord))
            time=toc()
            
            Time.append(time)
            print('Time'+str(time))
            ESP.patch('ledStatus',False)    
            self.Connect_and_Home()
            

            return Time,coord
        
    def SSAP_ALL_I(self,Robots,PR_pos,loops,tests):
        global Error
        Error=[ [] for _ in range(3)]
        
        zi=PR_pos[0]
        zf=PR_pos[1]
        
        ZL=np.linspace(zi,zf,tests)
        
        Error[0].append(ZL)
        Time=[]

        ESP.patch('ledStatus',False)   
        
        for R in Robots: #home everyone
        
            R.Connect_and_Home()
        i=0
            
        for R in Robots:
            print('Starting robot call')
            i+=1
            print('SLAP on'+ str(R))
            
        
            coord=(0,0,0)
            
            tic()
            for Etests in (range(0,tests)):
                
                self.SetTargets(ZL[Etests])
                    
                time, coord = R.SSAP( ZL[Etests] , loops, 1 )
                print(coord)
                #Everything to do after an alignment
                
                Error[i].append(coord)
                print(Error)
                time=toc()
                Time.append(time)
                print('Time'+str(time))
                #-----------------------------------------
                # if input('Do You Want To Continue? (y/n)') != 'y':
                #     ESP.patch('ledStatus',False)
                #     break
                
            print(0)
    
        return Time,Error
                
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
        
        TPP1.setParent(LED)
        TPP2.setParent(LED)
        
        # Define the position of the photoresistor
        zPR=-PR_pos
        #set pose of Hardware
        PR.setPose(Pose(0,0,zPR,0,0,0)) #photo Resistor
        CT.setPose(Pose(0,0,zPR/2,0,0,0)) #crystal Target point
        
        
        #Define offsets of robot 1
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
        
        XLT1=-np.sqrt(FD**2-YLT1**2) #x tool offset for robot 1
        XLT2=-np.sqrt(FD**2-YLT2**2)  #x tool offset for robot 2
        
        TPP1.setPose(Pose(XLT1,YLT1,zPR/2,90,90-bt1,-90))
        TPP2.setPose(Pose(XLT2,YLT2,zPR/2,90,90-bt2,-90))
        
        RB1.setPose(Pose(x1,y1,z1,a1,b1,c1))#set postion of robot 1 base
        RB2.setPose(Pose(x2,y2,z2,a2,b2,c2))#set postion of robot 2 base
        
        BI.setParent(TPP1)
        BI.setPose(Pose(-zPR/2,0,1000,0,0,0))
        
def printSLAP():
    print(
    '   _____ _               _____             _      _____ _____ _   _ __  __ ______ _   _ _______ \n'
    '  / ____| |        /\   |  __ \      /\   | |    |_   _/ ____| \ | |  \/  |  ____| \ | |__   __|\n'
    ' | (___ | |       /  \  | |__) |    /  \  | |      | || |  __|  \| | \  / | |__  |  \| |  | |   \n'
    '  \___ \| |      / /\ \ |  ___/    / /\ \ | |      | || | |_ | . ` | |\/| |  __| | . ` |  | |   \n'
    '  ____) | |____ / ____ \| |       / ____ \| |____ _| || |__| | |\  | |  | | |____| |\  |  | |   \n'
    ' |_____/|______/_/    \_\_|      /_/    \_\______|_____\_____|_| \_|_|  |_|______|_| \_|  |_|   \n'
    '                                                                                                \n'
                                                                                                    
    )




def printBAR():
    print(
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                 
                                                                                                                                                                                                                                     
     ' ______ ______ ______ ______ ______ ______ ______ ______ ______ ______ ______ ______ ______    \n'
     '|______|______|______|______|______|______|______|______|______|______|______|______|______| \n'
     '\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \/ /\ \n'
     ' >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  <  >  \n'
     '/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/\_\/_/ \n'
     '|______|______|______|______|______|______|______|______|______|______|______|______|______| \n'
     '\n'                                                                                                                                                            
    )
    
    
    
def printLOGO(): 
       print(
    '   _____ __    ___   ______   _  _____ ____________________  _  _____   __  \n'
    '  / ___// /   /   | / ____/  / |/ / _ /_  __/_  __/  _/ __ \/ |/ / _ | / / \n'
    '  \__ \/ /   / /| |/ /      /    / __ |/ /   / / _/ // /_/ /    / __ |/ /__ \n'
    ' ___/ / /___/ ___ / /___   /_/|_/_/_|_/_/___/_/_/___/\____/_/|_/_/ |_/____/__  ___  \n'
    '/____/_____/_/  |_\____/     / _ |/ ___/ ___/ __/ /  / __/ _ \/ _ /_  __/ __ \/ _ \ \n'
    '                            / __ / /__/ /__/ _// /__/ _// , _/ __ |/ / / /_/ / , _/ \n'
    ' ____________       ___    /_/_|_\___/\___/___/____/___/_/|_/_/_|_/_/__\____/_/|_|  \n'
    '|\____________\    |\__\   / /  / _ | / _ )/ __ \/ _ \/ _ /_  __/ __ \/ _ \ \/ /  \n'
    '\|____________|    \|__|  / /__/ __ |/ _  / /_/ / , _/ __ |/ / / /_/ / , _/\  /  \n'
    '                         /____/_/ |_/____/\____/_/|_/_/ |_/_/  \____/_/|_| /_/  \n'
                            
    )
                           


    

#%%
ESP=Webserver('192.168.0.99','values')
R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','LEDT1','Home1',ESP)
R2=Robot('2 Mecademic Meca500 R3','192.168.0.101','2 CrystalToolSample','TPP2','APP2','LEDT2','Home2',ESP)

#R1.SetTargets(300)
ESP.patch('ledStatus',False)
R1.Connect_and_Home()
R2.Connect_and_Home()
#%%
Time,EE=R1.SSAP(400,1,1)

#%%
printLOGO()
printBAR()
printSLAP()
printBAR()

PR_pos=[100,400]

Robots=[R1,R2]
Time,Error=R1.SSAP_ALL_I(Robots, PR_pos, 1, 3)

#%%
ESP=Webserver('192.168.0.99','values')
R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','LEDT1','Home1',ESP)
R2=Robot('2 Mecademic Meca500 R3','192.168.0.101','2 CrystalToolSample','TPP2','APP2','LEDT2','Home2',ESP)

#R1.SetTargets(300)
#ESP.patch('ledStatus',False)
#R1.Connect_and_Home()
#R2.Connect_and_Home()

E=[[100., 250., 400.],
  [(-0.43206826952948507, 0.13363120837748896, 4.0),
   (0.029686906117525898, -1.253417329681019, 3.0),
   (-0.615752255473664, -2.6427825445537203, 3.0)],
  [(1.2983878132051547, -3.0491394495704838, 3.0),
   (0.8371217487513968, -2.151956559525346, 2.0),
   (0.0, 0.0, 3.0)]]
 

i=0
ESP.patch("goalposLED", E[0][i])
coord=E[1][i]
#coord=(0,0,0)

R1.SetTargets(E[0][i])

R1.SafeJM(RDK.Item('TPP1'), coord)
ESP.patch('ledStatus',True)


#%%


# printLOGO()
# printBAR()
# printSLAP()
# printBAR()

# Robots=[R1,R2]
# loops=2
# tests=1
# Time=R1.SSAP_ALL(Robots,300,loops,tests)

#%%



#%%

# loops=3
# ESP=Webserver('192.168.0.99','values')
# R1=Robot('1 Mecademic Meca500 R3','192.168.0.100','1 CrystalToolSample','TPP1','APP1','LEDT1','Home1',ESP)
# R2=Robot('2 Mecademic Meca500 R3','192.168.0.101','2 CrystalToolSample','TPP2','APP2','LEDT2','Home2',ESP)

# ESP.patch('ledStatus',False)
# if ESP.get('zero')==0:#If gantry hasn't been zeroed zero it
#     ESP.Zero()
# ESP.EJog(150) #jog gantry to position for alignment
# ESP.get('zCurrent') #wait until the gantry move has stopped and can process a get before continuing

# R1.Connect_and_Home()
# R2.Connect_and_Home()

# Robots=[R1,R2]

# for R in Robots:
#     print(R)

#     coord=(0,0,0)

#     s_prom=10
#     s_width=2
#     s_height=50
    
#     ESP.patch('ledStatus',True)
#     for L in (range(1,loops+1)):
#         cur_spiral=Spiral(10,1/L,200,0,0)
#         spiral=cur_spiral.build()
        
#         x,y=spiral
#         plt.plot(x,y)
#         plt.axis('equal')
#         plt.show()
#         Lstep=0.5/L**2
#         Pstep=1/L**2
        
#         print('Loop'+str(L))
        
        
#         coord,s_center,Sdata=R.Spiral_scan(spiral,coord,s_prom,s_width,s_height)
#         coord,Ldata=R.Linear_scan(Sdata, Lstep, coord,s_center,s_prom,s_width,s_height)
#         coord=R.PushPull_scan(Ldata, Pstep, coord,3,s_prom,s_width,s_height)
        
#     ESP.patch('ledStatus',False)    
#     R.Connect_and_Home()
        