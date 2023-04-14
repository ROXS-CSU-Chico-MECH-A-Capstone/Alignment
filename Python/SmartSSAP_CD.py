
"""
Created on Thu Jan 19 18:43:18 2023

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
#******CLASSES AND CLASS FUNCTIONS****************************************************************************
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

class Photoresistor:
    def __init__(self,url) -> None:
        self.url=url
    def get(self):
        response = requests.get(self.url) #go to webserver under the url of the arduino/photoresistor
        if response.ok: #if we can connect
            value = response.text #screen scrape the web server
            pr=re.findall(r'\d+',value) #take the number off the webserver
            pri=-int(pr[0]) #reverse the sign of the number
            #print('Received value:', value)
            #print('Number', PR)
        else:
            print('Error:', response.status_code)  
          

        return pri
    
def SafeLM(joints,Final_Pose):
    if robot.MoveL_Test(joints,Final_Pose) == 0:         # linear move to the approach position
        robot.MoveL(Final_Pose)
    else :
        print("collision avoided")
        
def SafeJM(joints,Final_Pose):
    if robot.MoveJ_Test(joints,Final_Pose) == 0:         # linear move to the approach position
        robot.MoveJ(Final_Pose)
    else :
        print("collision avoided")
        




#def Exit_SSAP(Intensity,plateau_size,height):
 #   if len(find_peaks(Intensity,plateau_size=pleatau_size,height=height)[0])>0:
    #    break

#***SPIRAL SCANNING ALIGNMENT PROCEDURE**************************************************************************************************************
#___Initialization______________
xl=0 #starting offset from target for spiral
yl=0
zl=0
Loops=3 #number of loops to run SSAP


# Set up the IP address and port number of the robot
IPs = ['192.168.0.100','192.168.0.101']
PORT = [30002,30003]
R_NAMES = ['1 Mecademic Meca500 R3','2 Mecademic Meca500 R3']
RTOOLS = ['1 CrystalToolSample','2 CrystalToolSample']

#setCollisionActivePair()

robots={ 'Names':R_NAMES,'Tools':RTOOLS,'IP': IPs,'Port': PORT}
rdf = pd.DataFrame(robots)


for nrobot in range(0,len(rdf)-1):
    robot=RDK.Item(rdf['Names'][nrobot])
    IP=rdf['IP'][nrobot]
    Port=rdf['Port'][nrobot]
    robot.Connect(IP)
    print(robot)
    #robot.setConnectionParams(IP, Port, 'socket')
    tool=RDK.Item(rdf['Tools'][nrobot])
    robot.setPoseTool(tool) # Update the TCP
    robot.setSpeed(100)  # Set linear speed in mm/s
    robot.setJoints([0,0,0,0,0,0])      # set all robot axes to zero
    

    # Connect to the robot

print('Robots homed')



PR= RDK.Item('Photo Resistor',ITEM_TYPE_TARGET)
X,Y,Z,A,B,C=Pose_2_KUKA(PR.Pose())

# Initialize ESP Webserver 
IP="192.168.0.99"
ext="values"
ESP=Webserver(IP,ext)

#Move Detector to position
print('Attempting gantry move to', -Z)
ESP.patch("goalposLED",-Z)
ESP.get("goalposLED")


for nrobot in range(0,len(rdf)):
    print(nrobot)
    robot=RDK.Item(rdf['Names'][nrobot])

    TPP='TPP'+str(nrobot+1)
    APP='APP'+str(nrobot+1)
    Home='Home'+str(nrobot+1)
    print(TPP)
    target = RDK.Item(TPP)# retrieve the Target item
    Home=RDK.Item(Home)
    tool=RDK.Item(rdf['Tools'][nrobot])
    robot.setPoseTool(tool) # Update the TCP
    robot.setSpeed(100)  # Set linear speed in mm/s
  

    #robot.MoveJ(target)  
    
    rj=robot.Joints()
    Final = target.Pose()
    SafeJM(rj,Final)
    
    print('Robot target approached')
    
    ESP.patch("ledStatus",True)
    for L in range(1,Loops+1):
        if __name__=="__main__":
            print('Starting Spiral '+str(L))
    
            cur_spiral=Spiral(xl,yl,400) #setting params of the spiral we will be using
            revs=20/L**L
            points=200
            x,y=cur_spiral.build(20/L**2,2/L**2) #build our spiral x y coords
            step=.05/L #linear step size
            zstep=2/L
    
            Ix=-10
            Iy=10
            TestInt=Intensity(Ix,Iy,0,0,10,10) #set test intensity parameters     
            prominence=10
            width=1
            height=10
            plateau_size=2 ############################################
            IAS=10
            
            

        #initialize params
        X=[]
        Y=[]
        Z=[]
        Iv=[]
        Cv=[]
    
        for i in range(1,len(x)): #loop to get coord and corresponding intensity
            if len(find_peaks(Iv,prominence=[prominence],width=width,height=height)[0])<=0: #check if we found a peak
                if i==len(x):
                    print('No peaks found in spiral: Spiral centered')
                    coord=(xl,yl,zl)
                    
                    #robot.MoveJ(target.Pose()*transl(coord))
                    
                    
                    rj=robot.Joints()
                    Final = target.Pose()*transl(coord)
                    SafeJM(rj,Final)
                    
                    raise Exception('Alignment Acheived')
                    break
                coord=(x[i],y[i],zl)
                X.append(x[i])
                Y.append(y[i])
                Z.append(0)
                #Ic=PR.get() #get intensity for PR for online run
                Ic=ESP.get("PRInt")
                #Ic=TestInt.get(x[i],y[i]) #test intensity get
                Iv.append(np.array(Ic))
    
                #robot.MoveL(target.Pose()*transl(coord))
                
                rj=robot.Joints()
                Final = target.Pose()*transl(coord)
                SafeLM(rj,Final)
                
                robot.WaitFinished()
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
    
        if len(find_peaks(data['Intensity'],plateau_size=plateau_size,height=height)[0])>0: #check to see if we are at the peak
            PI=find_peaks(data['Intensity'],plateau_size=plateau_size,height=height)[0]
            xl=data['X'][PI]
            yl=data['Y'][PI]
            coord=(xl,yl,zl)
            print(coord)
            #robot.MoveJ(target.Pose()*transl(coord))
            
            rj=robot.Joints()
            Final = target.Pose()*transl(coord)
            SafeJM(rj,Final)
            
            print('SSAP Complete!')
            
            break
    
    
    
        CL=Linear(data,2) #current linear line definition
       
    
        peaks= find_peaks(data['Intensity'],height=height,prominence=[prominence],width=width) 
        IPeaks=peaks[0][0]
        IPeaks=Iv.index(max(Iv))
    
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
        
        print('Linear move '+str(L)+ ': STARTED')
        
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
            
            #robot.MoveL(target.Pose()*transl(coord))
            rj=robot.Joints()
            Final = target.Pose()*transl(coord)
            SafeLM(rj,Final)
    
            if len(find_peaks(-np.array(Il),plateau_size=4,prominence=30)[0])<=0:
                print('Linear moved past alignment')
                break
                
        print('Linear move '+str(L)+ ': COMPLETED')
        
        LLdata={'X':Xl,'Y':Yl,'I':Il}
        Ldata=pd.DataFrame(LLdata)
        
        Lmax_I=Il.index(max(Ldata['I']))
        xl=Xl[Lmax_I]
        yl=Yl[Lmax_I]
    
        fig = plt.figure()
        #plt.style.use('seaborn-notebook')
        ax = fig.add_subplot(111, projection='3d')
        ax.plot3D(data['X'], data['Y'], data['Intensity'])
        ax.plot3D(data['X'][IPeaks], data['Y'][IPeaks], data['Intensity'][IPeaks],"x")
        ax.plot3D(xl, yl, max(Ldata['I']),"x")
        ax.plot3D(Xl, Yl, Il,'o')
    
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        plt.show()
        
        if len(find_peaks(Ldata['I'],plateau_size=plateau_size,height=height)[0])>0: #check to see if we are at the peak
            PI=find_peaks(Ldata['I'],plateau_size=plateau_size,height=height)[0]
            xl=Ldata['X'][PI]
            yl=Ldata['Y'][PI]
            coord=(xl,yl,zl)
            
            #robot.MoveL(target.Pose()*transl(coord))
            rj=robot.Joints()
            Final = target.Pose()*transl(coord)
            SafeLM(rj,Final)
            print('SSAP Complete!')
            break
    
    
        #Push Pull Search
        IPP=[]
        ZLPP=[]
        
        for i in range(3):
            coord=(xl,yl,zl)
            ZLPP.append(zl)
            #robot.MoveJ(target.Pose()*transl(coord))
            
            rj=robot.Joints()
            Final = target.Pose()*transl(coord)
            SafeJM(rj,Final)
            
            PRPPT=[]
            for j in range(IAS):
                #PRPP=PR.get()
                PRPP=ESP.get("PRInt")
                PRPPT.append(PRPP)
            PPI=np.average(PRPPT)
            IPP.append(PPI)
    
            zl=zl+(zstep)
    
        plt.plot(ZLPP,IPP)
        plt.show()
    
        if IPP.index(max(IPP))>0:
            zstep=zstep
        else:
            zstep=-zstep
        ZLPP=[]
        IPP=[]
    
        coord=(xl,yl,zl)
        #robot.MoveJ(target.Pose()*transl(coord))
        
        rj=robot.Joints()
        Final = target.Pose()*transl(coord)
        SafeJM(rj,Final)        
        #maybe savitzky golay or moving average filter
        while len(find_peaks(IPP,prominence=5,height=height)[0])<=0 or PPI <= .9*max(IPP): #Push pull
            coord=(xl,yl,zl)
            ZLPP.append(zl)
            #robot.MoveL(target.Pose()*transl(coord))
            rj=robot.Joints()
            Final = target.Pose()*transl(coord)
            SafeLM(rj,Final)
            
            PRPPT=[]
            for j in range(IAS):
                #PRPP=PR.get()
                PRPP=ESP.get("PRInt")
                PRPPT.append(PRPP)
            PPI=np.average(PRPPT)
            IPP.append(PPI)
            zl=zl+(zstep)
    
    
        zl=ZLPP[IPP.index(max(IPP))]
        
     
        plt.plot(ZLPP,IPP)
        plt.show()
    

    coord=(xl,yl,zl)
    #robot.MoveL(target.Pose()*transl(coord))
    rj=robot.Joints()
    Final = target.Pose()*transl(coord)
    SafeJM(rj,Final)   
    
    ESP.patch("ledStatus",False)
    
    APP_RDK= RDK.Item(APP,ITEM_TYPE_TARGET)
    TPP_RDK= RDK.Item(TPP,ITEM_TYPE_TARGET)
    APP_RDK.setParent(TPP_RDK)
    APP_RDK.setPose(Pose(xl,yl,zl,0,0,0))
    
    rj=robot.Joints()
    Final = Home.Pose()
    SafeJM(rj,Final)
    

    print('***Linear Translation********************************************************')
    print(Ldata)
    print('_________SSAP Report:_________')
    print('Calculated peak is at',xl,yl)
    print('Linear step size is:',step)
