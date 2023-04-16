# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:00:55 2023

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
        

#_____________Robots go beep boop_______________________________________




# Set up the IP address and port number of the robot
IPs = ['192.168.0.100','192.168.0.101']
PORT = [30002,30003]
R_NAMES = ['1 Mecademic Meca500 R3','2 Mecademic Meca500 R3']
RTOOLS = ['1 CrystalToolSample','2 CrystalToolSample']

#setCollisionActivePair()

robots={ 'Names':R_NAMES,'Tools':RTOOLS,'IP': IPs,'Port': PORT}
rdf = pd.DataFrame(robots)


for nrobot in range(1,len(rdf)):
    print(rdf['Names'][nrobot])
    
    robot=RDK.Item(rdf['Names'][nrobot])
    IP=rdf['IP'][nrobot]
    robot.ConnectSafe(IP)
    APP='APP'+str(nrobot+1)

    APP=RDK.Item(APP)
    

    tool=RDK.Item(rdf['Tools'][nrobot])
    robot.setPoseTool(tool) # Update the TCP
    robot.setSpeed(-1,25)  # Set linear speed in mm/s
    robot.setSpeed(25)  # Set linear speed in mm/s
    
    rj=robot.Joints()
    Final = APP.Pose()
    SafeJM(rj,Final)

print('Robots sent to APPs')  