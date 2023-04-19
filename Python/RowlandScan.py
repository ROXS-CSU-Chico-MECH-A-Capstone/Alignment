# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:18:33 2023

@author: westj
"""

from robodk.robolink import *       # import the robolink library (bridge with RoboDK)
RDK = Robolink()   
from robodk.robomath import *   
import math
import numpy as np
import requests 
import json
import pandas as pd

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
        

  
print(0)
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




for nrobot in range(0,len(rdf)):
    #connect to robot
    print(rdf['Names'][nrobot])
    robot=RDK.Item(rdf['Names'][nrobot])
    IP=rdf['IP'][nrobot]
    robot.ConnectSafe(IP)

    #configure robot params
    tool=RDK.Item(rdf['Tools'][nrobot])

    #robot.setPoseTool(tool) # Update the TCP
    robot.setSpeed(-1,25)  # Set linear speed in mm/s
    robot.setSpeed(25)  # Set linear speed in mm/s

    #call move
    Home=RDK.Item('Home'+str(nrobot+1))
    Final = Home.Pose()
    SafeJM(robot.Joints(),Final)
    
print('Robots homed')  

#set target poses

X,Y,Z,A,B,C=Pose_2_KUKA(PR.Pose())

for nrobot in range(0,len(rdf)):
    print('LEDT'+str(nrobot+1))
    LEDT=RDK.Item('LEDT'+str(nrobot+1))
    TPP=RDK.Item('TPP'+str(nrobot+1))
    
    LEDT.setPose(Pose(0,0,0,0,90,0))
    tx,ty,tz,ta,tb,tc=Pose_2_KUKA(TPP.Pose())
    LEDT.setPose(LEDT.Pose()*rotx(np.deg2rad(-ta))) #set plane for rowland circle to be generated




#set scan parameters
zi=150
zf=400
points=50
PRZ=np.linspace(zi,zf,points)



ESP.patch("ledStatus",True)

print('Scan: Started')
for i in range(0,len(PRZ)):
    print(i)
    PR.setParent(LED)
    Z=-1*PRZ[i] #call coord for particular scan
    PR.setPose(Pose(0,0,Z,0,0,0))

    ESP.patch("speed","100")
    ESP.patch("goalposLED",abs(Z))
    ESP.get("goalposLED")

    theta=np.arccos(Z/(2*1000))
    R=Z/2*np.tan(theta) #distance from emitter detector axis


    #call robots
    for nrobot in range(0,len(rdf)):
        #connect to robot
        print(rdf['Names'][nrobot])
        robot=RDK.Item(rdf['Names'][nrobot])

        #configure robot params
        tool=RDK.Item(rdf['Tools'][nrobot]) #call tool
        robot.setPoseTool(tool) # Update the TCP
        robot.setSpeed(-1,25)  # Set linear speed in mm/s
        robot.setSpeed(25)  # Set linear speed in mm/s
        
        #call move
        TPP=RDK.Item('TPP'+str(nrobot+1)) #find TPP for robot
        APP=RDK.Item('APP'+str(nrobot+1)) #find APP for robot
        ex,ey,ez,ea,eb,ec=Pose_2_KUKA(APP.Pose()) #set error offsets
        
        coord=(-Z/2+ex,0+ey,-R+ez)#rowland circle scan position offset from LEDT
        LEDT=RDK.Item('LEDT'+str(nrobot+1),ITEM_TYPE_TARGET) #find correct LEDT
        Final = LEDT.Pose()*transl(coord)           #create pose
        SafeJM(robot.Joints(),Final)                #call move
        print('Robot Moved')
    
    
    
    
    
    
    
    coord=(-Z/2+ex,0+ey,-R+ez)

    #RCS.setPose(LEDT.Pose()*transl(coord))
    
    #robot.MoveJ(LEDT.Pose()*transl(coord))
    #pause(0.5)
ESP.patch("ledStatus",False)
print('Scan Complete')