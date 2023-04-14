# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:10:20 2023

@author: westj
"""
import pandas as pd



# Set up the IP address and port number of the robot
IPs = ['192.168.0.100','192.168.0.101']
PORT = [30002,30003]
R_NAMES = ['1 Mecademic Meca500 R3','2 Mecademic Meca500 R3']

robots={ 'Names':R_NAMES,'IP': IPs,'Port': PORT}
rdf = pd.DataFrame(robots)


print(rdf)
print(rdf['Names'][0])

for nrobot in range(0,len(rdf)):
    print(nrobot)
    print(rdf['Names'][nrobot])
    robot=rdf['Names'][nrobot]
    TTP='TTP'+str(nrobot+1)
    print(TTP)
    for i in range(0,len(rdf)):
        if i== nrobot:
            robot.setjoints([0,0,0,0,0,0])
#%%
rj=robot.Joints()
Final = target.Pose()*transl(0,m,0)


rj=robot.Joints()
Final = target.Pose()*transl(0,m,0)
def SafeLM(joints,Final_Pose):
    if robot.MoveL_Test(joints,Final_Pose) == 0:         # linear move to the approach position
        robot.MoveL(Final_Pose)
    else :
        print("collision avoided")
        
LCDetect(robot.Joints,target.Pose()*transl(0,m,0))

#%%


# Connect to RoboDK
RDK = Robolink()

# Get the robot item from RoboDK
robot = RDK.Item('ABB IRB 120', ITEM_TYPE_ROBOT)

# Set up the TCP/IP communication with the robot
robot.setConnectionParams(IP_ADDRESS, PORT, 'socket')
#%%
from robodk.robolink import *       # import the robolink library (bridge with RoboDK)
from robodk.robomath import * 
RDK = Robolink()  

robot'1 Mecademic Meca500 R3'
#%%
s='wow'+str(1)*9
print(s)
s
