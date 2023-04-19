# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 11:28:24 2023

@author: westj
"""

class Robot:
    def __init__(self,rdf,index) -> None:
        
        self.Name=self.rdf['RID'][self.index]
        self.IP=self.rdf['IP'][self.index]
        self.Tool=self.rdf['Tool'][self.index]
        self.TPP=self.rdf['TPP'][self.index]
        self.APP=self.rdf['APP'][self.index]
        self.Ex=self.rdf['Ex'][self.index]
        self.Ey=self.rdf['Ey'][self.index]
        self.Ez=self.rdf['Ez'][self.index]

#%%
print(rdf)
robot=Robot(rdf,0)