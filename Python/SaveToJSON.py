# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 20:36:52 2023

@author: westj
"""

import json
	
# python object(dictionary) to be dumped

x=1
y=2
z=3
a=4
b=5
c=6
dict1 ={
	"R1offsets": {
		"x": str(x),"y": str(y),"z": str(z),
        "a": str(a), "b": str(b), "c": str(c),

	},
    
	"R2offsets": {
		"x": str(x),"y": str(y),"z": str(z),
        "a": str(a), "b": str(b), "c": str(c),
	},
}
	
# the json file where the output must be stored
out_file = open("myfile.json", "w")
	
json.dump(dict1, out_file, indent = 6)
	
out_file.close()

#%%

import json
file = open("myfile.json", "r+")
data= json.load(file)
print(data['R2offsets']['x'])
