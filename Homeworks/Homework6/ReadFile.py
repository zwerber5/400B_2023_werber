#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 17:10:43 2023

@author: zachwerber
"""

import numpy as np
import astropy.units as u


#This function reads in datafiles in the 400B format
#Input: the name of the file to be read
#returns the time (in Myr), number of particles, and descriptions of each particle (type, position, and velocity)
def Read(filename):
    file = open(filename,'r') #opens file
    line1= file.readline() #load the first line
    label,value = line1.split() #break line into type and actual value
    time = float(value)*u.Myr #apply units and save the current time of the file
    
    line2 = file.readline() #load in the second line
    label,value = line1.split() #break line into type and actual value
    particle_count = float(value) #Save the total number of particles
    
    file.close() #close the file
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3) #load in the actual data from the file
    #labels: type, m,x,y,z,vx,vy,vz
    
    return time, particle_count, data


