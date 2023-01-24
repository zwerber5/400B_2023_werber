#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 22:35:40 2023

@author: zachwerber
"""

import numpy as np
import astropy.units as u
from ReadFile import Read

#This function gives the velocity and distance of a specific particle from a given file
#Inputs: filename to be read; the type of particle of interest; the number of the particle of interest
#returns: the distance in kpc, the velocity in km/s, and the mass in M_sun
def ParticleInfo(filename,particle_type,particle_number):
    time,particle_count,data = Read(filename) #get the data from the file as well as its current time and number of particles
    index = np.where(data['type']==particle_type) #get index of all particles of interest
    newdata=data[index] #make new data table of just the desired particle types
    
    
    #note the -1 on each is for indexing since index starts at 0
    desired_x=newdata['x'][particle_number-1] #get the x position of the desired particle
    desired_y=newdata['y'][particle_number-1] #get the y position of the desired particle
    desired_z=newdata['z'][particle_number-1] #get the z position of the desired particle
    
    distance = np.sqrt(desired_x**2 + desired_y**2 + desired_z**2) #get the 3D distance to the particle
    distance = np.around(distance,3)*u.kpc #round to 3 decimals and apply units
    
    desired_vx=newdata['vx'][particle_number-1] #get the x velocity of the desired particle
    desired_vy=newdata['vy'][particle_number-1] #get the y velocity of the desired particle
    desired_vz=newdata['vz'][particle_number-1] #get the z velocity of the desired particle
    
    velo = np.sqrt(desired_vx**2 + desired_vy**2 + desired_vz**2) #get the 3D distance to the particle
    velo = np.around(velo,3)*u.km / u.s #round to 3 decimals and apply units    
    
    mass = newdata['m'][particle_number]*1e10*u.M_sun #get the mass (in 1e10) of the desired particle
    
    return distance,velo,mass #return the desired distance, velocity, and mass

distance,velo,mass = ParticleInfo('MW_000.txt', 2, 100)
