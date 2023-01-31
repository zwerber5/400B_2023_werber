#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 17:21:56 2023

@author: zachwerber
"""

import numpy as np
import astropy.units as u
from ReadFile import Read

def Componentmass(filename, particle_type):
    '''
    This function determines the mass of all particles of a certain type in the given file/galaxy

    Inputs
    ----------
    filename : filethe name of the file to be read
        and determine the mass.
    particle_type : the desired particle type to determine the total mass of.

    Output:
    -------
    returns the total mass of all particles of a given type
    in units of 10^12 solar masses
    '''
    time,particle_count,data = Read(filename) # read in file
    
    index = np.where(data['type']==particle_type) #get index of all particles of interest
    newdata=data[index] #make new data table of just the desired particle types
    
    mass=newdata['m'] # make a list of the masses of each particle of the desired type
    totalmass=np.sum(mass)*u.Msun # sum over the the masses and apply units
    totalmass/=1e2 # convert to units of 10^12 Msun
    totalmass=np.around(totalmass,3) #round to 3 decimals
    
    return totalmass


