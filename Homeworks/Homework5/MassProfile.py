#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:10:31 2023

@author: zachwerber
"""

import numpy as np
import astropy.units as u
import astropy.table as tbl
from astropy.constants import G

from ReadFile import Read
from CenterofMass import CenterOfMass
import matplotlib.pyplot as plt
import matplotlib


class MassProfile:
# Class to define Mass Profile for a given galaxy and simulation snapshot

    def __init__(self, galaxy, snap):
        '''
        Class to calculate the mass profile of a given galaxy

        Parameters
        ----------
        galaxy : "str"
            Name of the galaxy for the mass profile to be determined
            
        snap : "int"
            Time in which the mass profile will be determined

        '''
        
        #add a string of the filenumber to the value '000'
        ilbl = '000'+str(snap)
        #remove all but last 3 digits
        ilbl = ilbl[-3:]
        self.filename= "%s_"%(galaxy) + ilbl +'.txt'
        
        #load in data
        self.time, self.total, self.data = Read(self.filename)
        
        self.m = self.data['m']
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
        self.gname = galaxy
    
    def MassEnclosed(self,ptype,rad):
        '''
        This function calculates the mass enclosed of a given radius from the center of mass
        determined by a given particle type

        Parameters
        ----------
        ptype : int
            This number indicates the particle type to be used for calculating the mass enclosed
        rad : array
            array of radii in which the mass enclosed will be calculated from

        Returns
        -------
        enclosed_M: "astropy quantity"
            The mass enclosed for a certain particle type in units of Msun

        '''
        
        #call center of mass class to get COM position
        COM = CenterOfMass(self.filename, ptype)
        COM_p = COM.COM_P(0.1)
        
        index= np.where(self.data['type'] == ptype)[0]
        
        enclosed_M= np.array([]) # create array to store enclosed mass for each radius
        new_m=self.m[index]
        
        #find the radial distance from the center of mass for each particle of the desired type
        COM_rad = np.sqrt((self.x[index] - COM_p[0])**2 +(self.y[index]- COM_p[1])**2+(self.z[index]- COM_p[2])**2)

    
        
        for r in rad:
            new_index = np.where(COM_rad<r) #find where the radial distance of the particle is within the given radius
            summed_m = new_m[new_index] #create an array of masses for the particles that are within the given radius
            enclosed = np.sum(summed_m) #sum all the desired particle masses
            enclosed_M = np.append(enclosed_M,enclosed*1e10)  #append the total particle mass within the radius
        
        return enclosed_M *u.Msun # return array with proper units

    def MassEnclosedTotal(self,rad):
        '''
        Finds the total mass enclosed for all particle types by

        Parameters
        ----------
        rad : array
            array of radii in which the mass enclosed will be calculated from.
            Finding the mass within these radii from the center of mass

        Returns
        -------
        enclosed_by_type : astropy quantity
            This is an array of the total mass of the galaxy within the given radius 

        '''
        enclosed_by_type= np.zeros(len(rad))*u.Msun # create an empty array to add the mass to
        if self.gname!= "M33": # check if you need to ignore bulge mass
            for i in range(1,4):
                mass_enclosed = self.MassEnclosed(i,rad) # get enclosed mass of a particle type for each radius
    
                enclosed_by_type += mass_enclosed #add the mass enclosed for each radius to the array
        else: # ignoring bulge mass
            for i in range(1,3):
                mass_enclosed = self.MassEnclosed(i,rad) # get enclosed mass of a particle type for each radius
    
                enclosed_by_type += mass_enclosed #add the mass enclosed for each radius to the array
        
        return enclosed_by_type
    
    def HernquistMass(self,radius,a,Mhalo):
        '''
        Determine the Mass profile using the Hernquist profile

        Parameters
        ----------
        radius : astropy quantity
            The radius (in kpc) where the enclosed mass will be calculated from
        a : astropy quantity
            The scale factor (units of kpc) for determining the hernquist mass and density profiles
        Mhalo : astropy quantity
            this is the halo mass (in Msun) of the galaxy

        Returns
        -------
        Astropy quantity
            This returns the mass in Msun for the given radius

        '''
        return (Mhalo*(radius**2))/((a+radius)**2) #Hernquist mass returned
    
    def CircularVelocity(self,ptype,rad):
        '''
        This function calculates the circular velocity for a given particle type 
        based on the enclosed mass for given radii

        Parameters
        ----------
        ptype : int
            This number indicates the particle type to be used for calculating the mass enclosed
        rad : array
            array of radii where we will get the mass enclosed from

        Returns
        -------
        astropy quantity
        The circular velocity (in km/s) for a given radius
            

        '''
        new_G=G.to(u.kpc*u.km**2/u.s**2/u.Msun) # convert gravitational constant to proper units
        total_enclosed_m=self.MassEnclosed(ptype,rad) #get the enclosed mass of the given particle type
        velos= np.sqrt(new_G*total_enclosed_m / rad) #standard circular velocity in km/s
        return np.around(velos,2)
    
    
    def CircularVelocityTotal(self,rad):
        '''
        This function calculates the circular velocity for a given radius based on the total enclosed mass

        Parameters
        ----------
        rad : array
            an array of radii that we will calculate the circular velocity for each value

        Returns
        -------
        astropy quantity
            the circular velocity for each given radius in units of km/s

        '''
        total_mass_by_rad = self.MassEnclosedTotal(rad) #get the total mass enclosed
        new_G=G.to(u.kpc*u.km**2/u.s**2/u.Msun) #convert gravitational constant to proper units
        total_velos= np.sqrt(new_G*total_mass_by_rad/rad) #standard circular velocity in km/s
        return np.around(total_velos,2)
        
        
    
    def HernquistVCirc(self,radius,a,Mhalo):
        '''
        This function calculates the circular velocity for a given radius based on 
        the enclosed mass given from the hernquist mass profile

        Parameters
        ----------
        radius : astropy quantity
            The radius (in kpc) where the enclosed mass will be calculated from
        a : astropy quantity
            The scale factor (units of kpc) for determining the hernquist mass and density profiles
        Mhalo : astropy quantity
            this is the halo mass (in Msun) of the galaxy

        Returns
        -------
        astropy quantity
            the circular velocity in km/s based on the hernquist mass profile

        '''
        new_G=G.to(u.kpc*u.km**2/u.s**2/u.Msun) #convert gravitational constant to proper units
        total_enclosed_m=self.HernquistMass(radius, a, Mhalo) #get the total mass enclosed
        velo_hernquist= np.sqrt(new_G*total_enclosed_m / radius) #standard circular velocity in km/s
        return np.around(velo_hernquist,2)


def plotting(r):
    '''
    This function plots the velocity and mass profiles for each galaxy

    Parameters
    ----------
    r : astropy quantity
        an array of radii away from the CoM to compute enclosed mass

    Returns
    -------
    None. saves each plot as a png

    '''
    MW=MassProfile('MW', 0) #initialize for Milky way
    
    #get mass profiles
    DM_mass=MW.MassEnclosed(1, r)
    disk_mass=MW.MassEnclosed(2, r)
    bulge_mass=MW.MassEnclosed(3, r)
    
    #plot mass profiles
    plt.plot(r,DM_mass,c='black',linestyle='--',label='Dark Matter')
    plt.plot(r,disk_mass,c='red',linestyle=':',label='Disk Mass')
    plt.plot(r,bulge_mass,c='blue',linestyle='-',label='Bulge Mass')
    
    #get and plot total mass on same graph
    total_mass=MW.MassEnclosedTotal(r)
    plt.plot(r,total_mass,c='green',linestyle='dashdot',label='Total Mass')
    
    #alter graph layout
    plt.semilogy()
    plt.xlabel('Radius (kpc)', fontsize=12)
    plt.ylabel(r'Mass Enclosed ($M_\odot$)', fontsize=12)
    plt.title("Mass Enclosed by Particle Type for the Milky Way")
    
    
    #get hernquist mass profile
    hern = np.zeros(len(DM_mass))*u.Msun #initialize the array
    index= np.where(MW.data['type']==1) #find all dark matter particles
    for i in range(len(DM_mass)): #loop through the radii
        hern[i]+= MW.HernquistMass(r[i],60*u.kpc, np.sum(MW.m[index])*1e10*u.Msun) #call hernquist mass with a=60
    plt.plot(r,hern,color='orange',label='Hernquist')
    
    plt.legend()
    plt.savefig('MW_mass')
    
    
    plt.clf() #clear figure to graph new data after saving
    
    #get circular velocity profiles
    DM_velo=MW.CircularVelocity(1, r)
    disk_velo=MW.CircularVelocity(2, r)
    bulge_velo=MW.CircularVelocity(3, r)
    
    #graph circular velocity profiles
    plt.plot(r,DM_velo,c='black',linestyle='--',label='Dark Matter')
    plt.plot(r,disk_velo,c='red',linestyle=':',label='Disk')
    plt.plot(r,bulge_velo,c='blue',linestyle='-',label='Bulge')
    
    #get and graph circular velocity for the total mass
    total_velo=MW.CircularVelocityTotal(r)
    plt.plot(r,total_velo,c='green',linestyle='dashdot',label='Total Mass')
    
    #plt.semilogy()
    
    
    plt.xlabel('Radius (kpc)', fontsize=12)
    plt.ylabel(r'Velocity (km/s)', fontsize=12)
    plt.title('Velocity for distances from the Milky Way center')
    
    #get velocity with hernquist mass profile
    hern = np.zeros(len(DM_velo))*u.km / u.s
    index= np.where(MW.data['type']==1)
    for i in range(len(DM_velo)):
        hern[i]+= MW.HernquistVCirc(r[i],60*u.kpc, np.sum(MW.m[index])*1e10*u.Msun)
    plt.plot(r,hern,color='orange',label='Hernquist')
    
    plt.legend()
    plt.savefig("MW_velocity")
    
    plt.clf()
    
    
    
    
    
    
    
    
    M33=MassProfile('M33', 0)#initialize for M33
    
    #get mass profiles
    DM_mass=M33.MassEnclosed(1, r)
    disk_mass=M33.MassEnclosed(2, r)
    #bulge_mass=M33.MassEnclosed(3, r)
    
    #graph mass profiles
    plt.plot(r,DM_mass,c='black',linestyle='--',label='Dark Matter')
    plt.plot(r,disk_mass,c='red',linestyle=':',label='Disk Mass')
    
    #get total mass profile
    total_mass=M33.MassEnclosedTotal(r)
    plt.plot(r,total_mass,c='green',linestyle='dashdot',label='Total Mass')
    
    plt.semilogy()
    
    
    plt.xlabel('Radius (kpc)', fontsize=12)
    plt.ylabel(r'Mass Enclosed ($M_\odot$)', fontsize=12)
    plt.title("Mass Enclosed by Particle Type for the M33")
    
    
    #get hernquist mass profile
    hern = np.zeros(len(DM_mass))*u.Msun
    index= np.where(M33.data['type']==1) # find all dark matter particles
    for i in range(len(DM_mass)):
        hern[i]+= M33.HernquistMass(r[i],25*u.kpc, np.sum(M33.m[index])*1e10*u.Msun) # find hernquist mass with a=25
    plt.plot(r,hern,color='orange',label='Hernquist')
    
    plt.legend()
    
    plt.savefig('M33_mass')
    
    plt.clf()
    
    #get circular velocity profiles
    DM_velo=M33.CircularVelocity(1, r)
    disk_velo=M33.CircularVelocity(2, r)
    #bulge_velo=M33.CircularVelocity(3, r)
    
    #graph circular velocity profiles
    plt.plot(r,DM_velo,c='black',linestyle='--',label='Dark Matter')
    plt.plot(r,disk_velo,c='red',linestyle=':',label='Disk')
    #plt.plot(r,bulge_velo,c='blue',linestyle='-',label='Bulge Mass')
    
    #find and graph total circular velocity based on the total mass
    total_velo=M33.CircularVelocityTotal(r)
    plt.plot(r,total_velo,c='green',linestyle='dashdot',label='Total Mass')
    
    #plt.semilogy()
    
    
    plt.xlabel('Radius (kpc)', fontsize=12)
    plt.ylabel(r'Velocity (km/s)', fontsize=12)
    plt.title('Velocity for distances from the M33 center')
    
    #find the circular velocity based on the hernquist mass profile
    hern = np.zeros(len(DM_velo))*u.km / u.s
    index= np.where(M33.data['type']==1)
    for i in range(len(DM_velo)):
        hern[i]+= M33.HernquistVCirc(r[i],25*u.kpc, np.sum(M33.m[index])*1e10*u.Msun)
    plt.plot(r,hern,color='orange',label='Hernquist')
    
    plt.legend()
    plt.savefig('M33_velo')
    
    plt.clf()
    
    
    
    
    
    
    
    
    
    M31=MassProfile('M31', 0) #initialize for M31
    
    #find the mass profiles
    DM_mass=M31.MassEnclosed(1, r)
    disk_mass=M31.MassEnclosed(2, r)
    bulge_mass=M31.MassEnclosed(3, r)
    
    #plot the mass profile
    plt.plot(r,DM_mass,c='black',linestyle='--',label='Dark Matter')
    plt.plot(r,disk_mass,c='red',linestyle=':',label='Disk Mass')
    plt.plot(r,bulge_mass,c='blue',linestyle='-',label='Bulge Mass')
    
    #find and plot the total mass profile
    total_mass=M31.MassEnclosedTotal(r)
    plt.plot(r,total_mass,c='green',linestyle='dashdot',label='Total Mass')
    
    plt.semilogy()
    
    
    plt.xlabel('Radius (kpc)', fontsize=12)
    plt.ylabel(r'Mass Enclosed ($M_\odot$)', fontsize=12)
    plt.title("Mass Enclosed by Particle Type for the M31")
    
    #find the hernquist mass
    hern = np.zeros(len(DM_mass))*u.Msun
    index= np.where(M31.data['type']==1)
    for i in range(len(DM_mass)):
        hern[i]+= M31.HernquistMass(r[i],60*u.kpc, np.sum(M31.m[index])*1e10*u.Msun)
    plt.plot(r,hern,color='orange',label='Hernquist')
    
    plt.legend()
    plt.savefig('M31_mass')
    
    plt.clf()
    
    
    
    #find the circular velocity profiles
    DM_velo=M31.CircularVelocity(1, r)
    disk_velo=M31.CircularVelocity(2, r)
    bulge_velo=M31.CircularVelocity(3, r)
    
    #graph the circular velocity profiles
    plt.plot(r,DM_velo,c='black',linestyle='--',label='Dark Matter')
    plt.plot(r,disk_velo,c='red',linestyle=':',label='Disk')
    plt.plot(r,bulge_velo,c='blue',linestyle='-',label='Bulge')
    
    total_velo=M31.CircularVelocityTotal(r)
    plt.plot(r,total_velo,c='green',linestyle='dashdot',label='Total Mass')
    
    #plt.semilogy()
    
    
    plt.xlabel('Radius (kpc)', fontsize=12)
    plt.ylabel(r'Velocity (km/s)', fontsize=12)
    plt.title('Velocity for distances from the M31 center')
    
    #find the circular velocity using the hernquist mass profile
    hern = np.zeros(len(DM_velo))*u.km / u.s
    index= np.where(M31.data['type']==1)
    for i in range(len(DM_velo)):
        hern[i]+= M31.HernquistVCirc(r[i],60*u.kpc, np.sum(M31.m[index])*1e10*u.Msun)
    plt.plot(r,hern,color='orange',label='Hernquist')
    
    plt.legend()
    plt.savefig('M31_velo')
    
    plt.clf()


if __name__ == "__main__":
    r = np.arange(0.25,30.5,1.5)*u.kpc
    
    plotting(r)