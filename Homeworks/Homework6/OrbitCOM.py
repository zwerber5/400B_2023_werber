

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
import sys

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterofMass_voldec import CenterOfMass



def OrbitCOM(galaxy,start,end,n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy: string
            name of the galaxy to loop over
        start: int
            the number of the first snapshot
        end: int
            the number of the last snapshot
        n: int
            intervals over which the COM will be returned
          
    outputs: 
        COM_p: astropy quantity
            the center of mass position at a certain time
        COM_v: astropy quantity
            the velocity vector of the center of mass position
        time: int
            the time that passed
    """
    
    # compose the filename for output
    fileout= 'Orbit_'+galaxy+'_'+str(start)+'_'+str(end)+'.txt'
    
    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    delta=.1
    voldec=2
    # for M33 that is stripped more, use different values for VolDec
    M33_delta=.1
    M33_voldec=4

    
    # generate the snapshot id sequence 
    snap_ids = np.arange(start,end+1,n)
    # it is always a good idea to also check if the input is eligible (not required)
    arraysize = np.size(snap_ids)
    if arraysize==0:
        print('This array is empty and thus not valid. Provide better values you doughnut')
        sys.exit()
    
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros([len(snap_ids),7])
    
    
    # a for loop 
    for  i,snap_id in enumerate(snap_ids): # loop over files
        # compose the data filename (be careful about the folder)
        ilbl = '000'+str(snap_id)
        #remove all but last 3 digits
        ilbl = ilbl[-3:]
        filename= './lowrestexts/'+galaxy+'_' + str(ilbl) +'.txt'
        
        
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        if galaxy!= 'M33':
            COM_p = COM.COM_P(delta,voldec)
            COM_v= COM.COM_V(COM_p[0], COM_p[1], COM_p[2])
        else:
            COM_p = COM.COM_P(M33_delta,M33_voldec)
            COM_v= COM.COM_V(COM_p[0], COM_p[1], COM_p[2])
            
       
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        time = COM.time / 1000
        orbit[i][0] =time.value
        orbit[i][1] =COM_p[0].value
        orbit[i][2] =COM_p[1].value
        orbit[i][3] =COM_p[2].value
        orbit[i][4] =COM_v[0].value
        orbit[i][5] =COM_v[1].value
        orbit[i][6] =COM_v[2].value
        # note that you can store 
        # a[i] = var1, *tuple(array1)

        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
#OrbitCOM('MW', 0, 800, 5)
#OrbitCOM('M31', 0, 800, 5)
#OrbitCOM('M33', 0, 800, 5)
# Note: This might take a little while - test your code with a smaller number of snapshots first! 




# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
data_MW = np.genfromtxt("Orbit_MW_0_800.txt",dtype=None,names=True,skip_header=0)
data_M31=np.genfromtxt("Orbit_M31_0_800.txt",dtype=None,names=True,skip_header=0)
data_M33=np.genfromtxt("Orbit_M33_0_800.txt",dtype=None,names=True,skip_header=0)




# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def VectorDifference(vector1,vector2):
    '''
    This is a generic function to compute the difference between 2 vectors and then find the magnitude of
    the difference vector i.e., the relative vector

    Parameters
    ----------
    vector1 : array
        array with the coordinates of the first vector
    vector2 : array
        array with the coordinates of the second vector

    Returns
    -------
    magnitude_difference: float
        the magnitude of the difference of the 2 vectors that were passed in

    '''
    dif = vector1-vector2 # difference of 2 vectors
    magnitude_difference = 0 # initialize the magnitude
    for x in dif: # loop through difference vector
        magnitude_difference+= x**2 # sum the squares
    magnitude_difference= np.sqrt(magnitude_difference) # square root of the sum of squars
    return magnitude_difference


# Determine the magnitude of the relative position and velocities 

MW_positions= np.zeros([len(data_MW['x']),3])
for i in range(0,len(data_MW['x'])):
    MW_positions[i] = [data_MW['x'][i],data_MW['y'][i],data_MW['z'][i]]

M31_positions= np.zeros([len(data_M31['x']),3])
for i in range(0,len(data_M31['x'])):
    M31_positions[i] = [data_M31['x'][i],data_M31['y'][i],data_M31['z'][i]]


M33_positions= np.zeros([len(data_M33['x']),3])
for i in range(0,len(data_M33['x'])):
    M33_positions[i] = [data_M33['x'][i],data_M33['y'][i],data_M33['z'][i]]





MW_velos= np.zeros([len(data_MW['vx']),3])
for i in range(0,len(data_MW['vx'])):
    MW_velos[i] = [data_MW['vx'][i],data_MW['vy'][i],data_MW['vz'][i]]

M31_velos= np.zeros([len(data_M31['vx']),3])
for i in range(0,len(data_M31['vx'])):
    M31_velos[i] = [data_M31['vx'][i],data_M31['vy'][i],data_M31['vz'][i]]


M33_velos= np.zeros([len(data_M33['vx']),3])
for i in range(0,len(data_M33['vx'])):
    M33_velos[i] = [data_M33['vx'][i],data_M33['vy'][i],data_M33['vz'][i]]

# of MW and M31
difference_positions_MW_M31 = np.array([])
difference_velos_MW_M31 = np.array([])
for i in range(0,len(data_MW['x'])):
    difference_positions_MW_M31 = np.append(difference_positions_MW_M31,VectorDifference(MW_positions[i], M31_positions[i]))
    difference_velos_MW_M31 = np.append(difference_velos_MW_M31,VectorDifference(MW_velos[i], M31_velos[i]))




# of M33 and M31
difference_positions_M33_M31 = np.array([]) #initialize arrays
difference_velos_M33_M31 = np.array([])
for i in range(0,len(data_M33['x'])):
    difference_positions_M33_M31 = np.append(difference_positions_M33_M31,VectorDifference(M33_positions[i], M31_positions[i]))
    difference_velos_M33_M31 = np.append(difference_velos_M33_M31,VectorDifference(M33_velos[i], M31_velos[i]))



# =============================================================================
# # Plot the Orbit of the galaxies 
# #################################
# plt.plot(data_MW['t'],difference_positions_M33_M31, color = 'red', linewidth =2, label = 'M33-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Distance (Kpc)')
# plt.title('Orbit of M31-M33')
# plt.legend()
# plt.savefig('Orbit_M31_M33')
# 
# plt.plot(data_MW['t'],difference_positions_MW_M31, color = 'red', linewidth =2, label = 'MW-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Distance (Kpc)')
# plt.title('Orbit of MW-M31')
# plt.legend()
# plt.savefig('Orbit_MW_M31')
# 
# 
# 
# 
# # Plot the orbital velocities of the galaxies 
# #################################
# plt.plot(data_MW['t'],difference_velos_M33_M31, color = 'red', linewidth =2, label = 'M33-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Velocity (km/s)')
# plt.title('Relative Velocity of M31-M33')
# plt.legend()
# plt.savefig('Velo_M31_M33')
# 
# plt.plot(data_MW['t'],difference_velos_MW_M31, color = 'red', linewidth =2, label = 'MW-M31')
# plt.xlabel('time (Gyr)')
# plt.ylabel('Velocity (km/s)')
# plt.title('Relative Velocity of MW-M31')
# plt.legend()
# plt.savefig('Velo_MW_M31')
# =============================================================================

'''
Questions

1) There will be 2 close encounters between MW and M33 before they merger in the third encounter
2) As the galaxies get closer, the velocity increases. Essentially, velocity appears inversely related to
    separation.
3) MW and M31 fully merge just after 6 Gyr from now
bonus) M33 appears to be decaying by ~ 10 kpc / Gyr just after 6 Gyr. If this continued, it would take ~7.5 Gyrs for M33 to 
    then merge with the MW-M31 remnant.
'''