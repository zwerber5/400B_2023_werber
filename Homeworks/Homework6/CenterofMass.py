# Homework 4
# Center of Mass Position and Velocity
# Solutions: G.Besla, R. Li, H. Foote


# remember this is just a template, you don't need to follow every step
# if you have your own method to solve the homework, it is totally fine



# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read




class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        # write your own code to compute the generic COM 
        #using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        a_com = np.sum(a*m)/np.sum(m) #np arrays multiplied together just multiply same index. Then sum over the product
        # ycomponent Center of mass
        b_com = np.sum(b*m)/np.sum(m)
        # zcomponent Center of mass
        c_com = np.sum(c*m)/np.sum(m)
        
        # return the 3 components separately
        return a_com, b_com, c_com
       
    
    
    def COM_P(self, delta):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.
        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2) # standard magnitude of 3D


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2) # standard magnitude of 3D

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(r_new < r_max)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)
            # compute the new 3D COM position
            # write your own code below
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2) # standard magnitude of 3D

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            # uncomment the following line if you want to check this                                                                                               
            #print ("CHANGE = ", change)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            r_max /= 2.0
            # check this.                                                                                              
            #print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            x_new = self.x - x_COM2
            y_new = self.y - y_COM2
            z_new = self.z - z_COM2
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2) # standard magnitude of 3D

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position                                                                                                                                                       
            #p_COM = np.array[x_COM, y_COM, z_COM]

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        x_COM=np.around(x_COM,2)
        y_COM=np.around(y_COM,2)
        z_COM=np.around(z_COM,2)
        r_COM = np.array([x_COM,y_COM,z_COM])*u.kpc
    
        return r_COM
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.
        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        # write your own code below
        xV = self.x*u.kpc-x_COM #x distance from center of mass
        yV = self.y*u.kpc-y_COM #y distance from center of mass
        zV = self.z*u.kpc-z_COM #z distance from center of mass
        rV = np.sqrt(xV**2+yV**2+zV**2) #actual particle distance from center of mass
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(rV < rv_max)
        
        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new =  self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        
        # create an array to store the COM velocity
        # write your own code below
        v_COM = np.array([vx_COM,vy_COM,vz_COM])
        

        # return the COM vector
        # set the correct units using astropy
        # round all values                                                                                        
        
        v_COM= np.around(v_COM,2)*u.km/u.s
        return v_COM
    

# ANSWERING QUESTIONS
#######################
#get MW COM data
# =============================================================================
# if __name__ == '__main__' : 
# 
#     # Create a Center of mass object for the MW, M31 and M33
#     # below is an example of using the class for MW
#     MW_COM = CenterOfMass("MW_000.txt", 2)
# 
#     MW_COM_p = MW_COM.COM_P(0.1)
#     print('COM position of Milky Way:',MW_COM_p)
#     MW_COM_v = MW_COM.COM_V(MW_COM_p[0],  MW_COM_p[1],  MW_COM_p[2])
#     print('COM velocity of Milky Way:',MW_COM_v)
# 
# #get M31 COM data
#     M31_COM = CenterOfMass("M31_000.txt", 2)
#     M31_COM_p = M31_COM.COM_P(0.1)
#     M31_COM_v = M31_COM.COM_V(M31_COM_p[0],  M31_COM_p[1],  M31_COM_p[2])
#     print('COM position of M31:',M31_COM_p)
#     print('COM velocity of M31:',M31_COM_v)
#     
# #get M33 data
#     M33_COM = CenterOfMass("M33_000.txt", 2)
#     M33_COM_p = M33_COM.COM_P(0.1)
#     M33_COM_v = M33_COM.COM_V(M33_COM_p[0],  M33_COM_p[1],  M33_COM_p[2])
#     print('COM position of M33:',M33_COM_p)
#     print('COM velocity of M33:',M33_COM_v)
# 
# 
# 
# MWto31_r=MW_COM_p - M31_COM_p #get difference in COM xyz
# MWto31_r=np.around(np.sqrt(np.sum(MWto31_r**2)),3) #sum the squares of xyz and then square root and round to 3
# print('Current magnitude of separation of COM of MW and M31:',MWto31_r) 
# 
# MWto31_v=MW_COM_v - M31_COM_v #get difference in COM vx,vy,vz
# MWto31_v=np.around(np.sqrt(np.sum(MWto31_v**2)),3) #sum the squares of vx,vy,vz and then square root and round to 3
# print('Current velocity magnitude of COM of MW and M31:',MWto31_v)
# 
# 
# M33to31_r=M33_COM_p - M31_COM_p #get difference in COM xyz
# M33to31_r=np.around(np.sqrt(np.sum(M33to31_r**2)),3) #sum the squares of xyz and then square root and round to 3
# print('Current magnitude of separation of COM of M31 and M33:',M33to31_r)
# 
# M33to31_v=M33_COM_v - M31_COM_v #get difference in COM vx,vy,vz
# M33to31_v=np.around(np.sqrt(np.sum(M33to31_v**2)),3) #sum the squares of vx,vy,vz and then square root and round to 3
# print('Current velocity magnitude of COM of M31 and M33:',M33to31_v)
# 
# 
# print('q4:')
# print('An iterative process is important since the galaxies could be colliding already. \n'
#       'In this case, the particles become difficult to distinguish who they belong to and \n'
#       'iterating to smaller bounds allows us to capture the proper particles to the proper \n'
#       ' galaxy to still track the changing of the center of mass position and velocity.')
# =============================================================================

    