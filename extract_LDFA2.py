#! /usr/bin/env python
import numpy as np
import sys
import time
from scipy.interpolate import UnivariateSpline
from sys import argv
from cube import cube

class ldfa:
    """ Usage: model = ldfa('cubefile')
                tensor = model.get_friction(position_of_atom)
                Where position of atom is a list e.g [1,1,1] in cartesian coordinates / Ang
    """
    def __init__(self,cubefile):
        """
        calculates ldfa friction for H on metal
        """
        h_bar = 1.0545718*1e-34 # kg*m**2/s
        m_e= 9.10938215*1e-31 #kg
        pi = 3.14159265359
        self.prefac = h_bar/1.66e-27/(1.0E-10)**2*(pi/3)**(1./3.)*4 

        self.c=cube()
        self.c.read(cubefile)


    def get_fits(self,phaseshift_filename):
        
        phaseshift_data = np.loadtxt(phaseshift_filename,comments='#')

        r_s = phaseshift_data[:,0]
        Q = phaseshift_data[:,-1]
        delta = phaseshift_data[:,1:-1]
        no_l = np.shape(delta)[1]

        fits = []
        for i in range(no_l):
            fits.append(UnivariateSpline(r_s,delta[:,i],k=3,s=0.025))

        self.fits = fits
        self.delta = delta

    def get_friction_coefficient(self,phaseshift_filename, position_vector):

        self.get_fits(phaseshift_filename)

        fits = self.fits
        delta = self.delta
        r = np.array(position_vector)
        no_l = np.shape(delta)[1]

        friction = 0.0
        density = self.c(r[0],r[1],r[2])
        if density<1e-7:
            friction = 0.0
        else:
            radius = (3/(4*np.pi*density))**(1./3.)
            if radius<10.:
                for l in range(no_l-1):
                    friction+= self.prefac* \
                            density**(2./3.)*(l+1)* \
                            np.sin(fits[l](radius)-fits[l+1](radius))**2
            else:
                pass

        friction = friction/1e12

        return friction

