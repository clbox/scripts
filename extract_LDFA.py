#! /usr/bin/env python
import numpy as np
import sys
import time
from scipy.interpolate import UnivariateSpline
from sys import argv
import ase
from cube import cube

def num_grad(c, x,y,z):
    """
    calculates finite difference gradient for a given density on a grid
    """
    gradient = np.zeros(3)
    d = 0.1
    dens = c(x,y,z)
    gradient[0] = np.abs((c(x+d,y,z)-c(x-d,y,z))/(2*d))
    gradient[1] = np.abs((c(x,y+d,z)-c(x,y-d,z))/(2*d))
    gradient[2] = np.abs((c(x,y,z+d)-c(x,y,z-d))/(2*d))
    gradient = gradient / dens
    return gradient

def num_grad2(c, x,y,z):
    """
    calculates finite difference gradient for a given density on a grid
    """
    gradient = np.zeros(3)
    d = 0.1
    dens = c(x,y,z)
    gradient[0] = np.abs((c(x+d,y,z)+c(x-d,y,z)-2.*c(x,y,z))/(d*d))
    gradient[1] = np.abs((c(x,y+d,z)+c(x,y-d,z)-2.*c(x,y,z))/(d*d))
    gradient[2] = np.abs((c(x,y,z+d)+c(x,y,z-d)-2.*c(x,y,z))/(d*d))
    gradient = gradient #/ dens**2
    return gradient

class ldfa:

    def __init__(self,cubefile):
        """
        calculates ldfa friction for H on metal
        """
        h_bar = 1.0545718*1e-34 # kg*m**2/s
        m_e= 9.10938215*1e-31 #kg
        pi = 3.14159265359
        self.prefac = h_bar/1.66e-27/(1.0E-10)**2*(pi/3)**(1./3.)*4 

        r_s = np.loadtxt('LDFA.txt',usecols=range(0,1))
        Q = np.loadtxt('LDFA.txt',usecols=range(7,8))
        delta = np.zeros([6,len(Q)])
        fits = []
        for i in range(6):
            delta[i] = np.loadtxt('LDFA.txt',usecols=range(i+1,i+2))
            # print delta[i]
            fits.append(UnivariateSpline(r_s,delta[i],k=3,s=0.025))

        self.fits = fits
        self.delta = delta

        self.c=cube()
        self.c.read(cubefile)

    def get_friction(self,pos):
        """
        returns friction tensor for H2 on Ag(111)
        """
        fits = self.fits
        delta = self.delta
        pi = 3.14159265359

        cell = np.array([
                [5.883128419472, 0.0 ,0.0],
                [-2.9415642097, 5.094938664989 ,0.0],
                [0.0, 0.0 ,37.0],])
        celli = np.linalg.inv(cell)
        gamma_list = []
        for p in [pos[0],pos[1]]:
            f = np.dot(p,celli)
            f[0] = f[0]%1.0
            f[1] = f[1]%1.0
            if f[0]>0.5278:
                f[0] -= 1.0
            if f[1]>0.5278:
                f[1] -= 1.0
            p = np.dot(f,cell)
            gamma = 0.0
            x,y,z = p
            density = self.c(x,y,z+7.210,silent=True)
            if density<1e-7:
                gamma = 0.0
            else:
                radius = (3/(4*pi*density))**(1./3.)
                if radius<10.:
                    for l in range(len(delta)-1):
                        gamma+= self.prefac* \
                                density**(2./3.)*(l+1)* \
                                np.sin(fits[l](radius)-fits[l+1](radius))**2
                else:
                    pass
            # print z, gamma/1e12 
            gamma = gamma/1e12
            gamma_list.append(gamma)

        gamma_list = np.array(gamma_list).repeat(3)
        friction = np.eye(6)*gamma_list

        return friction


    def get_friction_atom(self,pos):
        """
        returns friction tensor for H on Ag(111)
        """
        fits = self.fits
        delta = self.delta
        pi = 3.14159265359

        cell = np.array([
                [5.883128419472, 0.0 ,0.0],
                [-2.9415642097, 5.094938664989 ,0.0],
                [0.0, 0.0 ,37.0],])

                
        celli = np.linalg.inv(cell)
        gamma_list = []
        for p in [pos[0],pos[1]]:
            f = np.dot(p,celli)
            f[0] = f[0]%1.0
            f[1] = f[1]%1.0
            if f[0]>0.5278:
                f[0] -= 1.0
            if f[1]>0.5278:
                f[1] -= 1.0
            p = np.dot(f,cell)
            gamma = 0.0
            x,y,z = p
            density = self.c(x,y,z+7.210,silent=True)
            if density<1e-7:
                gamma = 0.0
            else:
                radius = (3/(4*pi*density))**(1./3.)
                if radius<10.:
                    for l in range(len(delta)-1):
                        gamma+= self.prefac* \
                                density**(2./3.)*(l+1)* \
                                np.sin(fits[l](radius)-fits[l+1](radius))**2
                else:
                    pass
            # print z, gamma/1e12 
            gamma = gamma/1e12
            gamma_list.append(gamma)

        gamma_list = np.array(gamma_list).repeat(3)
        friction = np.eye(6)*gamma_list

        return friction

