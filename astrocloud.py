#!/usr/bin/env python
r"""Package for implementing a Monte Carlo simulation of photons propagating through a scattering cloud.

"""

__author__ = "Jesse Weller"
__copyright__ = "Copyright 2019, Jesse Weller, All rights reserved"
__version__ = "1.0"
__email__ = "wellerj@oregonstate.edu"



import numpy as np
import random
import math
import matplotlib.pyplot as plt

class Photon:

    def __init__(self, phi, theta, position):
        self.phi = phi
        self.theta = theta
        self.position = position
        self.stepnum = 0

    def get_dir(self):
        """Gets Photon direction.
            Returns: Float,Float """
        return self.phi,self.theta

    def get_pos(self):
        """Gets Photon position.
            Returns: List [x,y,z]"""
        return self.position

    def get_numsteps(self):
        """Gets number of scattering steps for Photon.
            Returns: Int"""
        return self.stepnum

    def set_dir(self, new_phi,new_theta):
        """Sets new Photon direction in spherical coordinates.
            Returns: None"""
        self.phi = new_phi
        self.theta = new_theta

    def set_pos(self, new_position):
        """Sets new Photon position. new_position should be an array or list in the form [x,y,z]
        Returns: None"""
        self.position = new_position
        
    def scatter(self, length):
        """Executes Photon scattering step. Updates the Photon position and direction.
            Returns: None"""
        dx = math.cos(self.theta)*math.sin(self.phi)
        dy = math.sin(self.theta)*math.sin(self.phi)
        dz = math.cos(self.phi)
        if self.stepnum==0:
            if dz<0:
                print(dx,dy,dz)
                print('wrong direction!')        
        step = length*np.asarray([dx,dy,dz])
        self.position += step
        self.stepnum += 1



class Cloud:

    def __init__(self,
                radius,
                star_distance,
                observer_distance,
                observer_angle,
                scatter_coefficient,
                scatter_type = 'uniform'):
        
        self.r = radius # cloud radius
        self.dstar = star_distance # star distance from cloud
        self.dobs = observer_distance # observer distance from cloud
        self.angobs = observer_angle  # angle within which photons are considered to be seen by observer.
        self.asc = scatter_coefficient # cloud scattering coefficient

        self.scatter_type = scatter_type # scattering type key DEFAULT: 'uniform'


    def enter_photon(self):
        """Initiates a new photon. Chooses a random incident angle from which the initial position is calculated.

            Returns: Photon.
            """
        maxphi_star = math.atan2(self.r,self.dstar)

        # initial direction of travel
        phi_star = random.random()*maxphi_star
        theta_star = random.random()*maxphi_star
        if phi_star>(math.pi/2):
            print('wrong_direction!!!')

        phi = math.asin(self.dstar*math.cos(phi_star)*math.sin(phi_star)/self.r)
        theta = math.asin(self.dstar*math.cos(theta_star)*math.sin(theta_star)/self.r)

        # initial position
        x = -math.cos(theta)*math.sin(phi)
        y = -math.sin(theta)*math.sin(phi)
        z = -math.cos(phi)
        position = self.r*np.asarray([x,y,z])

        return Photon(phi_star,theta_star,position)

    def check_if_escaped(self, photon):
        """Checks if photon has escaped the cloud radius. Returns true if escaped. False otherwise.
             Returns: Boolean."""
        if np.sqrt(np.sum(np.square(photon.get_pos()))) > (self.r*1.0000000001):
            return True
        return False


    def photon_step(self, photon):
        """Executes a photon scattering step. Samples scattering length L, calls Photon.scatter to step photon in its current
            direction, then samples new random scattering direction based on scattering type. Sets new photon direction.
            Returns: None
            """
        rand = random.random()
        tau = -math.log(rand)   # randomly choose optical depth traveled
        L = tau/self.asc        # photon travels distance L before scattered
        
        photon.scatter(min(L,2*self.r))       # move photon another scatter step

        if not self.check_if_escaped(photon):            
            # choose new phi and theta
            if self.scatter_type == 'uniform':
                new_phi = random.random()*math.pi
                new_theta = random.random()*2*math.pi
            # if self.scatter_type == 'other':
                # Define 'other' sampling procedure here.
                        
            photon.set_dir(new_phi,new_theta)
        
        

    def hit_target(self, photon):
        """Checks if photon will be seen by the observer. Returns True if photon is within target area. False otherwise.
             Returns: Boolean."""
        target = np.asarray([0,0,self.dobs])
        pos = photon.get_pos()

        dpos = target-pos
        dpos = dpos/np.linalg.norm(dpos) # unit vector from exit pointing toward center of target
        x,y,z = dpos[0],dpos[1],dpos[2]
        phi_target = math.atan2(math.sqrt(x**2+y**2),z)
        theta_target = math.atan2(y,x)

        phi_target = (phi_target+2*math.pi) % (2*math.pi)
        if phi_target>(math.pi/2) or phi_target<0:
            print('wrong target phi!')
        theta_target = (theta_target+2*math.pi) % (2*math.pi)
        if theta_target > (2*math.pi):
            print('wrong target theta!')
        
        phi_photon = photon.get_dir()[0]
        theta_photon = photon.get_dir()[1]



        dphi = abs(phi_photon - phi_target)
        dth = abs(theta_photon - phi_target)

        if max(dphi,dth) < self.angobs:
            return True, max(dphi,dth)

        return False, max(dphi,dth)




def run_scatter(cloud, maxphi,R,Dstar,Dobs,Phobs,asc,num_photons):
    """Example scattering simulation. Returns list of simulation data.

    Returns: List.
            [maxphi: Maximum photon incident angle
            R: Cloud radius. 
            Dstar: Star distance from center of cloud.
            Dobs: Observer distance from center of cloud.
            Phobs: Angle within which photons are considered to be seen by observer.
            asc: Cloud scattering coefficient.
            hits: Number of photons seen by observer.
            stuck: Number of photons that don't escape cloud within 10**7 steps
            reflected: Number of photons that leave cloud on same side as they enter (z-position is negative)
            mean_phsteps: Mean photon scattering steps.
            phpos: List of each photon's exit position.
            angle_wrt_observer: List of angles each photon radiates with respect to observer]
    
    """
    phsteps = [] # list to record number of scattering steps
    phpos = [] # list to record each photon's exit position
    angle_wrt_observer = [] # list to record angle photon radiates wrt observer
    hits = 0 # number of photons seen by observer
    out = 0 # number of photons with initial positions outside of cloud (for testing purposes)
    stuck = 0 # number of photons that exceed 10**7 scattering steps within cloud (for testing purposes)
    reflected = 0 # number of photons that leave cloud on same side as they enter (z-position is negative)
    for photon in range(num_photons):
        if photon % max(10,(num_photons/100)) == 0:
            print('\nReflected: ', reflected, 'of', photon)
            print('Hit Target: ', hits, 'of', photon)
        ph = cloud.enter_photon()

        if cloud.check_if_escaped(ph):
            print('landed out')
            out+=1
            
        while not cloud.check_if_escaped(ph):
            if ph.get_numsteps() > 10**7:
                print('photon is stuck')
                print('position: ', ph.get_pos())
                print('dir: ', ph.get_dir())
                stuck += 1
                break
            cloud.photon_step(ph)
        hit_bool, angle = cloud.hit_target(ph)
        angle_wrt_observer.append(angle)
        phsteps.append(ph.get_numsteps())   # record number of steps
##        phpos.append(ph.get_pos())        # record photon exit position (produces a lot of data)
        if hit_bool:
            hits +=1
        if ph.get_pos()[2] < 0:
            reflected += 1
    
    print('\nLanded Out: ', out, 'of', num_photons)
    print('Stuck: ', stuck, 'of', num_photons)
    print('Reflected: ', reflected, 'of', num_photons)
    print('Hit Target: ', hits, 'of', num_photons)
    
    mean_phsteps = np.mean(phsteps) # calculate mean photon steps
    angle_wrt_observer = np.round((np.asarray(angle_wrt_observer)*180/math.pi*60)).astype(np.int) # convert to arcminutes
    return [maxphi,R,Dstar,Dobs,Phobs,asc,hits,stuck,reflected,mean_phsteps,phpos,angle_wrt_observer]


if __name__ == '__main__':
    """Example main function for running astrocloud simulation."""
    
    maxphi = np.pi/180*10 # cross section of cloud from star viewpoint

    R = 10**6 # cm
    Dstar = R/math.tan(maxphi)    # star_distance
    Dobs = Dstar    # observer_distance
    Phobs = np.pi/180/60*1   # observer_angle
    asc = 10**-6
    scatter_type = 'uniform'

    cloud = Cloud(R, Dstar, Dobs, Phobs, asc, scatter_type)

    num_photons = 10**5
    phsteps = []
    hits = 0 
    out = 0
    stuck = 0
    reflected = 0
    for photon in range(num_photons):
        if photon % max(10,(num_photons/100)) == 0:
            print('..'+str(photon),end = '')
        ph = cloud.enter_photon()
        
        if cloud.check_if_escaped(ph):
            print('landed out')
            out+=1
            
        while not cloud.check_if_escaped(ph):
            if ph.get_numsteps() > 10**7:
                print('photon is stuck')
                print('position: ', ph.get_pos())
                print('dir: ', ph.get_dir())
                stuck += 1
                break
            cloud.photon_step(ph)
        phsteps.append(ph.get_numsteps())
        if cloud.hit_target(ph):
            hits +=1
        if ph.get_pos()[2] < 0:
            reflected += 1

    # print results of simulation
    print('\nLanded Out: ', out, 'of', num_photons)
    print('Stuck in cloud: ', stuck, 'of', num_photons)
    print('Reflected: ', reflected, 'of', num_photons)
    print('Seen by observer: ', hits, 'of', num_photons)


    # create histogram of photon steps before escaping cloud.
    plt.figure()
    bins = np.arange(0,100,1)
    plt.hist(phsteps, bins=bins, align='left')
    plt.xlabel('Photon Steps')

    plt.show()


