astrocloud.py

NAME
    astrocloud - Module for implementing a Monte Carlo simulation of photons propagating through a scattering cloud.

CLASSES
    builtins.object
        Cloud
        Photon

    class Cloud(builtins.object)
     |  Methods defined here:
     |
     |  __init__(self, radius, star_distance, observer_distance, observer_angle, scatter_coefficient, scatter_type='uniform')
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  check_if_escaped(self, photon)
     |      Checks if photon has escaped the cloud radius. Returns true if escaped. False otherwise.
     |      Returns: Boolean.
     |
     |  enter_photon(self)
     |      Initiates a new photon. Chooses a random incident angle from which the initial position is calculated.
     |
     |      Returns: Photon.
     |
     |  hit_target(self, photon)
     |      Checks if photon will be seen by the observer. Returns True if photon is within target area. False otherwise.
     |      Returns: Boolean.
     |
     |  photon_step(self, photon)
     |      Executes a photon scattering step. Samples scattering length L, calls Photon.scatter to step photon in its current
     |      direction, then samples new random scattering direction based on scattering type. Sets new photon direction.     |      Returns: None
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

    class Photon(builtins.object)
     |  Methods defined here:
     |
     |  __init__(self, phi, theta, position)
     |      Initialize self.  See help(type(self)) for accurate signature.
     |
     |  get_dir(self)
     |      Gets Photon direction.
     |      Returns: Float,Float
     |
     |  get_numsteps(self)
     |      Gets number of scattering steps for Photon.
     |      Returns: Int
     |
     |  get_pos(self)
     |      Gets Photon position.
     |      Returns: List [x,y,z]
     |
     |  scatter(self, length)
     |      Executes Photon scattering step. Updates the Photon position and direction.
     |      Returns: None
     |
     |  set_dir(self, new_phi, new_theta)
     |      Sets new Photon direction in spherical coordinates.
     |      Returns: None
     |
     |  set_pos(self, new_position)
     |      Sets new Photon position. new_position should be an array or list in the form [x,y,z]
     |      Returns: None
     |
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |
     |  __dict__
     |      dictionary for instance variables (if defined)
     |
     |  __weakref__
     |      list of weak references to the object (if defined)

FUNCTIONS
    run_scatter(cloud, maxphi, R, Dstar, Dobs, Phobs, asc, num_photons)
        Example scattering simulation. Returns list of simulation data.

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

DATA
    __copyright__ = 'Copyright 2019, Jesse Weller, All rights reserved'
    __email__ = 'wellerj@oregonstate.edu'

VERSION
    1.0

AUTHOR
    Jesse Weller

FILE
    astrocloud.py
