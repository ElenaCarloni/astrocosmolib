import numpy as np
from scipy.integrate import quad 
from scipy.constants import c as speed_of_light
from pyACC.integrate.integration import Integrate

class CosmologicalDistances:
    def __init__(self, hubble_function, *cosmo_pars):
        self.hubble_function = hubble_function
        self.cosmo_pars = cosmo_pars

    def comoving_distance(self, z, *cosmo_pars, **integ_args):
        if len(cosmo_pars) != 0:
            integrand = lambda z: speed_of_light/1.e3 / self.hubble_function(z, *cosmo_pars)
        else:
            integrand = lambda z: speed_of_light/1.e3 / self.hubble_function(z, *self.cosmo_pars)
            
        return Integrate(integrand)("quad", 0.0, z, **integ_args)  

    def angular_distance(self, z, *cosmo_pars, **integ_args):
        return self.comoving_distance(z, *cosmo_pars, **integ_args) / (1 + z)
    
    def luminosity_distance(self, z, *cosmo_pars, **integ_args):
        return (1 + z) * self.comoving_distance(z, *cosmo_pars, **integ_args)

    def hubble_distance(self, z, *cosmo_pars, **integ_args):
        if len(cosmo_pars) !=0:
            return speed_of_light/1.e3 / self.hubble_function(z, *cosmo_pars)
        else:
            return speed_of_light/1.e3 / self.hubble_function(z, *self.cosmo_pars)
        
    def volume_distance(self, z, *cosmo_pars, **integ_args):
        if len(cosmo_pars) !=0:
            return ((speed_of_light / 1.e3 * z) * (self.comoving_distance(z, *cosmo_pars, **integ_args)**2 / self.hubble_function(z, *cosmo_pars)))**(1./3)
        else:
            return ((speed_of_light / 1.e3 * z) * (self.comoving_distance(z, *cosmo_pars, **integ_args)**2 / self.hubble_function(z, *self.cosmo_pars)))**(1./3)
        
    def transverse_comoving_distance(self, z, *cosmo_pars, **integ_args):
        if len(cosmo_pars) != 0:
            return self.comoving_distance(z, *cosmo_pars, **integ_args)
        else:
            return self.comoving_distance(z, *self.cosmo_pars, **integ_args)
    
    def isotropic_volume_distance(self, z, *cosmo_pars, **integ_args):
        if len(cosmo_pars) != 0:
            dm = self.transverse_comoving_distance(z, *cosmo_pars, **integ_args)
            dh = self.hubble_distance(z, *cosmo_pars, **integ_args)
        else:
            dm = self.transverse_comoving_distance(z, *self.cosmo_pars, **integ_args)
            dh = self.hubble_distance(z, *self.cosmo_pars, **integ_args)
        return (z * dm**2 * dh)**(1./3)