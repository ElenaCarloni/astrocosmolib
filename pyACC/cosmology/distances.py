import numpy as np
<<<<<<< Exam
from scipy.integrate import quad 
=======
from scipy.integrate import quad
>>>>>>> main
from scipy.constants import c as speed_of_light
from pyACC.integrate.integration import Integrate

class Distances:
    def __init__(self, hubble_function, h_units=True):
        self.hubble_function = hubble_function
        self.h_units = h_units


    def comoving_distance(self, z, quad_kwargs):
        inverse_E = lambda zz : 1.0 / self.hubble_function.E(zz)
        if isinstance(z, (list, np.ndarray)):
            integral = np.array([quad(inverse_E, 0, zi, **quad_kwargs)[0] for zi in z])
        else:
            integral = quad(inverse_E, 0, z, **quad_kwargs)[0]
        if self.h_units:
            return integral * speed_of_light/1.e3 / 100
        else:
            return integral * speed_of_light/1.e3 / self.hubble_function.H0


    def luminosity_distance(self, z, quad_kwargs):
        return (1 + z) * self.comoving_distance(z, quad_kwargs)
    

    def transverse_comoving_distance(self, z, quad_kwargs):
        dc = self.comoving_distance(z, quad_kwargs)
        if self.hubble_function.Omega_k == 0:
            return dc
        elif self.hubble_function.Omega_k > 0:
            return dc / np.sqrt(np.abs(self.hubble_function.Omega_k)) * np.sinh(np.abs(self.hubble_function.Omega_k))
        else:
            return dc / np.sqrt(np.abs(self.hubble_function.Omega_k)) * np.sin(np.abs(self.hubble_function.Omega_k))
        

    def angular_diameter_distance(self, z, quad_kwargs):
        return self.transverse_comoving_distance(z, quad_kwargs) / (1 + z)


    def hubble_distance(self, z):
        if self.h_units:
            return speed_of_light / 1.e3 / (100 * self.hubble_function.E(z))
        else:
            return speed_of_light / 1.e3 / (self.hubble_function.H(z))

    
    def isotropic_volume_distance(self, z, quad_kwargs):
        dm = self.transverse_comoving_distance(z, quad_kwargs)
        dh = self.hubble_distance(z)
        return ( z  * dm**2 * dh)**(1./3)
    



# This class is used for a specific exercise
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