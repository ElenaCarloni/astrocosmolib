from astropy.io import fits
import numpy as np
from pyACC.fits_wrapper import FitsManager  # Importa FitsManager se necessario

class ComovingDistance:
    def __init__(self, input_file, cosmology):
        "Costruttore del file"
        if isinstance(input_file, FitsManager):
            # Se input_file è un oggetto FitsManager, utilizza la sua hdulist
            self.hdulist = input_file.hdulist
        else:
            # Se è una stringa, apri il file FITS
            self.hdulist = fits.open(input_file)
        self.cosmology = cosmology

    def comoving_coordinates(self, redshift, alpha, delta):
        cd = self.cosmology.comoving_distance(redshift).value
        "Calcola le coordinate comoventi"
        x = cd * np.cos(np.radians(delta)) * np.cos(np.radians(alpha))
        y = cd * np.cos(np.radians(delta)) * np.sin(np.radians(alpha))
        z = cd * np.sin(np.radians(delta))
        return x, y, z
