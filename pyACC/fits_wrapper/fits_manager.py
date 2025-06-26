from astropy.io import fits
import numpy as np

from ..helpers.logger import Logger

class FitsManager:
    "class to merge fits_file using astropy.io.fits module"
    #define the constructure:
    def __init__(self, input_file):
        "Constructure of the file"
        self.input_file = input_file
        self.hdulist = fits.open(input_file)

    def get_hdu_count(self):
        "get the number of HDU in the fits file"
        return len(self.hdulist)

    def get_header(self, hdu_index):
        "get the header of the fits file"

        if hdu_index < 0 or hdu_index >= len(self.hdulist):
            raise ValueError("Invalid HDU index")

        return self.hdulist[hdu_index].header
    
    def get_data(self, hdu_index):
        "get the data of the fits file"
        return self.hdulist[hdu_index].data
    
    def close(self):
        "close the fits file"
        self.hdulist.close() 
