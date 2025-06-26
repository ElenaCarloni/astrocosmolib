import numpy as np
import Corrfunc 

class TwoPointIsotropic:

    def __init__(self):
        pass

    def set_separation(self, r_min, r_max, delta_r):
        self.r_min = r_min
        self.r_max = r_max
        self.delta_r = delta_r

        self.r_edges = np.arange(r_min, r_max + delta_r, delta_r)
        self.r_centers = 0.5 * (self.r_edges[:-1] + self.r_edges[1:])
        self.n_bins = len(self.r_centers)
    
    def compute_auto_pairs(self, xcoords, ycoords, zcoords, weights=None, periodic=True, boxsize=None, nthreads=1):
        if periodic and boxsize is None:
            raise ValueError("Box size must be provided for periodic boundary conditions")
        
        weight_type = None
        if weights is not None:
            weight_type = "pair_product"

        result = Corrfunc.theory.DD(True,
            nthreads=1,
            binfile=self.r_edges,
            X1=xcoords,
            Y1=ycoords,
            Z1=zcoords,
            weights1=weights,
            periodic=periodic,
            boxsize=boxsize,
            weight_type=weight_type
        )

        return result

    def compute_cross_pairs(self, x1coords, y1coords, z1coords, x2coords, y2coords, z2coords, weights1=None, weights2=None, periodic=True, boxsize=None, nthreads=1):
        if periodic and boxsize is None:
            raise ValueError("Box size must be provided for periodic boundary conditions")
        
        weight_type = None
        if (weights1 is not None) or (weights2 is not None):
            weight_type = "pair_product"

        result = Corrfunc.theory.DD(autocorr=False,
            nthreads=1,
            binfile=self.r_edges,
            X1=x1coords,
            Y1=y1coords,
            Z1=z1coords,
            weights1=weights1,
            X2=x2coords,
            Y2=y2coords,
            Z2=z2coords,
            weights2=weights2,
            periodic=periodic,
            boxsize=boxsize,
            weight_type=weight_type
        )

        return result

    def compute(self, 
                r_min,
                r_max,
                delta_r,
                x1coords,
                y1coords,
                z1coords,
                x2coords,
                y2coords,
                z2coords,
                weights1=None,
                weights2=None,
                periodic=True,
                boxsize=None,
                nthreads=1):
        
        self.set_separation(r_min, r_max, delta_r)

        self.DD = self.compute_auto_pairs(x1coords,
                                          y1coords,
                                          z1coords,
                                          weights=weights1,
                                          periodic=periodic,
                                          boxsize=boxsize,
                                          nthreads=nthreads)    

        self.RR = self.compute_auto_pairs(x2coords,
                                          y2coords,
                                          z2coords,
                                          weights=weights2,
                                          periodic=periodic,
                                          boxsize=boxsize,
                                          nthreads=nthreads)

        self.DR = self.compute_cross_pairs(x1coords,
                                           y1coords,
                                           z1coords,
                                           x2coords,
                                           y2coords,
                                           z2coords,
                                           weights1=weights1,
                                           weights2=weights2,
                                           periodic=periodic,
                                           boxsize=boxsize,
                                           nthreads=nthreads)

        DD = self.DD['npairs'].astype(np.float64)
        RR = self.RR['npairs'].astype(np.float64)
        DR = self.DR['npairs'].astype(np.float64)

        Nd = len(x1coords)**2
        Nr = len(x2coords)**2
        Ndr = len(x1coords) * len(x2coords)

        if (weights1 is not None): 
            DD *= self.DD['weightavg']
            Nd = np.sum(weights1)**2
        
        if (weights2 is not None):
            RR *= self.RR['weightavg']
            Nr = np.sum(weights2)**2

        if (weights1 is not None) and (weights2 is not None):
            DR *= self.DR['weightavg']
            Ndr = np.sum(weights1) * np.sum(weights2)

        self.xi = (DD/Nd + RR/Nr - 2*DR/Ndr) / (RR/Nr)