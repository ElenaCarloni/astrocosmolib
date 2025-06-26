import numpy as np

class LambdaCDM:
    def __init__(self, H0, Omega_m, Omega_Lambda, Omega_radiation=0.0):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_Lambda = Omega_Lambda
        self.Omega_radiation = Omega_radiation
        self.Omega_k = 1.0 - (self.Omega_m + self.Omega_Lambda + self.Omega_radiation)
    
    def E(self, z):
        return np.sqrt(
            self.Omega_m * (1 + z)**3 +
            self.Omega_radiation * (1 + z)**4 +
            self.Omega_k * (1 + z)**2 +
            self.Omega_Lambda
        ) 
    
    def H(self, z):
        return self.H0 * self.E(z)

class FlatLambdaCDM:
    def __init__(self, H0, Omega_m, Omega_radiation=0.0):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_radiation = Omega_radiation
        self.Omega_k = 0
        self.Omega_Lambda = 1 - self.Omega_m - self.Omega_radiation
    
    def E(self, z):
        return np.sqrt(
            self.Omega_m * (1 + z)**3 +
            self.Omega_radiation * (1 + z)**4 +
            self.Omega_Lambda
        )
    
    def H(self, z):
        return self.H0 * self.E(z)

class w0waFlatCDM:
    def __init__(self, H0, Omega_m, Omega_radiation=0.0, w0=-1, wa=0):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_DarkEnergy = 1 - Omega_m - Omega_radiation
        self.Omega_radiation = Omega_radiation
        self.Omega_k = 0
        self.w0 = w0
        self.wa = wa
        self.Omega_k = 0
    
    def E(self, z):
        return np.sqrt(
            self.Omega_m * (1 + z)**3 +
            self.Omega_radiation * (1 + z)**4 +
            self.Omega_DarkEnergy * (1 + z)**(3 * (1 + self.w0 + self.wa)) * np.exp(-3 * self.wa * z/(1 + z))
        )
    
    def H(self, z):
        return self.H0 * self.E(z)
    