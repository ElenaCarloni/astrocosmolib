import numpy as np
from scipy.integrate import quad

class Integrate:
    def __init__(self, integrand, **kwargs):
        self.integrand = lambda x : integrand(x, **kwargs)
    
    def quadrature(self, a, b, **kwargs):
        I, _ = quad(self.integrand, a, b, **kwargs)
        return I

    def trapezoidal(self, a, b, N):
        x = np.linspace(a, b, N+1)
        y = self.integrand(x)
        h = (b - a)/N
        I = h*(0.5*y[0] + 0.5*y[-1] + np.sum(y[1:-1]))
        return I

    def __call__(self, method, a, b, **kwargs):
        if method=='trapezoidal':
            return self.trapezoidal(a, b, **kwargs)
        elif method=='quad':
            return self.quadrature(a, b, **kwargs)
        else:
            raise ValueError('Unknown integration method. Use "trapezoidal" or "quad".')
        