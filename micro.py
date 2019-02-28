import numpy as np
from scipy.special import gamma


def calcgammadist(m,n,coef,ex,nu):
    '''Returns diameters and number for size distribution'''

    dvals = np.linspace(1,500,5000)*1e-6
    d=(m/(n*coef))**(1./ex)
    dn = (gamma(nu)/gamma(nu+ex))**(1./ex)*d
    func = (1/gamma(nu))*((dvals/dn)**ex)*(1/dn)*np.exp(-1*dvals/dn)
    return dvals, func

