import numpy as np
from scipy.special import gamma
from rachelutils.hdfload import getvar

def calcgammadist(m,n,coef,ex,nu):
    '''Returns diameters and number for size distribution'''

    dvals = np.linspace(1,500,5000)*1e-6
    d=(m/(n*coef))**(1./ex)
    dn = (gamma(nu)/gamma(nu+ex))**(1./ex)*d
    func = (1/gamma(nu))*((dvals/dn)**ex)*(1/dn)*np.exp(-1*dvals/dn)
    return dvals, func

def calcd(fil,var):
    nnames = {  'pristine':'pris_concen_kg',
                'cloud':'cloud_concen_mg',
                'aggregates':'agg_concen_kg',
                'snow':'snow_concen_kg',
                'drizzle':'drizzle_concen_mg',
                'hail':'hail_concen_kg',
                'rain':'rain_concen_kg',
                'graupel':'graup_concen_kg'
            }
    cfmas = { 'pristine': 110.8,
              'cloud': 524.,
              'aggregates': 0.496,
              'snow': 0.002739,
              'drizzle': 524.,
              'hail': 471.,
              'rain': 524.,
              'graupel': 157.   }

    pwmas = { 'cloud': 3.,
              'drizzle': 3.,
              'rain': 3.,
              'pristine': 2.91,
              'snow': 1.74,
              'aggregates': 2.4,
              'graupel': 3.,
              'hail': 3.  }
    mass = getvar(fil,var)
    num = getvar(fil,nnames[var])
    
    coef = cfmas[var]
    ex = pwmas[var]

    d=(mass/(num*coef))**(1./ex)
    d[num==0]=0.0
    return d
