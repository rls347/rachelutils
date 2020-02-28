import numpy as np
from scipy.special import gamma
from rachelutils.hdfload import getvar

def getdefault(var):
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
  
    nu = {    'cloud': 4.,
              'drizzle': 4.,
              'rain': 2.,
              'pristine': 2.,
              'snow': 2.,
              'aggregates': 2.,
              'graupel':2.,
              'hail': 2.  }
    return cfmas[var],pwmas[var],nu[var]

def calcgammadist(m,n,coef,ex,nu):
    '''Returns diameters and number for size distribution'''
    dvals = np.linspace(1,10000,5000)*1e-6
    with np.errstate(invalid='ignore',divide='ignore'):
        d=(m/(n*coef))**(1./ex)
        
    dn = (gamma(nu)/gamma(nu+ex))**(1./ex)*d
    func = (1/gamma(nu))*((dvals/dn)**(nu-1))*(1/dn)*np.exp(-1*dvals/dn)
    return dvals, func
    
def gammadist_diam(d,ex,nu):
    '''Returns diameters and number for size distribution given mean mass diameter'''

    #dvals = np.linspace(1,500,5000)*1e-6
    dvals = np.arange(1,5000)*1e-7
    d=(m/(n*coef))**(1./ex)

    dn = (gamma(nu)/gamma(nu+ex))**(1./ex)*d

    func = (1/gamma(nu))*((dvals/dn)**ex)*(1/dn)*np.exp(-1*dvals/dn)
    
    return dvals, func

def calcd(fil,var,x=None):
    nnames = {  'pristine':'pris_concen_kg',
                'cloud':'cloud_concen_mg',
                'aggregates':'agg_concen_kg',
                'snow':'snow_concen_kg',
                'drizzle':'drizzle_concen_mg',
                'hail':'hail_concen_kg',
                'rain':'rain_concen_kg',
                'graupel':'graup_concen_kg'
            }

    mass = getvar(fil,var)/1000.
    num = getvar(fil,nnames[var])
    if var == 'cloud' or var == 'drizzle':
        num=num*1000.*1000.
    
    coef, ex, nu = getdefault(var)
    if x is not None:
        coef=x
    
    with np.errstate(invalid='ignore',divide='ignore'):
        d=(mass/(num*coef))**(1./ex)
        d[num==0]=0.0
    
    return d
    
def getdefault(var):
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
  
    nu = {    'cloud': 4.,
              'drizzle': 4.,
              'rain': 2.,
              'pristine': 2.,
              'snow': 2.,
              'aggregates': 2.,
              'graupel':2.,
              'hail': 2.  }
    return cfmas[var],pwmas[var],nu[var]
              

def intsizedist(fil,var,coord=None):
    nnames = {  'pristine':'pris_concen_kg',
                'cloud':'cloud_concen_mg',
                'aggregates':'agg_concen_kg',
                'snow':'snow_concen_kg',
                'drizzle':'drizzle_concen_mg',
                'hail':'hail_concen_kg',
                'rain':'rain_concen_kg',
                'graupel':'graup_concen_kg'
            }
    #If no coord given, uses None, which just returns full array + a dimension
    
    
    mass = np.array(getvar(fil,var)[coord]/1000.)
    num = np.array(getvar(fil,nnames[var])[coord])
    
        
    if var == 'cloud' or var == 'drizzle':
        num=num*1000.*1000.
    coef, ex, nu = getdefault(var) 

    m=mass.flatten()
    n=num.flatten()
    tval = np.zeros_like(m)
    bval = np.zeros_like(m)
    for i in range(len(m)):
        if n[i] >0:
            dvals, func = calcgammadist(m[i],n[i],coef,ex,nu)
            diffr = np.zeros_like(dvals)
            diffr[1:] = np.diff(dvals)/2
            tval[i] = np.sum( (dvals/2.)**3 * diffr * func * 3.14 )
            bval[i] = np.sum( (dvals/2.)**2 * diffr * func * 3.14 )
    bot = np.reshape(bval,mass.shape)
    top = np.reshape(tval,mass.shape)

    return mass, num, bot,top

def getcloudtop(fil):
    q=getvar(fil,'total_cond')
    c=np.where(q>.01)
    q[c]=100.
    inds = q.shape[0]-1-np.argmax(q[::-1,:,:],0)
    inds[np.max(q,0)<.01]=0

    return inds

def intsizedist(fil,var,coord=None):
    nnames = {  'pristine':'pris_concen_kg',
                'cloud':'cloud_concen_mg',
                'aggregates':'agg_concen_kg',
                'snow':'snow_concen_kg',
                'drizzle':'drizzle_concen_mg',
                'hail':'hail_concen_kg',
                'rain':'rain_concen_kg',
                'graupel':'graup_concen_kg'
            }
    #If no coord given, uses None, which just returns full array + a dimension
    
    
    mass = np.array(getvar(fil,var)[coord]/1000.)
    num = np.array(getvar(fil,nnames[var])[coord])
    
        
    if var == 'cloud' or var == 'drizzle':
        num=num*1000.*1000.
    coef, ex, nu = getdefault(var) 

    m=mass.flatten()
    n=num.flatten()
    tval = np.zeros_like(m)
    bval = np.zeros_like(m)
    for i in range(len(m)):
        if n[i] >0:
            dvals, func = calcgammadist(m[i],n[i],coef,ex,nu)
            diffr = np.zeros_like(dvals)
            diffr[1:] = np.diff(dvals)/2
            tval[i] = np.sum( (dvals/2.)**3 * diffr * func * 3.14 )
            bval[i] = np.sum( (dvals/2.)**2 * diffr * func * 3.14 )
    bot = np.reshape(bval,mass.shape)
    top = np.reshape(tval,mass.shape)

    return mass, num, bot,top

