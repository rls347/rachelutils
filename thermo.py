import h5py as hdf
import numpy as np
from rachelutils.hdfload import getvar, getdz

G = 9.8
CP = 1004.
RD = 287.
EPS = .622
P0 = 1000.
LV = 2.5e6
CV = CP-RD
GAMMAD = G/CP
KAPPA = RD/CP

def satmixratio(tempk,press):
    '''expects temperature in kelvin, pressure in millibars'''
    tc = tempk-273.15
    es = 6.112 * np.exp( (17.67*tc)/(tc+243.5) )  #in mb
    rs = (EPS *es) / (press-es)
    return rs

def find_lcl_sfc(height,tempk,q,press):  
    '''expects height, temperature in kelvin, q in g/kg, pressure in mb'''
    tprev = tempk[1]
    qsfc = q[1]/1000.
    for z in range(2,len(q)):
        p=press[z]
        t = tprev - ((height[z]-height[z-1]) * GAMMAD)
        rs = satmixratio(t,p)
        tprev = t
        if rs <= qsfc:
            lcl = height[z]
            tlcl = t
            break
    return z, lcl, tlcl

def find_lcl_ml100(height,tempk,q,press):
    '''expects height, temperature in kelvin, q in g/kg, pressure in mb'''
    psfc = press[1]
    lev100 = np.argmin(np.abs(press-(psfc-100)))

    rho = (press*100)/(RD*tempk)
    meanz = np.mean(height[1:lev100+1]*rho[1:lev100+1])/np.mean(rho[1:lev100+1])
    meant = np.mean(tempk[1:lev100+1]*rho[1:lev100+1])/np.mean(rho[1:lev100+1])
    meanq = np.mean(q[1:lev100+1]*rho[1:lev100+1])/np.mean(rho[1:lev100+1])

    zstart = np.argmin(np.abs(height-meanz))
    tprev = meant
    qsfc = meanq/1000.
    for z in range(zstart,len(q)):
        p=press[z]
        t = tprev - ((height[z]-height[z-1]) * GAMMAD)
        rs = satmixratio(t,p)
        tprev = t
        if rs <= qsfc:
            lcl = height[z]
            tlcl = t
            break
    return z, lcl, tlcl

def find_lcl_mostunstable(height,tempk,q,press):
    '''expects height, temperature in kelvin, q in g/kg, pressure in mb'''
    rho = (press*100)/(RD*tempk)
    theta = tempk*(P0/press)**(KAPPA)
    thetae = theta * np.exp((LV*(q/1000.))/(CP*tempk))
    thetae = thetae[height<5000]

    zstart = np.argmax(thetae)
    print zstart
    if zstart ==0:
        zstart =1
    tprev = tempk[zstart] 
    qsfc = q[zstart]/1000.
    for z in range(zstart,len(q)):
        p=press[z]
        t = tprev - ((height[z]-height[z-1]) * GAMMAD)
        rs = satmixratio(t,p)
        tprev = t
        if rs <= qsfc:
            lcl = height[z]
            tlcl = t
            break
    return z, lcl, tlcl

def moistadiabat(height,press,zlcl,tlcl,tempk,q):
    parcel = np.zeros_like(press)
    parcel[zlcl]=tlcl
    for z in range(zlcl+1,len(press)):
        t=parcel[z-1]
        #qv=q[z-1]/1000.
        qv = satmixratio(t,press[z-1])
        lqrt = (LV*qv)/(RD*t)
        lqrt2 = ((LV**2)*qv*(EPS+qv)) / (RD*(t**2))
        gammaPA = G * ( ((1+qv)*(1+lqrt)) / (CP + (qv*CV) + lqrt2) )
        parcel[z]= parcel[z-1] - ((height[z]-height[z-1]) * gammaPA)
    for z in range(zlcl-1,-2,-1):
        first = parcel[z+1] + ( (height[z+1]-height[z]) * GAMMAD)
        parcel[z]=first
    return parcel


def capecalc(height,tempk,q,press,dz,flag):
    if flag == 'sfc':
        zlcl, lcl, tlcl = find_lcl_sfc(height,tempk,q,press)
    elif flag == 'mu':
        zlcl, lcl, tlcl = find_lcl_mostunstable(height,tempk,q,press)
    else:
        zlcl, lcl, tlcl = find_lcl_ml100(height,tempk,q,press)

    parcel = moistadiabat(height,press,zlcl,tlcl,tempk,q)
    parceltv = parcel*(1+(.61*satmixratio(parcel,press)))
    Tv = tempk*(1+(.61*(q/1000.)))
    posarea = np.where(parceltv>Tv)
    intvar = (parceltv-Tv)*(G/Tv)
    lfc = np.min(np.asarray(posarea))
    el = np.max(np.asarray(posarea))
    intvarpos=intvar[lfc:el+1]
    zpos = height[lfc:el+1]

#    cape1 = np.sum(intvar[posarea])
    cape = np.trapz(intvarpos, dx=np.diff(zpos))

    return cape



def get_cape(filename,flag,coords = None):
    try:
        fil = hdf.File(filename, 'r')
        height = getvar(fil, 'z_coords')
        tempk = getvar(fil, 'tempk')
        q = getvar(fil, 'vapor')
        press = getvar(fil, 'press')
        dz = getdz(fil)
        fil.close()
    except:
        print 'error reading variables'

    if type(coords) is tuple:
        c0,c1 = coords
        outvar = capecalc(height,tempk[:,c0,c1],q[:,c0,c1],press[:,c0,c1],dz,flag)
    else:
        outvar = np.zeros((tempk.shape[1],tempk.shape[2]))
        for c0 in range(outvar.shape[0]):
            for c1 in range(outvar.shape[1]):
                outvar[c0,c1] = capecalc(height,tempk[:,c0,c1],q[:,c0,c1],press[:,c0,c1],dz,flag)

    return outvar
