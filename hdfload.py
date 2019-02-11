import h5py as hdf
import numpy as np

def getvar(fil, varname):
	try:
		var = np.squeeze(fil[varname].value)
	except:
		filey = hdf.File(fil, 'r')
		var = np.squeeze(filey[varname].value)
		filey.close()

	return var

def getdz(fil):
    height = getvar(fil, 'z_coords')
    dz = np.zeros_like(height)
#    for i in range(1,len(dz)-1):
#        dz[i]=.5*(height[i+1]-height[i-1])
    
    dz[:-1]=np.diff(height)
    dz[-1]=dz[-2]
    return dz

def getrho(fil):
    press = getvar(fil, 'press')
    tempk = getvar(fil, 'tempk')
    rho = (press*100) / (287*tempk)
    return rho

def meanprof(fil, varname):
    if varname == 'rho':
        var = getrho(fil)
    else:
        var = getvar(fil,varname)
    varout = np.mean(np.mean(var,1),1)
    return varout
