import numpy as np
from os import system

def addone(var):
    '''add a term linearly to the end of an array'''
    varout = np.zeros(len(var)+1)
    varout[:-1] = var[:]
    num = var[-1] - (var[-2]-var[-1])
    if (round(num,4))<=0:  #if it's a possily negative number like u, v, do nothing, otherwise keep constant
        if np.min(var)>0:   
            num=var[-1]
    varout[-1]=num
    return varout
    

def makeramsin(press, tempk, rv, u, v, name):
    '''Makes a ramsin file with a sounding 
       Will make sure units make sense, though u and v must be components.
       Uses the name for the model directory as well as the ramsin filename. '''
    
    #probably take out??
    # if press[0] > 1000:
    #     press = np.delete(press,0)
    #     tempc = np.delete(tempc,0)
    #     rv = np.delete(rv,0)
    #     u = np.delete(u,0)
    #     v = np.delete(v,0)
    
    #fileout will be a temporary file, ramsout the final RAMSIN
    fileout = 'tmp-'+name+'.txt'
    ramsout = 'ramsin-'+name

    #check units
    if rv.max() < 1:
        rv = rv*1000.
    if tempk.max() <100:
        tempk = tempk + 273.15
    if press.max() > 1500:
        press = press/ 100. 

    #Stretch sounding high enough
    if press.min() < 60:
        toolow = np.where(press < 60)
        press = np.delete(press,toolow)
        tempk = np.delete(tempk,toolow)
        rv = np.delete(rv,toolow)
        u = np.delete(u,toolow)
        v = np.delete(v,toolow)

    for i in range(5):
        press = addone(press)
        tempk = addone(tempk)
        rv = addone(rv)
        u = addone(u)
        v = addone(v)
    
    out = open(fileout, 'w')
    out.write("   ! Sounding:  "+ name + "  \n")
    out.write('\n')
    out.write("   PS = ")
    for i, p in enumerate(press):
        out.write(str(p))
        out.write(", ")
        if (i> 1 and i % 6 == 0):
            out.write('\n')
            out.write("   ")
    out.write("\n \n   TS = ")
    for i, p in enumerate(tempk):
        out.write(str(round(p,2)))
        out.write(", ")
        if (i> 1 and i % 6 == 0):
            out.write('\n')
            out.write("   ")
    out.write("\n \n   RTS = ")
    for i, p in enumerate(rv):
        out.write(str(round(p,4)))
        out.write(", ")
        if (i> 1 and i % 6 == 0):
            out.write('\n')
            out.write("   ")
    out.write("\n \n   US = ")
    for i, p in enumerate(u):
        out.write(str(round(p,4)))
        out.write(", ")
        if (i> 1 and i % 6 == 0):
            out.write('\n')
            out.write("   ")
    out.write("\n \n   VS = ")
    for i, p in enumerate(v):
        out.write(str(round(p,4)))
        out.write(", ")
        if (i> 1 and i % 6 == 0):
            out.write('\n')
            out.write("   ")
    out.write("\n \n")

    out.close()
    
    catline = 'cat ramsin1 ' + fileout + ' ramsin2 > ' + ramsout
    system(catline)
    system('rm *txt')
    sedline = "sed -i '' 's/NAMENAME/"+name+"/g'  "  + ramsout
    system(sedline)
    sedline = "sed -i '' 's/XXXX/"+str(round(tempk[0],2))+"/g'  " + ramsout
    system(sedline)
