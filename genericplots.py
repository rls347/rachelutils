import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def timeseries(pcp, xs, outname, title, xlab, ylab):
    plt.plot(xs,pcp,linewidth=3)
    plt.xlabel(xlab, size=14)
    plt.ylabel(ylab, size=14)
    plt.title(title, size=20)
    plt.savefig(outname)
    plt.clf()

def movie(pcp, outname, title, height):
    numfiles = pcp.shape[0]
#    pcp = np.asarray(pcp)
    #pcp = np.swapaxes(pcp,1,2)
    print pcp.max(), pcp.min()
    xs = np.arange(pcp.shape[2])*.25
    if (pcp.shape[1]==pcp.shape[2]):
        ys = xs
    else:
        ys = height/1000.
    print xs.shape, ys.shape
    fig = plt.figure()
    z2 = pcp[0,:,:]
    print z2.shape
    
    cont2 = plt.contourf(xs,ys,z2,levels=np.linspace(np.min(pcp),np.max(pcp),20))
    c2 = plt.colorbar(cont2,label = title) 
#    plt.ylim(0,16)

    def animate(i):
        z2 = pcp[i,:,:]
        plt.clf()
        cont2 = plt.contourf(xs,ys,z2,levels=np.linspace(np.min(pcp),np.max(pcp),20))
        c2 = plt.colorbar(cont2,label = title)
#        plt.ylim(0,16)
        onetitle='Time ' + str(i)
        plt.title(onetitle)
        return cont2

    anim = animation.FuncAnimation(fig, animate, frames=numfiles)
    anim.save(outname+title+'.mp4')
    plt.clf()
