
import matplotlib.pylab as plt


import numpy as np
import pylab as py
from matplotlib.cm import plasma
from scipy.stats import sigmaclip
import astropy
import matplotlib 
from astropy.stats import sigma_clip


def getImage(filename):

    """

    Simple wrapper to read an image from file

    """
    
    image=astropy.io.fits.getdata(filename)

    return image
    
def to_fits(filename):

    #quick routine to convert raw to FITS

    a=np.fromfile(filename,dtype='uint16')
    image=a.reshape([5778,8960])

    pf.writeto(filename+".fits",image)
    print(image.min(),image.mean(),image.max())


def checkCentroids(xc,yc,cutrange,prefix,inter):

    """

    Quick plot of centroids to check results

    input

    xc,yc: centroid coordinates
    cutrange: limit the region of plotting if needed (for bad data)
    prefix: prefix for plots

    returns: plot, to screen and file

    """

    fig,ax = plt.subplots()

    #scatter plot
    
    ax.scatter(xc,yc)

    #set limit if desired
    if(cutrange==1):
        np.xlim([-8,344])
        np.ylim([-8,344])

    #plt.axes().set_aspect('equal', 'datalim')
    ax.set_aspect('equal')
    #display and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_checkpoints.png")

def checkMatched(xx,yy,xs,ys,prefix,inter):

    """

    quick plotting routine for measured centroids and pinhole coordinates

    input: 

    xx,yy mask coordiantes
    xs,ys: spot coordinates
    prefix: prefix for image files

    """

    
    fig,ax = plt.subplots()
    
    #scatter plot: centroids in circles, mask in red dots
    
    ax.scatter(xs,ys)
    ax.scatter(xx,yy,s=20,color='r')

    #save and show
    plt.savefig(prefix+"_checkpoints1.png")
    if(inter == 1):
        plt.show()

def plotDistortion(c,c1,pts1,pts2,diffx,diffy,fxs,fys,peaks,limit,prefix,units,inter):

    """

    Distortion plots

    input:

    c,c1,pts1,pts2,diffx,diffy: as returned by simpleDistortion
    fxs,fyx: fwhms
    peaks: peak values

    limit: upper limit for filtering out bad data (in c units)

    prefix: prefix for output files

    """

    #a quick kludge to filter out bad quality points in poorly focussed images
    
    if(limit > 0):
        ind=np.where((c < limit) & (fxs < 10) & (fys < 10) & (fxs > 0) & (fys > 0) & (peaks != 0))
    else:
        ind=np.arange(len(c))
        
    #quiver plot
    
    fig,ax = plt.subplots(facecolor='g')

    ax.quiver(pts1[0,ind,0],pts1[0,ind,1],-diffx[ind],-diffy[ind])
    plt.xlabel("x ("+units+")")
    plt.ylabel("y ("+units+")")
    plt.title("Distortion")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_distortion.png")

    #t1=np.partition(c.flatten(), -2)[-2]
    #t2=np.partition(c1.flatten(), -2)[-2]
    
    fig,ax = plt.subplots()
    #plot map in mm
    sc=ax.scatter(pts2[0,ind,0],pts2[0,ind,1],marker='s',c=c[ind],cmap='plasma',s=200)
    fig.colorbar(sc)
    plt.title("Distortion ("+units+")")
    plt.xlabel("x ("+units+")")
    plt.ylabel("y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_distortion_col.png")

    fig,ax = plt.subplots()
    #plot map in %
    sc=ax.scatter(pts2[0,ind,0],pts2[0,ind,1],marker='s',c=c1[ind]*4,cmap='plasma',s=200)
    fig.colorbar(sc)
    plt.title("Distortion (% of field size)")
    plt.xlabel("x ("+units+")")
    plt.ylabel("y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_distortion_col1.png")


def plotVal(xs,ys,val,limit,plotrange,titl,prefix,suffix,units,inter):

    """
    
    routine for scatter plot of a variable

    input

    xs,ys: coordinates

    val: variable to plot
    limit: upper limit to filter out of plots
    plotrange: colour range in [vmin,vmax] format, or None for default
    title: title for plot
    prefix: prefix for output files
    suffix: suffix for output files

    """

    
    #a quick kludge to filter out bad quality points in poorly focussed images
    if(limit > 0):
        ind=np.where((val < limit) & (val > 0) & (val > 0))
    else:
        ind=np.arange(len(val))

    #scatter plot, with or without ragne limit
    
    fig, axes = plt.subplots()
    
    if(plotrange != None):
        sc=axes.scatter(xs[ind],ys[ind],c=val[ind],marker="o",cmap='Purples',lw=0,s=20,vmin=plotrange[0],vmax=plotrange[1])
    
    else:
        sc=axes.scatter(xs[ind],ys[ind],c=val[ind],marker="o",cmap='Purples',lw=0,s=20)

    fig.colorbar(sc)
    plt.title(titl)
    plt.xlabel("X ("+units+")")
    plt.ylabel("Y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")

def checkPlots(files,inter):

    """

    TAkes a list of files and plots the mean and rms of the pixels by frame. 

    input: list of files

    output: plots

    """

    
    #text file for output
    nfiles=len(files)
    

    #set up variables
    av=[]
    rms=[]
    frame=np.arange(nfiles)+1

    #cycle through files
    for file in files:
        print(file)
        #read in image
        image=getImage(file)

        #calculate Stats
        rms.append(image.std())
        av.append(image.mean())

    #plot
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(frame,av,marker='d',linestyle="-")
    ax2.plot(frame,rms,marker='d',linestyle="-")
    ax1.set_title("Frame Mean")
    ax2.set_title("Frame RMS")
    ax2.set_xlabel("Frame")
    ax1.set_ylabel("Mean")
    ax2.set_ylabel("RMS")
    plt.tight_layout()
    if(inter == 1):
        plt.show()

def plotTransByFrame(fxFrameAv,fyFrameAv,peakFrameAv,sxAll,syAll,xdAll,ydAll,rotAll,prefix,inter):

    """

    plot the calculated transformations values and averages by frame number.
    takes output generated by getTransByFrame

    input:
    fxFrameAv,fyFrameAv: average FWHM by frame
    peakFrameAv: peak value by frame

    sxAll,syAll: scale in x and y direction
    xdAll,ydAll: trnaslation in x and y
    rotAll: rotation 

    output: plots

    """

    #get number of frames
    frames=np.arange(len(fxFrameAv))

    fig,ax=plt.subplots()
    ax.plot(frames,fxFrameAv,marker='d',linestyle="-")
    ax.plot(frames,fyFrameAv,marker='d',linestyle="-")
    plt.title("FWHM Average by Frame")
    plt.xlabel("Frame #")
    plt.ylabel("FHWM (pixels)")
    plt.savefig(prefix+"_fwhmbyframe.png")
    if(inter == 1):
        fig.show()

    fig,ax=plt.subplots()
    ax.plot(frames,peakFrameAv,marker='d',linestyle="-")
    plt.title("Peak Average by Frame")
    plt.xlabel("Frame #")
    plt.ylabel("Peak")
    plt.savefig(prefix+"_peakbyframe.png")

    if(inter == 1):
        fig.show()

    fig,ax=plt.subplots()
    ax.plot(frames,sxAll,marker='d',linestyle="-")
    ax.plot(frames,syAll,marker='d',linestyle="-")
    plt.title("Scale Average by Frame")
    plt.xlabel("Frame #")
    plt.ylabel("Scale")
    plt.savefig(prefix+"_scalebyframe.png")

    if(inter == 1):
        fig.show()
    
    fig,ax=plt.subplots()
    ax.plot(frames,xdAll,marker='d',linestyle="-")
    ax.plot(frames,ydAll,marker='d',linestyle="-")
    plt.title("Translation Average by Frame")
    plt.xlabel("Frame #")
    plt.ylabel("Translation (pixels)")
    plt.savefig(prefix+"_transbyframe.png")

    if(inter == 1):
        fig.show()

    fig,ax=plt.subplots()
    ax.plot(frames,rotAll,marker='d',linestyle="-")
    plt.title("Rotation Average by Frame")
    plt.xlabel("Frame #")
    plt.ylabel("Rotation (radians)")
    plt.savefig(prefix+"_rotbyframe.png")

    if(inter == 1):
         fig.show()

def plotImageStats(image,prefix,inter):


    back = sigma_clip(image, sigma=2, iters=2)
    backImage=back.mean()
    rmsImage=back.std()

    logbins = np.geomspace(image.min(), image.max(), 50)
    
    fig,ax = plt.subplots()
    ax.hist(image.flatten(),bins=logbins,histtype="step")
    print("here")
    plt.title("Histogram of Region of Interest")
    plt.xlabel("Flux Value")
    plt.ylabel("N")
    plt.yscale("log")
    plt.savefig(prefix+"_stats.png")

    if(inter == 1):
        fig.show()

    return backImage,rmsImage
