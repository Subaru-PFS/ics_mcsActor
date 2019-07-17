
import numpy as np
import numpy.ma as ma
from scipy.stats import sigmaclip
import matplotlib.pylab as plt


def pairPlot(xs,    ys,    val1,val2,           plotrange,titl,      prefix, suffix, valtype,   units,   nbins,inter,stitle=""):

    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)


    sc=axes[0].scatter(xs,ys,c=val1,marker="o",cmap='Purples',lw=0,s=20,vmin=plotrange[0],vmax=plotrange[1])
    axes[0].set_xlabel("X ("+units+")")
    axes[0].set_ylabel("Y ("+units+")")
    
    binsize=(plotrange[1]-plotrange[0])/nbins
    bins=np.arange(plotrange[0],plotrange[1]+binsize,binsize)

    hi=axes[1].hist(val2,bins=bins)

    axes[1].set_xlabel(valtype)
    axes[1].set_ylabel("N")

    plt.suptitle(valtype)

    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")


def plotTransFrame(fxFrameAv,fyFrameAv,peakFrameAv,sxAll,syAll,xdAll,ydAll,rotAll,prefix,inter,stitle=""):

    fig,axes=plt.subplots(1,2)

    ax[0].plot(frames,fxFrameAv,marker='d',linestyle="-",color="#1f77b4")
    ax[0].plot(frames,fxFrameAv,marker='d',linestyle="-",color="#ff7f0e")

    ax[0].set_xlabel("Frame #")
    ax[0].set_ylabel("FWHM (pixels)")


    ax[1].plot(frames,peakFrameAv,marker='d',linestyle="-")

    ax[0].set_xlabel("Frame #")
    ax[0].set_ylabel("FWHM (pixels)")



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
    #if(cutrange==1):
    #    np.xlim([-8,344])
    #    np.ylim([-8,344])

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

def plotDistortion(c,c1,pts1,pts2,x1,y1,diffx,diffy,fxs,fys,peaks,limit,prefix,units,inter,inst_scale,stitle=""):

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
    
    if(limit !=0):
        ind=np.where((c < limit) & (fxs < 10) & (fys < 10) & (fxs > 0) & (fys > 0) & (peaks != 0))
    else:
        ind=np.arange(len(c))
        
    #quiver plot
    
    fig,ax = plt.subplots(facecolor='g')

    x1=ma.masked_values(x1,0)
    y1=ma.masked_values(y1,0)

    print(x1.min(),y1.min())
   
    ax.quiver(pts1[:,ind,0].ravel(),pts1[:,ind,1].ravel(),-diffx[ind],-diffy[ind])
    #ax.quiver(x1,y1,-diffx[ind],-diffy[ind])
    plt.xlabel("x ("+units+")")
    plt.ylabel("y ("+units+")")
    plt.title("Distortion"+stitle)
    plt.axis('equal')
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_distortion.png")

    #t1=np.partition(c.flatten(), -2)[-2]
    #t2=np.partition(c1.flatten(), -2)[-2]
    
    fig,ax = plt.subplots()
    #plot map in mm
    #sc=ax.scatter(pts1[:,ind,0].ravel(),pts1[:,ind,1].ravel(),marker='o',c=c[ind],cmap='plasma',s=25)
    sc=ax.scatter(pts1[:,ind,0].ravel(),pts1[:,ind,1].ravel(),marker='o',c=c[ind]/inst_scale,cmap='plasma',s=25)
    fig.colorbar(sc)
    plt.title("Distortion ("+units+")"+stitle)
    plt.xlabel("x ("+units+")")
    plt.ylabel("y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_distortion_col.png")

    fig,ax = plt.subplots()
    #plot map in %
    sc=ax.scatter(pts1[:,ind,0].ravel(),pts1[:,ind,1].ravel(),marker='o',c=c1[ind]/inst_scale,cmap='plasma',s=25)
    fig.colorbar(sc)
    plt.title("Distortion (% of field size)"+stitle)
    plt.xlabel("x ("+units+")")
    plt.ylabel("y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_distortion_col1.png")


def plotVal(xs,ys,val,limit,plotrange,titl,prefix,suffix,units,inter,stitle=""):

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
    plt.title(titl+stitle)
    plt.xlabel("X ("+units+")")
    plt.ylabel("Y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")

def checkPlots(files,inter,stitle=""):

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
    ax1.set_title("Frame Mean"+stitle)
    ax2.set_title("Frame RMS"+stitle)
    ax2.set_xlabel("Frame")
    ax1.set_ylabel("Mean")
    ax2.set_ylabel("RMS")
    plt.tight_layout()
    if(inter == 1):
        plt.show()

def plotTransByFrame(fxFrameAv,fyFrameAv,peakFrameAv,sxAll,syAll,xdAll,ydAll,rotAll,prefix,inter,stitle=""):

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


    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)
    
    axes[0].plot(frames,fxFrameAv,marker='d',linestyle="-",color="#1f77b4")
    axes[0].plot(frames,fyFrameAv,marker='s',linestyle="-",color="#ff7f0e")
    axes[0].set_title("FWHM Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("FHWM (pixels)")

    axes[1].plot(frames,xdAll-xdAll.mean(),marker='d',linestyle="-",color="#1f77b4")
    axes[1].plot(frames,ydAll-ydAll.mean(),marker='s',linestyle="-",color="#ff7f0e")
    axes[1].set_title("Translation Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Translation (pixels)")
    
    plt.savefig(prefix+"_byframe1.png")
    if(inter == 1):
        fig.show()

    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)
        
    axes[0].plot(frames,peakFrameAv,marker='d',linestyle="-")
    axes[0].set_title("Peak Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("Peak")

    axes[1].plot(frames,peakFrameAv,marker='d',linestyle="-")
    axes[1].set_title("Back Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Back")
    
    plt.savefig(prefix+"_byframe2.png")
    if(inter == 1):
        fig.show()

    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)
 
    axes[0].plot(frames,sxAll,marker='d',linestyle="-")
    axes[0].plot(frames,syAll,marker='d',linestyle="-")
    axes[0].set_title("Scale Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("Scale")
    
    axes[1].plot(frames,rotAll,marker='d',linestyle="-")
    axes[1].set_title("Rotation Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Rotation (radians)")

    plt.savefig(prefix+"_byframe3.png")
    if(inter == 1):
         fig.show()

def plotFailed(ncent,nqual,ncentR,nqualR,inr,prefix,inter,stitle=""):

    """

    """

    
    #get number of frames
    frames=np.arange(len(ncent))
         
    fig,ax=plt.subplots()
    ax.plot(frames,ncent,marker='d',linestyle="-")
    plt.title("Number of Detected Points by Frame"+stitle)
    plt.xlabel("Frame #")
    plt.ylabel("Number of Detected Points")

    if(inr==180):
        plt.axhline(y=3444,linestyle="-",color="red")
    else:
        plt.axhline(y=3542,linestyle="-",color="red")

        plt.savefig(prefix+"_centbyframe.png")

    if(inter == 1):
         fig.show()
    
    fig,ax=plt.subplots()
    ax.plot(frames,nqual,marker='d',linestyle="-")
    plt.title("Number of Flagged Fits by Frame"+stitle)
    plt.xlabel("Frame #")
    plt.ylabel("Number of Flagged Fits")
    plt.savefig(prefix+"_flaggbyframe.png")

    if(inter == 1):
         fig.show()
    
         
def plotImageStats(image,prefix,inter,stitle=""):


    back = sigma_clip(image, sigma=2, iters=2)
    backImage=back.mean()
    rmsImage=back.std()

    logbins = np.geomspace(image.min(), image.max(), 50)
    
    fig,ax = plt.subplots()
    ax.hist(image.flatten(),bins=logbins,histtype="step")
    print("here")
    plt.title("Histogram of Region of Interest"+stitle)
    plt.xlabel("Flux Value")
    plt.ylabel("N")
    plt.yscale("log")
    plt.savefig(prefix+"_stats.png")

    if(inter == 1):
        fig.show()

    return backImage,rmsImage

def plotValHist(val,plotrange,titl,prefix,suffix,valtype,inter,nbins,stitle=""):

    #routine to plot a histogram of a variable

    fig,ax = plt.subplots()

    binsize=(plotrange[1]-plotrange[0])/nbins
    bins=np.arange(plotrange[0],plotrange[1]+binsize,binsize)
    ax.hist(val,bins=bins)
    plt.title(titl+stitle)
    plt.xlabel(valtype)
    plt.ylabel("N")
    plt.savefig(prefix+"_"+suffix+".png")

    if(inter == 1):
        fig.show()


def makeReg(x,y,outfile):

    """

    Dump a series of points to ds9 region file

    """

    ff=open(outfile,"w")
    for i in range(len(x)):
        print("circle point ",x[i],y[i],file=ff)
    ff.close()

def movieCutout(files,xc,yc,sz,mmin,mmax,prefix):

    """
    
    generate a series of cutouts around an xy point

    """
    
    fig,ax=plt.subplots()

    for i in range(len(files)):
        image=getImage(files[i])
        
        ax.imshow(image[xc-sz:xc+sz,yc-sz:yc+sz],vmin=mmin,vmax=mmax)
        plt.savefig(prefix+str(i).zfill(2)+".png")


    
