
import numpy as np
import numpy.ma as ma
from scipy.stats import sigmaclip
import matplotlib.pylab as plt


def pairPlot(xs, ys, val1, val2, plotrange, titl, prefix, suffix, valtype, units, nbins, inter, stitle=""):
    """

    plot a pair of plots, one showing the map, the other the histogram.

    Input: 
       xs,ys - coordinates of points
       val1 - value for the map, same dimensions as xs,ys
       val2 - value for the histogram
       plotrange - [min,max] for polotting range on both plots
       titl - title
       prefix - prefix for output files
       suffix - suffix for output files
       valtype - type of variable plotted (e.g., FWHM (x))
       units - unites (eg pixels)
       nbins - number of bins for histogram
       inter - interactive flag
       stitle - optional subbbtitle.

    REturns

       plots to png files, if inter=1 to screen

    """

    # set up plot
    fig, axes = plt.subplots(1, 2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    # scatter plot. size of points optimized for file/notebooks
    sc = axes[0].scatter(xs, ys, c=val1, marker="o", cmap='Purples',
                         lw=0, s=20, vmin=plotrange[0], vmax=plotrange[1])
    # label the axes
    axes[0].set_xlabel("X ("+units+")")
    axes[0].set_ylabel("Y ("+units+")")

    # calculate the bins
    binsize = (plotrange[1]-plotrange[0])/nbins
    bins = np.arange(plotrange[0], plotrange[1]+binsize, binsize)

    # histogram
    hi = axes[1].hist(val2, bins=bins)

    # labels
    axes[1].set_xlabel(valtype)
    axes[1].set_ylabel("N")

    plt.suptitle(valtype)

    #show and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")


def checkCentroids(xc, yc, cutrange, prefix, inter):
    """

    Quick plot of centroids to check results

    input

    xc,yc: centroid coordinates
    cutrange: limit the region of plotting if needed (for bad data)
    prefix: prefix for plots

    returns: plot, to screen and file

    """

    fig, ax = plt.subplots()

    # scatter plot

    ax.scatter(xc, yc)
    ax.set_aspect('equal')

    #display and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_checkpoints.png")


def checkMatched(xx, yy, xs, ys, prefix, inter):
    """

    quick plotting routine for measured centroids and pinhole coordinates

    input: 

    xx,yy mask coordiantes
    xs,ys: spot coordinates
    prefix: prefix for image files

    """

    fig, ax = plt.subplots()

    # scatter plot: centroids in circles, mask in red dots

    ax.scatter(xs, ys)
    ax.scatter(xx, yy, s=20, color='r')

    #save and show
    plt.savefig(prefix+"_checkpoints1.png")
    if(inter == 1):
        plt.show()


def plotVal(xs, ys, val, limit, plotrange, titl, prefix, suffix, units, inter, stitle=""):
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

    # a quick kludge to filter out bad quality points in poorly focussed images
    if(limit > 0):
        ind = np.where((val < limit) & (val > 0) & (val > 0))
    else:
        ind = np.arange(len(val))

    # scatter plot, with or without ragne limit

    fig, axes = plt.subplots()

    if(plotrange != None):
        sc = axes.scatter(xs[ind], ys[ind], c=val[ind], marker="o", cmap='Purples',
                          lw=0, s=20, vmin=plotrange[0], vmax=plotrange[1])

    else:
        sc = axes.scatter(xs[ind], ys[ind], c=val[ind], marker="o", cmap='Purples', lw=0, s=20)

    fig.colorbar(sc)
    plt.title(titl+stitle)
    plt.xlabel("X ("+units+")")
    plt.ylabel("Y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")


def checkPlots(files, inter, stitle=""):
    """

    TAkes a list of files and plots the mean and rms of the pixels by frame. 

    input: list of files

    output: plots

    """

    # text file for output
    nfiles = len(files)

    # set up variables
    av = []
    rms = []
    frame = np.arange(nfiles)+1

    # cycle through files
    for file in files:
        print(file)
        #read in image
        image = getImage(file)

        # calculate Stats
        rms.append(image.std())
        av.append(image.mean())

    # plot
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(frame, av, marker='d', linestyle="-")
    ax2.plot(frame, rms, marker='d', linestyle="-")
    ax1.set_title("Frame Mean"+stitle)
    ax2.set_title("Frame RMS"+stitle)
    ax2.set_xlabel("Frame")
    ax1.set_ylabel("Mean")
    ax2.set_ylabel("RMS")
    plt.tight_layout()
    if(inter == 1):
        plt.show()


def plotTransByFrame(fxFrameAv, fyFrameAv, peakFrameAv, sxAll, syAll, xdAll, ydAll, rotAll, prefix, inter, stitle=""):
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

    # get number of frames
    frames = np.arange(len(fxFrameAv))

    # first set - fwhms and translation (most useful)
    fig, axes = plt.subplots(1, 2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    axes[0].plot(frames, fxFrameAv, marker='d', linestyle="-", color="#1f77b4")
    axes[0].plot(frames, fyFrameAv, marker='s', linestyle="-", color="#ff7f0e")
    axes[0].set_title("FWHM Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("FHWM (pixels)")

    axes[1].plot(frames, xdAll-xdAll.mean(), marker='d', linestyle="-", color="#1f77b4")
    axes[1].plot(frames, ydAll-ydAll.mean(), marker='s', linestyle="-", color="#ff7f0e")
    axes[1].set_title("Translation Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Translation (pixels)")

    plt.savefig(prefix+"_byframe1.png")
    if(inter == 1):
        fig.show()

    # second set - peaks and bakcgrounds
    fig, axes = plt.subplots(1, 2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    axes[0].plot(frames, peakFrameAv, marker='d', linestyle="-")
    axes[0].set_title("Peak Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("Peak")

    axes[1].plot(frames, peakFrameAv, marker='d', linestyle="-")
    axes[1].set_title("Back Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Back")

    plt.savefig(prefix+"_byframe2.png")
    if(inter == 1):
        fig.show()

    # third set - scale and rotation
    fig, axes = plt.subplots(1, 2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    axes[0].plot(frames, sxAll, marker='d', linestyle="-")
    axes[0].plot(frames, syAll, marker='d', linestyle="-")
    axes[0].set_title("Scale Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("Scale")

    axes[1].plot(frames, rotAll, marker='d', linestyle="-")
    axes[1].set_title("Rotation Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Rotation (radians)")

    plt.savefig(prefix+"_byframe3.png")
    if(inter == 1):
        fig.show()

    # fourth set - nper frame


def plotImageStats(image, prefix, inter, stitle=""):
    """

    plot histogram of an image

    """

    back = sigmaclip(image, sigma=2, iters=2)
    backImage = back.mean()
    rmsImage = back.std()

    logbins = np.geomspace(image.min(), image.max(), 50)

    fig, ax = plt.subplots()
    ax.hist(image.flatten(), bins=logbins, histtype="step")
    print("here")
    plt.title("Histogram of Region of Interest"+stitle)
    plt.xlabel("Flux Value")
    plt.ylabel("N")
    plt.yscale("log")
    plt.savefig(prefix+"_stats.png")

    if(inter == 1):
        fig.show()

    return backImage, rmsImage


def plotValHist(val, plotrange, titl, prefix, suffix, valtype, inter, nbins, stitle=""):

    # routine to plot a histogram of a variable

    fig, ax = plt.subplots()

    binsize = (plotrange[1]-plotrange[0])/nbins
    bins = np.arange(plotrange[0], plotrange[1]+binsize, binsize)
    ax.hist(val, bins=bins)
    plt.title(titl+stitle)
    plt.xlabel(valtype)
    plt.ylabel("N")
    plt.savefig(prefix+"_"+suffix+".png")

    if(inter == 1):
        fig.show()


def makeReg(x, y, outfile):
    """

    Dump a series of points to ds9 region file

    """

    ff = open(outfile, "w")
    for i in range(len(x)):
        print("circle point ", x[i], y[i], file=ff)
    ff.close()


def movieCutout(files, xc, yc, sz, mmin, mmax, prefix):
    """

    generate a series of cutouts around an xy point

    """

    fig, ax = plt.subplots()

    for i in range(len(files)):
        image = getImage(files[i])

        ax.imshow(image[xc-sz:xc+sz, yc-sz:yc+sz], vmin=mmin, vmax=mmax)
        plt.savefig(prefix+str(i).zfill(2)+".png")


def quiverPlot(x, y, dx, dy):
    """

    plot a distortion map in quiver and colour format

    Input:

      x,y - positions
      dx,dy - difference from expected

    returns: plot to file nad screen


    """

    fig, ax = plt.subplots(1, 2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    ax[0].quiver(x, y, dx, dy)
    ax[0].set_xlabel("X (pixels)")
    ax[0].set_ylabel("Y (pixels)")

    # map shows the total deviation

    dist = np.sqrt((dx)**2+(dy)**2)
    sc = ax[1].scatter(x, y, c=dist)
    ax[1].set_xlabel("X (pixels)")
    fig.colorbar(sc, ax=ax[1])

    fig.suptitle("Distortion Map")

    fig.show()


def diagPlot(image, mCentroids, dx, dy):
    """

    plot the image wtih the distortion overplotted in images. 

    This needs a (working) functino to get the fibreID by clicking on it a point 

    Input: 
       iimage - image
       mCentroids - array with matched centroids
       dx,dy - distortion

    """

    fig, ax = plt.subplots()
    x = mCentroids[:, 1]
    y = mCentroids[:, 2]
    print(x.min(), x.max(), y.min(), y.max())
    ax.imshow(image, origin='lower')
    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min(), y.max()])
    ax.quiver(x, y, dx, dy, color='white', headlength=0, headaxislength=0)

    # this shoudl read the position, but doesn't work yet.
    def onM(event):
        xpos = event.xdata
        ypos = event.ydata
        print(xpos, ypos)
        fig.canvas.draw()

    fig.canvas.mpl_connect('button_press_event', onM)
    plt.show()


def checkThreshold(image, xrange, yrange):
    """

    quick polot to show the image and overplot the region for the threshold calculation

    """

    fig, ax = plt.subplots()
    ax.imshow(image)
    ax.scatter([yrange[0], yrange[0], yrange[1], yrange[1]], [xrange[0], xrange[1], xrange[0], xrange[1]])
    fig.show()
