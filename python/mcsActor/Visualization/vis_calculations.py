

import numpy as np
import numpy.ma as ma
import cv2
import astropy
from scipy.stats import sigmaclip
from scipy.optimize import curve_fit


# this is so that you can run the code in a MHS environment, or stand
# alone

try:
    import mcsActor.Visualization.vis_plotting as visplot
except:
    import vis_plotting as visplot

try:
    import mcsActor.Visualization.vis_coordinates as viscoords
except:
    import vis_coordinates as viscoords

try:
    import mcsActor.mpfitCentroid.centroid as centroid
except:
    import centroid as centroid


def getCentroids(image, fwhm, boxsize, thresh, rl, rh, sl, sh):
    """

    wrapper to run centroid command, and convert the output to arrays

    input: 

    image (numpy array)
    FWHM, boxsize, thresh, rl,rh,sl,sh: parameters for centroiding routine

    output:

    x,y: positions
    fx, fy: FWHMs
    back: background
    peak: peak values


    """

    # centroid
    a = centroid.centroid_only(image.astype('<i4'), fwhm, thresh, boxsize, 2, sl, sh, rl, rh, 0)
    centroids = np.frombuffer(a, dtype='<f8')
    npoint = len(centroids)//7
    centroids = np.reshape(centroids, (npoint, 7))

    # extract useful values

    x = centroids[:, 0]
    y = centroids[:, 1]

    fx = (centroids[:, 2])
    fy = (centroids[:, 3])
    back = (centroids[:, 4])
    peak = (centroids[:, 5])

    return x, y, fx, fy, back, peak


def getAllCentroids(files, fwhm, boxsize, thresh, rl, rh, sl, sh, outfile):
    """

    wrapper to do the centroiding for a set of files. Writes to a text
    output file.  Returns the last set of xy positions (used for the
    next routine.

    input: 

    files: list of FITS files
    outfile: name of output text file
    FWHM, boxsize, thresh, rl,rh,sl,sh: parameters for centroiding routine

    output: writes to a file with each line of the format
       framenum, x, y, fx, fy, peak, back

    """

    # text file for output
    ff = open(outfile, "w")
    nfiles = len(files)

    # cycle through files
    filenum = 0

    print(str(nfiles)+" Frames. Centroiding ", end="")
    for file in files:
        #read in image
        image = visplot.getImage(file)
        # centroid
        x, y, fx, fy, back, peak = getCentroids(image, fwhm, boxsize, thresh, rl, rh, sl, sh)
        print(str(filenum+1)+", found "+str(len(x))+" centroids, ", end="")

        # write results to file, labelling the file number for bookkeeping
        for i in range(len(x)):
            print(filenum, x[i], y[i], fx[i], fy[i], back[i], peak[i], file=ff)
        filenum += 1

    print()
    ff.close()
    return x, y


def simpleDistortion(xx, yy, xs, ys, removeScale):
    """

    Calculate the distortion map using opencv library

    input: 

    xx,yy mask  coordinates
    xs,ys: spot coordinates

    returns:

    c,c1: distortion magnitude in mm and % of field
    pts1,pts2: points with average removed
    diffx,diffy: difference between pts2 and pts with translation/rotation removed

    """

    # subtract average to align

    nn = len(xs)
    pts1 = np.zeros((1, nn, 2))
    pts2 = np.zeros((1, nn, 2))
    xav1 = xs.mean()
    yav1 = ys.mean()
    xav2 = xx.mean()
    yav2 = yy.mean()

    for i in range(len(xs)):
        # pts1[0,i,0]=xs[i]-xav1
        # pts1[0,i,1]=ys[i]-yav1
        # pts2[0,i,0]=xx[i]-xav2
        # pts2[0,i,1]=yy[i]-yav2
        pts1[0, i, 0] = xs[i]
        pts1[0, i, 1] = ys[i]
        pts2[0, i, 0] = xx[i]
        pts2[0, i, 1] = yy[i]

    # cv2 needs float32

    pts1 = np.float32(pts1)
    pts2 = np.float32(pts2)

    # get the affine transformation - translation, rotation, scaling
    # uses opencv library
    transformation = cv2.estimateRigidTransform(pts1, pts2, False)

    visplot.checkMatched(pts1[0, :, 0], pts1[0, :, 1], pts2[0, :, 0], pts2[0, :, 1], "test", 1)

    transformFull = transformation.copy()

    # calculate the scaling from the result
    sx = np.sqrt(transformation[0, 0]**2+transformation[0, 1]**2)
    sy = np.sqrt(transformation[1, 0]**2+transformation[1, 1]**2)

    print("sx=", sx, " sy=", sy)

    if(removeScale == 1):
        # remove the scaling from the transformation matrix
        transformation[0, 0] /= sx
        transformation[0, 1] /= sx
        transformation[1, 0] /= sy
        transformation[1, 1] /= sy

    # calculate the transformed points
    points = cv2.transform(pts1, transformation)
    visplot.checkMatched(points[0, :, 0], points[0, :, 1], pts2[0, :, 0], pts2[0, :, 1], "test", 1)

    # difference between expected and actual point positions
    diffx = (pts2[0, :, 0]-points[0, :, 0])
    diffy = (pts2[0, :, 1]-points[0, :, 1])

    # in mm and in % of field

    c = np.sqrt(diffx*diffx+diffy*diffy)
    c1 = np.sqrt(diffx*diffx+diffy*diffy)/336.*100

    pts1[:, :, 0] = ma.masked_values(pts1[:, :, 0], 0)
    pts1[:, :, 1] = ma.masked_values(pts1[:, :, 1], 0)
    pts2[:, :, 0] = ma.masked_values(pts2[:, :, 0], 0)
    pts2[:, :, 1] = ma.masked_values(pts2[:, :, 1], 0)
    pts1 = ma.masked_invalid(pts1)
    pts2 = ma.masked_invalid(pts2)

    return c, c1, pts1, pts2, diffx, diffy, transformFull


def simpleDistortionBack(xx, yy, xs, ys, x1, y1):
    """

    Calculate the distortion map using opencv library

    input: 

    xx,yy mask  coordinates
    xs,ys: spot coordinates

    returns:

    c,c1: distortion magnitude in mm and % of field
    pts1,pts2: points with average removed
    diffx,diffy: difference between pts2 and pts with translation/rotation removed

    """

    # subtract average to align

    nn = len(xs)
    pts1 = np.zeros((1, nn, 2))
    pts2 = np.zeros((1, nn, 2))
    xav1 = xs.mean()
    yav1 = ys.mean()
    xav2 = xx.mean()
    yav2 = yy.mean()

    for i in range(len(xs)):
        # pts1[0,i,0]=xs[i]-xav1
        # pts1[0,i,1]=ys[i]-yav1
        # pts2[0,i,0]=xx[i]-xav2
        # pts2[0,i,1]=yy[i]-yav2
        pts1[0, i, 0] = xs[i]
        pts1[0, i, 1] = ys[i]
        pts2[0, i, 0] = xx[i]
        pts2[0, i, 1] = yy[i]

    # cv2 needs float32

    pts1 = np.float32(pts1)
    pts2 = np.float32(pts2)

    # get the affine transformation - translation, rotation, scaling
    # uses opencv library

    transformation = cv2.estimateRigidTransform(pts1, pts2, False)
    transformFull = transformation.copy()

    # calculate the scaling from the result
    sx = np.sqrt(transformation[0, 0]**2+transformation[0, 1]**2)
    sy = np.sqrt(transformation[1, 0]**2+transformation[1, 1]**2)

    print("sx=", sx, " sy=", sy)

    # remove the scaling from the transformation matrix
    transformation[0, 0] /= sx
    transformation[0, 1] /= sx
    transformation[1, 0] /= sy
    transformation[1, 1] /= sy

    # calculate the transformed points
    points = cv2.transform(pts1, transformation)

    # difference between expected and actual point positions
    diffx = (pts2[0, :, 0]-points[0, :, 0])
    diffy = (pts2[0, :, 1]-points[0, :, 1])

    # in mm and in % of field

    c = np.sqrt(diffx*diffx+diffy*diffy)
    c1 = np.sqrt(diffx*diffx+diffy*diffy)/336.*100

    pts1[:, :, 0] = ma.masked_values(pts1[:, :, 0], 0)
    pts1[:, :, 1] = ma.masked_values(pts1[:, :, 1], 0)
    pts2[:, :, 0] = ma.masked_values(pts2[:, :, 0], 0)
    pts2[:, :, 1] = ma.masked_values(pts2[:, :, 1], 0)
    pts1 = ma.masked_invalid(pts1)
    pts2 = ma.masked_invalid(pts2)

    return c, c1, pts1, pts2, diffx, diffy, transformFull


def totalFlux(xs, ys, image):
    """

    Very simple estimate of total flux of a point. Sums up the values in a box,
    subracting average background of image. 

    input: 

    xs,ys: coords of spots (pixel)
    image: image array

    returns: array of fluxes

    """

    ind = np.where(image > 0)

    # get average background
    back, a, b = sigmaclip(np.ravel(image[ind]))
    back = back.mean()

    # initialize array

    flux = np.zeros((len(xs)))

    # cycle through the points, do a very simple aperture photometry
    # summing values in box

    for i in range(len(xc)):
        # for i in range(2):
        xx = np.round(xc[i])
        yy = np.round(yc[i])
        j1 = int(xx)
        i1 = int(yy)

        for ii in range(-4, 5):
            for jj in range(-4, 5):
                j1 = int(xx+ii)
                i1 = int(yy+jj)
                flux[i] = flux[i]+image[i1, j1]-back

    return flux


def getRMSStats(xArray, yArray, fxArray, fyArray, peakArray, backArray):
    """

    calculate the RMS and average statistics for a set of seeing measurements.

    input: 
    xArray, yArray: arrays of registered positions
    fxArray, fyArray, peakArray, backArray: spot parameters to go with above

    output:
    xAv,yAv: average x and y positions 
    fxAv,fyAv: average FWHMS
    peakAv,backAv: average background and peak
    rms: rms in pixel positions, after x-y shift is removed
    nMatch: number of measurements that went into each spot

    """

    # get dimensions
    npoints, nframes = xArray.shape

    # get RMS and Averages
    xAv = xArray.mean(axis=1)
    yAv = yArray.mean(axis=1)
    fxAv = fxArray.mean(axis=1)
    fyAv = fyArray.mean(axis=1)
    peakAv = peakArray.mean(axis=1)
    backAv = backArray.mean(axis=1)
    nMatch = peakArray.count(axis=1)

    # now calculate and subtract the average of each frame (to remove)
    # frame-to-frame shifts

    xArray1 = np.zeros((npoints, nframes))
    yArray1 = np.zeros((npoints, nframes))

    for i in range(nframes):
        xFrameOff = (xArray[:, i]-xAv)
        yFrameOff = (yArray[:, i]-yAv)

        xArray1[:, i] = xFrameOff - xFrameOff.mean()
        yArray1[:, i] = yFrameOff - yFrameOff.mean()

        # print(xFrameOff.mean(),yFrameOff.mean())
        # xArray1[:,i]=xArray[:,i]-xFrameAv
        # yArray1[:,i]=yArray[:,i]-yFrameAv

    # get distance of change in each frame

    dd = np.sqrt((xArray1)**2+(yArray1)**2)

    # calculate RMS
    rms = dd.std(axis=1)

    return xAv, yAv, fxAv, fyAv, peakAv, backAv, rms, nMatch


def makeFakeCentroids(nframes, rms, avShift, avRot, avScale):
    """

    testing and demonstration routine to calculate a grid of simulated
    positions with a given RMS and  shift in position.

    input:
    nframes: number of frames to generate
    rms: desired RMS in position (in pixels)
    avShift: desired RMS of shift in position (in pixels)
    avScale: desired RMS of change in scale (in percent)
    avRot: RMS of rotation (in radians)

    output: returns an analogue to output of matchAllPoints, plus the
            transformed mask coordinates

    xArray,yArray: array of psotions
    fxArray,fyArray: array of FWHMs
    backArray,peakArray: array of background/peak values
    xm,ym: mask coordinates, transformed to image


    """

    # this approximates real pixel coordinates

    xd = 4439.53976392
    yd = 2892.62786802
    theta = 0
    scale = 14.723529411764705

    # get the mask coordinates in mm
    x, y = viscoords.maskinMM(1)
    x = x-168
    y = y-168
    npoints = len(x)

    # set up the arrays
    xArray = np.zeros((npoints, nframes))
    yArray = np.zeros((npoints, nframes))
    fxArray = np.zeros((npoints, nframes))
    fyArray = np.zeros((npoints, nframes))
    peakArray = np.zeros((npoints, nframes))
    backArray = np.zeros((npoints, nframes))

    for i in range(nframes):

        # add a small shift, scale and rotation to each frame
        xshift = np.random.normal(0, avShift, 1)
        yshift = np.random.normal(0, avShift, 1)
        scalevar = scale*np.random.normal(1, avScale-1, 1)
        rotvar = theta + np.random.normal(0, avRot, 1)

        print(xshift, yshift, scalevar, rotvar)

        # transform the mask to the data
        xm, ym = viscoords.transformPoints(x, y, xd+xshift, yd+yshift, rotvar, scalevar)

        # and Gaussian noise for seeing (and in other variables)
        xArray[:, i] = xm+np.random.normal(0, rms, npoints)
        yArray[:, i] = ym+np.random.normal(0, rms, npoints)
        fxArray[:, i] = np.random.normal(3.5, 1, npoints)
        fyArray[:, i] = np.random.normal(4.2, 1.5, npoints)
        peakArray[:, i] = np.random.normal(1000, 200, npoints)
        backArray[:, i] = np.random.normal(100, 10, npoints)

    xArray = ma.masked_where(peakArray == 0, xArray)
    yArray = ma.masked_where(peakArray == 0, yArray)
    fxArray = ma.masked_where(peakArray == 0, fxArray)
    fyArray = ma.masked_where(peakArray == 0, fyArray)
    backArray = ma.masked_where(peakArray == 0, backArray)
    peakArray = ma.masked_where(peakArray == 0, peakArray)

    return xArray, yArray, fxArray, fyArray, backArray, peakArray, x, y


def getTransByFrame(xArray, yArray, fxArray, fyArray, peakArray, xm, ym):
    """

    Estimate the affine transformation on a frame by frame basis, compared to the 
    pinhole mask. 

    input:
    xArray,yArray: array of psotions
    fxArray,fyArray: array of FWHMs
    backArray,peakArray: array of background/peak values
    xm,ym: mask coordinates

    output: affine transformation results
    xdAll,ydAll: translations
    sxAll,syAll: scaling
    rotAll: rotation
    fxFrameAv,fyFrameAv,peakFrameAv: average in spot characteristics

    """

    # set up variables
    npoints, nframes = xArray.shape
    xdAll = []
    ydAll = []
    sxAll = []
    syAll = []
    fxFrameAv = []
    fyFrameAv = []
    peakFrameAv = []
    rotAll = []

    allTrans = []

    for i in range(nframes):

        # use CV2 library to calculate the affine transformation

        transform, xd, yd, sx, sy, rotation = viscoords.getTransform(xArray[:, i], yArray[:, i], xm, ym, 1)
        xdAll.append(xd)
        ydAll.append(yd)
        sxAll.append(sx)
        syAll.append(sy)
        rotAll.append(rotation)

        # calculate the average values
        fxFrameAv.append(fxArray[:, i].mean())
        fyFrameAv.append(fyArray[:, i].mean())
        peakFrameAv.append(peakArray[:, i].mean())

        allTrans.append(transform)

    # convert data to numpy arrays
    xdAll = np.array(xdAll)
    ydAll = np.array(ydAll)
    fxFrameAv = np.array(fxFrameAv)
    fyFrameAv = np.array(fyFrameAv)
    peakFrameAv = np.array(peakFrameAv)
    sxAll = np.array(sxAll)
    syAll = np.array(syAll)
    rotAll = np.array(rotAll)

    return xdAll, ydAll, sxAll, syAll, rotAll, fxFrameAv, fyFrameAv, peakFrameAv, allTrans


def rotFunc(coords, xt, yt, rot):

    n = int(len(coords))
    x = coords[0:int(n/2)]
    y = coords[int(n/2):n]
    yval = (x-xt)*np.cos(rot)-(y-yt)*np.sin(rot)+xt
    xval = (x-xt)*np.sin(rot)+(y-yt)*np.cos(rot)+yt

    return np.array([xval, yval]).ravel()


def rotationParams(xs1, ys1, xs2, ys2, xc, yc, rot):

    params = np.array([xc, yc, rot])

    x = np.array([xs1, ys1]).ravel()
    y = np.array([xs2, ys2]).ravel()
    p0 = np.array([xc, yc, rot])

    popt, pconv = curve_fit(rotFunc, x, y, p0)

    print(popt)
    return popt


def matchRot(x1, y1, x2, y2, xc, yc, theta):

    # rotate second set of points to match first, based on guess of centre

    x_adj, y_adj = viscoords.transformGeneral(x1, y1, 0, 0, theta, 1, xc, yc)

    # plt.clf()
    # plt.plot(x_adj,y_adj,'dg')
    # plt.plot(x2,y2,'ob')
    # plt.title("Check")
    # plt.show()

    # x_adj=(x1-xc)*np.cos(theta)-(y1-yc)*np.sin(theta)+xc
    # y_adj=(x1-xc)*np.sin(theta)+(y1-yc)*np.cos(theta)+yc

    npoints1 = len(x1)
    npoints2 = len(x2)

    xa = np.zeros((len(x1)))
    ya = np.zeros((len(x1)))
    xb = np.zeros((len(x1)))
    yb = np.zeros((len(x1)))
    xx = np.zeros((len(x1)))
    yy = np.zeros((len(y1)))
    dval = np.zeros((len(x1)))

    for i in range(npoints1):

        # dd=100
        # ind=0

        dd = (x2-x_adj[i])**2+(y2-y_adj[i])**2
        ind = np.argmin(dd)
        dval[i] = dd[ind]

        # for j in range(npoints2):
        #
        #    nd=(x2[j]-x_adj[i])**2+(y2[j]-y_adj[i])**2
        #    if(nd < dd):
        #        dd=nd
        #        ind=j

        # dval[i]=dd
        xa[i] = x1[i]
        ya[i] = y1[i]
        xb[i] = x2[ind]
        yb[i] = y2[ind]
        xx[i] = x_adj[i]
        yy[i] = y_adj[i]

    plt.clf()
    plt.plot(xx, xb, 'dg')
    plt.plot(yy, yb, 'ob')
    plt.title("Check")
    plt.show()

    # mask outliers/unmatched points
    print(dval)
    back, a, b = sigmaclip(dval)
    mn = back.mean()
    st = back.std()
    print("stats", mn, st, max(dval))

    ind = np.where(abs(dval-mn) > 1*st)
    ma.masked_where(abs(dval-mn) > 1*st, xa)
    ma.masked_where(abs(dval-mn) > 1*st, xb)
    ma.masked_where(abs(dval-mn) > 1*st, ya)
    ma.masked_where(abs(dval-mn) > 1*st, yb)

    visplot.checkMatched(x2, y2, x_adj, y_adj, "aaa", 0)

    transformation, xd, yd, sx, sy, rotation = viscoords.getTransform(xa, ya, xb, yb, 1)

    # xc1,yc1,theta1=rotationParams(xa,ya,xb,yb,xc,yc,theta)
    return transformation, xd, yd, sx, sy, rotation


def rotatePoints(x, y, xc, yc, theta):

    x_rot = (x-xc)*np.cos(theta)-(y-yc)*np.sin(theta)+xc
    y_rot = (x-xc)*np.sin(theta)+(y-yc)*np.cos(theta)+yc

    return x_rot, y_rot
