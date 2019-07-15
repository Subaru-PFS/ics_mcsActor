

import numpy as np
from astropy.io.fits import getdata
import matplotlib.pylab as plt
from scipy.stats import sigmaclip
import numpy.ma as ma



def dist(x1,y1,x2,y2):

    """Simple distance between two points (or point and array)"""
    
    dd=np.sqrt((x1-x2)**2+(y1-y2)**2)
    return dd


def findHomes(centroids,fibrePos,tol):

    """

    Do nearest neighbour matching on a set of centroids (ie, home position case).
    The routine is flexible with the number of paramters in the centroid array - they 
    will all be copied to the output array. 

    Input: 
       centroids: centroid array in the form
           id,x,y,other parameters.....
       fibrePos: expected positions, (Nx3) array with first column as index
       tol: maximum separation for match, in pixels

    Returns: sorted array with index in the first column

 
    """

    #get shape of output array
    n=len(fibrePos[:,0])
    nparm=centroids.shape[1]

    #create an array to hold the output
    outArray=np.zeros((n,nparm))

    #now do the matching
    for i in range(n):

        dd=dist(fibrePos[i,1],fibrePos[i,2],centroids[:,1],centroids[:,2])
        ind=np.argmin(dd)
        #found an acceptable match
        if(dd[ind] < tol):
            outArray[i,:]=centroids[ind,:]

        #otherwise set values to NaN
        else:
            outArray[i,:]=np.empty((nparm)).fill(np.nan)

        #and assign the index
        outArray[i,0]=fibrePos[i,0]

    #convert to masked array, masking unfound points
    outArray=ma.masked_where(outArray!=outArray,outArray)

    return outArray

def getThresh(image,threshMethod,threshSigma,threshFact,findSigma,centSigma,fibrePos=None):

    """

    Wrapper routine for Threshold calcualaton, toggling beteen calibrarion adn operation mode.

    input: 

    image: image
    threshMethod: =calib for calibration mode (must be rotated by multple of 90 degrees)
                  =fieldID for real operation

    threshSigma: threshold for sigma clipping
    threshFact: a fudge factor for the calib mode which takes into account summer over different sized
                axes for a square array. 
    findSigma, centSigma: multiple of RMS for finding and centroiding spots
    fibrePos = fibrePos array for fieldID, or None for calib

    Returns:
        thresh1,thresh2: the two thresholds
        xrange,yrange: outer bounds of region used for calculation. This is a circle for fieldID, 
        and a square for calib. 

    """

    #toggle between two methods
    
    if(threshMethod=='fieldID'):

        a=getAutoThresh(image,threshSigma,fibrePos)

        #widest points of a circle
        xrange=[fibrePos[:,1].mean()-1500,fibrePos[:,1].mean()+1500]
        yrange=[fibrePos[:,2].mean()-1500,fibrePos[:,2].mean()+1500]

    elif(threshMethod=='calib'):

        xrange,yrange=getRegion(image,threshSigma,threshFact)

        a=getManualThresh(image,xrange,yrange,threshSigma)

    #get the thresholds
    thresh1=a.mean()+a.std()*findSigma
    thresh2=a.mean()+a.std()*centSigma

    return thresh1,thresh2,xrange,yrange

def getManualThresh(image,xrange,yrange,sigmaThresh):

    """

    returns a sigma clipped array. 

    Input: 
       image: image
       xrange,yrange: bounds of region 
       sigmaThresh: threshold for sigma clipping

    Returns:
       clipped array

    """

    #sigma clip the image to get rid of the spots
    subIm=image[xrange[0]:xrange[1],yrange[0]:yrange[1]]
    a,b,c=sigmaclip(subIm.ravel(),sigmaThresh,sigmaThresh)

    return a


def getRegion(image,high,factor):

    """

    A function to manually find the region of the pinhole mask when it's not known in advance. 

    Input: 
        Image: image
        high: factor for sigma clipping
        factor: fudge factor to take into account that 

    """

    #first a sigmaclip
    im,a,b=sigmaclip(image,high=high)

    #find the boundaries of the region
    xlow,xhigh,ylow,yhigh=getBoundary(image,a,b,0)

    #and a second iteration
    im,a,b=sigmaclip(image[xlow:xhigh,ylow:yhigh],high=high)
    rrms=im.std()/factor
    xlow,xhigh,ylow,yhigh=getBoundary(image,a,b,rrms)

    return [xlow,xhigh],[ylow,yhigh]

def getBoundary(image,a,b,rms):

    """

    find the region that is occupied by the pinhole mask. 

    input:
 
        image
        a,b - upper and lower limits for sigma clipping
        rms - adjustment factor for short axis of image (default zero)

    output:

        x1,y1,x2,y2: limits of box for the pinhole mask structure. 

    """

    #set up the variables
    sz=image.shape

    prof1=np.zeros((sz[0]))
    prof0=np.zeros((sz[1]))

    
    #do a sigma clipped summing collapsing along each axis.
    #there's probabably a fancy python way to make this more
    #efficient
    
    for i in range(sz[0]):
        ind=np.where(image[i,:] < b)
        prof1[i]=np.mean(image[i,ind])
    for i in range(sz[1]):
        ind=np.where(image[:,i] < b)
        prof0[i]=np.mean(image[ind,i])

        
    #get the mean of each profile, with an adjustment for the short axis
    #of the image
    
    pm0=prof0.mean()
    pm1=prof1.mean()-rms

    #next step - move along the summed profile and find the step function
    #portion by looking for the place in which the step crossed the (adjusted)
    #mean value. 

    #start from the left and move right, then start from the right and move left.
    #if we reach the middle of the image without finding it, then we assume
    #the region is right up against the edge of the image. 
    
    ###########################################
    
    found=0
    i=0
    while((found==0) & (i < sz[0]/2)):
        if((prof1[i]-pm1)*(prof1[i+1]-pm1)<0):
            found=1
        i=i+1
        
    x1=i
        
    ###########################################3
    i=int(sz[0]-1)
    found=0
    while((found==0) & (i > sz[0]/2)):
        if((prof1[i]-pm1)*(prof1[i-1]-pm1)<0):
            found=1
        i=i-1

    if(found==1):
        x2=i
    else:
        x2=sz[0]-1
        
    ###########################################3

    #start at 1000 rather than zero in the long axis
    #to avoid variation in background flux at the edge of the image
    
    found=0
    i=1000
    
    while((found==0) & (i < sz[1]/2)):
        if((prof0[i]-pm0)*(prof0[i+1]-pm0)<0):
            found=1
        i=i+1
        
    y1=i
    
    ###########################################3
    i=sz[1]-1
    found=0
    while((found==0) & (i > 1)):
        if((prof0[i]-pm1)*(prof0[i-1]-pm0)<0):
            found=1
        i=i-1
        
    y2=i

    ###########################################3

    return x1,x2,y1,y2

def getAutoThresh(image,threshSigma,fibrePos):

       
    #approximate centre of image
    mpx=fibrePos[:,1].mean()
    mpy=fibrePos[:,2].mean()

    #find points in a circle near the centre
    xx,yy=np.meshgrid(np.arange(image.shape[0]),np.arange(image.shape[1]))
    
    dfromc=np.sqrt((xx-mpx)**2+(yy-mpy)**2)
    
    ind=np.where(dfromc < 1500)

    #sigma clip the 
    a,b,c=sigmaclip(image[ind],threshSigma,threshSigma)

    return a

            
def calc_R(x,y, xc, yc):
    """
    calculate the distance of each 2D points from the center (xc, yc)
    """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f(c, x, y):
    """
    calculate the algebraic distance between the data points
    and the mean circle centered at c=(xc, yc)
    """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()

def least_squares_circle(x,y):
    """
    Circle fit using least-squares solver.
    Inputs:

        - coords, list or numpy array with len>2 of the form:
        [
    [x_coord, y_coord],
    ...,
    [x_coord, y_coord]
    ]

    Outputs:

        - xc : x-coordinate of solution center (float)
        - yc : y-coordinate of solution center (float)
        - R : Radius of solution (float)
        - residu : MSE of solution against training data (float)
    """
    # coordinates of the barycenter

    #x = np.array([x[0] for x in coords])
    #y = np.array([x[1] for x in coords])
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x,y))
    xc, yc = center
    Ri       = calc_R(x, y, *center)
    R        = Ri.mean()
    residu   = np.sum((Ri - R)**2)
    return xc, yc, R, residu

def flatField(image,flat):

    return image/flat



    
