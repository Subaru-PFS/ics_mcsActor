
"""

Routines related to coordinates and calculations

"""

import sys
import numpy as np
import cv2
import numpy.ma as ma
import itertools
from scipy import optimize
sys.path.append("/Users/karr/Science/PFS/NewCode/Code/pfs_utils/python/pfs/utils/coordinates/")
import CoordTransp

import windowedCentroid as centroid

    
def getCorners(x,y):

    ds,inds=getOrientation(x,y)
    ind=np.argsort(ds)

    x0=x[inds[ind[3]]]
    x1=x[inds[ind[2]]]
    y0=y[inds[ind[3]]]
    y1=y[inds[ind[2]]]

    return x0,x1,y0,y1
    
def getOrientation(xlast,ylast):

    
    xm=xlast.mean()
    ym=ylast.mean()

    #find the four 'corners' by distance from the mean point

    #divide into quadrands
    ind1 = np.where( np.logical_or( xlast-xm > 0, ylast-ym>0) )[0]
    ind2 = np.where( np.logical_or( xlast-xm > 0, ylast-ym<0) )[0]
    ind3 = np.where( np.logical_or( xlast-xm < 0, ylast-ym>0) )[0]
    ind4 = np.where( np.logical_or( xlast-xm < 0, ylast-ym<0) )[0]
    d1=np.sqrt((xlast[ind1]-xm)**2+(ylast[ind1]-ym)**2)
    d2=np.sqrt((xlast[ind2]-xm)**2+(ylast[ind2]-ym)**2)
    d3=np.sqrt((xlast[ind3]-xm)**2+(ylast[ind3]-ym)**2)
    d4=np.sqrt((xlast[ind4]-xm)**2+(ylast[ind4]-ym)**2)

    #distances for each
    d=np.sqrt((xlast-xm)**2+(ylast-ym)**2)
    d1=d.copy()
    d2=d.copy()
    d3=d.copy()
    d4=d.copy()

    #mask irrelevant points
    d1[ind1]=0
    d2[ind2]=0
    d3[ind3]=0
    d4[ind4]=0

    #max distance
    dm1=d1.max()
    dm2=d2.max()
    dm3=d3.max()
    dm4=d4.max()

    #index thereof
    ind1=d1.argmax()
    ind2=d2.argmax()
    ind3=d3.argmax()
    ind4=d4.argmax()

    #now find the two largest. These will be the good corners
    
    ds=np.array([dm1,dm2,dm3,dm4])
    inds=np.array([ind1,ind2,ind3,ind4])

    ind=np.argsort(ds)

    #and the positions
    x1=xlast[inds[ind[3]]]
    y1=ylast[inds[ind[3]]]
    x2=xlast[inds[ind[2]]]
    y2=ylast[inds[ind[2]]]

    return ds,inds
    #now get which orientation
    
    #if(x1-x2 < 1000):
    #    scale=168*2/(y1-y2)
    #else:
    #    scale=168*2/(x1-x2)
    #
    #if (np.logical_and( x1-xm > 0, y1-ym>0)):
    #    xtrans=x1-336*scale
    #    ytrans=y1
    #elif (np.logical_and( x1-xm > 0, y1-ym<0)):
    #    xtrans=x1-336*scale
    #    ytrans=y1-336*scale
    #elif (np.logical_and( x1-xm < 0, y1-ym<0)):
    #    xtrans=x1-336*scale
    #    ytrans=y1-336*scale


def pointLine(a,b,x,y,delta,sortby):

    """

    find points within a given distance of a line.

    input:

    a,b: coefficients of line (y=ax+b)
    x,y: set of coordinates
    delta: maximum dstance

    returns:

    x,y: subset of points within distance delta of line

    """

    
    d=abs(a*x-y+b)/np.sqrt(a*a+1)

    ind=np.where(d < delta)
    xline=x[ind]
    yline=y[ind]

    if(sortby=='y'):
        ind=(-yline).argsort()
    if(sortby=='x'):
        ind=(-xline).argsort()

    xline=xline[ind]
    yline=yline[ind]
    
    return xline,yline

def getLine(x1,x2,y1,y2):

    """

    Get the coefficients of a line between two points

    """
    
    #get the coefficients of a line between two points
    
    a=(y2-y1)/(x2-x1)
    b=y1-a*x1

    return a,b

def getPerpLine(a,b,x,y):

    """
    
    given the coefficients of a line (y=ax+b), and a point
    return the coefficents of a perpendicular line through the point.
    sorts the results by the given axis
    
    """
    
    a1=-1/a
    b1=y-a1*x
    
    return a1,b1

def polyfit2d(x, y, z, order=3):

    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

    
def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros(len(x))
    for a, (i,j) in zip(m, ij):
        #print(a,i,j)
        z += a * x**i * y**j
    return z

    
def rotAngle(x0,x1,y0,y1):

    angle=-(np.arctan2(x0-x1,y0-y1)+np.pi/2)

    return angle

def fakeFiducial(x,y):

    scale=23/2
    
    #get corner positions
    x0,x1,y0,y1=getCorners(x,y)

    angle=rotAngle(x0,x1,y0,y1)
    xmean=(x1+x0)/2.-scale*8
    ymean=(y1+y0)/2.

    #get the line through the top
    a1,b1=getLine(x0,x1,y0,y1)

    #perpendicular line
    a,b=getPerpLine(a1,b1,xmean,ymean)
    
    #points close to the line
    xline1,yline1=pointLine(a,b,x,y,20,'y')
    xline2,yline2=pointLine(a1,b1,x,y,20,'x')

    a2=a1
    b2=b1+yline1[27]-ymean
    xline3,yline3=pointLine(a2,b2,x,y,30,'x')

    xline=np.concatenate((xline1,xline2,xline3))
    yline=np.concatenate((yline1,yline2,yline3))

    
    #xline=xline1
    #yline=yline1

    return xline,yline,angle
    
def fakeFiducialTrue(offx,offy,angle):

    nn=43

    xtrue1=np.zeros((nn))-8
    ytrue1=(np.arange(nn-1,-1,-1))*8-336/2

    xtrue2=np.arange(nn-1,-1,-1)*8-336/2
    ytrue2=(np.zeros((nn)))+336-336/2

    xtrue3=np.arange(nn-1,-1,-1)*8-336/2
    ytrue3=(np.zeros((nn))+(16)*8.)-336/2

    xtrue=np.concatenate((xtrue1,xtrue2,xtrue3))
    ytrue=np.concatenate((ytrue1,ytrue2,ytrue3))

    xtrue=xtrue+offx
    ytrue=ytrue+offy
     
    #xtrue=xtrue1+offx
    #ytrue=ytrue1+offy

    xt=np.cos(angle)*xtrue-np.sin(angle)*ytrue
    yt=np.sin(angle)*xtrue+np.cos(angle)*ytrue

    return xt,yt

def xytoR(x,y,xc,yc):

    r=np.sqrt((x-xc)**2+(y-yc)**2)
    theta=np.arctan2(y-yc,x-xc)
    
    return r,theta

def rtoXY(r,theta,xc,yc):

    x=r*np.cos(theta)+xc
    y=r*np.sin(theta)+yc
    
    return x,y
            
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


def matchAllPoints(centroids,xx,yy,tol):

    """

    takes a batch of centroids created by getAllCentroids, registers
    the positions for each frame to the mask, and returns a set of arrays for
    each variable, in which points that were not detected have been masked. 

    xx, yy are the mask hole positions, and should be transformed to match the
    image

    input: 
    centroidFile: file with centroid output
    xx,yy: transformed mask coordinates

    output: 
    xArray,yArray: array of registered positions
    fxArray,fyArray,backArray,peakArray: registered parameters for spots

    """

    #create arrays
    #centroids=np.loadtxt(centroidFile)

    #size of output array
    npoints=len(xx)
    nfiles=np.int(centroids[:,0].max())
    
    xArray=np.zeros((npoints,nfiles))
    yArray=np.zeros((npoints,nfiles))
    fxArray=np.zeros((npoints,nfiles))
    fyArray=np.zeros((npoints,nfiles))
    peakArray=np.zeros((npoints,nfiles))
    backArray=np.zeros((npoints,nfiles))
    qualArray=np.zeros((npoints,nfiles))

    #load the centroid data
    
    print(str(nfiles)+" frames. Matching ",end="")
    #for each image
    for i in range(nfiles):

        print(str(i+1)+", ",end="")

        #get spot positions, etc from a particular image

        ind=np.where(centroids[:,0]==i)

        ll=centroids[ind,0].shape[1]
        x=centroids[ind,1].reshape(ll)
        y=centroids[ind,2].reshape(ll)
        fx=centroids[ind,3].reshape(ll)
        fy=centroids[ind,4].reshape(ll)
        peak=centroids[ind,5].reshape(ll)
        back=centroids[ind,6].reshape(ll)
        qual=centroids[ind,7].reshape(ll)

        #nearest neighbour matching
        for j in range(npoints):

            dd=np.sqrt((xx[j]-x)**2+(yy[j]-y)**2)
            ind=np.where(dd==dd.min())

            #need to filter in case points are missing
            if(dd.min() < tol):
                xArray[j,i]=x[ind]
                yArray[j,i]=y[ind]
                fxArray[j,i]=fx[ind]
                fyArray[j,i]=fy[ind]
                peakArray[j,i]=peak[ind]
                backArray[j,i]=back[ind]
                qualArray[j,i]=qual[ind]
    print()
    #mask unfound values

    
    xArray=ma.masked_where(xArray <= 0 ,xArray)
    yArray=ma.masked_where(xArray <= 0 ,yArray)
    fxArray=ma.masked_where(xArray <= 0 ,fxArray)
    fyArray=ma.masked_where(xArray <= 0 ,fyArray)
    backArray=ma.masked_where(xArray <= 0 ,backArray)
    peakArray=ma.masked_where(xArray <= 0 ,peakArray)
    qualArray=ma.masked_where(xArray <= 0 ,qualArray)
    
    return xArray,yArray,fxArray,fyArray,backArray,peakArray,qualArray

    
def matchPoints(centroids,xp,yp,tol):

    """

    Do nearest neighbour matching on a set of centroids

    Input: 
       centroids: centroid array in the form
           x,y,other parameters.....
       xp,yp: array of points to match
       tol: maximum separation for match
 
    """

    
    n=len(xp)
    nparm=centroids.shape[1]
    
    outArray=np.zeros((n,nparm))
        
    for i in range(n):
        dd=dist(xp[i],y[i],x,y)
        ind=np.argmin(dd)

        #found an acceptable match
        if(dd[ind[0]] < tol):

            outArray[i,:]=centroids[ind[0],:]

        #otherwise set values to NaN
        else:
            outArray[i,:]=np.empty((nparms)).fill(np.nan)

    #convert to masked array, masking unfound points
    outArray=ma.masked_where(outArray==np.nan)

    return outArray

def transformPoints(rCenter,offset,zAngle,inRot,fine):

    xx,yy=maskinMM(fine)
    #xx,yy=readPositions()

    xx=xx+xoff
    yy=yy+yoff

    xylist=[xx,yy]
    
    xyout=CoordTranp.CoordinateTransform(xylist,zAngle,"pfi_mcs_wofe",cent=rCenter,inr=inRot)

    return xyout[:,0],xyout[:,1]
    
def matchPoints(centroids,xp,yp,tol):

    """

    Do nearest neighbour matching on a set of centroids

    Input: 
       centroids: centroid array in the form
           x,y,other parameters.....
       xp,yp: array of points to match
       tol: maximum separation for match
 
    """

    
    n=len(xp)
    nparm=centroids.shape[1]
    
    outArray=np.zeros((n,nparm))
        
    for i in range(n):
        dd=dist(xp[i],y[i],x,y)
        ind=np.argmin(dd)

        #found an acceptable match
        if(dd[ind[0]] < tol):

            outArray[i,:]=centroids[ind[0],:]

        #otherwise set values to NaN
        else:
            outArray[i,:]=np.empty((nparms)).fill(np.nan)

    #convert to masked array, masking unfound points
    outArray=ma.masked_where(outArray==np.nan)

    return outArray

def transformMask(rCenter,offset,zAngle,inRot,fine):

    xx,yy=maskinMM(fine)
    #xx,yy=readPositions()

    xx=xx+xoff
    yy=yy+yoff

    xylist=[xx,yy]
    
    xyout=CoordTranp.CoordinateTransform(xylist,zAngle,"pfi_mcs_wofe",cent=rCenter,inr=inRot)

    return xyout[:,0],xyout[:,1]

    
def getRMSStats(xArray,yArray,fxArray,fyArray,peakArray,backArray,xdAll,ydAll,sxAll,syAll,rotAll,xm,ym):

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
    

    #get dimensions
    npoints,nframes=xArray.shape
    
    #get RMS and Averages
    #xAv=xArray.mean(axis=1)
    #yAv=yArray.mean(axis=1)
    fxAv=fxArray.mean(axis=1)
    fyAv=fyArray.mean(axis=1)
    peakAv=peakArray.mean(axis=1)
    backAv=backArray.mean(axis=1)
    nMatch=xArray.count(axis=1)

    #now calculate and subtract the best transform between the mask and points.

    xArray1=np.zeros((npoints,nframes))
    yArray1=np.zeros((npoints,nframes))
    
    for i in range(nframes):
        xArray1[:,i], yArray1[:,i] = transformPointsNew(xArray[:,i],yArray[:,i],xdAll[i],ydAll[i],rotAll[i],sxAll[i],syAll[i])

    #get distance of change in each frame
    #xAv=xArray1.mean(axis=1)
    #yAv=yArray1.mean(axis=1)
    xArray1=ma.masked_where((xArray1 < 100) | (xArray.mask == True),xArray1)
    yArray1=ma.masked_where((yArray1 < 100) | (xArray.mask == True),yArray1)

    dd=np.zeros(xArray.shape)
    xd=np.zeros(xArray.shape)
    yd=np.zeros(xArray.shape)
    for i in range(xArray.shape[0]):
        dd[i,:]=np.sqrt((xArray1[i,:]-xm[i])**2+(yArray1[i,:]-ym[i])**2)
        xd[i,:]=np.sqrt((xArray1[i,:]-xm[i])**2)
        yd[i,:]=np.sqrt((yArray1[i,:]-ym[i])**2)
        
    #calculate RMS
    
    dd=ma.masked_where(((dd <=0) | (xArray.mask == True)), dd)
    xd=ma.masked_where(((dd <=0) | (xArray.mask == True)), xd)
    yd=ma.masked_where(((dd <=0) | (xArray.mask == True)), yd)
    xAv=xArray1.mean(axis=1)
    yAv=yArray1.mean(axis=1)
    rms=dd.std(axis=1)
    rmsX=xd.std(axis=1)
    rmsY=yd.std(axis=1)
    

    return xAv,yAv,fxAv,fyAv,peakAv,backAv,rms,nMatch,xArray1,yArray1,dd,rmsX,rmsY,xd,yd

def getTransByFrame(xArray,yArray,fxArray,fyArray,peakArray,xm,ym):

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

    #set up variables
    npoints,nframes=xArray.shape
    xdAll=[]
    ydAll=[]
    sxAll=[]
    syAll=[]
    fxFrameAv=[]
    fyFrameAv=[]
    peakFrameAv=[]
    rotAll=[]

    allTrans=[]

    print("Translating ",end="")

    for i in range(nframes):
        print(i+1,', ',end="")
        #use CV2 library to calculate the affine transformation

        transform,xd,yd,sx,sy,rotation=getTransform(xArray[:,i],yArray[:,i],xm,ym,1)
        xdAll.append(xd)
        ydAll.append(yd)
        sxAll.append(sx)
        syAll.append(sy)
        rotAll.append(rotation)

        #calculate the average values
        fxFrameAv.append(fxArray[:,i].mean())
        fyFrameAv.append(fyArray[:,i].mean())
        peakFrameAv.append(peakArray[:,i].mean())

        allTrans.append(transform)

    print()
    #convert data to numpy arrays
    xdAll=np.array(xdAll)
    ydAll=np.array(ydAll)
    fxFrameAv=np.array(fxFrameAv)
    fyFrameAv=np.array(fyFrameAv)
    peakFrameAv=np.array(peakFrameAv)
    sxAll=np.array(sxAll)
    syAll=np.array(syAll)
    rotAll=np.array(rotAll)
    
    return xdAll,ydAll,sxAll,syAll,rotAll,fxFrameAv,fyFrameAv,peakFrameAv,allTrans

def getTransform(x,y,xx,yy,getVals):

    """

    given two sets of registered points, estimate the rigid transformation
    this is in a separate routine mostly for bookkeeping purposes. 
    Returns transformation matrix, and if getVals == 1 returns the 
    extracted parameters (rotation, translation, scale) as well. 

    input:
    x,y: input positions
    xx,yy: transformed positions
    getVales: if ==1, return parameters too

    output: 
    transformation: matrix 
    xd,yd: translations
    sx,sy: scalings
    rotation: rotation (radians)

    """

    #turn data into right form
    pts1=np.zeros((1,len(x),2))
    pts2=np.zeros((1,len(x),2))

    pts1[0,:,0]=x
    pts1[0,:,1]=y

    pts2[0,:,0]=xx
    pts2[0,:,1]=yy

    #float32 is needed
    pts1=np.float32(pts1)
    pts2=np.float32(pts2)

    #calculate the transformation
    transformation = cv2.estimateRigidTransform(pts1, pts2, False)
    if(getVals == 0):
        return transformation
    
    if(getVals == 1):

        #extract the parameters

        sx=np.sqrt(transformation[0,0]**2+transformation[0,1]**2)
        sy=np.sqrt(transformation[1,0]**2+transformation[1,1]**2)

        xd=transformation[0,2]
        yd=transformation[1,2]

        rotation = np.arctan2(transformation[1,0]/sx,transformation[1,1]/sy)
        
        return transformation,xd,yd,sx,sy,rotation



def getFieldDefinition(fieldID):

    if(fieldID=='oct18'):
        fiducials=np.loadtxt('fiducials.dat',delimiter=',')
        scienceFibres=np.loadtxt('scienceFibres.dat',delimiter=',')

    return fiducials,scienceFibres

def cobrasToMCS(fiducials,scienceFibres,za,inr,rotCent,offset):

    xyin=np.array([np.concatenate([fiducials[:,1],scienceFibres[:,1]]),np.concatenate([fiducials[:,2],scienceFibres[:,2]])])

    xyin[0,:]
    xyin[1,:]

    xyout=CoordTransp.CoordinateTransform(xyin,za,'pfi_mcs_wofe',inr=-inr*180/np.pi,cent=rotCent)
    
    return xyout[0,:],xyout[1,:],xyout[2,:]

def transformPointsNew(x,y,xd,yd,theta,sx,sy):

    """
    Apply a rigid transformation to the mask (trans, scale, rot). Mostly bookkeeping
    stuff. 

    input:
    x,y: mask positions
    xd,yd: translation
    theta: rotation (radians)
    s: scale

    output: transformed x,y

    """
    
    #create transformation matrix
    matrix=np.zeros((2,3))
    matrix[0,0]=np.cos(theta)*sx
    matrix[0,1]=-np.sin(theta)*sy
    matrix[1,0]=np.sin(theta)*sx
    matrix[1,1]=np.cos(theta)*sy
    matrix[0,2]=xd
    matrix[1,2]=yd

    #bookkeeping for coordinate format
    pts=np.zeros((1,len(x),2))

    pts[0,:,0]=x
    pts[0,:,1]=y

    pts=np.float32(pts)

    #do the transform
    pts1=cv2.transform(pts,matrix)

    #more bookkeeping
    xx=pts1[0,:,0]
    yy=pts1[0,:,1]

    return xx,yy


