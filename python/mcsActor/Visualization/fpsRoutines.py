
import sys
import numpy as np
import cv2
import numpy.ma as ma
import itertools
from scipy import optimize
sys.path.append("/Users/karr/Science/PFS/NewCode/Code/pfs_utils/python/pfs/utils/coordinates/")
import CoordTransp



def getCorners(x,y):

    """
    
    Routine to find the two square corners of the set of pinhole mask spots. This is used when
    Doing the rotation calculation, and the backup routine for finding fiducial fibres. 

    Input: x,y - positions

    Returns: x0,x1,y0,y1 - x and y positions of points

    Note that this needs a fairly clean set of centroids to work. 

    """

    #call the geometry routine
    
    ds,inds=getOrientation(x,y)

    #sort by distance from centre and assign
    ind=np.argsort(ds)

    x0=x[inds[ind[3]]]
    x1=x[inds[ind[2]]]
    y0=y[inds[ind[3]]]
    y1=y[inds[ind[2]]]

    return x0,x1,y0,y1
    
def getOrientation(xlast,ylast):

    """Routine to fine the corners of a set of pinhole mask spots. The
    definition of 'corner' is the point in each quadrant that is
    farthest from the mean position. The square corners give the
    larger distances; the asymmetry of the mask breaks the degenerecy,
    allowing us to select a specific corner. 

    Input: x,y - positions

    Returns: arrays of indices and distances for the four cornesrs

    """
    
    #first, find the mean position. 
    xm=xlast.mean()
    ym=ylast.mean()

    #find the four 'corners' by distance from the mean point

    #indices for quadrants (all points not in a quadrant)
    ind1 = np.where( np.logical_or( xlast-xm > 0, ylast-ym>0) )[0]
    ind2 = np.where( np.logical_or( xlast-xm > 0, ylast-ym<0) )[0]
    ind3 = np.where( np.logical_or( xlast-xm < 0, ylast-ym>0) )[0]
    ind4 = np.where( np.logical_or( xlast-xm < 0, ylast-ym<0) )[0]

    #distances for each
    d=np.sqrt((xlast-xm)**2+(ylast-ym)**2)
    d1=d.copy()
    d2=d.copy()
    d3=d.copy()
    d4=d.copy()

    #set irrelevant points to zero
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
 
    return ds,inds

def getFieldDefinition(fieldID):

    """

    program to read in the set of cobras. Currently a dummy program
    reading from file. 

    """

    fiducials=np.loadtxt('fiducials.dat',delimiter=',')
    scienceFibres=np.loadtxt('scienceFibres.dat',delimiter=',')

    return fiducials,scienceFibres

def getFibrePos(fiducials,scienceFibres,za,inr,rotCent,offset):

    """
    
    Convert a set of fiducials and science fibres (in mask coordinates)
    into expected pixel positons on the image. 


    """

    #rotation adjustment
    inr=inr-180
    if(inr < 0):
        inr=inr+360
        
    #concatenate - MCS doesn't care which is fiducial and which is science

    xx=np.array([np.concatenate([fiducials[:,1],scienceFibres[:,1]])]).ravel()
    yy=np.array([np.concatenate([fiducials[:,2],scienceFibres[:,2]])]).ravel()

    #now offset to centre of rotation
    xx-=offset[0]
    yy-=offset[1]

    #correect input format
    xyin=np.array([xx,yy])

    #call the routine
    xyout=CoordTransp.CoordinateTransform(xyin,za,'pfi_mcs_wofe',inr=inr,cent=rotCent)
    
    return np.array([np.arange(len(xx)),xyout[0,:],xyout[1,:]]).T

def getDiff(centroids,fibrePos):

    dx=centroids[:,1]-fibrePos[:,1]
    dy=centroids[:,2]-fibrePos[:,2]

    return dx,dy
