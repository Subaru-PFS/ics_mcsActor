import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib
from importlib import reload
import sys
from pfs.utils import coordinates
import pickle

from pfs.utils import coordinates
from pfs.utils.coordinates import CoordTransp

import pandas as pd
from scipy.stats import sigmaclip
import centroid as centroid
#sys.path.insert(1, "/Users/karr/Science/PFS/pfsPackages/ics_cobraCharmer/python/ics/cobraCharmer/")
#sys.path.insert(1, "/Users/karr/software/mhs/products/DarwinX86/spt_operational_database/0.0.6/python/opdb-0.1-py3.7.egg/")
rootPath=os.path.join(os.environ['ICS_MHS_ROOT'])
dbPath=os.path.join(rootPath,"products/Linux64//spt_operational_database/0.0.6/python/opdb-0.1-py3.8.egg/")
sys.path.insert(1, dbPaths)
from opdb import models,utils,opdb
import pfi as pfi 
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import cv2
import copy

def readCobraGeometryFake():

    """
    creates fake cobra geometry for pinhole mask images, fo testing purposes. The centre is the mask position,
    the arm length is set to 4mm, the dot positions are offset and arbitrary. All "cobras" are good.
    """

    centrePos=np.loadtxt("/Users/karr/ics_mcsActor/python/mcsActor/Visualization/scienceFibres.dat",delimiter=",")
    centrePos[:,0]+=1
    armLength=np.repeat(3,centrePos.shape[0])
    nCobras=len(armLength)

    dotPos=np.zeros((centrePos.shape[0],4))
    dotPos[:,0]=centrePos[:,0]
    dotPos[:,1]=centrePos[:,1]+1.19
    dotPos[:,2]=centrePos[:,2]+1.19
    dotPos[:,3]=1.5

    goodIdx=np.copy(centrePos[:,0])-1
    return centrePos,armLength,dotPos,goodIdx

def readCobraGeometry(xmlFile,dotFile):

    """
    read cobra geometry from configuration file/inst_config

    dot positions from CSVfile at the moment, this will change

    The results will be return in whatever unit the input XML file is in
    """

    #geometry XML file

    pfic=pfi.PFI(fpgaHost='localhost',doConnect=False,logDir=None)
    aa=pfic.loadModel([pathlib.Path(xmlFile)])
    
    #first figure out the good cobras (bad positions are set to 0)
    centersAll=pfic.calibModel.centers
    goodIdx=np.array(np.where(centersAll.real != 0)).astype('int').ravel()
    
    #then extract the parameters for good fibres only
    centrePos=np.array([goodIdx+1,centersAll[goodIdx].real,centersAll[goodIdx].imag]).T
    armLength=(pfic.calibModel.L1[goodIdx]+pfic.calibModel.L2[goodIdx])
    
    #number of cobras
    nCobras=len(armLength)
#PFS_INSTDATA_DIR
    #at the moment read dots from CSV file
    dotData=pd.read_csv(dotFile,delimiter=",")
    dotPos=np.zeros((len(goodIdx),4))

    dotPos[:,0]=goodIdx+1
    dotPos[:,1]=dotData['x_tran'].values[goodIdx]
    dotPos[:,2]=dotData['y_tran'].values[goodIdx]
    dotPos[:,3]=dotData['r_tran'].values[goodIdx]

    with open('dumpGeom.pkl', 'wb') as output:
        pickle.dump([centrePos,armLength,dotPos,goodIdx],output)

    return centrePos,armLength,dotPos,goodIdx

def transformToMM(posPix,rotCent,offset,zenithAngle,insRot,fieldElement,pixScale=0):
    
    """
    transform coordinates from MCS pixels to MM distortion only.
    This is a wrapper for the pfs_utils.CoordTransp routine

    if pixScale !=0 a simple scaling is applied instead (for testing)

    Input: 
        posPix: input positions in nx3 shape, first column is an id number, next two coordinates
        rotCent: boresight, 2 element array
        offset: offset from centre of instrument
        zenithAngle, insRot: from FIRST header
        pixScale: set non zero to do a simple scaling instead

    Returns
       posMM: coorinates in pixels, same format as input

    """
 
    #for testing purposes assumes a straight pixel scaling
    if(pixScale!=0):
        return np.array([posPix[:,0],posPix[:,1]*pixScale,posPix[:,2]*pixScale]).T

    #rotation adjustment
    insRot=insRot-180
    if(insRot < 0):
        insRot=insRot+360

    xyout=CoordTransp.CoordinateTransform(posPix[:,1:3],zenithAngle,'mcs_pfi',inr=insRot,cent=rotCent)
    xyout[0,:]+=offset[0]
    xyout[1,:]+=offset[1]
    return np.array([posPix[:,0],xyout[0,:],xyout[1,:]]).T
     
def transformToPix(posMM,rotCent,offset,zenithAngle,insRot,fieldElement,pixScale=0):

    """
    
    transform coordinates form MM to MCS pixels distortion only
    This is a wrapper for the pfs_utils.CoordTransp routine


    If pixScale is non zero, a simple scaling factor is applied

    Input: 
        posMM: input positions in nx3 shape, first column is an id number, next two coordinates
        rotCent: boresight, 2 element array
        offset: offset from centre of instrument
        zenithAngle, insRot: from FIRST header
        pixScale: set non zero to do a simple scaling instead

    Returns
       posPix: coorinates in pixels, same format as input

    """
    
    #for testing purposes assumes a straight pixel scaling
    if(pixScale != 0):
        return np.array([posMM[:,0],posMM[:,1]/pixScale,posMM[:,2]/pixScale]).T

    if(fieldElement==True):
        fElem='pfi_mcs'
    else:
        fElem='pfi_mcs_wofe'
    #rotation adjustment
    insRot=insRot-180
    if(insRot < 0):
        insRot=insRot+360
        
    xyin=np.array([posMM[:,1]-offset[0],posMM[:,2]-offset[1]])
    #call the routine
    xyout=CoordTransp.CoordinateTransform(xyin,zenithAngle,fElem,inr=insRot,cent=rotCent).T

    return np.array([posMM[:,0],xyout[:,0],xyout[:,1]]).T


def calcAffineTransform(pos1,pos2):

    """

    given two sets of registered points, estimate the rigid transformation
    this is a wrapper for the cv2 routine
    Returns transformation matrix and extracted parameters (rotation, translation, scale)

    input:
    pos1 input positions in nx3 shape, first column is an id number, next two coordinates
    xx,yy: transformed positions
    getVales: if ==1, return parameters too

    output: 
    transformation: matrix 
    xd,yd: translations
    sx,sy: scalings
    rotation: rotation (radians)

    """

    sz=pos1.shape[0]
    #turn data into right form
    pts1=np.zeros((1,sz,2))
    pts2=np.zeros((1,sz,2))

    pts1[0,:,0]=pos1[:,1]
    pts1[0,:,1]=pos1[:,2]

    pts2[0,:,0]=pos2[:,1]
    pts2[0,:,1]=pos2[:,2]

    #float32 is needed
    pts1=np.float32(pts1)
    pts2=np.float32(pts2)

    #calculate the transformation
    #transformation = cv2.estimateRigidTransform(pts1, pts2, False)

    afCoeff,inlier=cv2.estimateAffinePartial2D(pts1, pts2)
    
    
    #extract the parameters


    sx=np.sqrt(afCoeff[0,0]**2+afCoeff[0,1]**2)
    sy=np.sqrt(afCoeff[1,0]**2+afCoeff[1,1]**2)
        
    xd=afCoeff[0,2]
    yd=afCoeff[1,2]
    
    rotation=np.arctan2(afCoeff[1,0]/np.sqrt(afCoeff[0,0]**2+afCoeff[0,1]**2),
                        afCoeff[1,1]/np.sqrt(afCoeff[1,0]**2+afCoeff[1,1]**2))


    return afCoeff,xd,yd,sx,sy,rotation

    
 
def makeAdjacentList(ff,armLength):

    """

    construct the list of adjacent cobras. This is used for fibre identifications

    input
      list of coordinates (nX2)
      armlengths (nx1)

    output:
      ordered list of adjacent cobras
    """

    adjacent=[]

    #use cKDTree for faster distances

    
    cobraTree=cKDTree(ff)

    for i in range(len(ff[:,0])):
        #list of adjacent centers

        #factor of 2.2 pads in case of variation in arm length from measurement to measurement
        ind1=cobraTree.query_ball_point(ff[i],armLength[i]*2.2)
        #remove the central value 
        ind2=cobraTree.query_ball_point(ff[i],1)
        ind1.remove(ind2[0])
        adjacent.append(ind1)
    return(adjacent)


def fibreID(points,targets,centers,arms,dotPos,adjacentCobras,goodIdx):

    """
    runs the fibreID code for fibre identification

    input
  
    points: nx9 array of spot positions. Later columns are shape info
    targets: nx3 array of expected target positions *AFTER THIS MOVE* (not the final target)
    centers: nx3 array of cobra center of motion
    arms: nx1 array of arm length (L1+L2)
    dotPos: nx4 array of dot positions. Final column is dot radius
    adjacentCobras: list of adjacent cobras
    goodIdx: index of good cobras, for bookkeeping

    returns

    cobraMatch: nx5 array of results, with fields cobraID, spotID, x position, y position, flags


    """
    
    
    anyChange=0
    nPoints=points.shape[0]
    nCobras=targets.shape[0]

    #set up variables
    aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch=prepWork(points,nPoints,nCobras,centers,arms,goodIdx)

    #first pass - assign cobra/spot pairs based on the spots poiint of view
    aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange=firstPassNew(aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange)

    #second pass - assign cobra/spot pairs based on the cobra point of view
    aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange=secondPassNew(aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,adjacentCobras,anyChange)

    #last pass - figure out the spots that can belong to more than one cobra, and things hidden by dots
    aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange=lastPassDist(aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,points,targets,centers,'point',anyChange)

    #some final tidying up
    aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange=secondPassNew(aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,adjacentCobras,anyChange)

    #turn the results into an array to be written to teh database
    cobraMatch=np.empty((nCobras,5))
    
    ii=0
    for i in range(int(goodIdx[-1]+2)):
        if(i in goodIdx):

            #cobras assigned to spots
            if(len(potPointMatch[ii])==1):
                cobraMatch[ii,0]=i+1
                cobraMatch[ii,1]=potPointMatch[ii][0]
                cobraMatch[ii,2]=points[potPointMatch[ii][0],1]
                cobraMatch[ii,3]=points[potPointMatch[ii][0],2]
                cobraMatch[ii,4]=0

            #cobras assigned to dots - set spot_id to -1 and set flag
            elif(len(potPointMatch[ii])==0):
                cobraMatch[ii,0]=i+1
                cobraMatch[ii,1]=-1
                cobraMatch[ii,2]=dotPos[ii,1]
                cobraMatch[ii,3]=dotPos[ii,2]
                cobraMatch[ii,4]=2
            ii=ii+1
            
    return cobraMatch
    
def nearestNeighbourMatching(points,targets,nTarg):

    """
    simple matching for fiducial fibres

    input
       points: set to match (nx3)
       targets: set to match to (nx3)
       nTarg: length of targets

    """

    #use cKDTree for speed
    pointTree=cKDTree(points[:,1:3])

    matchPoint=np.zeros((nTarg,3))
    for i in range(nTarg):
        dd,ii=pointTree.query(targets[i,1:3],k=1)
        matchPoint[i]=points[ii]

    return matchPoint
    
    
def prepWork(points,nPoints,nCobras,centers,arms,goodIdx,armFudge=1):

    """
    Create initial list of potential cobra/pooint matches
    and assigned/unasigned cobras. 

    input
    
    points: list of points in n x 3 array (ID, x, y)
    nPoints: number of points
    nCobras: number of cobras
    centers: centers of cobras
    arms: list of arm lenghts (l1+l2)
    armFucge: amount in pixels, by which to increae the arm length to take into account measurement uncertainties

    note that the cobra values aer assumed to be for good cobras (ie, goodIdx)

    returns
    unaCobras,aCobras: lists of assigned and unassigned cobras. union(unaCobras,aCobras)=all cobras
    unaPoints,aPoints: lists of assigned and unassigned points
    potCobraMatch: nPoint long list of list of cobras which could potentially match with each point
    potPointMatch: nCobra long list of list of points which could potentially match with each cobra

    Note that we can't assume the same arm lenght for all points (which would allow fo faster calculation methods)
    as there is enough variation to amke a difference

    """

    potCobraMatch=[]
    potPointMatch=[]
    unaCobras=list(range(nCobras))
    unaPoints=list(range(nPoints))
    aCobras=[]
    aPoints=[]
    dotCobras=[]

    #get the distnace between cobras and points. cdist is pretty fast, check total time
    D = cdist(points[:,1:3], centers[:,1:3])
                   

    #find the cobras which are within arm length of each point and add to the list
    print(D.shape,arms.shape)
    for i in range(nPoints):
        ind1=np.where(D[i,:] < (arms+armFudge))
        
        potCobraMatch.append(list(ind1[0]))

    #now the mirror - find the points which are within arm length of each cobra and add to the list
    for i in range(nCobras):
        ind1=np.where(D[:,i] < (arms[i]+armFudge))
        potPointMatch.append(list(ind1[0]))

    return aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch
    
def firstPassNew(aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange):


    """

    first run through the points, assigning anything that can be done logically from the point side.

    - assign points that are only reachable by one cobra
    - iterate, taking into account cobras that have been assigned (ie, newly freed cobras)
    - the bookkeeping is done via several lists

    input: (all created by prepWork)

    unaCobras,aCobras,dotCobras: lists of assigned, unassigned, and dot-associated cobras. union(unaCobras,aCobras,dotCobras)=all cobras
    unaPoints,aPoints: lists of assigned and unassigned points
    potCobraMatch: nPoint long list of list of cobras which could potentially match with each point
    potPointMatch: nCobra long list of list of points which could potentially match with each cobra


    returns:
      updated versions of the input variables
      anyChange: flag if any change has been made to the input variables (for iteration purposes)

    """

    change=1
    nLoop=0
    
    #loop until no more changes
    while(change==1):
        nLoop=nLoop+1
        change=0

        #pick out the single match cobras, ie case like potCobraMatch[iPoint]=[12] where there is
        #only one point the cobra can be matched with
       
        #go through unassigned points
        for iPoint in unaPoints:
            #if there's only one possible cobra match
            if(len(potCobraMatch[iPoint])==1):
                
                #get the cobra number
                iCob=potCobraMatch[iPoint][0]

                #note that a change has been made
                change=1
                anyChange=1

                #update assigned/unassigned lists
                unaPoints.remove(iPoint)
                aPoints.append(iPoint)

                #remove other points from teh point match list
                potPointMatch[iCob]=[iPoint]

                #this is for debugging and should be removed for final implementation
                try:
                    unaCobras.remove(iCob)
                    aCobras.append(iCob)
                except:
                    print("AAA",iCob)

                #and remove the newly assigned match from all other potential match lists
                for l in unaPoints:
                    if(iCob in potCobraMatch[l]):
                        potCobraMatch[l].remove(iCob)
                for l in unaCobras:
                    if(iPoint in potPointMatch[l]):
                        potPointMatch[l].remove(iPoint)

        #we then need to iterate, to take into account newly freed cobras
                    
    return aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange


def secondPassNew(aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,adjacentCobras,anyChange):


    """
    the second pass deals with things from the cobra's perpective. This is a little tricker due to dots;
    if a point oly matchest with a single cobra, it is unique, but a cobra matched with a single point
    is only uniq if all the surrounding points are assigned.

    input
      same as firstPass

    returns
      same as firstPass
    """

    #set the change variables
    change=1
    nLoop=0
    with open('dump.pkl', 'wb') as output:
        pickle.dump([unaCobras,potPointMatch,adjacentCobras],output)
    #loop until no more changes
    while(change==1):
        change=0

        #cycle through the unassigned cobras
        for iCobra in unaCobras:

            #no potential match for the cobra, point must be hidden by dot
            if(len(potPointMatch[iCobra])==0):

                #update variables
                unaCobras.remove(iCobra)
                dotCobras.append(iCobra)
                change=1
                anyChange=1
                for l in unaCobras:
                    if(iCobra in potPointMatch[l]):
                        potPointMatch[l].remove(iCobra)
            #if there is one potential match, check the surroudning cobras. If they are all assigned,
            #there can't be a dot involved, and we can assign the cobra-point pair
            
            if(len(potPointMatch[iCobra])==1):
                
                #check the assignment
                allAss=checkAdjacentCobras(iCobra,adjacentCobras,unaCobras)
                if(allAss==True):
                    change=1
                    anychange=1
                    #update the lists
                    unaCobras.remove(iCobra)
                    aCobras.append(iCobra)

                    for l in unaPoints:
                        if(iCob in potCobraMatch[l]):
                            potCobraMatch[l].remove(iCob)
                    for l in unaCobras:
                        if(iPoint in potPointMatch[l]):
                            potPointMatch[l].remove(iPoint)

            
    return aCobras,unaCobras,dotCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange


def lastPassDist(aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,points,targets,centers,dFrom,anyChange):

    """
    final pass to deal with tricky assignments (more than one potential cobra assignment, dots)

    this version has two modes for sorthing things out; distance from the expected target, or distance from a line
    connecting the cobra centre to the final target.

    input
      same as firstPass, plus 
         targets: expected targets *AFTER THIS MOVE* (not final target)
         centers: cobra centres
         dFrom: flag for method: 'target' for target, 'line' for line

    returns
         updated variables

    """

    #temporary list of unassigned cobras, to keep track
    tempUnaCobras=copy.deepcopy(unaCobras)
    tempUnaPoints=copy.deepcopy(unaPoints)

    #cycle through unassigned cobras

    #cycle through unassigned cobras
    for iCob in tempUnaCobras:
        elem=potPointMatch[iCob]

        #is there more than one match for the point
        if(len(elem) > 1):

            #two distance measures. "target" is from teh predicted intermediate target,
            #"line" the distance from a line between the cobra center and predicted target (ie, radial move)
            if(dFrom=="target"):
                dTar=np.sqrt((points[elem,1]-targets[iCob,1])**2+(points[elem,2]-targets[iCob,2])**2)
            elif(dFrom=="line"):
                dTar=distancePointLine(targets[iCob,1:2],centers[iCob,1:2],points[elem,1:2])
            else:
                print("dFrom invalid")
                return

            iPoint=elem[dTar.argmin()]

            #update variables
            unaCobras.remove(iCob)
            aCobras.append(iCob)
            unaPoints.remove(iPoint)
            aPoints.append(iPoint)

            for l in tempUnaPoints:
                if(iCob in potCobraMatch[l]):
                    potCobraMatch[l].remove(iCob)
            for l in tempUnaCobras:
                if(iPoint in potPointMatch[l]):
                    potPointMatch[l].remove(iPoint)

        #if only one point, we can update
        elif(len(elem) == 1):
            iPoint=elem[0]
            potPointMatch[iPoint]=[iPoint]
            unaCobras.remove(iCob)
            aCobras.append(iCob)
            unaPoints.remove(iPoint)
            aPoints.append(iPoint)

            for l in tempUnaPoints:
                if(iCob in potCobraMatch[l]):
                    potCobraMatch[l].remove(iCob)
            for l in tempUnaCobras:
                if(iPoint in potPointMatch[l]):
                    potPointMatch[l].remove(iPoint)

   
    return aCobras,unaCobras,aPoints,unaPoints,potCobraMatch,potPointMatch,anyChange
 
def applyAffineTransform(points,afCoeff):

    """
    apply an affine transform to a set of points

    this is a wrapper for the CV2 routines.

    Input: 
      points nx3 array of positions
      afCoeff: transformation matrix
    """
    
    pts=np.zeros((1,points.shape[0],2))

    pts[0,:,0]=points[:,1]
    pts[0,:,1]=points[:,2]

    #cv2 needs 32 bit floats.
    pts=np.float32(pts)

    #do the transform
    pts1=cv2.transform(pts,afCoeff)

    #more bookkeeping
    xx=pts1[0,:,0]
    yy=pts1[0,:,1]

    return np.array([points[:,0],xx,yy]).T

  
def getThresh(image,cobraPos,threshMethod,sigmaThresh,threshFact,findSigma,centSigma):

    """
    wrapper for getting threshold

    there are two modes. 

      calib is specifically for the case when the expected location of the spots on the camera
      is unknown, and therefor the illuminated region. The image statistics in the illuminated/
      unilluminated parts are different enough to skew the results. 

      It was designed for the pinhole mask, but appears to work for the PFI case.

      otherwise, the expected position of the cobras is used to figure out the illuminated region. 


    input: 
      image: image array
      cobraPos: positions of the cobras **in MCS pixels** relative to image
      sigmaFind, sigmaCent: sigma values for finding/centroiding
      sigmaThresh: sigma for sigmaclipping image
      threshMethod: flag for method

    returns
       threshFind: threshold for finding spots, in pixel values
       threshCent: threshold for centroiding spots, in pixel values

    """

    #engineering mode
    if(threshMethod=='calib'):
        xrange,yrange=getRegion(image,threshSigma,threshFact)
        a=getManualThresh(image,xrange,yrange,threshSigma)
    else:

        #get the centre
        mpx=cobraPos[:,1].mean()
        mpy=cobraPos[:,2].mean()
        #grid of pixels value, and distance from the center
        xx,yy=np.meshgrid(np.arange(image.shape[0]),np.arange(image.shape[1]))
        dfromc=np.sqrt((xx-mpx)**2+(yy-mpy)**2)

        #crop to a radius of 1500 for hte whole PFI, to select the illuminated region
    
        ind=np.where(dfromc < 1500)
        #sigma clip values
        a,b,c=sigmaclip(image[ind],sigmaThresh,sigmaThresh)

        #return the mean + sigma value
        threshFind=a.mean()+a.std()*findSigma
        threshCent=a.mean()+a.std()*centSigma

        return threshFind,threshCent


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

    
def checkAdjacentCobras(iCobra,adjacentCobras,unaCobras):

    """
    check if the cobras adjacent to the given cobra aer all assigned to spots

    Input
       iCobra: cobra index
       adjacentCobras: list of adjacent cobras (generated by adjacentCobras)
       unaCobras: list of unassigned cobras

    """

    
    allAss=True
    #print(iCobra,adjacentCobras[iCobra])
    for i in adjacentCobras[iCobra]:    
        if(i in unaCobras):
            allAss=False

    return allAss