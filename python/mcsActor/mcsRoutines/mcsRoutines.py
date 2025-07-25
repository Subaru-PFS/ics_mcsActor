import os
import numpy as np

from scipy.stats import sigmaclip
import mcsActor.mcsRoutines.dbRoutinesMCS as dbTools
import sep

#import ics.cobraCharmer.pfi as pfi

from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import cv2
import copy
import yaml
from scipy import optimize


def getCentroidParams(cmd, configuredCentParms):
    """Given the default configuration from pfs_instdata, update with any parameters in the command."""
    try:
        cmdKeys = cmd.cmd.keywords
    except:
        cmdKeys = []

    # returns just the values dictionary
    centParms = configuredCentParms.copy()

    if('fwhmx' in cmdKeys):
        centParms['fwhmx'] = cmd.cmd.keywords["fwhmx"].values[0]
    if('fwhmy' in cmdKeys):
        centParms['fwhmy'] = cmd.cmd.keywords["fwhmy"].values[0]

    if('boxFind' in cmdKeys):
        centParms['boxFind'] = cmd.cmd.keywords["boxFind"].values[0]
    if('boxCent' in cmdKeys):
        centParms['boxCent'] = cmd.cmd.keywords["boxCent"].values[0]

    if('findSigma' in cmdKeys):
        centParms['findSigma'] = cmd.cmd.keywords["findSigma"].values[0]
    if('centSigma' in cmdKeys):
        centParms['centSigma'] = cmd.cmd.keywords["centSigma"].values[0]
    if('threshSigma' in cmdKeys):
        centParms['threshSigma'] = cmd.cmd.keywords["threshSigma"].values[0]

    if('nmin' in cmdKeys):
        centParms['nmin'] = cmd.cmd.keywords["nmin"].values[0]
    if('maxIt' in cmdKeys):
        centParms['maxIt'] = cmd.cmd.keywords["maxIt"].values[0]

    if('aperture' in cmdKeys):
        centParms['aperture'] = cmd.cmd.keywords["aperture"].values[0]
    if('innerRad' in cmdKeys):
        centParms['innerRad'] = cmd.cmd.keywords["innerRad"].values[0]
    if('outerRad' in cmdKeys):
        centParms['outerRad'] = cmd.cmd.keywords["outerRad"].values[0]

    return centParms

def readFiducialMasks(fids):

    """
    read good/bad/outer ring fiducial information from yaml file
    """

    instPath = os.path.join(os.environ['PFS_INSTDATA_DIR'])

    fidFile = os.path.join(instPath,"data/pfi/fiducials/fiducialFiberFlags.yaml")
    with open(fidFile, 'r') as inFile:
        fiducialFlags = yaml.safe_load(inFile)

    fidsOuterRing = fids[fids.fiducialId.isin(fiducialFlags['outerRingIds'])]
    badFids = fiducialFlags['badFidIds']
    goodFids = list(set(fids['fiducialId'].values)-set(badFids))
    fidsGood = fids[fids.fiducialId.isin(goodFids)]

    return fidsOuterRing, fidsGood
    
def readCobraGeometry(des,dotData):
    """
    read cobra geometry from configuration file/inst_config

    dot positions from CSVfile at the moment, this will change

    The results will be return in whatever unit the input XML file is in
    """

    # first figure out the good cobras (bad positions are set to 0)
    centersAll = des.centers

    #get the list of good cobras
    cobs = des.findCobraByCobraIndex(np.arange(0, 2394))
    goodIdx = []
    #cycle through the indices, get the module/cobra number, and check its status
    for i in range(2394):
        cob = des.findCobraByCobraIndex([i])
        status = des.cobraStatus(cob[0][1], cob[0][0])
        if(status == 1):
            goodIdx.append(i)
    #return a numpy array for ease of later use
    goodIdx = np.array(goodIdx).astype('int')

    
    # then extract the parameters for good fibres only
    #centrePos = np.array([goodIdx+1, centersAll[goodIdx].real, centersAll[goodIdx].imag]).T
    #armLength = (des.L1[goodIdx]+des.L2[goodIdx])

    # get the centres and armlengths
    centrePos = np.array([np.arange(0,2394), centersAll.real,centersAll.imag]).T
    armLength = (des.L1 + des.L2)

    # find bad arm lengths and put fake values in it
    ind=np.where(np.any([armLength < 3.5,armLength > 6.5],axis=0))
    armLength[ind]=4.5

    goodIdx = np.arange(2394)
    
    # number of cobras
    nCobras = len(armLength)
    dotPos = np.zeros((len(goodIdx), 4))

    dotPos[:, 0] = dotData['spotId'].values[goodIdx]
    dotPos[:, 1] = dotData['x'].values[goodIdx]
    dotPos[:, 2] = dotData['y'].values[goodIdx]
    dotPos[:, 3] = dotData['r'].values[goodIdx]


    return centrePos, armLength, dotPos, goodIdx, des

def calcAffineTransform(pos1, pos2):
    """

    given two sets of registered points, estimate the rigid transformation
    this is a wrapper for the cv2 routine
    Returns transformation matrix and extracted parameters (rotation, translation, scale)

    input:
    pos1 input positions in nx3 shape, first column is an id number, next two coordinates
    xx, yy: transformed positions
    getVales: if == 1, return parameters too

    output: 
    transformation: matrix 
    xd, yd: translations
    sx, sy: scalings
    rotation: rotation (radians)

    """

    sz = pos1.shape[0]
    # turn data into right form
    pts1 = np.zeros((1, sz, 2))
    pts2 = np.zeros((1, sz, 2))

    pts1[0, :, 0] = pos1[:, 1]
    pts1[0, :, 1] = pos1[:, 2]

    pts2[0, :, 0] = pos2[:, 1]
    pts2[0, :, 1] = pos2[:, 2]

    #float32 is needed
    pts1 = np.float32(pts1)
    pts2 = np.float32(pts2)

    # calculate the transformation
    #transformation = cv2.estimateRigidTransform(pts1, pts2, False)

    afCoeff, inlier = cv2.estimateAffinePartial2D(pts1, pts2)

    # extract the parameters

    sx = np.sqrt(afCoeff[0, 0]**2+afCoeff[0, 1]**2)
    sy = np.sqrt(afCoeff[1, 0]**2+afCoeff[1, 1]**2)

    xd = afCoeff[0, 2]
    yd = afCoeff[1, 2]

    rotation = np.arctan2(afCoeff[1, 0]/np.sqrt(afCoeff[0, 0]**2+afCoeff[0, 1]**2), 
                          afCoeff[1, 1]/np.sqrt(afCoeff[1, 0]**2+afCoeff[1, 1]**2))

    return afCoeff, xd, yd, sx, sy, rotation


def makeAdjacentList(ff, armLength):
    """

    construct the list of adjacent cobras. This is used for fibre identifications

    input
      list of coordinates (nX2)
      armlengths (nx1)

    output:
      ordered list of adjacent cobras
    """

    adjacent = []

    # use cKDTree for faster distances
    # temproary fix for system issues, change back!!!

    cobraTree = cKDTree(np.ascontiguousarray(ff), copy_data = True)
    for i in range(len(ff[:, 0])):
        # list of adjacent centers
        dd = np.sqrt((ff[i, 1]-ff[:, 1])**2+(ff[i, 2]-ff[:, 2])**2)
        ind1 = np.where(np.all([dd > 0, dd < armLength[i]*2.2], axis = 0))
        adjacent.append(ind1)
        # factor of 2.2 pads in case of variation in arm length from measurement to measurement
        #ind1 = cobraTree.query_ball_point(ff[i, 1:3], armLength[i]*3)
        ##remove the central value
        #ind2 = cobraTree.query_ball_point(ff[i, 1:3], 1)
        #ind1.remove(ind2[0])
        #adjacent.append(ind1)
        
    return(adjacent)


def fibreId(centroids, centrePos, armLength, tarPos, fids, dotPos, goodIdx, adjacentCobras):

    
    centers = centrePos
    points = centroids
    nPoints = len(centroids)
    nCobras = len(centrePos)
    arms = armLength
    fidPos =  np.array([fids['fiducialId'],fids['x_mm'],fids['y_mm']]).T

    anyChange = 0
    targets = tarPos


    
    # these are the effective number of cobras (ie, goodIdx)
    nPoints = points.shape[0]
    nCobras = targets.shape[0]


    # set up variables
    aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod = prepWork(
        points, nPoints, nCobras, centers, arms, goodIdx, fidPos, armFudge=0.5)


    # first pass - assign cobra/spot pairs based on the spots poiint of view
    aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange = firstPass(
        aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange)


    # second pass - assign cobra/spot pairs based on the cobra point of view
    aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange = secondPass(
        aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, adjacentCobras, assignMethod, anyChange)


    # last pass - figure out the spots that can belong to more than one cobra, and things hidden by dots
    aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange = lastPassDist(
        aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, points, targets, centers, tarPos, 't', assignMethod, anyChange, goodIdx)


    # some final tidying up
    aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange = secondPass(
        aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, adjacentCobras, assignMethod, anyChange)


    # turn the results into an array to be written to teh database
    cobraMatch = np.empty((nCobras, 5), dtype='f4')

    # That is HORRIBLE and DANGEROUS: ids are INTs, and you 
    # basically never need to use doubles. Want something like the following
    # but that requires more effort should go i=on this ticket. I merely changed to 
    # 32-bit reals.
    #cobraMatch = np.empty((nCobras, 5), 
    #                      dtype=[('cobra_id', 'int32'), ('spot_id', 'int32'), 
    #                             ('x_mm', 'float32'), ('y_mm', 'float32'), 
    #                             ('flags', 'int32')])
    ii = 0
    flag = 0
    for i in range(int(goodIdx[-1]+2)):
        if(i in goodIdx):

            # cobras assigned to spots
            if(len(potPointMatch[ii]) == 1):
                cobraMatch[ii, 0] = i+1
                cobraMatch[ii, 1] = points[potPointMatch[ii][0],0]
                cobraMatch[ii, 2] = points[potPointMatch[ii][0], 1]
                cobraMatch[ii, 3] = points[potPointMatch[ii][0], 2]
                cobraMatch[ii, 4] = 0
                if(assignMethod[ii]==0):
                    cobraMatch[ii, 4] += 1

                
            # cobras assigned to dots - set spot_id to -1 and set flag
            elif(len(potPointMatch[ii]) == 0):
                cobraMatch[ii, 0] = i+1
                cobraMatch[ii, 1] = -1
                cobraMatch[ii, 2] = dotPos[ii, 1]
                cobraMatch[ii, 3] = dotPos[ii, 2]
                cobraMatch[ii, 4] = 2
                if(assignMethod[ii]==0):
                    cobraMatch[ii, 4] += 1

            # if there is more than one potential match for a cobra, something has gone badly wrong
            else:
                flag = flag+1
                #print(ii,i,potPointMatch[ii])
            ii = ii+1
            
    return cobraMatch, unaPoints, flag


def nearestNeighbourMatching(points, targets):
    
    """
    simple matching for fiducial fibres, for use in home position

    input
       points: set to match (nx3)
       targets: set to match to (nx3)
       nTarg: length of targets

    """

    nTarg = points.shape[0]

    # use cKDTree for speed
    pointTree = cKDTree(points[:, 1:3])
    matchPoint = np.zeros((nTarg, 3))
    
    for i in range(nTarg):
        dd, ii = pointTree.query(targets[i, 1:3], k = 1)
        ii = np.argmin(dd)
        matchPoint[i] = points[ii]

    return matchPoint

def nearestNeighbourMatchingBore(points, targets, unrot):
    
    """
    simple matching for fiducial fibres, for use in home position

    input
       points: set to match (nx3)
       targets: set to match to (nx3)
       nTarg: length of targets

    """

    nTarg = points.shape[0]

    # use cKDTree for speed
    pointTree = cKDTree(points[:, 1:3])
    matchPoint = np.zeros((len(targets), 3))
    
    for i in range(len(targets)):
        dd, ii = pointTree.query(targets[i, 1:3], k = 1)
        matchPoint[i] = unrot[ii]

    return matchPoint

def prepWork(points, nPoints, nCobras, centers, arms, goodIdx, fidPos, armFudge = 0.08):
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
    unaCobras, aCobras: lists of assigned and unassigned cobras. union(unaCobras, aCobras) = all cobras
    unaPoints, aPoints: lists of assigned and unassigned points
    potCobraMatch: nPoint long list of list of cobras which could potentially match with each point
    potPointMatch: nCobra long list of list of points which could potentially match with each cobra

    Note that we can't assume the same arm lenght for all points (which would allow fo faster calculation methods)
    as there is enough variation to amke a difference

    """

    potCobraMatch = []
    potPointMatch = []
    unaCobras = list(range(nCobras))
    unaPoints = list(range(nPoints))
    aCobras = []
    aPoints = []
    dotCobras = []
    fidPoints = []

    assignMethod=np.zeros(nCobras, dtype=np.int32)

    bPoints = []  # non real points (fids, stuck fibres)

    fileName = os.path.join(os.environ['ICS_MCSACTOR_DIR'],  'etc',  'stuck.txt')

    stuckPos = np.loadtxt(fileName)

    
    #first, quick positional matching to remove fiducial fibres from list of matchable points
    #note that the matching should return either 0 points (unilluminated fiducials) or
    #1 points (match) as the fiducial fibres don't have a patrol radius

    D = cdist(fidPos[:,1:3],points[:,1:3])
    for i in range(len(fidPos)):
        ind = np.where(D[i, :] < 1)
        #print("Fid Match", i, ind, len(ind[0]))
        if len(ind[0]) > 0:
            unaPoints.remove(ind[0][0])
            bPoints.append(ind[0][0])

    #and the same for stuck but illuminated cobras. This is currently a bit of a cludge, based
    #on empirical averages of positions
    fileName = os.path.join(os.environ['ICS_MCSACTOR_DIR'],  'etc',  'stuck.txt')
    stuckPos = np.loadtxt(fileName)

    #D = cdist(stuckPos[:,1:3], points[:,1:3])
    #for i in range(len(stuckPos)):
    #    ind = np.where(D[i, :] < 1)
    #    if len(ind[0]) > 0:
    #        unaPoints.remove(ind[0][0])
    #        bPoints.append(ind[0][0])
    # get the distnace between cobras and points. cdist is pretty fast, check total time
    D = cdist(points[:, 1:3], centers[:, 1:3])

    # find the cobras which are within arm length of each point and add to the list
    
    for i in range(nPoints):
        ind1 = np.where(D[i, :] < (arms+armFudge))

        potCobraMatch.append(list(ind1[0]))

    # now the mirror - find the points which are within arm length of each cobra and add to the list
    for i in range(nCobras):
        ind1 = np.where(D[:, i] < (arms[i]+armFudge))
        potPointMatch.append(list(ind1[0]))

    # now remove the non cobra points
    for iPoint in bPoints:
        for l in unaCobras:
            if(iPoint in potPointMatch[l]):
                potPointMatch[l].remove(iPoint)

    return aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod


def firstPass(aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange):
    """

    first run through the points, assigning anything that can be done logically from the point side.

    - assign points that are only reachable by one cobra
    - iterate, taking into account cobras that have been assigned (ie, newly freed cobras)
    - the bookkeeping is done via several lists

    input: (all created by prepWork)

    unaCobras, aCobras, dotCobras: lists of assigned, unassigned, and dot-associated cobras. union(unaCobras, aCobras, dotCobras) = all cobras
    unaPoints, aPoints: lists of assigned and unassigned points
    potCobraMatch: nPoint long list of list of cobras which could potentially match with each point
    potPointMatch: nCobra long list of list of points which could potentially match with each cobra


    returns:
      updated versions of the input variables
      anyChange: flag if any change has been made to the input variables (for iteration purposes)

    """

    change = 1
    nLoop = 0

    # loop until no more changes
    while(change == 1):
        nLoop = nLoop+1
        change = 0

        # pick out the single match cobras, ie case like potCobraMatch[iPoint] = [12] where there is
        # only one point the cobra can be matched with

        # go through unassigned points
        for iPoint in unaPoints:
            # if there's only one possible cobra match
            if(len(potCobraMatch[iPoint]) == 1):

                # get the cobra number
                iCob = potCobraMatch[iPoint][0]

                assignMethod[iCob]=1
                # note that a change has been made
                change = 1
                anyChange = 1

                # update assigned/unassigned lists
                unaPoints.remove(iPoint)
                aPoints.append(iPoint)

                # remove other points from teh point match list
                potPointMatch[iCob] = [iPoint]

                # this is for debugging and should be removed for final implementation
                try:
                    unaCobras.remove(iCob)
                    aCobras.append(iCob)
                except:
                    print("AAA", iCob)

                # and remove the newly assigned match from all other potential match lists
                for l in unaPoints:
                    if(iCob in potCobraMatch[l]):
                        potCobraMatch[l].remove(iCob)
                for l in unaCobras:
                    if(iPoint in potPointMatch[l]):
                        potPointMatch[l].remove(iPoint)

        # we then need to iterate, to take into account newly freed cobras

    return aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange


def secondPass(aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, adjacentCobras, assignMethod, anyChange):
    """
    the second pass deals with things from the cobra's perpective. This is a little tricker due to dots;
    if a point oly matchest with a single cobra, it is unique, but a cobra matched with a single point
    is only uniq if all the surrounding points are assigned.

    input
      same as firstPass

    returns
      same as firstPass
    """

    # set the change variables
    change = 1
    nLoop = 0

    # loop until no more changes
    while(change == 1):
        change = 0

        # cycle through the unassigned cobras
        for iCobra in unaCobras:

            # no potential match for the cobra, point must be hidden by dot
            if(len(potPointMatch[iCobra]) == 0):
                # update variables
                unaCobras.remove(iCobra)
                dotCobras.append(iCobra)
                assignMethod[iCobra]=1
                change = 1
                anyChange = 1
                for l in unaCobras:
                    if(iCobra in potPointMatch[l]):
                        potPointMatch[l].remove(iCobra)
            # if there is one potential match, check the surroudning cobras. If they are all assigned, 
            # there can't be a dot involved, and we can assign the cobra-point pair

            elif(len(potPointMatch[iCobra]) == 1):

                # check the assignment
                allAss = checkAdjacentCobras(iCobra, adjacentCobras, unaCobras)
                if(allAss == True):
                    iPoint = potPointMatch[iCobra][0]

                    change = 1
                    anychange = 1
                    # update the lists
                    unaCobras.remove(iCobra)
                    aCobras.append(iCobra)
                    assignMethod[iCobra]=1

                    try:
                        unaPoints.remove(iPoint)
                    except:
                        pass
                    aPoints.append(iPoint)

                    for l in unaPoints:
                        if(iCobra in potCobraMatch[l]):
                            potCobraMatch[l].remove(iCobra)
                    for l in unaCobras:
                        if(iPoint in potPointMatch[l]):
                            potPointMatch[l].remove(iPoint)
                    potCobraMatch[iPoint] = [iCobra] #!!
                    potPointMatch[iCobra] = [iPoint] #!!

    return aCobras, unaCobras, dotCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange


def lastPassDist(aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, points, targets, centers, prevPos, dFrom, assignMethod, anyChange, goodIdx):

    # temporary list of unassigned cobras, to keep track
    tempUnaCobras = copy.deepcopy(unaCobras)
    tempUnaPoints = copy.deepcopy(unaPoints)

    # repeat of first pass to account for newly singled points/cobras
    change = 1

    
    masterUnaPoints = np.copy(unaPoints)
    masterUnaCobras = np.copy(unaCobras)

    while(change == 1):
        change = 0
        for iPoint in unaPoints:
            elem = potCobraMatch[iPoint]

            if(len(elem) == 1):
                iCob = elem[0]
                if(iCob in unaCobras):
                    unaCobras.remove(iCob)
                aCobras.append(iCob)
                assignMethod[iCob]=1
            
                if(iPoint in unaPoints):
                    unaPoints.remove(iPoint)
                aPoints.append(iPoint)

                for l in tempUnaPoints:
                    if(iCob in potCobraMatch[l]):
                        potCobraMatch[l].remove(iCob)
                for l in tempUnaCobras:
                    if(iPoint in potPointMatch[l]):
                        potPointMatch[l].remove(iPoint)

                potCobraMatch[iPoint] = [iCob]
                potPointMatch[iCob] = [iPoint]
                change = 1


    nLoop=0
    change = 1

    tsum=0

    #get distances between unassigned points and unassigned cobras
    #sort, to find the lowest value
    D = cdist(points[unaPoints, 1:3], targets[unaCobras, 1:3])
    ind = np.unravel_index(np.argsort(D, axis = None), D.shape)
    


    nnn=0
    # turn inds in to lists so we can remove values. We crop the indices to the
    # 3* the number of unassigned points, which will give us all the points-cobras
    # distances that are w/i a patrol (one cobra can be in at most 3 patrol regions)
    
    ind1 = list(ind[0][0:3*len(unaPoints)])
    ind2 = list(ind[1][0:3*len(unaPoints)])

    while(change == 1):

        change = 0
        
        nLoop=nLoop+1
        
        #now find the smallest separation among the allowed values
        found = False
        i = 0
            
        while (not found and i < len(ind1)):

            # use the master values to keep track of the original indices. 
            iPoint = masterUnaPoints[ind1[i]]
            iCob = masterUnaCobras[ind2[i]]
        
            if(iPoint in potPointMatch[iCob]):
                if(iCob in potCobraMatch[iPoint]):
                   found = True
            else:
                pass
                    
            i = i+1
        if(found == True):
        
                change = 1
                #update variables
                if(iPoint in unaPoints):
                    unaPoints.remove(iPoint)
                aPoints.append(iPoint)
                if(iCob in unaCobras):
                    unaCobras.remove(iCob)
                aCobras.append(iCob)
                
                for l in tempUnaPoints:
                    if(iCob in potCobraMatch[l]):
                        potCobraMatch[l].remove(iCob)
                for l in tempUnaCobras:
                    if(iPoint in potPointMatch[l]):
                        potPointMatch[l].remove(iPoint)

                potCobraMatch[iPoint] = [iCob]
                potPointMatch[iCob] = [iPoint]

                ind1.remove(ind1[i-1])
                ind2.remove(ind2[i-1])


                
        #print()
        #get singled points to not break things
        nchange = 1
 
        while(nchange == 1):
            nchange = 0
            for iPoint in unaPoints:
                elem = potCobraMatch[iPoint]
        
                if(len(elem) == 1):
                    iCob = elem[0]
                    
                    if(iCob in unaCobras):
                        unaCobras.remove(iCob)
                        
                    aCobras.append(iCob)
                
                    if(iPoint in unaPoints):
                        unaPoints.remove(iPoint)
                    aPoints.append(iPoint)
        
                    for l in tempUnaPoints:
                        if(iCob in potCobraMatch[l]):
                            potCobraMatch[l].remove(iCob)
                    for l in tempUnaCobras:
                        if(iPoint in potPointMatch[l]):
                            potPointMatch[l].remove(iPoint)
        
                    potCobraMatch[iPoint] = [iCob]
                    potPointMatch[iCob] = [iPoint]
                    nchange = 1
    
    return aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, assignMethod, anyChange

def getThreshBench(image, boreSight, sigmaThresh, findSigma, centSigma):
  
    mpx = boreSight[0]
    mpy = boreSight[1]

    fact = 1;

    # grid of pixels value, and distance from the center
    xx, yy = np.meshgrid(np.arange(image.shape[0]), np.arange(image.shape[1]))
    dfromc = np.sqrt((xx-mpx)**2+(yy-mpy)**2)
    
    # crop to a radius of 1500 for hte whole PFI, to select the illuminated region

    ind = np.where(dfromc < 500)
    if(np.median(image[ind])==0):
        ind = np.where(image > 0)
        fact = 2;

    # sigma clip values
    a, b, c = sigmaclip(image[ind], sigmaThresh, sigmaThresh)

    # return the mean + sigma value
    threshFind = (a.mean()+a.std()*findSigma) / fact
    threshCent = (a.mean()+a.std()*centSigma) / fact

    return threshFind, threshCent, a.mean()

def getThresh(image, boreSight, sigmaThresh, findSigma, centSigma):
    """
    wrapper for getting threshold

    input: 
      image: image array
      cobraPos: positions of the cobras **in MCS pixels** relative to image
      sigmaFind, sigmaCent: sigma values for finding/centroiding
      sigmaThresh: sigma for sigmaclipping image
      threshMethod: flag for method

    returns
       threshFind: threshold for finding spots, in pixel values
       threshCent: threshold for centroiding spots, in pixel values
       avBack: average background value

    """

    # get the centre
    mpx = boreSight[0]
    mpy = boreSight[1]
    # grid of pixels value, and distance from the center
    xx, yy = np.meshgrid(np.arange(image.shape[0]), np.arange(image.shape[1]))
    dfromc = np.sqrt((xx-mpx)**2+(yy-mpy)**2)
    
    # crop to a radius of 1500 for hte whole PFI, to select the illuminated region
    
    ind = np.where(dfromc < 500)
    # sigma clip values
    a, b, c = sigmaclip(image[ind], sigmaThresh, sigmaThresh)

    # return the mean + sigma value
    threshFind = a.mean()+a.std()*findSigma
    threshCent = a.mean()+a.std()*centSigma

    return threshFind, threshCent, a.mean()

def checkAdjacentCobras(iCobra, adjacentCobras, unaCobras):
    """
    check if the cobras adjacent to the given cobra aer all assigned to spots

    Input
       iCobra: cobra index
       adjacentCobras: list of adjacent cobras (generated by adjacentCobras)
       unaCobras: list of unassigned cobras

    """

    allAss = True
    for i in adjacentCobras[int(iCobra)][0]:
        if(i in unaCobras):
            allAss = False

    return allAss

def distancePointLine(p1, p2, p):
    """ get closest distance from a line between two points and another point """

    y1 = p1[1]
    x1 = p1[0]
    y2 = p2[1]
    x2 = p2[0]
    x = p[0]
    y = p[1]

    d = abs((y2-y1)*x+(x2-x1)*y+x2*y1-y2*x1)/np.sqrt((y2-y1)**2+(x2-x1)**2)

    return p


def extractRotCent(afCoeff):    
    """extract centre of rotation from affine matrix"""
        
    #affine matrix is of the form
    #alpha   beta   (1-alpha)*xc-beta*yc
    #-beta   alpha  beta*xc-(1-alpha)*yc 
    #so we solve for xc and yc
    
    A = 1-afCoeff[0, 0]
    B = afCoeff[0, 1]

    xd = afCoeff[0, 2]
    yd = afCoeff[1, 2]

    yc = (yd-B/A*xd)/(B**2/A+A)
    xc = (xd+B*yc)/A

    return xc, yc

def f(c, x, y):
    """
    calculate the algebraic distance between the data points
    and the mean circle centered at c=(xc, yc)
    """
    Ri = calc_R(x, y, *c)
    return Ri - Ri.mean()


def least_squares_circle(x, y):
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

    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x, y))
    xc, yc = center
    Ri = calc_R(x, y, *center)
    R = Ri.mean()
    residu = np.sum((Ri - R)**2)
    return xc, yc, R, residu

def calc_R(x, y, xc, yc):
    
    """
    calculate the distance of each 2D points from the center (xc, yc)
    """
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def initialBoresight(db, frameIDs):

    """ 
    An initial approximation of the boresight from a sequence of spot measurements 

    We calculate the mean position of the spots for each frame, and then fit a circle
    to the sequence of means. As the rotation centre is slightly offset from the image centre,
    the mean positions will trace a circle around the centre of rotation.

    """
    
    xC = []
    yC = []
    #for each frame, pull the spots from database, and calculate mean position
    for frameId in frameIDs:
        # retrieve all spots
        points = dbTools.loadCentroidsFromDB(db, frameId)

        # get means
        xC.append(np.nanmean(points[:, 1]))
        yC.append(np.nanmean(points[:, 2]))

    # do the fit
    xCentre, yCentre, radius, residuals = least_squares_circle(xC, yC)
    
    return [xCentre, yCentre]

def refineBoresight(db, frameId1, frameId2, boresightEstimate):
    """
    Refine the boresight measurement to subpixel accuracy.

    - take two frames at different angles, 
    - rotate using the estimated centre and angle,
    - do nearest neighbour matching
    - calculate the affine transform
    - extract the centre of rotation

    Input: 
    frameId1, frameId2: mcs_frame_id for the two sets of points
    boresightEstimate: [xC,yC] returned by initialBoresight
    
    """

    # retrieve two sets of points
    points = dbTools.loadCentroidsFromDB(db, frameId1)
    points1 = dbTools.loadCentroidsFromDB(db, frameId2)

    points=points[~np.isnan(points).any(axis=1)]
    points1=points1[~np.isnan(points1).any(axis=1)]
    
    # and their instrumetn rotation
    zenithAngle1, insRot1 = dbTools.loadTelescopeParametersFromDB(db, frameId1)
    zenithAngle2, insRot2 = dbTools.loadTelescopeParametersFromDB(db, frameId2)

    thetaDiff = (insRot1 - insRot2) * np.pi/180

    x1 = points[:, 1]
    y1 = points[:, 2]

    # the estimated centre
    x0 = boresightEstimate[0]
    y0 = boresightEstimate[1]

    # rotate to match the second set of point, using theta differnce
    
    x2 = (x1-x0) * np.cos(thetaDiff) - (y1-y0) * np.sin(thetaDiff) + x0
    y2 = (x1-x0) * np.sin(thetaDiff) + (y1-y0) * np.cos(thetaDiff) + y0

    # some bookkeeping for the format expected by nearestNeighbourMatching

    #three sets of poitns, rotated, unrotated and target

    unRot = points[:, 0:3]
    source = np.array([x2, x2, y2]).T
    target = points1[:, 0:3]

    # do nearest neighbour matching on *transformed* values,
    # and return the *untransformed* values matched to the first set

    matchPoint = nearestNeighbourMatchingBore(source, target, unRot)
 
    #gethe affine transform
    afCoeff, xd, yd, sx, sy, rotation = calcAffineTransform(target[:, 0:3], matchPoint[:, 0:3])

    #and extract the centre of rotation from the matrix
    xc, yc = extractRotCent(afCoeff)
        
    return [xc, yc]


def calcBoresight(db, frameIds, pfsVisitId):

    """
    wrapper for boresight calculationg.

    Input: 
       db: database connection
       frameIds: list of mcs_frame_ids for the data set

    Output: 
       returns boresight
       write updated value to database

    """

    boresightEstimate = initialBoresight(db, frameIds)
    boresight = refineBoresight(db, frameIds[0], frameIds[1], boresightEstimate)

    dbTools.writeBoresightToDB(db, pfsVisitId, boresight)

    return boresight

def mcsPhotometry(image, xPos, yPos, centParms):

    """do aperture photometry on the image"""

    # subtract the background 
    bkg = sep.Background(image)
    image = image - bkg

    # call do the photometry
    flux, fluxerr, flag = sep.sum_circle(image, xPos, yPos, centParms['aperture'],err=bkg.globalrms, gain=1.0, bkgann=(centParms['innerAnn'],centParms['outerAnn']))
    return flux, fluxerr


    
