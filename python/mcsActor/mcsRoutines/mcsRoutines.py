import os
import numpy as np

from scipy.stats import sigmaclip
import mcsActor.mcsRoutines.dbRoutinesMCS as dbTools
import sep

#import ics.cobraCharmer.pfi as pfi

from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
import copy
from scipy import optimize


def getCentroidParams(cmd, configuredCentParms):
    
    """
    Given the default configuration from pfs_instdata, update with any parameters in the command.
    """
    
    try:
        cmdKeys = cmd.cmd.keywords
    except:
        cmdKeys = []

    # returns just the values dictionary
    centParms = configuredCentParms.copy()

    # run through the potential keywords, and update if necessary.
    # if a new keyword is added, update this part
    
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
    
def readCobraGeometry(des,dotData):
    """
    Given the cobra geometry object and dot information,  create the 
    variables used by fibreID. 

    Each variabe is a numpy array, rather than a pandas variable, 
    for computational speed
    """

    #get the list of good cobras
    
    cobs = des.findCobraByCobraIndex(np.arange(0, 2394))
    goodIdx = []
    # cycle through the indices, get the module/cobra number, and check its status
    for i in range(2394):
        cob = des.findCobraByCobraIndex([i])
        status = des.cobraStatus(cob[0][1], cob[0][0])
        if(status == 1):
            goodIdx.append(i)
            
    #return a numpy array for ease of later use
    goodIdx = np.array(goodIdx).astype('int')

    # get the centres and armlengths
    centersAll = des.centers
    centrePos = np.array([np.arange(0,2394), centersAll.real,centersAll.imag]).T
    armLength = (des.L1 + des.L2)

    # find bad arm lengths and put fake values in them
    ind=np.where(np.any([armLength < 3.5,armLength > 6.5],axis=0))
    armLength[ind]=4.5

    # this is a bit kludgy, but is kept for compatibility with the fibre identificaion code.
    # !!todo fix this, but after unit tests written
    goodIdx = np.arange(2394)
    
    nCobras = len(armLength)

    # and the dot positions and sizes
    
    dotPos = np.zeros((len(goodIdx), 4))

    dotPos[:, 0] = dotData['spotId'].values[goodIdx]
    dotPos[:, 1] = dotData['x'].values[goodIdx]
    dotPos[:, 2] = dotData['y'].values[goodIdx]
    dotPos[:, 3] = dotData['r'].values[goodIdx]


    return centrePos, armLength, dotPos, goodIdx, des


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

    cobraTree = cKDTree(np.ascontiguousarray(ff), copy_data = True)
    for i in range(len(ff[:, 0])):
        dd = np.sqrt((ff[i, 1]-ff[:, 1])**2+(ff[i, 2]-ff[:, 2])**2)
        ind1 = np.where(np.all([dd > 0, dd < armLength[i]*2.2], axis = 0))
        adjacent.append(ind1)
        # factor of 2.2 pads in case of variation in arm length from measurement to measurement
        
    return(adjacent)


def fibreId(centroids, centrePos, armLength, tarPos, fids, dotPos, goodIdx, adjacentCobras):

    """
    Main routine for fibre identification. 
    Sets up the variables, calls the sub-routines, sets flags, formats the results.
    """
    
    centers = centrePos
    points = centroids
    nPoints = len(centroids)
    nCobras = len(centrePos)
    arms = armLength
    fidPos =  np.array([fids['fiducialId'],fids['x_mm'],fids['y_mm']]).T

    anyChange = 0
    targets = tarPos

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
    simple matching for fiducial fibres, for use in home position. 
    Not currently used, but is occasionally useful for debugging, so left in. 

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

    Note that we can't assume the same arm length for all points 
    (which would allow fo faster calculation methods)
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
        if len(ind[0]) > 0:
            unaPoints.remove(ind[0][0])
            bPoints.append(ind[0][0])

    #and the same for stuck but illuminated cobras. This is currently a bit of a cludge, based
    #on empirical averages of positions
    
    fileName = os.path.join(os.environ['ICS_MCSACTOR_DIR'],  'etc',  'stuck.txt')
    stuckPos = np.loadtxt(fileName)

    # get the distnace between cobras and points. cdist is pretty fast
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
    the second pass deals with things from the cobra's perspective. This is a little tricker due to dots;
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

    """
    A final pass to finish up. This is the routine that makes decisions for cobras that have more
    than one possible match, that can't be distinguished by pure logic. 

    The deicison is made based on which spot is closest to the targetted position, cycling through
    the unassigned cobras globally from closest to most distant. 
    """

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
    
    # turn inds in to lists so we can remove values. We crop the indices to the
    # 3* the number of unassigned points, which will give us all the points-cobras
    # distances that are w/i a patrol (one cobra can be in at most 3 patrol regions)
    
    ind1 = list(ind[0][0:3*len(unaPoints)])
    ind2 = list(ind[1][0:3*len(unaPoints)])

    # cycle through until no more assignments are possible
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
        #and assign newly singled points to not break things
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

    """
    Threshold routine optimized for bench tests and the unit selector. 
    The main difference is that the noise is very low in the unit selector, 
    which can results in nonsense statistics when most of the background is 
    zero.  If the image median is 0, we calculate thresholds based only on 
    non zero pixels; otherwise the same as getThresh. 
    """

    
    mpx = boreSight[0]
    mpy = boreSight[1]

    fact = 1;

    # grid of pixels value, and distance from the center
    xx, yy = np.meshgrid(np.arange(image.shape[0]), np.arange(image.shape[1]))
    dfromc = np.sqrt((xx-mpx)**2+(yy-mpy)**2)
    
    # crop to a radius of 500 for hte whole PFI, to select the illuminated region

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
    wrapper for getting threshold. It uses the central region of the
    field of view, and calculates sigmaclipped means and standard 
    deviations. 

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
    check if the cobras adjacent to the given cobra are all assigned to spots

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


def mcsPhotometry(image, xPos, yPos, centParms):

    """do aperture photometry on the image"""

    # subtract the background 
    bkg = sep.Background(image)
    image = image - bkg

    # call do the photometry
    flux, fluxerr, flag = sep.sum_circle(image, xPos, yPos, centParms['aperture'],err=bkg.globalrms, gain=1.0, bkgann=(centParms['innerAnn'],centParms['outerAnn']))
    return flux, fluxerr


    
