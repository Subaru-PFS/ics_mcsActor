
"""
Routines for fibre identification

for coding purposes

potCobraMatch is a list BY COBRA with each element containing the indices of potential matching points.
     it is the same length as goodIdx (list of good cobras)
potPointMatch is a list BY POINT with each element containing the indices of potential matching cobras
     it is the same length as points (measures spots)

unaCobras, unaPoints are indices of cobras/points which have not been assigned
aCobras, aPoints are indices of cobras/points 
union(unaCobras,aCobras) = list of cobras
intersection(unaCobras,aCobras) = empty list

to get the final list, get the cobra number from goodIdx[i], and the
coordinates from points[potPointMatch[i],:]

"""

import pfi as pfi
from matplotlib import cm
import itertools
import copy
import seaborn as sns
import matplotlib.pyplot as plt
import centroid as centroid
from scipy.stats import sigmaclip
from scipy.spatial.distance import cdist
import centroid
import pathlib
from astropy.io.fits import getdata
import sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(1, "/Users/karr/Science/PFS/pfsPackages/ics_cobraCharmer/python/procedures/moduleTest/")
sys.path.insert(1, "/Users/karr/Science/PFS/pfsPackages/ics_cobraCharmer/python/ics/cobraCharmer/")
sys.path.insert(1, "/Users/karr/Science/PFS/cobraData/")


def secondPassBench(aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, points, targets):
    """
    second pass trhough the ids to find confused cobras; if found, match based on proximity to targets
    """

    # need copies because we're altering the list

    tempUnaCobras = copy.deepcopy(unaCobras)
    tempUnaPoints = copy.deepcopy(unaPoints)

    # cycle through the cobras
    for iCob in tempUnaCobras:
        elem = potPointMatch[iCob]

        # is there more than one match
        if(len(elem) > 1):
            # get distance from targets
            dTar = np.sqrt((points[elem, 1]-targets[iCob, 1])**2+(points[elem, 2]-targets[iCob, 2])**2)
            # choose the minimum
            iPoint = elem[dTar.argmin()]

            # now we go through the bookkeeping (note, must dump this in a subroutine)
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
        # now account for the points that have had the other element removed; same bookkeeping
        elif(len(elem) == 1):
            iPoint = elem[0]
            potPointMatch[iPoint] = [iPoint]
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

    return aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch


def check(aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, centers, goodIdx, arms, points):
    """
    debugging routine
    """
    ff1 = open("test.reg", "w")
    ff2 = open("miss.reg", "w")
    for i, elem in enumerate(potCobraMatch):
        elem = potCobraMatch[i]

        if(len(elem) == 0):
            pass
        elif(len(elem) > 1):
            # pass
            print("AA", i, elem)
            print("circle point ", points[i, 1], points[i, 2], "#color=blue", file=ff2)
            for e in elem:
                print("cross point ", centers[e, 1], centers[e, 2], "#color=blue", file=ff2)
        else:
            dd = np.sqrt((centers[elem[0], 1]-points[i, 1])**2+(centers[elem[0], 2]-points[i, 2])**2)
            print("cross point ", centers[elem[0], 1], centers[elem[0], 2], "#color=red", file=ff1)
            print("circle point ", points[i, 1], points[i, 2], "#color=green", file=ff1)
            print("line ", centers[elem[0], 1], centers[elem[0], 2],
                  points[i, 1], points[i, 2], "#color=yellow", file=ff1)
            if(dd > arms[elem[0]]+2):
                print("BB", i, elem, dd)
    for i, elem in enumerate(potPointMatch):
        elem = potPointMatch[i]
        if(len(elem) > 1):
            print("CC", i, elem)

    ff1.close()
    ff2.close()


def checkResults(aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch):

    # check to see if an element shows up more than once

    set = []
    for elem in potPointMatch:
        try:
            set.append(elem[0])
        except:
            pass
    aa = Counter(set)
    for elem in aa:
        if(aa[elem] > 1):
            print(elem, aa[elem])
    print("checking potCobraMatch uniq")
    set = []
    for elem in potCobraMatch:
        try:
            set.append(elem[0])
        except:
            pass
    aa = Counter(set)
    for elem in aa:
        if(aa[elem] == 0):
            print(elem, aa[elem])


def checkAdjacent(adjacentCobras):
    """
    quick sanity check from adjacent cobras routine
    """
    fig, ax = plt.subplots()
    for i in np.arange(100, 200, 10):
        ax.scatter(centers[i, 1], centers[i, 2], color='black')
        ax.scatter(centers[adjacentCobras[i], 1], centers[adjacentCobras[i], 2], color='blue')

    ax.set_aspect('equal')


def checkCentroids(points, nPoints):
    """
    quick sanity check for centroids
    """

    print(nPoints)
    fig, ax = plt.subplots()
    ax.scatter(points[:, 1], points[:, 2])


def checkSpiral(spiralX, spiralY):
    """
    quick sanity check for spirals
    """

    fig, ax = plt.subplots()

    palette = itertools.cycle(sns.color_palette())
    cval = np.arange(len(spiralX))
    ax.scatter(spiralX, spiralY, c=cval)


def standardParameters():
    """
    simple wrapper to return standard variables for centroiding
    """

    sigmaFind = 50  # sigma for finding spots
    sigmaCent = 15  # sigma for centroiding spots
    boxFind = 10  # box size for finding spots
    boxCent = 5  # box size for centroiding spots
    nMin = 8  # minimum number of pixels in spot
    nMax = 90  # maximum number of pixels in spot
    maxIt = 20  # maximum number of iterationss
    sigmaThresh = 4  # sigma for determining threshold
    fwhmX = 0  # fwhm for centroiding (=0 to fit)
    fwhmY = 0  # fwhm for centroiding (=0 to fit)

    return sigmaFind, sigmaCent, sigmaThresh, boxFind, boxCent, nMin, nMax, maxIt, fwhmX, fwhmY


def getThresh(image, cobraPos, sigmaFind, sigmaCent, sigmaThresh):
    """
    wrapper for getting threshold, knowing the approximate locations of the spots

    input: 
      image: image array
      cobraPos: positions of the cobras in MCS pixels relative to image
      sigmaFind, sigmaCent: sigma values for finding/centroiding
      sigmaThresh: sigma for sigmaclipping image

    returns
       threshFind: threshold for finding spots, in pixel values
       threshCent: threshold for centroiding spots, in pixel values

    """

    # find the centre of the cobra array
    mpx = cobraPos[:, 1].mean()
    mpy = cobraPos[:, 2].mean()
    # grid of pixels value, and distance from the center
    xx, yy = np.meshgrid(np.arange(image.shape[0]), np.arange(image.shape[1]))
    dfromc = np.sqrt((xx-mpx)**2+(yy-mpy)**2)

    # crop to a radius of 1500 for hte whole PFI, to select the illuminated region

    ind = np.where(dfromc < 1500)
    # sigma clip values
    a, b, c = sigmaclip(image[ind], sigmaThresh, sigmaThresh)

    # return the mean + sigma value
    threshFind = a.mean()+a.std()*sigmaFind
    threshCent = a.mean()+a.std()*sigmaCent

    return threshFind, threshCent


def makeSpirals(rmin, rmax):
    """

    Create a square spiral around a point, for finding arcs.

    Input: rmin, rmax: minimum and maximum radii, in pixels
    Output: spiralx,spiraly: x and y coords of spiral

    """

    x = 0
    y = 0
    d = 1
    m = 1

    spiralx = []
    spiraly = []

    finished = False
    for i in range((2*rmax+1)**2):

        while 2 * x * d < m:
            dd = np.sqrt(x*x+y*y)
            check = np.all([dd > rmin, dd < rmax])
            if(check):
                spiralx.append(x)
                spiraly.append(y)

                #print(x, y,dd)
            x = x + d

        while 2 * y * d < m:
            dd = np.sqrt(x*x+y*y)
            check = np.all([dd > rmin, dd < rmax])

            if(check):
                spiralx.append(x)
                spiraly.append(y)

                #print(x, y,dd)
            y = y + d

        d = -1 * d
        m = m + 1

    return spiralx, spiraly


def confusedFibres(potCobraMatch, aCobras, unaCobras, arc, spiralx, spiraly, points, centers, arcThresh, nCobras, goodIdx):
    """

    routine to handle the final confused fibre sections, which needs the arc images

    """

    # print(potCobraMatch)

    # cycle through the cobras
    for i in range(nCobras):

        # check for newly singled cobras
        if(np.all([len(potCobraMatch[i]) == 1, potCobraMatch[i][0] not in aCobras])):
            aCobras.append(potCobraMatch[i][0])
            tt = potCobraMatch[i][0]

            # if a cobra is newly singled, remove it from the other
            for l in range(i+1, nPoints):
                if(tt in potCobraMatch[l]):
                    potCobraMatch[l].remove(tt)
        if(len(potCobraMatch[i]) > 1):
            assigned = []
            # get point positions
            ddm = 600
            xp = points[i, 1]
            yp = points[i, 2]
            kk = 0

            # find arc
            notFound = True
            ii = 0
            while(np.all([notFound, ii < len(spiralx)])):
                xt = int(spiralx[ii]+yp)
                yt = int(spiraly[ii]+xp)

                if(arc[xt, yt] > arcThresh):
                    notFound = False
                ii = ii+1

            # now get distances from found to
            kk = 0
            toRemove = []
            for jj in potCobraMatch[i]:

                xc = centers[goodIdx[jj]].real
                yc = centers[goodIdx[jj]].imag

                dd = np.sqrt((xc-yt)**2+(yc-xt)**2)
                # print(jj,i)
                #print("  ",int(xp),int(yp))
                #print("  ",int(xc),int(yc))
                #print("  ",xt,yt)
                #print("  ",dd)
                if(dd < ddm):
                    ddm = dd
                    indm = jj
                kk = kk+1
            aCobras.append(indm)
            potCobraMatch[i] = [indm]
            for l in range(i+1, nCobras):
                for elem in toRemove:
                    if(indm in potCobraMatch[l]):
                        potCobraMatch[l].remove(indm)

    return potCobraMatch, aCobras, unaCobras


def prepWork(points, nPoints, nCobras, centers, arms, goodIdx, armFudge=1):
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

    potCobraMatch = []
    potPointMatch = []
    unaCobras = list(range(nCobras))
    unaPoints = list(range(nPoints))
    aCobras = []
    aPoints = []

    # get the distnace between cobras and points. cdist is pretty fast, check total time
    #

    D = cdist(points[:, 1:3], centers[:, 1:3])

    # find the cobras which are within arm length of each point and add to the list
    print(D.shape, arms.shape)
    for i in range(nPoints):
        ind1 = np.where(D[i, :] < (arms+armFudge))

        potCobraMatch.append(list(ind1[0]))

    #D = cdist(centers[:,1:3],points[:,1:3])

    # now the mirror - find the points which are within arm length of each cobra and add to the list
    for i in range(nCobras):
        ind1 = np.where(D[:, i] < (arms[i]+armFudge))
        potPointMatch.append(list(ind1[0]))

    return aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch


def firstPass(aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch):
    """

    first run through the points, assigning anything that can be done logically. 

    - assign points that are only reachable by one cobra
    - iterate, taking into account cobras that have been assigned (ie, newly cobras)
    - the bookkeeping is done via several lists

    input: (all created by prepWork)

    unaCobras,aCobras: lists of assigned and unassigned cobras. union(unaCobras,aCobras)=all cobras
    unaPoints,aPoints: lists of assigned and unassigned points
    potCobraMatch: nPoint long list of list of cobras which could potentially match with each point
    potPointMatch: nCobra long list of list of points which could potentially match with each cobra


    returns:
      updated versions of the input variables
      anyChange: flag if any change has been made to the input variables (for iteration purposes)

    """

    change = 1
    nLoop = 0
    anyChange = 0

    # loop until no more changes
    while(change == 1):
        nLoop = nLoop+1
        change = 0

        # pick out the single match cobras, ie case like potCobraMatch[iPoint]=[12] where there is
        # only one point the cobra can be matched with

        # go through unassigned points
        for iPoint in unaPoints:
            # if there's only one possible cobra match
            if(len(potCobraMatch[iPoint]) == 1):

                # get the cobra number
                iCob = potCobraMatch[iPoint][0]

                # note that a change has been made
                change = 1
                anyChange = 1

                # update assigned/unassigned lists
                unaPoints.remove(iPoint)
                aPoints.append(iPoint)

                # remove other points from teh point match list
                potPointMatch[iCob] = [iPoint]

                try:
                    unaCobras.remove(iCob)
                    aCobras.append(iCob)
                except:
                    print("AAA", iCob)

                # and remove the newly assigned match from all other match lists
                for l in unaPoints:
                    if(iCob in potCobraMatch[l]):
                        potCobraMatch[l].remove(iCob)
                for l in unaCobras:
                    if(iPoint in potPointMatch[l]):
                        potPointMatch[l].remove(iPoint)

    return aCobras, unaCobras, aPoints, unaPoints, potCobraMatch, potPointMatch, anyChange


def secondPass(potCobraMatch, aCobras, unaCobras, aPoints, unaPoints, dots):
    """
    second cycle checking for uniquely identified spots considering the dots

    this is a mirror of firstPass, iterating by cobra rather than point, but
    considering if a single match could be affected by a dot.

    """

    inDot = []
    change = 1
    anyChange = 0
    while(change == 1):
        change = 0
        for iCob in unaCob:
            # no possible matches for cobra
            if(len(potPointMatch[iCob][0]) == 0):
                inDot.append(iCob)

            # if there is once match, check if surrounding cobras are assigned.
            # if they are all matched, the spot for this cobra can't be under
            # a dot

            elif(len(potPointMatch[iCob][0]) == 1):
                allAss = checkAdjacent(iCob, aCobras)

                if(allAss == True):
                    unaCobras.remove(iCob)
                    aCobras.append(iCob)
                    iPoint = potPointMatch[iCob][0]
                    change = 1
                    anyChange = 1
                    try:
                        unaPoint.remove(iPoint)
                        aPoint.append(iPoint)
                    except:
                        pass
                    for l in unaPoints:
                        if(iCob in potCobraMatch[l]):
                            potCobraMatch[l].remove(iCob)
                    for l in unaCobras:
                        if(iPoint in potPointMatch[l]):
                            potPointMatch[l].remove(iPoint)


def checkAdjacent(iCob, aCobras, adjacentCobras):
    """
    For specified cobra, determine if all the adjacent
    cobras have been assigned. If so, return true
    """

    allAss = True
    for cNum in adjacentCobras[iCob]:
        if(cNum not in aCobras):
            allAss = False
            break
    return allAss


def centroidRegions(centers, centroids, nSpots, nCobras, regFile):
    """

    generate region files for visualization

    """

    ofile = open(regFile, "w")
    for i in range(nSpots):
        print("circle point ", centroids[i, 0], centroids[i, 1], file=ofile)


def getCentroids(data, threshFind, threshCent, boxFind, boxCent, nmin, nmax, maxIt):
    """

    wrapper for centroids

    """

    a = centroid.centroid_only(data.astype('<i4'), 0, 0, threshFind, threshCent,
                               boxFind, boxCent, nmin, nmax, maxIt, 0)
    centroids = np.frombuffer(a, dtype='<f8')
    centroids = np.reshape(centroids, (len(centroids)//7, 7))
    nSpots = centroids.shape[0]
    points = np.empty((nSpots, 3))
    points[:, 0] = np.arange(nSpots)
    points[:, 1:3] = centroids[:, 0:2]
    return points, nSpots


def readFiles(spotFile, arcFile):
    """

    read files from disk

    """

    data = getdata(spotFile)
    arc = getdata(arcFile)

    return data, arc


def getCalibBench(xmlFile, broken):
    """

    wrapper to get calibration information; centres, arm lengths, from XML file
    change to database as needed. 

    """

    # initialize adn read the file
    pfic = pfi.PFI(fpgaHost='localhost', doConnect=False,
                   logDir=None)
    aa = pfic.loadModel([pathlib.Path(xmlFile)])

    # get the arrays for all values (including brokens)
    armsAll = pfic.calibModel.L1+pfic.calibModel.L2
    nCobrasAll = len(armsAll)
    centersAll = np.array([np.arange(nCobrasAll), pfic.calibModel.centers.real,
                           pfic.calibModel.centers.imag]).T

    # make goodIdx variable and get number of cobras
    badIdx = np.array(broken, dtype='i4') - 1
    goodIdx = np.array([e for e in range(0, nCobrasAll) if e not in badIdx])

    # just the good ones
    nCobras = len(goodIdx)
    centers = centersAll[goodIdx, :]
    arms = armsAll[goodIdx]

    return pfic, centers, centersAll, goodIdx, arms, armsAll, nCobras, nCobrasAll


def nearTo(point, centers, arms):
    """

    return indices that are within a certain radius of a point

    """

    ind = []
    # for i in range(len(centers)):
    #    dd2=(point[0]-centers[i].real)**2+(point[1]-centers[i].imag)**2
    #    if(dd2 < arms[i]**2):
    #        ind.append(i)
    dd2 = (point[0]-centers.real)**2+(point[1]-centers.imag)**2
    ind = where(dd2 < arms**2)
    return ind


def makeCobraList(ff):
    """

    reorganize the cobras into an array

    """
    cobras = fi.makeCobraList(ff)

    cobraStruct = np.zeros((ff.shape[0], 6))
    cobraStruct[:, 0] = np.arange(1, ff.shape[0]+1).astype('int')
    cobraStruct[:, 1:3] = ff

    return cobraStruct


def makeAdjacentList(ff, armLength):
    """

    construct the list of adjacent cobras

    """

    adjacent = []

    cobraTree = cKDTree(ff)

    for i in range(len(ff[:, 0])):
        # list of adjacent centers
        ind1 = cobraTree.query_ball_point(ff[i], armLength[i]*2.2)
        # remove the central valueclose
        ind2 = cobraTree.query_ball_point(ff[i], 1)
        ind1.remove(ind2[0])
        adjacent.append(ind1)
    return(adjacent)


def testFirstPassSimple():

    ff = np.loadtxt("fiber_position_onchip.txt")
    nCobras = len(ff[:, 0])
    nPoints = nCobras
    centers = ff[:, 0]+ff[:, 1]*1j
    armLength = 55.35
    closest = 23

    goodIdx = np.arange(nCobras).astype('int')

    points, rightCob = simpleTestSet(armLength, closest, numHid=0)

    potCobraMatch, aCobras, unaCobas = firstPassNew(
        points, nPoints, nCobras, centers, armLength*np.ones((nCobras)), goodIdx)

    for i in range(nPoints):
        if(len(potCobraMatch[i] == 1)):
            if(potCobraMatch[i][0] != rightCob[i]):
                print(i, potCobraMatch[i][0], rightCob[i])


def simpleTestSet(armLength, closest, numHid=0, pU=5, nsize=0):
    """

    generate a simple test set, by randomly assigning points inside each cobra patrol region, and 
    filtering for points closer than the tip size.

    if numHid is non zero, numHid randomly chosen spots will be covered by the dots

    """

    # read in cobra centres
    ff = np.loadtxt("fiber_position_onchip.txt")
    if(nsize != 0):
        ff = ff[:nsize, :]

    # get the number of cobras
    nn = len(ff[:, 0])

    # rpoint=armLength*np.sqrt(np.random.uniform(0,1,nn))
    # tpoint=np.random.uniform(0,2*3.14,nn)

    # create the points array
    points = np.zeros((nn-numHid, 3))
    targets = np.zeros((nn, 3))
    points[:, 0] = np.arange(nn-numHid)
    targets[:, 0] = np.arange(nn)
    rightCob = []

    # generate a list of points to be hidden
    hidList = []
    i = 0
    while(i < numHid):
        r = np.random.randint(0, nn-1)
        if r not in hidList:
            hidList.append(r)
            i = i+1

    # cycle through the fibres
    nTot = 0
    for i in range(nn):
        if(i % 100 == 0):

            print(i, end=",")

        # if not hidden
        nota = 1
        # while loop to filter for spots that are too close together
        while(nota):
            # calculate parameters for random generation of positions by area

            rpoint = armLength*np.sqrt(np.random.uniform(0, 1))
            tpoint = np.random.uniform(0, 2*3.14)

            # get the random positions
            xi = ff[i, 0]+rpoint*np.cos(tpoint)
            yi = ff[i, 1]+rpoint*np.sin(tpoint)

            # check distances
            dd = np.sqrt((xi-points[:, 0])**2+(yi-points[:, 1])**2)
            if(dd.min() > closest):
                points[nTot, 1] = xi
                points[nTot, 2] = yi
                nota = 0
        nota = 1
        while(nota):
            offsetX = np.random.normal(loc=0, scale=pU)
            offsetY = np.random.normal(loc=0, scale=pU)
            xiT = xi+offsetX
            yiT = yi+offsetY
            dd = np.sqrt((xiT-targets[:, 1])**2+(yiT-targets[:, 2])**2)
            dC = np.sqrt((xiT-ff[i, 0])**2+(yiT-ff[i, 1])**2)
            if(np.all([dd.min() > closest, dC < armLength])):
                nota = 0
                targets[i, 1] = xiT
                targets[i, 2] = yiT

        if(i not in hidList):
            rightCob.append(i)
            points[nTot, 1] = xi
            points[nTot, 2] = yi

            nTot = nTot+1
    print()
    # ind=random.shuffle(range(nTot))
    return np.array(points), rightCob, hidList, targets


def idFibres(fibreIDs, fibreList, xc, yc, xh, yh):

    # cycle through measured positions and create list of potential ids

    dds = np.dot(xc, yc, xh, yh)


def makeTree():
    ff = np.loadtxt("fiber_position_onchip.txt")
    cobraTree = cKDTree(ff)
    return cobraTree


def potentialMatches(xc, yc, cobraTree, armLength):

    ind = cobraTree.query_ball_point((xc, yc), armLength)
    return ind


def pointLine(p1, p2):
    """ get equation of line between two points """

    slope = (p1[1]-p2[1])/(p1[0]-p2[0])
    intercept = p1[1]-slope*p1[0]

    return slope, intercept


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
