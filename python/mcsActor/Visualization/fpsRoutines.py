
import sys
import numpy as np
import cv2
import numpy.ma as ma
import itertools
from scipy import optimize
from pfs.utils.coordinates import CoordTransp
from pfs.utils.coordinates import DistortionCoefficients
from astropy.io import fits

import mcsActor.Visualization.visRoutines as vis


def applyAffineFPS(coord, trans):

    xf, yf = vis.transformPointsNew(coord[:, 1], coord[:, 2], trans['xTrans'],
                                    trans['yTrans'], trans['angle'], trans['xScale'], trans['yScale'])

    return np.array([coord[:, 0], xf, yf]).T


def getAffine(x1, x2, y1, y2):
    """

    Calculate and remove an affine transformation between two sets of
    matched points.

    input:

    x1,y1: coordinates of first set of data (the one that will be transformed)
    x2,y2: coordinates of second set of data (the reference set)

    returns:

    xf,yf transformed points
    trans transformation matrix used, stored in a dictionary

    """

    transform, xd, yd, sx, sy, rotation = vis.getTransform(x1, x2, y1, y2, 1)

    trans = {}
    trans['xTrans'] = xd
    trans['yTrans'] = yd
    trans['xScale'] = sx
    trans['yScale'] = sy
    trans['angle'] = rotation

    return trans


def removeAffine(x1, y1, x2, y2):
    """

    Calculate and remove an affine transformation between two sets of
    matched points.

    input:

    x1,y1: coordinates of first set of data (the one that will be transformed)
    x2,y2: coordinates of second set of data (the reference set)

    returns:

    xf,yf transformed points
    trans transformation matrix used, stored in a dictionary

    """

    transform, xd, yd, sx, sy, rotation = vis.getTransform(x1, x2, y1, y2, 1)

    trans = {}
    trans['xTrans'] = xd
    trans['yTrans'] = yd
    trans['xScale'] = sx
    trans['yScale'] = sy
    trans['angle'] = rotation

    xf, yf = vis.transformPointsNew(x1, y1, xd, yd, rotation, sx, sy)

    return xf, yf, trans


def getFieldDefinition(fieldID):
    """

    program to read in the set of cobras. Currently a dummy program
    reading from file. 

    """

    fiducials = np.loadtxt('fiducials.dat', delimiter=',')
    scienceFibres = np.loadtxt('scienceFibres.dat', delimiter=',')

    return fiducials, scienceFibres


def getInstConfig(fname):
    """

    Routine to retrieve configureation of telescope, including zenith angle and instrument rotation.

    """

    # PUT CODE HERE!!!

    hdu = fits.open(fname)
    hdr = hdu[0].header
    za = 90-hdr['Altitude']
    inr = hdr['inr-str']

    return za, inr


def getFibrePos(fiducials, scienceFibres, za, inr, rotCent, offset):
    """

    Convert a set of fiducials and science fibres (in mask coordinates)
    into expected pixel positons on the image. 


    """

    # rotation adjustment
    inr = inr-180
    if(inr < 0):
        inr = inr+360

    # concatenate - MCS doesn't care which is fiducial and which is science

    xx = np.array([np.concatenate([fiducials[:, 1], scienceFibres[:, 1]])]).ravel()
    yy = np.array([np.concatenate([fiducials[:, 2], scienceFibres[:, 2]])]).ravel()

    # now offset to centre of rotation
    # xx-=offset[0]
    # yy-=offset[1]

    # correect input format
    # xyin=np.array([xx,yy])

    xyinFid = np.array([fiducials[:, 1]-offset[0], fiducials[:, 2]-offset[1]])
    xyinSci = np.array([scienceFibres[:, 1]-offset[0], scienceFibres[:, 2]-offset[1]])

    # call the routine
    xyoutFid = CoordTransp.CoordinateTransform(xyinFid, za, 'pfi_mcs_wofe', inr=inr, cent=rotCent)
    xyoutSci = CoordTransp.CoordinateTransform(xyinSci, za, 'pfi_mcs_wofe', inr=inr, cent=rotCent)

    # reassemble into the right format and return

    return np.array([fiducials[:, 0], xyoutFid[0, :], xyoutFid[1, :]]).T, np.array([scienceFibres[:, 0], xyoutSci[0, :], xyoutSci[1, :]]).T


def getAllFibrePos(fiducials, scienceFibres, za, inr, rotCent, offset):
    """

    Convert a set of fiducials and science fibres (in mask coordinates)
    into expected pixel positons on the image. 


    """

    # rotation adjustment
    inr = inr-180
    if(inr < 0):
        inr = inr+360

    # concatenate - MCS doesn't care which is fiducial and which is science

    xx = np.array([np.concatenate([fiducials[:, 1], scienceFibres[:, 1]])]).ravel()
    yy = np.array([np.concatenate([fiducials[:, 2], scienceFibres[:, 2]])]).ravel()

    # now offset to centre of rotation
    xx -= offset[0]
    yy -= offset[1]

    # correect input format
    xyin = np.array([xx, yy])

    # call the routine
    xyout = CoordTransp.CoordinateTransform(xyin, za, 'pfi_mcs_wofe', inr=inr, cent=rotCent)

    return np.array([np.arange(len(xx)), xyout[0, :], xyout[1, :]]).T


def getDiff(centroids, fibrePos):

    dx = centroids[:, 1]-fibrePos[:, 1]
    dy = centroids[:, 2]-fibrePos[:, 2]

    return dx, dy


def calcRotationCentre(centroids, frameIDs):
    """ 

    Calculate the centre of rotation, given an array retrived by getAllCentroidsFromDB.

    Input
       centroids: Nx9 array of centroid information. Columns 2 and 3 are x and y positions
       frameIDs: list of frameIDs

    """

    xCorner = []
    yCorner = []

    # NEED TO RETRIEVE THE SET OF CENTROIDS FROM

    for i in frameIDs:
        ind = np.where(allCentroids[:, 0] == i)
        x = centroids[ind, 2].ravel()
        y = centroids[ind, 3].ravel()

        x0, x1, y0, y1 = getCorners(x, y)
        xCorner.append(x0)
        yCorner.append(y0)

    xCorner = np.array(xCorner)
    yCorner = np.array(yCorner)

    coords = [xCorner, yCorner]
    xc, yc, r, _ = mcs.least_squares_circle(xCorner, yCorner)

    return xc, yc


def calc_R(x, y, xc, yc):
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
    # coordinates of the barycenter

    #x = np.array([x[0] for x in coords])
    #y = np.array([x[1] for x in coords])
    x_m = np.mean(x)
    y_m = np.mean(y)
    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, args=(x, y))
    xc, yc = center
    Ri = calc_R(x, y, *center)
    R = Ri.mean()
    residu = np.sum((Ri - R)**2)
    return xc, yc, R, residu


def getCentroidsDB(conn, frameID, moveID):
    """

    retrieves a set of centroids from the database, for a sequence of frameIDs

    Input:  
       conn: database connection
       frameID: frameID

    """

    # put code here

    return centroids


def getAllCentroidsDB(conn, frameIDs):
    """

    retrieves a set of centroids from the database, for a sequence of frameIDs

    Input:  
       conn: database connection
       frameIDs: list of frame id numbers

    """

    # for commissioning with MCS

    moveId = 1

    # make a blank array for the centroid array
    centroids = np.array([])

    i = 0

    # cycle through each ID number
    for id in frameIDs:

        # SQL for getting a set of centroids
        #cmd_string = f"""select * from mcsEngTable where frameId={id} and moveId=1"""
        # cmd_string=""
        data = np.array([])
        n = 0
        with conn.cursor() as curs:
            curs.execute(cmd_string)
            rows = curs.fetchall()
            for idx, val in enumerate(rows):
                if idx == 0:
                    data = val
                if idx != 0:
                    data = np.vstack([data, val])
        conn.commit()

        # some data massaging into the right form.
        cen = data[:, 5:11]
        cen1 = np.zeros((cen.shape[0], 7))

        # add an index to the first number
        cen1[:, 0] = i

        # copy over the centroids
        cen1[:, 1:7] = cen

        # create master array
        if(i == 0):
            centroids = cen1
        else:
            centroids = np.concatenate((centroids, cen1), axis=0)

    return centroids
