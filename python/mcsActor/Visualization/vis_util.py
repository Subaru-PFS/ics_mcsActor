"""

Miscellaneous utility routines

"""


import numpy as np
from astropy.io.fits import getdata
import matplotlib.pylab as plt


def dbConnect(user='pfs', host='133.40.164.208', passwd='psfpass'):
    """

    Connect to the specified database. 

    Input: user, host, password

    Returns: connection

    """

    try:
        #conn = psycopg2.connect("dbname='fps' user='pfs' host='133.40.164.208' password=pfspass")
        conn = psycopg2.connect("dbname='fps' user='pfs' host='133.40.164.208' password="+passwd)
        print('text="Connected to FPS database."')
        return conn
    except:
        print('text="I am unable to connect to the database."')
        return None


def getFileName(frameID, sourceDir, fPref, dataType):
    """

    get a single file name

    """

    if(dataType == 'taichung'):
        fname = sourceDir+"/"+fPref+str(frameID).zfill(4)+".fits"
        prefix = fPref

    elif(dataType == 'pinhole'):
        fname = sourceDir+"/"+fPref+str(frameID).zfill(6)+"00.fits"
        prefix = "check_"+str(frameID).zfill(5)

    return fname, prefix


def getFileNames(frameId1, frameId2, frameSkip, sourceDir, fPref, dataType):
    """

    Generate a list of file names with complete path.

    the dataType flag toggles between different input sources
    (Taichung lab vs pinhole mask), and can be added to as needed.
    This will need to be updated when there are exposure + move ids.


    Input
       frameId1, frameId2 - first and last frames (inclusive)
       frameSkip - list of frames to be excluded
       sourceDir - full path of source directory
       fPref - file prefix (PFSC for MCS)
       dataType - 'pinhole' - commissioning run format (MCS, move id = 0)
                - 'taichung' - various lab data from testing

    Returns
       files = list of fileanmes (full path)
       prefix - prefix for output plot file names
       centroidFile - filename for saving centroids in local operation
       frameIDs - list of frameIds


    """

    # generate a list of frame IDs, delete any skipped ones

    frameIDs = list(np.arange(frameId1, frameId2+1))
    for ids in frameSkip:
        frameIDs.remove(ids)

    files = []

    # different input data, assemle the file names

    if(dataType == 'taichung'):
        for i in frameIDs:
            files.append(sourceDir+"/"+fPref+str(i).zfill(4)+".fits")

        prefix = fPref
        centroidFile = prefix+"_ncentroids.dat"

    elif(dataType == 'pinhole'):
        for i in frameIDs:
            files.append(sourceDir+"/"+fPref+str(i).zfill(6)+"00.fits")

        prefix = "see_"+str(frameId1).zfill(5)+"_"+str(frameId2).zfill(5)

    centroidFile = prefix+"_centroids.dat"

    return files, prefix, centroidFile, frameIDs


def loadInstParams(config):
    """

    load instrument parameters. Dummy function, update!

    """

    if(config == 'oct18'):
        rotCent = [4691.5, 3095.7]
        offset = [927, 964]

    return rotCent, offset


def getImage(filename):
    """

    Simple wrapper to read an image from file

    """

    image = getdata(filename)

    return image


def toFits(filename):
    """

    Quick routine to convert raw image to fits.

    """

    a = np.fromfile(filename, dtype='uint16')
    image = a.reshape([5778, 8960])

    pf.writeto(filename+".fits", image)
    print(image.min(), image.mean(), image.max())
