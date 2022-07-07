"""
routines for testing fibreID offline
"""

import mcsActor.windowedCentroid.centroid as centroid
import mcsActor.mcsRoutines.mcsRoutines as mcsTools
import mcsActor.mcsRoutines.dbRoutinesMCS as dbTools
from astropy.io import fits
import matplotlib.pylab as plt
from importlib import reload
from pfs.utils import butler
import pandas as pd
import pfs.utils.coordinates.transform as transformUtils
import numpy as np
import sys
from scipy.spatial.distance import cdist

def initFibreID(xmlFile,dotFile):
    centParms = mcsTools.getCentroidParams(None)
    centParms['nmax']=1e6

    centrePos, armLength, dotPos, goodIdx, des = mcsTools.readCobraGeometry(xmlFile,dotFile)
    adjacentCobras = mcsTools.makeAdjacentList(centrePos, armLength)
    fids = pd.read_csv("/Users/karr/Science/PFS/newestMHS/pfs_instdata/data/pfi/fiducial_positions.csv", comment='#')

    fidPos =  np.array([fids['fiducialId'],fids['x_mm'],fids['y_mm']]).T

    return centParms, centrePos, armLength, dotPos, goodIdx, des, adjacentCobras, fids, fidPos



def runMatchFile(frameID,dataPath,centParms,cameraName,fids,centrePos,armLength,dotPos,goodIdx,adjacentCobras,targetPos):


    # get fits file information and data

    centrePosPix = np.array([5033,3447]).T
    fName=dataPath+"PFSC0"+str(frameID)+".fits"
    hdu = fits.open(fName)
    hdr = hdu[0].header
    altitude = hdr['altitude']
    insrot = hdr['insrot']
    image=fits.getdata(fName)

    # get thresholds
    findThresh, centThresh, avBack = mcsTools.getThresh(
        image, centrePosPix,  centParms['threshSigma'], centParms['findSigma'], centParms['centSigma'])
    

    # centroid
    a = centroid.centroid_only(image.astype('<i4'),
                               centParms['fwhmx'], centParms['fwhmy'], findThresh, centThresh,
                               centParms['boxFind'], centParms['boxCent'],
                               centParms['nmin'], centParms['maxIt'], 0)
    centroids = np.frombuffer(a, dtype='<f8')
    centroids = np.reshape(centroids, (len(centroids)//7, 7))
    nSpots = centroids.shape[0]
    points = np.empty((nSpots, 8))
    points[:, 0] = np.arange(nSpots)
    points[:, 1:] = centroids[:, 0:]


    sz = points.shape
    mcsFrameId = frameID
    mcsFrameIDs = np.repeat(mcsFrameId, sz[0]).astype('int')

    # massage the format for the transforms
    
    frame = np.zeros((sz[0], 9))
    frame[:, 0] = mcsFrameIDs
    frame[:, 1:] = points
    # column names
    columns = ['mcs_frame_id', 'spot_id', 'mcs_center_x_pix', 'mcs_center_y_pix', 'mcs_second_moment_x_pix',
               'mcs_second_moment_y_pix', 'peakvalue', 'bgvalue', 'mcs_second_moment_xy_pix']
    mcsData = pd.DataFrame(frame, columns=columns)
    centroids = points

    # load the transformations
    pfiTransform = transformUtils.fromCameraName(cameraName, altitude=altitude, insrot=insrot)

    fidList=list(fids['fiducialId'].values)
    badFid = [1,32,34,61,68,75,88,89,2,4,33,36,37,65,66,67,68,69]

    for i in badFid:
        try:
            fidList.remove(i)
        except:
                pass
    outerRing = np.zeros(len(fids), dtype=bool)   
    goodFid = np.zeros(len(fids), dtype=bool)   
    for i in [29, 30, 31, 61, 62, 64, 93, 94, 95, 96]:
        outerRing[fids.fiducialId == i] = True
    for i in fidList:
        goodFid[fids.fiducialId == i] = True

    # match fiducials
    ffids0,dist0=pfiTransform.updateTransform(mcsData, fids[outerRing], matchRadius=8.0, nMatchMin=0.1)
    ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=4, nMatchMin=0.1)

    # transform
    
    mmCentroids=np.copy(centroids)
    mmCentroids[:,1], mmCentroids[:,2] = pfiTransform.mcsToPfi(centroids[:,1],centroids[:,2])


    
    # identify
    cobraMatch, unaPoints = mcsTools.fibreId(mmCentroids, centrePos, armLength, targetPos, fids, dotPos, goodIdx, adjacentCobras)
    fidPos =  np.array([fids['fiducialId'],fids['x_mm'],fids['y_mm']]).T
    fidMatch=np.zeros((len(fidPos),4))
    # get the matching for fiducials
    D = cdist(fidPos[:,1:3],mmCentroids[:,1:3])
    
    for i in range(len(fidPos)):
        ind = np.where(D[i, :] < 1)
        if len(ind[0]) > 0:
            fidMatch[i,0]=fidPos[i,0]
            fidMatch[i,1]=ind[0][0]
            fidMatch[i,2]=mmCentroids[ind[0][0],1]
            fidMatch[i,3]=mmCentroids[ind[0][0],2]
        else:
            fidMatch[i,0]=fidPos[i,0]
            fidMatch[i,1]=-1
            fidMatch[i,2]=-1
            fidMatch[i,3]=-1
    
    return cobraMatch,fidMatch,centroids

def runMatchDF(frameID,mcsData,expData,cameraName,fids,centrePos,armLength,dotPos,goodIdx,adjacentCobras,targetPos):


    # get fits file information and data

    altitude = expData['altitude'].values[0]
    insrot = expData['insrot'].values[0]

    # the format the matching wants
    
    centroids = np.array([mcsData['spot_id'].values,mcsData['mcs_center_x_pix'].values,mcsData['mcs_center_y_pix']]).T

    # load the transformations
    pfiTransform = transformUtils.fromCameraName(cameraName, altitude=altitude, insrot=insrot)

    fidList=list(fids['fiducialId'].values)
    badFid = [1,32,34,61,68,75,88,89,2,4,33,36,37,65,66,67,68,69]

    for i in badFid:
        try:
            fidList.remove(i)
        except:
                pass
    outerRing = np.zeros(len(fids), dtype=bool)   
    goodFid = np.zeros(len(fids), dtype=bool)   
    for i in [29, 30, 31, 61, 62, 64, 93, 94, 95, 96]:
        outerRing[fids.fiducialId == i] = True
    for i in fidList:
        goodFid[fids.fiducialId == i] = True

    # match fiducials
    ffids0,dist0=pfiTransform.updateTransform(mcsData, fids[outerRing], matchRadius=8.0, nMatchMin=0.1)
    ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=4,nMatchMin=0.1)

    # transform
    
    mmCentroids=np.copy(centroids)
    mmCentroids[:,1], mmCentroids[:,2] = pfiTransform.mcsToPfi(centroids[:,1],centroids[:,2])


    
    # identify
    cobraMatch, unaPoints = mcsTools.fibreId(mmCentroids, centrePos, armLength, targetPos, fids, dotPos, goodIdx, adjacentCobras)
    fidPos =  np.array([fids['fiducialId'],fids['x_mm'],fids['y_mm']]).T
    fidMatch=np.zeros((len(fidPos),4))
    # get the matching for fiducials
    D = cdist(fidPos[:,1:3],mmCentroids[:,1:3])
    
    for i in range(len(fidPos)):
        ind = np.where(D[i, :] < 1)
        if len(ind[0]) > 0:
            fidMatch[i,0]=fidPos[i,0]
            fidMatch[i,1]=ind[0][0]
            fidMatch[i,2]=mmCentroids[ind[0][0],1]
            fidMatch[i,3]=mmCentroids[ind[0][0],2]
        else:
            fidMatch[i,0]=fidPos[i,0]
            fidMatch[i,1]=-1
            fidMatch[i,2]=-1
            fidMatch[i,3]=-1
    
    return cobraMatch,fidMatch,centroids
