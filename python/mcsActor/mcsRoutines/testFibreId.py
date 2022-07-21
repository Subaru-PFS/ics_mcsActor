

import mcsActor.windowedCentroid.centroid as centroid
import mcsActor.mcsRoutines.mcsRoutines as mcsTools
import mcsActor.mcsRoutines.dbRoutinesMCS as dbTools
from astropy.io import fits
import matplotlib.pylab as plt
from importlib import reload
from pfs.utils import butler
import pandas as pd
#import pfs.utils.coordinates.transform as transformUtils
import numpy as np
import sys
sys.path.append("/Users/karr/software/mhs/products/Darwin/pfs_utils/6.1.5/python/pfs.utils-0.0.0-py3.8.egg/pfs/utils/coordinates/")
sys.path.append("/Users/karr/software/mhs/products/Darwin/pfs_utils/6.1.5/python/pfs.utils-0.0.0-py3.8.egg/pfs/utils/")
import transform as transformUtils
import butler


def initFibreID(xmlFile,dotFile):
    centParms = mcsTools.getCentroidParams(None)
    centParms['nmax']=1e6

    centrePos, armLength, dotPos, goodIdx, des = mcsTools.readCobraGeometry(xmlFile,dotFile)
    adjacentCobras = mcsTools.makeAdjacentList(centrePos, armLength)
    fids = pd.read_csv("/Users/karr/Science/PFS/newestMHS/pfs_instdata/data/pfi/fiducial_positions.csv", comment='#')

    fidPos =  np.array([fids['fiducialId'],fids['x_mm'],fids['y_mm']]).T

    return centParms, centrePos, armLength, dotPos, goodIdx, des, adjacentCobras, fids, fidPos


def checkTransform(dataPath,frameIds,camera_name,centParms, centrePos, armLength, dotPos, goodIdx, des, adjacentCobras, fids, fidPos):

    centrePosPix = [5033,3447]


    #centrePosPix = np.array([[5033,5033],[5033,5033],[3447,3447]]).T
    tarPos = centrePos
    prevPos = centrePos
    iter = 0
    #orphans = []
    for frameId in frameIds:
        print("frameId = ",frameId)
        fName=dataPath+"PFSC0"+str(frameId)+".fits"
        image=fits.getdata(fName)
        if(iter==0):
            findThresh, centThresh, avBack = mcsTools.getThresh(
            image, centrePosPix, centParms['threshSigma'], centParms['findSigma'], centParms['centSigma'])
            iter = 1

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
        mcsFrameId = frameId
        mcsFrameIDs = np.repeat(mcsFrameId, sz[0]).astype('int')
        # make a data frame
        frame = np.zeros((sz[0], 9))
        frame[:, 0] = mcsFrameIDs
        frame[:, 1:] = points
        # column names
        columns = ['mcs_frame_id', 'spot_id', 'mcs_center_x_pix', 'mcs_center_y_pix', 'mcs_second_moment_x_pix',
            'mcs_second_moment_y_pix', 'peakvalue', 'bgvalue', 'mcs_second_moment_xy_pix']
        mcsData = pd.DataFrame(frame, columns=columns)
        centroids = points
        
        pfiTransform = transformUtils.fromCameraName(cameraName, altitude=60, insrot=0)
        fids = pd.read_csv("/Users/karr/Science/PFS/newestMHS/pfs_instdata/data/pfi/fiducial_positions.csv", comment='#')

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
        
        ffids0,dist0=pfiTransform.updateTransform(mcsData, fids[outerRing], matchRadius=8.0, nMatchMin=0.1)
        ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=4.2,nMatchMin=0.1)
        ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=2,nMatchMin=0.1)
        #ffids2,dist2=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=3.0,nMatchMin=0.1)
        pfiTrans = pfiTransform

        centroids[:,1], centroids[:,2] = pfiTransform.mcsToPfi(centroids[:,1],centroids[:,2])

  
        xd=[]
        yd=[]
        xp=[]
        yp=[]
        d=[]

        xx,yy =  pfiTransform.mcsToPfi(mcsData['mcs_center_x_pix'].values,mcsData['mcs_center_y_pix'].values)
        for x,y in zip(fids['x_mm'].values,fids['y_mm'].values):
            dd=np.sqrt((x-xx)**2+(y-yy)**2).ravel()
            ind=np.argmin(dd)
            xd.append(x-xx[ind])
            yd.append(y-yy[ind])
            xp.append(x)
            yp.append(y)
            d.append(dd[ind])
        
        fig,ax=plt.subplots()
        keyLength = 0.1
        q = ax.quiver(xp,yp,xd,yd)
        fTitle = str(keyLength)+" mm"
        ax.quiverkey(q,-200,200,keyLength,fTitle,coordinates='data',color='red',labelcolor='red')
        ax.set_ylabel("y")
        ax.set_xlabel("x")
        ax.set_aspect('equal')
        plt.savefig("v2.png")

        fig,ax=plt.subplots()
        ax.hist(d,bins=np.arange(0,2,0.1))


def testFibreIDCSV(dataPath,frameIds,cameraName, df, centParms, centrePos, armLength, dotPos, goodIdx, des, adjacentCobras, fids, fidPos, calibModel):
                   

    targets = centrePos
    tarPos = centrePos
    prevPos = centrePos
    iter = 0
    #orphans = []
    #fig1,ax1=plt.subplots()
    #fig2,ax2=plt.subplots()
    for frameId in frameIds:
        
        fName=dataPath+"PFSC0"+str(frameId)+".fits"
        hdu = fits.open(fName)
        hdr = hdu[0].header
        altitude = hdr['altitude']
        insrot = hdr['insrot']


        mcsData=df.loc[df['mcs_frame_id']==frameId]

        centroids = np.array([mcsData['spot_id'].values,mcsData['mcs_center_x_pix'],mcsData['mcs_center_y_pix']]).T
        pfiTransform = transformUtils.fromCameraName(cameraName, altitude=altitude, insrot=insrot)
        #fix,ax=plt.subplots()
        
        #ax.scatter(centroids[:,1],centroids[:,2])
        #plt.savefig(str(int(frameId))+".png")

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
        
        ffids0,dist0=pfiTransform.updateTransform(mcsData, fids[outerRing], matchRadius=8.0, nMatchMin=0.1)
        ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=4.2,nMatchMin=0.1)
        #ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=2,nMatchMin=0.1)
        #ffids2,dist2=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=3.0,nMatchMin=0.1)
        pfiTrans = pfiTransform

        mmCentroids=np.copy(centroids)
        mmCentroids[:,1], mmCentroids[:,2] = pfiTransform.mcsToPfi(centroids[:,1],centroids[:,2])
        cobraMatch, unaPoints = mcsTools.fibreId(mmCentroids, centrePos, armLength, prevPos, fids, dotPos, goodIdx, adjacentCobras)
 
        #fig3,ax3=plt.subplots()
        #for i in range(len(goodIdx)):
        #    ax1.plot([cobraMatch[i,2],centrePos[i,1]],[cobraMatch[i,3],centrePos[i,2]])
        #ax2.scatter(mmCentroids[:,1],mmCentroids[:,2],color="tab:blue")
    
        #ax3.scatter(centroids[:,1],centroids[:,2])
        #diff = np.sqrt((cobraMatch[:,2]-c[:,1])**2+(cobraMatch[:,3]-prevPos[:,2])**2)
        ##ax3.scatter(centroids[unaPoints,1],centroids[unaPoints,2])
        #
        ##ax1.hist(diff.ravel())
        #ind = np.where(diff > 4)
        #print("D",cobraMatch[ind,1])
        prevPos = cobraMatch[:,[0,2,3]]
        #orphans.append(centroids[unaPoints,:])

        #ax1.scatter(centrePos[:,1],centrePos[:,2])
        #ax1.scatter(fids['x_mm'],fids['y_mm'],marker="d")
        #ax2.scatter(centrePos[:,1],centrePos[:,2],color="black")
        #ax2.scatter(fids['x_mm'],fids['y_mm'],marker="d",color="tab:orange")
        ##ax3.scatter(centrePos[:,1],centrePos[:,2])
        #
        #for i in range(len(centrePos)):
        #    circle=plt.Circle((centrePos[i,1],centrePos[i,2]),armLength[i],fill=False,color='black')
        #    circle1=plt.Circle((centrePos[i,1],centrePos[i,2]),armLength[i],fill=False,color='black')
        #    circle3=plt.Circle((centrePos[i,1],centrePos[i,2]),armLength[i],fill=False,color='black')
        #    a=ax1.add_artist(circle)
        #    a=ax2.add_artist(circle1)
        #    #a=ax3.add_artist(circle3)

    return cobraMatch,centroids, mmCentroids, unaPoints

def testFibreIDNearest(dataPath,frameIds,cameraName,centParms, centrePos, armLength, dotPos, goodIdx, des, adjacentCobras, fids, fidPos,calibModel):


    cobraMatches=[]
    allCentroids=[]
    centrePosPix = np.array([5033,3447]).T
    tarPos = centrePos
    prevPos = centrePos
    iter = 0
    #orphans = []
    fig1,ax1=plt.subplots()
    fig2,ax2=plt.subplots()
    for frameId in frameIds:
    
        print("frameId = ",frameId)
        fName=dataPath+"PFSC0"+str(frameId)+".fits"
        hdu = fits.open(fName)
        hdr = hdu[0].header
        altitude = hdr['altitude']
        insrot = hdr['insrot']
        image=fits.getdata(fName)
        if(iter==0):
            findThresh, centThresh, avBack = mcsTools.getThresh(
            image, centrePosPix,  centParms['threshSigma'], centParms['findSigma'], centParms['centSigma'])
            iter = 1

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
        mcsFrameId = frameId
        mcsFrameIDs = np.repeat(mcsFrameId, sz[0]).astype('int')
        # make a data frame
        frame = np.zeros((sz[0], 9))
        frame[:, 0] = mcsFrameIDs
        frame[:, 1:] = points
        # column names
        columns = ['mcs_frame_id', 'spot_id', 'mcs_center_x_pix', 'mcs_center_y_pix', 'mcs_second_moment_x_pix',
            'mcs_second_moment_y_pix', 'peakvalue', 'bgvalue', 'mcs_second_moment_xy_pix']
        mcsData = pd.DataFrame(frame, columns=columns)
        centroids = points
        print(altitude,insrot)
        pfiTransform = transformUtils.fromCameraName(cameraName, altitude=altitude, insrot=insrot)
        fids = pd.read_csv("/Users/karr/Science/PFS/newestMHS/pfs_instdata/data/pfi/fiducial_positions.csv", comment='#')

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
        
        ffids0,dist0=pfiTransform.updateTransform(mcsData, fids[outerRing], matchRadius=8.0, nMatchMin=0.1)
        ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=4.2,nMatchMin=0.1)
        ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=2,nMatchMin=0.1)
        #ffids2,dist2=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=3.0,nMatchMin=0.1)
        pfiTrans = pfiTransform

        centroids[:,1], centroids[:,2] = pfiTransform.mcsToPfi(centroids[:,1],centroids[:,2])
        cobraMatch, unaPoints = mcsTools.fibreId(centroids, centrePos, armLength, prevPos, fids, dotPos, goodIdx, adjacentCobras)
        targets = cobraMatch[:,[0,2,3]]

        #fig3,ax3=plt.subplots()
        for i in range(len(goodIdx)):
            ax1.plot([cobraMatch[i,2],calibModel.centers.real[goodIdx[i]]],[cobraMatch[i,3],calibModel.centers.imag[goodIdx[i]]])
        ax2.scatter(centroids[:,1],centroids[:,2],color="tab:blue")
    
        #ax3.scatter(centroids[:,1],centroids[:,2])
        #diff = np.sqrt((cobraMatch[:,2]-centrePos[:,1])**2+(cobraMatch[:,3]-centrePos[:,2])**2)
        ##ax3.scatter(centroids[unaPoints,1],centroids[unaPoints,2])
        #
        ##ax1.hist(diff.ravel())
        #ind = np.where(diff > 4)
        #print("D",cobraMatch[ind,1])
        prevPos = cobraMatch[:,[0,2,3]]
        #orphans.append(centroids[unaPoints,:])

        ax1.scatter(centrePos[:,1],centrePos[:,2])
        ax1.scatter(fids['x_mm'],fids['y_mm'],marker="d")
        ax2.scatter(centrePos[:,1],centrePos[:,2],color="black")
        ax2.scatter(fids['x_mm'],fids['y_mm'],marker="d",color="tab:orange")
        #ax3.scatter(centrePos[:,1],centrePos[:,2])

        for i in range(len(centrePos)):
            circle=plt.Circle((centrePos[i,1],centrePos[i,2]),armLength[i],fill=False,color='black')
            circle1=plt.Circle((centrePos[i,1],centrePos[i,2]),armLength[i],fill=False,color='black')
            circle3=plt.Circle((centrePos[i,1],centrePos[i,2]),armLength[i],fill=False,color='black')
            a=ax1.add_artist(circle)
            a=ax2.add_artist(circle1)
            #a=ax3.add_artist(circle3)

        cobraMatches.append(cobraMatch)
        allCentroids.append(centroids)
    return cobraMatches,allCentroids
