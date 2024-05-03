
import pathlib
from importlib import reload 
import matplotlib.pyplot as plt
import numpy as np
import os, fnmatch
import astropy.io.fits as pyfits
import pandas as pd
from opdb import opdb
from photutils.detection import DAOStarFinder
import time
from ics.cobraCharmer.cobraCoach import visDianosticPlot

reload(visDianosticPlot)

import psycopg2
from sqlalchemy import create_engine

from pfs.utils.fiberids import FiberIds
from pfs.utils import butler 

import pfs.utils.coordinates.transform as transformUtils
from pfs.utils.coordinates.transform import PfiTransform
from pfs.utils.coordinates.transform import matchIds

import yaml

from mcsActor.mcsRoutines import speedCentroid
import mcsActor.windowedCentroid.centroid as centroid

newXml = pathlib.Path('/software/mhs/products/Linux64/pfs_instdata/1.7.36/data/pfi/modules/ALL/ALL.xml')


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

def solveForTransform(mcsData, fids, camera, insrot, altitude = None,
                      nsigma=0, alphaRot=2, makePlots=False):
    
    if camera == 'canon50M':
    
        pt = transformUtils.fromCameraName('canon50M',
            altitude=altitude,
            insrot=insrot)
    else:
        pt = transformUtils.fromCameraName(camera, altitude=90, insrot=insrot, nsigma=0, alphaRot=alphaRot)
    
    outerRingIds, fidsGood = readFiducialMasks(fids)
    #outerRingIds = [29, 30, 31, 61, 62, 64, 93, 94, 95, 96]
    #fidsOuterRing = fids[fids.fiducialId.isin(outerRingIds)]

    fig = None
    ifig = 1  # first figure to use
    if makePlots:
        fig = plt.figure(ifig); plt.clf()
    pt.updateTransform(mcsData, outerRingIds, matchRadius=10.0, nMatchMin=0.1, fig=fig)
    pt.mcsDistort.setArgs(pt.mcsDistort.getArgs())

    pt.nsigma = nsigma
    pt.alphaRot = 0

    for i in range(2):
        if makePlots:
            fig = plt.figure(ifig); plt.clf()
        fid, dmin = pt.updateTransform(mcsData, fids, matchRadius=4.2, nMatchMin=0.1, fig=fig)
        
    return pt, fid, dmin


def getMCSdata(frameNum):
    fids = butler.Butler().get('fiducials')

    conn_str = "postgresql://pfs@db-ics/opdb"
    engine = create_engine(conn_str)

    with engine.connect() as conn:
        mcsData = pd.read_sql(f'''
            SELECT DISTINCT 
                spot_id, mcs_center_x_pix, mcs_center_y_pix
            FROM mcs_data
            WHERE
              mcs_frame_id = {frameNum}
            -- limit 10
            ''', conn)

    mcs_x = mcsData['mcs_center_x_pix'].to_numpy()          
    mcs_y = mcsData['mcs_center_y_pix'].to_numpy()

    x_fid_mm = fids.x_mm.to_numpy()
    y_fid_mm = fids.y_mm.to_numpy()
    fiducialId = fids.fiducialId.to_numpy()

    return mcsData, fids

def diffMCS(visit0, visit1):
    conn = psycopg2.connect("dbname='opdb' host='db-ics' port=5432 user='pfs'") 
    engine = create_engine('postgresql+psycopg2://', creator=lambda: conn)


    mcsData0 = pd.read_sql(f'''
        SELECT DISTINCT 
            *
        FROM cobra_match
        WHERE
            cobra_match.mcs_frame_id = %(pfs_visit_id)s
        ORDER BY
            cobra_match.cobra_id ASC
        --limit 10
    ''', engine, params={'pfs_visit_id': visit0})
    
    mcsData1 = pd.read_sql(f'''
        SELECT DISTINCT 
            *
        FROM cobra_match
        WHERE
            cobra_match.mcs_frame_id = %(pfs_visit_id)s
        ORDER BY
            cobra_match.cobra_id ASC
        --limit 10
    ''', engine, params={'pfs_visit_id': visit1})

    xx=mcsData0['pfi_center_x_mm']- mcsData1['pfi_center_x_mm']
    yy=mcsData0['pfi_center_y_mm']- mcsData1['pfi_center_y_mm']

    return xx.values,yy.values


visit0 = 10749000
visit1 = 10749001

def visCobraMatchDiff(visit0, visit1):
    xx,yy, = diffMCS(visit0, visit1)

    vis=visDianosticPlot.VisDianosticPlot(xml=newXml)
    vis.visCreateNewPlot(f'Cobra Convergence {visit0} - {visit1}', 
                                'x', 'Y',size=(6, 6))
    ax=plt.gca()
    #ax.set_aspect('auto')
    x = vis.calibModel.centers.real
    y = vis.calibModel.centers.imag

    vectorLength = 0.01
    q=ax.quiver(x,y, 
        xx, yy, color='red',units='xy')


    ax.quiverkey(q, X=0.2, Y=0.95, U=vectorLength,
        label=f'length = {vectorLength} mm', labelpos='E')

    plt.show()

def getCobraPixel(visit0, cobraIdx):
    mcsData, fids = getMCSdata(visit0)
    vis=visDianosticPlot.VisDianosticPlot(xml=newXml)

    db=opdb.OpDB(hostname='db-ics', port=5432,dbname='opdb',
                        username='pfs')
    teleInfo = db.bulkSelect('mcs_exposure','select altitude, insrot from mcs_exposure where '
            f'mcs_frame_id = {visit0}')
    cobraMatch = db.bulkSelect('cobra_match','select cobra_id, pfi_center_x_mm, pfi_center_y_mm from cobra_match where '
            f'mcs_frame_id = {visit0}')
    db.close()

    pt, fid, dmin =solveForTransform(mcsData, fids,'canon50M', teleInfo['insrot'].values[0],
                                    altitude=teleInfo['altitude'].values[0], nsigma=0, alphaRot=0, makePlots=False)

    xcenter, ycenter = pt.pfiToMcs([vis.calibModel.centers[cobraIdx].real],
                                [vis.calibModel.centers[cobraIdx].imag]) 
    
    return xcenter, ycenter

def runCentroidDAO(visit):
    image = pyfits.open(pathlib.Path(f'/data/MCS/{visDianosticPlot.findRunDir(int(visit/100))}/data/PFSC{visit}.fits'))
    data = image[1].data
    m, s = np.mean(data), np.std(data)

    t0 = time.time()
    daofind = DAOStarFinder(fwhm=3.0, threshold=7.*s)  
    centroids = daofind(data - m) 
    t1 = time.time()

    print(f'Time = {t1 - t0}')

    return centroids

def runCentroidSEP(visit):
    image = pyfits.open(pathlib.Path(f'/data/MCS/{visDianosticPlot.findRunDir(int(visit/100))}/data/PFSC{visit}.fits'))

    t0 = time.time()
    spCenMT = speedCentroid.speedCentroid(image[1].data)
    spCenMT.cores = 8
    spCenMT.runCentroidMP()
    spCenMT.arrangeCentroid()
    t1 = time.time()

    print(f'Time = {t1 - t0}')

    return spCenMT.centroids
    
def runCentroid(visit):
    image = pyfits.open(pathlib.Path(f'/data/MCS/{visDianosticPlot.findRunDir(int(visit/100))}/data/PFSC{visit}.fits'))

    centParms= {'boxCent': 6, 'boxFind': 10, 'centSigma': 15, 'findSigma': 25,
        'fwhmx': 2, 'fwhmy': 2, 'matchRad': 20, 'maxIt': 20, 'nmax': 180, 'nmin': 8, 
        'threshFact': 4, 'threshMode': 'full', 'threshSigma': 4,
        'activeX': 2665, 'activeY': 2665}
    
    #findThresh =1975.7806658274599
    findThresh =3002.7806658274599
    #centThresh = 1902.3914938148669
    centThresh = 3002.3914938148669
    t0 = time.time()
    a = centroid.centroid_only(image[1].data.astype('<i4'),
           centParms['fwhmx'], centParms['fwhmy'], findThresh, centThresh,
           centParms['boxFind'], centParms['boxCent'],
           centParms['nmin'], centParms['maxIt'], 0)

    centroids = np.frombuffer(a, dtype='<f8')
    centroids = np.reshape(centroids, (len(centroids)//7, 7))
    t1 = time.time()

    print(f'Time = {t1 - t0}')

    return centroids


def showMcsDataImage(visit0, visit1, cobraIdx, 
                    centroid0, centroid1, boxSize=10, method='win'):
    image0 = pyfits.open(pathlib.Path(f'/data/MCS/{visDianosticPlot.findRunDir(int(visit0/100))}/data/PFSC{visit0}.fits'))
    image1 = pyfits.open(pathlib.Path(f'/data/MCS/{visDianosticPlot.findRunDir(int(visit1/100))}/data/PFSC{visit1}.fits'))


    mcsData, fids = getMCSdata(visit0)
    mcsData2, fids = getMCSdata(visit1)
    
    xcenter, ycenter = getCobraPixel(visit0, cobraIdx)

    vis=visDianosticPlot.VisDianosticPlot(xml=newXml)
    data = image0[1].data-image1[1].data
    data = image0[1].data
    m, s = np.mean(data), np.std(data)

    vis.visCreateNewPlot(f'Image  {visit0}', 
                             'X', 'Y',size=(8, 8))

    ax = plt.gca()

    im = ax.imshow(data, interpolation='nearest', 
                cmap='gray', vmin=m-s, vmax=m+3*s, origin='lower')
    ax.scatter(mcsData['mcs_center_x_pix'].values,mcsData['mcs_center_y_pix'].values,label=f'DB {visit0}')
    ax.scatter(mcsData2['mcs_center_x_pix'].values,mcsData2['mcs_center_y_pix'].values,label=f'DB {visit1}')

    if method == 'win':
        ax.scatter(centroid0[:,0],centroid0[:,1],marker='x',label=f'win {visit0}')
        ax.scatter(centroid1[:,0],centroid1[:,1],marker='x',label=f'win {visit1}')

    if method == 'sep':
    #ax.scatter(centroids[:,0],centroids[:,1],marker='x',label=f'Manual')
        ax.scatter(centroid0['x'],centroid0['y'],marker='x',label=f'sep {visit0}')
        ax.scatter(centroid1['x'],centroid1['y'],marker='x',label=f'sep {visit1}')

    if method == 'dao':
        ax.scatter(centroid0['xcentroid'].value,centroid0['ycentroid'].value,marker='x',label=f'dao {visit0}')
        ax.scatter(centroid1['xcentroid'].value,centroid1['ycentroid'].value,marker='x',label=f'dao {visit1}')
    ax.legend()

    vis.visSetAxesLimits([xcenter-boxSize,xcenter+boxSize],[ycenter-boxSize,ycenter+boxSize])

def main():
    runCentroid(10749001)


if __name__:
    main()