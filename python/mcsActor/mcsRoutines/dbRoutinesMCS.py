"""
Routines to read and write information to DB for ics_mcsActor

Most of this is wrappers to the opdb routines, creating the SQL query,
or turning the data into the correct pandas dataframe for insertion.

tables involved

mcs_data (write)
mcs_exposure (write, read)
mcs_boresight (read, write in calibration mode)
mcs_pfi_transformation (write)
cobra_target (read)
cobra_match (write)
fiducial_fiber_geometry (read)

pfs_visit (must be populated)
fiducial_fiber (must be populated)


"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import sys

import pandas as pd
from scipy.stats import sigmaclip
import copy
from opdb import opdb


def connectToDB(hostname='',port='',dbname='opdb_api_test',username='pfs',passwd=None):
    
    db = opdb.OpDB(hostname, port, dbname, username, passwd)
    db.connect()
    
    return db

def loadCobraMatchFromDB(db,frameId):

    """
    read the cobra_match table information, convert to x+ij format, return positions and flags
    """

    
    sql=f'SELECT cobra_match.cobra_id,cobra_match.pfi_center_x_mm,cobra_match.pfi_center_x_mm,cobra_match.flags from cobra_match where cobra_match.mcs_frame_id={frameId}'

    df=db.fetch_query(sql)

    positions=df['pfi_center_x_mm']+df['pfi_center_y_mm']*1j
    flags=df['flags']

    return positions,flags


def loadTelescopeParametersFromDB(db,frameId):

    sql=f'SELECT mcs_exposure.insrot,mcs_exposure.altitude FROM mcs_exposure WHERE mcs_exposure.mcs_frame_id={frameId}'
    df=db.fetch_query(sql)
    zenithAngle=90-df['altitude'][0]
    insRot=df['insrot'][0]

    return zenithAngle,insRot

def loadBoresightFromDB(db,pfsVisitId):
    
    """
    read boresight informatino from database
    table = mcs_boresight
    """


    sql=f'SELECT mcs_boresight.mcs_boresight_x_pix,mcs_boresight.mcs_boresight_y_pix FROM mcs_boresight ORDER BY calculated_at DESC LIMIT 1'
    df=db.fetch_query(sql)
    return [df['mcs_boresight_x_pix'][0],df['mcs_boresight_y_pix'][0]]


def loadFiducialsFromDB(db):

    """
    load fiducial fibre positions from the DB
    table=fiducial_fiber_geometry

    returns an  Nx3 array with the columns equal to the fiducial ID, x and y position
    """

    sql=f'SELECT fiducial_fiber_geometry.fiducial_fiber_id , fiducial_fiber_geometry.ff_center_on_pfi_x_mm ,fiducial_fiber_geometry.ff_center_on_pfi_y_mm from fiducial_fiber_geometry'

    df=db.fetch_query(sql)

    
    return np.array([df['fiducial_fiber_id'],df['ff_center_on_pfi_x_mm'],df['ff_center_on_pfi_y_mm']]).T
    

def loadTargetsFromDB(db,frameId):

    """
    load the intermediate target positions from the database
    table=cobra_target
    
    returns an Nx3 array with the columns equal to the cobra ID, x and y position
    """

    visitId = frameId // 100
    iteration = frameId % 100
    
    sql=f'SELECT cobra_target.cobra_id,cobra_target.pfi_nominal_x_mm,cobra_target.pfi_nominal_y_mm FROM cobra_target WHERE cobra_target.iteration={iteration} and cobra_target.pfs_visit_id={visitId}'
    df=db.fetch_query(sql)

    return np.array([df['cobra_id'],df['pfi_nominal_x_mm'],df['pfi_nominal_y_mm']]).T


def writeCentroidsToDB(db,centroids,mcsFrameId):

    """
    write the centroids to the database
    table=mcs_data
    variables=spot_id,mcs_center_x_pix,mcs_center_y_pix
              mcs_second_moment_x_pix,mcs_second_moment_y_pix,
              mcs_second_moment_xy_pix,bgvalue,peakvalue
    """


    #get size of array
    sz=centroids.shape
    
    #create array of frameIDs (same for all spots)
    mcsFrameIDs=np.repeat(mcsFrameId,sz[0]).astype('int')
    #make a data frame
    frame=np.zeros((sz[0],9))
    frame[:,0]=mcsFrameIDs
    frame[:,1:]=centroids
    frame[:,8]=np.zeros((sz[0]))
    #column names
    columns=['mcs_frame_id','spot_id','mcs_center_x_pix','mcs_center_y_pix','mcs_second_moment_x_pix','mcs_second_moment_y_pix','mcs_second_moment_xy_pix','bgvalue','peakvalue']
    df=pd.DataFrame(frame,columns=columns)
    db.insert("mcs_data",df)

def writeMatchesToDB(db,cobraMatch,mcsFrameId):

    """
    write the centroids to the database
    table=mcs_data
    variables=spot_id,mcs_center_x_pix,mcs_center_y_pix
              mcs_second_moment_x_pix,mcs_second_moment_y_pix,
              mcs_second_moment_xy_pix,bgvalue,peakvalue
    """

    np.save("cobraMatch.npy",cobraMatch)
    #get size of array
    sz=cobraMatch.shape
    #create array of frameIDs (same for all spots)
    mcsFrameIDs=np.repeat(mcsFrameId,sz[0]).astype('int')
    visitIds=np.repeat(mcsFrameId // 100,sz[0]).astype('int')
    iterations=np.repeat(mcsFrameId % 100,sz[0]).astype('int')
    #make a data frame
    frame=np.zeros((sz[0],8))
    frame[:,2]=mcsFrameIDs
    frame[:,0]=visitIds
    frame[:,1]=iterations
    frame[:,3:]=cobraMatch

    ind=np.where(cobraMatch[:,1]!=-1)
    np.save("frame.npy",frame)
    #column names
    columns=['pfs_visit_id','iteration','mcs_frame_id','cobra_id','spot_id','pfi_center_x_mm','pfi_center_y_mm','flags']
    df=pd.DataFrame(frame[ind[0],:],columns=columns)
    db.insert("cobra_match",df)

    columns=['pfs_visit_id','iteration','mcs_frame_id','cobra_id','pfi_center_x_mm','pfi_center_y_mm','flags']
    ind=np.where(cobraMatch[:,1]==-1)
    if(len(ind[0])>0):
       ff=frame[ind[0],:]
       df=pd.DataFrame(ff[:,[0,1,2,3,5,6,7]],columns=columns)

       db.insert("cobra_match",df)


def writeAffineToDB(db,afCoeff,frameId):

    """
    write the affine transformation to DB
    """

    sx=np.sqrt(afCoeff[0,0]**2+afCoeff[0,1]**2)
    sy=np.sqrt(afCoeff[1,0]**2+afCoeff[1,1]**2)
        
    xd=afCoeff[0,2]
    yd=afCoeff[1,2]
    
    rotation=np.arctan2(afCoeff[1,0]/np.sqrt(afCoeff[0,0]**2+afCoeff[0,1]**2),
                        afCoeff[1,1]/np.sqrt(afCoeff[1,0]**2+afCoeff[1,1]**2))

    df = pd.DataFrame({'mcs_frame_id':[frameId],'x_trans':[xd],'y_trans':[yd],'x_scale':[sx],'y_scale':[sy],'angle':[rotation]})
    db.insert('mcs_pfi_transformation',df)
