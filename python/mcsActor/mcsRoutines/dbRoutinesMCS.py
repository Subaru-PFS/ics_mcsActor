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
import io 
import logging

import pandas as pd
from scipy.stats import sigmaclip
import copy
from opdb import opdb
from datetime import datetime, timezone

def connectToDB(hostname='', port='', dbname='opdb', username='pfs'):

    db = opdb.OpDB(hostname=hostname, port=port,
                   dbname=dbname,
                   username=username)
    db.connect()

    return db


def loadCobraMatchFromDB(db, frameId):
    """
    read the cobra_match table information, convert to x+ij format, return positions and flags
    """

    sql = f'SELECT cobra_match.cobra_id,cobra_match.pfi_center_x_mm,cobra_match.pfi_center_x_mm,cobra_match.flags from cobra_match where cobra_match.mcs_frame_id={frameId}'

    df = db.fetch_query(sql)

    positions = df['pfi_center_x_mm']+df['pfi_center_y_mm']*1j
    flags = df['flags']

    return positions, flags


def loadTelescopeParametersFromDB(db, frameId):

    sql = f'SELECT mcs_exposure.insrot,mcs_exposure.altitude FROM mcs_exposure WHERE mcs_exposure.mcs_frame_id={frameId}'
    df = db.fetch_query(sql)

    if df['altitude'][0] < -99:
        zenithAngle = 90
    else:
        zenithAngle = 90-df['altitude'][0]

    insRot = df['insrot'][0]

    return zenithAngle, insRot
    
def loadBoresightFromDB(db, pfsVisitId):
    """
    read boresight informatino from database
    table = mcs_boresight
    """
    sql = f'''SELECT * FROM mcs_boresight ORDER BY calculated_at DESC FETCH FIRST ROW ONLY'''
    #sql=f'SELECT mcs_boresight.mcs_boresight_x_pix,mcs_boresight.mcs_boresight_y_pix from mcs_boresight where mcs_boresight.pfs_visit_id={pfsVisitId}'
    df = db.fetch_query(sql)
    return [df['mcs_boresight_x_pix'][0], df['mcs_boresight_y_pix'][0]]


def loadCentroidsFromDB(db, mcsFrameId):
    """ retrieve a set of centroids from database and return as a numpy array"""
    
    sql = f'select mcs_data.spot_id, mcs_data.mcs_center_x_pix, mcs_data.mcs_center_y_pix from mcs_data where mcs_data.mcs_frame_id={mcsFrameId}'
    df = db.fetch_query(sql)
    return df.to_numpy()
    
def loadFiducialsFromDB(db):
    """
    load fiducial fibre positions from the DB
    table=fiducial_fiber_geometry

    returns an  Nx3 array with the columns equal to the fiducial ID, x and y position
    """

    sql = f'SELECT fiducial_fiber_geometry.fiducial_fiber_id , fiducial_fiber_geometry.ff_center_on_pfi_x_mm ,fiducial_fiber_geometry.ff_center_on_pfi_y_mm from fiducial_fiber_geometry'

    df = db.fetch_query(sql)

    return np.array([df['fiducial_fiber_id'], df['ff_center_on_pfi_x_mm'], df['ff_center_on_pfi_y_mm']]).T


def loadTargetsFromDB(db, frameId):
    """
    load the intermediate target positions from the database
    table=cobra_target

    returns an Nx3 array with the columns equal to the cobra ID, x and y position
    """

    visitId = frameId // 100
    iteration = frameId % 100

    sql = f'SELECT cobra_target.cobra_id,cobra_target.pfi_nominal_x_mm,cobra_target.pfi_nominal_y_mm FROM cobra_target WHERE cobra_target.iteration={iteration} and cobra_target.pfs_visit_id={visitId}'
    df = db.bulkSelect('cobra_target',sql)

    return np.array([df['cobra_id'], df['pfi_nominal_x_mm'], df['pfi_nominal_y_mm']]).T

def writeTransformToDB(db, frameId, pfiTransform, cameraName):
    """
    write transformation coefficients to database
    """
    trans=pfiTransform.mcsDistort.getArgs()
    res = db.session.execute('select * FROM "mcs_pfi_transformation" where false')
    colnames = tuple(res.keys())
    realcolnames = colnames[0:]
    line = '%d,%f,%f,%f,%e,%e,%f,%s' % (frameId, trans[0].astype('float64'),
           trans[1], trans[2], trans[3],trans[4],
           pfiTransform.alphaRot, 'canon50M')
                                                                        
    buf = io.StringIO()
    buf.write(line)
    buf.seek(0, 0)
    _writeData('mcs_pfi_transformation', realcolnames, buf)


    data = {'mcs_frame_id': [frameId],
            'x0': [trans[0]],
            'y0': [trans[1]],
            'dscale': [trans[2]],
            'scale2': [trans[3]],
            'theta': [trans[4]],
            'alpha_rot': [pfiTransform.alphaRot],
            'camera_name': [cameraName]}
    df = pd.DataFrame(data=data)
    return df['mcs_frame_id'].values,df['x0'].values,df['y0'].values,df['dscale'].values, \
        df['scale2'].values,df['theta'].values,df['alpha_rot'].values,df['camera_name'].values
    
    
def _writeData(tableName, columnNames, dataBuf):
    """Wrap a direct COPY_FROM via sqlalchemy. """
    columns = ','.join('"{}"'.format(k) for k in columnNames)
    sql = 'COPY {} ({}) FROM STDIN WITH CSV'.format(
        tableName, columns)
    try:
        db=opdb.OpDB(hostname='db-ics', port=5432,dbname='opdb',
                        username='pfs')
        session = db.session
        with session.connection().connection.cursor() as cursor:
            cursor.copy_expert(sql, dataBuf)
            cursor.close()
        session.execute('commit')
    except Exception as e:
        logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S")
        logger = logging.getLogger('mcscmd')
        logger.setLevel(logging.INFO)
        logger.warn(f"failed to write with {sql}: {e}")

def writeTargetToDB(db, frameId, target, mpos):
    visitId = frameId // 100
    iteration = frameId % 100
    
    # To-Do here we need a better implementation.
    data = {'pfs_visit_id': np.repeat(visitId,2394),     
            'iteration' : np.repeat(iteration,2394),   
            'cobra_id':np.arange(2394)+1,
            
    }

    df = pd.DataFrame(data=data)
    db.insert("cobra_target", df)

def writeFakeTargetToDB(db, centers, frameId):

    """
    make sure there is a target value if the target d
    """
    visitId = frameId // 100
    iteration = frameId % 100

    nCob = 2394
    
    # To-Do here we need a better implementation.
    data = {'pfs_visit_id': np.repeat(visitId,nCob).astype('int'),
            'iteration' : np.repeat(iteration,nCob).astype('int'),
            'cobra_id':np.arange(2394).astype('int')+1,
            'pfs_config_id':np.repeat(0,2394).astype('int'),
            'pfi_nominal_x_mm':centers.real,
            'pfi_nominal_y_mm':centers.imag,
            'pfi_target_x_mm': centers.real,
            'pfi_target_y_mm':centers.imag,
            'cobra_mortor_model_id_theta': np.repeat(0,2394).astype('int'),
            'motor_target_theta':np.repeat(0,2394),
            'motor_num_step_theta':np.repeat(0,2394),
            'motor_on_time_theta':np.repeat(0,2394),
            'cobra_mortor_model_id_phi': np.repeat(0,2394).astype('int'),
            'motor_target_phi':np.repeat(0,2394),
            'motor_num_step_phi':np.repeat(0,2394),
            'motor_on_time_phi':np.repeat(0,2394),
            'flags':np.repeat(0,2394).astype('int')
    }

    df = pd.DataFrame(data=data)
    db.insert("cobra_target", df)

def writeCobraCenterToDB(db, frameId, centers, mpos):
    visitId = frameId // 100
    iteration = frameId % 100

    targetTable = {'pfs_visit_id':np.repeat(visitId,2394),
                    'iteration':np.repeat(iteration,2394),
                    'cobra_id':np.arange(2394)+1,
                    'pfs_config_id':np.repeat(0,2394),
                    'pfi_nominal_x_mm':centers.real,
                    'pfi_nominal_y_mm':centers.imag,
                    'pfi_target_x_mm': mpos.real,
                    'pfi_target_y_mm':mpos.imag,
                    'motor_target_theta':np.repeat(0,2394),
                    'motor_target_phi':np.repeat(0,2394),
            }

    df = pd.DataFrame(data=targetTable)
    db.bulkInsert("cobra_target", df)

def writeBoresightToDB(db, pfsVisitId, boresight):
    """ write boresight to database with current timestamp """
    
    dt = datetime.now(timezone.utc)
   
    df = pd.DataFrame({'pfs_visit_id': [pfsVisitId], 'mcs_boresight_x_pix': [boresight[0]], 'mcs_boresight_y_pix': [boresight[1]],
                       'calculated_at': [dt]})
    db.bulkInsert('mcs_boresight', df)


def writeCentroidsToDB(db, centroids, mcsFrameId):
    """
    write the centroids to the database
    table=mcs_data
    variables=spot_id,mcs_center_x_pix,mcs_center_y_pix
              mcs_second_moment_x_pix,mcs_second_moment_y_pix,
              mcs_second_moment_xy_pix,bgvalue,peakvalue
    """

    # get size of array
    sz = centroids.shape

    # create array of frameIDs (same for all spots)
    mcsFrameIDs = np.repeat(mcsFrameId, sz[0]).astype('int')
    # make a data frame
    frame = np.zeros((sz[0], 9))
    frame[:, 0] = mcsFrameIDs
    frame[:, 1:] = centroids
    # column names
    columns = ['mcs_frame_id', 'spot_id', 'mcs_center_x_pix', 'mcs_center_y_pix', 'mcs_second_moment_x_pix',
               'mcs_second_moment_y_pix', 'peakvalue', 'bgvalue', 'mcs_second_moment_xy_pix']

    df = pd.DataFrame(frame, columns=columns)
    db.insert("mcs_data", df)
    

def readMatchFromDB(db, mcsFrameId):
    
    match = db.bulkSelect('cobra_match',f'select * from cobra_match where '
                      'mcs_frame_id = {mcsFrameId}')
    return match


def writeMatchesToDB(db, cobraMatch, mcsFrameId):
    """
    write the centroids to the database
    table=mcs_data
    variables=spot_id,mcs_center_x_pix,mcs_center_y_pix
              mcs_second_moment_x_pix,mcs_second_moment_y_pix,
              mcs_second_moment_xy_pix,bgvalue,peakvalue
    """

    np.save("cobraMatch.npy", cobraMatch)
    # get size of array
    sz = cobraMatch.shape
    # create array of frameIDs (same for all spots)
    mcsFrameIDs = np.repeat(mcsFrameId, sz[0]).astype('int')
    visitIds = np.repeat(mcsFrameId // 100, sz[0]).astype('int')
    iterations = np.repeat(mcsFrameId % 100, sz[0]).astype('int')
    # make a data frame
    frame = np.zeros((sz[0], 8))
    frame[:, 2] = mcsFrameIDs
    frame[:, 0] = visitIds
    frame[:, 1] = iterations
    frame[:, 3:] = cobraMatch

    columns = ['pfs_visit_id', 'iteration', 'mcs_frame_id', 'cobra_id',
               'spot_id', 'pfi_center_x_mm', 'pfi_center_y_mm', 'flags']

    visitId = mcsFrameId // 100
    iteration = mcsFrameId % 100
    targetTable = {'pfs_visit_id':np.repeat(visitId,2394).astype('int'),
                    'iteration':np.repeat(iteration,2394).astype('int'),
                    'cobra_id':cobraMatch[:,0].astype('int'),
                    'mcs_frame_id':np.repeat(mcsFrameId,2394).astype('int'),
                    'spot_id':cobraMatch[:,1].astype('int'),
                    'pfi_center_x_mm':cobraMatch[:,2],
                    'pfi_center_y_mm':cobraMatch[:,3],
                    'flags':cobraMatch[:,4].astype('int')
                   }
    df = pd.DataFrame(data=targetTable)                   
    db.insert("cobra_match", df)

    #ind = np.where(cobraMatch[:, 1] != -1)
    #np.save("frame.npy", frame)
    # column names
    #columns = ['pfs_visit_id', 'iteration', 'mcs_frame_id', 'cobra_id',
    #           'spot_id', 'pfi_center_x_mm', 'pfi_center_y_mm', 'flags']
    #df = pd.DataFrame(frame[ind[0], :], columns=columns)
    #db.insert("cobra_match", df)

    #columns = ['pfs_visit_id', 'iteration', 'mcs_frame_id',
    #           'cobra_id', 'pfi_center_x_mm', 'pfi_center_y_mm', 'flags']
    #ind = np.where(cobraMatch[:, 1] == -1)
    #if(len(ind[0]) > 0):
    #    ff = frame[ind[0], :]
    #    df = pd.DataFrame(ff[:, [0, 1, 2, 3, 5, 6, 7]], columns=columns)

    #    db.insert("cobra_match", df)


def writeAffineToDB(db, afCoeff, frameId):
    """
    write the affine transformation to DB
    """

    sx = np.sqrt(afCoeff[0, 0]**2+afCoeff[0, 1]**2)
    sy = np.sqrt(afCoeff[1, 0]**2+afCoeff[1, 1]**2)

    xd = afCoeff[0, 2]
    yd = afCoeff[1, 2]

    rotation = np.arctan2(afCoeff[1, 0]/np.sqrt(afCoeff[0, 0]**2+afCoeff[0, 1]**2),
                          afCoeff[1, 1]/np.sqrt(afCoeff[1, 0]**2+afCoeff[1, 1]**2))

    df = pd.DataFrame({'mcs_frame_id': [frameId], 'x_trans': [xd], 'y_trans': [
                      yd], 'x_scale': [sx], 'y_scale': [sy], 'angle': [rotation]})
    db.insert('mcs_pfi_transformation', df)

def writeFidToDB(db, ffid,  mcs_frame_id):

    """
    write the fiducial fibre matches to db.

    only writes fibres with matches
    """

    # get indices of matched FF
    ind=np.where(ffid != -1)
    ffids=ffid[ind]

    # generate the dataframe
    
    pfs_visit_id = mcs_frame_id // 100
    iteration = mcs_frame_id % 100
    sz = len(ffids)
    frame = np.zeros((sz, 6))
    frame[:,0] = np.repeat(pfs_visit_id,sz)
    frame[:,1] = np.repeat(iteration,sz)
    frame[:,2] = np.repeat(mcs_frame_id,sz)
    frame[:,3] = ffids
    frame[:,4] = ind[0]
    frame[:,5] = np.repeat(0,sz)
    columns = ['pfs_visit_id','iteration','mcs_frame_id', 'fiducial_fiber_id', 'spot_id', 'flags']

    df = pd.DataFrame(frame, columns=columns)

    #insert
    db.bulkInsert("fiducial_fiber_match", df)
 
