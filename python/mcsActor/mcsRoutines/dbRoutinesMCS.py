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
import pandas as pd
import numpy as np
import io 
import logging

from opdb import opdb
from datetime import datetime, timezone
from sqlalchemy import text as sqlText


logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S")
logger = logging.getLogger('mcscmd')
logger.setLevel(logging.INFO)

def connectToDB(hostname='', port='', dbname='opdb', username='pfs'):

    """
    create a connection to the DB
    """

    db = opdb.OpDB(hostname=hostname, port=port,
                   dbname=dbname,
                   username=username)
    db.connect()

    return db


def loadTelescopeParametersFromDB(db, frameId):

    """
    read zenithAngle and insrot from mcs_exposure. 
    """

    
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
    read boresight information from database
    table = mcs_boresight
    """
    
    sql = f'''SELECT * FROM mcs_boresight ORDER BY calculated_at DESC FETCH FIRST ROW ONLY'''
    df = db.fetch_query(sql)
    
    return [df['mcs_boresight_x_pix'][0], df['mcs_boresight_y_pix'][0]]


def loadTargetsFromDB(db, frameId):
    """
    load the intermediate target positions from the database
    table=cobra_target

    returns an Nx3 array with the columns equal to the cobra ID, x and y position
    """

    visitId = frameId // 100
    iteration = frameId % 100

    sql = f'SELECT cobra_target.cobra_id,cobra_target.pfi_target_x_mm,cobra_target.pfi_target_y_mm FROM cobra_target WHERE cobra_target.iteration={iteration} and cobra_target.pfs_visit_id={visitId}'
    df = db.bulkSelect('cobra_target',sql)

    return np.array([df['cobra_id'], df['pfi_target_x_mm'], df['pfi_target_y_mm']]).T

def writeTransformToDB(db, frameId, pfiTransform, cameraName):
    """
    write transformation coefficients to database
    """
    pfs_visit_id = frameId // 100
    iteration = frameId % 100

    # just recording mcs_boresight for the first iteration.
    if iteration == 0:
        opDB.insert('mcs_boresight',
                    pfs_visit_id=pfs_visit_id,
                    mcs_boresight_x_pix=pfiTransform.mcs_boresight_x_pix,
                    mcs_boresight_y_pix=pfiTransform.mcs_boresight_y_pix,
                    calculated_at='now')

    trans=pfiTransform.mcsDistort.getArgs()
    res = db.session.execute(sqlText('select * FROM "mcs_pfi_transformation" where false'))
    colnames = tuple(res.keys())
    realcolnames = colnames[0:]
    line = '%d,%f,%f,%f,%e,%e,%f,%s' % (frameId, trans[0].astype('float64'),
           trans[1], trans[3], trans[4],trans[2],
           pfiTransform.alphaRot, cameraName)
                                                                        
    buf = io.StringIO()
    buf.write(line)
    buf.seek(0, 0)
    _writeData(db, 'mcs_pfi_transformation', realcolnames, buf)


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
    
    
def _writeData(db, tableName, columnNames, dataBuf):
    """Wrap a direct COPY_FROM via sqlalchemy. """

    columns = ','.join('"{}"'.format(k) for k in columnNames)
    sql = 'COPY {} ({}) FROM STDIN WITH CSV'.format(
        tableName, columns)
    try:
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

def writeFakeMoveToDB(db, frameId):
    """
    make sure there is a target value if the target doesn't exist
    """
    visitId = frameId // 100
    iteration = frameId % 100

    nCob = 2394
    
    # To-Do here we need a better implementation.
    data = {'pfs_visit_id': np.repeat(visitId,nCob).astype('int'),
            'iteration' : np.repeat(iteration,nCob).astype('int'),
            'cobra_id':np.arange(2394).astype('int')+1,
            #'pfs_config_id':np.repeat(0,2394).astype('int'),
            #'pfi_nominal_x_mm':centers.real,
            #'pfi_nominal_y_mm':centers.imag,
            #'pfi_target_x_mm': centers.real,
            #'pfi_target_y_mm':centers.imag,
            'cobra_mortor_model_id_theta': np.repeat(0,2394).astype('int'),
            'motor_target_theta':np.repeat(0,2394),
            'motor_num_step_theta':np.repeat(0,2394),
            'motor_on_time_theta':np.repeat(0,2394),
            'cobra_mortor_model_id_phi': np.repeat(0,2394).astype('int'),
            'motor_target_phi':np.repeat(0,2394),
            'motor_num_step_phi':np.repeat(0,2394),
            'motor_on_time_phi':np.repeat(0,2394),
            'flags':np.repeat(0,2394).astype('int'),
    }

    df = pd.DataFrame(data=data)
    db.insert("cobra_move", df)



def writeFakeTargetToDB(db, centers, frameId):

    """
    make sure there is a target value if the target doesn't exist
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
            #'cobra_mortor_model_id_theta': np.repeat(0,2394).astype('int'),
            #'motor_target_theta':np.repeat(0,2394),
            #'motor_num_step_theta':np.repeat(0,2394),
            #'motor_on_time_theta':np.repeat(0,2394),
            #'cobra_mortor_model_id_phi': np.repeat(0,2394).astype('int'),
            #'motor_target_phi':np.repeat(0,2394),
            #'motor_num_step_phi':np.repeat(0,2394),
            #'motor_on_time_phi':np.repeat(0,2394),
            'flags':np.repeat(0,2394).astype('int')
    }

    df = pd.DataFrame(data=data)
    db.insert("cobra_target", df)


def writeMatchesToDB(db, cobraMatch, mcsFrameId):
    """
    write the centroids to the database
    table=mcs_data
    variables=spot_id,mcs_center_x_pix,mcs_center_y_pix
              mcs_second_moment_x_pix,mcs_second_moment_y_pix,
              mcs_second_moment_xy_pix,bgvalue,peakvalue
    """

    #np.save("cobraMatch.npy", cobraMatch)
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

def writeFidToDB(db, ffid, mcsData, mcs_frame_id, fids):
    """
    Write the fiducial table to DB, combining with matched MCS data.
    Each fiducial will have its matched spot info if available.
    """
    

    # Prepare output columns
    nFids = len(fids)
    pfs_visit_id = mcs_frame_id // 100
    iteration = mcs_frame_id % 100

    # Prepare arrays for matched data
    spot_id = np.full(nFids, -1, dtype=int)
    pfi_center_x_mm = np.full(nFids, np.nan)
    pfi_center_y_mm = np.full(nFids, np.nan)

    # For each fiducial, check if it was matched in ffid
    for i, fid in enumerate(fids['fiducialId']):
        # Find the index in ffid where this fiducial was matched
        idx = np.where(ffid == fid)[0]
        if len(idx) > 0:
            # Take the first match (or handle multiple matches as needed)
            spot_id[i] = mcsData['spot_id'].iloc[idx[0]]
            pfi_center_x_mm[i] = mcsData['pfi_center_x_mm'].iloc[idx[0]]
            pfi_center_y_mm[i] = mcsData['pfi_center_y_mm'].iloc[idx[0]]

    logging.info(f"Writing {len(fids)} fiducials to DB for frame {mcs_frame_id}")
    
    # Build DataFrame for DB insertion
    df = pd.DataFrame({
        'pfs_visit_id': np.repeat(pfs_visit_id, nFids),
        'iteration': np.repeat(iteration, nFids),
        'mcs_frame_id': np.repeat(mcs_frame_id, nFids),
        'fiducial_fiber_id': fids['fiducialId'],
        'spot_id': spot_id,
        'flags': np.repeat(0, nFids),  # or use your own flags
        'pfi_center_x_mm': pfi_center_x_mm,
        'pfi_center_y_mm': pfi_center_y_mm,
        'match_mask': fids['match_mask'] if 'match_mask' in fids.columns else np.repeat(0, nFids),
        'fiducial_tweaked_x_mm': fids['fiducial_tweaked_x_mm'] if 'fiducial_tweaked_x_mm' in fids.columns else np.repeat(np.nan, nFids),
        'fiducial_tweaked_y_mm': fids['fiducial_tweaked_y_mm'] if 'fiducial_tweaked_y_mm' in fids.columns else np.repeat(np.nan, nFids)
    })
    logging.info(f"fiducial_fiber_match DataFrame shape: {df.shape}")
    logging.info(f"fiducial_tweaked_x_mm: {df['fiducial_tweaked_x_mm'].values}")
    logging.info(f"fiducial_tweaked_y_mm: {df['fiducial_tweaked_y_mm'].values}")
    
    db.insert("fiducial_fiber_match", df)
