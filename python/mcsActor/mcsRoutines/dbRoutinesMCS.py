
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
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from pfs.utils.database import opdb

logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
                    datefmt="%Y-%m-%dT%H:%M:%S")
logger = logging.getLogger('mcscmd')
logger.setLevel(logging.INFO)


def connectToDB(hostname=None, port=None, dbname=None, username=None):

    """ connect to database """
    
    raise NotImplementedError()

    db = opdb.OpDB(hostname=hostname, port=port,
                   dbname=dbname,
                   username=username)
    return db


def loadTelescopeParametersFromDB(db, frameId):

    """ load telescope parameters from database """
    
    sql = f'SELECT mcs_exposure.insrot,mcs_exposure.altitude FROM mcs_exposure WHERE mcs_exposure.mcs_frame_id={frameId}'
    df = db.query_dataframe(sql)

    if df['altitude'][0] < -99:
        zenithAngle = 90
    else:
        zenithAngle = 90 - df['altitude'][0]

    insRot = df['insrot'][0]

    return zenithAngle, insRot


def loadTargetsFromDB(db, frameId):
    """
    load the intermediate target positions from the database
    table=cobra_target

    returns an Nx3 array with the columns equal to the cobra ID, x and y position
    """

    visitId = frameId // 100
    iteration = frameId % 100

    sql = f'SELECT cobra_target.cobra_id,cobra_target.pfi_target_x_mm,cobra_target.pfi_target_y_mm FROM cobra_target WHERE cobra_target.iteration={iteration} and cobra_target.pfs_visit_id={visitId}'
    df = db.query_dataframe(sql)

    return np.array([df['cobra_id'], df['pfi_target_x_mm'], df['pfi_target_y_mm']]).T


def writeTransformToDB(db, frameId, pfiTransform, cameraName, doCloseTransaction=False):
    """Write [x0,y0,theta,dscale,scale2,alpha_rot,camera_name] for one frame into mcs_pfi_transformation.
    If doCloseTransaction=True and a transaction is already active, commit here; else delegate to db.insert().
    """
    pfs_visit_id = frameId // 100
    iteration = frameId % 100

    # just recording mcs_boresight for the first iteration.
    if iteration == 0:
        db.insert_kw('mcs_boresight',
                     pfs_visit_id=pfs_visit_id,
                     mcs_boresight_x_pix=float(pfiTransform.mcs_boresight_x_pix),
                     mcs_boresight_y_pix=float(pfiTransform.mcs_boresight_y_pix),
                     calculated_at='now')

    mcsDistortCols = ['x0', 'y0', 'theta', 'dscale', 'scale2']

    df = pd.DataFrame(pfiTransform.mcsDistort.getArgs().reshape(1, len(mcsDistortCols)), columns=mcsDistortCols)
    df['mcs_frame_id'] = frameId
    df['alpha_rot'] = float(pfiTransform.alphaRot)
    df['camera_name'] = cameraName

    db.insert_dataframe('mcs_pfi_transformation', df)

    return (df['mcs_frame_id'].values, df['x0'].values, df['y0'].values, df['dscale'].values,
            df['scale2'].values, df['theta'].values, df['alpha_rot'].values, df['camera_name'].values)


def writeFakeMoveToDB(db, frameId):
    """
    make sure there is a target value if the target database is empty
    """
    visitId = frameId // 100
    iteration = frameId % 100

    nCob = 2394

    # To-Do here we need a better implementation.
    data = {'pfs_visit_id': np.repeat(visitId, nCob).astype('int'),
            'iteration': np.repeat(iteration, nCob).astype('int'),
            'cobra_id': np.arange(2394).astype('int') + 1,
            # 'pfs_config_id':np.repeat(0,2394).astype('int'),
            # 'pfi_nominal_x_mm':centers.real,
            # 'pfi_nominal_y_mm':centers.imag,
            # 'pfi_target_x_mm': centers.real,
            # 'pfi_target_y_mm':centers.imag,
            'cobra_motor_model_id_theta': np.repeat(0, 2394).astype('int'),
            'motor_target_theta': np.repeat(0, 2394),
            'motor_num_step_theta': np.repeat(0, 2394),
            'motor_on_time_theta': np.repeat(0, 2394),
            'cobra_motor_model_id_phi': np.repeat(0, 2394).astype('int'),
            'motor_target_phi': np.repeat(0, 2394),
            'motor_num_step_phi': np.repeat(0, 2394),
            'motor_on_time_phi': np.repeat(0, 2394),
            'flags': np.repeat(0, 2394).astype('int'),
            }

    df = pd.DataFrame(data=data)
    db.insert_dataframe("cobra_move", df)


def writeFakeTargetToDB(db, centers, frameId):
    """
    make sure there is a target value if the target database is empty
    """
    visitId = frameId // 100
    iteration = frameId % 100

    nCob = 2394

    # To-Do here we need a better implementation.
    data = {'pfs_visit_id': np.repeat(visitId, nCob).astype('int'),
            'iteration': np.repeat(iteration, nCob).astype('int'),
            'cobra_id': np.arange(nCob).astype('int') + 1,
            'pfs_config_id': np.repeat(0, nCob).astype('int'),
            'pfi_nominal_x_mm': centers.real,
            'pfi_nominal_y_mm': centers.imag,
            'pfi_target_x_mm': centers.real,
            'pfi_target_y_mm': centers.imag,
            'flags': np.repeat(0, nCob).astype('int')
            }

    df = pd.DataFrame(data=data)
    db.insert_dataframe("cobra_target", df)

def writeMatchesToDB(db, cobraMatch, mcsFrameId):
    """
    write the centroids to the database
    table=mcs_data
    variables=spot_id,mcs_center_x_pix,mcs_center_y_pix
              mcs_second_moment_x_pix,mcs_second_moment_y_pix,
              mcs_second_moment_xy_pix,bgvalue,peakvalue
    """

    visitId = mcsFrameId // 100
    iteration = mcsFrameId % 100
    targetTable = {'pfs_visit_id': np.repeat(visitId, 2394).astype('int'),
                   'iteration': np.repeat(iteration, 2394).astype('int'),
                   'cobra_id': cobraMatch[:, 0].astype('int'),
                   'mcs_frame_id': np.repeat(mcsFrameId, 2394).astype('int'),
                   'spot_id': cobraMatch[:, 1].astype('int'),
                   'pfi_center_x_mm': cobraMatch[:, 2],
                   'pfi_center_y_mm': cobraMatch[:, 3],
                   'flags': cobraMatch[:, 4].astype('int')
                   }
    df = pd.DataFrame(data=targetTable)
    db.insert_dataframe("cobra_match", df)

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
        'fiducial_tweaked_x_mm': fids[
            'fiducial_tweaked_x_mm'] if 'fiducial_tweaked_x_mm' in fids.columns else np.repeat(np.nan, nFids),
        'fiducial_tweaked_y_mm': fids[
            'fiducial_tweaked_y_mm'] if 'fiducial_tweaked_y_mm' in fids.columns else np.repeat(np.nan, nFids)
    })
    logging.info(f"fiducial_fiber_match DataFrame shape: {df.shape}")
    logging.info(f"fiducial_tweaked_x_mm: {df['fiducial_tweaked_x_mm'].values}")
    logging.info(f"fiducial_tweaked_y_mm: {df['fiducial_tweaked_y_mm'].values}")

    db.insert_dataframe("fiducial_fiber_match", df)
