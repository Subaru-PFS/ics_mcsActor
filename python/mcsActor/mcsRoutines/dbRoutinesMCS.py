
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
    raise NotImplementedError()

    db = opdb.OpDB(hostname=hostname, port=port,
                   dbname=dbname,
                   username=username)
    return db

def loadCobraMatchFromDB(db, frameId):
    """
    read the cobra_match table information, convert to x+ij format, return positions and flags
    """

    sql = f'SELECT cobra_match.cobra_id,cobra_match.pfi_center_x_mm,cobra_match.pfi_center_x_mm,cobra_match.flags from cobra_match where cobra_match.mcs_frame_id={frameId}'

    df = db.query_dataframe(sql)

    positions = df['pfi_center_x_mm'] + df['pfi_center_y_mm'] * 1j
    flags = df['flags']

    return positions, flags


def loadTelescopeParametersFromDB(db, frameId):

    sql = f'SELECT mcs_exposure.insrot,mcs_exposure.altitude FROM mcs_exposure WHERE mcs_exposure.mcs_frame_id={frameId}'
    df = db.query_dataframe(sql)

    if df['altitude'][0] < -99:
        zenithAngle = 90
    else:
        zenithAngle = 90 - df['altitude'][0]

    insRot = df['insrot'][0]

    return zenithAngle, insRot


def loadBoresightFromDB(db, pfsVisitId):
    """
    read boresight informatino from database
    table = mcs_boresight
    """
    sql = f'''SELECT * FROM mcs_boresight ORDER BY calculated_at DESC FETCH FIRST ROW ONLY'''
    # sql=f'SELECT mcs_boresight.mcs_boresight_x_pix,mcs_boresight.mcs_boresight_y_pix from mcs_boresight where mcs_boresight.pfs_visit_id={pfsVisitId}'
    df = db.query_dataframe(sql)
    return [df['mcs_boresight_x_pix'][0], df['mcs_boresight_y_pix'][0]]


def loadCentroidsFromDB(db, mcsFrameId):
    """ retrieve a set of centroids from database and return as a numpy array"""

    sql = f'select mcs_data.spot_id, mcs_data.mcs_center_x_pix, mcs_data.mcs_center_y_pix from mcs_data where mcs_data.mcs_frame_id={mcsFrameId}'
    df = db.query_dataframe(sql)
    return df.to_numpy()


def loadFiducialsFromDB(db):
    """
    load fiducial fibre positions from the DB
    table=fiducial_fiber_geometry

    returns an  Nx3 array with the columns equal to the fiducial ID, x and y position
    """

    sql = f'SELECT fiducial_fiber_geometry.fiducial_fiber_id , fiducial_fiber_geometry.ff_center_on_pfi_x_mm ,fiducial_fiber_geometry.ff_center_on_pfi_y_mm from fiducial_fiber_geometry'

    df = db.query_dataframe(sql)

    return np.array([df['fiducial_fiber_id'], df['ff_center_on_pfi_x_mm'], df['ff_center_on_pfi_y_mm']]).T


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


def _writeData(db, tableName, columnNames, dataBuf):
    """
    COPY FROM STDIN (CSV) into a table using the SQLAlchemy Session.

    Parameters
    ----------
    db : opdb.opdb.OpDB
        Handle exposing `db.session` (SQLAlchemy Session).
    tableName : str
        Target table name, optionally schema-qualified (e.g. 'schema.table').
    columnNames : iterable[str]
        Columns to load, in the exact order they appear in the CSV stream.
    dataBuf : io.StringIO
        Text file-like object positioned anywhere; will be rewound.

    Behavior
    --------
    - Reuses an active transaction on the Session; otherwise opens one and commits.
    - Rolls back on error and re-raises the exception.
    - Does not close the Session; only the DBAPI cursor is context-managed.

    Raises
    ------
    Exception
        Any error from COPY or the database driver; the Session is rolled back.
    """

    raise NotImplementedError()

    columns = ','.join(f'"{c}"' for c in columnNames)
    sql = f'COPY {tableName} ({columns}) FROM STDIN WITH CSV'
    session = db.session
    # make sure to rewind the buffer
    dataBuf.seek(0)

    def doCopy():
        conn = session.connection()  # SQLA Connection (bound to current txn)
        raw = conn.connection  # DBAPI connection (psycopg2)
        with raw.cursor() as cur:
            cur.copy_expert(sql, dataBuf)  # context manager closes cursor

    try:
        if session.in_transaction():
            doCopy()  # reuse caller's transaction
        else:
            with session.begin():
                doCopy()  # commit on success

        return True

    except Exception as e:
        if session.in_transaction():
            session.rollback()

        logger.warning(f"failed to write with {sql}: {e}")

        return False


def writeTargetToDB(db, frameId, target, mpos):
    visitId = frameId // 100
    iteration = frameId % 100

    # To-Do here we need a better implementation.
    data = {'pfs_visit_id': np.repeat(visitId, 2394),
            'iteration': np.repeat(iteration, 2394),
            'cobra_id': np.arange(2394) + 1,

            }

    df = pd.DataFrame(data=data)
    db.insert_dataframe("cobra_target", df)


def writeFakeMoveToDB(db, frameId):
    """
    make sure there is a target value if the target d
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
    make sure there is a target value if the target d
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

def writeBoresightToDB(db, pfsVisitId, boresight):
    """ write boresight to database with current timestamp """

    dt = datetime.now(timezone.utc)

    db.insert_kw('mcs_boresight',
                 pfs_visit_id=pfsVisitId,
                 mcs_boresight_x_pix=boresight[0],
                 mcs_boresight_y_pix=boresight[1],
                 calculated_at=dt)

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
    db.insert_dataframe("mcs_data", df)


def readMatchFromDB(db, mcsFrameId):
    match = db.query_dataframe(f'select * from cobra_match where '
                               f'mcs_frame_id = {mcsFrameId}')
    return match


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

def writeAffineToDB(db, afCoeff, frameId):
    """
    write the affine transformation to DB
    """

    sx = np.sqrt(afCoeff[0, 0] ** 2 + afCoeff[0, 1] ** 2)
    sy = np.sqrt(afCoeff[1, 0] ** 2 + afCoeff[1, 1] ** 2)

    xd = afCoeff[0, 2]
    yd = afCoeff[1, 2]

    rotation = np.arctan2(afCoeff[1, 0] / np.sqrt(afCoeff[0, 0] ** 2 + afCoeff[0, 1] ** 2),
                          afCoeff[1, 1] / np.sqrt(afCoeff[1, 0] ** 2 + afCoeff[1, 1] ** 2))

    df = pd.DataFrame({'mcs_frame_id': [frameId],
                       'x_trans': [xd],
                       'y_trans': [yd],
                       'x_scale': [sx],
                       'y_scale': [sy],
                       'angle': [rotation]})
    db.insert_dataframe('mcs_pfi_transformation', df)


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
