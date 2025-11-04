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
import io
import logging
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from opdb import opdb

logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
                    datefmt="%Y-%m-%dT%H:%M:%S")
logger = logging.getLogger('mcscmd')
logger.setLevel(logging.INFO)


def connectToDB(hostname='', port='', dbname='opdb', username='pfs'):
    return opdb.OpDB(hostname=hostname, port=port, dbname=dbname, username=username)  # already connect in init


def loadCobraMatchFromDB(db, frameId):
    sql = (
        "SELECT cobra_id, pfi_center_x_mm, pfi_center_y_mm, flags "
        f"FROM cobra_match WHERE mcs_frame_id={int(frameId)}"
    )
    df = db.fetch_query(sql)
    positions = df['pfi_center_x_mm'] + df['pfi_center_y_mm'] * 1j
    return positions, df['flags']


def loadTelescopeParametersFromDB(db, frameId):
    sql = f'SELECT mcs_exposure.insrot,mcs_exposure.altitude FROM mcs_exposure WHERE mcs_exposure.mcs_frame_id={frameId}'
    df = db.fetch_query(sql)

    if df['altitude'][0] < -99:
        zenithAngle = 90
    else:
        zenithAngle = 90 - df['altitude'][0]

    insRot = df['insrot'][0]

    return zenithAngle, insRot


def loadBoresightFromDB(db, pfsVisitId):
    """Return the most recently calculated boresight (x_pix, y_pix)."""
    sql = "SELECT mcs_boresight_x_pix, mcs_boresight_y_pix " \
          "FROM mcs_boresight ORDER BY calculated_at DESC FETCH FIRST ROW ONLY"
    df = db.fetch_query(sql)
    return [df['mcs_boresight_x_pix'][0], df['mcs_boresight_y_pix'][0]]


def loadCentroidsFromDB(db, mcsFrameId):
    """ retrieve a set of centroids from database and return as a numpy array"""

    sql = 'select mcs_data.spot_id, mcs_data.mcs_center_x_pix, mcs_data.mcs_center_y_pix from mcs_data where ' \
          f'mcs_data.mcs_frame_id={mcsFrameId}'
    df = db.fetch_query(sql)
    return df.to_numpy()


def loadFiducialsFromDB(db):
    """
    load fiducial fibre positions from the DB
    table=fiducial_fiber_geometry

    returns an  Nx3 array with the columns equal to the fiducial ID, x and y position
    """

    sql = f'SELECT fiducial_fiber_geometry.fiducial_fiber_id,fiducial_fiber_geometry.ff_center_on_pfi_x_mm,' \
          'fiducial_fiber_geometry.ff_center_on_pfi_y_mm from fiducial_fiber_geometry'

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

    sql = 'SELECT cobra_target.cobra_id,cobra_target.pfi_target_x_mm,cobra_target.pfi_target_y_mm ' \
          f'FROM cobra_target WHERE cobra_target.iteration={iteration} and cobra_target.pfs_visit_id={visitId}'
    df = db.bulkSelect('cobra_target', sql)

    return np.array([df['cobra_id'], df['pfi_target_x_mm'], df['pfi_target_y_mm']]).T


def writeTransformToDB(db, frameId, pfiTransform, cameraName):
    """Write [x0,y0,theta,dscale,scale2,alpha_rot,camera_name] for one frame into mcs_pfi_transformation.
    If doCloseTransaction=True and a transaction is already active, commit here; else delegate to db.insert().
    """
    pfs_visit_id = frameId // 100
    iteration = frameId % 100

    # just recording mcs_boresight for the first iteration.
    if iteration == 0:
        if iteration == 0:
            db.insert('mcs_boresight', pd.DataFrame({
                'pfs_visit_id': [pfs_visit_id],
                'mcs_boresight_x_pix': [pfiTransform.mcs_boresight_x_pix],
                'mcs_boresight_y_pix': [pfiTransform.mcs_boresight_y_pix],
                'calculated_at': ['now'],
            }))

    mcsDistortCols = ['x0', 'y0', 'theta', 'dscale', 'scale2']

    df = pd.DataFrame(pfiTransform.mcsDistort.getArgs().reshape(1, len(mcsDistortCols)), columns=mcsDistortCols)
    df['mcs_frame_id'] = frameId
    df['alpha_rot'] = float(pfiTransform.alphaRot)
    df['camera_name'] = cameraName

    db.insert('mcs_pfi_transformation', df)

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
    columns = ','.join(f'"{c}"' for c in columnNames)
    sql = f'COPY {tableName} ({columns}) FROM STDIN WITH CSV'

    # make sure to rewind the buffer
    dataBuf.seek(0)

    try:
        with db.engine.begin() as conn:
            with conn.connection.cursor() as cur:  # DBAPI cursor
                cur.copy_expert(sql, dataBuf)
        return True

    except Exception as e:
        logger.warning(f"failed to write with {sql}: {e}")
        return False


def _readData(db, sql):
    dataBuf = io.StringIO()
    with db.engine.begin() as conn:
        with conn.connection.cursor() as cur:
            cur.copy_expert(sql, dataBuf)  # sql is a plain string
    dataBuf.seek(0)
    return dataBuf


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
            'cobra_mortor_model_id_theta': np.repeat(0, 2394).astype('int'),
            'motor_target_theta': np.repeat(0, 2394),
            'motor_num_step_theta': np.repeat(0, 2394),
            'motor_on_time_theta': np.repeat(0, 2394),
            'cobra_mortor_model_id_phi': np.repeat(0, 2394).astype('int'),
            'motor_target_phi': np.repeat(0, 2394),
            'motor_num_step_phi': np.repeat(0, 2394),
            'motor_on_time_phi': np.repeat(0, 2394),
            'flags': np.repeat(0, 2394).astype('int'),
            }

    df = pd.DataFrame(data=data)
    db.insert("cobra_move", df)


def writeFakeTargetToDB(db, centers, frameId):
    """
    Bulk-COPY a minimal target set into cobra_target using cobra centers.
    """
    visitId = int(frameId // 100)
    iteration = int(frameId % 100)
    nCob = 2394

    columns = [
        'pfs_visit_id',
        'iteration',
        'cobra_id',
        'pfs_config_id',
        'pfi_nominal_x_mm',
        'pfi_nominal_y_mm',
        'pfi_target_x_mm',
        'pfi_target_y_mm',
        'flags',
    ]

    df = pd.DataFrame({
        'pfs_visit_id': np.repeat(visitId, nCob).astype(int),
        'iteration': np.repeat(iteration, nCob).astype(int),
        'cobra_id': (np.arange(nCob) + 1).astype(int),
        'pfs_config_id': np.repeat(0, nCob).astype(int),
        'pfi_nominal_x_mm': centers.real,
        'pfi_nominal_y_mm': centers.imag,
        'pfi_target_x_mm': centers.real,
        'pfi_target_y_mm': centers.imag,
        'flags': np.repeat(0, nCob).astype(int),
    })

    buf = io.StringIO()
    df.to_csv(buf, index=False, header=False)
    buf.seek(0)
    _writeData(db, 'cobra_target', columns, buf)


def writeCobraCenterToDB(db, frameId, centers, mpos):
    visitId = frameId // 100
    iteration = frameId % 100

    targetTable = {'pfs_visit_id': np.repeat(visitId, 2394),
                   'iteration': np.repeat(iteration, 2394),
                   'cobra_id': np.arange(2394) + 1,
                   'pfs_config_id': np.repeat(0, 2394),
                   'pfi_nominal_x_mm': centers.real,
                   'pfi_nominal_y_mm': centers.imag,
                   'pfi_target_x_mm': mpos.real,
                   'pfi_target_y_mm': mpos.imag,
                   'motor_target_theta': np.repeat(0, 2394),
                   'motor_target_phi': np.repeat(0, 2394),
                   }

    df = pd.DataFrame(data=targetTable)
    db.bulkInsert("cobra_target", df)


def writeBoresightToDB(db, pfsVisitId, boresight):
    """ write boresight to database with current timestamp """

    dt = datetime.now(timezone.utc)

    df = pd.DataFrame(
        {'pfs_visit_id': [pfsVisitId], 'mcs_boresight_x_pix': [boresight[0]], 'mcs_boresight_y_pix': [boresight[1]],
         'calculated_at': [dt]})
    db.bulkInsert('mcs_boresight', df)


def writeCentroidsToDB(db, centroids, mcsFrameId):
    """
    Bulk-COPY insert into mcs_data (atomic, no Session state).
    Expects `centroids` as Nx8 or Nx? array with columns:
      [spot_id, x, y, x2, y2, xy, peak, bg] (your current order).
    """
    n = centroids.shape[0]
    if n == 0:
        return

    # Map your array -> table columns explicitly
    columns = [
        'mcs_frame_id',
        'spot_id',
        'mcs_center_x_pix',
        'mcs_center_y_pix',
        'mcs_second_moment_x_pix',
        'mcs_second_moment_y_pix',
        'mcs_second_moment_xy_pix',
        'peakvalue',
        'bgvalue',
    ]

    df = pd.DataFrame({
        'mcs_frame_id': np.repeat(int(mcsFrameId), n),
        'spot_id': centroids[:, 0].astype(int),
        'mcs_center_x_pix': centroids[:, 1],
        'mcs_center_y_pix': centroids[:, 2],
        'mcs_second_moment_x_pix': centroids[:, 3],
        'mcs_second_moment_y_pix': centroids[:, 4],
        'mcs_second_moment_xy_pix': centroids[:, 5],
        'peakvalue': centroids[:, 6],
        'bgvalue': centroids[:, 7],
    })

    buf = io.StringIO()
    df.to_csv(buf, index=False, header=False)
    buf.seek(0)
    _writeData(db, 'mcs_data', columns, buf)


def readMatchFromDB(db, mcsFrameId):
    sql = f"select * from cobra_match where mcs_frame_id = {int(mcsFrameId)}"
    return db.bulkSelect('cobra_match', sql)


def writeMatchesToDB(db, cobraMatch, mcsFrameId):
    """
    Bulk-COPY insert into cobra_match.
    `cobraMatch` columns: [cobra_id, spot_id, pfi_x_mm, pfi_y_mm, flags]
    """
    n = cobraMatch.shape[0]
    if n == 0:
        return

    visitId = int(mcsFrameId // 100)
    iteration = int(mcsFrameId % 100)

    columns = [
        'pfs_visit_id',
        'iteration',
        'mcs_frame_id',
        'cobra_id',
        'spot_id',
        'pfi_center_x_mm',
        'pfi_center_y_mm',
        'flags',
    ]

    df = pd.DataFrame({
        'pfs_visit_id': np.repeat(visitId, n),
        'iteration': np.repeat(iteration, n),
        'mcs_frame_id': np.repeat(int(mcsFrameId), n),
        'cobra_id': cobraMatch[:, 0].astype(int),
        'spot_id': cobraMatch[:, 1].astype(int),
        'pfi_center_x_mm': cobraMatch[:, 2],
        'pfi_center_y_mm': cobraMatch[:, 3],
        'flags': cobraMatch[:, 4].astype(int),
    })

    buf = io.StringIO()
    df.to_csv(buf, index=False, header=False)
    buf.seek(0)
    _writeData(db, 'cobra_match', columns, buf)


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

    df = pd.DataFrame({'mcs_frame_id': [frameId], 'x_trans': [xd], 'y_trans': [
        yd], 'x_scale': [sx], 'y_scale': [sy], 'angle': [rotation]})
    db.insert('mcs_pfi_transformation', df)


# --- new helper: write mcs_exposure row via COPY ---
def writeMcsExposure(db, telescope_info: dict, camera_name: str, centroid_method: str):
    """
    COPY one row into mcs_exposure.

    Parameters
    ----------
    db : opdb.OpDB
    telescope_info : dict
        Keys: frameid, visitid, exptime(ms), altitude, azimuth, instrot, adc_pa,
              dome_temperature, dome_pressure, dome_humidity,
              outside_temperature, outside_pressure, outside_humidity,
              mcs_cover_temperature, mcs_m1_temperature, starttime (ISO string)
    camera_name : str
    centroid_method : str
        'win' -> 'windowed', anything else -> 'sep'
    """
    columns = [
        'mcs_frame_id',
        'pfs_visit_id',
        'exptime',
        'altitude',
        'azimuth',
        'insrot',
        'adc_pa',
        'dome_temperature',
        'dome_pressure',
        'dome_humidity',
        'outside_temperature',
        'outside_pressure',
        'outside_humidity',
        'mcs_cover_temperature',
        'mcs_m1_temperature',
        'taken_at',
        'taken_in_hst_at',
        'mcs_camera_id',
        'measurement_algorithm',
        'version_actor',
        'version_instdata',
    ]

    if camera_name == 'canon_50m':
        mcs_camera_id = 0
    elif camera_name == 'rmod_71m':
        mcs_camera_id = 1
    else:
        mcs_camera_id = -1

    measurement_algorithm = 'windowed' if centroid_method == 'win' else 'sep'

    df = pd.DataFrame([{
        'mcs_frame_id': int(telescope_info['frameid']),
        'pfs_visit_id': int(telescope_info['visitid']),
        'exptime': float(telescope_info['exptime']) / 1000.0,
        'altitude': float(telescope_info['altitude']),
        'azimuth': float(telescope_info['azimuth']),
        'insrot': float(telescope_info['instrot']),
        'adc_pa': float(telescope_info['adc_pa']),
        'dome_temperature': float(telescope_info['dome_temperature']),
        'dome_pressure': float(telescope_info['dome_pressure']),
        'dome_humidity': float(telescope_info['dome_humidity']),
        'outside_temperature': float(telescope_info['outside_temperature']),
        'outside_pressure': float(telescope_info['outside_pressure']),
        'outside_humidity': float(telescope_info['outside_humidity']),
        'mcs_cover_temperature': float(telescope_info['mcs_cover_temperature']),
        'mcs_m1_temperature': float(telescope_info['mcs_m1_temperature']),
        'taken_at': telescope_info['starttime'],
        'taken_in_hst_at': telescope_info['starttime'],
        'mcs_camera_id': int(mcs_camera_id),
        'measurement_algorithm': measurement_algorithm,
        'version_actor': '0',
        'version_instdata': '0',
    }])

    buf = io.StringIO()
    df.to_csv(buf, index=False, header=False)
    buf.seek(0)
    _writeData(db, 'mcs_exposure', columns, buf)


# --- new helper: write cobra_target rows via COPY ---
def writeCobraTarget(db, df_targets: pd.DataFrame):
    """
    COPY many rows into cobra_target. Expects exactly these columns in df:
      ['pfs_visit_id','iteration','cobra_id','pfs_config_id',
       'pfi_nominal_x_mm','pfi_nominal_y_mm','pfi_target_x_mm',
       'pfi_target_y_mm','flags']
    """
    required = [
        'pfs_visit_id', 'iteration', 'cobra_id', 'pfs_config_id',
        'pfi_nominal_x_mm', 'pfi_nominal_y_mm',
        'pfi_target_x_mm', 'pfi_target_y_mm', 'flags',
    ]
    missing = [c for c in required if c not in df_targets.columns]
    if missing:
        raise ValueError(f'write_cobra_target: missing columns: {missing}')

    buf = io.StringIO()
    df_targets[required].to_csv(buf, index=False, header=False)
    buf.seek(0)
    _writeData(db, 'cobra_target', required, buf)


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

    db.insert("fiducial_fiber_match", df)
