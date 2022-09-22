#!/usr/bin/env python

import cv2
import logging
import numpy as np
from pfs.utils.coordinates import CoordTransp
from pfs.utils import coordinates
import pandas as pd
import time
import datetime
import io
import pathlib
import queue
import threading
import sep
import sys

import os
import astropy.io.fits as pyfits
import astropy
import fitsio

import opscore.protocols.keys as keys
import opscore.protocols.types as types

#from actorcore.utility import fits as fitsUtils
#from actorcore.utility import timecards

from ics.utils.fits import mhs as fitsUtils
from ics.utils.fits import timecards

from opscore.utility.qstr import qstr

import pfs.utils.coordinates.transform as transformUtils
import pfs.utils.coordinates.MakeWCSCards as pfsWcs
from scipy.spatial import cKDTree

import psycopg2
import psycopg2.extras
from xml.etree.ElementTree import dump
from procedures.moduleTest import calculation

import mcsActor.windowedCentroid.centroid as centroid
import mcsActor.mcsRoutines.mcsRoutines as mcsTools
import mcsActor.mcsRoutines.dbRoutinesMCS as dbTools
from pfs.utils import butler

from opdb import opdb

from importlib import reload
reload(dbTools)
reload(mcsTools)

reload(CoordTransp)

class McsCmd(object):

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor
        self.expTime = 1000
        self.newTable = None
        self.simulationPath = None
        self._connectionToDB = None
        self._conn = None
        self._db = None

        self.findThresh = None
        self.centThresh = None
        self.cMethod = 'win'
        self.fMethod = 'target'

        # self.setCentroidParams(None)
        self.adjacentCobras = None
        self.geometrySet = False
        self.geomFile = None
        self.dotFile = None
        self.fidsGood = None
        self.fidsOuterRing = None

        logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S")
        self.logger = logging.getLogger('mcscmd')
        self.logger.setLevel(logging.INFO)
        self.butler = butler.Butler()

        # Declare the commands we implement. When the actor is started
        # these are registered with the parser, which will call the
        # associated methods when matched. The callbacks will be
        # passed a single argument, the parsed and typed command.
        #
        self.vocab = [
            ('ping', '', self.ping),
            ('status', '', self.status),
            ('expose', '@(bias|test) [<frameId>]', self.expose),
            ('expose', '@(dark|flat) <expTime> [<frameId>]', self.expose),
            ('expose',
                'object <expTime> [<frameId>] [@noCentroid] [@doCentroid] [@doFibreID] [@simDot] '
                '[@newField]', self.expose),
            ('runCentroid', '[@newTable]', self.runCentroid),
            #('runFibreID', '[@newTable]', self.runFibreID),
            ('reconnect', '', self.reconnect),
            ('resetThreshold', '', self.resetThreshold),
            ('setCentroidParams', '[<fwhmx>] [<fwhmy>] [<boxFind>] [<boxCent>] [<nmin>] [<maxIt>] [<centSigma>] [<findSigma>]',
             self.setCentroidParams),
            ('calcThresh', '[<threshSigma>] [<threshFact>]', self.calcThresh),
            ('simulate', '<path>', self.simulateOn),
            ('simulate', 'off', self.simulateOff),
            ('switchCMethod', '<cMethod>', self.switchCMethod),
            ('switchFMethod', '<fMethod>', self.switchFMethod),
            ('resetGeometry', '', self.resetGeometry),
            ('resetGeometryFile', '<geomFile>', self.resetGeometryFile)
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time, seconds"),
                                        keys.Key("expType", types.String(), help="The exposure type"),
                                        keys.Key("cameraName", types.String(), help="The camera used for MCS"),
                                        keys.Key("frameId", types.Int(), help="exposure frameID"),
                                        keys.Key("filename", types.String(), help="exposure filename"),
                                        keys.Key("path", types.String(), help="Simulated image directory"),
                                        keys.Key("getArc", types.Int(), help="flag for arc image"),
                                        keys.Key("fwhmx", types.Float(), help="X fwhm for centroid routine"),
                                        keys.Key("fwhmy", types.Float(), help="Y fwhm for centroid routine"),
                                        keys.Key("boxFind", types.Int(), help="box size for finding spots"),
                                        keys.Key("boxCent", types.Int(),
                                                 help="box size for centroiding spots"),
                                        keys.Key("nmin", types.Int(),
                                                 help="minimum number of points for spot"),
                                        keys.Key("maxIt", types.Int(),
                                                 help="maximum number of iterations for centroiding"),
                                        keys.Key("findSigma", types.Float(),
                                                 help="threshhold for finding spots"),
                                        keys.Key("centSigma", types.Float(),
                                                 help="threshhold for calculating moments of spots"),
                                        keys.Key("threshSigma", types.Float(),
                                                 help="threshhold calculating background level"),
                                        keys.Key("threshFact", types.Float(),
                                                 help="factor for engineering threshold measurements"),
                                        keys.Key("matchRad", types.Int(),
                                                 help="radius in pixels for matching positions"),
                                        keys.Key("geomFile", types.String(), help="file for geometry"),
                                        keys.Key("dotFile", types.String(), help="file for dot information"),
                                        keys.Key("fieldID", types.String(),
                                                 help="fieldID for getting instrument parameters"),
                                        keys.Key("cMethod", types.String(),
                                                 help="method for centroiding (of 'win', 'cent', default 'win')"),
                                        keys.Key("fMethod", types.String(),
                                                 help="method for fibreId (of 'target', 'previous', default 'target')")
                                        )

    def connectToDB(self, cmd=None):
        """connect to the database if not already connected"""

        if self._db is not None:
            return self._db

        if cmd is None:
            cmd = self.actor.bcast

        try:
            config = self.actor.config
            hostname = config.get('db', 'hostname', fallback='db-ics')
            dbname = config.get('db', 'dbname', fallback='opdb')
            port = config.get('db', 'port', fallback=5432)
            username = config.get('db', 'username', fallback='pfs')
        except Exception as e:
            raise RuntimeError(f'failed to load opdb configuration: {e}')

        cmd.diag(f'text="connecting to db, {username}@{hostname}:{port}, db={dbname}"')
        try:
            db = dbTools.connectToDB(hostname=hostname,
                                     port=port,
                                     dbname=dbname,
                                     username=username)
        except Exception as e:
            raise RuntimeError(f"unable to connect to the database: {e}")

        if cmd is not None:
            cmd.inform('text="Connected to Database"')

        self._db = db
        return self._db
    

    def ping(self, cmd):
        """Query the actor for liveness/happiness."""
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd, self.expTime)

        cmd.finish("text='Present and (probably) well'")

    def reconnect(self, cmd):
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd, self.expTime)

        cmdKeys = cmd.cmd.keywords
        cmdKeys['cameraName'] = self.actor.cameraName

        cmd.finish('text="Camera connected!"')

    def status(self, cmd):
        """Report status and version; obtain and send current data"""

        self.actor.sendVersionKey(cmd)
        self.actor.camera.sendStatusKeys(cmd)
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd, self.expTime)

        cmd.inform(f'text="MCS camera present! camera name = {self.actor.cameraName}"')
        cmd.finish()

    def simulateOff(self, cmd):
        """turn off simulation mode"""

        self.simulationPath = None
        cmd.finish('text="set simulation path to %s"' % str(self.simulationPath))

    def simulateOn(self, cmd):
        """set path for simulated images"""

        cmdKeys = cmd.cmd.keywords

        path = cmdKeys['path'].values[0]

        if not os.path.isdir(path):
            cmd.fail('text="path %s is not a directory"' % (path))
            return

        self.simulationPath = (path, 0, '')
        cmd.finish('text="set simulation path to %s"' % str(self.simulationPath))

    def getNextSimulationImage(self, cmd):
        import glob

        path, idx, lastname = self.simulationPath
        cmd.inform('text="frameIds1 = %s." '%{path})
        files = sorted(glob.glob(os.path.join(path, 'PFSC*.fits')))

        cmd.debug('text="%i of %i files in %s..."' % (idx, len(files), path))
        if len(files) == 0:
            raise RuntimeError(f"no .fits files in {path}")

        if idx+1 > len(files):
            # I don't think this is what we want: when the data set
            # has been read we should *stop*, not loop back.
            # idx = 0
            raise RuntimeError(f"no more .fits files in {path}")

        imagePath = files[idx]
        cmd.inform('text="frameIds2 = %s." '%{imagePath})
        image, hdr = pyfits.getdata(imagePath, 0, header=True)

        self.simulationPath = (path, idx+1, imagePath)
        imagePath = pathlib.Path(imagePath)
        frameId = int(imagePath.stem[4:], base=10)
        self.visitId = frameId // 100

        try:
            fwfile = sorted(glob.glob(os.path.join(path, 'thetaFW.npy')))
            rvfile = sorted(glob.glob(os.path.join(path, 'thetaRV.npy')))
            fw = np.load(fwfile[0])
            rv = np.load(rvfile[0])
            pos_array = np.append(fw, rv, axis=2)
            targets = pos_array[:, 0, idx]
        except:
            targets = np.zeros(2394)+np.zeros(2394)*1j

        cmd.inform('text="returning simulation file %s"' % (imagePath))
        return imagePath, image, targets

    def requestNextFileIds(self, cmd, frameId):
        """ Return a queue which will eventually contain our fileIda. """

        q = queue.Queue()
        cmd.inform('text="frameIds = %s." '%{frameId})

        def worker(q=q, cmd=cmd):
            filename = self.getNextFileIds(cmd, frameId)
            q.put(filename)
            q.task_done()

        task = threading.Thread(target=worker,
                                name='visitFetcher', daemon=True)
        task.start()

        return q

    def getNextFileIds(self, cmd, frameId):
        """ Fetch next image filename components

        Args
        ----
        cmd : MHS Command
          who to report to
        frameId : `int`
          The full, 8-digit frame number. Or None if we should fetch a new one.

        Returns
        -------
        fileIds : dictionary of ids sufficient for butler.getPath('mcsFile')

        """

        if frameId is None:
            ret = self.actor.cmdr.call(actor='gen2', cmdStr='getVisit caller=mcs', timeLim=15.0)
            if ret.didFail:
                raise RuntimeError("getNextFilename failed getting a visit number in 15s!")
            visit = self.actor.models['gen2'].keyVarDict['visit'].valueList[0]
            subframeId = 0
            cmd.inform(f'text="gen2.getVisit = {visit}"')
        else:
            visit = frameId // 100
            subframeId = frameId % 100

            t0 = time.time()
            ret = self.actor.cmdr.call(actor='gen2', cmdStr='updateTelStatus caller=mcs', timeLim=5.0)
            if ret.didFail:
                raise RuntimeError("getNextFilename failed updating telesceop status in 15s!")
            t1 = time.time()
            if t1-t0 >= 2:
                cmd.warn(f'text="it took {t1-t0:0.2f}s to update telescope status"')

        idDict = dict(pfsDay=self.butler.getPfsDay(),
                      visit=visit,
                      frame=subframeId,
                      frameId=frameId)

        return idDict

    def _getInstHeader(self, cmd):
        """ Gather FITS cards from all actors we are interested in. """

        modelNames = set(self.actor.models.keys())
        cmd.debug('text="fetching MHS cards for %s..."' % (modelNames))
        cards = fitsUtils.gatherHeaderCards(cmd, self.actor,
                                            modelNames=modelNames, shortNames=True)
        cmd.debug('text="fetched %d MHS cards..."' % (len(cards)))

        # Until we convert to fitsio, convert cards to pyfits
        pycards = []
        for c in cards:
            if isinstance(c, str):
                pcard = 'COMMENT', c
            else:
                pcard = c['name'], c['value'], c.get('comment', '')
            pycards.append(pcard)
            cmd.debug('text=%s' % (qstr("fetched card: %s" % (str(pcard)))))

        if self.simulationPath is not None:
            pcard = 'W_MCSMNM', self.simulationPath[-1], 'Simulated file path'
            pycards.append(pcard)

        return pycards

    def _constructHeader(self, cmd, filename, expType, expTime, expStart):
        if expType == 'bias':
            expTime = 0.0

        # Subaru rules
        if expType == 'object':
            expType = 'acquisition'

        hdr = pyfits.Header()

        try:
            detectorTemp = self.actor.models['meb'].keyVarDict['temps'].valueList[1]
        except Exception as e:
            cmd.warn('text="FAILED to fetch MEB cards: %s"' % (e))
            detectorTemp = str(np.nan)

        try:
            config = pfsConfig.getConfigDict('mcs')
            detectorId = config['serial']
            gain = config['gain']
        except Exception as e:
            cmd.warn('text="FAILED to fetch config cards: %s"' % (e))
            detectorId = 'MCS cam#1'
            gain = 2.24

        hdr.append(('DETECTOR', detectorId, 'Name of the detector/CCD'))
        hdr.append(('GAIN', gain, '[e-/ADU] AD conversion factor'))
        hdr.append(('DET-TMP', detectorTemp, '[degC] Detector temperature'))

        try:
            expTimeSeconds = expTime / 1000.0
            expStart = astropy.time.Time(expStart, format='unix', scale='utc')
            timeCards = timecards.TimeCards(startTime=expStart)
            timeCards.end(expTime=expTimeSeconds)
            cards = timeCards.getCards()

            # MCS is still pyfits based. So convert cards formats, grr.
            hdr.append(('EXPTIME', expTimeSeconds, '[s] commanded exposure time'))
            pyCards = [(c['name'], c['value'], c['comment']) for c in cards]
            hdr.extend(pyCards)
        except Exception as e:
            cmd.warn(f'text="FAILED to gather time cards: {e}"')

        try:
            instCards = self._getInstHeader(cmd)
            hdr.add_comment('Subaru Device Dependent Header Block for PFS-MCS')
            hdr.extend(instCards, bottom=True)
        except Exception as e:
            cmd.warn(f'text="FAILED to gather instrument cards: {e}"')

        return hdr

    def _makeImageHeader(self, cmd):
        """ Create a complete WCS header.

        Notes
        ----
        - We need to get image center from config file.
        - So far, just spit out MCS-to-PFI linear conversion. Later, add SIP terms, and MCS-to-sky.
        - Needs to be pulled out into actorcore or pfs_utils.

        Returns
        -------
        hdr : pyfits.Header instance.

        """

        gen2Keys = self.actor.models['gen2'].keyVarDict

        imageCenter = (8960//2, 5778//2)
        rot = gen2Keys['tel_rot'][1]
        alt = gen2Keys['tel_axes'][1]
        az = gen2Keys['tel_axes'][0]

        cmd.debug(f'text="center={imageCenter} axes={az},{alt},{rot}"')

        hdr = None
        wcs, sip = pfsWcs.WCSParameters('mcs_pfi', imageCenter, rot, alt, az)
        hdr = wcs.to_header()

        return hdr

    def _doExpose(self, cmd, expTime, expType, frameId, mask=None):
        """ Take an exposure and save it to disk. """

        fileIdsQ = self.requestNextFileIds(cmd, frameId)
        cmd.diag(f'text="new exposure"')
        expStart = time.time()
        if self.simulationPath is None:
            filename = '/tmp/scratchFile'
            image = self.actor.camera.expose(cmd, expTime, expType, filename, doCopy=False)
        else:
            imagePath, image, target = self.getNextSimulationImage(cmd)
        cmd.inform(f'text="done: image shape = {image.shape}"')

        try:
            fileIds = fileIdsQ.get(timeout=15.0)
        except queue.Empty:
            cmd.warn('text="failed to get a new filename in time"')

        frameId = fileIds['frameId']
        filename = self.butler.getPath('mcsFile', **fileIds)
        pathdir = filename.parent
        if not os.path.isdir(pathdir):
            os.makedirs(pathdir, 0o755, exist_ok=True)
        cmd.inform(f'text="newpath={filename}"')

        # Now, after getting the filename, get predicted locations
        if self.simulationPath is not None:
            self._writeExpectTarget(cmd, frameId, target)
        else:
            pass

        if mask is not None:
            cmd.inform(f'text="mask image shape: {mask.shape} type:{mask.dtype}"')
            cmd.inform(f'text="image shape: {image.shape} type:{image.dtype}"')
            image = image*mask.astype('uint16')

        hdr = self._constructHeader(cmd, filename, expType, expTime, expStart)
        cmd.diag(f'text="hdr done: {len(hdr)}"')
        phdu = pyfits.PrimaryHDU(header=hdr)

        try:
            imgHdr = self._makeImageHeader(cmd)
        except Exception as e:
            cmd.warn(f'text="FAILED to generate WCS header: {e}"')
            imgHdr = pyfits.Header()
        imgHdu = pyfits.CompImageHDU(image.astype('uint16'), name='IMAGE', compression_type='RICE_1')

        imgHdu.header.extend(imgHdr)
        hduList = pyfits.HDUList([phdu, imgHdu])

        # Patch core FITS card comments to match Subaru requirements.
        imgHdr = imgHdu.header
        imgHdr.set('BZERO', comment='Real=fits-value*BSCALE+BZERO')
        imgHdr.set('BSCALE', comment='Real=fits-value*BSCALE+BZERO')
        imgHdr.append(('BUNIT', 'ADU', 'Unit of original pixel values'))
        imgHdr.append(('BLANK', 32767, 'Value used for NULL pixels'))  # 12-bit camera
        imgHdr.append(('BIN-FCT1', 1, '[pixel] Binning factor of X axis'))
        imgHdr.append(('BIN-FCT2', 1, '[pixel] Binning factor of X axis'))

        pHdr = phdu.header
        pHdr.set('EXTEND', comment='Presence of FITS Extension')

        hduList.writeto(filename, checksum=True, overwrite=True)
        cmd.inform(f'text="write image to filename={filename}"')

        cmd.inform(f'mcsFileIds={fileIds["pfsDay"]},{fileIds["visit"]},{fileIds["frame"]}')
        return filename, image

    def resetThreshold(self, cmd):
        """reset the thresholds"""

        self.findThresh = None
        self.centThresh = None
        cmd.finish('Centroid threshold=d=reset')

    def expose(self, cmd):
        """
        call the sequence of steps for an exposure 
         - take image (or read from disk)
         - if noCentroids not set, do cetnroiding, write to database
         - if doFibreID set, do transformations, fibre identification, write to DB
        """

        cmdKeys = cmd.cmd.keywords

        # Switch from default no centroids to default do centroids
        noCentroidArg = 'noCentroid' in cmdKeys
        if noCentroidArg:
            doCentroid = False
        else:
            doCentroid = True

        doFibreID = 'doFibreID' in cmdKeys
        newField = 'newField' in cmdKeys

        simDot = 'simDot' in cmdKeys
        if simDot:
            dotmask = self._makeDotMask(cmd)
        else:
            dotmask = None

        cmd.inform(f'text="doCentroid= {doCentroid} doFibreID = {doFibreID}')
        
        # get frame ID if explicitly set, otherise reset
        expType = cmdKeys[0].name
        if 'frameId' in cmdKeys:
            frameId = cmdKeys['frameId'].values[0]
        else:
            frameId = None

        # set exposure time
        if expType in ('bias', 'test'):
            expTime = self.expTime
        else:
            expTime = cmd.cmd.keywords["expTime"].values[0] * 1000

        if (expTime != self.expTime):
            self.actor.camera.setExposureTime(cmd, expTime)

        cmd.inform('text="Exposure time now is %d ms." ' % (expTime))
        filename, image = self._doExpose(cmd, expTime, expType, frameId, mask=dotmask)

        if frameId is None:
            filename = pathlib.Path(filename)
            frameId = int(filename.stem[4:], base=10)

        self.handleTelescopeGeometry(cmd, filename, frameId, expTime)

        # set visitID
        self.visitId = frameId // 100
        self.actor.image = image

        # load telescope values from the DB
        cmd.inform(f'text="loading telescope parameters for frame={frameId}"')
        zenithAngle, insRot = dbTools.loadTelescopeParametersFromDB(self._db, int(frameId))

        # get the geometry if it hasn't been loaded yet
        cmd.inform('text="loading geometry"')
        self.getGeometry(cmd)

        # if the centroid flag is set
        if doCentroid:
            
            # connect to DB
            db = self.connectToDB(cmd)

            cmd.inform('text="Setting centroid parameters." ')
            self.setCentroidParams(cmd)

            if self.findThresh is None or frameId % 100 == 0:
                cmd.inform('text="Calculating threshold." ')
                self.calcThresh(cmd, frameId, zenithAngle, insRot, self.centParms)

            cmd.inform('text="Running centroid on current image" ')

            # switch for different centroid methods. Call with switchCMethod
            t1 = time.time()
            
            if self.actor.cameraName == 'rmod_71m':
                self.cMethod = 'sep'
                cmd.inform(f'text="Bench camera RMOD-71M is used. Using SEP" ')

            if(self.cMethod == 'sep'):
                cmd.inform(f'text="Using SExtractor for centroid" ')
                self.runCentroidSEP(cmd)
            else:
                self.runCentroid(cmd, self.centParms)
            
            t2 = time.time()
            cmd.inform(f'text="Centroids done in {t2-t1} second" ')

            # dumpCentroidtoDB
            self.dumpCentroidtoDB(cmd, frameId)
           
            cmd.inform('text="Sending centroid data to database" ')

        # do the fibre identification
        if doFibreID:
        
            cmd.inform('text="zenith angle=%s"'%(zenithAngle))
            cmd.inform('text="instrument rotation=%s"'%(insRot))
            
            # Get last two degits of frameID
            iterNum = frameId % 100
            if iterNum == 0:
                newField = True
                cmd.inform(f'text="New field because iterNum = {iterNum} "')
            else:
                newField = False

            enableEasyID=False

            if enableEasyID:
                #if newField:
                self.establishTransform(cmd, 90-zenithAngle, insRot, frameId)
                
                self.easyFiberID(cmd, frameId)

            else:

                self.establishTransform(cmd, 90-zenithAngle, insRot, frameId)
                if(self.adjacentCobras == None):
                    self.adjacentCobras = mcsTools.makeAdjacentList(self.centrePos, self.armLength)
                    cmd.inform(f'text="made adjacent lists"')

                # fibreID
                self.fibreID(cmd, frameId, zenithAngle, insRot)

        cmd.inform(f'frameId={frameId}; filename={filename}')


        cmd.finish('exposureState=done')

    def dumpCentroidtoDB(self, cmd, frameId):
        """Connect to database and return json string to an attribute."""

        conn = self.connectToDB()
        cmd.diag(f'text="dumping centroids to db {conn}"')

        # The section should be removed, replaced by the createTables command.
        if self.newTable:
            self._makeTables(conn, doDrop=True)
            cmd.inform('text="Centroid table created. "')
        else:
            #self._makeTables(conn, doDrop=False)
            cmd.inform('text="Attaching centroid to exsiting table. "')

        buf = self._writeCentroids(self.centroids, frameId, 1, conn)

        cmd.inform('text="Centroids of exposure ID %08d dumped."' % (frameId))

    def switchFMethod(self, cmd):
        cmdKeys = cmd.cmd.keywords
        self.fMethod = cmdKeys['fMethod'].values[0]
        cmd.inform(f'text="fMethod = {self.fMethod}"')
        cmd.finish('switchFMethod=done')

    def switchCMethod(self, cmd):
        cmdKeys = cmd.cmd.keywords
        self.cMethod = cmdKeys['cMethod'].values[0]
        cmd.inform(f'text="cMethod = {self.cMethod}"')
        cmd.finish('switchCMethod=done')

    def resetGeometry():
        """
        reset the geometry flag. Next call of getGeometry will reload parameters.
        """

        self.geometrySet = False

    def getGeometry(self, cmd):

        db = self.connectToDB(cmd)

        cmd.inform(f'text="getting geometry"')

        if(self.geometrySet == True):
            cmd.inform('text="geometry is already set"')
            return

        # boresight centre in pixels
        self.rotCent = dbTools.loadBoresightFromDB(db, int(self.visitId))
        cmd.inform(f'text="boresight={self.rotCent[0]},{self.rotCent[1]}"')

        # read xmlFile
        instPath = os.path.join(os.environ['PFS_INSTDATA_DIR'])
        #if(self.geomFile == None):
        #    self.geomFile = os.path.join(instPath, 'data/pfi/modules/ALL/ALL_final_20210920_mm.xml')
        if(self.dotFile == None):
            self.dotFile = os.path.join(
                instPath, "data/pfi/dot/black_dots_mm.csv")

        pfi = self.butler.get("moduleXml", moduleName="ALL", version="")
        dots = self.butler.get("black_dots", moduleName="ALL", version="")

        cmd.inform(f'text="loading XML from butler"')
        cmd.inform(f'text="loading DOT location from butler"')
        self.centrePos, self.armLength, self.dotPos, self.goodIdx, self.calibModel = mcsTools.readCobraGeometry(
            pfi, dots)
        
        fids = self.butler.get('fiducials')
        self.outerRingIds, self.badFidIds = mcsTools.readFiducialMasks(fids)
        cmd.inform('text="cobra geometry read"')
        self.geometrySet = True

    def establishTransform(self, cmd, altitude, insrot, frameID):

        # Read fiducial and spot geometry
        fids = self.butler.get('fiducials')

        # need for fibreID later
        self.fids = fids

        db=opdb.OpDB(hostname='db-ics', port=5432,
                   dbname='opdb',
                   username='pfs')
        mcsData = db.bulkSelect('mcs_data',f'select spot_id, mcs_center_x_pix, mcs_center_y_pix '
                f'from mcs_data where mcs_frame_id = {frameID}')

        if 'rmod' in self.actor.cameraName.lower():
            altitude = 90.0
            insrot = 0

        self.logger.info(f'Camera name: {self.actor.cameraName}')
        cmd.inform(f'text="camera name: {self.actor.cameraName} altitude using {altitude}"')
        cmd.inform(f'text="camera name: {self.actor.cameraName} altitude using {insrot}"')


        pfiTransform = transformUtils.fromCameraName(self.actor.cameraName, 
            altitude=altitude, insrot=insrot,nsigma=0, alphaRot=1)

        # these values are now read via mcsToolds.readFiducialMasks

        # set the good fiducials and outer ring fiducials if not yet set
        if(self.fidsGood == None):
            self.fidsOuterRing, self.fidsGood = mcsTools.readFiducialMasks(fids)
        
        #outerRingIds = [29, 30, 31, 61, 62, 64, 93, 94, 95, 96]
        #fidsOuterRing = fids[fids.fiducialId.isin(outerRingIds)]
        #badFids = [1,32,34,61,68,75,88,89,2,4,33,36,37,65,66,67,68,69]
        #goodFids = list(set(fids['fiducialId'].values)-set(badFids))
        #fidsGood = fids[fids.fiducialId.isin(goodFids)]
        
        pfiTransform.updateTransform(mcsData, self.fidsOuterRing, matchRadius=8.0, nMatchMin=0.1)

        nsigma = 0
        pfiTransform.nsigma = nsigma
        pfiTransform.alphaRot = 0

        for i in range(2):
            ffid, dist = pfiTransform.updateTransform(mcsData, self.fidsGood, matchRadius=4.2,nMatchMin=0.1)
        #pfiTransform.updateTransform(mcsData, fids, matchRadius=2.0)

        dbToools.writeTransformToDB(db, frameId, pfiTransform, self.actor.cameraName)
        #dbTools.writeFidToDB(db, ffid, frameID)
        
        self.pfiTrans = pfiTransform


        cmd.inform(f'text="PFI transformation method built"')


    def setDotFile(self, cmd):

        self.geomFile = cmd.cmd.keywords["dotFile"].values[0]
        cmd.inform(f'text="geometry file set to {self.dotFile}"')

    def resetGeometryFile(self, cmd):

        self.geomFile = None
        cmd.inform(f'text="geometry file set to {self.geomFile}"')

    def setGeometryFile(self, cmd):

        self.geomFile = cmd.cmd.keywords["geomFile"].values[0]
        cmd.inform(f'text="geometry file set to {self.geomFile}"')

    def easyFiberID(self, cmd, frameId):
        reload(calculation)
        
        db = self.connectToDB(cmd)
        cmd.inform(f'text="Running easy FibreID"')

        # sorting the spot id, so that the ID returned by matching should 
        # keep the same.
        mcsData = db.bulkSelect('mcs_data','select * from mcs_data where '
                f'mcs_frame_id = {frameId}').sort_values(by=['spot_id'])

        renames = dict(mcs_frame_id='mcsId',
                       spot_id='fiberId',
                       mcs_center_x_pix='centroidx',
                       mcs_center_y_pix='centroidy')

        mcsData.rename(columns=renames, inplace=True)
        df=mcsData.loc[mcsData['fiberId'] > 0]
        
        
        '''
            OK, here we transform all mcs_data to pfi
        '''
        x_mm, y_mm = self.pfiTrans.mcsToPfi(df['centroidx'].values,df['centroidy'].values)
        
        pos=x_mm+y_mm*(1j)
        
        centers=self.calibModel.centers
        cmd.inform(f'text="cobra centers = {centers}" ')

        cmd.inform(f'text="Loading arm-length = {self.calibModel.L1}" ')

        # It should be one spot in patrol region. 
        target = calculation.lazyIdentification(centers, pos, 
            radii=self.calibModel.L1+self.calibModel.L2)
        
        cmd.inform(f'text="Total ID target = {target}" ')


        cmd.inform(f'text="Loading DOT file for missing IDs" ')

        dotData = pd.read_csv(self.dotFile, delimiter = ",", header=0)
        dotCenter = dotData['x'].values+dotData['y'].values*1j

        mpos = np.zeros(len(target), dtype=complex)
        for n, k in enumerate(target):
            if k < 0:
                # If the target failed to match, we think it is highly possible in dot
                mpos[n] = dotCenter[n]
            else:
                mpos[n] = pos[k]
        indx = np.where(target < 0 )

        visitId = frameId // 100
        iteration = frameId % 100
        cobraTarget = db.bulkSelect('cobra_target','select pfs_visit_id from cobra_target where '
            f'(pfs_visit_id = {visitId}) AND iteration = {iteration}').reset_index()
        db.close()

        db = self.connectToDB(cmd)
        if len(cobraTarget) == 0:
            cmd.inform(f'text="Fall back using cobra centers as target." ')
            dbTools.writeFakeTargetToDB(db, self.calibModel.centers, int(frameId))
    
        cobraMatch = np.zeros((2394, 5))
        cobraMatch[:,0] = np.arange(2394)+1 
        cobraMatch[:,1] = target+1
        cobraMatch[:,2] = mpos.real
        cobraMatch[:,3] = mpos.imag
        cobraMatch[indx,4] = 1
        
        dbTools.writeMatchesToDB(db, cobraMatch, int(frameId))


    def fibreID(self, cmd, frameId, zenithAngle, insRot):

        db = self.connectToDB(cmd)

        # if the iteration is the first, the previous position is the home position
        if(frameId % 100 == 0):
            self.prevPos = self.centrePos

        self.mmCentroids = np.copy(self.centroids)
        #transform the coordinates to mm in place
        self.mmCentroids[:,1], self.mmCentroids[:,2] = self.pfiTrans.mcsToPfi(self.centroids[:,1],self.centroids[:,2])

        # if the method is target, load from database, otherwise the target = previous position
        if(self.fMethod == 'target'):
            # load target positions
            tarPos = dbTools.loadTargetsFromDB(db, int(frameId))
            db.close()
            
            cmd.inform(f'text="loaded {len(tarPos)} targets from DB"')
            if(len(tarPos)==0):
                db = self.connectToDB(cmd)
                dbTools.writeFakeTargetToDB(db, self.calibModel.centers, int(frameId))
                db.close()

                visitId = frameId // 100
                iteration = frameId % 100
                cmd.inform(f'text="Fall back using cobra centers as target." ')
                cmd.inform(f'text="Writing minimal information to target database."')

                tarPos = self.prevPos
        else:
            db.close()
            
            tarPos = dbTools.loadTargetsFromDB(db, int(frameId))
            if(len(tarPos)==0):
                db = self.connectToDB(cmd)
                dbTools.writeFakeTargetToDB(db, self.calibModel.centers, int(frameId))
                db.close()
            tarPos = self.prevPos
            
        # do the identification
        cmd.inform(f'text="Starting Fiber ID"')
        t0 = time.time()
        cobraMatch, unaPoints = mcsTools.fibreId(self.mmCentroids, self.centrePos, self.armLength, tarPos,
                                      self.fids, self.dotPos, self.goodIdx, self.adjacentCobras)
        t1 = time.time()
        cmd.inform(f'text="Fiber ID finished in {t1-t0:0.2f}s"')


        db = self.connectToDB(cmd)
        dbTools.writeMatchesToDB(db, cobraMatch, int(frameId))
        db.close()
        cmd.inform(f'text="wrote matched cobras to database"')

        # save the values to the previous position
        self.prevPos = cobraMatch[:, [0, 2, 3]]

    def handleTelescopeGeometry(self, cmd, filename, frameId, expTime):
        if self.simulationPath is None:
            # We are live: use Gen2 telescope info.
            gen2Model = self.actor.models['gen2'].keyVarDict

            axes = gen2Model['tel_axes'].getValue()
            az, alt, *_ = axes

            rot = gen2Model['tel_rot'].getValue()
            posAngle, instrot = rot
            startTime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        else:

            # We are reading images from disk: get the geometry from the headers.
            simPath = str(filename)
            cmd.inform(f'text="filename= {simPath}"')

            simHdr = fitsio.read_header(str(simPath), 0)
            cmd.inform('text="loaded telescope info from %s"' % (simPath))

            az = simHdr.get('AZIMUTH', -9998.0)
            alt = simHdr.get('ALTITUDE', -9998.0)

            if az is None:
                az = -9998.0
            if alt is None:
                alt = -9998.0

            expTime = simHdr.get('EXPTIME', -9998.0)
            instrot = simHdr.get('INR-STR', -9998.0)

            # Redefine instrot to be 0.5 since we know this fact from ASRD test.
            instrot = 0.5
            startTime = simHdr.get('UTC-STR', None)

            if startTime is None:
                ctime = os.stat(filename).st_ctime
                startTime = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(ctime))
                cmd.warn(f'text="no start card in {simPath}, using file date: {startTime}"')

        visitId = frameId // 100
        cmd.inform(f'text="frame={frameId} visit={visitId} az={az} alt={alt} instrot={instrot}"')
        # Packing information into data structure
        telescopeInfo = {'frameid': frameId,
                         'visitid': visitId,
                         'starttime': startTime,
                         'exptime': expTime,
                         'altitude': alt,
                         'azimuth': az,
                         'instrot': instrot}

        self._writeTelescopeInfo(cmd, telescopeInfo)

    def calcThresh(self, cmd, frameId, zenithAngle, insRot, centParms):
        """  Calculate thresholds for finding/centroiding from image 
        """        
        image = self.actor.image

        self.findThresh, self.centThresh, self.avBack = mcsTools.getThresh(
            image, self.rotCent, self.centParms['threshSigma'], self.centParms['findSigma'], self.centParms['centSigma'])

        a1 = self.centParms['threshSigma']
        a2 = self.centParms['findSigma']
        a3 = self.centParms['centSigma']
        cmd.inform(f'text="findThresh ={self.findThresh:.2}, centThresh = {self.centThresh:.2}"')
        cmd.inform(f'text="threshSigma={a1} findSigma={a2} centSigma={a3}"')

    def setCentroidParams(self, cmd):
        """
        top level routine for setting centroid parameters. REads the defaults from teh config fil,e
        then changes any specified in the keywords argument. 

        """

        self.centParms = mcsTools.getCentroidParams(cmd)
        self.centParms['findSigma']=50
        self.centParms['centSigma']=50

        self.logger.info(f'centParms: {self.centParms}')

    def runCentroidSEP(self, cmd):
        cmdKeys = cmd.cmd.keywords
        self.newTable = "newTable" in cmdKeys

        cmd.debug('text="newTable value = %s"' % (self.newTable))

        image = self.actor.image

        cmd.inform(f'state="measuring cached image: {image.shape}"')
        
        # Run centroid

        bkg = sep.Background(image.astype(float), bw=64, bh=64, fw=3, fh=3)
        bkg_image = bkg.back()

        data_sub = image - bkg

        #sigma = np.std(data_sub)
        centroids = sep.extract(data_sub.astype(float), 20 , err=bkg.globalrms,
            filter_type='conv', minarea=10)
        
        
        npoint = centroids.shape[0]
        tCentroids = np.zeros((npoint, 8))

        tCentroids[:, 0] = np.arange(npoint)+1
        tCentroids[:, 1] = centroids['x']
        tCentroids[:, 2] = centroids['y']
        tCentroids[:, 3] = centroids['x2']
        tCentroids[:, 4] = centroids['y2']
        tCentroids[:, 5] = centroids['xy']
        tCentroids[:, 6] = centroids['thresh']
        tCentroids[:, 7] = centroids['peak']

        self.centroids = tCentroids
        self.nCentroid = len(centroids)
        cmd.inform('text="%d centroids"' % (len(centroids)))
        cmd.inform('state="centroids measured"')

    def runCentroid(self, cmd, centParms):
        """ Measure centroids on the last acquired image. """

        cmdKeys = cmd.cmd.keywords
        self.newTable = "newTable" in cmdKeys

        cmd.debug('text="newTable value = %s"' % (self.newTable))

        image = self.actor.image

        cmd.inform(f'text="centroidParms: findThresh {self.findThresh},centThresh {self.centThresh}"')
        cmd.inform(f'centroidParms: nmin {centParms["nmin"]}, maxIt {centParms["maxIt"]} boxFind {centParms["boxFind"]} boxCent {centParms["boxCent"]}')

        
        cmd.inform(f'state="measuring cached image: {image.shape}"')


        # get hte active region, based on the current boresight
        #xmin = self.rotCent[0] + centParms['activeX']
        #ymin = self.rotCent[1] + centParms['activeY']
        #xmax = self.rotCent[0] + centParms['activeX']
        #ymax = self.rotCent[1] + centParms['activeY']
        
        # crop the region to the active image section
        
        #subImage = image[xmin:xmax,ymin:ymax]
        a = centroid.centroid_only(image.astype('<i4'),
                                   centParms['fwhmx'], centParms['fwhmy'], self.findThresh, self.centThresh,
                                   centParms['boxFind'], centParms['boxCent'],
                                   centParms['nmin'], centParms['maxIt'], 0)

        centroids = np.frombuffer(a, dtype='<f8')
        centroids = np.reshape(centroids, (len(centroids)//7, 7))
                         

        # adjust the coordinates back to global values
        #centroids[:,1] += xmin
        #centroids[:,2] += ymin
        
        # TEMPORARY FILTER STUFF FOR STRAY LIGHT DURING COMMISSIONING RUN, NEED TO FIX THIS
        # MORE ELEGANTLY ONCE DB SCHEMA HAS FLAGS COLUMN IN MCS_DATA??!!!

        # get rid of overly large points
        #ind=np.where(centroids[:,2] > 15)
        #centroids=centroids[ind].squeeze()

        nSpots = centroids.shape[0]
        points = np.empty((nSpots, 8))

        # ADD A PLUS 1 TO MATCH THE OTHER CENTROIDING AND STOP CAUSING INDEXING ERRORS
        points[:, 0] = np.arange(nSpots)+1
        points[:, 1:] = centroids[:, 0:]

        points[:,-1]=np.repeat(self.avBack,len(points))

        # Swap last two fields
        #points[:,[-2,-3]] = points[:,[-3,-2]]
        points[:,[-1,-2]] = points[:,[-2,-1]]
        points[:,[-1,-3]] = points[:,[-3,-1]]

        self.centroids = points

        self.nCentroid = len(points)
        cmd.inform('text="%d centroids"' % (len(centroids)))
        cmd.inform('state="centroids measured"')

    def _writeTelescopeInfo(self, cmd, telescopeInfo, conn=None):

        # Let the database handle the primary key
        db = self.connectToDB(cmd)
        res = db.session.execute('select * FROM "mcs_exposure" where false')
        colnames = tuple(res.keys())
        realcolnames = colnames[0:]

        """
          TODO: Those are the fake values for making PFI to work now, adding actual code later
        """
        adc_pa = 0
        dome_temperature = 5
        dome_pressure = 101
        dome_humidity = 0.3
        outside_temperature = 5
        outside_pressure = 101
        outside_humidity = 0.3
        mcs_cover_temperature = 5
        mcs_m1_temperature = 6
        taken_at = "'"+telescopeInfo['starttime']+"'"
        taken_in_hst_at = "'"+telescopeInfo['starttime']+"'"
        if self.actor.cameraName == 'canon_50m':
            mcs_camera_id = 0
        if self.actor.cameraName == 'rmod_71m':
            mcs_camera_id = 1
        if self.cMethod == 'win':
            measurement_algorithm = "'windowed'"
        if self.cMethod == 'sep':
            measurement_algorithm = "'sep'"
        version_actor = '0'
        version_instdata = '0'

        line = '%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s,%s,%d,%s,%s,%s' % (telescopeInfo['frameid'],
                                                                       telescopeInfo['visitid'],
                                                                       telescopeInfo['exptime']/1000.0,
                                                                       telescopeInfo['altitude'],
                                                                       telescopeInfo['azimuth'],
                                                                       telescopeInfo['instrot'],
                                                                       adc_pa, dome_temperature,
                                                                       dome_pressure, dome_humidity,
                                                                       outside_temperature, outside_pressure,
                                                                       outside_humidity,
                                                                       mcs_cover_temperature,mcs_m1_temperature,
                                                                       taken_at, taken_in_hst_at,
                                                                       mcs_camera_id,
                                                                             measurement_algorithm,
                                                                             version_actor, version_instdata)
                                                                        

        buf = io.StringIO()
        buf.write(line)
        buf.seek(0, 0)

        self._writeData('mcs_exposure', realcolnames, buf)
        buf.seek(0, 0)

        cmd.inform('text="Telescope information for frame %s populated."' % (telescopeInfo['frameid']))

        return buf

    def _writeData(self, tableName, columnNames, dataBuf):
        """Wrap a direct COPY_FROM via sqlalchemy. """

        columns = ','.join('"{}"'.format(k) for k in columnNames)
        sql = 'COPY {} ({}) FROM STDIN WITH CSV'.format(
            tableName, columns)

        #try:
        #    db = self.connectToDB(None)
        #    with db.engine.connect() as conn:
        #        with conn.connection.cursor() as cursor:
        #            cursor.copy_expert(sql, dataBuf)
        #except Exception as e:
        #    self.logger.warn(f"failed to write with {sql}: {e}")
        
        try:
            db = self.connectToDB(None)
            session = db.session
            with session.connection().connection.cursor() as cursor:
                cursor.copy_expert(sql, dataBuf)
                cursor.close()
            session.execute('commit')
        except Exception as e:
            self.logger.warn(f"failed to write with {sql}: {e}")

    def _readData(self, sql):
        """Wrap a direct COPY_TO via sqlalchemy. """

        dataBuf = io.StringIO()

        try:
            db = self.connectToDB(None)
            session = db.session
            with session.connection().connection.cursor() as cursor:
                cursor.copy_expert(sql, dataBuf)
            dataBuf.seek(0, 0)
            return dataBuf
        except Exception as e:
            self.logger.warn(f"failed to read with {sql}: {e}")

    def _writeCentroids(self, centArr, frameId, moveId, conn=None):
        """ Write all measurements for a given (frameId, moveId) """

        # if too many centroids have been returned, only save the first 5000
        
        nItems = len(centArr)
        if(nItems > 5000):
            nItems = 5000
            
        buf = io.StringIO()
        for l_i in range(nItems):
            line = '%d,%d,%f,%f,%f,%f,%f,%f,%f\n' % (frameId, l_i+1, centArr[l_i,1], 
                                        centArr[l_i,2], centArr[l_i,3], centArr[l_i,4], 
                                        centArr[l_i,5], centArr[l_i,6], centArr[l_i,7])
            
            buf.write(line)
        
        line = line = '%d,%d,%f,%f,%f,%f,%f,%f,%f\n' % (frameId, -1, np.nan, 
                                        np.nan, np.nan, np.nan, 
                                        np.nan, np.nan, np.nan)
        buf.write(line)
        buf.seek(0, 0)


        # Let the database handle the primary key
        db = self.connectToDB(None)
        colnames = db.session.execute('select * FROM "mcs_data" where false')
        realcolnames = tuple(colnames.keys())[0:]


        self._writeData('mcs_data', realcolnames, buf)
        buf.seek(0, 0)
        return buf

    def _readCentroids(self, conn, frameId, moveId):
        """ Read all measurements for a given (frameId, moveId)"""

        if conn is None:
            conn = self.connectToDB()

        if self.simulationPath is None:
            cmd = """copy (select * from mcsPerFiber where frameId={frameId} and moveId={moveId}) to stdout delimiter ',' """
            buf = self._readData(cmd)

            # Skip the frameId, etc. columns.
            arr = np.genfromtxt(buf, dtype='f4',
                                delimiter=',', usecols=range(4, 24))
        else:
            cmd = f"""copy (select * from 'mcsData' where frameId={frameId} and moveId={moveId}) to stdout delimiter ',' """
            buf = self._readData(cmd)

            # Skip the frameId, etc. columns.
            arr = np.genfromtxt(buf, dtype='f4',
                                delimiter=',', usecols=range(4, 8))

        return arr

    def _makeDotMask(self, cmd):
        """ Make dot mask for simulation or proessing purpose)"""

        # read dot location
        dotfile = '/home/pfs/mhs/devel/pfs_instdata/data/pfi/dot/dot_measurements_20210428_el90_rot+00_ave.csv'
        dotpos = pd.read_csv(dotfile)

        unit_height = 30
        unit_weight = 30
        r_mean = np.around(np.mean(dotpos['r_tran'].values)).astype('int')
        cmd.inform('text="Making dot image"')

        mask = self._create_circular_mask(unit_height, unit_weight, radius=10)

        masked_img = np.zeros([unit_height, unit_weight])+1

        masked_img[mask] = 0
        # To-do: change here for the image size
        dotmask = np.zeros([7096, 10000])+1

        for i in range(len(dotpos)):
            xstart = np.around(dotpos['x_tran'].values[i]-(unit_weight/2)).astype('int')

            ystart = np.around(dotpos['y_tran'].values[i]-(unit_height/2)).astype('int')

            dotmask[ystart:ystart+unit_height, xstart:xstart+unit_weight] = masked_img

        return dotmask

    def _create_circular_mask(self, h, w, center=None, radius=None):

        if center is None:  # use the middle of the image
            center = (int(w/2), int(h/2))
        if radius is None:  # use the smallest distance between the center and image walls
            radius = min(center[0], center[1], w-center[0], h-center[1])

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

        mask = dist_from_center <= radius
        return mask

    def _writeExpectTarget(self, cmd, frameId, targets):
        '''
        Write the expect target to databse.
        '''

        mcs_f3c_model = {'camMatrix': np.array([[8.36047142e+03, 0.00000000e+00, 5.01140651e+03],
                                                [0.00000000e+00, 8.75498927e+03, 3.53572485e+03],
                                                [0.00000000e+00, 0.00000000e+00, 1.00000000e+00]]),
                         'camDistor': np.array([[0.04236098, -0.11608623, 0.0012094, -0.00031588, 0.08034843]]),
                         'camRotVec': np.array([[-0.02458225],
                                                [-0.06933687],
                                                [-0.01775009]]),
                         'camTranVec': np.array([[-72269.86663875],
                                                 [-48804.84781669],
                                                 [112015.99005353]])}

        cobra_obj = np.array([targets.real, targets.imag, np.zeros(len(targets))]).T

        imgpoints2, _ = cv2.projectPoints(cobra_obj.astype(np.float32),
                                          mcs_f3c_model['camRotVec'], mcs_f3c_model['camTranVec'],
                                          mcs_f3c_model['camMatrix'], mcs_f3c_model['camDistor'])
        imgarr2 = imgpoints2[:, 0, :]
        pfi_x = imgarr2[:, 0]
        pfi_y = imgarr2[:, 1]

        db = self.connectToDB(None)
        colnames = db.session.execute('select * FROM "cobra_target" where false')
        realcolnames = tuple(colnames.keys())[0:]

        buf = io.StringIO()

        for i in range(len(pfi_x)):
            line = '%d,%d,%d,%d, %f,%f,%f,%f,%d,%d,%d,%f,%d,%d,%d,%f,%d\n'%(frameId/100, 1, i+1, 99, pfi_x[i], pfi_y[i],
                                                                            pfi_x[i], pfi_y[i], 1, 1, 99, 99, 1, 1, 99, 99, 0)
            buf.write(line)
        buf.seek(0, 0)

        self._writeData('cobra_target', realcolnames, buf)
        buf.seek(0, 0)
