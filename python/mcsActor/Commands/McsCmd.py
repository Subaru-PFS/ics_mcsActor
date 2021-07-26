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

import os
import astropy.io.fits as pyfits
import fitsio

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from actorcore.utility import fits as fitsUtils
from opscore.utility.qstr import qstr

import pfs.utils.coordinates.MakeWCSCards as pfsWcs
from scipy.spatial import cKDTree

import psycopg2
import psycopg2.extras
from xml.etree.ElementTree import dump

import mcsActor.windowedCentroid.centroid as centroid
import mcsActor.mcsRoutines.mcsRoutines as mcsToolsNew
import mcsActor.mcsRoutines.dbRoutinesMCS as dbTools

from importlib import reload
reload(dbTools)
reload(mcsToolsNew)
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

        # self.setCentroidParams(None)
        self.adjacentCobras = None
        self.geometrySet = False
        self.geomFile = None
        self.dotFile = None
        self.fibreMode = 'full'

        logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S")
        self.logger = logging.getLogger('mcscmd')
        self.logger.setLevel(logging.INFO)

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
             'object <expTime> [<frameId>] [@noCentroid] [@doCentroid] [@doFibreID] [@simDot]', self.expose),
            ('runCentroid', '[@newTable]', self.runCentroid),
            #('runFibreID', '[@newTable]', self.runFibreID),
            ('reconnect', '', self.reconnect),
            ('resetThreshold', '', self.resetThreshold),
            ('setCentroidParams', '[<fwhmx>] [<fwhmy>] [<boxFind>] [<boxCent>] [<nmin>] [<nmax>] [<maxIt>]',
             self.setCentroidParams),
            ('calcThresh', '[<threshMethod>] [<threshSigma>] [<threshFact>]', self.calcThresh),
            ('simulate', '<path>', self.simulateOn),
            ('simulate', 'off', self.simulateOff),
            ('switchFibreMode', '<fibreMode>', self.switchFibreMode),
            ('resetGeometry', '', self.resetGeometry),
            ('resetGeometryFile', '<geomFile>', self.resetGeometryFile)
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time, seconds"),
                                        keys.Key("expType", types.String(), help="The exposure type"),
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
                                        keys.Key("nmax", types.Int(), help="max number of points for spot"),
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
                                        keys.Key("threshMethod", types.String(),
                                                 help="method for thresholding"),
                                        keys.Key("threshMode", types.Float(), help="mode for threshold"),
                                        keys.Key("geomFile", types.String(), help="file for geometry"),
                                        keys.Key("dotFile", types.String(), help="file for dot information"),
                                        keys.Key("fieldID", types.String(),
                                                 help="fieldID for getting instrument parameters"),
                                        keys.Key("fibreMode", types.String(),
                                                 help="flag for testing different inputs")
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
        self.actor.connectCamera(cmd, camera='rmod')
        self.actor.camera.setExposureTime(cmd, self.expTime)

        cmd.finish('text="Camera connected!"')

    def status(self, cmd):
        """Report status and version; obtain and send current data"""

        self.actor.sendVersionKey(cmd)
        self.actor.camera.sendStatusKeys(cmd)
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd, self.expTime)

        cmd.inform('text="MCS camera present!"')
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

    def requestNextFilename(self, cmd, frameId):
        """ Return a queue which will eventually contain a filename. """

        q = queue.Queue()
        cmd.inform('text="frameIds = %s." '%{frameId})

        def worker(q=q, cmd=cmd):
            filename = self.getNextFilename(cmd, frameId)
            q.put(filename)
            q.task_done()

        task = threading.Thread(target=worker,
                                name='visitFetcher', daemon=True)
        task.start()

        return q

    def getNextFilename(self, cmd, frameId):
        """ Fetch next image filename. 

        In real life, we will instantiate a Subaru-compliant image pathname generating object.  

        """

        if frameId is None:
            ret = self.actor.cmdr.call(actor='gen2', cmdStr='getVisit caller=mcs', timeLim=15.0)
            if ret.didFail:
                raise RuntimeError("getNextFilename failed getting a visit number in 15s!")
            visit = self.actor.models['gen2'].keyVarDict['visit'].valueList[0]
            frameId = visit * 100

            cmd.inform(f'text="gen2.getVisit = {visit}"')
        else:
            if False:           # Need to make this async for safety.
                ret = self.actor.cmdr.call(actor='gen2', cmdStr='updateTelStatus caller=mcs', timeLim=15.0)
                if ret.didFail:
                    raise RuntimeError("getNextFilename failed getting a visit number in 15s!")

        # CPL -- switch to some ics.utils function
        path = os.path.join('/data/raw', time.strftime('%Y-%m-%d', time.gmtime()), 'mcs')
        path = os.path.expandvars(os.path.expanduser(path))
        if not os.path.isdir(path):
            os.makedirs(path, 0o755, exist_ok=True)

        newpath = os.path.join(path, 'PFSC%08d.fits' % (frameId))
        cmd.inform(f'text="newpath={newpath}"')

        return newpath

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

    def _constructHeader(self, cmd, filename, expType, expTime):
        if expType == 'bias':
            expTime = 0.0

        # Subaru rules
        if expType == 'object':
            expType = 'acquisition'

        frameId = os.path.splitext(os.path.basename(filename))[0]

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
            instCards = self._getInstHeader(cmd)
            hdr.add_comment('Subaru Device Dependent Header Block for PFS-MCS')
            hdr.extend(instCards, bottom=True)
        except Exception as e:
            cmd.warn(f'text="FAILED to gather instrument cards: {e}"')

        return hdr

    def _writeCentroidsToDB(self, cmd, frameId):
        """write centroids to database"""

        db = self.connectToDB(cmd)
        dbTools.writeCentroidsToDB(db, self.centroids, int(frameId))
        cmd.inform(f'text="Centroids of exposure ID {frameId} populated"')

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

        nameQ = self.requestNextFilename(cmd, frameId)
        cmd.diag(f'text="new exposure"')
        if self.simulationPath is None:
            filename = 'scratchFile'
            image = self.actor.camera.expose(cmd, expTime, expType, filename, doCopy=False)
        else:
            imagePath, image, target = self.getNextSimulationImage(cmd)
        cmd.inform(f'text="done: image shape = {image.shape}"')

        try:
            filename = nameQ.get(timeout=5.0)
        except queue.Empty:
            cmd.warn('text="failed to get a new filename in time"')

        frameId = int(pathlib.Path(filename).stem[4:], base=10)

        # Now, after getting the filename, get predicted locations
        if self.simulationPath is not None:
            self._writeExpectTarget(cmd, frameId, target)
        else:
            pass

        self.logger.info(f'read filename: {filename}')

        if mask is not None:
            cmd.inform(f'text="mask image shape: {mask.shape} type:{mask.dtype}"')
            cmd.inform(f'text="image shape: {image.shape} type:{image.dtype}"')
            image = image*mask.astype('uint16')

        hdr = self._constructHeader(cmd, filename, expType, expTime)
        cmd.diag(f'text="hdr done: {len(hdr)}"')
        phdu = pyfits.PrimaryHDU(header=hdr)

        try:
            imgHdr = self._makeImageHeader(cmd)
        except Exception as e:
            cmd.warn(f'text="FAILED to generate WCS header: {e}"')
            imgHdr = pyfits.Header()
        imgHdu = pyfits.CompImageHDU(image.astype('uint16'), name='IMAGE', compression_type='RICE_1')

        #imgHdu = pyfits.PrimaryHDU(image.astype('uint16'))
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

        hduList.writeto(filename, checksum=False, overwrite=True)
        cmd.inform(f'filename={filename}')
        cmd.inform(f'text="write image to filename={filename}"')

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

        simDotArg = 'simDot' in cmdKeys
        if simDotArg:
            simDot = True
            dotmask = self._makeDotMask(cmd)
        else:
            simDot = False

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
        if simDot is True:
            filename, image = self._doExpose(cmd, expTime, expType, frameId, mask=dotmask)
        else:
            filename, image = self._doExpose(cmd, expTime, expType, frameId)

        if frameId is None:
            filename = pathlib.Path(filename)
            frameId = int(filename.stem[4:], base=10)

        self.handleTelescopeGeometry(cmd, filename, frameId, expTime)

        # set visitID
        self.visitId = frameId // 100
        self.actor.image = image
        cmd.inform(f'frameId={frameId}')
        cmd.inform(f'filename={filename}')

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

            if self.findThresh is None:
                cmd.inform('text="Calculating threshold." ')
                self.calcThresh(cmd, frameId, zenithAngle, insRot, self.centParms)

            cmd.inform('text="Running centroid on current image" ')
            self.runCentroidSEP(cmd)

            # self.runCentroid(cmd,self.centParms)

            # if the threshold has changed, recalculate: this was added during commissioning run
            if (self.nCentroid < 2000):
                self.calcThresh(cmd, frameId, zenithAngle, insRot, self.centParms)
                self.runCentroid(cmd, self.centParms)

            cmd.inform('text="Sending centroid data to database" ')
            self.dumpCentroidtoDB(cmd, frameId)

        # do the fibre identification
        if doFibreID:

            cmd.inform('text="zenith angle=%s"'%(zenithAngle))
            cmd.inform('text="instrument rotation=%s"'%(insRot))

            # read FF from the database, get list of adjacent fibres if they haven't been calculated yet.
            if(self.adjacentCobras == None):
                adjacentCobras = mcsToolsNew.makeAdjacentList(self.centrePos[:, 1:3], self.armLength)
                cmd.inform(f'text="made adjacent lists"')
                self.fidPos = dbTools.loadFiducialsFromDB(db)
                cmd.inform(f'text="loaded fiducial fibres"')

            # transform centroids to MM
            self.transformations(cmd, frameId, zenithAngle, insRot)

            # fibreID
            self.fibreID(cmd, frameId, zenithAngle, insRot)

        cmd.inform('text="filename=%s"'%(filename))

        # if doFibreID:
        #self.runFibreID(cmd, doFinish=False)
        #self.dumpCentroidtoDB(cmd, frameId)

        cmd.finish('exposureState=done')

    def dumpCentroidtoDB(self, cmd, frameId):
        """Connect to database and return json string to an attribute."""

        conn = self.conn
        cmd.diag(f'text="dumping centroids to db {conn}"')

        # The section should be removed, replaced by the createTables command.
        if self.newTable:
            self._makeTables(conn, doDrop=True)
            cmd.inform('text="Centroid table created. "')
        else:
            #self._makeTables(conn, doDrop=False)
            cmd.inform('text="Attaching centroid to exsiting table. "')

        buf = self._writeCentroids(self.centroids, 1, frameId, 1, conn)

        cmd.inform('text="Centroids of exposure ID %08d dumped."' % (frameId))

    def switchFibreMode(self, cmd):
        cmdKeys = cmd.cmd.keywords
        self.fibreMode = cmdKeys['fibreMode'].values[0]
        cmd.inform(f'text="fibreMode = {self.fibreMode}"')
        cmd.finish('switchFibreMode=done')

    def resetGeometry():
        """
        reset the geometry flag. Next call of getGeometry will reload parameters.
        """

        self.geometrySet = False

    def getGeometry(self, cmd):

        db = self.connectToDB(cmd)
        cmd.inform(f'text="getting geometry for mode {self.fibreMode}"')

        # three modes here - asrd is in pixels, full is in mm, and comm is in mm with fake cobra geometry becuase
        # there are no cobras

        if(self.geometrySet == True):
            cmd.inform('text="geometry is already set"')
            return

        if(self.fibreMode == "asrd"):

            """
            asrd mode, all in pixels, no fiducial fibres
            """

            #not used
            self.offset = [0, 0]
            cmd.inform(f'text="offset={self.offset[0]},{self.offset[1]}"')

            # boresight centre in pixels
            self.rotCent = dbTools.loadBoresightFromDB(db, int(self.visitId))
            cmd.inform(f'text="boresight={self.rotCent[0]},{self.rotCent[1]}"')

            # read xmlFile
            #where is inst_data

            instPath = os.path.join(os.environ['PFS_INSTDATA_DIR'])
            if(self.geomFile == None):
                self.geomFile = os.path.join(instPath, "data/pfi/modules/ALL/ALL_final.xml")
            if(self.dotFile == None):
                self.dotFile = os.path.join(
                    instPath, "data/pfi/dot/dot_measurements_20210428_el30_rot+00_ave.csv")

            cmd.inform(f'text="{instPath} {self.geomFile} {self.dotFile}"')

            # xmlFile="/Users/karr/Science/PFS/cobraData/Full2D/20210219_002/output/ALL_new.xml"
            # dotFile="/Users/karr/software/mhs/products/DarwinX86/pfs_instdata/1.0.1/data/pfi/dot_measurements_20210428_el90_rot+00_ave.csv"
            cmd.inform(f'text="reading geometry from {self.geomFile} {self.dotFile}"')

            self.centrePos, self.armLength, self.dotPos, self.goodIdx = mcsToolsNew.readCobraGeometry(
                self.geomFile, self.dotFile)
            cmd.inform('text="cobra geometry read"')

        elif(self.fibreMode == "full"):

            # check this value
            self.offset = [0, 0]
            cmd.inform(f'text="offset={self.offset[0]},{self.offset[1]}"')

            # boresight centre in pixels
            self.rotCent = dbTools.loadBoresightFromDB(db, int(self.visitId))
            cmd.inform(f'text="boresight={self.rotCent[0]},{self.rotCent[1]}"')

            # read xmlFile
            instPath = os.path.join(os.environ['PFS_INSTDATA_DIR'])
            if(self.geomFile == None):
                self.geomFile = os.path.join(instPath, "data/pfi/modules/ALL/ALL_final.xml")
            if(self.dotFile == None):
                self.dotFile = os.path.join(
                    instPath, "data/pfi/dot/dot_measurements_20210428_el30_rot+00_ave.csv")

            cmd.inform(f'text="reading geometry from {self.geomFile} {self.dotFile}"')
            self.centrePos, self.armLength, self.dotPos, self.goodIdx = mcsToolsNew.readCobraGeometry(
                self.geomFile, self.dotFile)
            cmd.inform('text="cobra geometry read"')

        elif(self.fibreMode == "comm"):

            """
            commissioning mode, full version with fake fiducial fibres, no cobra movemnet, fake arms
            """

            # offset of pinhole mask from center of instrument
            self.offset = [0, -85]

            cmd.inform(f'text="offset={self.offset[0]},{self.offset[1]}"')

            # boresight centre, in pixels
            self.rotCent = dbTools.loadBoresightFromDB(db, int(self.visitId))

            cmd.inform(f'text="boresight={self.rotCent[0]},{self.rotCent[1]}"')

            # get fake geometry
            self.centrePos, self.armLength, self.dotPos, self.goodIdx = mcsToolsNew.readCobraGeometryFake()
            cmd.inform('text="cobra geometry read"')

        cmd.inform('text="fiducial fibres read"')

    def setDotFile(self, cmd):

        self.geomFile = cmd.cmd.keywords["dotFile"].values[0]
        cmd.inform(f'text="geometry file set to {self.dotFile}"')

    def resetGeometryFile(self, cmd):

        self.geomFile = None
        cmd.inform(f'text="geometry file set to {self.geomFile}"')

    def setGeometryFile(self, cmd):

        self.geomFile = cmd.cmd.keywords["geomFile"].values[0]
        cmd.inform(f'text="geometry file set to {self.geomFile}"')

    def transformations(self, cmd, frameId, zenithAngle, insRot):

        # two caes here, the full mm version, and the asrd pixel version, with no ff
        cmd.inform(f'text="fibreMode {self.fibreMode}"')

        if(self.fibreMode == 'comm'):
            fieldElement = False
        if(self.fibreMode == 'full'):
            fieldElement = True

        if(self.fibreMode in ('comm', 'full')):
            db = self.connectToDB(cmd)
            cmd.inform(f'text="tests centroid shape {self.centroids[:,1:3].shape}"')
            centroidsMM = mcsToolsNew.transformToMM(
                self.centroids, self.rotCent, self.offset, zenithAngle, fieldElement, insRot, pixScale=0)
            np.save("centroidsMM.npy", centroidsMM)

            # match FF to the transformed centroids
            nFid = len(self.fidPos[:, 0])
            matchPoint = mcsToolsNew.nearestNeighbourMatching(centroidsMM, self.fidPos, nFid)
            cmd.inform(f'text="fiducial fibres matched"')

            afCoeff, xd, yd, sx, sy, rotation = mcsToolsNew.calcAffineTransform(matchPoint, self.fidPos)
            dbTools.writeAffineToDB(db, afCoeff, int(frameId))

            cmd.inform(f'text="transform calculated {xd} {yd} {sx} {sy} {rotation}"')

            self.centroidsMMTrans = mcsToolsNew.applyAffineTransform(centroidsMM, afCoeff)
            cmd.inform(f'text="affine applied to centroids"')

        elif(self.fibreMode == 'asrd'):

            self.centroidsMMTrans = self.centroids

    def fibreID(self, cmd, frameId, zenithAngle, insRot):

        # if the iteration is the first, the previous position is the home position
        if(frameId % 100 == 1):
            self.prevPos = self.centrePos

        db = self.connectToDB(cmd)

        tarPos = dbTools.loadTargetsFromDB(db, int(frameId))
        cmd.inform(f'text="loaded target postitions from DB"')
        cobraMatch = mcsToolsNew.fibreID(self.centroidsMMTrans, tarPos, self.centrePos,
                                         self.armLength, self.dotPos, self.adjacentCobras, self.goodIdx)
        cmd.inform(f'text="identified fibres"')
        dbTools.writeMatchesToDB(db, cobraMatch, int(frameId))

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

        3 methods: fieldID uses the known system geometry to figure out the right region
                   calib is for calibration when the system is not known, calculates from the image characteristics
                   direct sets teh values manually (backup method)


        """
        image = self.actor.image
        db = self.connectToDB(cmd)

        cmd.inform(f'text="loading telescope parameters for frame={frameId} with fibreMode={self.fibreMode}'
                   ' at z={zenithAngle) rot={insRot}"')

        # zenithAngle,insRot=dbTools.loadTelescopeParametersFromDB(db,int(frameId))
        #cmd.diag(f'text="zenithAngle={zenithAngle}, insRot={insRot}"')

        # different transforms for different setups: with and w/o field elements
        if(self.fibreMode == 'full'):
            centrePosPix = mcsToolsNew.transformToPix(
                self.centrePos, self.rotCent, self.offset, zenithAngle, insRot, fieldElement=True, pixScale=0)
        elif(self.fibreMode == 'comm'):
            centrePosPix = mcsToolsNew.transformToPix(
                self.centrePos, self.rotCent, self.offset, zenithAngle, insRot, fieldElement=False, pixScale=0)
        elif(self.fibreMode == 'asrd'):
            centrePosPix = self.centrePos

        np.save("cpos.npy", centrePosPix)
        self.findThresh, self.centThresh = mcsToolsNew.getThresh(
            image, centrePosPix, 'full', self.centParms['threshSigma'], self.centParms['findSigma'], self.centParms['centSigma'])

        a1 = self.centParms['threshSigma']
        a2 = self.centParms['findSigma']
        a3 = self.centParms['centSigma']
        cmd.inform(f'text="findThresh ={self.findThresh:.2}, centThresh = {self.centThresh:.2}"')
        cmd.inform(f'text="{a1} {a2} {a3}"')

    def setCentroidParams(self, cmd):
        """
        top level routine for setting centroid parameters. REads the defaults from teh config fil,e
        then changes any specified in the keywords argument. 

        """

        self.centParms = mcsToolsNew.getCentroidParams(cmd)

    def runCentroidSEP(self, cmd):
        cmdKeys = cmd.cmd.keywords
        self.newTable = "newTable" in cmdKeys

        cmd.debug('text="newTable value = %s"' % (self.newTable))

        image = self.actor.image

        cmd.inform(f'state="measuring cached image: {image.shape}"')
        centroids = sep.extract(image.astype(float), 300)
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

        cmd.inform(f'state="measuring cached image: {image.shape}"')
        a = centroid.centroid_only(image.astype('<i4'),
                                   centParms['fwhmx'], centParms['fwhmy'], self.findThresh, self.centThresh,
                                   centParms['boxFind'], centParms['boxCent'],
                                   centParms['nmin'], centParms['nmax'], centParms['maxIt'], 0)

        centroids = np.frombuffer(a, dtype='<f8')
        centroids = np.reshape(centroids, (len(centroids)//7, 7))
        nSpots = centroids.shape[0]
        points = np.empty((nSpots, 8))
        points[:, 0] = np.arange(nSpots)
        points[:, 1:] = centroids[:, 0:]

        self.centroids = points

        self.nCentroid = len(points)
        cmd.inform('text="%d centroids"' % (len(centroids)))
        cmd.inform('state="centroids measured"')

    def _writeTelescopeInfo(self, cmd, telescopeInfo, conn=None):

        # Let the database handle the primary key
        db = self.connectToDB(cmd)
        res = db.session.execute('select * FROM "mcs_exposure" where false')
        colnames = res.keys()
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
        taken_at = telescopeInfo['starttime']
        taken_in_hst_at = telescopeInfo['starttime']
        line = '%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s,%s' % (telescopeInfo['frameid'],
                                                                       telescopeInfo['visitid'],
                                                                       telescopeInfo['exptime']/1000.0,
                                                                       telescopeInfo['altitude'],
                                                                       telescopeInfo['azimuth'],
                                                                       telescopeInfo['instrot'],
                                                                       adc_pa, dome_temperature,
                                                                       dome_pressure, dome_humidity,
                                                                       outside_temperature, outside_pressure,
                                                                       outside_humidity,
                                                                       mcs_cover_temperature,
                                                                       mcs_m1_temperature,
                                                                       taken_at, taken_in_hst_at)

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
        db = self.connectToDB(None)
        session = db.session
        with session.connection().connection.cursor() as cursor:
            cursor.copy_expert(sql, dataBuf)
            cursor.close()
        session.execute('commit')

    def _readData(self, sql):
        """Wrap a direct COPY_TO via sqlalchemy. """

        dataBuf = io.StringIO()
        db = self.connectToDB(None)
        session = db.session
        with session.connection().connection.cursor() as cursor:
            cursor.copy_expert(sql, dataBuf)
        dataBuf.seek(0, 0)
        return dataBuf

    def _writeCentroids(self, centArr, nextRowId, frameId, moveId, conn=None):
        """ Write all measurements for a given (frameId, moveId) """

        now = datetime.datetime.now()
        now.strftime("%Y-%m-%d %H:%M:%S")

        # Save measurements to a CSV buffer
        measBuf = io.StringIO()

        data = np.insert(centArr, 5, 0, axis=1)

        np.savetxt(measBuf, data[:, 1:8], delimiter=',', fmt='%0.6g')
        measBuf.seek(0, 0)

        # Let the database handle the primary key
        db = self.connectToDB(None)
        colnames = db.session.execute('select * FROM "mcs_data" where false')
        realcolnames = colnames.keys()[0:]

        colname = []
        for i in realcolnames:
            x = '"'+i+'"'
            colname.append(x)
        buf = io.StringIO()
        for l_i in range(len(centArr)):
            # line = '%s,%d,%d,%d,%s' % (now.strftime("%Y-%m-%d %H:%M:%S"),
            #                           frameId, moveId, l_i+1, measBuf.readline())
            line = '%d,%d,%s' % (frameId, l_i+1, measBuf.readline())
            buf.write(line)
        buf.seek(0, 0)

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
        realcolnames = colnames.keys()[0:]

        buf = io.StringIO()

        for i in range(len(pfi_x)):
            line = '%d,%d,%d,%d, %f,%f,%f,%f,%d,%d,%d,%f,%d,%d,%d,%f,%d\n'%(frameId/100, 1, i+1, 99, pfi_x[i], pfi_y[i],
                                                                            pfi_x[i], pfi_y[i], 1, 1, 99, 99, 1, 1, 99, 99, 0)
            buf.write(line)
        buf.seek(0, 0)

        self._writeData('cobra_target', realcolnames, buf)
        buf.seek(0, 0)
