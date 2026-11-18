"""What the actor records of a frame: detector, fiducial LEDs and bia in the header, software releases, the
centroid parameters it measures with."""
import logging
import pathlib
import queue
import tempfile
import types
import unittest
import unittest.mock

import numpy as np

import mcsActor.Commands.McsCmd as McsCmd
import test_findSpots

CONFIG = dict(centroidParams=dict(fwhmx=4, fwhmy=9, boxFind=14, boxCent=14, nmin=10, maxIt=20, centSigma=25,
                                  findSigma=25, aperture=10, innerAnn=20, outerAnn=24),
              detector=dict(name='MCS cam#1', gain=2.24))


class FakeCmd(test_findSpots.FakeCmd):
    """Keeps the replies of a command, diagnostics included."""

    def diag(self, text):
        self.replies.append(('diag', text))


def keyVar(*values):
    """Stand-in for a keyVar holding values."""
    return types.SimpleNamespace(valueList=list(values))


def commandSet(**models):
    """An McsCmd without its constructor, with the config above and the models given."""
    owner = McsCmd.McsCmd.__new__(McsCmd.McsCmd)
    owner.actor = types.SimpleNamespace(actorConfig=CONFIG, models=models)
    owner.centOverrides = dict()
    owner.centParms = McsCmd.centroidParameters(CONFIG['centroidParams'], owner.centOverrides)
    owner.detector = McsCmd.configuredDetector(CONFIG)
    owner.logger = logging.getLogger('test')
    return owner


def finished(cmd):
    return [text for level, text in cmd.replies if level == 'finish']


class DetectorTestCase(unittest.TestCase):

    def testConfiguredDetector(self):
        """Name and gain come from the detector section; without a name or a numeric gain the config is refused."""
        self.assertEqual(McsCmd.configuredDetector(CONFIG), dict(name='MCS cam#1', gain=2.24))
        for config in ({}, dict(detector=None), dict(detector=dict(name='MCS cam#1')),
                       dict(detector=dict(gain=2.24)), dict(detector=dict(name='MCS cam#1', gain='high'))):
            with self.assertRaises(ValueError, msg=config):
                McsCmd.configuredDetector(config)

    def testHeaderDetector(self):
        """DETECTOR and GAIN are those of the config, and the cards of the exposure follow the others."""
        owner = commandSet()
        exposureCards = [('OBS-MOD', 'SPEC_IMAG', 'mode'), ('W_BIAPOW', 40, 'bia power'),
                         ('W_FIDDTY', 6.9, '[%] fiducial LED duty cycle')]
        hdr = McsCmd.McsCmd._constructHeader(owner, FakeCmd(), 'file', 'object', 4800.0, 1.7e9, 15503500,
                                             exposureCards)

        self.assertEqual((hdr['DETECTOR'], hdr['GAIN']), ('MCS cam#1', 2.24))
        self.assertEqual((hdr['W_BIAPOW'], hdr['W_FIDDTY']), (40, 6.9))
        self.assertLess(list(hdr.keys()).index('GAIN'), list(hdr.keys()).index('W_BIAPOW'))


class LedCardsTestCase(unittest.TestCase):

    def testLedCards(self):
        """The duty cycle and period the fiducial LEDs output, in % and us."""
        owner = commandSet(peb=types.SimpleNamespace(keyVarDict=dict(dutycycle=keyVar(6.9, 1.0, 10.3),
                                                                      ledperiod=keyVar(100, 100, 100000))))
        cmd = FakeCmd()
        cards = McsCmd.McsCmd._ledCards(owner, cmd)

        self.assertEqual([card[:2] for card in cards], [('W_FIDDTY', 6.9), ('W_FIDPER', 100.0)])
        self.assertEqual(cmd.replies, [])

    def testUnknownLeds(self):
        """Values peb has not reported, or a missing peb model, give -9999; only the missing model warns."""
        owner = commandSet(peb=types.SimpleNamespace(keyVarDict=dict(dutycycle=keyVar(None, None, None),
                                                                      ledperiod=keyVar(None, None, None))))
        cmd = FakeCmd()
        self.assertEqual([card[1] for card in McsCmd.McsCmd._ledCards(owner, cmd)], [-9999.0, -9999.0])
        self.assertEqual(cmd.replies, [])

        owner, cmd = commandSet(), FakeCmd()
        self.assertEqual([card[1] for card in McsCmd.McsCmd._ledCards(owner, cmd)], [-9999.0, -9999.0])
        self.assertEqual([level for level, _ in cmd.replies], ['warn', 'warn'])


class ExposureCardsTestCase(unittest.TestCase):

    def testCardsGatheredWhileBiaIsOn(self):
        """The cards of the other actors are gathered after the camera has exposed and before the bia goes off."""
        events = []
        owner = commandSet()
        owner.actor.biaControl = types.SimpleNamespace(switchOn=lambda cmd, expType: events.append('bia on'),
                                                       switchOff=lambda cmd, expType: events.append('bia off'))
        frame = np.zeros((4, 6), dtype='<u2')
        owner.actor.camera = types.SimpleNamespace(expose=lambda *args, **kwargs: events.append('expose') or frame)
        owner._exposureCards = lambda cmd: events.append('cards') or ['cards']
        fileIds = queue.Queue()
        fileIds.put(dict(frameId=15503500))
        owner.requestNextFileIds = lambda cmd, frameId: fileIds
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        owner.butler = types.SimpleNamespace(getPath=lambda *args, **kwargs: pathlib.Path(directory.name, 'f.fits'))
        headers = []
        owner._constructHeader = lambda *args: headers.append(args) or 'header'

        _, hdr, image = McsCmd.McsCmd._doExpose(owner, FakeCmd(), 4800.0, 'object', 15503500)

        self.assertEqual(events, ['bia on', 'expose', 'cards', 'bia off'])
        self.assertEqual((hdr, headers[0][-1]), ('header', ['cards']))
        self.assertIs(image, frame)

    def testBiaOffWhenTheCameraFails(self):
        """A failed exposure still turns the bia off."""
        events = []
        owner = commandSet()
        owner.actor.biaControl = types.SimpleNamespace(switchOn=lambda cmd, expType: events.append('bia on'),
                                                       switchOff=lambda cmd, expType: events.append('bia off'))

        def expose(*args, **kwargs):
            raise RuntimeError('camera failed')

        owner.actor.camera = types.SimpleNamespace(expose=expose)
        owner.requestNextFileIds = lambda cmd, frameId: queue.Queue()
        with self.assertRaises(RuntimeError):
            McsCmd.McsCmd._doExpose(owner, FakeCmd(), 4800.0, 'object', 15503500)
        self.assertEqual(events, ['bia on', 'bia off'])


class ReleaseTestCase(unittest.TestCase):

    def testUnreleasedVersions(self):
        """Tags as git describe and EUPS print them are releases; anything past, beside or without a tag is not."""
        released = ('1.10.10', '7.4.46', 'w.2026.46', '1.8.119')
        unreleased = ('1.10.10-1-g8ba5dec', '1.10.10-dirty', '7.4.46-3-gabc1234-dirty', '8ba5dec', 'unknown',
                      'tickets/INSTRM-3109', 'LOCAL:/home/pfs/ics_mcsActor', '', None)
        versions = {f'r{k}': version for k, version in enumerate(released)}
        versions.update({f'u{k}': version for k, version in enumerate(unreleased)})

        self.assertEqual(McsCmd.unreleasedVersions(versions),
                         {f'u{k}': version for k, version in enumerate(unreleased)})

    def testWarning(self):
        """A frame warns naming the products that are not releases, and is silent when all are."""
        owner, cmd = commandSet(), FakeCmd()
        owner.softwareVersions = dict(version_actor='1.10.10-1-g8ba5dec', version_instdata='1.8.119',
                                      version_pfs_utils=None)
        McsCmd.McsCmd.warnUnreleased(owner, cmd)

        self.assertEqual(len(cmd.replies), 1)
        level, text = cmd.replies[0]
        self.assertEqual(level, 'warn')
        self.assertIn('version_actor=1.10.10-1-g8ba5dec', text)
        self.assertIn('version_pfs_utils=None', text)
        self.assertNotIn('instdata', text)

        owner.softwareVersions = dict(version_actor='1.10.11', version_instdata='1.8.119', version_pfs_utils='7.4.46')
        cmd = FakeCmd()
        McsCmd.McsCmd.warnUnreleased(owner, cmd)
        self.assertEqual(cmd.replies, [])

    def testUnreadableVersion(self):
        """A version that cannot be read is None."""
        with unittest.mock.patch.object(McsCmd, 'getVersion', side_effect=RuntimeError('no version')):
            self.assertIsNone(McsCmd.softwareVersion('pfs_instdata'))
        with unittest.mock.patch.object(McsCmd, 'getVersion', return_value='1.8.119'):
            self.assertEqual(McsCmd.softwareVersion('pfs_instdata'), '1.8.119')


class CentroidParamsTestCase(unittest.TestCase):

    def exposure(self, owner):
        """The parameters an exposure measures with."""
        return McsCmd.centroidParameters(owner.actor.actorConfig['centroidParams'], owner.centOverrides)

    def testOverridesLast(self):
        """Overrides finish the command and hold for the exposures after; the other parameters stay configured."""
        owner, cmd = commandSet(), FakeCmd(fwhmx=12.0, boxCent=20)
        McsCmd.McsCmd.setCentroidParams(owner, cmd)

        self.assertEqual(len(finished(cmd)), 1)
        self.assertIn('fwhmx=12.0', finished(cmd)[0])
        for _ in range(2):
            self.assertEqual(self.exposure(owner), dict(CONFIG['centroidParams'], fwhmx=12.0, boxCent=20))
        self.assertEqual(owner.centParms, self.exposure(owner))
        self.assertEqual(CONFIG['centroidParams']['fwhmx'], 4)

    def testApertureOverrides(self):
        """setApertureParams overrides the aperture and the annulus radii, keeping earlier overrides."""
        owner = commandSet()
        McsCmd.McsCmd.setCentroidParams(owner, FakeCmd(findSigma=30))
        cmd = FakeCmd(aperture=8.0, innerRad=18.0, outerRad=22.0)
        McsCmd.McsCmd.setApertureParams(owner, cmd)

        self.assertEqual(len(finished(cmd)), 1)
        self.assertEqual(self.exposure(owner), dict(CONFIG['centroidParams'], findSigma=30, aperture=8.0,
                                                    innerAnn=18.0, outerAnn=22.0))

    def testReset(self):
        """reset drops every override."""
        owner = commandSet()
        McsCmd.McsCmd.setCentroidParams(owner, FakeCmd(fwhmx=12.0))
        McsCmd.McsCmd.setApertureParams(owner, FakeCmd(aperture=8.0))
        cmd = FakeCmd(reset=None)
        McsCmd.McsCmd.resetCentroidParams(owner, cmd)

        self.assertEqual(self.exposure(owner), CONFIG['centroidParams'])
        self.assertEqual(owner.centParms, CONFIG['centroidParams'])
        self.assertIn('none', finished(cmd)[0])

    def testOtherKeywordsIgnored(self):
        """Keywords that are not centroid parameters override nothing."""
        owner = commandSet()
        McsCmd.McsCmd.setCentroidParams(owner, FakeCmd(expTime=4.8, threshFact=3))

        self.assertEqual(self.exposure(owner), CONFIG['centroidParams'])


if __name__ == '__main__':
    unittest.main()
