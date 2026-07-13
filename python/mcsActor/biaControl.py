from builtins import object


class BiaControl(object):
    """Drive the sps bia illuminator around an MCS exposure.

    The bia is only actually driven when this has been activated with the mcs biaControl command *and* the
    canon camera is loaded, which is how we know that we are on the telescope and not in a lab.
    """
    # camera which is only ever used on the telescope.
    telescopeCamera = 'canon_50m'
    # exposure types for which the cobras need to be illuminated.
    biaExposureTypes = ('object', 'test')

    switchOnTimeLim = 30
    switchOffTimeLim = 30

    def __init__(self, actor):
        self.actor = actor
        self.activated = False

    @property
    def onTelescope(self):
        return self.actor.cameraName == BiaControl.telescopeCamera

    def declare(self, cmd, activated):
        """Activate or deactivate the bia control."""
        self.activated = activated and self.onTelescope
        self.genKeys(cmd)

    def genKeys(self, cmd):
        """Generate biaControl keyword."""
        cmd.inform('biaControl=%s' % ('on' if self.activated else 'off'))

    def doDriveBia(self, expType):
        """Whether the bia should be driven for that exposure."""
        return self.activated and expType in BiaControl.biaExposureTypes

    def switchOn(self, cmd, expType):
        """Switch bia on, note that sps only returns when the bia is actually emitting light.

        :raise: RuntimeError if sps failed to turn the bia on, the exposure is not worth taking in that case.
        """
        if not self.doDriveBia(expType):
            return

        cmdVar = self.actor.cmdr.call(actor='sps', cmdStr='bia on', forUserCmd=cmd,
                                      timeLim=BiaControl.switchOnTimeLim)
        if cmdVar.didFail:
            raise RuntimeError('failed to turn the bia on, cobras would not be illuminated...')

    def switchOff(self, cmd, expType):
        """Switch bia off, called as soon as the camera exposure is done."""
        if not self.doDriveBia(expType):
            return

        cmdVar = self.actor.cmdr.call(actor='sps', cmdStr='bia off', forUserCmd=cmd,
                                      timeLim=BiaControl.switchOffTimeLim)
        if cmdVar.didFail:
            cmd.warn('text="failed to turn the bia off !"')
