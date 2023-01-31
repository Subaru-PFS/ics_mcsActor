#!/usr/bin/env python

from importlib import reload
import os
import pathlib
import socket

from opscore.utility import sdss3logging
import actorcore.ICC

class Mcs(actorcore.ICC.ICC):
    def __init__(self, name, productName=None, configFile=None, debugLevel=30):
        # This sets up the connections to/from the hub, the logger, and the twisted reactor.
        #
        super().__init__(name,
                         productName=productName,
                         modelNames=('gen2', 'meb', 'mcs', 'peb', 'sps'),
                         configFile=configFile)

        self.connectCamera(self.bcast)

        self.cameraName = None

    def connectCamera(self, cmd, camera=None, doFinish=True):
        # For the mcsActor at Subaru, we have two possible cameras:
        #  - on the pfic host from ASRD, use the rmod. Else use the canon.
        if camera is None:
            if socket.gethostname() == 'us-mcs':
                camera = 'rmod_71m'
            else:
                camera = 'canon_50m'

        configDir = pathlib.Path(os.environ['ICS_MCSACTOR_DIR'], 'etc')

        try:
            if camera == 'rmod_71m':
                import rmodCamera
                reload(rmodCamera)

                configPath = pathlib.Path(configDir, 'illunis-71mp.cfg')
                self.camera = rmodCamera.rmodCamera()
                self.camera.initialCamera(cmd, configPath=configPath)
                cmd.inform('text="loaded RMOD-71M camera"')
            if camera == 'canon_50m':
                import mcsCamera
                reload(mcsCamera)

                configPath = pathlib.Path(configDir, 'canon-8960x5778.cfg')
                self.camera = mcsCamera.mcsCamera()
                self.camera.initialCamera(cmd, configPath=configPath)
                cmd.inform('text="loaded real MCS camera"')
        except Exception as e:
            cmd.warn('text="failed to load MCS camera, loading fakeCamera: %s"' % (e))

            import fakeCamera

            self.camera = fakeCamera.FakeCamera()
            self.camera.initialCamera(cmd)

        self.camera.sendStatusKeys(cmd)
        self.cameraName = camera
        cmd.inform(f'text="camera name = {camera}"')
# To work


def main():
    mcs = Mcs('mcs', productName='mcsActor')
    mcs.run()


if __name__ == '__main__':
    main()
