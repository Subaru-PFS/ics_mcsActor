#!/usr/bin/env python

from __future__ import absolute_import
import actorcore.ICC
import aitCamera
import mcsCamera
import fakeCamera
from past.builtins import reload

class Mcs(actorcore.ICC.ICC):
    def __init__(self, name, productName=None, configFile=None, debugLevel=30):
        # This sets up the connections to/from the hub, the logger, and the twisted reactor.
        #
        super().__init__(name, 
                         productName=productName,
                         modelNames=('gen2','meb','mcs'),
                         configFile=configFile)

        # We will actually use a allocator with "global" sequencing
        self.exposureID = 0
        
        self.connectCamera(self.bcast)
        
    def connectCamera(self, cmd, doFinish=True):

        reload(mcsCamera)
        try:
            self.camera = mcsCamera.mcsCamera()
            self.camera.initialCamera(cmd)
            cmd.inform('text="loaded real MCS camera"')
        except Exception as e:
            cmd.warn('text="failed to load MCS camera, loading fakeCamera: %s"' % (e))
            self.camera = fakeCamera.FakeCamera()
            self.camera.initialCamera(cmd)

        self.camera.sendStatusKeys(cmd)
#
# To work

def main():
    mcs = Mcs('mcs', productName='mcsActor')
    mcs.run()

if __name__ == '__main__':
    main()

    
    
