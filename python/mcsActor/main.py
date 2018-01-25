#!/usr/bin/env python

import actorcore.Actor
import aitCamera
import fakeCamera
from past.builtins import reload

class Mcs(actorcore.Actor.Actor):
    def __init__(self, name, productName=None, configFile=None, debugLevel=30):
        # This sets up the connections to/from the hub, the logger, and the twisted reactor.
        #
        actorcore.Actor.Actor.__init__(self, name, 
                                       productName=productName, 
                                       configFile=configFile)

        # We will actually use a allocator with "global" sequencing
        self.exposureID = 0
        
        self.connectCamera(self.bcast)
        
    def connectCamera(self, cmd, doFinish=True):
        reload(aitCamera)
        self.camera = aitCamera.aitCamera()
        #self.camera = fakeCamera.FakeCamera()
        self.camera.sendStatusKeys(cmd)
        self.camera.initialCamera(cmd)
#
# To work

def main():
    mcs = Mcs('mcs', productName='mcsActor')
    mcs.run()

if __name__ == '__main__':
    main()

    
    
