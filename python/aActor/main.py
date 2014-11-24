#!/usr/bin/env python

from actorcore.Actor import Actor
import fakeCamera

class AActor(Actor):
    def __init__(self, name, productName=None, configFile=None, debugLevel=30):
        # This sets up the connections to/from the hub, the logger, and the twisted reactor.
        #
        Actor.__init__(self, name, 
                       productName=productName, 
                       configFile=configFile)

        # We will actually use a allocator with "global" sequencing
        self.exposureID = 0
        
        self.connectCamera(self.bcast)
        
    def connectCamera(self, cmd, doFinish=True):
        reload(fakeCamera)
        self.camera = fakeCamera.FakeCamera()
        self.camera.sendStatusKeys(cmd)

#
# To work

def main():
    actor = AActor('a', productName='aActor')
    actor.run()

if __name__ == '__main__':
    main()

    
    
