#!/usr/bin/env python

import actorcore.Actor
import fakeCamera

class Mcs(actorcore.Actor.Actor):
    def __init__(self, name, productName=None, configFile=None, debugLevel=30):
        # This sets up the connections to/from the hub, the logger, and the twisted reactor.
        #
        actorcore.Actor.Actor.__init__(self, name, 
                                       productName=productName, 
                                       configFile=configFile)


    def connectCamera(self):
        
        self.camera = 
#
# To work

def main():
    mcs = Mcs('mcs', productName='mcsActor')
    mcs.run()

if __name__ == '__main__':
    main()

    
    
