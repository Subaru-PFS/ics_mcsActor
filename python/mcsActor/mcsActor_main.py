#!/usr/bin/env python

import opscore.utility.sdss3logging

import opscore.actor.model
import actorcore.Actor

class Mcs(actorcore.Actor.Actor):
    def __init__(self, name, productName=None, configFile=None, debugLevel=30):
        # This sets up the connections to/from the hub, the logger, and the twisted reactor.
        #
        actorcore.Actor.Actor.__init__(self, name, productName=productName, configFile=configFile)
        #self.logger.setLevel(debugLevel)
        
        # Explicitly declare which actors we need to know about, so we can access their keywords.
        #
        self.models = {}
        for actor in []:
            self.models[actor] = opscore.actor.model.Model(actor)

#
# To work

if __name__ == '__main__':
    mcs = Mcs('mcs', 'mcsActor')
    mcs.run()
    
