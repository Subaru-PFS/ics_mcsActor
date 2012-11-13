#!/usr/bin/env python

import time

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr

class McsCmd(object):

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor

        # Declare the commands we implement. When the actor is started
        # these are registered with the parser, which will call the
        # associated methods when matched. The callbacks will be
        # passed a single argument, the parsed and typed command.
        #
        self.vocab = [
            ('ping', '', self.ping),
            ('status', '', self.status),
            ('expose', '<expTime>', self.expose),
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time"),
                                        )


    def ping(self, cmd):
        """Query the actor for liveness/happiness."""

        cmd.finish("text='Present and (probably) well'")

    def status(self, cmd):
        """Report status and version; obtain and send current data"""

        self.actor.sendVersionKey(cmd)

        keyStrings = ['text="nothing to say, really"']
        keyMsg = '; '.join(keyStrings)

        cmd.inform(keyMsg)
        cmd.diag('text="still nothing to say"')
        cmd.finish()

    def expose(self, cmd):
        """ Do something pointless. """

        expTime = cmd.cmd.keywords["expTime"].values[0]

        cmd.inform('state="exposing"')
        if expTime > 0:
            time.sleep(expTime + 0.5)

        cmd.inform('state="measuring"')
        centroids = ["%4.4f,%4.4f" % (0.5*x,0.5*x+0.25) for x in range(2400)]
        cmd.inform('state="measured"')
        cmd.inform("centroids=%s" % ", ".join(centroids))
       
        cmd.finish('state="done"')

