#!/usr/bin/env python

import base64
import numpy
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
            ('centroid', '<expTime>', self.centroid),
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

    def _doExpose(self, cmd):
        """ Take an exposure. """
        
        expTime = cmd.cmd.keywords["expTime"].values[0]

        cmd.inform('state="exposing"')
        if expTime > 0:
            time.sleep(expTime + 0.5)
        image = None
        
        filename = "fakeFilename.fits"
        cmd.inform("filename=%s" % (filename))

        return filename, image
            
    def expose(self, cmd):
        """ Take an exposure. Does not centroid. """

        filename, image = self._doExpose(cmd)
        cmd.finish('state="done"')

    def _encodeArray(self, array):
        """ Temporarily wrap a binary array transfer encoding. """

        # Actually, we want dtype,naxis,axNlen,base64(array)
        return base64.b64encode(array.tostring())
        
    def centroid(self, cmd):
        """ Take an exposure and measure centroids. """

        filename, image = self._doExpose(cmd)

        # The encoding scheme is temporary, and will be replaced with 
        cmd.inform('state="measuring"')
        centroids = numpy.random.random(4800).astype('f4').reshape(2400,2)

        centroidsStr = self._encodeArray(centroids)
        cmd.inform('state="measured"; centroids=%s' % (centroidsStr))

        cmd.finish('state="done"')
       

