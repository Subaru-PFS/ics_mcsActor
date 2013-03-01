import numpy
import time
import camera

class FakeCamera(camera.Camera):
    def __init__(self):
        # Should read from config file....
        self.imageSize = (4096, 4096)
        self.biasLevel = 100
        self.readNoise = 5.0

    def _readoutTime(self):
        return 0.5

    def _exposureOverheadTime(self):
        return 0.1

    def expose(self, cmd, expTime, expType):
        """ Generate an 'exposure' image. We don't have an actual camera, so generate some 
        plausible image. 

        Args:
           cmd     - a Command object to report to. Ignored if None.
           expTime - the fake exposure time. 
           expType - ("bias", "dark", "object", "flat", "test")
           
        Returns:
           - the image.

        Keys:
           exposureState
        """

        if not expType:
            expType = 'bias'
        if cmd:
            cmd.inform('exposureState="exposing"')
        if expType not in ('bias', 'dark') and expTime >= 0:
            time.sleep(expTime + self._exposureOverheadTime)

        if cmd:
            cmd.inform('exposureState="reading"')
        image = numpy.random.normal(self.biasLevel, 
                                    scale=self.readNoise, 
                                    size=self.imageSize).astype('u2')

        return image
        
        
