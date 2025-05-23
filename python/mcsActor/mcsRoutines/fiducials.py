import pandas as pd


class Fiducials(pd.DataFrame):
    """
    A class to represent fiducials data as a pandas DataFrame.

    Attributes
    ----------
    BROKEN : int
        Flag value indicating a broken fiducial.
    BAD : int
        Flag value indicating a bad fiducial.
    OUTER_RING : int
        Flag value indicating a fiducial in the outer ring.

    Methods
    -------
    read(butler)
        Retrieve fiducials data using a Butler instance.
    goodMask
        A property to return a boolean mask for good fiducials (not broken or bad).
    outerRingMask
        A property to return a boolean mask for fiducials in the outer ring.
    """
    BROKEN = 1
    BAD = 2
    OUTER_RING = 4
    SIGMACLIP = 8

    @classmethod
    def read(cls, butler):
        """
        Retrieve fiducials data using a Butler instance.

        Parameters
        ----------
        butler : pfs.utils.butler.Butler
            A Butler instance to retrieve fiducials data.

        Returns
        -------
        Fiducials
            A Fiducials instance containing the fiducials data.
        """
        # Retrieve fiducials data using butler
        return cls(butler.get('fiducials'))

    @property
    def goodMask(self):
        """
        A boolean mask for good fiducials (not broken or bad).

        Returns
        -------
        pandas.Series
            A boolean pandas Series indicating good fiducials.
        """
        return (self['flag'] & (self.BROKEN | self.BAD)) == 0

    @property
    def outerRingMask(self):
        """
        A boolean mask for fiducials in the outer ring.

        Returns
        -------
        pandas.Series
            A boolean pandas Series indicating fiducials in the outer ring.
        """
        return (self['flag'] & self.OUTER_RING) != 0
