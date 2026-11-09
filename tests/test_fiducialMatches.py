"""Fiducial matches recorded over the transform passes, and written to fiducial_fiber_match."""
import unittest

import numpy as np
import pandas as pd

import mcsActor.mcsRoutines.dbRoutinesMCS as dbTools
import mcsActor.mcsRoutines.fiducials as fiducials
import mcsActor.mcsRoutines.mcsRoutines as mcsRoutines


class FakeDB:
    """Keeps what is inserted, per table."""

    def __init__(self):
        self.inserted = {}

    def insert_dataframe(self, table, df):
        self.inserted[table] = df


def passes(*ffids, nFids=4):
    """fidMask and fidSpot after recording each pass's ffid with bits 1, 2, 4..."""
    fidMask = np.zeros(nFids, dtype=int)
    fidSpot = np.full(nFids, -1, dtype=int)
    for i, ffid in enumerate(ffids):
        mcsRoutines.recordFiducialMatches(np.array(ffid), 2 ** i, fidMask, fidSpot)
    return fidMask, fidSpot


class FakeButler:
    """Hands out one fiducial table."""

    def __init__(self, table):
        self.table = table

    def get(self, name):
        return self.table.copy()


class ReadFiducialsTestCase(unittest.TestCase):

    def testRowsFollowFiducialId(self):
        """A table stored out of order is read in fiducialId order."""
        table = pd.DataFrame(dict(fiducialId=[3, 1, 2], x_mm=[30., 10., 20.], flag=[0, 4, 1]))
        fids = fiducials.Fiducials.read(FakeButler(table))

        self.assertEqual(fids.fiducialId.tolist(), [1, 2, 3])
        self.assertEqual(fids.x_mm.tolist(), [10., 20., 30.])
        self.assertEqual(fids.index.tolist(), [0, 1, 2])
        self.assertEqual(fids.outerRingMask.tolist(), [True, False, False])


class RecordFiducialMatchesTestCase(unittest.TestCase):

    def testMatchedEveryPass(self):
        """A fiducial matched by every pass keeps its spot and every bit."""
        fidMask, fidSpot = passes([-1, 2, -1], [-1, 2, -1], [-1, 2, -1])
        self.assertEqual(fidMask.tolist(), [0, 7, 0, 0])
        self.assertEqual(fidSpot.tolist(), [-1, 1, -1, -1])

    def testClippedKeepsEarlierSpot(self):
        """A fiducial the last pass clips keeps the spot of the pass that matched it."""
        fidMask, fidSpot = passes([-1, 2, 3], [-1, 2, 3], [-1, 2, -1])
        self.assertEqual(fidMask.tolist(), [0, 7, 3, 0])
        self.assertEqual(fidSpot.tolist(), [-1, 1, 2, -1])

    def testLaterPassMovesSpot(self):
        """A later pass matching a fiducial to another spot replaces the earlier one."""
        fidMask, fidSpot = passes([1, -1, -1], [-1, 1, -1])
        self.assertEqual(fidMask.tolist(), [3, 0, 0, 0])
        self.assertEqual(fidSpot.tolist(), [1, -1, -1, -1])

    def testSpotGivenToAnotherFiducial(self):
        """A spot a later pass gives to another fiducial is taken from the earlier one."""
        fidMask, fidSpot = passes([3, -1], [4, -1])
        self.assertEqual(fidMask.tolist(), [0, 0, 1, 2])
        self.assertEqual(fidSpot.tolist(), [-1, -1, -1, 0])

    def testNoMatch(self):
        """A pass matching nothing changes nothing."""
        fidMask, fidSpot = passes([2, -1], [-1, -1])
        self.assertEqual(fidMask.tolist(), [0, 1, 0, 0])
        self.assertEqual(fidSpot.tolist(), [-1, 0, -1, -1])


class WriteFidToDBTestCase(unittest.TestCase):

    def testEveryMatchedFiducialHasItsSpot(self):
        """Clipped fiducials are written with their spot; unmatched ones with -1 and NaN."""
        fids = pd.DataFrame(dict(fiducialId=[1, 2, 3], match_mask=[7, 3, 0],
                                 fiducial_tweaked_x_mm=[1., 2., 3.], fiducial_tweaked_y_mm=[4., 5., 6.]))
        mcsData = pd.DataFrame(dict(spot_id=[10, 11, 12], pfi_center_x_mm=[.1, .2, .3],
                                    pfi_center_y_mm=[.4, .5, .6]))
        db = FakeDB()
        dbTools.writeFidToDB(db, np.array([2, 0, -1]), mcsData, 12345601, fids)

        df = db.inserted['fiducial_fiber_match']
        self.assertEqual(df.spot_id.tolist(), [12, 10, -1])
        self.assertEqual(df.match_mask.tolist(), [7, 3, 0])
        np.testing.assert_array_equal(df.pfi_center_x_mm, [.3, .1, np.nan])
        np.testing.assert_array_equal(df.pfi_center_y_mm, [.6, .4, np.nan])
        self.assertEqual(df.pfs_visit_id.tolist(), [123456] * 3)
        self.assertEqual(df.iteration.tolist(), [1] * 3)


if __name__ == '__main__':
    unittest.main()
