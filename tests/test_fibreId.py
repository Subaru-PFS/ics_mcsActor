"""Spot-to-cobra matching on small synthetic cobra layouts."""
import unittest

import numpy as np

import mcsActor.mcsRoutines.mcsRoutines as mcsRoutines

ARM = 4.7  # mm, L1+L2 of every synthetic cobra


def layout(centres):
    """Geometry arguments of fibreId for cobras at `centres` (mm), all with arm length ARM."""
    nCobra = len(centres)
    centrePos = np.column_stack([np.arange(nCobra), centres])
    armLength = np.full(nCobra, ARM)
    adjacent = mcsRoutines.makeAdjacentList(centrePos, armLength)
    fids = dict(fiducialId=np.array([1]), x_mm=np.array([500.]), y_mm=np.array([500.]))
    return centrePos, armLength, fids, np.zeros((nCobra, 4)), np.arange(nCobra), adjacent


def match(centres, targets, spots):
    """cobra index -> (spot_id, flags) from fibreId; spot k has spot_id 100+k."""
    centrePos, armLength, fids, dotPos, goodIdx, adjacent = layout(np.asarray(centres, dtype=float))
    points = np.column_stack([100 + np.arange(len(spots)), np.asarray(spots, dtype=float)])
    tarPos = np.column_stack([np.arange(len(targets)), np.asarray(targets, dtype=float)])
    cobraMatch, _, nAmbiguous = mcsRoutines.fibreId(points, centrePos, armLength, tarPos, fids, dotPos,
                                                    goodIdx, adjacent, 'target')
    assert nAmbiguous == 0
    return {i: (int(row[1]), int(row[4])) for i, row in enumerate(cobraMatch)}


def unassigned(nPoints, nCobras, potCobraMatch, potPointMatch):
    """lastPassDist's list arguments with every point and cobra unassigned."""
    return [], list(range(nCobras)), [], list(range(nPoints)), potCobraMatch, potPointMatch


class FibreIdTestCase(unittest.TestCase):

    def testSpotsAtTargets(self):
        """Every cobra finds the spot on its target."""
        rng = np.random.default_rng(0)
        centres = [(8. * i + 4. * (j % 2), 6.93 * j) for j in range(3) for i in range(4)]
        radius, angle = rng.uniform(0.5, ARM, len(centres)), rng.uniform(0, 2 * np.pi, len(centres))
        targets = np.array(centres) + np.column_stack([radius * np.cos(angle), radius * np.sin(angle)])
        spots = targets + rng.normal(0, 0.005, targets.shape)

        self.assertEqual({i: s for i, (s, _) in match(centres, targets, spots).items()},
                         {i: 100 + i for i in range(len(centres))})

    def testOvershootKeptByItsCobra(self):
        """A spot 0.15 mm beyond its cobra's arm length still goes to that cobra, not the neighbour."""
        centres = [(0., 0.), (8., 0.)]
        targets = [(ARM, 0.), (11., 0.)]
        spots = [(ARM + 0.15, 0.), (11., 0.)]

        self.assertEqual({i: s for i, (s, _) in match(centres, targets, spots).items()}, {0: 100, 1: 101})

    def testDotCleanupKeepsOtherCandidates(self):
        """A cobra with no candidate is removed from the points' candidate lists, and nothing else is."""
        adjacent = [(np.array([1, 2]),), (np.array([0, 2]),), (np.array([0, 1]),)]
        potPointMatch = [[1, 2], [], [2]]
        potCobraMatch = [[], [0], [0, 2]]
        result = mcsRoutines.secondPass([], [0, 1, 2], [], [], [1, 2], potCobraMatch, potPointMatch,
                                        adjacent, np.zeros(3, dtype=int), 0)
        dotCobras, potPointMatch = result[2], result[6]

        self.assertEqual(dotCobras, [1])
        self.assertEqual(potPointMatch[0], [1, 2])

    def testLastPassAfterSingles(self):
        """Pairs taken by distance resolve to the right points once a single-candidate point is assigned."""
        points = np.array([[0, 0., 0.], [1, 2., 0.], [2, 1., 0.]])
        targets = np.array([[1, 0., 0.], [2, 1., 0.], [3, 2., 0.]])
        lists = unassigned(3, 3, [[0], [1, 2], [1, 2]], [[0], [1, 2], [1, 2]])
        result = mcsRoutines.lastPassDist(*lists, points, targets, None, None, 't', np.zeros(3, dtype=int), 0,
                                          np.arange(3))

        self.assertEqual(result[5], [[0], [2], [1]])

    def testLastPassWithNothingLeft(self):
        """lastPassDist accepts empty point or cobra lists."""
        points = np.array([[0, 0., 0.], [1, 1., 0.]])
        for nPoints, nCobras in ((0, 2), (2, 0), (0, 0)):
            lists = unassigned(nPoints, nCobras, [[0], [1]], [[0], [1]])
            mcsRoutines.lastPassDist(*lists, points, points, None, None, 't', np.zeros(2, dtype=int), 0,
                                     np.arange(2))


if __name__ == '__main__':
    unittest.main()
