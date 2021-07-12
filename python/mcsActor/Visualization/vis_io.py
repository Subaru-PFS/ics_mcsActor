
import numpy as np
import os
import yaml
import glob


def incrementFileName(prefix, suffix):

    # find the next increment for a file name - highest index + 1

    flist = glob.glob('{0}_*_{1}.yaml'.format(prefix, suffix))
    flist.sort()

    if(flist == []):
        i = 0
    else:
        lastFile = flist[-1]
        print(lastFile)
        aa = lastFile.split("_")
        print(aa)
        i = np.int(aa[1])+1
    return '{0}_{1:02d}_{2}.yaml'.format(prefix, i, suffix)


def getLastFile(prefix, suffix):

    flist = glob.glob('{0}_*_{1}.yaml'.format(prefix, suffix))
    flist.sort()

    return flist[-1]


def imageQualitySave(prefix, filename, xs, ys, fxs, fys, peaks, c, diffx, diffy, pts):

    outStruct = {}

    outStruct['filename'] = filename
    outStruct['x'] = xs
    outStruct['y'] = ys
    outStruct['fx'] = fxs
    outStruct['fy'] = fys
    outStruct['peak'] = peaks
    outStruct['dist'] = c
    outStruct['dist_x'] = diffx
    outStruct['dist_y'] = diffy
    outStruct['pts'] = pts

    outFile = incrementFileName(prefix, "image")

    yaml.dump(outStruct, open(outFile, "w"))


def imageQualityLoad(prefix, i):

    if(i < 0):
        outFile = getLastFile(prefix, "image")
    else:
        outFile = '{0}_{1:02d}_image.yaml'.format(prefix, i)
    print(outFile)
    outStruct = yaml.load(open(outFile))

    return outStruct


def seeingSave(prefix, filenames, frameIDs, xArray, yArray, fxArray, fyArray, peakArray, backArray, xm, ym, xAv, yAv, fxAv, fyAv, peakAv, backAv, rmsVal, nMatch, fxFrameAv, fyFrameAv, peakFrameAv, allTrans):

    outStruct = {}

    outStruct['filenames'] = filenames
    outStruct['frameIDs'] = frameIDs
    outStruct['xArray'] = xArray
    outStruct['yArray'] = yArray
    outStruct['fxArray'] = fxArray
    outStruct['fyArray'] = fyArray
    outStruct['peakArray'] = peakArray
    outStruct['backArray'] = backArray
    outStruct['xm'] = xm
    outStruct['ym'] = ym
    outStruct['xAv'] = xAv
    outStruct['yAv'] = yAv
    outStruct['fxAv'] = fxAv
    outStruct['fyAv'] = fyAv
    outStruct['peakAv'] = peakAv
    outStruct['backAv'] = backAv
    outStruct['rmsVal'] = rmsVal
    outStruct['nMatch'] = nMatch
    outStruct['fxFrameAv'] = fxFrameAv
    outStruct['fyFrameAv'] = fyFrameAv
    outStruct['peakFrameAv'] = peakFrameAv
    outStruct['allTrans'] = allTrans
    outStruct['rms'] = rmsVal.mean()

    outFile = incrementFileName(prefix, "seeing")

    yaml.dump(outStruct, open(outFile, "w"))


def seeingLoad(prefix, i):
    if(i < 0):
        outFile = getLastFile(prefix, "seeing")
    else:
        outFile = '{0}_{1:02d}_seeing.yaml'.format(prefix, i)

    outStruct = yaml.load(open(outFile))

    return outStruct
