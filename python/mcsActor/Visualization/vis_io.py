
import numpy as np
import os

def incrementFileName(prefix):

    #find the next increment for a file name

    i = 0
    while os.path.exists('{0}_{1:0d})_saved.dat'.format(prefix,i)):
        i += 1

    return '{0}_{1:0d}_saved.dat'.format(name,i)
        

def imageQualitySave(prefix,filename,cenParams,cenOut,distOut):

    outFile=incrementFileName(prefix)
    ff=open(outfile,"w")

    npoints=len(cenOut[0])

    
    print("#Input Filename")
    print("#FWHM,boxsize,thresh,rl,rh,sl,sh")
    print("#npoints x  [x,y,fx,fy,peak,back,c,c1,diffx,diffy,pts1,pts2"])
    print(filename)
    print(cenParams[0],cenParams[1],cenParams[2],cenParams[3],cenParams[4],cenParams[5],cenParams[6],cenParams[7])
    for i in range(npoints):

        for j in range(6):
            print(cenOut[j][i])
        for j in range(6):
            print(cenOut[j][i])

    print()

    ff.close()

def seeingSave(prefix,filename,cenParams,cenOut,distOut):


    #nframes, npoints
    
    #registered centroids
    xArray,yArray,fxArray,fyArray,backArray,peakArray

    #averages by frame
    xAv,yAv,fxAv,fyAv,peakAv,backAv,rmsVal,nMatch

    #averages by frame
    xdAll,ydAll,sxAll,syAll,rotAll,fxFrameAv,fyFrameAv,peakFrameAv

    print("#Input Filenames")
    print("#FWHM,boxsize,thresh,rl,rh,sl,sh")


def rotationSave():

    
