from cpython cimport array
import array
from cython cimport view
import astropy.io.fits as pyfits
import numpy as np
cimport numpy as np
import ctypes

cdef extern from "centroid_types.h" nogil:
        struct centroids:
                double x;
                double y;
                double fx;
                double fy;
                double back;
                double peak;
                double qual;


cdef extern from "centroid_types.h" nogil:
        struct fibreid:
                double xp;     
                double yp;               
                double xt;                
                double yt;       
                double xc;       
                double yc;       
                double x;       
                double y;       
                double peak;
                double back;
                double fx;
                double fy;
                int qual;
                int idnum;       
    

cdef extern from "fibreid.h" nogil:
     fibreid *centroid_coarse(int *image, int *arc_image,double *homes,int n_x,int n_y,int hmin, double fwhm, int boxsize,int nhome,  int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh,int verbose)

cdef extern from "fibreid.h" nogil:
     fibreid *centroid_fine(int *image, double *homes,double *xp,double *yp,int n_x,int n_y,int hmin, double fwhm, int boxsize,int nhome,  int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh,int verbose);

cdef extern from "centroid.h" nogil:
     centroids *centroid(int *image, int n_x, int n_y, int hmin, double fwhm,int boxsize,int *npoint, int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh,int verbose)

def centroid_fine_call(np.ndarray[int, ndim=2, mode="c"] image,np.ndarray[double, ndim=1, mode="c"] homes, np.ndarray[double, ndim=1, mode="c"] xp,np.ndarray[double, ndim=1, mode="c"] yp, double fwhm, int hmin, int boxsize, int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh,int verbose):

    cdef char *cp
    cdef fibreid *vals

    with nogil:
        vals=centroid_fine(<int *>image.data,<double *>homes.data,<double *>xp.data,<double *>yp.data,image.shape[1],image.shape[0],hmin,fwhm,boxsize,homes.shape[0]/2,fittype,sharpLow,sharpHigh,roundLow,roundHigh,verbose)


    cp = <char *> vals

    return cp[:sizeof(fibreid) * homes.shape[0]/2]


def centroid_coarse_call(np.ndarray[int, ndim=2, mode="c"] image,np.ndarray[int, ndim=2, mode="c"] arcimage,np.ndarray[double, ndim=1, mode="c"] homes, double fwhm, int hmin, int boxsize, int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh,int verbose):

    cdef char *cp
    cdef fibreid *vals

    with nogil:
        vals=centroid_coarse(<int *>image.data,<int *>arcimage.data,<double *> homes.data,image.shape[1],image.shape[0],hmin,fwhm,boxsize,homes.shape[0]/2,fittype,sharpLow,sharpHigh,roundLow,roundHigh,verbose)


    cp = <char *> vals

    return cp[:sizeof(fibreid) * homes.shape[0]/2]

def centroid_only(np.ndarray[int, ndim=2, mode="c"] image, double fwhm, int hmin, int boxsize, int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh,int verbose):

    """

    Take an image, return the centroids
        
    """ 

    cdef centroids *val
    cdef char *cp
    cdef int npoint[1]

    #The function call
    with nogil:
         vals=centroid(<int *>image.data,image.shape[1],image.shape[0],hmin,fwhm, boxsize, npoint, fittype, sharpLow, sharpHigh, roundLow, roundHigh, verbose)

    #convert the output into a buffer string
    cp = <char *> vals

    #and into a numpy record array
    
    result=np.frombuffer(cp[:sizeof(centroids) * npoint[0]],dtype=[('x','<f8'),('y','<f8'),('fx','<f8'),('fy','<f8'),('back','<f8'),('peak','<f8'),('qual','<f8')])

    return result
