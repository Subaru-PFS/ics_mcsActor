#ifndef spotFinder_H_   /* Include guard */
#define spotFinder_H_

/* status of findSpots, in place of a number of spots */
#define SPOTS_BAD_IMAGE -1   /* no pixel, or more rows than a column histogram can count */
#define SPOTS_FLOODED -2     /* more pixels above the threshold than MAX_BRIGHT_PIXELS */
#define SPOTS_TOO_MANY -3    /* more spots than MAX_SPOTS */
#define SPOTS_NO_MEMORY -4   /* an allocation failed */

typedef struct spotParams
{
  double findSigma;  /* detection threshold, in units of the noise of the frame */
  int nmin;          /* pixels above the threshold a spot holds at least */
  int boxFind;       /* spots closer than boxFind + 1 pixels are taken as one */
  int boxCent;       /* radius of the window of the centroid, pixels */
  double fwhmx;      /* full width at half maximum of the window weighting x, pixels */
  double fwhmy;      /* same for y */
  int maxIt;         /* iterations of the windowed centroid at most */
  double aperture;   /* radius the flux is summed in, pixels */
  double gain;       /* e-/ADU */
  int nThreads;      /* threads used at most */
} spotParams;

typedef struct spot
{
  double x;        /* windowed centroid, pixels; a pixel is centred on its indices */
  double y;
  double x2;       /* second moments of the pixels above the threshold, pixels^2 */
  double y2;
  double xy;
  double peak;     /* highest pixel, ADU above the level of its column */
  double flux;     /* sum within the aperture, ADU above the column levels */
  double fluxerr;  /* statistical error of the flux, ADU */
  double background;  /* background subtracted at the centroid: column level and deviation of the smooth background, ADU */
} spot;

typedef struct frameLevels
{
  double background;  /* median of the column levels, the smooth background's median included, ADU */
  double noise;       /* standard deviation of a pixel, ADU */
  double threshold;   /* detection threshold applied, ADU above the column levels */
} frameLevels;

int columnLevels(const int *image, int n_x, int n_y, int nThreads, double *levels, double *noise);
double apertureSum(const int *image, int n_x, int n_y, double x, double y, double radius, int *nPixels);
spot *findSpots(int *image, int n_x, int n_y, const spotParams *params, frameLevels *levels, int *nSpots);

#endif
