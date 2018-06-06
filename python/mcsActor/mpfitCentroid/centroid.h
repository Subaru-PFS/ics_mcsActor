

#ifndef centroid_H_   /* Include guard */
#define centroid_H_


int *read_fits_int(char *filename, int *nx, int *ny);

int *convol_sep(int *image, int *imagemask, double *kernel, int n_x, int n_y, int nbox,int verbose);  /* function declaration */

int gaussfunc1d(int m, int n,double *p, double *dy, double **dvec, void *vars);

int gaussfunc2d(int m, int n, double *p, double *dy, double **dvec, void *vars);

int *make_mask(int* image,  int nbox, long n_x, long n_y, double mthresh,int verbose);

int *convol_sep(int *image, int *imagemask, double *kernel, int n_x, int n_y, int nbox,int verbose);

struct cand_point *getfit_2d(int *image,int hpoint, int n_x,int n_y, int ix, int iy,double fwhm,int boxsize);

struct cand_point *getfit_1d(int *image,int hpoint, int n_x,int n_y, int ix, int iy,double fwhm,int boxsize);

struct cand_point *getfit_simple(int *image,int hpoint, int n_x,int n_y, int ix, int iy,double fwhm,int boxsize);

int input_values(int *nbox,double **gx,double **c1,int **mask, int *pixels,int *nhalf,double fwhm,int verbose);

double calculate_sharp(int *image,int *mask, int nhalf, int ix, int iy, int n_x, int n_y,int pixels,int nbox);

double calculate_round(int *image,int *mask, int nhalf, int ix, int iy, int n_x, int n_y,int pixels,int nbox,double *c1);

int *get_data(char *filename,long *n_x,long *n_y);

struct centroids *centroid(int *image, int n_x, int n_y, int hmin, double fwhm,int boxsize,int *np, int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh,int verbose);  /* function declaration */

#endif // centroid_H_

