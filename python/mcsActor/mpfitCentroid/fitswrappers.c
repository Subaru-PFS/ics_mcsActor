#include "fitsio.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int write_fits_simple_int(char filename[],long nx,long ny,int *image)
{

  /* Simple routine to dump a long integer array to a FITS file.*/

  //FITS VARIABLES
  fitsfile *outfptr;
  int status = 0;  
  
  long naxes[2];  
  int naxis=2;
  int bitpix = SHORT_IMG;

  long npix=nx*ny;

  naxes[0]=nx;
  naxes[1]=ny;

  //CREATE FILE
  fits_create_file(&outfptr, filename, &status);

  //SET UP IMAGE VALUES
  fits_create_img(outfptr,bitpix,naxis,naxes,&status);

  //WRITE
  fits_write_img(outfptr, TINT, 1, npix,image,&status);

  //CLOSE
  fits_close_file(outfptr, &status);
  
  return 0;
}
/*--------------------------------------------------------*/

int write_fits_simple_long(char filename[],long nx,long ny,long *image)
{

  /* Simple routine to dump a long integer array to a FITS file.*/

  //FITS VARIABLES
  fitsfile *outfptr;
  int status = 0;  
  
  long naxes[2];  
  int naxis=2;
  int bitpix = LONG_IMG;

  long npix=nx*ny;

  naxes[0]=nx;
  naxes[1]=ny;

  //CREATE FILE
  fits_create_file(&outfptr, filename, &status);

  //SET UP IMAGE VALUES
  fits_create_img(outfptr,bitpix,naxis,naxes,&status);

  //WRITE
  fits_write_img(outfptr, TLONG, 1, npix,image,&status);

  //CLOSE
  fits_close_file(outfptr, &status);
  
  return 0;
}
/*--------------------------------------------------------*/

int write_fits_simple_longlong(char filename[],long nx,long ny,long *image)
{

  /* Simple routine to dump a long integer array to a FITS file.*/

  //FITS VARIABLES
  fitsfile *outfptr;
  int status = 0;  
  
  long naxes[2];  
  int naxis=2;
  int bitpix = LONGLONG_IMG;

  long npix=nx*ny;

  naxes[0]=nx;
  naxes[1]=ny;

  //CREATE FILE
  fits_create_file(&outfptr, filename, &status);

  //SET UP IMAGE VALUES
  fits_create_img(outfptr,bitpix,naxis,naxes,&status);

  //WRITE
  fits_write_img(outfptr, TLONGLONG, 1, npix,image,&status);

  //CLOSE
  fits_close_file(outfptr, &status);
  
  return 0;
}
/*--------------------------------------------------------*/

int *read_fits_int(char *filename, int *nx, int *ny)
{

  /*routine to read a fits file*/


  fitsfile *afptr;                    //POINTER TO FITS OBJECT
  int status = 0,nfound;      //FITS READ VARIABLES. STATUS MUST = 0
  long fpix[2];                       //FIRST PIXELS OF SUBREGION [X,Y]
  long lpix[2];                       //LAST PIXELS OF SUBREGION [X,Y]
  long naxes[2],npixels;              //AXES DIMENSIONS, NUMBER OF PIXELS

  long inc[2]={1,1};                  //VARIABLE FOR READING IN FITS FILE

  //OPEN FITS FILE
  fits_open_image(&afptr, filename, READONLY, &status);
  //GET SIZE
  fits_read_keys_lng(afptr, "NAXIS", 1, 2, naxes, &nfound, &status);
  npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
  
  //SET THE PIXELS TO ONLY READ THE CENTRAL PART OF THE IMAGE. CURRENTLY HARD
  //CODED. ??SET AS A VARIABLE LATER

  fpix[0]=1;
  fpix[1]=1;
  lpix[0]=naxes[0];
  lpix[1]=naxes[1];

  //GET THE DIMENSIONS OF THE SUBREGION
  (*nx)=naxes[0];
  (*ny)=naxes[1];
  
  //ALLOCATE MEMORY FOR THE IMAGE
  int *image;
  image = malloc((*nx)*(*ny)*sizeof(int));
  
  //READ IN FITS FILE AND CLOSE
  fits_read_subset(afptr, TINT,fpix,lpix,inc,NULL,image,NULL,&status);
  fits_close_file(afptr, &status);

  return image;
}
/*--------------------------------------------------------------------*/
long *read_fits_long(char *filename, long *nx, long *ny)
{

  /*routine to read a fits file*/


  fitsfile *afptr;                    //POINTER TO FITS OBJECT
  int status = 0,nfound;      //FITS READ VARIABLES. STATUS MUST = 0
  long fpix[2];                       //FIRST PIXELS OF SUBREGION [X,Y]
  long lpix[2];                       //LAST PIXELS OF SUBREGION [X,Y]
  long naxes[2],npixels;              //AXES DIMENSIONS, NUMBER OF PIXELS

  long inc[2]={1,1};                  //VARIABLE FOR READING IN FITS FILE

  //OPEN FITS FILE
  fits_open_image(&afptr, filename, READONLY, &status);

  //GET SIZE
  fits_read_keys_lng(afptr, "NAXIS", 1, 2, naxes, &nfound, &status);
  npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
  
  //SET THE PIXELS TO ONLY READ THE CENTRAL PART OF THE IMAGE. CURRENTLY HARD
  //CODED. ??SET AS A VARIABLE LATER

  fpix[0]=1;
  fpix[1]=1;
  lpix[0]=naxes[0];
  lpix[1]=naxes[1];

  //GET THE DIMENSIONS OF THE SUBREGION
  (*nx)=naxes[0];
  (*ny)=naxes[1];
  
  //ALLOCATE MEMORY FOR THE IMAGE
  long *image;
  image = malloc((*nx)*(*ny)*sizeof(long));
  
  //READ IN FITS FILE AND CLOSE
  fits_read_subset(afptr, TLONG,fpix,lpix,inc,NULL,image,NULL,&status);
  fits_close_file(afptr, &status);

  return image;
}
/*--------------------------------------------------------------------*/
