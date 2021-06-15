/*
 * guidersqeuence.c
 *
 *  Created on: Feb 9, 2015
 *      Author: chyan
 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

//#include "sgc/sgc.h"
#include "fitsio.h"
#include "edtinc.h"
#include "pciload.h" /* for strip_newline function */

int WriteFitsImage(char *filename, int height, int width, u_char * image_p)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj;
    long  fpixel, nelements, exposure;
    unsigned short *array=NULL;

    /* initialize FITS image parameters */
    //char filename[] = "atestfil.fit";             /* name for new FITS file */
    int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
    long naxis    =   2;  /* 2-dimensional image                            */
    long naxes[2] = { width, height };   /* image is 640 pixels wide by 512 rows */

    /* allocate memory for the whole image */
    array = (unsigned short *)malloc( naxes[0] * naxes[1]
                                        * sizeof( unsigned short ) );

    /* initialize pointers to the start of each row of the image */
    for( ii=1; ii<naxes[1]; ii++ )
      array[ii] = array[ii-1] + naxes[0];

    remove(filename);               /* Delete old file if it already exists */

    status = 0;         /* initialize status before calling fitsio routines */

    fits_create_file(&fptr, filename, &status); /* create new FITS file */
    if (status != 0) {
    	fprintf(stderr, "Error: (%s:%s:%d) can not get create image in disk "
    			".\n", __FILE__, __func__, __LINE__);
    	exit(1);
    }

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    fits_create_img(fptr,  bitpix, naxis, naxes, &status);
    if (status != 0) {
    	fprintf(stderr, "Error: (%s:%s:%d) can not get create image in disk "
    			".\n", __FILE__, __func__, __LINE__);
        exit(1);
    }
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* initialize the values in the image with a linear ramp function */
    /*for (jj = 0; jj < naxes[1]; jj++)
    {   for (ii = 0; ii < naxes[0]; ii++)
        {
            //sarray[jj][ii]=image_p[jj*naxes[1]+ii];
        }
    }
*/
    memcpy(array, image_p, nelements*sizeof(unsigned short));

    for (ii=0;ii<nelements;ii++){
    	if (array[ii] > 65535) array[ii]=65535;
    }


    fpixel = 1;                               /* first pixel to write      */
    /* write the array of unsigned integers to the FITS file */
    fits_write_img(fptr, TUSHORT, fpixel, nelements, (unsigned short *)array, &status);
    //fits_write_img(fptr, TUSHORT, fpixel, nelements, (unsigned short *)image_p, &status);
    if (status != 0) {
    	fprintf(stderr, "Error: (%s:%s:%d) can not get close image in disk "
    			".\n", __FILE__, __func__, __LINE__);
    	        exit(1);
    }

    free( array );  /* free previously allocated memory */

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    exposure = 1500.;
    fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Total Exposure Time", &status);
    if (status != 0) {
		fprintf(stderr, "Error: (%s:%s:%d) can not get close image in disk "
				".\n", __FILE__, __func__, __LINE__);
		exit(1);
    }

    fits_close_file(fptr, &status);                /* close the file */
    if (status != 0) {
    	fprintf(stderr, "Error: (%s:%s:%d) can not get close image in disk "
    			".\n", __FILE__, __func__, __LINE__);
        	exit(1);
     }

    return(0);
}


/*
 * Get the current clock timestamp
 */
static double
getClockTime(void) {

   struct timespec ts;

   if (clock_gettime(CLOCK_REALTIME, &ts) == -1) {
      //cfht_logv(CFHT_MAIN, CFHT_LOGONLY,
		//"(%s:%d) unable to get clock timestamp : %s (errno=%d)",
		//__FILE__, __LINE__, strerror(errno), errno);
      return 0;
   }
   return (double)(ts.tv_sec + ts.tv_nsec / 1000000000.0);
}

void delay(unsigned int mseconds)
{
    clock_t goal = mseconds + clock();
    while (goal > clock());
}


/* Print out the proper program usage syntax */
static void
printUsageSyntax(char *prgname) {
   fprintf(stderr,
	   "Star guider sequence.\n"
	   "Usage: %s <INPUT> <OUTPUT> [options...]\n"
		"	-h, --help   display help message\n"
		"	-f, --file   name of FITS file to be saved.\n"
		"	-l, --loops  number of exposure loops.\n"
		"	-v, --verbose  turn on verbose.\n"
		, prgname);

}


int main(int argc, char *argv[]){
	int    opt;
	int    unit = 0;
	int    channel = 0 ;
	int    s_height,s_width,s_depth;
	int    imagesize;
	int    verbose=0,loops=1;
	int    i,ii;
	int    numbufs = 4;
	int    started;
	int    timeouts =0 , last_timeouts = 0;
	int    recovering_timeout = FALSE;

	EdtDev *pdv_p = NULL;

	u_char **bufs;
    u_char *image_p=NULL;

    char   *file=NULL;
	char   *unitstr = "0";
	char   edt_devname[256];
    char   errstr[64];
    char   string[256];

    double shutter_ts,start_ts,end_ts,save_ts;
    double dtime;
	/** Check the total number of the arguments */
	struct option longopts[] = {
         {"loops" ,0, NULL, 'l'},
	     {"file" ,0, NULL, 'f'},
		 {"verbose",0, NULL, 'v'},
		 {"help", 0, NULL, 'h'},
		 {0,0,0,0}};

	while((opt = getopt_long(argc, argv, "f:l:vh",
	   longopts, NULL))  != -1){
	      switch(opt) {
	         case 'l':
	               loops = atoi(optarg);
	               break;
	         case 'f':
	               file = optarg;
	               break;
	         case 'v':
	               verbose = 1;
	               break;
	         case 'h':
	               printUsageSyntax(argv[0]);
	               exit(EXIT_FAILURE);
	               break;
	         case '?':
	               printUsageSyntax(argv[0]);
	               exit(EXIT_FAILURE);
	               break;
	      }
	}


	/** Print the usage syntax if there is no input */
	if (argc < 2 ) {
		printUsageSyntax(argv[0]);
		return EXIT_FAILURE;
	}

	if (file == NULL){
		fprintf(stderr, "Warning: (%s:%s:%d) there is no FITS file name specified, use \"exposureXX.fits\"."
		"\n", __FILE__, __func__, __LINE__);
		file="exposure";
	}

	if (loops > 1){
		if (verbose) printf("take %i exposures.\n",loops);
	}

	/* Start to establish EDT connection */
	unit = edt_parse_unit_channel(unitstr, edt_devname, "pdv", &channel);

	/*
	* pdv_open_channel is just pdv_open with an extra argument, the channel,
	* which is normally 0 unless you're running multiple cameras (2nd base
	* mode on a PCI DV C-Link or dasy-chained RCI boxes off a single PCI
	* FOI)
	*/

	if ((pdv_p = pdv_open_channel(edt_devname, unit, channel)) == NULL){
		fprintf(stderr, "Error:pdv_open(%s%d_%d)", edt_devname, unit, channel);
		pdv_perror(errstr);
		return EXIT_FAILURE;
	}

	pdv_flush_fifo(pdv_p);

	s_height=pdv_get_height(pdv_p);
	s_width=pdv_get_width(pdv_p);
	s_depth = pdv_get_depth(pdv_p);
	imagesize = pdv_get_imagesize(pdv_p);
	
	//image_p=pdv_alloc(pdv_image_size(pdv_p));

	if (verbose) printf("Image size --> Height = %i Width= %i\n", s_height, s_width);

	if (s_height<1 && s_width<1){
		fprintf(stderr, "Error: (%s:%s:%d) image size incorrect. "
		"Check ROI setting.\n", __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;

	}

		/* The number of buffers is limited only by the amount of host memory available,
		* up to approximately 3.5GBytes (or less, depending on other OS use of the low
		* 3.5 GB of memory). Each buffer has a certain amount of overhead, so setting
		* a large number, even if the images are small, is not recommended. Four is
		* the recommended number: at any time, one buffer is being read in, one buffer
		* is being read out, one is being set up for DMA, and one is in reserve in case
		* of overlap. Additional buffers may be necessary with very fast cameras;
		* 32 will almost always smooth out any problems with really fast cameras, and
		* if the system can't keep up with 64 buffers allocated, there may be other problems.
		*
		*/
	for (i=0; i<loops; i++){
		pdv_flush_fifo(pdv_p);

		pdv_multibuf(pdv_p, 4);
		if (pdv_p->dd_p->force_single){
			pdv_start_image(pdv_p);
			started = 1;
		}else{
			pdv_start_images(pdv_p, numbufs);
			started = numbufs;
		}


		if (verbose) fprintf(stdout,"Setting time-out limit %i\n", 100000);
		pdv_set_timeout(pdv_p, 100000);

		(void) edt_dtime();		/* init time for check */
	
		printf("getting image %d\r", i + 1);
        fflush(stdout);
		image_p=pdv_alloc(pdv_image_size(pdv_p));
		
		image_p = pdv_wait_image(pdv_p);
		if (i < loops - started){
             pdv_start_image(pdv_p);
        }
		dtime = edt_dtime();

        timeouts = pdv_timeouts(pdv_p);

         /*
          * check for timeouts or data overruns -- timeouts occur when data
          * is lost, camera isn't hooked up, etc, and application programs
          * should always check for them. data overruns usually occur as a
          * result of a timeout but should be checked for separately since
          * ROI can sometimes mask timeouts
          */
        if (timeouts > last_timeouts){
             /*
              * pdv_timeout_cleanup helps recover gracefully after a timeout,
              * particularly if multiple buffers were prestarted
              */
             pdv_timeout_restart(pdv_p, TRUE);
             last_timeouts = timeouts;
             recovering_timeout = TRUE;
             printf("\ntimeout....\n");
        } else if (recovering_timeout){
             pdv_timeout_restart(pdv_p, TRUE);
             recovering_timeout = FALSE;
             printf("\nrestarted....\n");
        }

		if (loops == 1){
			sprintf(string,"%s",file);
		} else {
			sprintf(string,"%s%04i%s",file,i+1,".fits");
		}
        WriteFitsImage(string, s_height, s_width,image_p);
		if (verbose) fprintf(stdout,"filename saved as %s\n", string);
		
		if (verbose) printf("%f frames/sec\n", dtime);


	}
    
	pdv_close(pdv_p);

    
	printf("done\n");

	return EXIT_SUCCESS;


}

