/*
 * guidersqeuence.c
 *
 *  Created on: Feb 9, 2015
 *      Author: chyan
 */


#include <linux/ioctl.h>
#include <sys/ioctl.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <getopt.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


#include "fitsio.h"
#include "edtinc.h"
#include "pciload.h" /* for strip_newline function */

#define GP0_OUTPUT_ENABLE	_IOW(0x81,0,int)
#define GP0_SET_VALUE		_IOW(0x81,1,int)
#define GP0_GET_VALUE		_IOR(0x81,2,int)

#define S_CLOSE  0
#define S_OPEN   1

/* 
    This value defines the ratio for detecting the shifted frames.  If we see very high
	standard deviation in first 100 pixels, that means these pixels are inserted into 
	the frames.  This cause frame shift in X direction.
*/ 
#define BADRATIO 6

float getStddev(float data[])
{
    float sum = 0.0, mean, standardDeviation = 0.0;
    int i;
    for(i=0; i<10; ++i)
    {
        sum += data[i];
    }
    mean = sum/10;
    for(i=0; i<10; ++i)
        standardDeviation += pow(data[i] - mean, 2);
    return sqrt(standardDeviation/10);
}


int setShutterStatus(int status)
{
	int i, outputData;
	int res;

	i = open("/dev/ttyUSB0", O_RDWR);
	outputData=status;
	res = ioctl(i, GP0_SET_VALUE,&outputData);
	close(i);
	return(0);
}

struct shutterArgs {
    int time;
    int verbose;
};


int OpenShutterTime(void* shutterArgs)
{

	int i, outputData;
	int res;
	int time, verbose;
	int delay = 700;

	time = ((struct shutterArgs*)shutterArgs)->time;
	verbose = ((struct shutterArgs*)shutterArgs)->verbose;

	i = open("/dev/ttyUSB0", O_RDWR);

	/* Adding a delay so that the image readout will start before shutter is opened */
	usleep(delay*1000);
	
	if (verbose){
		printf("Opening shutter for %d msec.\n",time);
	}
	
	outputData=S_OPEN;
	res = ioctl(i, GP0_SET_VALUE,&outputData);
	usleep(time*1000);
	outputData=S_CLOSE;

	res = ioctl(i, GP0_SET_VALUE,&outputData);

	close(i);

	return(0);
}

int WriteFitsImage(char *filename, int height, int width, u_char * image_p, int exptime)
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii;
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
    	fprintf(stderr, "Error: (%s:%s:%d) can not get create image %s in disk "
    			".\n", __FILE__, __func__, __LINE__,filenmae);
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
    exposure = exptime;
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
      
      return 0;
   }
   return (double)(ts.tv_sec + ts.tv_nsec / 1000000000.0);
}

void delay(unsigned int mseconds)
{
    clock_t goal = mseconds + clock();
    while (goal > clock());
}

void strupp(char* lower)
{
	while (*lower = toupper(*lower)) lower++;
}


/* Print out the proper program usage syntax */
static void
printUsageSyntax(char *prgname) {
   fprintf(stderr,
	   "Star guider sequence.\n"
	   "Usage: %s <INPUT> <OUTPUT> [options...]\n"
		"	-h, --help   display help message\n"
		"	-f, --file   name of FITS file to be saved.\n"
		"	-e, --etype  exposure type [flat|object|dark|object].\n"
		"	-t, --exptime  shutter time.\n"
		"	-c, --coadd  produce a co-added image.\n"
		"	-v, --verbose  turn on verbose.\n"
		, prgname);

}


int main(int argc, char *argv[]){
	pthread_t id;

	int    opt;
	int    unit = 0;
	int    channel = 0 ;
	int    s_height,s_width,s_depth;
	int    imagesize;
	int    verbose=0, loops;
	int    i,ii;
	int    exptime = 0;
	int    ret;
	int    coadd = 0;
	int    flag = 0;

	EdtDev *pdv_p = NULL;

	u_char **bufs;
    u_char *image_p=NULL;

	float  *stddev_img = NULL;
	float  *stddev_sec = NULL;

	u_char *coaddframe=NULL;
	float  nloops = 0;

    char   *file=NULL;
	char   *unitstr = "0";
	char   edt_devname[256];
    char   errstr[64];
    char   string[256];
    char   *etype=NULL;

	char etype_list[5][20]={"DARK","BIAS","FLAT","OBJECT","TEST"};

	struct shutterArgs *args = (struct shutterArgs *)malloc(sizeof(struct shutterArgs));
	//FILE *fp;

    double shutter_ts,start_ts,save_ts;
    double dtime;
	/** Check the total number of the arguments */
	struct option longopts[] = {
	     {"file" ,0, NULL, 'f'},
	     {"exptime" ,0, NULL, 't'},
	     {"etype" ,0, NULL, 'e'},
	     {"coadd" ,0, NULL, 'c'},
		 {"verbose",0, NULL, 'v'},
		 {"help", 0, NULL, 'h'},
		 {0,0,0,0}};

	while((opt = getopt_long(argc, argv, "e:f:l:t:vhc",
	   longopts, NULL))  != -1){
	      switch(opt) {
	         case 'e':
	               etype = optarg;
	               break;
	         case 't':
	               exptime = atoi(optarg);
	               break;
	         case 'f':
	               file = optarg;
	               break;
	         case 'v':
	               verbose = 1;
	               break;
			 case 'c':
	               coadd = 1;
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

	if (etype == NULL){
		fprintf(stderr, "Warning: (%s:%s:%d) ETYPE is not set, use OBJECT."
		"\n", __FILE__, __func__, __LINE__);
		etype=etype_list[3];
	}
	if (file == NULL){
		fprintf(stderr, "Warning: (%s:%s:%d) there is no FITS file name specified, use \"exposureXX.fits\"."
		"\n", __FILE__, __func__, __LINE__);
		file="exposure";
	}
	
	/* Determin the total loops needed if exposure time is set */
	nloops=ceil(((float)exptime/800.0))+1;
	loops = (int)nloops;
	if (loops > 1){
		if (verbose) printf("take %i exposures.\n",loops);
	}

	strupp(etype);
	for (i=0;i<5;i++){
		if (strcmp(etype, etype_list[i]) == 0){
			flag++;
		}
	}
	if (flag == 0){
		fprintf(stderr, "Error: (%s:%s:%d) ETYPE setting incorrect. "
				"Check ETYPE setting.\n", __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}


	if (strcmp(etype, "DARK") == 0 || strcmp(etype, "BIAS") == 0){
		if (verbose) {
			printf("Exposure type = %s .\n",etype);
			printf("Setting shutter time = 0\n");
		}
		exptime=0;
	} else {
		if (verbose) printf("Exposure type = %s .\n",etype);
	}

	//return EXIT_SUCCESS;


	/* Start to establish EDT connection */
	unit = edt_parse_unit_channel(unitstr, edt_devname, "pdv", &channel);
	
	/*
	 *   Waiting for shutter
 	 */
	shutter_ts = getClockTime();
	//sleep(1);

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

	s_height=pdv_get_height(pdv_p);
	s_width=pdv_get_width(pdv_p);
    s_depth = pdv_get_depth(pdv_p);
    imagesize = pdv_get_imagesize(pdv_p);
	
	image_p=pdv_alloc(pdv_image_size(pdv_p));

	if (verbose) printf("Image size --> Height = %i Width= %i\n", s_height, s_width);
	if (verbose) printf("Total pixels = %i\n",pdv_image_size(pdv_p));

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
	pdv_multibuf(pdv_p, 4);
	//pdv_start_image(pdv_p);
	
	bufs = (u_char **)malloc(loops * sizeof(u_char *));

	// Allocated the array for standard deviation check.
	stddev_img = (float *)malloc(100 * sizeof(float));
	stddev_sec = (float *)malloc(100 * sizeof(float));

    for (i=0; i<loops; i++){
		if ((bufs[i] = edt_alloc(imagesize)) == NULL){
	    	printf("buffer allocation FAILED (probably too many images specified)\n");
	    	exit(1);
		}
    }
	args->time = exptime;
	args->verbose = verbose;

    if (exptime > 0){
    	ret=pthread_create(&id,NULL,(void *) OpenShutterTime, (void *)args);
    }

	(void) edt_dtime();		/* init time for check */
	pdv_start_images(pdv_p, loops);

    for (i=0; i<loops; i++){
		image_p = pdv_wait_image(pdv_p);
		
		/*   Try to detecting the image shifting by calculating the stand deviation of 
		 *    the fisrt 100 pixel. 
		 */
		for (ii=0; ii<200; ii+=2){
			stddev_sec[ii/2] = image_p[ii] | (image_p[ii+1] << 8);
			stddev_img[ii/2] = image_p[200+ii] | (image_p[200+ii+1] << 8);
		}

		if (getStddev(stddev_sec)/getStddev(stddev_img) > BADRATIO ){
			fprintf(stderr, "Warning: (%s:%s:%d) Image shift detected. "
				"Re-issue exposure command.\n", __FILE__, __func__, __LINE__);
			image_p = pdv_wait_image(pdv_p);
		} 

		memcpy(bufs[i], image_p, imagesize);
    }


    dtime = edt_dtime();
	if (verbose){
    	printf("Image reading finished with %f frames/sec\n",(double) (loops) / dtime);
	}
	
	/* Finishing the shutter thread */
    if (exptime > 0){
    		pthread_join(id,NULL);
    }
	/* Stacking frames if the flag is set */
	if (coadd){
		if (verbose) printf("Coadding all frames.\n");

		u_short *coaddshorts; /* alias. */


		coaddframe = (u_char *)calloc(imagesize, sizeof(u_char));
		coaddshorts= (u_short *)coaddframe;

		//sprintf(string,"%s","coadd.fits");
		sprintf(string,"%s",file);

		for(i=0;i<nloops;i++) {
			for (ii=0;ii<imagesize;ii+=2){
			  u_short pixel;

			  pixel = bufs[i][ii] | (bufs[i][ii+1] << 8);
			  coaddshorts[ii/2] += pixel;
			}
		}

		WriteFitsImage(string, s_height, s_width, coaddframe, exptime);
	} else {

		i=0;
		while(loops) {
			start_ts = getClockTime();

			sprintf(string,"%s%04i%s",file,i+1,".fits");

			/*process and/or display image previously acquired here*/
			WriteFitsImage(string, s_height, s_width,bufs[i], 800);

			if (verbose){
				save_ts=getClockTime();;
				fprintf(stdout,"%02i Image saving runtime = %f\n",i+1, save_ts-start_ts);
			}

			loops--;
			i++;

		}
	}


	/* Free imaeg blocks */
	for (i=0; i<loops; i++){
		edt_free(bufs[i]);
	}

	//pdv_free(image_p);
	free(coaddframe);
	free(bufs);
	free(stddev_img);
	free(stddev_sec);
	pdv_close(pdv_p);

	if (verbose){
		printf("Exposure sequence is done\n");
	}
	return EXIT_SUCCESS;

}

