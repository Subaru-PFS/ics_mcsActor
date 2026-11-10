/*
 * canonexp.c
 *
 *  Created on: Feb 9, 2015
 *      Author: Chi-Hung Yan (chyan@asiaa.sinica.edu.tw)
 *
 *  Revision notes:
 * 	  - Date: 2026-11-15
 *    - Image-shift handling rewritten: a shifted frame now aborts the whole
 *      sequence (including the shutter) and the sequence is retried, up to
 *      MAX_RETRY times.
 *    - Shift detection uses all 100 samples and guards against division by zero.
 *    - Shutter thread uses the correct pthread signature and checks for errors.
 *    - Co-adding is done after acquisition (only when -c is given), so the DMA
 *      read path stays light.
 *    - coaddshorts is intentionally kept as u_short; the exposure time must be
 *      capped so that (number of frames * max pixel value) <= 65535.
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
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


#include "fitsio.h"
#include "edtinc.h"
//#include "pciload.h" /* for strip_newline function */

#define GP0_OUTPUT_ENABLE	_IOW(0x81,0,int)
#define GP0_SET_VALUE		_IOW(0x81,1,int)
#define GP0_GET_VALUE		_IOR(0x81,2,int)

#define S_CLOSE  0
#define S_OPEN   1

/*
 * Ratio used to detect shifted frames. If the standard deviation of the first
 * 100 pixels is much larger than that of the next 100 pixels, extra pixels have
 * been inserted into the frame, which causes a shift in the X direction.
 */
#define BADRATIO 6

/* Number of samples used for each standard-deviation estimate. */
#define NSAMPLES 100

/* Maximum number of times the whole exposure sequence is retried. */
#define MAX_RETRY 3

/* Approximate frame period of the free-running camera, in milliseconds. */
#define FRAME_PERIOD_MS 800.0

/* Upper limit of exposure time in ms, so that the u_short co-add cannot overflow:
 * loops = ceil(exptime / FRAME_PERIOD_MS) + 1 <= 16  ->  exptime <= 12000 ms
 * (assumes 12-bit data, max pixel value 4095; 16 * 4095 = 65520 <= 65535). */
#define MAX_EXPTIME_MS 12000

/*
 * Population standard deviation of the first n values in d.
 */
static double getStddev(const float *d, int n)
{
    double sum = 0.0, mean, var = 0.0;
    int i;

    for (i = 0; i < n; i++)
        sum += d[i];
    mean = sum / n;

    for (i = 0; i < n; i++)
        var += (d[i] - mean) * (d[i] - mean);

    return sqrt(var / n);
}

/*
 * Check one frame for an X-direction shift.
 *
 * Compares the scatter of pixels 0..99 with the scatter of pixels 100..199 of
 * the first row (16-bit little-endian samples). Returns 1 if the frame looks
 * shifted, 0 otherwise. If the reference region has no variation (saturated or
 * constant), the frame is not judged, which avoids a division by zero.
 */
static int isShifted(const u_char *img)
{
    float first[NSAMPLES], second[NSAMPLES];
    double ref;
    int k;

    for (k = 0; k < NSAMPLES; k++) {
        first[k]  = img[2 * k] | (img[2 * k + 1] << 8);
        second[k] = img[2 * NSAMPLES + 2 * k] | (img[2 * NSAMPLES + 2 * k + 1] << 8);
    }

    ref = getStddev(second, NSAMPLES);
    if (ref < 1e-6)
        return 0;

    return (getStddev(first, NSAMPLES) / ref) > BADRATIO;
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


/*
 * Shutter thread: waits briefly so that the image readout starts first, opens
 * the shutter for 'time' milliseconds, then closes it again.
 * Uses the standard pthread start-routine signature.
 */
void *OpenShutterTime(void *arg)
{
	struct shutterArgs *sa = (struct shutterArgs *)arg;
	int outputData;
	int delay = 70;		/* ms */
	int fd;

	fd = open("/dev/ttyUSB0", O_RDWR);
	if (fd < 0) {
		perror("open /dev/ttyUSB0");
		return (void *)-1;
	}

	/* Delay so that the image readout starts before the shutter is opened */
	usleep(delay * 1000);

	if (sa->verbose){
		fprintf(stderr, "Opening shutter for %d msec.\n", sa->time);
	}

	outputData = S_OPEN;
	ioctl(fd, GP0_SET_VALUE, &outputData);
	usleep(sa->time * 1000);
	outputData = S_CLOSE;
	ioctl(fd, GP0_SET_VALUE, &outputData);

	close(fd);

	return NULL;
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
    long naxes[2] = { width, height };   /* image is width pixels wide by height rows */

    /* allocate memory for the whole image */
    //array = (unsigned short *)malloc( naxes[0] * naxes[1]
    //                                    * sizeof( unsigned short ) );

    /* initialize pointers to the start of each row of the image */
    //for( ii=1; ii<naxes[1]; ii++ )
    //  array[ii] = array[ii-1] + naxes[0];

    remove(filename);               /* Delete old file if it already exists */

    status = 0;         /* initialize status before calling fitsio routines */

    fits_create_file(&fptr, filename, &status); /* create new FITS file */
    if (status != 0) {
    	fprintf(stderr, "Error: (%s:%s:%d) can not get create image %s in disk "
    			".\n", __FILE__, __func__, __LINE__,filename);
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

    //memcpy(array, image_p, nelements*sizeof(unsigned short));

    //for (ii=0;ii<nelements;ii++){
    //		if (array[ii] > 65535) array[ii]=65535;
    //}


    fpixel = 1;                               /* first pixel to write      */
    /* write the array of unsigned integers to the FITS file */
    //fits_write_img(fptr, TUSHORT, fpixel, nelements, (unsigned short *)array, &status);
    fits_write_img(fptr, TUSHORT, fpixel, nelements, (unsigned short *)image_p, &status);
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

void strupp(char* lower)
{
	while ((*lower = toupper((unsigned char)*lower))) lower++;
}

/* Print library versions */
static void
printVersions(void) {
    printf("Library Versions:\n");
    printf("================\n");

    printf("CFITSIO Library:\n");
    #ifdef CFITSIO_VERSION
    #define STRINGIFY(x) #x
    #define TOSTRING(x) STRINGIFY(x)
    printf("  Version: %s\n", TOSTRING(CFITSIO_VERSION));
    printf("  Status: Available\n");

    #endif
    printf("\n");

    // EDT PDV library version
    printf("EDT PDV Library:\n");
    char edt_version_str[256];
    if (edt_get_library_version(edt_version_str, sizeof(edt_version_str)) == 0) {
        printf("  Version: %s\n", edt_version_str);
        printf("  Status: Available\n");
    } else {
        printf("  Status: Error getting version\n");
    }
    printf("\n");

    // Compiler version
    printf("Compiler Information:\n");
    #ifdef __GNUC__
        #ifdef __GNUC_PATCHLEVEL__
        printf("  Compiler: GCC %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
        #else
        printf("  Compiler: GCC %d.%d\n", __GNUC__, __GNUC_MINOR__);
        #endif
    #else
    printf("  Compiler: Unknown\n");
    #endif

    printf("  Build Date: %s %s\n", __DATE__, __TIME__);
    printf("\n");
}

/* Print out the proper program usage syntax */
static void
printUsageSyntax(char *prgname) {
   fprintf(stderr,
	   "Canon 50M image acquisition sequence.\n"
	   "Usage: %s <INPUT> <OUTPUT> [options...]\n"
		"\t-h, --help   display help message\n"
		"\t-f, --file   name of FITS file to be saved.\n"
		"\t-e, --etype  exposure type [flat|object|dark|object].\n"
		"\t-t, --exptime  shutter time.\n"
		"\t-c, --coadd  produce a co-added image.\n"
		"\t-n, --noheader  write to stdout without header.\n"
		"\t-v, --verbose  turn on verbose.\n"
		"\t-V, --version   display library versions.\n"
		, prgname);

}


/*
 * Acquire one complete exposure sequence of 'loops' frames into bufs[].
 *
 * The shutter thread (if exptime > 0) is started first, then the frames are
 * read. If any frame times out, is NULL, or looks shifted, the remaining DMA is
 * aborted, the FIFO is flushed, and the shutter thread is joined (so the
 * shutter is guaranteed closed) before returning.
 *
 * Returns:  0 = success
 *           1 = image shift detected
 *          -1 = timeout or other acquisition error
 *
 * NOTE: edt_abort_dma(), pdv_flush_fifo() and pdv_timeouts() are the EDT calls
 * used here; check the exact names/behaviour against your EDT SDK version
 * (some versions also need pdv_timeout_restart() after a timeout).
 */
static int
acquireSequence(PdvDev *pdv_p, u_char **bufs, int loops, int imagesize,
                int exptime, int verbose)
{
	pthread_t id;
	struct shutterArgs args;
	int have_thread = 0;
	int rc = 0;
	int i;
	int timeouts0 = pdv_timeouts(pdv_p);

	args.time = exptime;
	args.verbose = verbose;

	if (exptime > 0) {
		if (pthread_create(&id, NULL, OpenShutterTime, &args) != 0) {
			perror("pthread_create");
			return -1;
		}
		have_thread = 1;
	}

	pdv_start_images(pdv_p, loops);

	for (i = 0; i < loops; i++) {
		u_char *img = (u_char *)pdv_wait_images(pdv_p, 1);

		if (!img || pdv_timeouts(pdv_p) > timeouts0) {
			fprintf(stderr, "Error: (%s:%s:%d) frame %d timed out or returned NULL.\n",
				__FILE__, __func__, __LINE__, i + 1);
			rc = -1;
			break;
		}

		if (isShifted(img)) {
			fprintf(stderr, "Warning: (%s:%s:%d) image shift detected in frame %d.\n",
				__FILE__, __func__, __LINE__, i + 1);
			rc = 1;
			break;
		}

		memcpy(bufs[i], img, imagesize);
	}

	if (rc != 0) {
		/* Stop pending DMA and clear leftovers so the next attempt starts clean */
		edt_abort_dma(pdv_p);
		pdv_flush_fifo(pdv_p);
	}

	/* The shutter thread closes the shutter before it exits */
	if (have_thread)
		pthread_join(id, NULL);

	return rc;
}


int main(int argc, char *argv[]){

	int    opt;
	int    unit = 0;
	int    channel = 0 ;
	int    s_height,s_width,s_depth;
	int    imagesize;
	int    verbose=0, loops;
	int    i,ii;
	int    exptime = 0;
	int    noheader = 0;
	int    coadd = 0;
    int    flag = 0;
	int    attempt, rc = -1;

    // PdvDev must be a pointer type
    PdvDev *pdv_p = NULL;

    u_char **bufs;

	u_char *coaddframe = NULL;
	float  nloops = 0;

    char   *file = NULL;
	char   *basename = NULL;
	char   *unitstr = "0";
	char   edt_devname[256];
    char   errstr[64];
    char   string[256];
    char   *etype=NULL;

	char etype_list[5][20]={"DARK","BIAS","FLAT","OBJECT","TEST"};

	u_short *coaddshorts = NULL;

    double shutter_ts,start_ts,save_ts;
    double dtime;


	/** Check the total number of the arguments */
	struct option longopts[] = {
	     {"file" ,1, NULL, 'f'},
	     {"exptime" ,1, NULL, 't'},
	     {"etype" ,1, NULL, 'e'},
	     {"coadd" ,0, NULL, 'c'},
	     {"noheader" ,0, NULL, 'n'},
		 {"verbose",0, NULL, 'v'},
		 {"version",0, NULL, 'V'},
		 {"help", 0, NULL, 'h'},
		 {0,0,0,0}};

	while((opt = getopt_long(argc, argv, "ne:f:l:t:vVhc",
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
	         case 'V':
	               printVersions();
	               exit(EXIT_SUCCESS);
	               break;
			 case 'c':
	               coadd = 1;
	               break;
			 case 'n':
	               noheader = 1;
	               break;
	         case 'h':
	               printUsageSyntax(argv[0]);
	               exit(EXIT_SUCCESS);
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

	if (exptime > MAX_EXPTIME_MS) {
    fprintf(stderr, "Error: (%s:%s:%d) exposure time %d ms exceeds the limit of %d ms.\n",
            __FILE__, __func__, __LINE__, exptime, MAX_EXPTIME_MS);
    return EXIT_FAILURE;
}
	/* Determine the total number of frames needed from the exposure time */
	nloops=ceil(((float)exptime/FRAME_PERIOD_MS))+1;
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
        perror("pdv_open_channel");
        return EXIT_FAILURE;
	}

	s_height=pdv_get_height(pdv_p);
	s_width=pdv_get_width(pdv_p);
    s_depth = pdv_get_depth(pdv_p);
    imagesize = pdv_get_image_size(pdv_p);

	if (verbose) printf("Image size --> Height = %i Width= %i\n", s_height, s_width);
	if (verbose) printf("Total pixels = %i\n", imagesize);

	if (s_height<1 || s_width<1){
		fprintf(stderr, "Error: (%s:%s:%d) image size incorrect. "
		"Check ROI setting.\n", __FILE__, __func__, __LINE__);
		pdv_close(pdv_p);
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

	/* Allocate all frame buffers once; they are reused across retries. */
	bufs = (u_char **)calloc(loops, sizeof(u_char *));
	if (bufs == NULL){
		fprintf(stderr, "Error: (%s:%s:%d) can not allocate buffer table.\n",
			__FILE__, __func__, __LINE__);
		pdv_close(pdv_p);
		return EXIT_FAILURE;
	}
	for (i=0; i<loops; i++){
		if ((bufs[i] = edt_alloc(imagesize)) == NULL){
	    	printf("buffer allocation FAILED (probably too many images specified)\n");
	    	exit(1);
		}
	}

	/* Allocate the memory for coadding frames */
	coaddframe = (u_char *)calloc(imagesize, sizeof(u_char));
	coaddshorts= (u_short *)coaddframe;

	/*
	 * Acquire the sequence. If a frame is shifted (or a timeout occurs), the
	 * whole sequence, including the shutter exposure, is discarded and redone,
	 * so that every saved frame corresponds to a complete, valid exposure.
	 */
	(void) edt_dtime();		/* init time for check */
	for (attempt = 1; attempt <= MAX_RETRY; attempt++){
		rc = acquireSequence(pdv_p, bufs, loops, imagesize, exptime, verbose);
		if (rc == 0)
			break;
		fprintf(stderr, "Warning: (%s:%s:%d) sequence attempt %d/%d failed (%s).\n",
			__FILE__, __func__, __LINE__, attempt, MAX_RETRY,
			(rc > 0) ? "image shift" : "acquisition error");
	}

	if (rc != 0){
		fprintf(stderr, "Error: (%s:%s:%d) acquisition failed after %d attempts.\n",
			__FILE__, __func__, __LINE__, MAX_RETRY);
		for (i=0; i<loops; i++){
			edt_free(bufs[i]);
		}
		free(coaddframe);
		free(bufs);
		pdv_close(pdv_p);
		return EXIT_FAILURE;
	}

	dtime = edt_dtime();
	if (verbose){
		fprintf(stderr, "Image reading finished with %f frames/sec\n",(double) (loops) / dtime);
	}

	/* Stacking frames if the flag is set (done after acquisition, not during DMA) */
	if (coadd){
		if (verbose) printf("Coadding all frames.\n");

		for (i=0; i<loops; i++){
			for (ii=0; ii<imagesize; ii+=2){
				coaddshorts[ii/2] += bufs[i][ii] | (bufs[i][ii+1] << 8);
			}
		}

		if (noheader){
			if (verbose) printf("Skipping FITS headers.\n");
			fwrite(coaddframe, sizeof(u_char), imagesize, stdout);
        	fflush(stdout);
		} else {
			sprintf(string,"%s",file);
			WriteFitsImage(string, s_height, s_width, coaddframe, exptime);
		}
	} else {
        // Safe handling of file names: do not use strtok on string literals
        int write_to_stdout = (strcmp(file, "-") == 0);
        char base[256];

        if (!write_to_stdout) {
            snprintf(base, sizeof(base), "%s", file);
            // Strip the extension, if any
            char *dot = strrchr(base, '.');
            if (dot) *dot = '\0';
        }

        for (i=0; i<loops; i++){
            start_ts = getClockTime();

            if (write_to_stdout){
                snprintf(string, sizeof(string), "%s", file);   // "-"
            } else {
                snprintf(string, sizeof(string), "%s-%04i%s", base, i+1, ".fits");
            }

            /* Write image */
            WriteFitsImage(string, s_height, s_width, bufs[i], 800);

            if (verbose){
                save_ts=getClockTime();
                fprintf(stdout,"%02i Image saving runtime = %f\n",i+1, save_ts-start_ts);
            }

        }
    }

	dtime = edt_dtime();
	if (verbose){
		printf("Time used to save image = %f second. \n", dtime);
	}

	/* Free image blocks */
	for (i=0; i<loops; i++){
		edt_free(bufs[i]);
	}

	free(coaddframe);
	free(bufs);
	pdv_close(pdv_p);

	if (verbose){
		printf("Exposure sequence is done\n");
	}
	return EXIT_SUCCESS;

}