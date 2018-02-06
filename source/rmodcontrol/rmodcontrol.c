/*
 * guidersqeuence.c
 *
 *  Created on: Feb 9, 2015
 *      Author: chyan
 */

#define SERIALTIMEOUT  6
#define SERBUFSIZE    512
#define UNIT       0
#define BAUD       9600
#define CHANNEL    0

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "fitsio.h"
#include "edtinc.h"
#include "pciload.h" /* for strip_newline function */

void
printify(u_char *buf, u_char *dest, int n, int show_nonprint)
{
    int     has_lf = FALSE;
    int     i, j;

    for (i = 0; i < n; i++)
    {
        if (buf[i] == 10)
            has_lf = TRUE;
    }

    for (i = 0, j = 0; i < n; i++)
    {

        if (buf[i] == 13)
        {
            if (!has_lf)
            dest[j++] = '\n';
        }
        else if (buf[i] == 10 || (buf[i] >= ' ' && buf[i] < 127))
            dest[j++] = buf[i];
        else if (show_nonprint)
        {
            switch (buf[i])
            {
            case 0x06:
            strcpy((char *) dest + j, "<ACK>");
            break;
            case 0x15:
            strcpy((char *) dest + j, "<NAK>");
            break;
            case 0x02:
            strcpy((char *) dest + j, "<STX>");
            break;
            case 0x03:
            strcpy((char *) dest + j, "<ETX>");
            break;
            default:
            sprintf((char *) dest + j, "<%02x>", buf[i]);
            }

            j = strlen((char *) dest);

        }

        dest[j] = 0;
    }
}

char *sltrim(char *str) {

   char *p;

   if (str == NULL)
      return NULL;

   for (p = str; (*p && isspace(*p)); p++)
      ;

   return p;
}

/* Trim off trailing whitespace in a string */

char *
srtrim(char *str) {

   int i;

   if (str != NULL) {
      for (i = strlen(str) - 1; (i >= 0 && isspace(str[i])); i--)
         ;
      str[++i] = '\0';
   }

   return str;
}

/* Trim off all leading and trailing whitespace from a string */

char * strim(char *str) {

        return srtrim(sltrim(str));
}

char** str_split(char* a_str, const char a_delim)
{
    char** result    = 0;
    size_t count     = 0;
    char* tmp        = a_str;
    char* last_comma = 0;
    char delim[2];
    delim[0] = a_delim;
    delim[1] = 0;
    //printf("%s\n",a_str);
    /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (a_str + strlen(a_str) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(a_str, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) = strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
    }

    return result;
}

void print_ascii_string(char *buf)
{
#if 0
    int     i = 0;
    char   *p = buf;
#endif
    char    tmpbuf[SERBUFSIZE];

    printify((u_char *) buf, (u_char *) tmpbuf, strlen(buf), 1);

#if 0
    while (*p)
    {
        if ((*p == '\r' || *p == '\n')
            && (p != buf) && (*(p - 1) != '\n') && (*(p + 1) != '\n'))
            tmpbuf[i++] = '\n';
        else
            tmpbuf[i++] = *p;
        ++p;
    }
    if (tmpbuf[i - 1] != '\n')
        tmpbuf[i++] = '\n';
    tmpbuf[i] = '\0';

#endif

    fputs(tmpbuf, stdout);
    fflush(stdout);
}



long PdvSerialWriteRead(char *ibuf_p, int verbose, char **reps){
    int   timeout = SERIALTIMEOUT;
    int   ret,i;
    int   length=0;
    char  temp[256]="";
    char  buf[SERBUFSIZE+1];
    u_char  hbuf[SERBUFSIZE];

    u_char  lastbyte,waitc;

	EdtDev *ed;


    /* open a handle to the device     */
    ed = pdv_open_channel(EDT_INTERFACE, UNIT, CHANNEL);
    if (ed == NULL)
    {
        pdv_perror(EDT_INTERFACE);
        return -1;
    }

    pdv_serial_read_enable(ed);

    /** Getting timeout value from EDT card configuration */
    if (timeout < 1) timeout = ed->dd_p->serial_timeout;
    if (verbose) printf("serial timeout %d\n", timeout);

    /** Setting baud rate */
    pdv_set_baud(ed, BAUD);


    /* flush any junk */
    (void) pdv_serial_read(ed, buf, SERBUFSIZE);

    u_int val;
	i = 0;
	sprintf(buf, "%s\r", ibuf_p);
    
    if (verbose) printf("writing <%s>\n", ibuf_p);
    pdv_serial_command(ed, buf);
	/*strip_newline(ibuf_p);
	if (verbose) printf("buf=%s\n",ibuf_p);
	while (*ibuf_p)
	{
		while ((*ibuf_p == ' ') || (*ibuf_p == '\t'))
		++ibuf_p;

		if (*ibuf_p == '\0')
		break;

		if (sscanf(ibuf_p, "%x", &val) != 1){
			printf("error reading input byte %d\n",i);
			i = 0;
			break;
		}

		if (val > 0xff){
			printf("Hex string format error -- expect hex bytes separated by spaces, e.g. '00 a0 ff ...' \n");
			i = 0;
			break;
		}
		else hbuf[i++] = val;

		while ((*ibuf_p != ' ') && (*ibuf_p != '\t') && (*ibuf_p != '\0'))
			++ibuf_p;
	}
*/	/*
	 * using pdv_serial_binary_command instead of
	 * pdv_serial_write because it prepends a 'c' if FOI
	 */
	if (i)
		if (pdv_serial_binary_command(ed, (char *) hbuf, i) !=0){
			fprintf(stderr, "Error: (%s:%s:%d) can not send serial binary "
			"command to camera.\n", __FILE__, __func__, __LINE__);
			return EXIT_FAILURE;
		}

    /*
     * serial_timeout comes from the config file (or -t override flag in
     * this app), or if not present defaults to 500 unless readonly
     * defaults to 60000
     */

	pdv_serial_wait(ed, timeout, 64);

	do{
	 	if (verbose) printf("read serial.  \n");
		ret = pdv_serial_read(ed, buf, SERBUFSIZE);
	    if (verbose)
	        printf("read returned %d\n", ret);

	    if (*buf)
	        lastbyte = (u_char)buf[strlen(buf)-1];

	    if (ret != 0)
	    {
	        buf[ret + 1] = 0;
	        
	        //print_ascii_string(buf);
	        length += ret;
	    }

	    if (ed->devid == PDVFOI_ID)
	        ret = pdv_serial_wait(ed, 500, 0);
	    else if (pdv_get_waitchar(ed, &waitc) && (lastbyte == waitc))
	        ret = 0; /* jump out if waitchar is enabled/received */
	    else ret = pdv_serial_wait(ed, 500, 64);

	}while (ret > 0);

    pdv_close(ed);

    *reps=strim(buf);
    if (verbose) printf("reps = %s. \n",*reps);

	return(0);
}

/*
 *    Get the exposure time in the unit of ms
 */
int GetRmodExptime(int verbose, unsigned long* count){
	int status;
	int num,low,up;
	char expvaluehex[5];

	char start[2];
	char checksum[3];

	char lowbyte[3];
	char upbyte[3];


	char *reps;
	char *string;


	string = "{r0200000000}";
	status = PdvSerialWriteRead(string, verbose, &reps);
	if (verbose) printf("string = %s \n", &reps[0]);
	

	/* Make sure the response is correct */
	memcpy(start, &reps[0], 1 );
	start[1] = '\0';
	if (strcmp(start,"!") != 0){
    	fprintf(stderr, "Error: (%s:%s:%d) camera response is not correct"
    			".\n", __FILE__, __func__, __LINE__);
    	return(1);
    }

    /* Make sure the checksum is correct */
    memcpy(checksum, &reps[11], 2 );
    checksum[2]= '\0';
    num = (int)strtol(checksum, NULL, 16);

    memcpy(upbyte, &reps[7], 2 );
    upbyte[2]= '\0';
    up = (int)strtol(upbyte, NULL, 16);

	memcpy(lowbyte, &reps[9], 2 );
    lowbyte[2]= '\0';
    low = (int)strtol(lowbyte, NULL, 16);

    if (256-low-up != num){
    	fprintf(stderr, "Error: (%s:%s:%d) checksum is not correct"
    			".\n", __FILE__, __func__, __LINE__);
    	return(1);
    }

	/* extracting exposure time in Hex format */
	memcpy(expvaluehex, &reps[7], 4 );
	expvaluehex[4] = '\0';
	num = (int)strtol(expvaluehex, NULL, 16);

 
	if (verbose)  printf("exposure value(hex) = %s \n", expvaluehex);

	

	*count = (int)strtol(expvaluehex, NULL, 16);
	if (verbose) printf("count = %lu \n", *count);

	return(0);

}



/*
 *    Get the image depth
 */
int SetRmodExptime(int verbose, unsigned int count){
	int status;
	int num,low,up;
	char valuehex[5];

	char start[2];
	char checksum[3];

	char lowbyte[3];
	char upbyte[3];


	char reps[256];
	char *string;


	string = "{w0200";
	sprintf(&valuehex[0], "%04x", count);
	if (verbose)  printf(" value(hex) = %s \n", valuehex);

	memcpy(upbyte, &valuehex[0], 2 );
    upbyte[2]= '\0';
    up = (int)strtol(upbyte, NULL, 16);

	memcpy(lowbyte, &valuehex[2], 2 );
    lowbyte[2]= '\0';
    low = (int)strtol(lowbyte, NULL, 16);

	sprintf(&checksum[0], "%02x", 256-low-up);
	sprintf(&reps[0],"%s%s%s%s",string,valuehex,checksum,"}");

    if (verbose) printf("sending command = %s\n",reps);

	status = PdvSerialWriteRead(reps, verbose, &string);
	if (verbose) printf("string = %s \n", &string[0]);

	/* Make sure the response is correct */
	memcpy(start, &string[0], 1 );
	start[1] = '\0';
	if (strcmp(start,"!") != 0){
    	fprintf(stderr, "Error: (%s:%s:%d) camera response is not correct"
    			".\n", __FILE__, __func__, __LINE__);
    	return(1);
    }

	return(0);

}



/*
 *    Get the exposure time in the unit of ms
 */
int GetRmodMode(int verbose, int* count){
	int status;
	int num,low,up;
	char valuehex[5];

	char start[2];
	char checksum[3];

	char lowbyte[3];
	char upbyte[3];


	char *reps;
	char *string;


	string = "{r04000001ff}";
	status = PdvSerialWriteRead(string, verbose, &reps);
	if (verbose) printf("string = %s \n", &reps[0]);
	

	/* Make sure the response is correct */
	memcpy(start, &reps[0], 1 );
	start[1] = '\0';
	if (strcmp(start,"!") != 0){
    	fprintf(stderr, "Error: (%s:%s:%d) camera response is not correct"
    			".\n", __FILE__, __func__, __LINE__);
    	return(1);
    }

    /* Make sure the checksum is correct */
    memcpy(checksum, &reps[11], 2 );
    checksum[2]= '\0';
    num = (int)strtol(checksum, NULL, 16);

    memcpy(upbyte, &reps[7], 2 );
    upbyte[2]= '\0';
    up = (int)strtol(upbyte, NULL, 16);

	memcpy(lowbyte, &reps[9], 2 );
    lowbyte[2]= '\0';
    low = (int)strtol(lowbyte, NULL, 16);

    if (256-low-up != num){
    	fprintf(stderr, "Error: (%s:%s:%d) checksum is not correct"
    			".\n", __FILE__, __func__, __LINE__);
    	return(1);
    }

	/* extracting exposure time in Hex format */
	memcpy(valuehex, &reps[7], 4 );
	valuehex[4] = '\0';
	num = (int)strtol(valuehex, NULL, 16);

 
	if (verbose)  printf(" value(hex) = %s \n", valuehex);

	

	*count = (int)strtol(valuehex, NULL, 16);
	if (verbose) printf("count = %d \n", *count);

	return(0);

}


/*
 *    Get the camera temperture
 */
int GetRmodTemp(int verbose, int* count){
	int status;
	int num,low,up;
	char valuehex[5];

	char start[2];
	char checksum[3];

	char lowbyte[3];
	char upbyte[3];


	char *reps;
	char *string;


	string = "{r0407000000}";
	status = PdvSerialWriteRead(string, verbose, &reps);
	if (verbose) printf("string = %s \n", &reps[0]);
	

	/* Make sure the response is correct */
	memcpy(start, &reps[0], 1 );
	start[1] = '\0';
	if (strcmp(start,"!") != 0){
    	fprintf(stderr, "Error: (%s:%s:%d) camera response is not correct"
    			".\n", __FILE__, __func__, __LINE__);
    	return(1);
    }

    /* Make sure the checksum is correct */
    memcpy(checksum, &reps[11], 2 );
    checksum[2]= '\0';
    num = (int)strtol(checksum, NULL, 16);

    memcpy(upbyte, &reps[7], 2 );
    upbyte[2]= '\0';
    up = (int)strtol(upbyte, NULL, 16);

	memcpy(lowbyte, &reps[9], 2 );
    lowbyte[2]= '\0';
    low = (int)strtol(lowbyte, NULL, 16);

    if (256-low-up != num){
    	fprintf(stderr, "Error: (%s:%s:%d) checksum is not correct"
    			".\n", __FILE__, __func__, __LINE__);
    	return(1);
    }

	/* extracting exposure time in Hex format */
	memcpy(valuehex, &reps[7], 4 );
	valuehex[4] = '\0';
	num = (int)strtol(valuehex, NULL, 16);

 
	if (verbose)  printf(" value(hex) = %s \n", valuehex);

	

	*count = (int)strtol(valuehex, NULL, 16);
	if (verbose) printf("count = %d \n", *count);

	return(0);

}

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




/* Print out the proper program usage syntax */
static void
printUsageSyntax(char *prgname) {
   fprintf(stderr,
	   "RMOD-71MP Control Agent \n"
	   "Usage: %s [options...]\n"
		"	-h, --help   display help message\n"
		"	-e, --etime  set a new exposure time (ms).\n"
		"	-l, --list   list system information.\n"
		"	-v, --verbose  turn on verbose.\n"
		, prgname);

}


int main(int argc, char *argv[]){
	int    opt;
	int    unit = 0;
	int    channel = 0 ;
	int    mode,temp;
	int    verbose=0,list=0;
	int    i,ii;
    int    timeouts = 0;
    int    status;


	EdtDev *pdv_p = NULL;

	u_char **bufs;
    u_char *image_p=NULL;

	char   *unitstr = "0";
	char   edt_devname[256];
    char   errstr[64];
    char   string[256];

    unsigned long etime=0, exptime;


    double shutter_ts,start_ts,end_ts,save_ts;
    double dtime;
	/** Check the total number of the arguments */
	struct option longopts[] = {
         {"list" ,0, NULL, 'l'},
	     {"etime" ,0, NULL, 'e'},
		 {"verbose",0, NULL, 'v'},
		 {"help", 0, NULL, 'h'},
		 {0,0,0,0}};

	while((opt = getopt_long(argc, argv, "e:lvh",
	   longopts, NULL))  != -1){
	      switch(opt) {
	         case 'l':
	               list = 1;
	               break;
	         case 'e':
	               etime=atoi(optarg);
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

	/* Start to establish EDT connection to update time out information */
	unit = edt_parse_unit_channel(unitstr, edt_devname, "pdv", &channel);
	/*
     * pdv_open_channel is just pdv_open with an extra argument, the channel,
     * which is normally 0 unless you're running multiple cameras (2nd base
     * mode on a PCI DV C-Link or dasy-chained RCI boxes off a single PCI
     * FOI)
     */
	if ((pdv_p = pdv_open_channel(edt_devname, unit, channel)) == NULL){
    	fprintf(stderr, "Error:pdv_open(%s%d_%d)", edt_devname, UNIT, CHANNEL);
        pdv_perror(errstr);
        return EXIT_FAILURE;
	}

	if (list == 1){
		printf("Listing camera information.\n");

		status=GetRmodMode(verbose,&mode);
		if (status != 0){
			fprintf(stderr, "Error: (%s:%s:%d) not able to get camera "
			"mode.\n", __FILE__, __func__, __LINE__);
			return EXIT_FAILURE;
		} else {
			if (mode == 1) printf("---> Camera is now under medium mode. \n");
			if (mode == 0) printf("---> Camera is now under base mode. \n");

		}	
		status=GetRmodTemp(verbose,&temp);
		if (status != 0){
			fprintf(stderr, "Error: (%s:%s:%d) not able to get camera "
			"mode.\n", __FILE__, __func__, __LINE__);
			return EXIT_FAILURE;
		} else {
			printf("---> Camera temperture = %d degC.\n", temp);
		}	

		status=GetRmodExptime(verbose,&exptime);
		if (status != 0){
			fprintf(stderr, "Error: (%s:%s:%d) not able to get camera "
			"exposure time.\n", __FILE__, __func__, __LINE__);
			return EXIT_FAILURE;
		} else {

			printf("---> Camera exposure time = %ld ms. \n",exptime);
		}	
	}

	if (etime != 0){
		if (verbose) printf("Set camera exposure time = %ld ms. \n",etime);

		status=SetRmodExptime(verbose,(unsigned int)etime);
		if (status != 0){
			fprintf(stderr, "Error: (%s:%s:%d) not able to get camera "
			"mode.\n", __FILE__, __func__, __LINE__);
			return EXIT_FAILURE;
		}

		status=GetRmodExptime(verbose,&etime);
		if (status != 0){
			fprintf(stderr, "Error: (%s:%s:%d) not able to get camera "
			"exposure time.\n", __FILE__, __func__, __LINE__);
			return EXIT_FAILURE;
		} else {

			printf("---> Camera exposure time = %ld ms. \n",etime);
		}	


	}


	return EXIT_SUCCESS;

}

