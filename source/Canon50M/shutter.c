#include <linux/ioctl.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>

#include <getopt.h>

#define GP0_OUTPUT_ENABLE	_IOW(0x81,0,int)
#define GP0_SET_VALUE		_IOW(0x81,1,int)
#define GP0_GET_VALUE		_IOR(0x81,2,int)

#define S_CLOSE  0
#define S_OPEN   1

/* Print out the proper program usage syntax */
static void
printUsageSyntax(char *prgname) {
   fprintf(stderr,
	   "Metrology camera shutter control.\n"
	   "Usage: %s [options...]\n"
       "	-o, --open     Open shutter.\n"
       "	-c, --close    Close shutter.\n"
       "	-t, --time     Set timimg for shutter (ms).\n"
	   "	-h, --help     display help message.\n"
	   "	-v, --verbose  turn on verbose.\n"
		, prgname);

}


int main(int argc, char* argv[]){
	int opt;
	int i,verbose=0;
	int shutter=S_CLOSE;
	int time=0;

	int iOUTPUT_ENABLE=1;
	int outputData;
	int res;

	/** Check the total number of the arguments */
	struct option longopts[] = {
		 {"open" ,0, NULL, 'o'},
		 {"close" ,0, NULL, 'c'},
		 {"time" ,1, NULL, 't'},
		 {"verbose",0, NULL, 'v'},
		 {"help", 0, NULL, 'h'},
		 {0,0,0,0}};

	while((opt = getopt_long(argc, argv, "oct:vh",
	   longopts, NULL))  != -1){
	      switch(opt) {
	         case 'o':
	               shutter = S_OPEN;
	               break;
	         case 'c':
	               shutter = S_CLOSE;
	               break;
	         case 't':
	               time = atoi(optarg);
	               shutter = -1;
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

	if (time !=0 && time < 10 ) {
		printf("Number too small. The unit of time is ms!!!\n");
		return EXIT_FAILURE;
	}

	i = open("/dev/ttyUSB0", O_RDWR);
	if (i==-1) {
		printf("Cannot open device!!!\n");
		return -1;
	}

	/* Enable GPIO_1 to output mode*/
	res = ioctl(i, GP0_OUTPUT_ENABLE,&iOUTPUT_ENABLE);
	if (res==-1) {
		printf("Unable to access shutter, check device.\n");
		return EXIT_FAILURE;	
	}

	if (shutter == S_OPEN){
		printf("Opening shutter.\n");
		outputData=S_OPEN;
		res = ioctl(i, GP0_SET_VALUE,&outputData);
	}

	if (shutter == S_CLOSE){
		printf("closing shutter.\n");
		outputData=S_CLOSE;
		res = ioctl(i, GP0_SET_VALUE,&outputData);

	}

	if (time != 0){
		printf("Opening shutter for %d msec.\n",time);
		outputData=1;
		res = ioctl(i, GP0_SET_VALUE,&outputData);
		usleep(time*1000);
		outputData=0;

		res = ioctl(i, GP0_SET_VALUE,&outputData);

	}

	close(i);
	return 0;
}

