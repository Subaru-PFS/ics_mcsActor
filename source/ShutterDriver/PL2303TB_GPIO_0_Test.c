#include <linux/ioctl.h>
#include <stdio.h>
#include <errno.h> 
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>

#define GP0_OUTPUT_ENABLE	_IOW(0x81,0,int)
#define GP0_SET_VALUE		_IOW(0x81,1,int)
#define GP0_GET_VALUE		_IOR(0x81,2,int)

int main(int argc, char* argv[])
{
	int i;
	int data;
	int inputData=0;
	int iOUTPUT_ENABLE=0;
	int iGPIO_Number=0;
	int outputData=0;
	char nextTest;
	//unsigned short buf;
	int res;
	
	i = open("/dev/ttyUSB0", O_RDWR);
	if (i==-1) {
		printf("Cannot open device!!!\n");
		return -1;
	}
	
	
	
	while (1)
       {
	    printf("Q1: PLease make a choice: GPIO_number[0]:\n");
	    scanf("%d",&iGPIO_Number);
	    printf("Q2: PLease make a choice: I/O mode[0,1],  '0': input(disable), '1': output(enable)\n");
	    scanf("%d",&iOUTPUT_ENABLE);

            switch (iGPIO_Number) { // S1
		case 0:
	          	res = ioctl(i, GP0_OUTPUT_ENABLE,&iOUTPUT_ENABLE);
			  break;		
		default:
		      	res = ioctl(i, GP0_OUTPUT_ENABLE,&iOUTPUT_ENABLE);
			 break;	
	    } //edn // S1
		
	   if (res==-1) {
		  printf("Ans1&2: ioctl error GPIO_%d OUTPUT_ENABLE, error#=%d\n",iGPIO_Number, errno);
		 break;	
	   }
	   else {
		  printf("Ans1&2: Set GPIO_%d_OUTPUT_ENABLE to %d\n",iGPIO_Number, iOUTPUT_ENABLE);
	   }

        if(iOUTPUT_ENABLE==1){  //S3
		// GPIO output mode, Set H or L
               printf("GPIO_%d is output mode,\n",iGPIO_Number);
	       printf("Q3: PLease make a choice: OUTPUT_Value[0,1] ,1 is High level, 0 is low level.\n");
	       scanf("%d",&outputData);

               if(outputData>=2) outputData =1;
		
	       switch (iGPIO_Number) {
		case 0:
	        	  res = ioctl(i, GP0_SET_VALUE,&outputData);
			  break;
		default:
			 res = ioctl(i, GP0_SET_VALUE,&outputData);
			 break;
	       }//end switch (iGPIO_Number)

            if (res==-1) {
		     printf("Ans3: ioctl error GPIO_%d SET_VALUE, error#=%d\n",iGPIO_Number, errno);
			  break;
	        }	
	        else {
		     printf("Ans3: Set GPIO_%d SET_VALUE to %d\n",iGPIO_Number, outputData);
			} //end if

        } //end if(iOUTPUT_ENABLE==1){ //s3
        else{ //S4

         //  GPIO input mode, Get Data
	    printf("GPIO_%d is input mode, Getting Data......\n",iGPIO_Number);
	     
	    switch (iGPIO_Number) {
		case 0:
	        	  res = ioctl(i, GP0_GET_VALUE,&inputData);
			  break;	
		default:
			  res = ioctl(i, GP0_GET_VALUE,&inputData);
			 break;
	        }
            
		if (res==-1) {
		       printf("ioctl error GPIO_%d_GET_VALUE, error#=%d\n",iGPIO_Number, errno);
			   break;
	        }	
	        else {
		      data =inputData;
		     printf("Get GPIO_%d_GET_VALUE, return=%d\n", iGPIO_Number,data);
		}

         }// end  //  GPIO input mode, Get Data
	
	    printf("Q4: Do you want to continue testing GPIO funciton: Yes[y] or No[n]\n");
	    scanf("%s",&nextTest);
		
	    /*CODE TO EXIT THE LOOP GOES HERE*/
        if (nextTest == 'n') break;
	} //end while (1)
	close(i);
     return -1;
}
