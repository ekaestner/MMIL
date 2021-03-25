/********************************************************
 * Copyright (C) 2006 University of California, San Diego
 * ADNI project program: Image Converter: DICOM to MGH
 * File Name:	dcm2mgh.cc 
 * Developer:  	Sumiko Abe
 * Date: 	2006, Jan, 10
*********************************************************/

#include "dctk.h"
#include "ofstdinc.h"
#include <limits.h>

#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h> 
#include "mgh.h"

#define DEBUG 0 
#define TEMP_STRING_LENGTH 256

int writeMgh(char *mghFileName,  mgh_header *mghInfor, bool endianInfor);
int dcm2mgh(char *dcmListFileName, char *outputMghFilePath, char *apiTxtFileName, char *createdTxtFileName,char *patientID, char *UCSDTrackingNum);


int main(int argc, char *argv[])
{
	if(argc != 7)
        {
                printf("\t#####Usage: %s <DICOM LIST> <OUTPUT_MGH_PATH> <API_TXT_FILE> <CREATED_TXT_FILE> <PATIENT ID> <UCSD TRACKING NUM>\n", argv[0]);
                exit(0);
        }
	
	int ret = dcm2mgh(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);

	if( ret == 1 )
		return 1;
	else	
		return 11;
}

