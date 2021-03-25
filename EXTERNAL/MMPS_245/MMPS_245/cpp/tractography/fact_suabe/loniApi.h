/********************************************************
 * Copyright (C) 2006 University of California, San Diego
 * ADNI project program: Image Converter
 * File Name:   loniApi.h 
 * Developer:   Sumiko Abe
 * Date:        2006, Feb, 20 
*********************************************************/
#include <sys/types.h>

#ifndef _LONI_API_
#define _LONI_API_

typedef struct loni_api
{
	char patientID[128];
	char loniExamID[64];
	char loniSeriesID[64];
	char loniScanDate[64];
	char loniScanDescription[128];
	char MR_strength[64];
	char coilType[64];
} loni_api;

typedef struct created_api
{
	char dicomExamNum[128];
	char dicomSeriesNum[64];
	char dicomGradStrength[128];
	char dicomScanDescription[128];
} created_api;


#endif /* _MGHHEADER_ */
