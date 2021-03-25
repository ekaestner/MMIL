/********************************************************
 * Copyright (C) 2006 University of California, San Diego
 * ADNI project program: RAS coordinator system creation
 * File Name:   getMghCoordinator.cc 
 * Developer:   Sumiko Abe
 * Date:        2006, Apr, 14
*********************************************************/

#include "dctk.h"
#include "ofstdinc.h"
#include <limits.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include "mgh.h"
#include "loniApi.h"

#define DEBUG 0
#define DEBUG2 1 
#define TEMP_STRING_LENGTH 2048

dicomInfo*  read_ADNI_DICOMInfo(char **dcmOrderedFileName, int fnum, char *dcm2mghErrorLogFile);
mgh_header* convert_ADNI_DICOM2RAS(dicomInfo *dicomInfoStruct);
float* readImageFloatArray(char **orderedDICOMList, int nx, int ny, int nz);

float *compute_UA(float *U, float *A, int N, int M, int P);
float *compute_UT(float *U, int N, int M);
float normal_v(float *v, int n);
float dot_UV(float *U, float *V, int N);

dicomInfo*  read_ADNI_DICOMInfo(char **dcmOrderedFileName, int fnum, char *dcm2mghErrorLogFile)
{
	int i, j, k;
	float firstPosition[3], lastPosition[3];
	float secondPosition[3];
	float ImageOrientationPatient[6];
	float PixelSpacing;
	int columns, rows;
	const char *Manufacturer;
	dicomInfo *dicomInforStruct;
	float SliceThickness;	

	dicomInforStruct = (dicomInfo*)calloc(1, sizeof(dicomInforStruct));
	{
		DcmFileFormat * dfile = new DcmFileFormat;
		DcmItem * dataset = dfile->getDataset();
		DcmTagKey key;
		OFCondition cond = dfile->loadFile(dcmOrderedFileName[0]);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " <<dcmOrderedFileName[0]  << " loading error 0: " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s loading error : %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;

        	}
		
        	Uint16 rows = 0;
        	Uint16 columns = 0;
        	Sint16 slices = 0;

        	key.set(0x0028, 0x0010);        // number of row
        	cond = dataset->findAndGetUint16(key, rows, 0, 0);
        	if(cond.bad()) {
               		FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                        	printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not get image rows in DICOM file <"<< dcmOrderedFileName[0] << ">. Error follows: " << cond.text() << endl;
                	fprintf(err_fp,"0003, DICOM ==> MGH CONVERTER ERROR: Can not get image rows in DICOM file %s.  Error follows: %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;

        	}
		dicomInforStruct->rows 	= 	rows;

	        key.set(0x0028, 0x0011);        // number of column 
       	 	cond = dataset->findAndGetUint16(key, columns, 0, 0);
        	if(cond.bad()) {
               	 	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                        	printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not get image columns in DICOM file <"<< dcmOrderedFileName[0] << ">. Error follows: " << (char *)cond.text() << endl;
                	fprintf(err_fp,"0004, DICOM ==> MGH CONVERTER ERROR: Can not get image columns in DICOM file %s.  Error follows: %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
        	}
		dicomInforStruct->columns 	= 	columns;

		const char *charSliceThickness;	
		key.set(0x0018, 0x0050); 
		cond = 	dataset->findAndGetString(key, charSliceThickness);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " <<dcmOrderedFileName[0]  << " can not get SliceThickness: " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s can not get SliceThickness: %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
		SliceThickness = (float)atof(charSliceThickness);
		
		const char *charPixelSpacing;
		key.set(0x0028, 0x0030); 
		cond = 	dataset->findAndGetString(key, charPixelSpacing);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[0] << " can not get PixelSpacing: " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s can not get PixelSpacing: %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
		PixelSpacing = (float)atof(charPixelSpacing);

		key.set(0x0008, 0x0070); 
		cond = 	dataset->findAndGetString(key, Manufacturer);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[0] << "can not get Manufacturer : " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s can not get Manufacturer: %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}

		const char *firstImagePositionPatinet;	
		key.set(0x0020, 0x0032); 
		cond = 	dataset->findAndGetString(key, firstImagePositionPatinet);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[0] << "can not get ImagePositionPatinet : " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s can not get ImagePositionPatinet: %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
		char tempchar[3][64];
		for(i=0; i<3; i++)
			tempchar[i][0] = NULL;
		j = 0; k = 0;
		int ll = strlen(firstImagePositionPatinet);

		for(i=0; i < ll;  i++)
		{
			if(firstImagePositionPatinet[i] != '\\' )
			{
				tempchar[j][k] = firstImagePositionPatinet[i];
				k++;
			}
			else
			{
				tempchar[j][k] = '\0';
				j ++;
				k = 0;
			}
		}
		tempchar[2][k] = '\0';
		for(i = 0; i<3; i++)
			firstPosition[i] = (float)atof(tempchar[i]);

		const char *charImageOrientationPatient;
		key.set(0x0020, 0x0037); 
		cond = 	dataset->findAndGetString(key, charImageOrientationPatient);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[0] << " can not get ImageOrientationPatient : " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s can not get ImageOrientationPatinet: %s.", dcmOrderedFileName[0], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
#if DEBUG2
cout << "charImageOrientationPatient = " << charImageOrientationPatient<<endl;
getchar();
#endif

		char tempchar2[6][64];
		for(i=0; i<6; i++)
			tempchar2[i][0] = NULL;
		j = 0; k = 0;
		ll = strlen(charImageOrientationPatient);
		for(i=0; i < ll;  i++)
		{
			if ( charImageOrientationPatient[i] != '\\')
			{
				tempchar2[j][k] = charImageOrientationPatient[i];
				k++;
			}
			else
			{
				tempchar2[j][k] = '\0';
				j++;
				k = 0;
			}

		}
		tempchar2[5][k] = '\0';
		for(i = 0; i<6; i++)
			ImageOrientationPatient[i] = (float)atof(tempchar2[i]);
		
		
		
	}
	
	//read second dicom image
	{
		DcmFileFormat * dfile = new DcmFileFormat;
		DcmItem * dataset = dfile->getDataset();
		DcmTagKey key;
		OFCondition cond = dfile->loadFile(dcmOrderedFileName[1]);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[1] << " loading error : " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s loading error: %s.", dcmOrderedFileName[1], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
		const char *secondImagePositionPatient = NULL;	
		key.set(0x0020, 0x0032); 
		cond = dataset->findAndGetString(key, secondImagePositionPatient);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[1] << "can not get ImagePositionPatient : " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s can not get ImagePositionPatient: %s.", dcmOrderedFileName[1], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
		char tempchar[3][64];
		for(i=0; i<3; i++)
			tempchar[i][0] = NULL;
		j = 0; k = 0;
		int ll = strlen(secondImagePositionPatient);

		for(i=0; i < ll;  i++)
		{
			if(secondImagePositionPatient[i] != '\\' )
			{
				tempchar[j][k] = secondImagePositionPatient[i];
				k++;
			}
			else
			{
				tempchar[j][k] = '\0';
				j ++;
				k = 0;
			}
		}
		tempchar[2][k] = '\0';
		for(i = 0; i<3; i++)
			secondPosition[i] = (float)atof(tempchar[i]);
	}

	//read last dicom image
	{
		DcmFileFormat * dfile = new DcmFileFormat;
		DcmItem * dataset = dfile->getDataset();
		DcmTagKey key;
		OFCondition cond = dfile->loadFile(dcmOrderedFileName[fnum-1]);

		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile[fnum - 1]);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[fnum-1] << " loading error : " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s loading error: %s.", dcmOrderedFileName[fnum-1], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
		const char *lastImagePositionPatinet = NULL;	
		key.set(0x0020, 0x0032); 
		cond = dataset->findAndGetString(key, lastImagePositionPatinet);
		if(!cond.good()) {
                	FILE *err_fp;
                	if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                       		printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return dicomInforStruct;
                	}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << dcmOrderedFileName[fnum-1] << " can not get ImagePostionPatient : " << cond.text() << endl;
                	fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s can not get ImagePositionPatient: %s.", dcmOrderedFileName[fnum-1], cond.text());
                	fclose(err_fp);
			return dicomInforStruct;
		}
		char tempchar[3][128];	
		for(i=0; i<3; i++)
		tempchar[i][0] = NULL;
		j = 0; k = 0;
		int ll = strlen(lastImagePositionPatinet);

		for(i=0; i < ll;  i++)
		{
			if(lastImagePositionPatinet[i] != '\\' )
			{
				tempchar[j][k] = lastImagePositionPatinet[i];
				k++;
			}
			else
			{
				tempchar[j][k] = '\0';
				j ++;
				k = 0;
			}
		}
		tempchar[2][k] = '\0';
		for(i = 0; i<3; i++)
		{
			lastPosition[i] = (float)atof(tempchar[i]);
		}
	}

	dicomInforStruct->SliceThickness	=	SliceThickness;
	dicomInforStruct->PixelSpacing[0] 	= 	PixelSpacing;
	dicomInforStruct->PixelSpacing[1] 	= 	PixelSpacing;
	for(i=0; i<3; i++)
	{
		dicomInforStruct->firstPosition[i]	= 	firstPosition[i];
		dicomInforStruct->secondPosition[i]	= 	secondPosition[i];
		dicomInforStruct->lastPosition[i]	= 	lastPosition[i];
	}
	for(i=0; i<6; i++)
	{
		dicomInforStruct->ImageOrientationPatient[i] = ImageOrientationPatient[i];
	}

	

	dicomInforStruct->slices	= 	fnum;
	dicomInforStruct->frames	= 	1;

	return dicomInforStruct;
}

	

mgh_header* convert_ADNI_DICOM2RAS(dicomInfo *dicomInfoStruct)
{
	int i, j, k;

	mgh_header *outputMgh;   // struct subsequently used for .mgh file writing
	outputMgh = (mgh_header *)calloc (1, sizeof(mgh_header));

        float M[4][4];
        dicomInfo *pdI=dicomInfoStruct;

	int nx = dicomInfoStruct->columns;
	int ny = dicomInfoStruct->rows;
	int nz = dicomInfoStruct->slices;
	int nt = dicomInfoStruct->frames;
	int ns = nx * ny;
	int np = ns * nz;
	int nv = np * nt;

	outputMgh->floatImageArray = (float*) calloc(nv, sizeof(float));
	for(i = 0; i<nv; i++)
		outputMgh->floatImageArray[i] = dicomInfoStruct->vol[i];

	outputMgh->v = 1;
	outputMgh->width = dicomInfoStruct->columns;
	outputMgh->height = dicomInfoStruct->rows;
	outputMgh->depth = dicomInfoStruct->slices;
	outputMgh->nframes = dicomInfoStruct->frames;
        outputMgh->type = MRI_FLOAT;
        outputMgh->dof = 1;
        outputMgh->ras_good_flag = 1;

        // Initialize a vox2LPH transform matrix
        for(i=0; i<4; i++)
           for (j=0; j<4; j++)
               M[i][j]=0.0;
        M[3][3]=1.0;

#if DEBUG2
	cout <<"first slice position :\t";
	for(i=0; i<3; i++)
		cout <<"\t      "<<  dicomInfoStruct->firstPosition[i];
	cout <<endl;
	cout <<"last slice position :\t" ;
	for(i=0; i<3; i++)
		cout <<"\t      "<<  dicomInfoStruct->lastPosition[i];
	cout <<endl;
#endif

	for(i = 0; i<3; i++)
        {
		M[i][0] =  pdI->PixelSpacing[1]*pdI->ImageOrientationPatient[3+i];
                M[i][1] =  pdI->PixelSpacing[0]*pdI->ImageOrientationPatient[i];
        }

        for(i=0; i<3; i++)
                M[i][2] = (pdI->lastPosition[i] - pdI->firstPosition[i])/(float)(pdI->slices-1);

        for(i=0; i<3; i++)
                M[i][3] = pdI->firstPosition[i];


// Given M is vox2LPH transform matrix, the following forces slice order into the orientation given below.
        float standardSagittal[3] = {1,  0, 0};   // want R->L sag. slices
	float standardCoronal[3]  = {0, -1, 0};   // want P->A cor. slices
	float standardAxial[3]    = {0,  0, 1};   // want I->S axial slices
	float scalarProduct = 0;
	int sliceReverseFlag = 0;
	float sliceDirect[3];
	for(i=0; i<3; i++)
		sliceDirect[i] = M[i][2];

//        for(i=0; i<3; i++)
//                cout <<"\t      "<<  sliceDirect[i];
//        cout <<endl; 

	if ( fabsf((scalarProduct = dot_UV(sliceDirect, standardSagittal, 3) )) > 0.5 )
	{
		sliceReverseFlag = scalarProduct < 0 ? 1 : 0;
	}
	else if ( fabsf((scalarProduct = dot_UV(sliceDirect, standardCoronal, 3) )) > 0.5 )
	{
		sliceReverseFlag = scalarProduct < 0 ? 1 : 0;
	}
	else if ( fabsf((scalarProduct = dot_UV(sliceDirect, standardAxial, 3) )) > 0.5 )
	{
		sliceReverseFlag = scalarProduct < 0 ? 1 : 0;
	}
	else 
	{
		cout << "WARNING: Could not find sliceReverseFlag, set sliceReverseFlag = 0" << endl;
	}

// If sliceReverseFlag = 1, then 
//    a) reverse the order of slices in the volume 
//    b) multiply the inter-slice direction column of M by -1
//    c) use ImagePositionPatient of original LAST slice as the offset coords in M.

        outputMgh->floatImageArray = (float*) calloc(nv, sizeof(float));
        float *tempVolume = (float*) calloc(nv, sizeof(float));

        if( sliceReverseFlag == 1)
        {
                M[0][2] = -M[0][2];
                M[1][2] = -M[1][2];
                M[2][2] = -M[2][2];

                for(i=0; i<3; i++)
                    M[i][3] = pdI->lastPosition[i];
 
                for(i=0; i<nz; i++)
                        for(j=0; j<ns; j++)
                                tempVolume[i*ns+j] = dicomInfoStruct->vol[(nz-1-i)*ns + j];
        }
        else
        {
                for(i=0; i<nv; i++)
                        tempVolume[i] =  dicomInfoStruct->vol[i];
        }
        for(i=0; i<nz; i++)
                for(j=0; j<ny; j++)
                        for(k=0; k<nx; k++)
                                outputMgh->floatImageArray[i*ns+j*nx+k] = tempVolume[i*ns +k*nx+j];

        free(tempVolume);

// convert M from vox2LPH to vox2RAS (negate first two rows of M)
	for(i=0; i<4; i++)
        {
		M[0][i] = -M[0][i];
                M[1][i] = -M[1][i];
        }

// Populate the remaining parts of outputMgh in preparation for .mgh file writing
        outputMgh->xyz_size[0] = pdI->PixelSpacing[1];
        outputMgh->xyz_size[1] = pdI->PixelSpacing[0];

        float tmp_vector[3];
        for(i = 0; i<3; i++)
                tmp_vector[i] = M[i][2];
        outputMgh->xyz_size[2] = normal_v(tmp_vector, 3);

        for(i=0; i<3; i++)
                outputMgh->x_ras[i] = M[i][0]/outputMgh->xyz_size[0];
        for(i=0; i<3; i++)
                outputMgh->y_ras[i] = M[i][1]/outputMgh->xyz_size[1];
        for(i=0; i<3; i++)
                outputMgh->z_ras[i] = M[i][2]/outputMgh->xyz_size[2];

        float c[3];
        c[0] = ((float)outputMgh->width)/2.0;
        c[1] = ((float)outputMgh->height)/2.0;
        c[2] = ((float)outputMgh->depth)/2.0;

        outputMgh->c_ras[0] = M[0][3] + M[0][0]*c[0] + M[0][1]*c[1] + M[0][2]*c[2];
        outputMgh->c_ras[1] = M[1][3] + M[1][0]*c[0] + M[1][1]*c[1] + M[1][2]*c[2];
        outputMgh->c_ras[2] = M[2][3] + M[2][0]*c[0] + M[2][1]*c[1] + M[2][2]*c[2];
        
        return outputMgh;
}

//read volume data from ordered dicom files
float* readImageFloatArray(char **orderedDICOMList, int nx, int ny, int nz)
{

        int i, j, k;
        int ns, np;
        short *imageShortArray;
        float *imageFloatArray;
        short *oneSliceImageShortArray;
        char temp_str[TEMP_STRING_LENGTH];

        ns = nx * ny;
        np = ns * nz;

        oneSliceImageShortArray = (short *)calloc(ns, sizeof(short));
        imageShortArray = (short *)calloc(np, sizeof(short));
        imageFloatArray = (float *)calloc(np, sizeof(float));

        for(i=0; i<nz; i++)
        {
                FILE *dicomfp;
                if ( (dicomfp = fopen(orderedDICOMList[i], "r")) == NULL) {
                        printf("Could not open %s, exiting ...\n", orderedDICOMList[i]);
                        return NULL;
                }
		fseek(dicomfp, -ns*sizeof(short), SEEK_END);
		fread(oneSliceImageShortArray, sizeof(short), ns, dicomfp);
		fclose(dicomfp);

		for(j=0; j<ns; j++)
			imageShortArray[i*ns +j] = oneSliceImageShortArray[j];

        }
        for(i=0; i<np; i++)
		imageFloatArray[i] = (float)imageShortArray[i];
/*
cout << "nx = " << nx <<", ny = " << ny << endl;
FILE *fpp;
fpp = fopen("./testnew.raw", "w");
fwrite(imageFloatArray, nx*ny*nz, sizeof(float), fpp);
fclose(fpp);
cout << "wrote!" << endl;
getchar();
*/

	free(oneSliceImageShortArray);
	free(imageShortArray);
        return (imageFloatArray);
}

