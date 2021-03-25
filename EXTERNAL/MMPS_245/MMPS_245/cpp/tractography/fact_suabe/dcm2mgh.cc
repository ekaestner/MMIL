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
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h> 
#include "mgh.h"
#include "loniApi.h"

#define DEBUG  0 
#define DEBUG2 0 
#define TEMP_STRING_LENGTH 2048

int writeMgh(char *mghFileName,  mgh_header *mghInfor, bool endianInfor);
int dcm2mgh(char *dcmListFileName, char *outputMghFilePath, char *apiTxtFileName, char *createdTxtFileName, char *patientID, char *UCSDTrackingNum);
char *createMghFileName(char *apiTxtFileName, char *createdTxtFileName, char *originalDicomDescription, char *outputMghFilePath,  
			char *errorFile, char *magnetStrength, int seriesNum, int B1_flag, int *scanType, char *seriesInstanseUID);
struct dicomInfo*  read_ADNI_DICOMInfo(char **dcmOrderedFileName, int fnum, char *dcm2mghErrorLogFile);
char **orderDICOMSlice(char *dcmListFileName, int fnum, char *dcm2mghErrorLogFile);
struct mgh_header* convert_ADNI_DICOM2RAS(struct dicomInfo *dicomInfoStruct);
float* readImageFloatArray(char **orderedDICOMList, int nx, int ny, int nz);


int dcm2mgh(char *dcmListFileName, char *outputMghFilePath, char *apiTxtFileName, char *createdTxtFileName, char *patientID, char *UCSDTrackingNum)
{
	int i, j, k, m, n, l, ll;
	int nx, ny, nz, nt;
	int np, ns, nv;

	FILE *fp1;
	FILE *fp2;


	if ( (fp1 = fopen(dcmListFileName, "r")) == NULL) {
		printf("Could not open %s, exiting ...\n", dcmListFileName);
		return(11);
	}
	/* Count the total number of files */
	int num_files = 0;
	char temp_str1[TEMP_STRING_LENGTH];
	while (fgets(temp_str1, TEMP_STRING_LENGTH, fp1) != NULL) 
		num_files++;
	printf("\t#####Processing %d files\n", num_files);

	char filename[TEMP_STRING_LENGTH];
        fseek(fp1, 0, SEEK_SET);		 /* Go back to the beginning of the file */
	fscanf(fp1, "%s\n", filename);		 /* Read first header name */
	fclose(fp1);

	char dcm2mghErrorLogFile[1024];
	sprintf(dcm2mghErrorLogFile, "%s/error.%s.%s.log\0", outputMghFilePath, patientID, UCSDTrackingNum);

        // verify a template dictionary is loaded
        if (!dcmDataDict.isDictionaryLoaded()) {
		FILE *err_fp;
		if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
			return(11);
		}
                cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM dictionary was not loaded!  Set DCMDICTPATH" << endl;
                fprintf(err_fp,"0001, DICOM ==> MGH CONVERTER ERROR: DICOM dictionary loading Error: can not find DICOM dictionary.");
                fclose(err_fp);
                return (11);
        }

	if ( (fp1 = fopen(dcmListFileName, "r")) == NULL) {
		printf("Could not open %s, exiting ...\n", dcmListFileName);
		return (11);
	}
	DcmFileFormat * dfile = new DcmFileFormat;
        DcmItem * dataset = dfile->getDataset();
	DcmTagKey key;

	//for debugging
	//cout << "filename : " << filename << endl;
       	//obtain direct link to the dataset

       	OFCondition cond = dfile->loadFile(filename);
       	if(!cond.good()) {
		FILE *err_fp;
		if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
			return(11);
		}
                cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << filename << " loading error : " << cond.text() << endl;
                fprintf(err_fp,"0002, DICOM ==> MGH CONVERTER ERROR: DICOM file %s loading error : %s.", filename, cond.text());
                fclose(err_fp);
              	return(11); 
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
       			 return(11); 
               }
               cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not get image rows in DICOM file <"<< filename << ">. Error follows: " << cond.text() << endl;
               fprintf(err_fp,"0003, DICOM ==> MGH CONVERTER ERROR: Can not get image rows in DICOM file %s.  Error follows: %s.", filename, cond.text());
               fclose(err_fp);
               return(11); 
        }

	key.set(0x0028, 0x0011);        // number of column 
        cond = dataset->findAndGetUint16(key, columns, 0, 0);
        if(cond.bad()) {
                FILE *err_fp;
                if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
                        printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
              		return(11); 
                }
               cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not get image columns in DICOM file <"<< filename << ">. Error follows: " << cond.text() << endl;
               fprintf(err_fp,"0004, DICOM ==> MGH CONVERTER ERROR: Can not get image columns in DICOM file %s.  Error follows: %s.", filename, cond.text());
               fclose(err_fp);
               return(11); 
        }

	key.set(0x0020,0x0011);     // series number
	Sint32 subSeriesNum = -1;
	cond = dataset->findAndGetSint32(key, subSeriesNum);
	if(cond.bad()) {
		FILE *err_fp;
		if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
			return(11);
		}
                cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not get image series number in DICOM file "<< filename << ". Error follows: " << cond.text() << endl;
                fprintf(err_fp,"0004, DICOM ==> MGH CONVERTER ERROR: Can not get image series number in DICOM file %s.  Error follows: %s.", filename, cond.text());
                fclose(err_fp);
		return(10);
	}

	const char * seriesDescription = NULL;
	key.set(0x0008,0x103e);     // series description 
	cond = dataset->findAndGetString(key, seriesDescription);
	if(cond.bad()) {
		FILE *err_fp;
		if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
			return(11);
		}
                cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not get image series description in DICOM file <"<< filename << ">. Error follows: " << cond.text() << endl;
                fprintf(err_fp,"0005, DICOM ==> MGH CONVERTER ERROR: Can not get image series description in DICOM file %s.  Error follows: %s.\n", filename, cond.text());
                fclose(err_fp);
		return (11);
	}

	const char *magnetStrength;
	float floatMagnetStrength;
	key.set(0x0018, 0x0087);        // Magnetic Field Strength 
	cond = dataset->findAndGetString(key, magnetStrength);
	if(cond.bad()) {
		cout << "\t#####ERROR: DICOM file %s couldn't get magnet strength. Error follows:" ;
		cout << cond.text();
		FILE *err_fp;
		if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
			return(10);
		}
		fprintf(err_fp,"0017, DICOM DB IMPORT ERROR: %s, can not get magnet strength. Error follows: %s", filename, cond.text());
		fclose(err_fp);
		return(11);
	}
	floatMagnetStrength = (float)atof(magnetStrength);
	if(floatMagnetStrength >= 1000.0)
		floatMagnetStrength = floatMagnetStrength/10000.0;

	const char *imageType;
	key.set(0x0008, 0x0008);        // Magnetic Field Strength 
	cond = dataset->findAndGetString(key, imageType);
	if(cond.bad()) {
		cout << "\t#####ERROR: DICOM file %s couldn't get image type. Error follows:" ;
		cout << cond.text();
		FILE *err_fp;
		if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
			return(10);
		}
		fprintf(err_fp,"0018, DICOM DB IMPORT ERROR: %s, can not get image type. Error follows: %s", filename, cond.text());
		fclose(err_fp);
		return(11);
	}
	const char *seriesInstanceUID;
	key.set(0x0020, 0x000e);        // series instanse UID
	cond = dataset->findAndGetString(key, seriesInstanceUID);
	if(cond.bad()) {
		cout << "\t#####ERROR: DICOM file %s couldn't get seriesInstanceUID. Error follows:" ;
		cout << cond.text();
		FILE *err_fp;
		if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
			return(10);
		}
		fprintf(err_fp,"0018, DICOM DB IMPORT ERROR: %s, can not get seriesInstnaceUID. Error follows: %s", filename, cond.text());
		fclose(err_fp);
		return(11);
	}
	fclose(fp1);

	floatMagnetStrength = floatMagnetStrength*10.0;
	floatMagnetStrength = ceil((double)floatMagnetStrength)/10.0;
	float diffMagnet1 = fabs(floatMagnetStrength - 1.5);
	float diffMagnet2 = fabs(floatMagnetStrength - 3.0);

	char charMagStrength[128];
	if(diffMagnet1 < diffMagnet2)
		sprintf(charMagStrength, "%.1f\0", 1.5);
	else
		sprintf(charMagStrength, "%.1f\0", 3.0);


	char* originalSeriesDescription;
	originalSeriesDescription = strdup((char*)seriesDescription);

	int B1_flag = 0;
	if( (strstr(imageType, "ORIGINAL\\PRIMARY\\M\\ND") != NULL) || (strstr(imageType, "ORIGINAL\\PRIMARY\\OTHER") != NULL))
		B1_flag = 1;
	else 
		B1_flag = 2;
		
	int scanType = 0;

	char *MghFileName = createMghFileName(apiTxtFileName, createdTxtFileName, originalSeriesDescription, 
					      outputMghFilePath, dcm2mghErrorLogFile, charMagStrength, 
					      subSeriesNum, B1_flag, &scanType, (char *)seriesInstanceUID);
	if (MghFileName == NULL)
		return(8);

	char outputMghFileName[512];
	sprintf(outputMghFileName, "%s/%s\0", outputMghFilePath, MghFileName);

	//get dicom image slice order by image instance number;
	char **dcmOrderedFileName;
	dcmOrderedFileName = orderDICOMSlice(dcmListFileName, num_files, dcm2mghErrorLogFile);

	//get dicom header informationr;
	struct dicomInfo *dicomInfoStruct;
	dicomInfoStruct = read_ADNI_DICOMInfo(dcmOrderedFileName, num_files, dcm2mghErrorLogFile);
	
	//get image volume data 
	dicomInfoStruct->vol = readImageFloatArray(dcmOrderedFileName, columns, rows, num_files);

#if DEBUG
	cout << "SliceThickness = 	" << dicomInfoStruct->SliceThickness << endl;
	cout << "PixelSpacing[0] = 	" << dicomInfoStruct->PixelSpacing[0] << endl;
	cout << "PixelSpacing[1] = 	" << dicomInfoStruct->PixelSpacing[1] << endl;
	cout << "firstPosition[1] = 	" << dicomInfoStruct->firstPosition[0] << endl;
	cout << "firstPosition[2] = 	" << dicomInfoStruct->firstPosition[1] << endl;
	cout << "firstPosition[3] = 	" << dicomInfoStruct->firstPosition[2] << endl;
	cout << "secondPosition[1] = 	" << dicomInfoStruct->secondPosition[0] << endl;
	cout << "secondPosition[2] = 	" << dicomInfoStruct->secondPosition[1] << endl;
	cout << "secondPosition[3] = 	" << dicomInfoStruct->secondPosition[2] << endl;
	cout << "lastPosition[1] = 	" << dicomInfoStruct->lastPosition[0] << endl;
	cout << "lastPosition[2] = 	" << dicomInfoStruct->lastPosition[1] << endl;
	cout << "lastPosition[3] = 	" << dicomInfoStruct->lastPosition[2] << endl;
	for(i= 0; i<6; i++)
		cout<< "imageOrientationPation = " <<dicomInfoStruct->imageOrientationPation[i] << endl;
	cout << "scanType = 	" << dicomInfoStruct->scanType << endl;
	cout << "columns = 	" << dicomInfoStruct->columns << endl;
	cout << "rows = 	" << dicomInfoStruct->rows << endl;
	cout << "slices = 	" << dicomInfoStruct->slices << endl;
	cout << "frames = 	" << dicomInfoStruct->frames << endl;
#endif
	

	//get Mvox2mgh matrix and image volum;
	struct mgh_header *outputMgh;
	outputMgh = convert_ADNI_DICOM2RAS(dicomInfoStruct);  // convert vox=>mgh
#if DEBUG2
	cout << "outputMgh.v = 	" <<outputMgh->v << endl;
	cout << "outputMgh.width = 	" <<outputMgh->width << endl;
	cout << "outputMgh.height = 	" <<outputMgh->height << endl;
	cout << "outputMgh.depth = 	" <<outputMgh->depth << endl;
	cout << "outputMgh.nframe = 	" <<outputMgh->nframes<< endl;
	cout << "outputMgh.type = 	" <<outputMgh->type << endl;
	cout << "outputMgh.dof = 	" <<outputMgh->dof << endl;
	cout << "outputMgh.ras_good_flag = 	" <<outputMgh->ras_good_flag << endl;
#endif
	
	writeMgh(outputMghFileName, outputMgh, 1);

	return (1);
}

char **orderDICOMSlice(char *dcmListFileName, int num_files, char *dcm2mghErrorLogFile)
{


	char **dcmOrderedFileName;
	char **tmpOrderedFileName;
	int i, j, k;

	dcmOrderedFileName = (char **)calloc(num_files, sizeof(char *));
	tmpOrderedFileName = (char **)calloc(num_files, sizeof(char *));
	for(i=0; i<num_files; i++)
	{
		dcmOrderedFileName[i] = (char *)calloc(4096, sizeof(char));
		tmpOrderedFileName[i] = (char *)calloc(4096, sizeof(char));
	}
		

	int *dcmFileOrder = (int *)calloc(num_files, sizeof(int));
	int *tempOrder = (int *)calloc(num_files, sizeof(int));
	int *fileNameLen = (int *)calloc(num_files, sizeof(int));
	FILE *fp1;


	if ( (fp1 = fopen(dcmListFileName, "r")) == NULL) {
		printf("Could not open %s, exiting ...\n", dcmListFileName);
		return (NULL);
	}

	for(i=0; i<num_files; i++)
	{
		char *cc;
		char temp_str = NULL;
		Sint32 instance_number = 0;

		k= 0; j = 0;
		int currentPosition = ftell(fp1);
		for(k = 0; temp_str != '\n' || feof(fp1) ; k++)
		{
			temp_str = fgetc(fp1);
			//cout << temp_str
			j++;
		}
		//getchar();
                fseek(fp1, currentPosition, SEEK_SET);

		cc = (char *)calloc(j+1, sizeof(char));
		
		for(k=0; k<j; k++)
		{
			fread(&temp_str, sizeof(char), 1, fp1);
			cc[k] = temp_str;
		}
		cc[j-1] = '\0';
		//dcmOrderedFileName[i] = strdup(cc);
		sprintf(dcmOrderedFileName[i], "%s\0", cc);

       		//obtain direct link to the dataset1
		DcmFileFormat * dfile1 = new DcmFileFormat;
       		DcmItem * dataset1 = dfile1->getDataset();
		DcmTagKey key1;
        	OFCondition cond1 = dfile1->loadFile(cc);
		if(!cond1.good()) 
		{
			FILE *err_fp;
			if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
				printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return(NULL);
			}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR: DICOM file " << cc << " loading error : " << cond1.text() << endl;
                	fprintf(err_fp,"0007, DICOM ==> MGH CONVERTER ERROR: DICOM file %s loading error : %s.", cc, cond1.text());
                	fclose(err_fp);
               		return (NULL);
       		}

		key1.set(0x0020, 0x0013); 	// number of instance
		cond1 = dataset1->findAndGetSint32(key1, instance_number, 0, 0);
		if(cond1.bad()) {
			FILE *err_fp;
			if ( (err_fp = fopen(dcm2mghErrorLogFile, "a")) == NULL) {
				printf("Could not open %s, exiting ...\n", dcm2mghErrorLogFile);
				return(NULL);
			}
                	cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not get image instance number in DICOM file <"<< cc << ">. Error follows: " << cond1.text() << endl;
                	fprintf(err_fp,"0008, DICOM ==> MGH CONVERTER ERROR: Can not get image instance number in DICOM file %s.  Error follows: %s.\n", cc, cond1.text());
                	fclose(err_fp);
			return (NULL);
			continue;
		}

		tempOrder[i] = instance_number;
		free(cc);
	}

	int temp_int;
	for(i = 0; i < num_files; i++)
	{
		for(j = i; j < num_files; j++)
		{
			if(i == j)
				continue;
			if(tempOrder[j] < tempOrder[i])
			{
				char *temp_char;

				temp_int = tempOrder[j] ;
				tempOrder[j] = tempOrder[i] ;
				tempOrder[i] = temp_int;

				temp_char = strdup(dcmOrderedFileName[j]);
				//strcpy(dcmOrderedFileName[j], dcmOrderedFileName[i]);
				sprintf(dcmOrderedFileName[j], "%s\0", dcmOrderedFileName[i]);
				sprintf(dcmOrderedFileName[i], "%s\0", temp_char);	
			}
				
		}
	}

	fclose(fp1);

	return(dcmOrderedFileName);

}

char *createMghFileName(char *apiTxtFileName, char *createdTxtFileName,
			char *originalDicomDescription, char *outputMghFilePath, 
			char *errorFile, char *magnetStrength, 
			int seriesNum, int B1_flag, int *scanType, char *seriesInstanceUID)
{
	int i, j, k;
	FILE *afp, *cfp;
	char *apiString;
	
	loni_api loniApi;
	created_api dicomApi;

#if DEBUG
	cout << "originalDicomDescription = " << originalDicomDescription << endl;
	cout << "b1_flag = " << B1_flag << endl;
#endif

	if ( (afp = fopen(apiTxtFileName, "r")) == NULL) {
		FILE *err_fp;
		if ( (err_fp = fopen(errorFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", errorFile);
			return NULL;
		}
               	cout << "\t#####DICOM ==> MGH CONVERTER ERROR : Can not open API text file <"<< apiTxtFileName << ">. " << endl;
               	fprintf(err_fp,"0009, DICOM ==> MGH CONVERTER ERROR: Can not open API text file <%s>. \n", apiTxtFileName);
               	fclose(err_fp);
		return NULL;
	}
	// remove space from original series description
	{
		char *p;
		p= strchr(originalDicomDescription, ' ');
		while( p != NULL)
		{
			*(p++) = '_' ;
			p= strchr(originalDicomDescription, ' ');
		}
	}
	// change upper case to lower case
	int strLen = strlen(originalDicomDescription);	
	char *tmp_descreption = (char *)calloc(strLen, sizeof(char));
	for(j = 0; j<strLen ; j++)
	{
		if(isupper(originalDicomDescription[j]) != 0)
			tmp_descreption[j] = tolower(originalDicomDescription[j]);
		else
			tmp_descreption[j] = originalDicomDescription[j];
	}
	char tempDicomApi[3][128];
	char tempLoniApi[7][128];
	sprintf((char*)originalDicomDescription, "%s", tmp_descreption);
	{
		char tmp_char[256];
		char *p, *q = ", ";
		int breakFlag = 0;
		while (fgets(tmp_char, TEMP_STRING_LENGTH, afp) != NULL && !feof(afp)) 
		{
			if(breakFlag == 1)
				break;
			int ii=0;
		
			/**/
			/**/

			while(tmp_char[ii] != '\n' && !feof(afp))
			{
				ii++;
			}
			if(tmp_char[ii-1] == '\n')
				tmp_char[ii-1] = '\0';

			// remove space from original series description
			/*
			{
				char *p;
				p= strchr(tmp_char, ' ');
				while( p != NULL)
				{
					*(p++) = '_' ;
					p= strchr(tmp_char, ' ');
				}
			}
			*/
			
			// change upper case to lower case
			int strLen = ii+1;	
			char *tmp_descreption2 ;
			tmp_descreption2= (char *)calloc(strLen, sizeof(char));
			
			for(j = 0; j<strLen ; j++)
			{
				/*
				if(isupper(tmp_char[j]) != 0)
					tmp_descreption2[j] = tolower(tmp_char[j]);
				else
				*/	
					tmp_descreption2[j] = tmp_char[j];
			}
	
#if DEBUG	
			cout << "tmp_char = " << tmp_char << endl;
			cout << "originalDicomDescription = " << originalDicomDescription << endl;
			cout << "tmp_char2 = " << tmp_descreption2 << endl;
			getchar();
#endif
					
			char *searchRst;
			apiString = strdup(tmp_descreption2);
			/*
			cout << "apiString = " << apiString ;
			getchar();
			*/
			{
				int jj = 0;
				for (p = strtok(apiString, q); p!=NULL; p=strtok(NULL, q))
				{
					sprintf(tempLoniApi[jj], "%s\0", p);
					//printf("jj = %d, TempLoniApi = %s\n", jj, tempLoniApi[jj]);
					jj++;
				}
				int txtReadSeriesNum = atoi(tempLoniApi[4]);
				/**
				printf("originalDicomDescription =	%s\n",  originalDicomDescription);
				printf("tempLoniApi[4] =	%s\n",  tempLoniApi[4]);
				printf("txtReadSeriesNum =	%d, seriesNum = %d\n",  txtReadSeriesNum, seriesNum);
				**/
				if( (strstr( tempLoniApi[3], seriesInstanceUID)) != NULL )
				{
					sprintf(loniApi.patientID, 		"%s\0", tempLoniApi[0]);
					sprintf(loniApi.loniExamID, 		"%s\0", tempLoniApi[1]);
					sprintf(loniApi.loniSeriesID, 		"%s\0", tempLoniApi[2]);
					sprintf(loniApi.loniScanDate, 		"%s\0", tempLoniApi[3]);
					sprintf(loniApi.loniScanDescription, 	"%s\0", tempLoniApi[4]);
					sprintf(loniApi.MR_strength, 		"%s\0", tempLoniApi[5]);
					sprintf(loniApi.coilType, 		"%s\0", tempLoniApi[6]);

#if DEBUG
					printf("lenthAPI = %d, lenthORI = %d \n", strlen(tempLoniApi[4]), strlen(originalDicomDescription));
					printf("loniApi.patientID = 		%s\n",  loniApi.patientID);
					printf("loniApi.loniExamID =		%s\n",  loniApi.loniExamID);
					printf("loniApi.loniSeriesID =		%s\n",  loniApi.loniSeriesID);
					printf("loniApi.loniScanDate =		%s\n",  loniApi.loniScanDate);
					printf("loniApi.loniScanDescription =	%s\n",  loniApi.loniScanDescription);
					printf("loniApi.MR_strength =		%s\n",  loniApi.MR_strength);
					printf("loniApi.coilType =%s\n",  loniApi.coilType);
					getchar();
#endif
					char tempchar[512];
					int nn;	
					/*
					for ( nn = 0; loniApi.loniExamID[nn] != '\0'; nn++)
					{
						tempchar[nn] = loniApi.loniExamID[nn];
					}
					sprintf(loniApi.loniExamID, "%s\0", tempchar);
					*/
					char tempchar2[512];
					nn = 0;
					int temp_flag = 0;
					while (!temp_flag)
					{
						if(loniApi.coilType[nn] == '\n')
						{
							loniApi.coilType[nn] = '\0';
							temp_flag = 1;
							break;
						}
						else
						{
							nn ++;
							continue;
						}
					}
		
				}
				/*
				else 
				{
					FILE *err_fp;
					if ( (err_fp = fopen(errorFile, "a")) == NULL) {
						printf("Could not open %s, exiting ...\n", errorFile);
						return NULL;
					}
                			cout<<"\t#####DICOM ==> MGH CONVERTER WARNING : Can not find matching feature for mgh file name : DICOM series: <"<<seriesNum << ">."<<endl;
                			fprintf(err_fp,"0008,  DICOM ==> MGH CONVERTER WARNING : Can not find matching feature pattern for mgh file name : DICOM series <%d>, %s.\n", 
						seriesNum, originalDicomDescription);
                			fclose(	err_fp);
					return NULL;
				}
				*/
			}
			
			searchRst = strstr(apiString, originalDicomDescription);		
			/**
			cout << "apiString	= " << apiString<< endl;		
			cout << "original	= " << originalDicomDescription<< endl;		
			printf("searchRst	= %s", searchRst);
			getchar(); 
			**/
		}
	}
	fclose(afp);

	char tmp_patientID[256];	
	for(j = 0; loniApi.patientID[j] != '\0'; j++)
	{
		if(islower(loniApi.patientID[j]) != 0)
			loniApi.patientID[j] = toupper(loniApi.patientID[j]);
	}

	char *mghFileName;
	char tmpMghFileName[256];
	//mghFileName = strdup(tmpMghFileName);
	//return(mghFileName);
	if( (strstr( loniApi.coilType, "B1-H")) != NULL )
	{
		if(B1_flag == 1)
			sprintf(tmpMghFileName, "B1-H_%s_%s_%s_mag.mgh\0", magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
		else
			sprintf(tmpMghFileName, "B1-H_%s_%s_%s_phi.mgh\0", magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
		mghFileName = strdup(tmpMghFileName);
#if DEBUG
		cout << "loniApi.coilType = " << loniApi.coilType << ", mghFileName = " << mghFileName << endl;
		cout << "loniApi.loniExamID = " << loniApi.loniExamID << ", mghFileName = " << mghFileName << endl;
		getchar();
#endif
		return(mghFileName);
	}
	else if( (strstr( loniApi.coilType, "B1-B")) != NULL )
	{
		if(B1_flag == 1)
			sprintf(tmpMghFileName, "B1-B_%s_%s_%s_mag.mgh\0", magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
		else
			sprintf(tmpMghFileName, "B1-B_%s_%s_%s_phi.mgh\0", magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
		mghFileName = strdup(tmpMghFileName);
#if DEBUG
		cout << "loniApi.coilType = " << loniApi.coilType << ", mghFileName = " << mghFileName << endl;
		cout << "loniApi.loniExamID = " << loniApi.loniExamID << ", mghFileName = " << mghFileName << endl;
		getchar();
#endif
		return(mghFileName);
	}
	else if( (strstr( loniApi.coilType, "MPR")) != NULL && (strstr( loniApi.coilType, "MPR-R")) == NULL && (strstr( loniApi.coilType, "MPR-C")) == NULL)
	{
		sprintf(tmpMghFileName, "MPR_%s_%s_%s.mgh\0", magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
		mghFileName = strdup(tmpMghFileName);
#if DEBUG
		cout << "loniApi.coilType = " << loniApi.coilType << ", mghFileName = " << mghFileName << endl;
		cout << "loniApi.loniExamID = " << loniApi.loniExamID << ", mghFileName = " << mghFileName << endl;
		getchar();
#endif
		return(mghFileName);
	}
	else if( (strstr( loniApi.coilType, "MPR-R")) != NULL)
	{
		sprintf(tmpMghFileName, "MPR-R_%s_%s_%s.mgh\0", magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
		mghFileName = strdup(tmpMghFileName);
#if DEBUG
		cout << "loniApi.coilType = " << loniApi.coilType << ", mghFileName = " << mghFileName << endl;
		cout << "loniApi.loniExamID = " << loniApi.loniExamID << ", mghFileName = " << mghFileName << endl;
		getchar();
#endif
		return(mghFileName);
	}
	else if( (strstr( loniApi.coilType, "MPR-C")) != NULL)
	{
		sprintf(tmpMghFileName, "MPR-C_%s_%s_%s.mgh\0", magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
		mghFileName = strdup(tmpMghFileName);
		(*scanType) = 1;
#if DEBUG
		cout << "MPR-C " <<endl;
		getchar();
		cout << "loniApi.coilType = " << loniApi.coilType << ", mghFileName = " << mghFileName << endl;
		cout << "loniApi.loniExamID = " << loniApi.loniExamID << ", mghFileName = " << mghFileName << endl;
		getchar();
#endif
		return(mghFileName);
	}
	else 
	{
		FILE *err_fp;
		if ( (err_fp = fopen(errorFile, "a")) == NULL) {
			printf("Could not open %s, exiting ...\n", errorFile);
			return NULL;
		}
       		cout<<"\t#####DICOM ==> MGH CONVERTER WARNING : Can not find matching feature for mgh file name : DICOM series: <"<<seriesNum << ">."<<endl;
       		fprintf(err_fp,"0008,  DICOM ==> MGH CONVERTER WARNING : Can not find matching feature pattern for mgh file name : DICOM series %d, %s.\n", 
				seriesNum, originalDicomDescription);
		fclose(	err_fp);
		return NULL;
	}

	/*
	cout << "tmpMghFileName = " << tmpMghFileName << endl;
	cout << "loniApi.loniExamID = " << loniApi.loniExamID << endl;
	cout << "loniApi.loniSeriesID = " << loniApi.loniSeriesID << endl;
	getchar();
	*/
	
	//sprintf(tmpMghFileName, "%s_%s%s%s.mgh\0", tmpMghFileName, magnetStrength, loniApi.loniExamID, loniApi.loniSeriesID);
	//mghFileName = strdup(tmpMghFileName);

	/*	
	printf("mghFileName = %s\n", mghFileName);
	*/
}
