/********************************************************
 * Copyright (C) 2006 University of California, San Diego
 * ADNI project program: Image Converter: DICOM to MGH
 * File Name:	readmgh.cc 
 * Developer:  	Sumiko Abe
 * Date: 	2006, Jan, 12
*********************************************************/

#define INCLUDE_CSTDLIB
#define INCLUDE_CSTRING

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h> 

#include <limits.h>
#include "dctk.h"
#include "ofstdinc.h"
#include "mgh.h"


#define NRANSI
#define DEBUG 0
#define TEMP_STRING_LENGTH 256


int writeMgh(char *mghFileName, mgh_header *mghInfor, bool endianInfor);

/******************************************/
/*** swap functions: 			***/
/*** little endian ==> big endian 	***/
/*** big endian ==> little endian 	***/
/******************************************/
void swap( char* p1, char* p2 ); // short type swap
void bswap( short& n );		 // short type swap
void bswap( int& n );		 // integer type swap
void bswap( long& n );		 // float type swap
void bswap( float& n );		 // float type swap
void bswap( double& n );	 // double type swap

int mgh_bswap( mgh_header *mghHeader );

int writeMgh(char *mghFileName,  mgh_header *mghInfor, bool endianInfor)
{

	int i, j, k;
	int tmp_int[7];
	short tmp_short;
	int nv ;
	int v =	mghInfor->v;
	int width = mghInfor->width;
	int height = mghInfor->height;
	int depth = mghInfor->depth;
	int nframes = mghInfor->nframes;
	int type = mghInfor->type;
	int dog = mghInfor->dof;
	int ras_good_flag = mghInfor->ras_good_flag;

	// check the endian type of local machine	
	int n=1;
	int endianFlag=-1; 	// big endian: 1; little endian: 0; initial: -1
	if( *(char*)&n == 1 ) 
	{
		endianFlag=0;
		cout << "\t#####Little endian machine." << endl; 
	} 
	else 	
	{
		endianFlag=1;
		cout << "Big endian machine." << endl;
	}
	
	cout << "\t#####Writing a MGH file : "<< mghFileName << endl ;
	FILE *fp;
	if ( (fp = fopen(mghFileName, "w")) == NULL) {
                printf("Could not open %s, exiting ...\n", mghFileName);
                return (-1);
        }	


	//cout << "endianFlag =" << endianFlag << ", endianInfor: " << endianInfor << endl;	
	if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
		mgh_bswap(mghInfor);

/*
	if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
	{
		cout << "111" << endl;
		bswap(mghInfor->v);
		bswap(mghInfor->width);
		bswap(mghInfor->height);
		bswap(mghInfor->depth);
		bswap(mghInfor->nframes);
		bswap(mghInfor->type);
		bswap(mghInfor->dof);
		bswap(mghInfor->ras_good_flag);
	}
*/
	fwrite( &mghInfor->v, sizeof(int), 1, fp);
	fwrite( &mghInfor->width, sizeof(int), 1, fp);
	fwrite( &mghInfor->height, sizeof(int), 1, fp);
	fwrite( &mghInfor->depth, sizeof(int), 1, fp);
	fwrite( &mghInfor->nframes, sizeof(int), 1, fp);
	fwrite( &mghInfor->type, sizeof(int), 1, fp);
	fwrite( &mghInfor->dof, sizeof(int), 1, fp);

	int unused_space_size = UNUSED_SPACE_SIZE-2;
/*
	cout << "float size = " << sizeof(float) << endl;
	getchar();
*/
	if(ras_good_flag == 1)
	{
/*
	cout << "unused_space_size1 = " << unused_space_size << endl;
	cout << "USED_SPACE_SIZE = " << USED_SPACE_SIZE << endl;
	getchar();
*/
		int tmp_space_size = USED_SPACE_SIZE;
		unused_space_size = unused_space_size - tmp_space_size;
/*
	cout << "unused_space_size1 = " << unused_space_size << endl;
	cout << "USED_SPACE_SIZE = " << USED_SPACE_SIZE << endl;
	getchar();
*/
		fwrite( &mghInfor->ras_good_flag, sizeof(short), 1, fp);
		fwrite( mghInfor->xyz_size, sizeof(float), 3, fp);
		fwrite( mghInfor->x_ras, sizeof(float), 3, fp);
		fwrite( mghInfor->y_ras, sizeof(float), 3, fp);
		fwrite( mghInfor->z_ras, sizeof(float), 3, fp);
		fwrite( mghInfor->c_ras, sizeof(float), 3, fp);
	}
	else
		fwrite( 0, sizeof(short), 1, fp);
		//fwrite( &mghInfor->ras_good_flag, sizeof(short), 1, fp);

	nv = width * height * depth * nframes;
	int slices =  depth * nframes;
	int ns = width * height;

	char unused_space_char[unused_space_size];
/**
	cout << "unused_space_size2 = " << unused_space_size << endl;
	getchar();
**/

	for(i=0; i<unused_space_size; i++)
		unused_space_char[i] = 0;
	fwrite( unused_space_char, sizeof(char), unused_space_size, fp);
	
/*


	float *oneSliceImageArray[slices];

	for(i=0; i<depth*nframes; i++)
	{
		oneSliceImageArray[i] = (float *)calloc(ns, sizeof(float));
		for(j = 0; j< ns; j++)
			oneSliceImageArray[i][j] = mghInfor->floatImageArray[i*ns + j];
		fwrite(oneSliceImageArray[i], sizeof(float), ns, fp);
	}

	
***/
	if(type ==MRI_UCHAR)
	{
		cout << "\t#####Image Data Type: unsigned char" << endl;
		fwrite(mghInfor->charImageArray, sizeof(unsigned char), nv, fp);
	}
	else if(type ==MRI_INT)
	{
		cout << "\t#####Image Data Type: int" << endl;
		fwrite(mghInfor->intImageArray, sizeof(int), nv, fp);
	}
	else if(type ==MRI_LONG)
	{
		cout << "\t#####Image Data Type: long" << endl;
		fwrite(mghInfor->longImageArray, sizeof(long), nv, fp);
	}
	else if(type ==MRI_FLOAT)
	{
		cout << "\t#####Image Data Type: float" << endl;
		fwrite(mghInfor->floatImageArray, sizeof(float), nv, fp);
			
	}
	else if(type ==MRI_SHORT)
	{
		cout << "\t#####Image Data Type: short" << endl;
		fwrite(mghInfor->shortImageArray, sizeof(short), nv, fp);
	}
	else
	{
		
		cout << "ERROR: Type is not defined, Exit!" << endl;
		return (-2);
	}

	char emptychar = 0;
	for (i=0; i<16; i++)
		fwrite(&emptychar, sizeof(char), 1, fp);
		

	fclose(fp);


/* for debug*
	if ( (fp = fopen("test6.raw", "w")) == NULL) {
                printf("Could not open %s, exiting ...\n", mghFileName);
                return (-1);
        }	
	fseek(fp, 0, SEEK_END);
	fwrite(mghInfor->floatImageArray, sizeof(float), nv, fp);
	fclose(fp);

	
	cout << "in writting fucntion: "<< endl;
	cout << "\tv = "<<mghInfor->v << endl;
	cout << "\twidth = " << mghInfor->width << endl;
	cout << "\theight = " << mghInfor->height << endl;
	cout << "\tdepth = " << mghInfor->depth << endl;
	cout << "\tnframes = " << mghInfor->nframes << endl;
	cout << "\ttype = " << mghInfor->type << endl;
	cout << "\tdof = " << mghInfor->dof << endl;

	cout << "\tv = "<<v << endl;
	cout << "\twidth = " << width << endl;
	cout << "\theight = " << height << endl;
	cout << "\tdepth = " << depth << endl;
	cout << "\tnframes = " << nframes << endl;
	cout << "\ttype = " << type << endl;

	getchar();
**/
	return (0);
}

