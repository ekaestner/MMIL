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


int read_mgh(char *mghFileName, mgh_header *mghInfor, int endianInfor);

//swap functions: little endian ==> big endian or big endian ==> little endian 
void swap( char* p1, char* p2 ); // short type swap
void bswap( short& n );		 // short type swap
void bswap( int& n );		 // integer type swap
void bswap( long& n );		 // float type swap
void bswap( float& n );		 // float type swap
void bswap( double& n );	 // double type swap

int mgh_bswap( mgh_header *mghHeader );

int read_mgh(char *mghFileName,  mgh_header *mghInfor, bool endianInfor)
{

	int i, j, k;
	int tmp_int[7];
	short tmp_short;

	// check the endian type of local machine       
	int n=1;
	int endianFlag=-1;      // big endian: 1; little endian: 0; initial: -1
	if( *(char*)&n == 1 )
	{
		endianFlag=0;
		//cout << "\t#####Little endian machine." << endl; 
	} 
	else
	{
		endianFlag=1;
		//cout << "Big endian machine." << endl;
	}

	FILE *fp;

	cout << "\n\t####Reading a MGH file 1: "<< mghFileName << endl ;

	if ( (fp = fopen(mghFileName, "r")) == NULL) {
                printf("Could not open %s, exiting ...\n", mghFileName);
                return (-1);
        }	

	fread( tmp_int, sizeof(int), 7, fp);
	if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
		for ( i=0; i<=6; i++)
		{
			if (endianInfor == 1)
				bswap(tmp_int[i]);
			// for debug
			//cout << "tmp_int["<<i<<"]:" <<tmp_int[i]<< endl;
		}

	fread(&tmp_short , sizeof(short), 1, fp);
	if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
		bswap(tmp_short);



	mghInfor->v = tmp_int[0];
	mghInfor->width = tmp_int[1];
	mghInfor->height = tmp_int[2];
	mghInfor->depth = tmp_int[3];
	mghInfor->nframes = tmp_int[4];
	mghInfor->type = tmp_int[5];
	mghInfor->dof = tmp_int[6];
	mghInfor->ras_good_flag = tmp_short;


	int unused_space_size = UNUSED_SPACE_SIZE-2;
	if(mghInfor->ras_good_flag)
	{
		unused_space_size = unused_space_size - USED_SPACE_SIZE;
		fread( mghInfor->xyz_size, sizeof(float), 3, fp);
		fread( mghInfor->x_ras, sizeof(float), 3, fp);
		fread( mghInfor->y_ras, sizeof(float), 3, fp);
		fread( mghInfor->z_ras, sizeof(float), 3, fp);
		fread( mghInfor->c_ras, sizeof(float), 3, fp);
		if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
		{

			for(i=0; i<3; i++)
			{
				
				bswap(mghInfor->xyz_size[i]);
				bswap(mghInfor->x_ras[i]);
				bswap(mghInfor->y_ras[i]);
				bswap(mghInfor->z_ras[i]);
				bswap(mghInfor->c_ras[i]);
			}
		}
	}

	fseek(fp, unused_space_size, SEEK_CUR) ;

	int nv = mghInfor->width * mghInfor->height * mghInfor->depth *mghInfor->nframes;

		if(mghInfor->type ==MRI_UCHAR)
		{
			//cout << "Type: unsigned char" << endl;
			mghInfor->charImageArray = (unsigned char *)calloc(nv, sizeof(unsigned char)); 
			fread(mghInfor->charImageArray, sizeof(unsigned char), nv, fp);
		}
		else if(mghInfor->type ==MRI_INT)
		{
			//cout << "Type: int" << endl;
			mghInfor->intImageArray = (int *)calloc(nv, sizeof(int)); 
			fread(mghInfor->intImageArray, sizeof(int), nv, fp);
			if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
				for(i=0; i<nv; i++)
					bswap(mghInfor->intImageArray[i]);
		}
		else if(mghInfor->type ==MRI_LONG)
		{
			//cout << "Type: long" << endl;
			mghInfor->longImageArray = (long *)calloc(nv, sizeof(long)); 
			fread(mghInfor->longImageArray, sizeof(long), nv, fp);
			if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
				for(i=0; i<nv; i++)
					bswap(mghInfor->longImageArray[i]);
		}
		else if(mghInfor->type ==MRI_FLOAT)
		{
			//cout << "Type: float" << endl;
			mghInfor->floatImageArray = (float *)calloc(nv, sizeof(float)); 
			fread(mghInfor->floatImageArray, sizeof(float), nv, fp);
			if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
				for(i=0; i<nv; i++)
					bswap(mghInfor->floatImageArray[i]);

		}
		else if(mghInfor->type ==MRI_SHORT)
		{
			//cout << "Type: short" << endl;
			mghInfor->shortImageArray = (short *)calloc(nv, sizeof(short)); 
			fread(mghInfor->shortImageArray, sizeof(short), nv, fp);
			if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
				for(i=0; i<nv; i++)
					bswap(mghInfor->shortImageArray[i]);
		}
		else
		{
			cout << "Type is not defined!, Exit!" << endl;
			return (-2);
		}

	fclose(fp);

	////////////////////////////////////////////////////////////////////
	////// little endian machine for reading a big endian data    //////
	////// or big endian machine for reading a littla endian data //////
	////////////////////////////////////////////////////////////////////
	/*
	{
		cout << "byte swapped!" << endl;
		mgh_bswap(mghInfor);
		getchar();
	}
	*/

/***** for debug *****

	cout << "v = "<<mghInfor->v << endl;
	cout << "width = " << mghInfor->width << endl;
	cout << "height = " << mghInfor->height << endl;
	cout << "depth = " << mghInfor->depth << endl;
	cout << "nframes = " << mghInfor->nframes << endl;
	cout << "type = " << mghInfor->type << endl;
	cout << "dof = " << mghInfor->dof << endl;

	getchar();

********************/
}

void swap( char* p1, char* p2 ){
  char t = *p1; *p1 = *p2; *p2 = t;
}

void bswap( short& n ){
  swap( 0+(char*)&n, 1+(char*)&n );
}

void bswap( int& n ){
  swap( 0+(char*)&n, 3+(char*)&n );
  swap( 1+(char*)&n, 2+(char*)&n );
}

void bswap( long& n ){
  swap( 0+(char*)&n, 3+(char*)&n );
  swap( 1+(char*)&n, 2+(char*)&n );
}

void bswap( float& n ){
  swap( 0+(char*)&n, 3+(char*)&n );
  swap( 1+(char*)&n, 2+(char*)&n );
}

void bswap( double& n ){
  swap( 0+(char*)&n, 7+(char*)&n );
  swap( 1+(char*)&n, 6+(char*)&n );
  swap( 2+(char*)&n, 5+(char*)&n );
  swap( 3+(char*)&n, 4+(char*)&n );
}


int mgh_bswap( mgh_header *mghHeader )
{
	int i, nv;
	
	nv = mghHeader->width * mghHeader->height * mghHeader->depth * mghHeader->nframes;
/*
	cout << "width : "<<  mghHeader->width<< endl;
	cout << "height : "<<  mghHeader->height<< endl;
	cout << "depth : "<<  mghHeader->depth<< endl;
	cout << "nframes : "<<  mghHeader->nframes<< endl;	
	cout << "type : "<<  mghHeader->type<< endl;	
	getchar();
*/
	if(mghHeader->type ==MRI_INT)
	{
		//cout << "int" << endl;
		for(i = 0; i < nv; i++)
			bswap(mghHeader->intImageArray[i]);
	}
	else if(mghHeader->type ==MRI_LONG)
	{
		//cout << "long" << endl;
		for(i = 0; i < nv; i++)
			bswap(mghHeader->longImageArray[i]);
	}
	else if(mghHeader->type ==MRI_FLOAT)
	{
		//cout << "float" << endl;
		for(i = 0; i < nv; i++)
		{
			
			bswap(mghHeader->floatImageArray[i]);	
		}
	}
	else if(mghHeader->type ==MRI_SHORT)
	{
		//cout << "short" << endl;
		for(i = 0; i < nv; i++)
			bswap(mghHeader->shortImageArray[i]);
	}
	else
	{
		//cout << "ERROR: Type is not defined, Exit!" << endl;
		return (-2);
	}

	if(mghHeader->ras_good_flag == 1)
		for(i=0; i<3; i++)
		{
			bswap(mghHeader->xyz_size[i]);
			bswap(mghHeader->x_ras[i]);
			bswap(mghHeader->y_ras[i]);
			bswap(mghHeader->z_ras[i]);
			bswap(mghHeader->c_ras[i]);
		}


	bswap(mghHeader->v);
	bswap(mghHeader->width);
	bswap(mghHeader->height);
	bswap(mghHeader->depth);
	bswap(mghHeader->nframes);
	bswap(mghHeader->type);
	bswap(mghHeader->dof);
	bswap(mghHeader->ras_good_flag);

	return(1);
}
