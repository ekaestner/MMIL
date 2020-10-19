/********************************************************
 * Copyright (C) 2006 University of California, San Diego
 * ADNI project program: Image Converter
 * File Name:   mgh.h 
 * Developer:   Sumiko Abe
 * Date:        2006, Jan, 12
*********************************************************/
#include <sys/types.h>

#ifndef _MGH_HEADER_
#define _MGH_HEADER_

#define UNUSED_SPACE_SIZE 256
#define USED_SPACE_SIZE  3*4+4*3*4
#define MRI_UCHAR   0 
#define MRI_INT     1 
#define MRI_LONG    2 
#define MRI_FLOAT   3 
#define MRI_SHORT   4 
#define MRI_BITMAP  5 
#define TEMP_STRING_LENGTH 2048

typedef struct mgh_header
{
	int v;
	int width;
	int height;
	int depth;
        int nframes;
	int type;
	int dof;
	short ras_good_flag;

	float xyz_size[3];
	float x_ras[3];
	float y_ras[3];
	float z_ras[3];
	float c_ras[3];

	unsigned char *charImageArray;
	short *shortImageArray;
	int *intImageArray;
	float *floatImageArray;
	long *longImageArray;

	float xFOV;
	float yFOV;
	
} mgh_header;

typedef struct dicomInfo
{
	float SliceThickness;
	float PixelSpacing[2];
	float firstPosition[3];
	float secondPosition[3];
	float lastPosition[3];
	float imageOrientationPation[6];
	int   scanType; // 0: axial, 1: coronal, 2: saggital
	int   columns;
	int   rows;
	int   slices;
	int   frames;
	float *vol;
} dicomInfo;

typedef struct Mvox2mgh
{
	float Mmatrix[16];
	int nx;
	int ny;
	int nz;
	int nt;
	float *vol;
} Mvox2ras;


#endif /* _MGHHEADER_ */
