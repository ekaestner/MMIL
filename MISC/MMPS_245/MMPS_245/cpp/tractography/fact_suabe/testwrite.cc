#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <time.h>

#include "fact.h"

#include "dctk.h"
#include "ofstdinc.h"
#include "mgh.h"
//#include "mri.h"

#define max(a,b) (((a)>(b))?(a):(b)) 
#define min(a,b) (((a)<(b))?(a):(b)) 
#define MAX_INPUTARGV_NUM 100 
#define INPUTARGV_TYPE_NUM 35


// todo: fiber color as option

// main entry
//================================================================================================
// Name:	main()
// Purpose: main entry of the program
//
// Usage:   Raw2Mgh  mgh_input_file raw_input_file output_mgh_file
//
//================================================================================================

int Fiber2Mgh(FACT_MAP stFactMap, mgh_header *mghinfo);
int read_mgh(char *mghFileName,  mgh_header *mghInfor, bool endianInfor);
int writeMgh(char *mghFileName, mgh_header *mghInfor, bool endianInfor);
int writeRaw(FACT_MAP stFactMap);
FACT_MAP getFactmap(int argc, char* argv[] );
int getFactPara(int argc, char* argv[] , FACT_PARA *pstFactPara);
int readFileList(char *RoiListFile, char RoiFiles[MAX_ROI_FILE_Nr][TEMP_STRING_LENGTH]);
int readRoiOperList(char *RoiOperListFile, int* RoiOperations);
char** createOptionArray( );
int  createOptionFlag(int argc,  char *argv[], optionStruct option_st[INPUTARGV_TYPE_NUM]);
float setDefaultStartFA();
float setDefaultStopFA();
float setDefaultTurnAngleDeg();
bool setDefaultByteSwap();
int setDefaultFiberLenMin();
bool setDefaultFlipVecX();
bool setDefaultFlipVecY();
bool setDefaultFlipVecZ();

int ReadMghImgFiles(FACT_MAP *pstFactMap);
int main(int argc, char* argv[]) 
{
	printf("test mgh writing function: write_mgh( ) \n");
	printf(" By Sumiko Abe, May 2007\n");
	//-----------------------------------------------------------
	// checking: command line parameter
	//-----------------------------------------------------------
	//if (argc < 13) { DispErrorMsg(DM_ERR_Usage); return DM_ERR_Usage; }
	if (argc != 4  )
	{ 
		DispErrorMsg(DM_ERR_Usage); 
		return DM_ERR_Usage; 
	}

	mgh_header mghIninfo; // = (mgh_header *)calloc(1, sizeof(mgh_header));
	mgh_header mghOutinfo; // = (mgh_header *)calloc(1, sizeof(mgh_header));
	char*	szFilePathIn = argv[1];	        // input file name
	char*	szFilePathOut = argv[3];	// output file name

	int ret = read_mgh(szFilePathIn,  &mghIninfo, 1);

	FILE *fp;

	fp = fopen ( argv[2], "r");
	if ( fp == NULL)
	{
		DispErrorMsg(DM_ERR_FileOpen); 
		return DM_ERR_Usage; 
	}
	fseek(fp, 0, SEEK_END);

	int filesize = ftell(fp);
	int ImageSize = filesize / sizeof(float);
	cout << "file length = " << filesize << endl;
	cout << "image size  = " << ImageSize << endl;
	getchar();

	

	rewind(fp);
	fread(mghIninfo.floatImageArray, sizeof(float),  ImageSize, fp ) ;
	for (int i=0; i < ImageSize; i++)
		if(mghIninfo.floatImageArray[i] != 0)
		{
			mghIninfo.floatImageArray[i] = 256.0;
		}

	int ns = mghIninfo.width * mghIninfo.height;
	int nv = ns * mghIninfo.depth;

	float *tempfloat = mghIninfo.floatImageArray;
	for (int i=0; i < mghIninfo.depth; i++)
		for(int j = 0; j < ns; j++)
			 mghIninfo.floatImageArray[(mghIninfo.depth-(i+1))*ns +j] = tempfloat[i*ns + j];
	

	fclose(fp);

	writeMgh(szFilePathOut, &mghIninfo, 1);

/**
	cout << "============================================== "<< endl;
	cout << "nx = " << stFactMap.pstFactPara->nImgWidth << endl;
	cout << "ny = " << stFactMap.pstFactPara->nImgHeight << endl;
	cout << "nz = " << stFactMap.pstFactPara->nImgSlices << endl;
	cout << "orientation = " << stFactMap.pstFactPara->nImgOrientation << endl;
	cout << "ImgSequencing = " <<stFactMap.pstFactPara->nImgSequencing << endl;
	cout << "FovWidth = " << stFactMap.pstFactPara->fFovWidth << endl;
	cout << "FovHeight = " << stFactMap.pstFactPara->fFovHeight << endl;
	cout << "pixel w = " << stFactMap.pstFactPara->fPixelSizeWidth << endl;
	cout << "pixel h = " << stFactMap.pstFactPara->fPixelSizeHeight << endl;
	cout << "pixel thick = " << stFactMap.pstFactPara->fSliceThickness << endl;
	cout << "byteswap = " << stFactMap.pstFactPara->bSwapBytes << endl;
	cout << "Flip x = " << stFactMap.pstFactPara->bFlipVecX << endl;
	cout << "Flip y = " << stFactMap.pstFactPara->bFlipVecY<< endl;
	cout << "Flip z = " << stFactMap.pstFactPara->bFlipVecZ<< endl;
	cout << "Start FA = " << stFactMap.pstFactPara->fStartFA<< endl;
	cout << "Stop FA = " <<stFactMap.pstFactPara->fStopFA<< endl;
	cout << "TurnDeg = " << stFactMap.pstFactPara->fTurnDeg<< endl;
	cout << "============================================== "<< endl;
	getchar();
**/

}

//================================================================================================
// Name:	FreeMemDpfPara()
// Purpose: release memory used by FACT_PARA structure
//
// Entry:	pstFactPara			-- pointer to the FACT_PARA structure
//
// Return:  N/A
//================================================================================================
void FreeMemFactPara(FACT_PARA*	pstFactPara)
{
	return;

	// release memory 
	if (pstFactPara->szImgVecFile != NULL) free(pstFactPara->szImgVecFile);
	if (pstFactPara->szImgAniFile != NULL) free(pstFactPara->szImgAniFile);

	if (pstFactPara->szFiberAllFile != NULL)	free(pstFactPara->szFiberAllFile);
	if (pstFactPara->szFiberSelFile != NULL)	free(pstFactPara->szFiberSelFile);
	if (pstFactPara->szFiberAllTxtFile != NULL) free(pstFactPara->szFiberAllTxtFile);
	if (pstFactPara->szFiberSelTxtFile != NULL) free(pstFactPara->szFiberSelTxtFile);
	if (pstFactPara->szFiberSelMghFile != NULL) free(pstFactPara->szFiberSelMghFile);

	for (int i=0; i<MAX_ROI_FILE_Nr; i++)
		if (pstFactPara->szBinRoiFile[i] != NULL) free(pstFactPara->szBinRoiFile[i]);
}

//================================================================================================
// Name:	FreeMemDtiMap()
// Purpose: release memory used by FACT_MAP structure
//
// Entry:	pstDtiMap			-- pointer to the FACT_MAP structure
//
// Return:  N/A
//================================================================================================
void FreeMemFactMap(FACT_MAP*	pstFactMap)
{
	// original images
	if (pstFactMap->pfAniImg != NULL) free(pstFactMap->pfAniImg);
	if (pstFactMap->pxyzVecImg != NULL) free(pstFactMap->pxyzVecImg);

	// fiber data
	long i;
	if (pstFactMap->pstFiberAry != NULL) {
		for (i=0; i<pstFactMap->lFiberNrTotal; i++) free(pstFactMap->pstFiberAry[i].pxyzChain);
		free(pstFactMap->pstFiberAry);
	}

	// fiber counter of each voxel
	if (pstFactMap->pnFiberCntAry != NULL) free(pstFactMap->pnFiberCntAry);

	// fiber data index of each voxel
	if (pstFactMap->pnFiberIdxAry != NULL) {
		long nImgBlkSize = pstFactMap->pstFactPara->nImgWidth *
						   pstFactMap->pstFactPara->nImgHeight *
						   pstFactMap->pstFactPara->nImgSlices;

		for (i=0; i<nImgBlkSize; i++) 
			if (pstFactMap->pnFiberIdxAry[i] != NULL) free(pstFactMap->pnFiberIdxAry[i]);

		free(pstFactMap->pnFiberIdxAry);
	}
}

//================================================================================================
// Name:	DispErrorMsg()
// Purpose: Display error message
//
// Entry:	DM_ERR_Code = DM_ERR_Unknow	-- error code, defined in ErrorCode.h
//			pPara0..pPara3=NULL			-- optional parameter pointers
//
// Return:  true, if successful, otherwise, return false;
//
// Note:	pPara0..pPara3 can be used for passing parameters into this function. for example,
//			we can pass nImgWidth and nImgHeight by:
//
//			DispErrorMsg(DM_ERR_ImgDim,  &nImgWidth, &nImgHeight)
//================================================================================================
bool DispErrorMsg(int DM_ERR_Code, void* pPara0, void* pPara1, void* pPara2, void* pPara3)
{
	switch(DM_ERR_Code) {
	case DM_ERR_Usage:		// parameter: none
		//printf(" Usage: fact paraFile\n");
		cout << " Program Usage: raw2mgh MGH_INPUT_FILE RAW_INPUT_FILE OUTPUT_FILENAME\n" << endl;
		cout << "\t\tCheck your input arguments, Exit!"<< endl;

		cout << endl;
		cout << endl;
		break;

	case DM_ERR_FileOpen:	// para0 = file name
		printf( "File open error: %s\n", pPara0);
		break;

	case DM_ERR_FileRead:	// para0 = file name [,para1=ErrMsg]
		if (!pPara1)
			printf( "File read error: %s\n", pPara0);
		else 
			printf( "File read error: %s %s\n", pPara0, pPara1);
		break;

	case DM_ERR_MemAlloc:	// para0 = file name
		printf( "Memory allocation error: %s\n", pPara0);
		break;

	default:
		printf(" Run time error\n");
		break;
	}
	return true;
}

// implementation file for Read Parameter files
//================================================================================================
// Name:	ReadLn()
// Purpose:	Read from text file, inFile, until newline or EOF, or buffer is full.
//
// Entry:	inFile			-- Pointer to FILE structure
//			pBuf			-- pointer to the input buffer
//			nBufLen > 0		-- length of the buffer
//
// Return:	number of bytes have read, the last byte of the buffer is set as \0. the file pointer
//			is updated
//================================================================================================
int ReadLn(FILE *inFile, char* pBuf, int nBufLen)
{
	char	cIn = '\0';
	int		nPos = 0;
	size_t	nRead = 0;

	while ((nPos < nBufLen) && ((nRead=fread(&cIn, 1, sizeof(char), inFile)) == 1) ) {
		pBuf[nPos++] = cIn;
		if (cIn == '\n') break;
    }
	pBuf[nPos] = 0;
	return nPos;
}

//================================================================================================
// Name:	ReadFactParaFile()
// Purpose: Parsing fiber tracking parameter file
//
// Entry:	szParaFileName	-- A string that is the path to the desired file
//			pstFactPara		-- pointer to the FACT_PARA structure defined in Fact.h
//
// Return:  DM_PROC_Ok, if successful, the FACT_PARA data structure is filled in.
//				otherwise, error code defined in Fact.h
//================================================================================================
int ReadFactParaFile(char* szParaFileName, FACT_PARA* pstFactPara)
{
	//-------------------------------------------
	// initializing FACT_PARA structure: pstDpfPara
	//-------------------------------------------
	pstFactPara->nImgWidth = -1;			// image dimension
	pstFactPara->nImgHeight = -1;
	pstFactPara->nImgSlices = -1;

	pstFactPara->nImgOrientation = 0;			// image dimension
	pstFactPara->nImgSequencing = 0;

	pstFactPara->fFovWidth = -1;			// FOV
	pstFactPara->fFovHeight = -1;
	pstFactPara->fPixelSizeWidth = -1;		// Pixel size
	pstFactPara->fPixelSizeHeight = -1;
	pstFactPara->fSliceThickness = -1;

	pstFactPara->szImgVecFile = NULL;		// input file name
	pstFactPara->szImgAniFile = NULL;
	pstFactPara->bSwapBytes = false;

	pstFactPara->bFlipVecX = false;			//flip x-y-z?
	pstFactPara->bFlipVecY = false;
	pstFactPara->bFlipVecZ = false;

	pstFactPara->fStartFA = -1;			// tracking threshold values
	pstFactPara->fStopFA = -1;
	pstFactPara->fTurnDeg = -1;
	pstFactPara->nFiberLenMin = -1;

	// output files
	pstFactPara->szFiberAllFile = NULL;
	pstFactPara->szFiberSelFile = NULL;
	pstFactPara->szFiberAllTxtFile = NULL;
	pstFactPara->szFiberSelTxtFile = NULL;
	pstFactPara->szFiberSelMghFile = NULL;

	// ROI fies
	for (int i=0; i<MAX_ROI_FILE_Nr; i++)	pstFactPara->szBinRoiFile[i] = NULL;
	//-------------------------------------------
	// Open parameter file, and read it
	//-------------------------------------------
	FILE*	inFile;
	if( (inFile  = fopen( szParaFileName, "rt" )) == NULL ) {
		DispErrorMsg(DM_ERR_FileOpen, (void*)szParaFileName);
		return	DM_ERR_FileOpen;
	}

	//-------------------------------------------
	// parsing the parameter file now
	//-------------------------------------------
	char	szLnBuf[512];
	
	// find "Begin" of the file
	char	szTmpBuf[512];
	bool	bBeginOfPara = false;					// iteration control flag
	while (!bBeginOfPara) {
		if (ReadLn(inFile, szLnBuf, 512) <= 0) {
			DispErrorMsg(DM_ERR_FileRead, (void*)szParaFileName, (void*)" fail to find the keyword: Begin");
			fclose(inFile);
			return	DM_ERR_FileRead;
		}

		// format scanning...
		sscanf(szLnBuf,"%s",szTmpBuf);			// get the 1st string of this line 
		if ((strcmp(szTmpBuf,"Begin:") == 0) || (strcmp(szTmpBuf,"Begin") == 0)) bBeginOfPara = true;
	}

	// go ahead
	int		nRoiFileCnt=0;
	int		nRoiOpCnt=0;
	bool	bEof = false;				// iteration loop controller
	while (!bEof) {
		if (ReadLn(inFile, szLnBuf, 512) <= 0) {
			DispErrorMsg(DM_ERR_FileRead, (void*)szParaFileName, (void*)" unexpected EOF");
			fclose(inFile);
			return	DM_ERR_FileRead;
		}

		// format scanning...
		sscanf(szLnBuf,"%s",szTmpBuf);			// get the 1st string of this line 
		// image dimensions
		if	(strcmp(szTmpBuf,"ImageWidth:") == 0)		sscanf(szLnBuf,"%s %d\n", szTmpBuf, &(pstFactPara->nImgWidth));		//ImageWidth
		else if (strcmp(szTmpBuf,"ImageHeight:") == 0)		sscanf(szLnBuf,"%s %d\n", szTmpBuf, &(pstFactPara->nImgHeight));		//ImageHeight
		else if (strcmp(szTmpBuf,"ImageSlices:") == 0)		sscanf(szLnBuf,"%s %d\n", szTmpBuf, &(pstFactPara->nImgSlices));		//ImageSlices

		// image orientation
		else if (strcmp(szTmpBuf,"ImageOrientation:") == 0)		sscanf(szLnBuf,"%s %d\n", szTmpBuf, &(pstFactPara->nImgOrientation));		//ImageOrientation: 0=coronal, 1=axial, 2=sagittal
		else if (strcmp(szTmpBuf,"ImageSequencing:") == 0)		sscanf(szLnBuf,"%s %d\n", szTmpBuf, &(pstFactPara->nImgSequencing));		//Image Sequencing:	0=normal, 1=flipped?

		// pixel size or FOV
		else if (strcmp(szTmpBuf,"PixelSize(X):") == 0)	sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fPixelSizeWidth));	//Pixel size(X)
		else if (strcmp(szTmpBuf,"PixelSize(Y):") == 0)	sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fPixelSizeHeight));	//Pixel Size(Y)
		else if (strcmp(szTmpBuf,"FieldOfView(X):") == 0)	sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fFovWidth));			//FieldOfView(X)
		else if (strcmp(szTmpBuf,"FieldOfView(Y):") == 0)	sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fFovHeight));			//FieldOfView(Y)
		else if (strcmp(szTmpBuf,"SliceThickness:") == 0)	sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fSliceThickness));	//SliceThickness

		// user-defined threshold for fiber tracking
		else if (strcmp(szTmpBuf,"TrackingStartFA:") == 0)	sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fStartFA));	//tracking start FA
		else if (strcmp(szTmpBuf,"TrackingEndFA:") == 0)	sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fStopFA));	//tracking stop FA
		else if (strcmp(szTmpBuf,"FiberTurningAngle:")==0) sscanf(szLnBuf,"%s %f\n", szTmpBuf, &(pstFactPara->fTurnDeg));	//tracking turning angle
		else if (strcmp(szTmpBuf,"MinimumFiberLen:")==0)   sscanf(szLnBuf,"%s %d\n", szTmpBuf, &(pstFactPara->nFiberLenMin));	//minimum fiber length
		
		// data manipulation
		else if (strcmp(szTmpBuf,"SwapBytes:") == 0) {
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);		//flag of swap bytes
			if (strcmp(szTmpBuf,"Yes") == 0) pstFactPara->bSwapBytes = true;
		}

		else if (strcmp(szTmpBuf,"FlipEigenVecX:") == 0) {
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);		//flip X-component?
			if (strcmp(szTmpBuf,"Yes") == 0) pstFactPara->bFlipVecX = true;
		}

		else if (strcmp(szTmpBuf,"FlipEigenVecY:") == 0) {
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);		//flip Y-component?
			if (strcmp(szTmpBuf,"Yes") == 0) pstFactPara->bFlipVecY = true;
		}
		
		else if (strcmp(szTmpBuf,"FlipEigenVecZ:") == 0) {
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);		//flip Z-component?
			if (strcmp(szTmpBuf,"Yes") == 0) pstFactPara->bFlipVecZ = true;
		}

		// input file name
		else if (strcmp(szTmpBuf,"InFileVecImg:") == 0) {		// Vector image file
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			pstFactPara->szImgVecFile  = (char*)malloc(strlen(szTmpBuf)+1);	//(strlen + '\0')
			strcpy(pstFactPara->szImgVecFile, szTmpBuf);
		}

		else if (strcmp(szTmpBuf,"InFileAniImg:") == 0) {		// Vector image file
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			pstFactPara->szImgAniFile  = (char*)malloc(strlen(szTmpBuf)+1);	//(strlen + '\0')
			strcpy(pstFactPara->szImgAniFile, szTmpBuf);
		}

		// output file name
		else if (strcmp(szTmpBuf,"OutFileFiberAll:") == 0) {		// output: all fibers
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			pstFactPara->szFiberAllFile  = (char*)malloc(strlen(szTmpBuf)+1);	//(strlen + '\0')
			strcpy(pstFactPara->szFiberAllFile, szTmpBuf);
		}
		else if (strcmp(szTmpBuf,"OutFileFiberSel:") == 0) {		// output: selected fibers
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			pstFactPara->szFiberSelFile  = (char*)malloc(strlen(szTmpBuf)+1);	//(strlen + '\0')
			strcpy(pstFactPara->szFiberSelFile, szTmpBuf);
		}
		else if (strcmp(szTmpBuf,"OutTxtFileFiberAll:") == 0) {		// text file output: all fibers
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			pstFactPara->szFiberAllTxtFile  = (char*)malloc(strlen(szTmpBuf)+1);	//(strlen + '\0')
			strcpy(pstFactPara->szFiberAllTxtFile, szTmpBuf);
		}
		else if (strcmp(szTmpBuf,"OutTxtFileFiberSel:") == 0) {		// text file output: selected fibers
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			pstFactPara->szFiberSelTxtFile  = (char*)malloc(strlen(szTmpBuf)+1);	//(strlen + '\0')
			strcpy(pstFactPara->szFiberSelTxtFile, szTmpBuf);
		}

		// binary ROI files
		else if (strcmp(szTmpBuf,"BinRoiFile:") == 0) {		// ROI file, for fiber selection
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			pstFactPara->szBinRoiFile[nRoiFileCnt]  = (char*)malloc(strlen(szTmpBuf)+1);	//(strlen + '\0')
			strcpy(pstFactPara->szBinRoiFile[nRoiFileCnt++], szTmpBuf);
		}

		else if (strcmp(szTmpBuf,"Operation:") == 0) {			// ROI operation by roi-file
			sscanf(szLnBuf,"%s %s\n", szTmpBuf, szTmpBuf);
			if (strcmp(szTmpBuf,"OR") == 0) pstFactPara->nBinRoiOp[nRoiOpCnt++] = ROI_OpOr;
			else if (strcmp(szTmpBuf,"AND") == 0) pstFactPara->nBinRoiOp[nRoiOpCnt++] = ROI_OpAnd;
			else if (strcmp(szTmpBuf,"NOT") == 0) pstFactPara->nBinRoiOp[nRoiOpCnt++] = ROI_OpNot;
			else if (strcmp(szTmpBuf,"CUT_OR") == 0) pstFactPara->nBinRoiOp[nRoiOpCnt++] = ROI_OpCutOr;
			else if (strcmp(szTmpBuf,"CUT_AND") == 0) pstFactPara->nBinRoiOp[nRoiOpCnt++] = ROI_OpCutAnd;
		}

		// end of parameters
		else if ((strcmp(szTmpBuf,"End:") == 0) || (strcmp(szTmpBuf,"End") == 0)) bEof = true;		// End
		else {
			//this line is meanless, skip it off
		}
	}
	fclose(inFile);				// ok, DPF parsing finished, close the file

	//--------------------------------------------------------------
	// Validating the parameters, are they all ready? 
	//--------------------------------------------------------------
	// REC file(s)
	if ( (pstFactPara->szImgVecFile == NULL) || (pstFactPara->szImgAniFile == NULL) ) {
		DispErrorMsg(DM_ERR_FileRead, (void*)szParaFileName, (void*)" fail to find vector or anisotropy image file name");
		return	DM_ERR_FileRead;
	}

	// image dimension
	if ((pstFactPara->nImgWidth < 0) || (pstFactPara->nImgHeight < 0) || (pstFactPara->nImgSlices < 0)) {
		DispErrorMsg(DM_ERR_FileRead, (void*)szParaFileName, (void*)" forget image dimensions (ImageWidth/Height/Slices)?");
		return DM_ERR_FileRead;
	}
	
	// FOV or Pixel Size
	if (((pstFactPara->fFovWidth < 0) && (pstFactPara->fPixelSizeWidth < 0)) || 
		((pstFactPara->fFovHeight < 0) && (pstFactPara->fPixelSizeHeight < 0)) ||
		( pstFactPara->fSliceThickness <= 0) ) {
		DispErrorMsg(DM_ERR_FileRead, (void*)szParaFileName, (void*)" forget FOV or pixel-size (FieldOfView/PixelSize/SliceThickness)?");
		return DM_ERR_FileRead;
	}

	// tracking threshold
	if ((pstFactPara->fStartFA < 0) || (pstFactPara->fStopFA < 0) || (pstFactPara->fTurnDeg < 0)) {
		DispErrorMsg(DM_ERR_FileRead, (void*)szParaFileName, (void*)" forget tracking threshold values (TrackingStartFA/StopFA/TractTurningAngle)?");
		return DM_ERR_FileRead;
	}
	
	//--------------------------------------------------------------
	// regulate parameters
	//--------------------------------------------------------------
	// FOV
	if ((pstFactPara->fFovWidth < 0) && (pstFactPara->fPixelSizeWidth > 0)) 	// your input is pixel size
		pstFactPara->fFovWidth = pstFactPara->nImgWidth * pstFactPara->fPixelSizeWidth;

	if ((pstFactPara->fFovHeight < 0) && (pstFactPara->fPixelSizeHeight > 0)) 	// your input is pixel size
		pstFactPara->fFovHeight = pstFactPara->nImgHeight * pstFactPara->fPixelSizeHeight;

	return DM_PROC_Ok;
}

// implementation file for ReadMghImgFiles.cpp
//================================================================================================
// Name:	ReadMghImgFiles()
// Purpose: Read Mgh image files (vector and FA map)
//
// Entry:	pstFactMap		-- pointer to the FACt_MAP structure that is defined in the fact.h
//
// Return:  DM_PROC_Ok, if successful. 	otherwise, error code defined in fact.h
//================================================================================================
int	ReadMghImgFiles(FACT_MAP*	pstFactMap)
{
	//------------------------------------------------------------
	// Check image input format type: if MGH format, the files
	// have been read and return OK. 
	//------------------------------------------------------------

	//------------------------------------------------------------
	// initialize DTI_MAP data structure: pstDtiMap
	//------------------------------------------------------------
	pstFactMap->pfAniImg = NULL;							// original FA map image pointer
	pstFactMap->pxyzVecImg = NULL;							// VEC image pointer

	pstFactMap->pstFiberAry = NULL;							// fiber data array
	pstFactMap->pnFiberIdxAry = NULL;						// fiber data index array
	pstFactMap->pnFiberCntAry = NULL;						// fiber counter in each voxel

	mgh_header mghFAinfo; // = (mgh_header *)calloc(1, sizeof(mgh_header));
	mgh_header mghV0info;  //= (mgh_header *)calloc(1, sizeof(mgh_header));
	char*	szFilePath = pstFactMap->pstFactPara->szMghAniFile;	// file name
	char*	szFilePath2 = pstFactMap->pstFactPara->szMghVecFile;	// file name
	int ret ;
	ret = read_mgh(szFilePath,  &mghFAinfo, 1);
	ret = read_mgh(szFilePath2, &mghV0info, 1);
	int i, j, k;

	if (mghFAinfo.x_ras[0] != 0 ) 				//axial ras matrix
	{
		pstFactMap->pstFactPara->nImgOrientation = 1;
	}
		
	else if ( mghFAinfo.x_ras[1] != 0 ) 			//coronal ras matrix
	{
		 pstFactMap->pstFactPara->nImgOrientation = 2;
		
	}
	else							//sagittal ras matrix
	{
		 pstFactMap->pstFactPara->nImgOrientation = 0;
	}

	//------------------------------------------------------------
	// get DWI parameters from DPF_PARA data structure
	//------------------------------------------------------------
	int		nImgWidth = pstFactMap->pstFactPara->nImgWidth = mghFAinfo.width;			// image dimensions
	int		nImgHeight = pstFactMap->pstFactPara->nImgHeight = mghFAinfo.height;
	int		nImgSlices = pstFactMap->pstFactPara->nImgSlices = mghFAinfo.depth;
	pstFactMap->pstFactPara->fFovWidth = mghFAinfo.xFOV;
	pstFactMap->pstFactPara->fFovHeight = mghFAinfo.yFOV;
	pstFactMap->pstFactPara->fPixelSizeWidth = mghFAinfo.xyz_size[0];
	pstFactMap->pstFactPara->fPixelSizeHeight = mghFAinfo.xyz_size[1];
	pstFactMap->pstFactPara->fSliceThickness = mghFAinfo.xyz_size[2];

	pstFactMap->pstFactPara->bFlipVecZ = 1;

	bool	bSwapBytes = pstFactMap->pstFactPara->bSwapBytes;		// swap byte? 
	
	bool	bFlipVecX = pstFactMap->pstFactPara->bFlipVecX;		// flip?
	bool	bFlipVecY = pstFactMap->pstFactPara->bFlipVecY;		// flip?
	bool	bFlipVecZ = pstFactMap->pstFactPara->bFlipVecZ;		// flip?

	//--------------------------------------------------------------------------------
	// read FA map
	//--------------------------------------------------------------------------------
	float*  pfAniImg, tmpAniImg;
	XYZ_TRIPLE*  pxyzVecImg , tmpVecImg;
	FILE* inFile;

	long	lwImgBlkSize = (long)nImgWidth * nImgHeight * nImgSlices;
	long	ns = (long)nImgWidth * nImgHeight;
	pxyzVecImg = pstFactMap->pxyzVecImg = (XYZ_TRIPLE*)malloc(sizeof(XYZ_TRIPLE) * lwImgBlkSize);
	pfAniImg = pstFactMap->pfAniImg = (float*)malloc(sizeof(float) * lwImgBlkSize);
	
	/*
	for ( i=0; i<lwImgBlkSize; i++)
	{
		pfAniImg[i] =  mghFAinfo.floatImageArray[i];
		pxyzVecImg[i].x = mghV0info.floatImageArray[i];
		pxyzVecImg[i].y =  mghV0info.floatImageArray[i+lwImgBlkSize];
		pxyzVecImg[i].z = mghV0info.floatImageArray[i+2*lwImgBlkSize];
	}
	*/
	

	for(i=lwImgBlkSize; i>0; i--)
	{
		int slices = i/ns+1;
		int pixels = i%ns;

		if(slices >= nImgSlices)
			continue;
		pfAniImg[i] = pstFactMap->pfAniImg[i] = mghFAinfo.floatImageArray[(nImgSlices-slices)*ns + pixels ];
		pxyzVecImg[i].x = pstFactMap->pxyzVecImg[i].x = mghV0info.floatImageArray[(nImgSlices-slices)*ns + pixels ];
		pxyzVecImg[i].y = pstFactMap->pxyzVecImg[i].y = mghV0info.floatImageArray[(nImgSlices-slices)*ns + pixels+lwImgBlkSize];
		pxyzVecImg[i].z = pstFactMap->pxyzVecImg[i].z = mghV0info.floatImageArray[(nImgSlices-slices)*ns + pixels+2*lwImgBlkSize];
	}

	

/*
	FILE *testfp;
	testfp = fopen("./testfa_2.raw","w");
	fwrite( mghFAinfo.floatImageArray, sizeof(float),  lwImgBlkSize, testfp);
	fclose(testfp); 
	testfp = fopen("./testfa_2-1.raw","w");
	fwrite( pfAniImg, sizeof(float),  lwImgBlkSize, testfp);
	fclose(testfp); 

	testfp = fopen("./testv0_2.raw","w");
	fwrite(pxyzVecImg, sizeof(XYZ_TRIPLE),  lwImgBlkSize,testfp);
	fclose(testfp); 
*/

	//------------------------------------------------------------
	// swap bytes if applicable
	//------------------------------------------------------------
	if (bSwapBytes == true) {
		SwapFloatImgBytes(pfAniImg, lwImgBlkSize);
		SwapFloatImgBytes((float*)pxyzVecImg, lwImgBlkSize*3);
	}


	//------------------------------------------------------------
	// flip (x,y,z) if applicable
	//------------------------------------------------------------
	if (bFlipVecX == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].x = -pxyzVecImg[i].x;
	if (bFlipVecY == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].y = -pxyzVecImg[i].y;
	if (bFlipVecZ == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].z = -pxyzVecImg[i].z;

	// adjusting the eigen-vector, considering the pixel size, so that the pixel
	// position can be regarded as interger-coordinate value (normalized pixel size as 1x1x1),
	// otherwise, the image pixel is acturally located in non-interger position that makes 
	// the intersection points of the straight line and cubic surface difficult to find out
	float fPixSizeW = pstFactMap->pstFactPara->fFovWidth / pstFactMap->pstFactPara->nImgWidth;
	float fPixSizeH = pstFactMap->pstFactPara->fFovHeight / pstFactMap->pstFactPara->nImgHeight;
	float fPixSizeD =  pstFactMap->pstFactPara->fSliceThickness;
	for (i=0; i<lwImgBlkSize; i++) {
		pxyzVecImg[i].x /= fPixSizeW;
		pxyzVecImg[i].y /= fPixSizeH;
		pxyzVecImg[i].z /= fPixSizeD;
	}


	return DM_PROC_Ok;
}

// implementation file for ReadImgFiles.cpp
//================================================================================================
// Name:	ReadImgFiles()
// Purpose: Read image files (vector and FA map)
//
// Entry:	pstFactMap		-- pointer to the FACt_MAP structure that is defined in the fact.h
//
// Return:  DM_PROC_Ok, if successful. 	otherwise, error code defined in fact.h
//================================================================================================
int	ReadImgFiles(FACT_MAP*	pstFactMap)
{
	//------------------------------------------------------------
	// Check image input format type: if MGH format, the files
	// have been read and return OK. 
	//------------------------------------------------------------
	if(pstFactMap->pstFactPara->mghFileInput[0] == 1 && pstFactMap->pstFactPara->mghFileInput[1] == 1 )
		return DM_PROC_Ok;

	//------------------------------------------------------------
	// initialize DTI_MAP data structure: pstDtiMap
	//------------------------------------------------------------
	pstFactMap->pfAniImg = NULL;							// original FA map image pointer
	pstFactMap->pxyzVecImg = NULL;							// VEC image pointer

	pstFactMap->pstFiberAry = NULL;							// fiber data array
	pstFactMap->pnFiberIdxAry = NULL;						// fiber data index array
	pstFactMap->pnFiberCntAry = NULL;						// fiber counter in each voxel

	//------------------------------------------------------------
	// get DWI parameters from DPF_PARA data structure
	//------------------------------------------------------------
	int		nImgWidth = pstFactMap->pstFactPara->nImgWidth;			// image dimensions
	int		nImgHeight = pstFactMap->pstFactPara->nImgHeight;
	int		nImgSlices = pstFactMap->pstFactPara->nImgSlices;

	bool	bSwapBytes = pstFactMap->pstFactPara->bSwapBytes;		// swap byte? 
	
	bool	bFlipVecX = pstFactMap->pstFactPara->bFlipVecX;		// flip?
	bool	bFlipVecY = pstFactMap->pstFactPara->bFlipVecY;		// flip?
	bool	bFlipVecZ = pstFactMap->pstFactPara->bFlipVecZ;		// flip?

	//--------------------------------------------------------------------------------
	// read FA map
	//--------------------------------------------------------------------------------
	char*	szFilePath = pstFactMap->pstFactPara->szImgAniFile;	// file name
	float*  pfAniImg;
	FILE* inFile;

	long	lwImgBlkSize = (long)nImgWidth * nImgHeight * nImgSlices;
	XYZ_TRIPLE*  pxyzVecImg ;

	// Open file
	printf("  ....read anisotropy image file....\n");
	if( (inFile = fopen( szFilePath, "rb" )) == NULL ) {
		DispErrorMsg(DM_ERR_FileOpen, (void*)szFilePath);
	return DM_ERR_FileOpen;
	}

	// is the file-size >= image-size?
	fseek( inFile, 0,  SEEK_END);
	if (ftell(inFile) < (long)(lwImgBlkSize * sizeof(float))) {
  		DispErrorMsg(DM_ERR_FileRead, (void*)szFilePath, (void*)": file size is smaller than expected");
		fclose(inFile);
		return DM_ERR_FileRead;
	}

	// allocate memory for input data
	pfAniImg = pstFactMap->pfAniImg = (float*)malloc(sizeof(float) * lwImgBlkSize);
	if (!pfAniImg) {
		DispErrorMsg(DM_ERR_MemAlloc, (void*)": anisotropy images");
		fclose(inFile);
		return	DM_ERR_MemAlloc;
	}
	memset(pfAniImg,0,sizeof(float)*lwImgBlkSize);	// clear memory

	// read file
	fseek( inFile, 0,  SEEK_SET);
	if ( fread(pfAniImg, sizeof(float),  lwImgBlkSize, inFile ) != (size_t)lwImgBlkSize) {
		DispErrorMsg(DM_ERR_FileRead, (void*)szFilePath);
		fclose(inFile);
		free(pfAniImg); pstFactMap->pfAniImg = NULL;
		return DM_ERR_FileRead;
	}

	//close file
	fclose(inFile);

/*
	FILE *testfp;
	testfp = fopen("./testfa.raw","w");
	fwrite(pfAniImg, sizeof(float),  lwImgBlkSize,testfp);
	fclose(testfp); 
*/
	

	
	//--------------------------------------------------------------------------------
	// read Vector map
	//--------------------------------------------------------------------------------
	szFilePath = pstFactMap->pstFactPara->szImgVecFile;	// file name
	printf("  ....read vector image file....\n");

	// Open file
	if( (inFile = fopen( szFilePath, "rb" )) == NULL ) {
		DispErrorMsg(DM_ERR_FileOpen, (void*)szFilePath);
		free(pfAniImg); pstFactMap->pfAniImg = NULL;
		return DM_ERR_FileOpen;
	}

	// is the file-size >= image-size?
	fseek( inFile, 0,  SEEK_END);
	if (ftell(inFile) < (long)(lwImgBlkSize * sizeof(XYZ_TRIPLE))) {
  		DispErrorMsg(DM_ERR_FileRead, (void*)szFilePath, (void*)": file size is smaller than expected");
		free(pfAniImg); pstFactMap->pfAniImg = NULL;
		fclose(inFile);
		return DM_ERR_FileRead;
	}

	// allocate memory for input data
	pxyzVecImg = pstFactMap->pxyzVecImg = (XYZ_TRIPLE*)malloc(sizeof(XYZ_TRIPLE) * lwImgBlkSize);
	if (!pxyzVecImg) {
		DispErrorMsg(DM_ERR_MemAlloc, (void*)" vector image");
		free(pfAniImg); pstFactMap->pfAniImg = NULL;
		fclose(inFile);
		return	DM_ERR_MemAlloc;
	}
	memset(pxyzVecImg,0,sizeof(XYZ_TRIPLE)*lwImgBlkSize);	// clear memory

	// read file
	fseek( inFile, 0,  SEEK_SET);
	if ( fread(pxyzVecImg, sizeof(XYZ_TRIPLE),  lwImgBlkSize, inFile ) != (size_t)lwImgBlkSize) {
		DispErrorMsg(DM_ERR_FileRead, (void*)szFilePath);
		fclose(inFile);
		free(pfAniImg); pstFactMap->pfAniImg = NULL;
		free(pxyzVecImg); pstFactMap->pxyzVecImg = NULL;
		return DM_ERR_FileRead;
	}

	//close file
	fclose(inFile);
/*
	testfp = fopen("./testv0.raw","w");
	fwrite(pxyzVecImg, sizeof(XYZ_TRIPLE),  lwImgBlkSize,testfp);
	fclose(testfp); 
*/

	//------------------------------------------------------------
	// swap bytes if applicable
	//------------------------------------------------------------
	if (bSwapBytes == true) {
		SwapFloatImgBytes(pfAniImg, lwImgBlkSize);
		SwapFloatImgBytes((float*)pxyzVecImg, lwImgBlkSize*3);
	}

	//------------------------------------------------------------
	// flip (x,y,z) if applicable
	//------------------------------------------------------------
	long i;
	if (bFlipVecX == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].x = -pxyzVecImg[i].x;
	if (bFlipVecY == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].y = -pxyzVecImg[i].y;
	if (bFlipVecZ == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].z = -pxyzVecImg[i].z;

	// adjusting the eigen-vector, considering the pixel size, so that the pixel
	// position can be regarded as interger-coordinate value (normalized pixel size as 1x1x1),
	// otherwise, the image pixel is acturally located in non-interger position that makes 
	// the intersection points of the straight line and cubic surface difficult to find out
	float fPixSizeW = pstFactMap->pstFactPara->fFovWidth / pstFactMap->pstFactPara->nImgWidth;
	float fPixSizeH = pstFactMap->pstFactPara->fFovHeight / pstFactMap->pstFactPara->nImgHeight;
	float fPixSizeD =  pstFactMap->pstFactPara->fSliceThickness;
	for (i=0; i<lwImgBlkSize; i++) {
		pxyzVecImg[i].x /= fPixSizeW;
		pxyzVecImg[i].y /= fPixSizeH;
		pxyzVecImg[i].z /= fPixSizeD;
	}

	return DM_PROC_Ok;
}

//================================================================================================
// Name:	SwapFloatImgBytes()
// Purpose: Swap float image bytes
//
// Entry:	pfImg			-- input float image pointer 
//			lwImgSize		-- image size = (Width x Height x Slices [x Blocks]
//
// Return:  n/a
//
// Note:	A tricky method, bitwise-exclusive-OR (^), is used to swap two bytes
//			suppose the two bytes are A and B,
//			then:	A^B^B results in A. 
//			and:	A^A^B results in B. 
//			so:		A^=B; A^=B; B^=A; will do the byte swap
//
//			the diagram below shows the status changing with these "XOR" operations
//
//			Operation   Status A         Status B
//				----	    A		         B					// initial status
//				A^=B;      A^B               B
//              B^=A       A^B               A (=B^(A^B))
//				A^=B        B (=(A^B)^A)     A					// finished
//================================================================================================
void SwapFloatImgBytes(float* pfImg, long lwImgSize)
{
	unsigned char*	pbImg = (unsigned char*) pfImg;
	long	lwImgSizeByte = (long)sizeof(float)*lwImgSize;

	for (long j=0; j<lwImgSizeByte; j+=4) {
		pbImg[j] ^= pbImg[j+3];			// swap byte 0 and 3
		pbImg[j+3] ^= pbImg[j];
		pbImg[j] ^= pbImg[j+3];

		pbImg[j+1] ^= pbImg[j+2];		// swap byte 1 and 2
		pbImg[j+2] ^= pbImg[j+1];
		pbImg[j+1] ^= pbImg[j+2];
	}
}


// Get coordinators of ROI pointers, 2007-01-25 by sumiko ucsd
roiInfo FindRoiMask(FACT_MAP* pstFactMap, unsigned char *ROI)
{
	int i, j, k;

	roiInfo RoiInfo;
	
	int nImgWidth = pstFactMap->pstFactPara->nImgWidth;
	int nImgHeight = pstFactMap->pstFactPara->nImgHeight;
	int nImgSlices = pstFactMap->pstFactPara->nImgSlices;
	
	int nv = nImgWidth * nImgHeight * nImgSlices;

	RoiInfo.ROI_size = 0;	
	for (i = 0; i<nv; i++)
		if(ROI[i] == 1)
			RoiInfo.ROI_size ++;

	RoiInfo.ROI_points = (int *)calloc(nv, sizeof(int));
	j = 0;
	for(i = 0; i<nv; i++)
	{
		if(ROI[i] == 1)
		{
			RoiInfo.ROI_points[j] = i;
			j++;
		}
	}
	FILE *testFile;  
	short* testroi = (short *)calloc(nv, sizeof(short));
	for(i = 0; i< RoiInfo.ROI_size; i++)
		testroi[RoiInfo.ROI_points[i]] = 256;
	
	if( (testFile = fopen("/home/suabe/work/tractography/morisrc/new/testROI.raw", "w") )== NULL)
	{
		printf("Can not open file, exit...");
		exit(0);
	} 
	fwrite(testroi, sizeof(short), nv, testFile);
	fclose(testFile);


	return RoiInfo;
}

// Tracking the pointers in ROI, 2007-01-25 by sumiko ucsd
int RoiFiberTrackingFACT(FACT_MAP*	pstFactMap, roiInfo RoiInfo)
{
	// get image files for FACT
	float*		pfImgFa = pstFactMap->pfAniImg;
	XYZ_TRIPLE*	pxyzImgVec = pstFactMap->pxyzVecImg;
	
	// get image dimensions for FACT
	int	nImgWidth = pstFactMap->pstFactPara->nImgWidth;
	int	nImgHeight = pstFactMap->pstFactPara->nImgHeight;
	int	nImgSlices = pstFactMap->pstFactPara->nImgSlices;
	int 	ns = nImgWidth * nImgHeight;
	int 	nv = ns * nImgSlices;

	// get tracking threshold for FACT
	float	fStartFA = pstFactMap->pstFactPara->fStartFA;
	float	fStopFA = pstFactMap->pstFactPara->fStopFA;
	float	fTurnDeg = pstFactMap->pstFactPara->fTurnDeg;
	int	nFiberLenMin = pstFactMap->pstFactPara->nFiberLenMin;

	// initialize fiber data array;
	int		nAryGrowBy=100;
	long	lFiberArySize = nAryGrowBy;
	long	lFiberCnt = 0;
	FIBER*	pstFiberAry = (FIBER*)malloc(sizeof(FIBER)*lFiberArySize);

	if (pstFiberAry == NULL) {
		DispErrorMsg(DM_ERR_MemAlloc, (void*)" fiber data array");
		return	DM_ERR_MemAlloc;
	}

	// go over all voxels of the image volumn
	long	lFiberLenMax = 0;
	double	dFiberLenMean = 0;
	double	dFiberLenSqrMean = 0;
	int i, j, k;
	for (int i=0; i<RoiInfo.ROI_size; i++) 
	{
		// skip off lower-FA point
		if (pfImgFa[RoiInfo.ROI_points[i]] < fStartFA)  continue;

		// some how, the eigen vector may equal to zero, although it is nonsense. but it does exist 
		// in some users synthetic dataset, so I added following line.
		if ((fabs(pxyzImgVec[RoiInfo.ROI_points[i]].x) < 1E-3) && 
			(fabs(pxyzImgVec[RoiInfo.ROI_points[i]].y) < 1E-3) && 
			(fabs(pxyzImgVec[RoiInfo.ROI_points[i]].z) < 1E-3) )  continue;
		// begin for fiber tracking
		float z_location = (float)(RoiInfo.ROI_points[i] / ns);
		float y_location = (float)(RoiInfo.ROI_points[i] % ns / nImgWidth);
		float x_location = (float)(RoiInfo.ROI_points[i] % ns % nImgWidth);
		XYZ_TRIPLE	xyzSeed = {x_location, y_location, z_location};
		int	nFiberLen;
		XYZ_TRIPLE* pxyzFiberChain = GetFiberChainFact(pfImgFa, pxyzImgVec, 
						nImgWidth, nImgHeight, nImgSlices,
						xyzSeed, fStopFA, fTurnDeg, nFiberLenMin, 
						&nFiberLen);

		if (pxyzFiberChain == NULL) continue;
		// set fiber data structure
		pstFiberAry[lFiberCnt].pxyzChain = pxyzFiberChain;
		pstFiberAry[lFiberCnt].nLength = nFiberLen;

		pstFiberAry[lFiberCnt].nSelStatus = 0;		// not selected = default value
		pstFiberAry[lFiberCnt].nSelCnt = 0;			// number of being selected by ROI
		pstFiberAry[lFiberCnt].xyzPtSeed = xyzSeed;	// seed point

		// update fiber counter, and incresing FiberArray size if necessary...
		lFiberCnt ++;
		if (lFiberCnt >= lFiberArySize) {
			lFiberArySize += nAryGrowBy;
			void* pTmp = realloc(pstFiberAry, sizeof(FIBER)*lFiberArySize);

			if (pTmp != NULL) pstFiberAry = (FIBER*)pTmp;
			else {
				DispErrorMsg(DM_ERR_MemAlloc, (void*)" expanding fiber data array");
				for (int i=0; i<lFiberCnt; i++) free(pstFiberAry[i].pxyzChain);
					free(pstFiberAry);
				return	DM_ERR_MemAlloc;
			}
		}
		// update fiber statistics
		lFiberLenMax = max(lFiberLenMax, nFiberLen);
		dFiberLenMean += nFiberLen;
		dFiberLenSqrMean += nFiberLen*nFiberLen;
	}

	// reform the fiber array
	if (lFiberCnt > 0) pstFiberAry = (FIBER*)realloc(pstFiberAry, sizeof(FIBER)*lFiberCnt);
	else {
		free(pstFiberAry);
		pstFiberAry = NULL;
	}

	// set the fiber array
	pstFiberAry = (FIBER*)realloc(pstFiberAry, sizeof(FIBER)*lFiberCnt);

	// fiber array statistics
	pstFactMap->pstFiberAry = pstFiberAry;
	pstFactMap->lFiberNrTotal = lFiberCnt;
	pstFactMap->lFiberLenMax = lFiberLenMax;
	if (lFiberCnt > 0) {
		pstFactMap->fFiberLenMean = (float)(dFiberLenMean/lFiberCnt);
		pstFactMap->fFiberLenStdDev = (float)(dFiberLenSqrMean/lFiberCnt - pstFactMap->fFiberLenMean*pstFactMap->fFiberLenMean);
	}
	return DM_PROC_Ok;
}


// implementation for fiber tracking
//=================================================================================
// Name:	FiberTrackingFACT()
// Purpose:	Fiber tracking for all voxels (brute-force) using FACT method, 
//
// Entry:	pstFactMap		-- data structure
//
// Return:  DM_PROC_Ok, if successful. 	and FACT_MAP structure members, pstFiberAry,
//			lFiberNr and lFiberLenMean, are ready. 
//			otherwise, error code defined in fact.h
//=================================================================================
int FiberTrackingFACT(FACT_MAP*	pstFactMap)
{
	// get image files for FACT
	float*		pfImgFa = pstFactMap->pfAniImg;
	XYZ_TRIPLE*	pxyzImgVec = pstFactMap->pxyzVecImg;

	// get image dimensions for FACT
	int	nImgWidth = pstFactMap->pstFactPara->nImgWidth;
	int	nImgHeight = pstFactMap->pstFactPara->nImgHeight;
	int	nImgSlices = pstFactMap->pstFactPara->nImgSlices;

	// get tracking threshold for FACT
	float	fStartFA = pstFactMap->pstFactPara->fStartFA;
	float	fStopFA = pstFactMap->pstFactPara->fStopFA;
	float	fTurnDeg = pstFactMap->pstFactPara->fTurnDeg;
	int		nFiberLenMin = pstFactMap->pstFactPara->nFiberLenMin;
/*
cout << "+++++++++++++++++++++++++++++" << endl;

	cout << nImgWidth << endl;
	cout << nImgHeight <<endl;
	cout << nImgSlices << endl;

	// get tracking threshold for FACT
	cout << fStartFA << endl;
	cout <<	fStopFA << endl; 
	cout <<	fTurnDeg << endl;
	cout << nFiberLenMin << endl;
cout << "+++++++++++++++++++++++++++++" << endl;

int nv = nImgWidth * nImgHeight * nImgSlices;
FILE *fp;
fp = fopen("./testfa3.raw", "w");
fwrite(pfImgFa, sizeof(float), nv, fp);
fclose(fp);
fp = fopen("./testv03.raw", "w");
fwrite(pxyzImgVec, sizeof(XYZ_TRIPLE), nv, fp);
fclose(fp);
*/
	// initialize fiber data array;
	int		nAryGrowBy=100;
	long	lFiberArySize = nAryGrowBy;
	long	lFiberCnt = 0;
	FIBER*	pstFiberAry = (FIBER*)malloc(sizeof(FIBER)*lFiberArySize);

	if (pstFiberAry == NULL) {
		DispErrorMsg(DM_ERR_MemAlloc, (void*)" fiber data array");
		return	DM_ERR_MemAlloc;
	}

	// go over all voxels of the image volumn
	long	lFiberLenMax = 0;
	double	dFiberLenMean = 0;
	double	dFiberLenSqrMean = 0;
	for (int nZ=1; nZ<nImgSlices-1; nZ++) {
		printf("\t\t.... tracking slice %d of %d ....\n", nZ, nImgSlices);
		long dwSliceOff = nZ * (long)nImgHeight*nImgWidth;

		for (int nY=1; nY<nImgHeight-1; nY++) {
		/*
		printf("ny ....%d\n", nY);
		getchar();
		*/
			long dwSliceLineOff = dwSliceOff + (nY * nImgWidth);

			for (int nX=1; nX<nImgWidth-1; nX++) {
				long	dwImgPixOff = dwSliceLineOff + nX;

				// skip off lower-FA point
				if (pfImgFa[dwImgPixOff] < fStartFA)  continue;

				// some how, the eigen vector may equal to zero, although it is nonsense. but it does exist 
				// in some users synthetic dataset, so I added following line.
				if ((fabs(pxyzImgVec[dwImgPixOff].x) < 1E-3) && 
					(fabs(pxyzImgVec[dwImgPixOff].y) < 1E-3) && 
					(fabs(pxyzImgVec[dwImgPixOff].z) < 1E-3) )  continue;

				// begin for fiber tracking
				XYZ_TRIPLE	xyzSeed = {(float)nX, (float)nY, (float)nZ};
				int			nFiberLen;

				XYZ_TRIPLE* pxyzFiberChain = GetFiberChainFact(pfImgFa, pxyzImgVec, 
										nImgWidth, nImgHeight, nImgSlices,
										xyzSeed, fStopFA, fTurnDeg, nFiberLenMin, 
										&nFiberLen);

				if (pxyzFiberChain == NULL) continue;

				// set fiber data structure
				pstFiberAry[lFiberCnt].pxyzChain = pxyzFiberChain;
				pstFiberAry[lFiberCnt].nLength = nFiberLen;

				pstFiberAry[lFiberCnt].nSelStatus = 0;		// not selected = default value
				pstFiberAry[lFiberCnt].nSelCnt = 0;			// number of being selected by ROI
				pstFiberAry[lFiberCnt].xyzPtSeed = xyzSeed;	// seed point

				// update fiber counter, and incresing FiberArray size if necessary...
				lFiberCnt ++;
				if (lFiberCnt >= lFiberArySize) {
					lFiberArySize += nAryGrowBy;
					void* pTmp = realloc(pstFiberAry, sizeof(FIBER)*lFiberArySize);

					if (pTmp != NULL) pstFiberAry = (FIBER*)pTmp;
					else {
						DispErrorMsg(DM_ERR_MemAlloc, (void*)" expanding fiber data array");
						for (int i=0; i<lFiberCnt; i++) free(pstFiberAry[i].pxyzChain);
						free(pstFiberAry);
						return	DM_ERR_MemAlloc;
					}
				}

				// update fiber statistics
				lFiberLenMax = max(lFiberLenMax, nFiberLen);
				dFiberLenMean += nFiberLen;
				dFiberLenSqrMean += nFiberLen*nFiberLen;

			}	// end of nX
		}	// end of nY
	}	// end of nZ

	// reform the fiber array
	if (lFiberCnt > 0) pstFiberAry = (FIBER*)realloc(pstFiberAry, sizeof(FIBER)*lFiberCnt);
	else {
		free(pstFiberAry);
		pstFiberAry = NULL;
	}

	// set the fiber array
	pstFiberAry = (FIBER*)realloc(pstFiberAry, sizeof(FIBER)*lFiberCnt);

	// fiber array statistics
	pstFactMap->pstFiberAry = pstFiberAry;
	pstFactMap->lFiberNrTotal = lFiberCnt;
	pstFactMap->lFiberLenMax = lFiberLenMax;
	if (lFiberCnt > 0) {
		pstFactMap->fFiberLenMean = (float)(dFiberLenMean/lFiberCnt);
		pstFactMap->fFiberLenStdDev = (float)(dFiberLenSqrMean/lFiberCnt - pstFactMap->fFiberLenMean*pstFactMap->fFiberLenMean);
	}
	return DM_PROC_Ok;
}

//=================================================================================
// Name: GetFiberChainFact()
// Purpose: do fiber tracking using FACT method
//
// Entry:	pfImgFA			-- pointer to FA image
//			pxyzImgVec		-- pointer to vector image
//			nImgWidth/Height/Slices		-- image dimensions
//			xyzSeedPt		-- start point for fiber tracking
//			fStopFA			-- FA threshold when tracking will stop
//			fTurnDeg		-- tract turning threshold when tracking stops
//			nFiberLenMin	-- minimum fiber chain length
//			pnFiberLen		-- pointer to the fiber length, a return value
//
// Return:  pointer to the fiber chain (could be NULL) 
//			and its length (could be 0), pnFiberLen;
//
//=================================================================================
XYZ_TRIPLE*	GetFiberChainFact(float* pfImgFA, XYZ_TRIPLE* pxyzImgVec, 
							  int nImgWidth, int nImgHeight, int nImgSlices,
							  XYZ_TRIPLE xyzSeedPt, 
							  float fStopFA,		// FA threshold where tracking stops
							  float fTurnDeg,		// fiber tract turning angle
							  int nFiberLenMin,		// minimum length of the fiber chain
							  int* pnFiberLen		// return value
							  )
{
	// forward tracking, get a half-fiber chain
	XYZ_TRIPLE xyzFiberPtF[4196];
	int	nFiberLenF = GetHalfChainFact(pfImgFA, pxyzImgVec,
									  nImgWidth, nImgHeight, nImgSlices,
									  xyzSeedPt,fStopFA, fTurnDeg, nFiberLenMin,
									  1, xyzFiberPtF);

	// backward tracking, get a half-fiber chain
	XYZ_TRIPLE	xyzFiberPtB[4196];	// backward tracking result
	int	nFiberLenB = GetHalfChainFact(pfImgFA, pxyzImgVec,
									  nImgWidth, nImgHeight, nImgSlices,
									  xyzSeedPt, fStopFA, fTurnDeg, nFiberLenMin,
									  -1, xyzFiberPtB);

	//-----------------------------------------------------------
	// combine them together
	//-----------------------------------------------------------
	// if whole fiber is too short, forget it.
	int nFiberPtLen = nFiberLenF + nFiberLenB -1;
	if (nFiberPtLen < nFiberLenMin) return NULL;

	// OK, combining two fiber segment together to a temprary array
	XYZ_TRIPLE* pxyzFiberPtWhole = (XYZ_TRIPLE*)malloc(sizeof(XYZ_TRIPLE)*nFiberPtLen);
	if (pxyzFiberPtWhole == NULL) {
		DispErrorMsg(DM_ERR_MemAlloc, (void*)" fiber tracking, combine two section together");
		return NULL;
	}

	int i, j;
	for (j=0, i=nFiberLenF-1; i>=0; j++, i--)
		pxyzFiberPtWhole[j] = xyzFiberPtF[i];

	for (i=1; i<nFiberLenB; j++, i++)
		pxyzFiberPtWhole[j] = xyzFiberPtB[i];

	//-----------------------------------------------------------
	// post processing: smoothing, interpolation and so on..
	//-----------------------------------------------------------

	//-------------------------------------------------------------------------
	// do nothing, return fiber data chain as uncooked one
	//-------------------------------------------------------------------------
	*pnFiberLen = nFiberPtLen;
	return pxyzFiberPtWhole;
}

//=================================================================================
// Name: GetHalfChainFact()
// Purpose: do one-dimensional fiber tracking using FACT method
//
// Entry:	pfImgFA			-- pointer to FA image
//			pxyzImgVec		-- pointer to vector image
//			nImgWidth/Height/Slices		-- image dimensions
//			xyzSeedPt		-- start point for fiber tracking
//			fStopFA			-- FA threshold when tracking will stop
//			fTurnDeg		-- tract turning threshold when tracking stops
//			nFiberLenMin	-- minimum fiber chain length
//			nTrackDir		-- tracking direction, 1=forward, -1=backward
//			pxyzFiberPt		-- fiber data array, return array
//
// Return:  number of fiber data points and the fiber data array, pnFiberLen;
//
//=================================================================================
float CosAlphaVec(XYZ_TRIPLE vecA, XYZ_TRIPLE vecB)
{
	float fAbsA = (float)sqrt(vecA.x*vecA.x + vecA.y*vecA.y + vecA.z*vecA.z);
	float fAbsB = (float)sqrt(vecB.x*vecB.x + vecB.y*vecB.y + vecB.z*vecB.z);
	return (vecA.x*vecB.x + vecA.y*vecB.y + vecA.z*vecB.z) / (fAbsA*fAbsB);
}

int GetHalfChainFact(float* pfImgFA, XYZ_TRIPLE* pxyzImgVec,
					 int nImgWidth, int nImgHeight, int nImgSlices,
					 XYZ_TRIPLE xyzSeedPt,float fStopFA,float fTurnDeg, int nFiberLenMin,
					 int nTrackDir, XYZ_TRIPLE xyzFiberPt[]		// return value
					 )
{
	//-----------------------------------------------------------
	// Now, preparing for fiber tracking
	//-----------------------------------------------------------
	const float fStopCosAngle = (float)cos(fTurnDeg*3.1415926/180);
	const float fTinyStep = 1E-3f;		// tracking step, the tiny one
	const int   nTinyStepMax = 5;		// the maximun number of non-intrupted tiny steps


	// pixel coordinate, and pixel offset
	int nPixelX = (int)xyzSeedPt.x;
	int nPixelY = (int)xyzSeedPt.y;
	int nPixelZ = (int)xyzSeedPt.z;

	long	dwCurPixelOff = (long)nPixelX +
			        (long)nPixelY * nImgWidth +
			        (long)nPixelZ * nImgWidth * nImgHeight;

	// eigen xyztor at current pixel
	XYZ_TRIPLE	xyzCurEigenVec = pxyzImgVec[dwCurPixelOff];
	if (nTrackDir <= -1) {		// backward tracking?
		xyzCurEigenVec.x = -xyzCurEigenVec.x;
		xyzCurEigenVec.y = -xyzCurEigenVec.y;
		xyzCurEigenVec.z = -xyzCurEigenVec.z;
	}

	// record seed point as the start point
	int		nFiberChainIdx = 0;								// fiber chain index

	xyzFiberPt[nFiberChainIdx].x = xyzSeedPt.x + 0.5f;
	xyzFiberPt[nFiberChainIdx].y = xyzSeedPt.y + 0.5f;
	xyzFiberPt[nFiberChainIdx].z = xyzSeedPt.z + 0.5f;

	// do tracking now
	bool	bTrackingOK = false;
 	int		nTinyStepNr = 0;		// number of continuous-small-steps 
	while (!bTrackingOK) {
		// some how, the exgin vector may equal to zero, although it is nonsense, it does exist in some
		// synthetic dataset, so I added following line.
		if ((fabs(xyzCurEigenVec.x) < 1E-3) && 
			(fabs(xyzCurEigenVec.y) < 1E-3) && 
			(fabs(xyzCurEigenVec.z) < 1E-3) )  {
			bTrackingOK = true;
			continue;
		}

		// get next fiber coordinate
		float fStepLen;
		XYZ_TRIPLE	xyzFiberPtNew = GetLineCubeInterPt(xyzFiberPt[nFiberChainIdx], xyzCurEigenVec, &fStepLen);
		if (fStepLen < fTinyStep) nTinyStepNr ++;
		else nTinyStepNr = 0;

		// cast float to integer coornate and verify its validaty
		nPixelX = (int)(xyzFiberPtNew.x + xyzCurEigenVec.x/1000);  // make sure the new point is the new one after (int)
		nPixelY = (int)(xyzFiberPtNew.y + xyzCurEigenVec.y/1000);
		nPixelZ = (int)(xyzFiberPtNew.z + xyzCurEigenVec.z/1000);

		if ( (nPixelX < 0) || (nPixelX >= nImgWidth) ||
			 (nPixelY < 0) || (nPixelY >= nImgHeight) ||
			 (nPixelZ < 0) || (nPixelZ >= nImgSlices) ) {

			bTrackingOK = true;
			continue;
		}

		xyzFiberPt[++nFiberChainIdx] = xyzFiberPtNew;
		dwCurPixelOff = (long)nPixelX +
				(long)nPixelY * nImgWidth +
				(long)nPixelZ * nImgWidth * nImgHeight;

//		printf( "FiberPt[%3d]=(%f, %f, %f), Vec[%3d,%3d,%2d]=(%f, %f, %f), step=%f, TinyNr=%d\n", 
//				nFiberChainIdx, xyzFiberPtNew.x, xyzFiberPtNew.y, xyzFiberPtNew.z, 
//				nPixelX, nPixelY, nPixelZ, 
//				pxyzImgVec[dwCurPixelOff].x, pxyzImgVec[dwCurPixelOff].y, pxyzImgVec[dwCurPixelOff].z,
//				fStepLen, nTinyStepNr);

		// test: is the inner product between two eigen-xyztors too smaller? or FA too small?
		XYZ_TRIPLE xyzNewVec = pxyzImgVec[dwCurPixelOff];
		float fCosAlphaVec = CosAlphaVec(xyzNewVec, xyzCurEigenVec);
		if ( (pfImgFA[dwCurPixelOff] < fStopFA) ||
			(fabs(fCosAlphaVec) < fStopCosAngle) ) {

			bTrackingOK = true;
			continue;
		}

		// OK, preparing for next iteration. 1st, update new eigen vector value, if necessary
		if (fCosAlphaVec < 0) {	
			xyzNewVec.x = -xyzNewVec.x; 
			xyzNewVec.y = -xyzNewVec.y; 
			xyzNewVec.z = -xyzNewVec.z; 
		}

		// adjust the direction if too many small steps occured
		if (nTinyStepNr > nTinyStepMax) {
			// method 1: the new vector is the averaging of xyzNewVec and xyzCurEigenVec, and 
			//			then tracking from xyzFiberPtNew
//			xyzNewVec.x = (xyzNewVec.x + xyzCurEigenVec.x)/2;
//			xyzNewVec.y = (xyzNewVec.y + xyzCurEigenVec.x)/2;
//			xyzNewVec.z = (xyzNewVec.z + xyzCurEigenVec.x)/2;

			// method 2: the new vector is xyzCurEigenVec, and then
			//			re-tracking from the last fiber point using modified vector
			xyzNewVec = xyzCurEigenVec;  nFiberChainIdx --;

			if ((fabs(xyzNewVec.x) > fabs(xyzNewVec.y)) && 
				(fabs(xyzNewVec.x) > fabs(xyzNewVec.z))) {
				xyzNewVec.y = xyzNewVec.z = 0;
			}
			else if ((fabs(xyzNewVec.y) > fabs(xyzNewVec.x)) && 
					 (fabs(xyzNewVec.y) > fabs(xyzNewVec.z))) {
				xyzNewVec.x = xyzNewVec.z = 0;
			}
			else if ((fabs(xyzNewVec.z) > fabs(xyzNewVec.x)) && 
					 (fabs(xyzNewVec.z) > fabs(xyzNewVec.y))) {
				xyzNewVec.x = xyzNewVec.y = 0;
			}
			nTinyStepNr = 0;
		}

		// update current eigen vector
		xyzCurEigenVec = xyzNewVec;

	} // end of while 

	return	++nFiberChainIdx;
}

//=================================================================================
// Name:	GetLineCubeInterPt()
// Purpose: find the insetction point of the straight line on the cube's surface 
//			from a point within the cube and direction of the line
//			called by GetHalfChainFact() for sub-pixel fiber tracking
//
// Entry:	vecPt0		-- start point of the straight line.
//			vecDir		-- direction number of the line;
//			pfStepLen	-- step length, return value
//
// Return:	a 3D-vector point that is the intersection point of the line and cube's surface
//
// Ref: Mathematic Handbook (Chinese) pp339-341
// Main Idea: the straight line is expressed as parameter equations. that is:
//		x=x0 + p*t;    y=y0 + q*t;    and  z=z0 + r*t;
//		where, (p,q,r) is the direction number of the line (or eigen vector in my case)
//			t the parameter (line segment length)
//		given a cube which enclose the point (x0,y0,z0), we can find parameter t that
//		makes the line just intersected with the nesrest-face of the cube.
//		we can imagine that the six faces of the cube are single variable
//		plains like: x=floor(x0), x=ceil(x0) and so on. So, we can find out
//		parameter t for every plain that makes the line touched with it. and
//		choose the minimum t as the result, since this t makes the line intersected
//		with the nearest-face of the cube.
//
//		of course, there are some tricks must be considered for this algorithm:
//		1) if the p,q or r is too small, forget it. since it will make t becomes
//			infinite
//		2) if the given point (x0,y0,z0) lies on the cube-face, we should be careful
//			to avoid the final result unchanged at all since the floor() or ceil()
//			could be the value of the point-component itself. in other words, the
//			the destination point (intersection point on the cube-face) should be
//			reasonable apart from the start point.
//
//			we have two method to deal with this case:
//			1: before using this method, always guarantee the start point is
//				located inside the cube by add a small step to the original start point
//				that is the method I'm going to use.
//			2: testing if the ceil() or floor() output is too close to the start point,
//				if it is so, make increment or decrement to the expected cubd-face before
//				you do others;
//				that is the method in GetLineCubeInterPt_General();
//	Note: May 26, 2005
//			using method 1, we take a risk that a small step (10% of the vector modules) 
//				towards into the inside of the pixel may left it if the fiber line is very 
//				short in this pixel. although it is not a big deal, the fiber-line will 
//				intercepted with the beyond one. anyway, it is not so good.
//			
//			by carefully choose the float to integer trunctation functions (ceil, floor, and int) 
//				it is possible to calculate the interception point without a shout step at all
//				please see my modified version. below
//========================================================================
XYZ_TRIPLE GetLineCubeInterPt(XYZ_TRIPLE xyzPt0, XYZ_TRIPLE xyzDir, float* pfStepLen)
{
	// declare some variables
	const float fVerySmall = 0.01f;
	double fT = 1e20;

	// initialize fT (from x-component)
	if (xyzDir.x > fVerySmall) fT = ((int)(xyzPt0.x+1) - xyzPt0.x)/xyzDir.x;
	else if (xyzDir.x < -fVerySmall) fT = (ceil(xyzPt0.x-1) - xyzPt0.x)/xyzDir.x;

	// update fT (from y-component)
	if (xyzDir.y > fVerySmall) fT = min(((int)(xyzPt0.y+1) - xyzPt0.y)/xyzDir.y, fT);
	else if (xyzDir.y < -fVerySmall) fT = min((ceil(xyzPt0.y-1) - xyzPt0.y)/xyzDir.y, fT);

	// update fT (from z-component)
	if (xyzDir.z > fVerySmall) fT = min(((int)(xyzPt0.z+1) - xyzPt0.z)/xyzDir.z, fT);
	else if (xyzDir.z < -fVerySmall) fT = min((ceil(xyzPt0.z-1) - xyzPt0.z)/xyzDir.z, fT);

	// return final result
	XYZ_TRIPLE xyzPtEnd;
	xyzPtEnd.x = (float)(xyzPt0.x + xyzDir.x*fT);
	xyzPtEnd.y = (float)(xyzPt0.y + xyzDir.y*fT);
	xyzPtEnd.z = (float)(xyzPt0.z + xyzDir.z*fT);

	*pfStepLen = (float)fT;
	return xyzPtEnd;
}

//=================================================================================
// Name:	CreateFiberIndexAry()
// Purpose:	As the name said. fiber index will be used for fiber selection
//
// Entry:	pstFactMap		-- data structure
//
// Return:  DM_PROC_Ok, if successful. 	and FACT_MAP structure members, pnFiberIdxAry,
//			and pnFiberCntAry are ready. 
//			otherwise, error code defined in fact.h
//=================================================================================
int CreateFiberIndexAry(FACT_MAP*	pstFactMap)
{
	// get image dimensions for FACT
	int	nImgWidth = pstFactMap->pstFactPara->nImgWidth;
	int	nImgHeight = pstFactMap->pstFactPara->nImgHeight;
	int	nImgSlices = pstFactMap->pstFactPara->nImgSlices;

	// get fiber data array;
	long	lFiberNrTotal = pstFactMap->lFiberNrTotal;
	FIBER*	pstFiberAry = pstFactMap->pstFiberAry;

	int i,j;
	//-----------------------------------------------------------------
	// create an fiber counter matrix, get fiber-# for each voxel,
	//-----------------------------------------------------------------
	// allocate memory for fiber counter of each voxel
	int* pnFiberCntAry = (int*)malloc(sizeof(int) * (long)nImgWidth*nImgHeight*nImgSlices);
	if (pnFiberCntAry == NULL) {
		DispErrorMsg(DM_ERR_MemAlloc, (void*)" fiber counter matrix");
		return	DM_ERR_MemAlloc;
	}
	memset(pnFiberCntAry, 0, sizeof(int) * (long)nImgWidth*nImgHeight*nImgSlices);

	// build Fiber-Count matrix
	for (i=0; i<lFiberNrTotal; i++) {	// for all fibers
		XYZ_TRIPLE* pxyzFiberChain = pstFiberAry[i].pxyzChain;
		int			nFiberLen = pstFiberAry[i].nLength;

		DWORD	dwPix0 = (DWORD)pxyzFiberChain[0].x +
							(DWORD)pxyzFiberChain[0].y*nImgWidth +
							(DWORD)pxyzFiberChain[0].z*nImgWidth*nImgHeight;

		pnFiberCntAry[dwPix0] ++;	// fiber count at this pixel +1

		DWORD	dwPix1;
		for (j=1; j<nFiberLen; j++) {		// for all remaining fiber point
			dwPix1 = (DWORD)pxyzFiberChain[j].x +
					   (DWORD)pxyzFiberChain[j].y*nImgWidth +
					   (DWORD)pxyzFiberChain[j].z*nImgWidth*nImgHeight;
			if (dwPix1 == dwPix0) continue;	// a new pixel position?
			dwPix0 = dwPix1;
			pnFiberCntAry[dwPix0] ++;	// fiber count at this pixel +1
		}
	}

	//-----------------------------------------------------------------
	// create an fiber data index array for each voxel,
	//-----------------------------------------------------------------
	// 1. allocate memory for fiber index array
	int** pnFiberIdxAry = (int**)malloc(sizeof(int*) * (long)nImgWidth*nImgHeight*nImgSlices);
	if (pnFiberCntAry == NULL) {
		DispErrorMsg(DM_ERR_MemAlloc, (void*)" fiber index array");
		free(pnFiberCntAry);
		return	DM_ERR_MemAlloc;
	}
	memset(pnFiberIdxAry, NULL, sizeof(int*) * (long)nImgWidth*nImgHeight*nImgSlices);

	// 2. allocate/initialize each element of the pnFiberIdxAry
	for (i=0; i<nImgWidth*nImgHeight*nImgSlices; i++)
		if (pnFiberCntAry[i] <= 0) pnFiberIdxAry[i] = NULL;
		else {
			pnFiberIdxAry[i] = (int*)malloc(sizeof(int) * pnFiberCntAry[i]);
			if (pnFiberIdxAry[i] == NULL) {
				DispErrorMsg(DM_ERR_MemAlloc, (void*)" element of the fiber index array");
				free(pnFiberCntAry);
				for (j=0; j<i; j++) free(pnFiberIdxAry[i]);
				free(pnFiberIdxAry);
				return	DM_ERR_MemAlloc;
			}
		}

	// 3. set fiber-index matrix
	memset(pnFiberCntAry, 0, sizeof(int) * (long)nImgWidth*nImgHeight*nImgSlices);
	for (i=0; i<lFiberNrTotal; i++) {	// for all fibers
		XYZ_TRIPLE* pxyzFiberChain = pstFiberAry[i].pxyzChain;
		int			nFiberLen = pstFiberAry[i].nLength;

		DWORD	dwPix0 = (DWORD)pxyzFiberChain[0].x +
							(DWORD)pxyzFiberChain[0].y*nImgWidth +
							(DWORD)pxyzFiberChain[0].z*nImgWidth*nImgHeight;

		int*	pnFiberIdx = pnFiberIdxAry[dwPix0];			// fiber-index of this pixel
		int		nFiberNrCur = pnFiberCntAry[dwPix0];		// current fiber-# of this pixel
		pnFiberIdx[nFiberNrCur] = i;							// add a fiber index to this fiber-index array

		pnFiberCntAry[dwPix0] ++;							// update current fiber-# of this pixel

		DWORD	dwPix1;
		for (j=1; j<nFiberLen; j++) {		// for all remaining fiber point
			dwPix1 = (DWORD)pxyzFiberChain[j].x +
					   (DWORD)pxyzFiberChain[j].y*nImgWidth +
					   (DWORD)pxyzFiberChain[j].z*nImgWidth*nImgHeight;
			if (dwPix1 == dwPix0) continue;	// a new pixel position?

			dwPix0 = dwPix1;

			pnFiberIdx = pnFiberIdxAry[dwPix0];		// fiber-index of this pixel
			nFiberNrCur = pnFiberCntAry[dwPix0];	// current fiber-# of this pixel
			pnFiberIdx[nFiberNrCur] = i;			// add index of fiber-record to this array

			pnFiberCntAry[dwPix0] ++;				// update current fiber-# of this pixel
		}
	}

	// assign fiber index array.
	pstFactMap->pnFiberIdxAry = pnFiberIdxAry;
	pstFactMap->pnFiberCntAry = pnFiberCntAry;

	return DM_PROC_Ok;
}

// fiber selection
//=================================================================================
// Name:	FiberSelectByRoiMap()
// Purpose:	Fiber selection via ROI maps
//
// Entry:	pstFactMap		-- data structure
//
// Return:  DM_PROC_Ok, if successful. 	selected fibers are marked as "Selected"
//			otherwise, error code defined in fact.h
//=================================================================================
int FiberSelectByRoiMap(FACT_MAP*	pstFactMap)
{
	if (pstFactMap->pstFactPara->szBinRoiFile[0] == NULL) {
		printf("no ROI files");
		return	DM_ERR_Unknow;
	}

	// get image dimensions
	int	nImgWidth = pstFactMap->pstFactPara->nImgWidth;
	int	nImgHeight = pstFactMap->pstFactPara->nImgHeight;
	int	nImgSlices = pstFactMap->pstFactPara->nImgSlices;
	long	lwImgBlkSize = (long)nImgWidth*nImgHeight*nImgSlices;

	// get fiber data
	FIBER*	pstFiberAry = pstFactMap->pstFiberAry;		// fiber data array
	int		nFiberArySize = pstFactMap->lFiberNrTotal;
	int**	pnFiberIdxAry = pstFactMap->pnFiberIdxAry;	// fiber data index of each voxel, 
	int*	pnFiberCntAry = pstFactMap->pnFiberCntAry;	// number of fibers through each voxel

	// get ROI files array and ROI operation code
	char**	szRoiFile = pstFactMap->pstFactPara->szBinRoiFile;
	int*	nRoiOp = pstFactMap->pstFactPara->nBinRoiOp;
	int*	nRoiFileType = pstFactMap->pstFactPara->szRoiFileType;

	//-------------------------------------------
	// do fiber labeling by binary ROI files
	//-------------------------------------------
	for (int i=0; i<MAX_ROI_FILE_Nr; i++) {
		// allocate memory for input data
		if (szRoiFile[i] == NULL) break;	// no more ROI files

		unsigned char*  pByteImg = (unsigned char*)malloc(sizeof(unsigned char) * lwImgBlkSize);
		if (!pByteImg) {
			DispErrorMsg(DM_ERR_MemAlloc, (void*)szRoiFile[i]);
			return	DM_ERR_MemAlloc;
		}
		memset(pByteImg,0,sizeof(unsigned char)*lwImgBlkSize);	// clear memory

		if(nRoiFileType[i] == 0)
		{

			// Open raw file
			FILE* inFile;
			if( (inFile = fopen( szRoiFile[i], "rb" )) == NULL ) {
				DispErrorMsg(DM_ERR_FileOpen, (void*)szRoiFile[i]);
				fclose(inFile);
				return DM_ERR_FileOpen;
			}

			// read file
			if ( fread(pByteImg, sizeof(unsigned char),  lwImgBlkSize, inFile ) != (size_t)lwImgBlkSize) {
				DispErrorMsg(DM_ERR_FileRead, (void*)szRoiFile[i]);
				fclose(inFile);
				free(pByteImg);
				return DM_ERR_FileRead;
			}
			fclose(inFile);

		}
		else if(nRoiFileType[i] == 1)
		{
			// Open mgh file
			mgh_header RoiMghInfo;
        		int ret = read_mgh(szRoiFile[i], &RoiMghInfo, 1);
			if(RoiMghInfo.type == MRI_UCHAR)
				for(int j=0; j<lwImgBlkSize; j++)
					pByteImg[j] = RoiMghInfo.charImageArray[j];
			else if(RoiMghInfo.type == MRI_INT)
				for(int j=0; j<lwImgBlkSize; j++)
					pByteImg[j] = (unsigned char)RoiMghInfo.intImageArray[j];
			else if(RoiMghInfo.type == MRI_LONG)
				for(int j=0; j<lwImgBlkSize; j++)
					pByteImg[j] = (unsigned char)RoiMghInfo.longImageArray[j];
			else if(RoiMghInfo.type == MRI_FLOAT)
				for(int j=0; j<lwImgBlkSize; j++)
					pByteImg[j] = (unsigned char)RoiMghInfo.floatImageArray[j];
			else if(RoiMghInfo.type == MRI_SHORT)
				for(int j=0; j<lwImgBlkSize; j++)
					pByteImg[j] = (unsigned char)RoiMghInfo.shortImageArray[j];
			
		}
		// labeling the fibers
		LabelFiberViaRoiImg(pstFiberAry,nFiberArySize, pnFiberIdxAry, pnFiberCntAry,
						pByteImg, nImgWidth, nImgHeight, nImgSlices,
						nRoiOp[i]								// ROI-operator
						);
		// release the byte image

		free(pByteImg);
	}
	return DM_PROC_Ok;
}

//=================================================================================
// Name:	LabelFiberViaRoiImg()
// Purpose:	Fiber selection via ROI map, called by FiberSelectByRoiMap()
//
// Entry:	stFiberAry[]		-- fiber data array, its size is (nImgW x nImgH x nImgS)
//			nFiberIdxAry[]		-- fiber data index of each voxel, 
//			nFiberCntAry[]		-- number of fibers through each voxel
//			nFiberNr			-- total number of filbers (or size of stFiberAry);
//
//			pByteImg			-- binary ROI file
//			nImgWidth/Height/Slices		-- image dimension
//			nRoiOp				-- ROI operation
//
// Return:  DM_PROC_Ok, if successful. 	selected fibers are marked as "Selected"
//			otherwise, error code defined in fact.h
//=================================================================================
void LabelFiberViaRoiImg(FIBER	stFiberAry[], int nFiberArySize, int*	nFiberIdxAry[], int	nFiberCntAry[], 
						unsigned char* pByteImg,						// binary ROI file
						int nImgWidth, int nImgHeight, int nImgSlices,		// image dimension
						int nRoiOp							// ROI operation
						)
{
	long dwFactOff = 0;
	int	xx,yy,zz, i,j;
	for (zz = 0; zz<nImgSlices; zz++) {
		long dwFactSliceOff = (long)zz * nImgWidth * nImgHeight;

		for (yy = 0; yy<nImgHeight; yy++) {
			long dwFactRowOff = (long)yy * nImgWidth;

			for (xx = 0; xx<nImgWidth; xx++, dwFactOff++) {

				// testing: ROI point?
				if (pByteImg[dwFactOff] == 0) continue;

				// Labeling fibers passint through this pixel
				long	dwFactVoxelOff = dwFactSliceOff + dwFactRowOff + xx;

				// no fibers, get off of here
				if (nFiberCntAry[dwFactVoxelOff] <= 0) continue;

				// labeling every fibers
				int* nFiberIdx = nFiberIdxAry[dwFactVoxelOff];			// get fiber index of this pixel
				for (i=0; i<nFiberCntAry[dwFactVoxelOff]; i++) {	// for all fibers through this pixel, do..

					stFiberAry[nFiberIdx[i]].xyzPtROI.x = (float)xx;			// remember the ROI point
					stFiberAry[nFiberIdx[i]].xyzPtROI.y = (float)yy;
					stFiberAry[nFiberIdx[i]].xyzPtROI.z = (float)zz;

					switch (nRoiOp) {
					case ROI_OpOr:
						stFiberAry[nFiberIdx[i]].nSelStatus = 1;	// = selected (OR operation) 
						stFiberAry[nFiberIdx[i]].nSelBeginIdx = 0;	// index start 0
						stFiberAry[nFiberIdx[i]].nSelEndIdx = stFiberAry[nFiberIdx[i]].nLength-1;	// end point of selection
						break;
					case ROI_OpNot:
						stFiberAry[nFiberIdx[i]].nSelStatus = 0;	// =delected (NOT operation)
						break;
					case ROI_OpAnd:		
						if (stFiberAry[nFiberIdx[i]].nSelStatus == 1) 	// been seleceted before?
							stFiberAry[nFiberIdx[i]].nSelStatus = 2;		// pre-selected, mark it as 2 (AND operation) 
						break;

					case ROI_OpCutOr:			// ROI_Cut0, similar to ROI_OR operation
						stFiberAry[nFiberIdx[i]].nSelStatus = 0x81;		// = selected (OR operation) 

						// find the hit-point of the fiber with current ROI-Pt, remember it as "SelStart"
						{
							XYZ_TRIPLE*	pxyzChain = stFiberAry[nFiberIdx[i]].pxyzChain;
							int			nChainLen = stFiberAry[nFiberIdx[i]].nLength;
							stFiberAry[nFiberIdx[i]].nSelBeginIdx = 0;			// initial value
							for (j=0; j<nChainLen; j++) {
								if (((long)pxyzChain[j].x + 
									 (long)pxyzChain[j].y*nImgWidth +
									 (long)pxyzChain[j].z*nImgWidth*nImgHeight) == dwFactOff) {
									stFiberAry[nFiberIdx[i]].nSelBeginIdx = j;
									break;
								}
							}
						}
						break;

					case ROI_OpCutAnd:			// ROI_Cut1, similar to ROI_AND operation
						 if (stFiberAry[nFiberIdx[i]].nSelStatus == 0x81) {	// been seleceted before?

							stFiberAry[nFiberIdx[i]].nSelStatus = 0x82;	// pre-selected (AND operation) 

							// find the hit-point of the fiber with current ROI-Pt, remember it as "SelEnd"
							{
								XYZ_TRIPLE*	pxyzChain = stFiberAry[nFiberIdx[i]].pxyzChain;
								int			nChainLen = stFiberAry[nFiberIdx[i]].nLength;
								stFiberAry[nFiberIdx[i]].nSelEndIdx = stFiberAry[nFiberIdx[i]].nLength-1;			// initial value
								for (j=0; j<nChainLen; j++) {
									if (((long)pxyzChain[j].x + 
										 (long)pxyzChain[j].y*nImgWidth +
										 (long)pxyzChain[j].z*nImgWidth*nImgHeight) == dwFactOff) {

										// note: nSelStart should be less then nSelEnd
										if (stFiberAry[nFiberIdx[i]].nSelBeginIdx <= j)	
											stFiberAry[nFiberIdx[i]].nSelEndIdx = j;
										else {
											stFiberAry[nFiberIdx[i]].nSelEndIdx = stFiberAry[nFiberIdx[i]].nSelBeginIdx;
											stFiberAry[nFiberIdx[i]].nSelBeginIdx = j;
										}
										break;
									}
								}
							}
						 }
						 break;

					default:
						break;
					}	// end of switch(RoiOP)
				}	// end of for (i=.. every fibers
			}	// end of for (xx=.. image width
		}	// end of for (yy=.. image height
	}	// end of for (zz=.. image slices


	// re-labeling the FiberSelect status for roi-AND operation, after all pixel in ROI has been processed.
	if (nRoiOp == ROI_OpAnd) {
		for (i=0; i<nFiberArySize; i++) {		// for all fiber chains do..
			if (stFiberAry[i].nSelStatus == 2)	{		// marked as OK for and-operation
				stFiberAry[i].nSelStatus = 1;			// set it as selected
				stFiberAry[i].nSelBeginIdx = 0;			// index start 0
				stFiberAry[i].nSelEndIdx = stFiberAry[i].nLength-1;	// end point of selection
			}
			else
				stFiberAry[i].nSelStatus = 0;	// otherwise, set it as deselected
		}
	}

	// re-labeling the FiberSelect status for roi-Cut-1 operation, after all pixel in ROI has been processed.
	else if (nRoiOp == ROI_OpCutAnd) {
		for (i=0; i<nFiberArySize; i++) {		// for all fiber chains do..

			if (stFiberAry[i].nSelStatus == 0x81)	// de-select fibers excluded in ROI-CutAnd, although it was selected by ROI-CutOR
				stFiberAry[i].nSelStatus = 0;
		}
	}
}

// fiber save
//=================================================================================
// Name:	FiberSave()
// Purpose:	Save fibers as data file
//
// Entry:	pstFactMap		-- data structure
//
// Return:  DM_PROC_Ok, if successful. 	otherwise, error code defined in fact.h
//=================================================================================
int FiberSave(FACT_MAP*	pstFactMap, bool bSaveSel)
{
	// define a file header for my fiber-data
	struct {
		char	sFiberFileTag[8];	// file tag = FiberDat

		int		nFiberNr;			// total number of fibers
		int		nFiberLenMax;		// max-length of fibers
		float	fFiberLenMean;		// mean-length of fibers

		int		nImgWidth;			// image dimension
		int		nImgHeight;
		int		nImgSlices;

		float	fPixelSizeWidth;	// voxel size
		float	fPixelSizeHeight;
		float	fSliceThickness;

//		SLICE_ORI	enumSliceOrientation;	// slice orientation
//		SLICE_SEQ	enumSliceSequencing;	// slice sequencing

		int			enumSliceOrientation;	// slice orientation
		int			enumSliceSequencing;	// slice sequencing

//		char	sVersion[8];		// version number
	} stFiberFileHeader;

	// get/create some parameters for file header
	int		nFiberNr = pstFactMap->lFiberNrTotal;
	int		nFiberLenMax = pstFactMap->lFiberLenMax;	// max-length of fibers
	float	fFiberLenMean = pstFactMap->fFiberLenMean;	// mean-length of fibers

	if (nFiberNr <= 0) {
		printf("No fibers in this dataset\n");
		return DM_ERR_Unknow;
	}

	int i,j;
	if (bSaveSel == true) {			// save selected fibers
		// counting: selected filbers
		nFiberNr = 0;
		nFiberLenMax = 0;
		fFiberLenMean = 0;
		for (i=0; i<pstFactMap->lFiberNrTotal; i++) {
			if (pstFactMap->pstFiberAry[i].nSelStatus != 0) {
				nFiberNr ++;
				nFiberLenMax = max(nFiberLenMax, pstFactMap->pstFiberAry[i].nLength);
				fFiberLenMean += pstFactMap->pstFiberAry[i].nLength;
			}
		}

		if (nFiberNr <= 0) {
			printf("No selected fibers in this dataset\n");
			return DM_ERR_Unknow;
		}

		fFiberLenMean /= nFiberNr;
	}
		
	// Create the fiber-data file
	char* szFileName = pstFactMap->pstFactPara->szFiberAllFile;
	if (bSaveSel == true) szFileName = pstFactMap->pstFactPara->szFiberSelFile;
	FILE* outFile = fopen(szFileName, "wb" );
	if (!outFile) {
		DispErrorMsg(DM_ERR_FileOpen, (void*)szFileName);
		return	DM_ERR_FileOpen;
	}

	// initializing some parameters for fiber-data header
	stFiberFileHeader.sFiberFileTag[0] = 'F';	// file tag = FiberDat
	stFiberFileHeader.sFiberFileTag[1] = 'i';
	stFiberFileHeader.sFiberFileTag[2] = 'b';
	stFiberFileHeader.sFiberFileTag[3] = 'e';
	stFiberFileHeader.sFiberFileTag[4] = 'r';
	stFiberFileHeader.sFiberFileTag[5] = 'D';
	stFiberFileHeader.sFiberFileTag[6] = 'a';
	stFiberFileHeader.sFiberFileTag[7] = 't';

	stFiberFileHeader.nImgWidth = pstFactMap->pstFactPara->nImgWidth;			// image dimension
	stFiberFileHeader.nImgHeight = pstFactMap->pstFactPara->nImgHeight;
	stFiberFileHeader.nImgSlices = pstFactMap->pstFactPara->nImgSlices;

	stFiberFileHeader.fPixelSizeWidth = pstFactMap->pstFactPara->fFovWidth / pstFactMap->pstFactPara->nImgWidth;			// voxel size
	stFiberFileHeader.fPixelSizeHeight = pstFactMap->pstFactPara->fFovHeight / pstFactMap->pstFactPara->nImgHeight;
	stFiberFileHeader.fSliceThickness = pstFactMap->pstFactPara->fSliceThickness;

	stFiberFileHeader.enumSliceOrientation = pstFactMap->pstFactPara->nImgOrientation;		// slice orientation
	stFiberFileHeader.enumSliceSequencing = pstFactMap->pstFactPara->nImgSequencing;		// slice sequencing

//	strncpy(stFiberFileHeader.sVersion, "2005.000", 8);							// version = YYYY.MM

	stFiberFileHeader.nFiberNr = nFiberNr;						// number of fibers
	stFiberFileHeader.nFiberLenMax = nFiberLenMax;				// max-length of fibers
	stFiberFileHeader.fFiberLenMean = fFiberLenMean;			// mean-length of fibers

	// save file header
	fwrite(&stFiberFileHeader, sizeof(stFiberFileHeader), 1, outFile);	// save file header

	// save fibers one-by-one, starts from file_offset = 128;
	fseek(outFile, 128, SEEK_SET);

	FIBER*		pstFiberAry = pstFactMap->pstFiberAry;
	unsigned char	bFiberSel = 0;
	if (bSaveSel == true) bFiberSel=1;

//	RGB_TRIPLE	rgbFiberClr = {255, 0, 0};		// define a fiber color
	RGB_TRIPLE	rgbFiberClr;					// generate random fiber color
	int			nFiberStart=0;
	int			nFiberLen;

	srand( (unsigned)time( NULL ) );
	float	fRandFactor = (float)255/RAND_MAX;

	for (i=0; i<pstFactMap->lFiberNrTotal; i++) {

		if (bSaveSel == true) {		// save selece fibers
			if (pstFiberAry[i].nSelStatus == 0) continue;
			if (pstFiberAry[i].nSelStatus == 0x82) {				// fiber cutting
				nFiberLen = pstFiberAry[i].nSelEndIdx - pstFiberAry[i].nSelBeginIdx +1;
				fwrite(&nFiberLen, sizeof(int), 1, outFile);				// save length
			}
			else {
				fwrite(&pstFiberAry[i].nLength, sizeof(int), 1, outFile);	// save length
			}

			fwrite(&bFiberSel, sizeof(unsigned char), 1, outFile);					// save sel-status

//			fwrite(&pstFiberAry[i].rgbColor, sizeof(RGB_TRIPLE));			// save color
			rgbFiberClr.r = int(255.0*rand() / RAND_MAX);
			rgbFiberClr.g = int(255.0*rand() / RAND_MAX);
			rgbFiberClr.b = int(255.0*rand() / RAND_MAX);
			fwrite(&rgbFiberClr, sizeof(RGB_TRIPLE), 1, outFile);			// save color

			// data fiber begin/end index (=0..fiber length
			if (pstFiberAry[i].nSelStatus == 0x82) {				// fiber cutting
				fwrite(&nFiberStart, sizeof(int), 1, outFile);			// save start index
				nFiberLen --;
				fwrite(&(nFiberLen), sizeof(int), 1, outFile);			// save end index
				nFiberLen ++;
			}
			else {
				fwrite(&pstFiberAry[i].nSelBeginIdx, sizeof(int), 1, outFile);	// save start index
				fwrite(&pstFiberAry[i].nSelEndIdx, sizeof(int), 1, outFile);	// save end index
			}

			// data fiber data
			if (pstFiberAry[i].nSelStatus == 0x82) {				// fiber cutting
				for (j=pstFiberAry[i].nSelBeginIdx; j<=pstFiberAry[i].nSelEndIdx; j++)
					fwrite(&pstFiberAry[i].pxyzChain[j], sizeof(XYZ_TRIPLE), 1, outFile);
			}
			else {
				fwrite(pstFiberAry[i].pxyzChain, sizeof(XYZ_TRIPLE), pstFiberAry[i].nLength, outFile);
				FILE *fp;
				fp = fopen ( "./points.txt", "a");
				for(j=0; j< pstFiberAry[i].nLength; j++)
					fprintf(fp, "(x, y, z) = (%3.1f, \t%3.1f, \t%3.1f)\n", pstFiberAry[i].pxyzChain[j].x, 
											       pstFiberAry[i].pxyzChain[j].y,
											       pstFiberAry[i].pxyzChain[j].z);
				fclose(fp);
				
				/* for test by sumiko
				XYZ_TRIPLE *tempChain = (XYZ_TRIPLE *)calloc( pstFiberAry[i].nLength, sizeof(XYZ_TRIPLE));
				for(j=0; j< pstFiberAry[i].nLength; j++)
				{

					tempChain[j].x=pstFiberAry[i].pxyzChain[j].x;
					tempChain[j].y=pstFiberAry[i].pxyzChain[j].y;
					tempChain[j].z=pstFiberAry[i].pxyzChain[j].z;
				}
				fwrite(tempChain, sizeof(XYZ_TRIPLE), pstFiberAry[i].nLength, outFile);
				free(tempChain);
				*/
					
			}
		}
		else {
			fwrite(&pstFiberAry[i].nLength, sizeof(int), 1, outFile);		// save length

	//		fwrite(&pstFiberAry[i].nSelStatus, sizeof(unsigned char), 1, outFile);	// save sel-status
			fwrite(&bFiberSel, sizeof(unsigned char), 1, outFile);					// save sel-status

	//		fwrite(&pstFiberAry[i].rgbColor, sizeof(RGB_TRIPLE));			// save color
			rgbFiberClr.r = int(255.0*rand() / RAND_MAX);
			rgbFiberClr.g = int(255.0*rand() / RAND_MAX);
			rgbFiberClr.b = int(255.0*rand() / RAND_MAX);
			fwrite(&rgbFiberClr, sizeof(RGB_TRIPLE), 1, outFile);			// save color

			fwrite(&pstFiberAry[i].nSelBeginIdx, sizeof(int), 1, outFile);	// save begin-end index
			fwrite(&pstFiberAry[i].nSelEndIdx, sizeof(int), 1, outFile);

			// save fiber data
			fwrite(pstFiberAry[i].pxyzChain, sizeof(XYZ_TRIPLE), pstFiberAry[i].nLength, outFile);
		}
	} // end of "for (i=0..
	fclose(outFile);

	return DM_PROC_Ok;
}

//=================================================================================
// Name:	FiberSaveTxt()
// Purpose:	Save fibers as a plain-text file
//
// Entry:	pstFactMap		-- data structure
//			bSaveSel		-- =TRUE if save selected fibers,  =FALSE if save all
//
// Return:  DM_PROC_Ok, if successful. 	otherwise, error code defined in fact.h
//=================================================================================
int FiberSaveTxt(FACT_MAP*	pstFactMap, bool bSaveSel)
{
	// get/create some parameters for file header
	int		nFiberNr = pstFactMap->lFiberNrTotal;
	int		nFiberLenMax = pstFactMap->lFiberLenMax;	// max-length of fibers
	float	fFiberLenMean = pstFactMap->fFiberLenMean;	// mean-length of fibers

	if (nFiberNr <= 0) {
		printf("No fibers in this dataset\n");
		return DM_ERR_Unknow;
	}

	int i,j;
	if (bSaveSel == true) {			// save selected fibers
		// counting: selected filbers
		nFiberNr = 0;
		nFiberLenMax = 0;
		fFiberLenMean = 0;
		for (i=0; i<pstFactMap->lFiberNrTotal; i++) {
			if (pstFactMap->pstFiberAry[i].nSelStatus != 0) {
				nFiberNr ++;
				nFiberLenMax = max(nFiberLenMax, pstFactMap->pstFiberAry[i].nLength);
				fFiberLenMean += pstFactMap->pstFiberAry[i].nLength;
			}
		}

		if (nFiberNr <= 0) {
			printf("No selected fibers in this dataset\n");
			return DM_ERR_Unknow;
		}

		fFiberLenMean /= nFiberNr;
	}
		
	// create a file
	char* szFileName = pstFactMap->pstFactPara->szFiberAllTxtFile;
	if (bSaveSel == true) szFileName = pstFactMap->pstFactPara->szFiberSelTxtFile;

	FILE* outFile = fopen(szFileName, "wt" );
	if (outFile == NULL) {
		DispErrorMsg(DM_ERR_FileOpen, (void*)szFileName);
		return	DM_ERR_FileOpen;
	}

	// write information header
	if (bSaveSel == true)
		fprintf(outFile, "Fiber tracts (selected ones) by fiber tracking toolkit: fact\n\n");	
	else
		fprintf(outFile, "Fiber tracts (all fibers) by fiber tracking toolkit: fact\n\n");	

	// image dimensions
	fprintf(outFile, "Image Width:\t  %d\n", pstFactMap->pstFactPara->nImgWidth);
	fprintf(outFile, "Image Height:\t  %d\n",pstFactMap->pstFactPara->nImgHeight);
	fprintf(outFile, "Image Slices:\t  %d\n", pstFactMap->pstFactPara->nImgSlices);
	fprintf(outFile, "\n");

	// pixel size
	fprintf(outFile, "Pixel Size (Width):\t  %f\n", pstFactMap->pstFactPara->fFovWidth / pstFactMap->pstFactPara->nImgWidth);
	fprintf(outFile, "Pixel Size (Height):\t  %f\n", pstFactMap->pstFactPara->fFovHeight / pstFactMap->pstFactPara->nImgHeight);
	fprintf(outFile, "Pixel Size (Thickness):\t  %f\n", pstFactMap->pstFactPara->fSliceThickness);
	fprintf(outFile, "\n");

	// save fiber counter
	fprintf(outFile, "Number of the fibers:\t  %d\n", nFiberNr);
	fprintf(outFile, "Max-length of the fibers:\t  %d\n", nFiberLenMax);
	fprintf(outFile, "Mean-length of the fibers:\t  %f\n", fFiberLenMean);
	fprintf(outFile, "\n");

	// Fiber coordinates,  one by one
	fprintf(outFile, "Unit of the coordinates:	Image Matrix\n");
	fprintf(outFile, "Format:  Fiber_Idx,  Fiber_Len, Seed_Pt[, ROI_Pt]\n");
	fprintf(outFile, "         fiber coordinates (x, y, z)\n\n");


	// Ok, output fiber data
	FIBER*		pstFiberAry = pstFactMap->pstFiberAry;
	XYZ_TRIPLE	xyzPt;			// a temp variabel
	int			nFiberChainStart, nFiberChainEnd;

	int		nFiberCnt = 0;		// fiber cunter
	int		nItemsPerLine = 0;	// items per line
	for (i=0; i<pstFactMap->lFiberNrTotal; i++) {
		if ((bSaveSel == true) && (pstFiberAry[i].nSelStatus == 0)) continue;

		// determining start/end pointer of the fiber chain
		nFiberChainStart = 0;
		nFiberChainEnd = pstFiberAry[i].nLength;
		if ((bSaveSel == true) && (pstFiberAry[i].nSelStatus == 0x82)) { // for fiber-cut-by-ROIs;
			nFiberChainStart = pstFiberAry[i].nSelBeginIdx;
			nFiberChainEnd = pstFiberAry[i].nSelEndIdx+1;
		}

		// print fiber index and length
		fprintf(outFile, "\nFiber:%8d \t Length:%8d\t", ++nFiberCnt, nFiberChainEnd-nFiberChainStart); 

		// output seed point and roi point
		if (bSaveSel == true)
			fprintf(outFile, "Seed: (%5.1f %5.1f %5.1f)  ROI: (%5.1f %5.1f %5.1f)\n",
						pstFiberAry[i].xyzPtSeed.x+0.5, pstFiberAry[i].xyzPtSeed.y+0.5, pstFiberAry[i].xyzPtSeed.z+0.5,
						pstFiberAry[i].xyzPtROI.x,  pstFiberAry[i].xyzPtROI.y,  pstFiberAry[i].xyzPtROI.z);
		else 
			fprintf(outFile, "Seed: (%5.1f %5.1f %5.1f)\n",
						pstFiberAry[i].xyzPtSeed.x+0.5, pstFiberAry[i].xyzPtSeed.y+0.5, pstFiberAry[i].xyzPtSeed.z+0.5);

		// output fiber coordinate 
		nItemsPerLine = 0;	// remember how many items have been print out
		for (j=nFiberChainStart; j<nFiberChainEnd; j++) {

			// output formated data in Image-Matrix coordinate
			xyzPt = pstFiberAry[i].pxyzChain[j];
			fprintf(outFile, "%13.4f %13.4f %13.4f",xyzPt.x, xyzPt.y, xyzPt.z);

			// print line-space
			fprintf(outFile, "\n"); continue;

			// or: print line-space if necessary
			nItemsPerLine ++;
			if ((nItemsPerLine % 2) == 0) {		// two vertices per line
				fprintf(outFile, "\n");
				nItemsPerLine = 0;
			}
		}
	}
	fprintf(outFile, "\n");
	fclose(outFile);
	return DM_PROC_Ok;
}


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

	cout << "\n\t####Reading a MGH file: "<< mghFileName << endl ;

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

	mghInfor->xFOV = mghInfor->xyz_size[0] *  mghInfor->width;
	mghInfor->yFOV = mghInfor->xyz_size[1] *  mghInfor->height;

	//fseek(fp, unused_space_size, SEEK_CUR) ;
	//cout << "current position : "<< ftell(fp) << endl ;

	int nv = mghInfor->width * mghInfor->height * mghInfor->depth *mghInfor->nframes;

		if(mghInfor->type ==MRI_UCHAR)
		{
			//cout << "Type: unsigned char" << endl;
			fseek(fp, 0, SEEK_END);
			int fileLen = ftell(fp);
			int headerSize = fileLen - nv*sizeof(unsigned char);
			mghInfor->charImageArray = (unsigned char *)calloc(nv, sizeof(unsigned char)); 
			fseek(fp, headerSize-16, SEEK_SET) ;
			fread(mghInfor->charImageArray, sizeof(unsigned char), nv, fp);
		}
		else if(mghInfor->type ==MRI_INT)
		{
			//cout << "Type: int" << endl;
			fseek(fp, 0, SEEK_END);
			int fileLen = ftell(fp);
			int headerSize = fileLen - nv*sizeof(int);
			mghInfor->intImageArray = (int *)calloc(nv, sizeof(int)); 
			fseek(fp, headerSize-16, SEEK_SET) ;
			fread(mghInfor->intImageArray, sizeof(int), nv, fp);
			if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
				for(i=0; i<nv; i++)
					bswap(mghInfor->intImageArray[i]);
		}
		else if(mghInfor->type ==MRI_LONG)
		{
			//cout << "Type: long" << endl;
			fseek(fp, 0, SEEK_END);
			int fileLen = ftell(fp);
			int headerSize = fileLen - nv*sizeof(MRI_LONG);
			mghInfor->longImageArray = (long *)calloc(nv, sizeof(long)); 
			fseek(fp, headerSize-16, SEEK_SET) ;
			fread(mghInfor->longImageArray, sizeof(long), nv, fp);
			if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
				for(i=0; i<nv; i++)
					bswap(mghInfor->longImageArray[i]);
		}
		else if(mghInfor->type ==MRI_FLOAT)
		{
			//cout << "Type: float" << endl;
			fseek(fp, 0, SEEK_END) ;
			int fileLen = ftell(fp);
			int headerSize = fileLen - nv*sizeof(MRI_FLOAT);
			mghInfor->floatImageArray = (float *)calloc(nv, sizeof(float)); 
			//fseek(fp, -nv*sizeof(MRI_FLOAT), SEEK_END) ;
			fseek(fp, headerSize-16, SEEK_SET) ;
			fread(mghInfor->floatImageArray, sizeof(float), nv, fp);
			if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0))
				for(i=0; i<nv; i++)
					bswap(mghInfor->floatImageArray[i]);

		}
		else if(mghInfor->type ==MRI_SHORT)
		{
			//cout << "Type: short" << endl;
			fseek(fp, 0, SEEK_END) ;
			int fileLen = ftell(fp);
			int headerSize = fileLen - nv*sizeof(short);
			mghInfor->shortImageArray = (short *)calloc(nv, sizeof(short)); 
			fseek(fp, headerSize-16, SEEK_SET) ;
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
/*
fopen("./testfa4.raw", "w");
fwrite(mghInfor->floatImageArray, sizeof(float), nv, fp);
fclose(fp);
*/

	return 1;

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

//get FACR_MAP structure
FACT_MAP getFactmap(int argc, char* argv[] )
{

	int i, j, k;
	mgh_header FAmghInfo;
	mgh_header V0mghInfo;


	// input file is a parameter text file
	FACT_MAP *stFactMap = (FACT_MAP*)calloc(1, sizeof(FACT_MAP));
	stFactMap->pstFactPara = (FACT_PARA*)calloc(1, sizeof(FACT_PARA));

	char *outputAllFiberFileBin;
	char *outputSelFiberFileBin;
	char *outputAllFiberFileTxt;
	char *outputSelFiberFileTxt;
	char *outputSelFiberMghFile;

	char **inputRoiMghOR;
	char **inputRoiMghAND;
	char **inputRoiMghNOT;
	char **inputRoiMghCUT_OR;
	char **inputRoiMghCUT_AND;


	int roiNum = 0;
	int mghFileInput[2];
	mghFileInput[0] = 0;
	mghFileInput[1] = 0;

	optionStruct option_st[INPUTARGV_TYPE_NUM];
	int res = createOptionFlag( argc,  argv, option_st);
	if (argc == 3 && strcmp(argv[1], "-iParatxt") == 0)
	{

		int nRetCode = ReadFactParaFile(argv[2], stFactMap->pstFactPara);
        	if (nRetCode != DM_PROC_Ok) {
			DispErrorMsg(DM_ERR_FileRead, (void*)argv[3], (void*)"check parameter file.");
                	FreeMemFactPara (stFactMap->pstFactPara);
			cout << "read txt file failed. " << endl; 
                	return *stFactMap;
        	}
		else 
		{
			printf("read image files....\n");
			nRetCode = ReadImgFiles(stFactMap);
			if (nRetCode != DM_PROC_Ok) {
				FreeMemFactMap (stFactMap);
				FreeMemFactPara (stFactMap->pstFactPara);
				DispErrorMsg(DM_ERR_FileRead, (void*)argv[3], (void*)"check image file.");
				return *stFactMap;
			}
			cout << "text file read completed. " << endl; 
			getchar();
			return *stFactMap;
		}
	}
	else
	{ 

		for(j=0; j<argc; j++)
		{
			if( strcmp(argv[j], "-iFAmgh") == 0 )
				mghFileInput[0] = 1;
			if( strcmp(argv[j], "-iV0mgh") == 0 )
				mghFileInput[1] = 1;
		}
		stFactMap->pstFactPara->mghFileInput[0] = mghFileInput[0];
		stFactMap->pstFactPara->mghFileInput[1] = mghFileInput[1];
		for(i=0; i<argc; i++)
		{
			if(strcmp(argv[i], "-iFAmgh") == 0)
			{
				stFactMap->pstFactPara->szMghAniFile = strdup( argv[i+1] );
        			int ret = read_mgh(argv[i+1], &FAmghInfo, 1);
				if (ret == -1)
				{
					cout << "FA MGH file error: check mgh file format." << endl;
					return *stFactMap;
				}
			}
			if(strcmp(argv[i], "-iV0mgh") == 0)
			{
				stFactMap->pstFactPara->szMghVecFile = strdup( argv[i+1] );
        			int ret = read_mgh(argv[i+1], &V0mghInfo, 1);
				if (ret == -1)
				{
					cout << "FA MGH file error: check mgh file format." << endl;
					return *stFactMap;
				}
			}
			if(strcmp(argv[i], "-iFAraw") == 0)
				stFactMap->pstFactPara->szImgAniFile = strdup( argv[i+1] );
			if(strcmp(argv[i], "-iV0raw") == 0)
				stFactMap->pstFactPara->szImgVecFile = strdup( argv[i+1] );
			if(strcmp(argv[i], "-iStartFA") == 0)
				stFactMap->pstFactPara->fStartFA = (float) atof( argv[i+1] );
			if(strcmp(argv[i], "-iStopFA") == 0)
				stFactMap->pstFactPara->fStopFA = (float) atof( argv[i+1] );
			if(strcmp(argv[i], "-iTurnAngleDeg") == 0)
				stFactMap->pstFactPara->fTurnDeg = (float) atof( argv[i+1] );
			if(strcmp(argv[i], "-iFlipEigenX") == 0)
				stFactMap->pstFactPara->bFlipVecX= (bool) atoi( argv[i+1] );
			if(strcmp(argv[i], "-iFlipEigenY") == 0)
				stFactMap->pstFactPara->bFlipVecY= (bool) atoi( argv[i+1] );
			if(strcmp(argv[i], "-iFlipEigenZ") == 0)
				stFactMap->pstFactPara->bFlipVecY= (bool) atoi( argv[i+1] );
			if(strcmp(argv[i], "-iRoimghOR") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum] = strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum] = 1;
				stFactMap->pstFactPara->nBinRoiOp[roiNum] = ROI_OpOr;
				roiNum ++;
			}	
			if(strcmp(argv[i],"-iRoimghAND") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 1;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpAnd;
				roiNum ++;
			}
			if(strcmp(argv[i], "-iRoimghNOT") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 1;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpNot;
				roiNum ++;
			}	
			if(strcmp(argv[i], "-iRoimghCUT_OR") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 1;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpCutOr;
				roiNum ++;
			}
			if(strcmp(argv[i], "-iRoimghCUT_AND") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 1;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpCutAnd;
				roiNum ++;
			}	
			if(strcmp(argv[i], "-iRoirawOR") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 0;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpOr;
				roiNum ++;
			}	
			if(strcmp(argv[i], "-iRoirawAND") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 0;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpAnd;
				roiNum ++;
			}	
			if(strcmp(argv[i],  "-iRoirawNOT") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 0;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpNot;
				roiNum ++;
			}	
			if(strcmp(argv[i], "-iRoirawCUT_OR") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 0;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpCutOr;
				roiNum ++;
			}	
			if(strcmp(argv[i], "-iRoirawCUT_AND") == 0)
			{
				stFactMap->pstFactPara->szBinRoiFile[roiNum]= strdup( argv[i+1] );
				stFactMap->pstFactPara->szRoiFileType[roiNum]= 0;
				stFactMap->pstFactPara->nBinRoiOp[roiNum]= ROI_OpCutAnd;
				roiNum ++;
			}	
			if( mghFileInput[0] == 0)
			{
					if(strcmp(argv[i], "-iImgDims") == 0)
					{
						stFactMap->pstFactPara->nImgWidth= atoi(argv[i+1] );
						stFactMap->pstFactPara->nImgHeight= atoi(argv[i+2] );
						stFactMap->pstFactPara->nImgSlices= atoi(argv[i+3] );
					}	
					if(strcmp(argv[i],  "-iVoxSize") == 0)
					{
						stFactMap->pstFactPara->fPixelSizeWidth= atof(argv[i+1] );
						stFactMap->pstFactPara->fPixelSizeHeight= atoi(argv[i+2] );
					}
					if(strcmp(argv[i],  "-iSliceThickness") == 0)
					{	
						stFactMap->pstFactPara->fSliceThickness= atoi(argv[i+1] );
					}	
					if(strcmp(argv[i], "-iFOVSize") == 0)
					{
						stFactMap->pstFactPara->fFovWidth= atof(argv[i+1] );
						stFactMap->pstFactPara->fFovHeight= atoi(argv[i+2] );
					}	
					if(strcmp(argv[i], "-iImageOrientation") == 0)
						stFactMap->pstFactPara->nImgOrientation= atoi(argv[i+1] ); // coronal: 0, axial: 1, sagittal: 2
		
			}
			if(strcmp(argv[i], "-oAllFiberFile") == 0)
				stFactMap->pstFactPara->szFiberAllFile= strdup(argv[i+1] ); 
			if(strcmp(argv[i], "-oSelFiberFile") == 0)
				stFactMap->pstFactPara->szFiberSelFile= strdup(argv[i+1] ); 
			if(strcmp(argv[i], "-oAllFiberTxtFile") == 0)
				stFactMap->pstFactPara->szFiberAllTxtFile= strdup(argv[i+1] ); 
			if(strcmp(argv[i], "-oSelFiberTxtFile") == 0)
				stFactMap->pstFactPara->szFiberSelTxtFile= strdup(argv[i+1] ); 
			if(strcmp(argv[i], "-oSelFiberMghFile") == 0)
				stFactMap->pstFactPara->szFiberSelMghFile= strdup(argv[i+1] ); 
			if(strcmp(argv[i], "-oSelFiberVolFile") == 0)
				stFactMap->pstFactPara->szFiberSelVolFile= strdup(argv[i+1] ); 
			if(strcmp(argv[i],  "-iByteSwapFlag") == 0)
				stFactMap->pstFactPara->bSwapBytes = (bool) atoi( argv[i+1] );
			if(strcmp(argv[i], "-iFiberLenMin") == 0)
					stFactMap->pstFactPara->nFiberLenMin =  atoi( argv[i+1] );
		}
	}
/*
	getchar();

	cout << "x = " <<  FAmghInfo.width << endl; 
	cout << "y = " <<  FAmghInfo.height << endl; 
	cout << "z = " <<  FAmghInfo.depth << endl; 
	cout << "xFOV = " <<  FAmghInfo.xFOV << endl; 
	cout << "yFOV = " <<  FAmghInfo.yFOV << endl; 
	getchar();
*/
	cout << "stFactMap->pstFactPara->mghFileInput[0] = " << stFactMap->pstFactPara->mghFileInput[0] << endl;
	cout << "stFactMap->pstFactPara->mghFileInput[1] = " << stFactMap->pstFactPara->mghFileInput[1] << endl;
	cout << "roiNumber = " <<roiNum << endl;
	getchar();

	int nv;
	if(  stFactMap->pstFactPara->mghFileInput[0] == 1 )
	{
		nv = FAmghInfo.width * FAmghInfo.height * FAmghInfo.depth;
		stFactMap->pstFactPara->nImgWidth 	= FAmghInfo.width;
		stFactMap->pstFactPara->nImgHeight	= FAmghInfo.height;
		stFactMap->pstFactPara->nImgSlices	= FAmghInfo.depth;

		stFactMap->pstFactPara->fPixelSizeWidth	= FAmghInfo.xyz_size[0];
		stFactMap->pstFactPara->fPixelSizeHeight= FAmghInfo.xyz_size[1];
		stFactMap->pstFactPara->fSliceThickness	= FAmghInfo.xyz_size[2];

		stFactMap->pstFactPara->fFovWidth	= FAmghInfo.xFOV;
		stFactMap->pstFactPara->fFovHeight	= FAmghInfo.yFOV;

		//ImageOrientation: 0=coronal, 1=axial, 2=sagittal
		if (FAmghInfo.x_ras[0] != 0 ) 				//coronal ras matrix
			 stFactMap->pstFactPara->nImgOrientation = 0;
		else if (FAmghInfo.x_ras[1] != 0 ) 			//axial ras matrix
			 stFactMap->pstFactPara->nImgOrientation = 1;
		else							//sagittal ras matrix
			stFactMap->pstFactPara->nImgOrientation = 2;
		if(FAmghInfo.type != MRI_FLOAT)
		{
			cout << "Input anisotropy image or tensor vector image is not float type. Stop program..." << endl;
			return *stFactMap;
		}
		stFactMap->pfAniImg = (float *)calloc(nv, sizeof(float));
		for (i=0; i<nv; i++)
			stFactMap->pfAniImg[i] = FAmghInfo.floatImageArray[i];

		FILE *testfp;
		testfp = fopen("fa.raw","w");
		fwrite( stFactMap->pfAniImg, sizeof(float), nv, testfp); 
		fclose(testfp);
		cout << "fa.raw output!" << endl;
		getchar();
	}	
	if(  stFactMap->pstFactPara->mghFileInput[1] == 1 )
	{
		nv = V0mghInfo.width * V0mghInfo.height * V0mghInfo.depth;
		stFactMap->pstFactPara->nImgWidth 	= V0mghInfo.width;
		stFactMap->pstFactPara->nImgHeight	= V0mghInfo.height;
		stFactMap->pstFactPara->nImgSlices	= V0mghInfo.depth;

		stFactMap->pstFactPara->fPixelSizeWidth	= V0mghInfo.xyz_size[0];
		stFactMap->pstFactPara->fPixelSizeHeight= V0mghInfo.xyz_size[1];
		stFactMap->pstFactPara->fSliceThickness	= V0mghInfo.xyz_size[2];

		stFactMap->pstFactPara->fFovWidth	= V0mghInfo.xFOV;
		stFactMap->pstFactPara->fFovHeight	= V0mghInfo.yFOV;

		//ImageOrientation: 0=coronal, 1=axial, 2=sagittal
		if (V0mghInfo.x_ras[0] != 0 ) 				//coronal ras matrix
			 stFactMap->pstFactPara->nImgOrientation = 0;
		else if (V0mghInfo.x_ras[1] != 0 ) 			//axial ras matrix
			 stFactMap->pstFactPara->nImgOrientation = 1;
		else							//sagittal ras matrix
			stFactMap->pstFactPara->nImgOrientation = 2;

		if(V0mghInfo.type != MRI_FLOAT)
		{
			cout << "Input anisotropy image or tensor vector image is not float type. Stop program..." << endl;
			return *stFactMap;
		}
		stFactMap->pxyzVecImg = (XYZ_TRIPLE *)calloc(nv, sizeof(XYZ_TRIPLE));
		float* vec_x = (float *)calloc(nv, sizeof(float));
		float* vec_y = (float *)calloc(nv, sizeof(float));
		float* vec_z = (float *)calloc(nv, sizeof(float));
		

		for (i=0; i<nv; i++)
		{
			stFactMap->pxyzVecImg[i].x = V0mghInfo.floatImageArray[i];
			stFactMap->pxyzVecImg[i].y = V0mghInfo.floatImageArray[i+nv];
			stFactMap->pxyzVecImg[i].z = V0mghInfo.floatImageArray[i+2*nv];
			vec_x[i] = stFactMap->pxyzVecImg[i].x;
			vec_y[i] = stFactMap->pxyzVecImg[i].y;
			vec_z[i] = stFactMap->pxyzVecImg[i].z;
		}
		FILE *testfp;
		testfp = fopen("v0.raw","w");
		fwrite( vec_x, sizeof(float), nv, testfp); 
		fwrite( vec_y, sizeof(float), nv, testfp); 
		fwrite( vec_z, sizeof(float), nv, testfp); 
		//fwrite( stFactMap->pxyzVecImg.z, sizeof(float), nv, testfp); 
		fclose(testfp);
		cout << "v0.raw output!" << endl;
		getchar();
	}

	if(option_st[3].optionflag == 0)
		stFactMap->pstFactPara->fStartFA = setDefaultStartFA();
	if(option_st[4].optionflag == 0)
		stFactMap->pstFactPara->fStopFA = setDefaultStopFA();
	if(option_st[5].optionflag == 0)
		stFactMap->pstFactPara->fTurnDeg = setDefaultTurnAngleDeg();
	if(option_st[6].optionflag == 0)
		stFactMap->pstFactPara->bFlipVecX = setDefaultFlipVecX();
	if(option_st[7].optionflag == 0)
		stFactMap->pstFactPara->bFlipVecY = setDefaultFlipVecY();

	if(option_st[8].optionflag == 0)
		stFactMap->pstFactPara->bFlipVecZ = setDefaultFlipVecZ();
	else 
		stFactMap->pstFactPara->bFlipVecZ = 1;

	if(option_st[30].optionflag == 0) 
		stFactMap->pstFactPara->bSwapBytes = setDefaultByteSwap();
	else
		stFactMap->pstFactPara->bSwapBytes = 1;

	if(option_st[31].optionflag == 0)
		stFactMap->pstFactPara->nFiberLenMin = setDefaultFiberLenMin();

	cout << "startFA: " << stFactMap->pstFactPara->fStartFA << endl;
	cout << "stopFA: " << stFactMap->pstFactPara->fStopFA  << endl;
	cout << "turnAng: " << stFactMap->pstFactPara->fTurnDeg << endl;
	cout << "flipx: " <<  stFactMap-> pstFactPara->bFlipVecX << endl; 
	cout << "flipy: " <<  stFactMap-> pstFactPara->bFlipVecY << endl; 
	cout << "flipz: " <<  stFactMap-> pstFactPara->bFlipVecZ << endl; 
	cout << "x = " <<  stFactMap->pstFactPara->nImgWidth << endl; 
	cout << "y = " <<  stFactMap->pstFactPara->nImgHeight << endl; 
	cout << "z = " <<  stFactMap->pstFactPara->nImgSlices << endl; 
	getchar();

	if( (option_st[1].optionflag == 0 && option_st[2].optionflag == 0) &&
	    (option_st[19].optionflag == 0 && option_st[20].optionflag == 0)) 
	{
		cout << "Input file error, we need FA and Vector data." << endl ;
		cout << "Porgram stopped." << endl << endl;
		return *stFactMap;
	}
	else if( (option_st[1].optionflag == 0 && option_st[2].optionflag == 0) ) 
	{
		if( option_st[21].optionflag == 0 || option_st[22].optionflag == 0 || 
		    option_st[23].optionflag == 0 || option_st[24].optionflag == 0 )
		{
			cout << "Input arguments errors on -iImgDims -iVoxSize -iFOVsize iImageOrientation." << endl ;
			cout << "Porgram stopped." << endl << endl;
			return *stFactMap;
		}
		/*
		else if( (option_st[19].optionflag == 1 && option_st[20].optionflag == 1) &&
			 (option_st[1 ].optionflag == 0 && option_st[2 ].optionflag == 0))
		{
			cout << "read *.dat data" << endl;
			stFactMap->pfAniImg = (float *) calloc (nv, sizeof(float));
			stFactMap->pxyzVecImg = (XYZ_TRIPLE *) calloc (nv, sizeof(XYZ_TRIPLE));

			FILE *inFile;
			inFile = fopen(stFactMap->pstFactPara->szImgAniFile, "r");
			fread(stFactMap->pfAniImg, sizeof(float),  nv, inFile );
			fclose(inFile);

			inFile = fopen(stFactMap->pstFactPara->szImgVecFile, "r");
			fread(stFactMap->pxyzVecImg, sizeof(XYZ_TRIPLE),  nv, inFile );
			fclose(inFile);
			cout << "szImgVecFile = " << stFactMap->pstFactPara->szImgVecFile << endl; 
			cout << "szImgAniFile = " << stFactMap->pstFactPara->szImgAniFile << endl; 
			getchar();
		}
		*/
		
	}
	cout << "flipping..."<<endl;
	getchar();
	if ( stFactMap->pstFactPara->bFlipVecX == true) for (i=0; i<nv; i++)  stFactMap->pxyzVecImg[i].x = -stFactMap->pxyzVecImg[i].x;
	if ( stFactMap->pstFactPara->bFlipVecY == true) for (i=0; i<nv; i++)  stFactMap->pxyzVecImg[i].y = -stFactMap->pxyzVecImg[i].y;
	if ( stFactMap->pstFactPara->bFlipVecZ == true) for (i=0; i<nv; i++)  stFactMap->pxyzVecImg[i].z = -stFactMap->pxyzVecImg[i].z;
	float fPixSizeW = stFactMap->pstFactPara->fFovWidth  / stFactMap->pstFactPara->nImgWidth;
	float fPixSizeH = stFactMap->pstFactPara->fFovHeight / stFactMap->pstFactPara->nImgHeight;
	float fPixSizeD = stFactMap->pstFactPara->fSliceThickness;
	for (i=0; i<nv; i++) {
		stFactMap->pxyzVecImg[i].x /= fPixSizeW;
		stFactMap->pxyzVecImg[i].y /= fPixSizeH;
		stFactMap->pxyzVecImg[i].z /= fPixSizeD;
	}

	
/*
	stFactMap->pstFactPara->bFlipVecX ;
	if ( stFactMap->pstFactPara->bFlipVecX == true) for (i=0; i<nv; i++)  stFactMap->pxyzVecImg[i].x = -stFactMap->pxyzVecImg[i].x;
	if ( stFactMap->pstFactPara->bFlipVecY == true) for (i=0; i<nv; i++)  stFactMap->pxyzVecImg[i].y = -stFactMap->pxyzVecImg[i].y;
	if ( stFactMap->pstFactPara->bFlipVecZ == true) for (i=0; i<nv; i++)  stFactMap->pxyzVecImg[i].z = -stFactMap->pxyzVecImg[i].y;
	if (bFlipVecY == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].y = -pxyzVecImg[i].y;
	if (bFlipVecZ == true) for (i=0; i<lwImgBlkSize; i++)  pxyzVecImg[i].z = -pxyzVecImg[i].z;
	stFactMap->pstFactPara->bFlipVecX= (bool) atoi( argv[i+1] );
	for(i=0; i<numRoiFiles; i++)
	{
		cout << "RoiFile[" << i << "] : " << RoiFiles[i] ;
		cout << ";   RoiOperation[" << i << "] : " << RoiOperations[i] << endl;
	}
	getchar();
	*/

	return *stFactMap;
}

float setDefaultStartFA()
{
	float startFA = 0.2;
	cout <<"\tDefault start FA = 0.2 ." << endl;
	return startFA;
	
}
float setDefaultStopFA()
{
	float stopFA = 0.2;
	cout <<"\tDefault stot FA = 0.2 ." << endl;
	return stopFA;
}
float setDefaultTurnAngleDeg()
{
	float turnAngleDeg = 20;
	cout <<"\tDefault turn angle degree = 20 deg" << endl;
	return turnAngleDeg;
}
bool setDefaultByteSwap()
{
	bool byteSwapFlag = 0;
	cout <<"\tDefault byte swap: 0." << endl;
	return byteSwapFlag;
}
int setDefaultFiberLenMin()
{
	int fiberLenMin = 7;
	cout <<"\tDefault setting of minimum fiber length is 7." << endl;
	return fiberLenMin;
}
bool setDefaultFlipVecX()
{
	int flip = 0;
	cout <<"\tDefault setting of vector flip in X direction: 0." << endl;
	return flip;
}
bool setDefaultFlipVecY()
{
	int flip = 0;
	cout <<"\tDefault setting of vector flip in Y direction: 0." << endl;
	return flip;
}
bool setDefaultFlipVecZ()
{
	int flip = 0;
	cout <<"\tDefault setting of vector flip in Z direction: 0." << endl;
	return flip;
}
	

int readFileList(char *RoiListFile, char RoiFiles[MAX_ROI_FILE_Nr][TEMP_STRING_LENGTH])
{
	int i, j;
	FILE *fp;
	if ( (fp = fopen(RoiListFile, "r")) == NULL) {
                printf("Could not open %s, exiting ...\n", RoiListFile);
                return(-1);
        }
	int num_files = 0;
	char temp_str1[TEMP_STRING_LENGTH];
	while (fgets(temp_str1, TEMP_STRING_LENGTH, fp) != NULL)
                num_files++;
        printf("\t#####Processing %d files\n", num_files);

	char filename[TEMP_STRING_LENGTH];
	fseek(fp, 0, SEEK_SET);

	for(i=0; i<num_files; i++)
	{
		fscanf(fp, "%s", filename);
		//cout << "filename = " << filename << endl;
		strcpy(RoiFiles[i], filename);
	}
	fclose(fp);

	return num_files;
}

int readRoiOperList(char *RoiOperListFile, int *RoiOperations)
{
	int i, j, k;
	FILE *fp;
	if ( (fp = fopen(RoiOperListFile, "r")) == NULL) {
                printf("Could not open %s, exiting ...\n", RoiOperListFile);
                return(-1);
        }
	int num_files = 0;
	char temp_str1[TEMP_STRING_LENGTH];
	while (fgets(temp_str1, TEMP_STRING_LENGTH, fp) != NULL)
                num_files++;
        printf("\t#####Processing %d files\n", num_files);

	char char_operation[TEMP_STRING_LENGTH];
	fseek(fp, 0, SEEK_SET);
	for(i=0; i<num_files; i++)
	{
		fscanf(fp, "%s", char_operation);
		RoiOperations[i] = atoi(char_operation);
		//cout << "ope " << i << " = " << RoiOperations[i] << endl;
	}

	return num_files;
}

char** createOptionArray( )
{
	int i;
	char **optionArray;

	optionArray = (char **)calloc(INPUTARGV_TYPE_NUM,sizeof(char*));

	for(i=0; i<INPUTARGV_TYPE_NUM; i++)
		optionArray[i] = (char *)calloc(32,sizeof(char));

	strcpy(optionArray[0], "-iParatxt");
	strcpy(optionArray[1], "-iFAmgh");
	strcpy(optionArray[2], "-iV0mgh");
	strcpy(optionArray[3], "-iStartFA");
	strcpy(optionArray[4], "-iStopFA");
	strcpy(optionArray[5], "-iTurnAngleDeg");
	strcpy(optionArray[6], "-iFlipEigenX");
	strcpy(optionArray[7], "-iFlipEigenY");
	strcpy(optionArray[8], "-iFlipEigenZ");
	strcpy(optionArray[9], "-iRoimghOR");
	strcpy(optionArray[10], "-iRoimghAND");
	strcpy(optionArray[11], "-iRoimghNOT");
	strcpy(optionArray[12], "-iRoimghCUT_OR");
	strcpy(optionArray[13], "-iRoimghCUT_AND");
	strcpy(optionArray[14], "-oAllFiberFile");
	strcpy(optionArray[15], "-oSelFiberFile");
	strcpy(optionArray[16], "-oAllFiberTxtFile");
	strcpy(optionArray[17], "-oSelFiberTxtFile");
	strcpy(optionArray[18], "-oSelFiberMghFile");
	strcpy(optionArray[19], "-iFAraw");
	strcpy(optionArray[20], "-iV0raw");
	strcpy(optionArray[21], "-iImgDims");
	strcpy(optionArray[22], "-iVoxSize");
	strcpy(optionArray[23], "-iFOVSize");
	strcpy(optionArray[24], "-iImageOrientation");
	strcpy(optionArray[25], "-iRoirawOR");
	strcpy(optionArray[26], "-iRoirawAND");
	strcpy(optionArray[27], "-iRoirawNOT");
	strcpy(optionArray[28], "-iRoirawCUT_OR");
	strcpy(optionArray[29], "-iRoirawCUT_AND");
	strcpy(optionArray[30], "-iByteSwapFlag");
	strcpy(optionArray[31], "-iFiberLenMin");
	strcpy(optionArray[32], "-iImgSequencing");
	strcpy(optionArray[33], "-iSliceThickness");
	strcpy(optionArray[34], "-oSelFiberVolFile");


	return optionArray;

}

int createOptionFlag( int argc,  char *argv[], optionStruct option_st[INPUTARGV_TYPE_NUM])
{
	//optionStruct* option_st = (optionStruct *) calloc ( INPUTARGV_TYPE_NUM, sizeof(optionStruct));
	//optionStruct option_st[INPUTARGV_TYPE_NUM];
	int i, j, k;
	int kk;
	char** optionArray = createOptionArray();
	char *tempchar;

/*
	for (j = 0; j<INPUTARGV_TYPE_NUM; j++)
		cout << optionArray[j] << endl;
	getchar();
	
	for(i = 0; i< argc; i++)
		for (j = 0; j<INPUTARGV_TYPE_NUM; j++)
			cout << "optionArray[" << j  << "] = " <<optionArray[j] << endl;
	getchar();
*/

	for (j = 0; j<INPUTARGV_TYPE_NUM; j++)
	{
		option_st[j].optionflag = 0;
		option_st[j].optiontype = 0;
	}

	for(i = 0; i< argc; i++)
	{
		if(argv[i][0] != '-')
			continue;
		else 
			for (j = 0; j<INPUTARGV_TYPE_NUM; j++)
			{
				if(strcmp(argv[i], optionArray[j]) == 0)
				{
					option_st[j].optionflag = 1;
					option_st[j].optiontype = i;
				}
			}
	} 
/*
	for (j = 0; j<INPUTARGV_TYPE_NUM; j++)
	{
		cout << j << ",  flag = " << option_st[j].optionflag << endl;
		getchar();
	}
*/
	return 1;

}

int Fiber2Mgh(FACT_MAP stFactMap, mgh_header *mghInfo)
{
	int i, j, k;

	if(stFactMap.pstFactPara->szMghAniFile == NULL)
	{
		cout << "No input mgh file, check input files." << endl;
		return -1;
	}
	int ret = read_mgh(stFactMap.pstFactPara->szMghAniFile, mghInfo, 1);
	if (ret  == -1)
	{
		cout << "Mgh file input error: check file name!" << endl;
		return -1;
	}

	mghInfo->nframes = 1;
	mghInfo->type = 1;

	int ns = stFactMap.pstFactPara->nImgWidth * stFactMap.pstFactPara->nImgHeight;
	int nv = ns * stFactMap.pstFactPara->nImgSlices;

	int *fiberVolume = (int *)calloc(nv, sizeof(int));
	for (i = 0; i<stFactMap.lFiberNrTotal; i++)
	{
		//cout << "fiberNumber = " << i << ", nLength = " << stFactMap.pstFiberAry[i].nLength << endl; 
		if(stFactMap.pstFiberAry[i].nSelStatus != 0)
		{
			for (j = 0; j<stFactMap.pstFiberAry[i].nLength; j++)
			{
				int x = (int)(stFactMap.pstFiberAry[i].pxyzChain[j].x+0.5);
				int y = (int)(stFactMap.pstFiberAry[i].pxyzChain[j].y+0.5);
				int z = (int)(stFactMap.pstFiberAry[i].pxyzChain[j].z+0.5);
				int location = x + y*mghInfo->width + z*ns;
				fiberVolume[location] = 256;
			}
		}	
	}

	mghInfo->intImageArray= (int *)calloc(nv, sizeof(int));
	for(i=0; i<nv; i++)
		mghInfo->intImageArray[i] = fiberVolume[i];

/**	
	FILE *fp;

	fp = fopen("volume.raw", "w");
	fwrite(fiberVolume, sizeof(int), nv, fp);
	fclose(fp);
	cout << "write fiber volumn!" << endl;
	getchar();
**/
	return 1;

}

int getFactPara(int argc, char* argv[], FACT_PARA *pstFactPara )
{

	//-------------------------------------------
	// initialiZIng FACT_PARA structure: pstDpfPara
	//-------------------------------------------
	
	//pstFactPara = (FACT_PARA *)calloc (1, sizeof( FACT_PARA ));

	pstFactPara->nImgWidth = -1;			// image dimension
	pstFactPara->nImgHeight = -1;
	pstFactPara->nImgSlices = -1;

	pstFactPara->nImgOrientation = 0;			// image dimension
	pstFactPara->nImgSequencing = 0;

	pstFactPara->fFovWidth = -1;			// FOV
	pstFactPara->fFovHeight = -1;
	pstFactPara->fPixelSizeWidth = -1;		// Pixel size
	pstFactPara->fPixelSizeHeight = -1;
	pstFactPara->fSliceThickness = -1;

	pstFactPara->szImgVecFile = NULL;		// input file name
	pstFactPara->szImgAniFile = NULL;
	pstFactPara->szMghVecFile = NULL;		// input file name
	pstFactPara->szMghAniFile = NULL;
	pstFactPara->bSwapBytes = false;

	pstFactPara->bFlipVecX = false;			//flip x-y-z?
	pstFactPara->bFlipVecY = false;
	pstFactPara->bFlipVecZ = false;

	pstFactPara->fStartFA = -1;			// tracking threshold values
	pstFactPara->fStopFA = -1;
	pstFactPara->fTurnDeg = -1;
	pstFactPara->nFiberLenMin = -1;

	// output files
	pstFactPara->szFiberAllFile = NULL;
	pstFactPara->szFiberSelFile = NULL;
	pstFactPara->szFiberAllTxtFile = NULL;
	pstFactPara->szFiberSelTxtFile = NULL;
	pstFactPara->szFiberSelMghFile = NULL;

	// ROI fies
	for (int i=0; i<MAX_ROI_FILE_Nr; i++)	pstFactPara->szBinRoiFile[i] = NULL;

	//-------------------------------------------
	// parsing the parameter file now
	//-------------------------------------------
	char	szLnBuf[512];
	
	int		nRoiFileCnt=0;
	int		nRoiOpCnt=0;
	int 		i, j, k;

	for(i=0; i<argc; i++)
	{
		//cout <<  "i = " << i << "; " << argv[i];
		//getchar();
		// text parameter file input
		if	(strcmp(argv[i],"-iParatxt") == 0)
		{
			int nRetCode = ReadFactParaFile(argv[2], pstFactPara);		//Read text parameter file
			if (nRetCode != DM_PROC_Ok) 
			{
                        	DispErrorMsg(DM_ERR_FileRead, (void*)argv[2], (void*)"check parameter file.");
                        	FreeMemFactPara (pstFactPara);
                        	cout << "read txt file failed. " << endl;
                        	getchar();
                        	return DM_ERR_FileRead;
                	}
                }

		// image dimensions
		else if	(strcmp(argv[i],"-iImgDims") == 0)
		{
			pstFactPara->nImgWidth = atoi(argv[i+1] );		//ImageWidth
			pstFactPara->nImgHeight = atoi(argv[i+2] );		//ImageHeight
			pstFactPara->nImgSlices = atoi(argv[i+3] );		//ImageSlices
		}

		// image orientation
		else if (strcmp(argv[i],"-iImageOrientation") == 0)
			pstFactPara->nImgOrientation = atoi(argv[i+1]) ;	//ImageOrientation: 0=coronal, 1=axial, 2=sagittal
		else if (strcmp(argv[i],"-iImageSequencing") == 0)	
			pstFactPara->nImgSequencing = atoi(argv[i+1]);		//Image Sequencing:	0=normal, 1=flipped?

		// pixel size or FOV
		else if (strcmp(argv[i],"-iVoxSize") == 0)	
		{
			pstFactPara->fPixelSizeWidth	= atof(argv[i+1]);	//Pixel size(X)
			pstFactPara->fPixelSizeHeight 	= atof(argv[i+2]);	//Pixel size(Y)
		}
		else if(strcmp(argv[i], "-iSliceThickness") == 0)
		{
			pstFactPara->fSliceThickness 	= atof(argv[i+1]);	//Slice thickness
		}
		else if (strcmp(argv[i],"-iFOVSize") == 0)	
		{
			pstFactPara->fFovWidth = atof(argv[i+1]);		//FieldOfView(X)	
			pstFactPara->fFovHeight = atof(argv[i+2]);		//FieldOfView(X)	
		}

		// user-defined threshold for fiber tracking
		else if (strcmp(argv[i],"-iStartFA") == 0)	
			pstFactPara->fStartFA = atof(argv[i+1]);	//tracking start FA
		else if (strcmp(argv[i],"-iStopFA") == 0)
			pstFactPara->fStopFA = atof(argv[i+1]);		//tracking stop FA
		else if (strcmp(argv[i],"-iTurnAngleDeg")==0) 
			pstFactPara->fTurnDeg = atof(argv[i+1]);	//tracking turning angle
		else if (strcmp(argv[i],"-iFiberLenMin")==0)   
			pstFactPara->nFiberLenMin = atof(argv[i+1]);	//minimum fiber length
		
		// data manipulation
		else if (strcmp(argv[i],"-iByteSwapFlag") == 0) {
			pstFactPara->bSwapBytes = (bool)atoi(argv[i+1]);
		}

		else if (strcmp(argv[i],"-iFlipEigenX") == 0) {
			pstFactPara->bFlipVecX = (bool)atoi(argv[i+1]);
		}

		else if (strcmp(argv[i],"-iFlipEigenY") == 0) {
			pstFactPara->bFlipVecY = (bool)atoi(argv[i+1]);
		}

		else if (strcmp(argv[i],"-iFlipEigenZ") == 0) {
			pstFactPara->bFlipVecZ = (bool)atoi(argv[i+1]);
		}
		
		// input file name
		else if (strcmp(argv[i],"-iFAraw") == 0) {		// FA image file
			pstFactPara->szImgAniFile = strdup( argv[i+1] );
		}

		else if (strcmp(argv[i],"-iV0raw") == 0) {		// Vector image file
			pstFactPara->szImgVecFile = strdup( argv[i+1] );
		}

		else if (strcmp(argv[i],"-iFAmgh") == 0) {		// FA mgh format image file
			pstFactPara->szMghAniFile = strdup( argv[i+1] );
		}

		else if (strcmp(argv[i],"-iV0mgh") == 0) {		// Vector mgh format image file
			pstFactPara->szMghVecFile = strdup( argv[i+1] );
		}


		// output file name
		else if (strcmp(argv[i], "-oAllFiberFile") == 0) {		// output: all fibers
			pstFactPara->szFiberAllFile= strdup( argv[i+1] );
		}
		else if (strcmp(argv[i],"-oSelFiberFile") == 0) {		// output: selected fibers
			pstFactPara->szFiberSelFile= strdup(argv[i+1] );
		}
		else if (strcmp(argv[i],"-oAllFiberTxtFile") == 0) {		// text file output: all fibers
			pstFactPara->szFiberAllTxtFile= strdup(argv[i+1] );
		}
		else if (strcmp(argv[i], "-oSelFiberTxtFile") == 0) {		// text file output: selected fibers
			pstFactPara->szFiberSelTxtFile= strdup(argv[i+1] );
		}
		else if (strcmp(argv[i], "-oSelFiberMghFile") == 0) {		// text file output: selected fibers
			pstFactPara->szFiberSelMghFile = strdup(argv[i+1] );
		}
		else if (strcmp(argv[i], "-oSelFiberVolFile") == 0) {		// text file output: selected fibers
			pstFactPara->szFiberSelVolFile = strdup(argv[i+1] );
		}

		// binary ROI files
		else if (strcmp(argv[i],"-iRoimghOR") == 0) {		// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  1 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpOr ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoimghAND") == 0) {		// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  1 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpAnd ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoimghNOT") == 0) {		// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  1 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpNot ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoimghCUT_OR") == 0) {	// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  1 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpCutOr ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoimghCUT_AND") == 0) {	// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  1 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpCutAnd ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoirawOR") == 0) {	// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  0 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpOr ;
			nRoiFileCnt++;
		}
		else if (strcmp(argv[i],"-iRoirawAND") == 0) {		// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  0 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpAnd ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoirawNOT") == 0) {		// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] = 0;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpNot ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoirawCUT_OR") == 0) {	// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  0 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpCutOr ;
			nRoiFileCnt++;
		}

		else if (strcmp(argv[i],"-iRoirawCUT_AND") == 0) {	// ROI file, for fiber selection
			pstFactPara->szBinRoiFile[nRoiFileCnt] =  strdup( argv[i+1] ); 
			pstFactPara->szRoiFileType[nRoiFileCnt] =  0 ;
			pstFactPara->nBinRoiOp[nRoiFileCnt] =  ROI_OpCutAnd ;
			nRoiFileCnt++;
		}
	}

	//--------------------------------------------------------------
	// Validating the parameters, are they all ready? 
	//--------------------------------------------------------------
	// REC file(s)
	if ( ((pstFactPara->szImgVecFile == NULL) || (pstFactPara->szImgAniFile == NULL)) &&
	     ((pstFactPara->szMghVecFile == NULL) && (pstFactPara->szMghAniFile == NULL))) {
		DispErrorMsg(DM_ERR_FileRead, (void*)pstFactPara->szImgVecFile, (void*)" fail to find  vector or anisotropy image file name");
		return	DM_ERR_FileRead;
	}
	if (((pstFactPara->szMghVecFile == NULL) || (pstFactPara->szMghAniFile == NULL)) &&
	    ((pstFactPara->szImgVecFile == NULL) && (pstFactPara->szImgAniFile == NULL))) {
		DispErrorMsg(DM_ERR_FileRead, (void*)pstFactPara->szImgVecFile, (void*)" fail to find mhg vector or anisotropy mgh image file name");
		return	DM_ERR_FileRead;
	}

	// image dimension
	if (  ((pstFactPara->szImgVecFile != NULL) && (pstFactPara->szImgAniFile != NULL)) &&
	      ((pstFactPara->nImgWidth < 0) || (pstFactPara->nImgHeight < 0) || (pstFactPara->nImgSlices < 0))) {
		DispErrorMsg(DM_ERR_FileRead, (void*)NULL, (void*)" forget image dimensions (ImageWidth/Height/Slices)?");
		return DM_ERR_FileRead;
	}
	
	// FOV or Pixel Size
	if (  ((pstFactPara->szImgVecFile != NULL) && (pstFactPara->szImgAniFile != NULL)) &&
		(((pstFactPara->fFovWidth < 0) && (pstFactPara->fPixelSizeWidth < 0)) || 
		((pstFactPara->fFovHeight < 0) && (pstFactPara->fPixelSizeHeight < 0)) ||
		( pstFactPara->fSliceThickness <= 0) )) {
		DispErrorMsg(DM_ERR_FileRead, (void*)NULL, (void*)" forget FOV or pixel-size (FieldOfView/PixelSize/SliceThickness)?");
		return DM_ERR_FileRead;
	}

	// tracking threshold
	if ((pstFactPara->fStartFA < 0) || (pstFactPara->fStopFA < 0) || (pstFactPara->fTurnDeg < 0)) {
		DispErrorMsg(DM_ERR_FileRead, (void*)NULL, (void*)" forget tracking threshold values (TrackingStartFA/StopFA/TractTurningAngle)?");
		return DM_ERR_FileRead;
	}
	
	//--------------------------------------------------------------
	// regulate parameters
	//--------------------------------------------------------------
	// FOV
	if ((pstFactPara->fFovWidth < 0) && (pstFactPara->fPixelSizeWidth > 0)) 	// your input is pixel size
		pstFactPara->fFovWidth = pstFactPara->nImgWidth * pstFactPara->fPixelSizeWidth;

	if ((pstFactPara->fFovHeight < 0) && (pstFactPara->fPixelSizeHeight > 0)) 	// your input is pixel size
		pstFactPara->fFovHeight = pstFactPara->nImgHeight * pstFactPara->fPixelSizeHeight;

	// calculate voxel size from FOV
	if (  ((pstFactPara->szImgVecFile != NULL) && (pstFactPara->szImgAniFile != NULL)) &&
	      ((pstFactPara->fFovWidth 		> 0) && (pstFactPara->nImgWidth 	> 0)) && 
	      ((pstFactPara->fFovHeight 	> 0) && (pstFactPara->nImgHeight 	> 0)) &&
	      ((pstFactPara->fPixelSizeWidth 	< 0) && (pstFactPara->fPixelSizeHeight 	< 0)) )
	{
		pstFactPara->fPixelSizeWidth = pstFactPara->fFovWidth / pstFactPara->nImgWidth;
		pstFactPara->fPixelSizeHeight = pstFactPara->fFovHeight / pstFactPara->nImgHeight;
	}

	// calculate FOV from voxel size	
	if (  ((pstFactPara->szImgVecFile != NULL) && (pstFactPara->szImgAniFile != NULL)) &&
	      ((pstFactPara->fPixelSizeWidth 	> 0) && (pstFactPara->nImgHeight 	> 0)) &&
	      ((pstFactPara->fPixelSizeHeight 	> 0) && (pstFactPara->nImgHeight 	> 0)) &&
	      ((pstFactPara->fFovWidth 		< 0) && (pstFactPara->fFovHeight 	< 0)))  
	{
		pstFactPara->fFovWidth  = pstFactPara->fPixelSizeWidth  * pstFactPara->nImgWidth;
		pstFactPara->fFovHeight = pstFactPara->fPixelSizeHeight * pstFactPara->nImgHeight;
	}

	return DM_PROC_Ok;
}
int writeRaw(FACT_MAP stFactMap)
{
	int i, j, k;

	int ns = stFactMap.pstFactPara->nImgWidth * stFactMap.pstFactPara->nImgHeight;
	int nv = ns * stFactMap.pstFactPara->nImgSlices;

	int *fiberVolume = (int *)calloc(nv, sizeof(int));
	for (i = 0; i<stFactMap.lFiberNrTotal; i++)
	{
		//cout << "fiberNumber = " << i << ", nLength = " << stFactMap.pstFiberAry[i].nLength << endl; 
		if(stFactMap.pstFiberAry[i].nSelStatus != 0)
		{
			for (j = 0; j<stFactMap.pstFiberAry[i].nLength; j++)
			{
				int x = (int)(stFactMap.pstFiberAry[i].pxyzChain[j].x+0.5);
				int y = (int)(stFactMap.pstFiberAry[i].pxyzChain[j].y+0.5);
				int z = (int)(stFactMap.pstFiberAry[i].pxyzChain[j].z+0.5);
				int location = x + y*stFactMap.pstFactPara->nImgWidth + z*ns;
				fiberVolume[location] = 256;
			}
		}	
	}


	FILE *fp;

	if( (fp = fopen(stFactMap.pstFactPara->szFiberSelVolFile, "w")) == NULL)
	{
		cout << "Can not open volume file: " << stFactMap.pstFactPara->szFiberSelVolFile << endl;
		return 0;
	}
	fwrite(fiberVolume, sizeof(int), nv, fp);
	fclose(fp);

	return 1;
}
