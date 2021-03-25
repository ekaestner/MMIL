//==================================================================
// Constant definition
//==================================================================
// maximum number of ROI files
#define		MAX_ROI_FILE_Nr		16

// debug compiling instructive
#define		_MyDebug
//#undef		_MyDebug

#define DEFAULT_STARTFA      0.2
#define DEFAULT_STOPFA       0.2
#define DEFAULT_TURNDEG      50
#define DEFAULT_FIBERMINLEN  7
#define DEFAULT_FIBERMAXLEN  4196

#define MAX_INPUTARGV_NUM 100 
#define INPUTARGV_TYPE_NUM 36

//==================================================================
// ROI operation code
//==================================================================
#define		ROI_OpOr			0
#define		ROI_OpAnd			1
#define		ROI_OpXor			2
#define		ROI_OpNot			3
#define		ROI_OpCutOr			4
#define		ROI_OpCutAnd			5

//==================================================================
// error code
//==================================================================
#define		DM_PROC_Ok			0						// the procedure is successful

#define		DM_ERR_CodeBase		0x7000
#define		DM_ERR_Usage		DM_ERR_CodeBase + 0		// Para = NULL

// General file I/O errors
#define		DM_ERR_FileOpen		DM_ERR_CodeBase + 1		// Para0 = FileName
#define		DM_ERR_FileRead		DM_ERR_CodeBase + 2		// Para0 = FileName [,para0=ErrMsg]
#define		DM_ERR_MemAlloc		DM_ERR_CodeBase + 5		// Para0 = ErrMsg

#define		DM_ERR_Unknow		DM_ERR_CodeBase + 700

//==================================================================
// macro functions
//==================================================================
//#define max(a,b) (((a)>(b))?(a):(b)) 
//#define min(a,b) (((a)<(b))?(a):(b)) 

//==================================================================
// new data type
//==================================================================
//#define		BYTE	unsigned char
#define		WORD	unsigned short
#define		DWORD	unsigned int

// vector data type
typedef struct tagXYZ_TRIPLE {
	float	x, y, z;
} XYZ_TRIPLE;

// R-G-B data type
typedef struct tagRGB_TRIPLE {
	unsigned char	r, g, b;
} RGB_TRIPLE;

// ROI data information, 2007-01-25 by sumiko, UCSD
typedef struct roiInfo {
	int *ROI_points;
	int ROI_size;
} roiInfo;
//==================================================================
// Data structure
//==================================================================
// fiber data
typedef struct tagFIBER {
	int			nLength;			// fiber length
	XYZ_TRIPLE*	pxyzChain;			// pointer of fiber data chain
	unsigned char		nSelStatus;			// fiber has been selected (>= 1) or deselected (== 0).
	int			nSelBeginIdx;		// start index of the selected fiber-segment
	int			nSelEndIdx;			// end index of the selected fiber-segment
	int			nSelCnt;			// number of being selected by ROI
	XYZ_TRIPLE	xyzPtSeed;			// seed point
	XYZ_TRIPLE	xyzPtROI;			// ROI point where the fiber been selected
} FIBER;

// parameter
typedef struct tagFACT_Para {
	// image diemnsion
	int		nImgWidth;			// image diemnsion
	int		nImgHeight;
	int		nImgSlices;

	int		nImgOrientation;
	int		nImgSequencing;

	// FOV or pixel size
	float	fFovWidth;			// FOV
	float	fFovHeight;
	float	fPixelSizeWidth;	// Pixel size
	float	fPixelSizeHeight;
	float	fSliceThickness;

	// input file name
	int	mghFileInput[2];
	char	*szImgVecFile;		// principal vector image 
	char	*szImgAniFile;		// anisotropy image 
	char	*szMghVecFile;		// MGH format principal vector image 
	char	*szMghAniFile;		// MGH format anisotropy image 

	bool	bSwapBytes;			
	bool	bFlipVecX, bFlipVecY, bFlipVecZ;

	// tracking threshold values
	float	fStartFA;
	float	fStopFA;
	float	fTurnDeg;
	int	nFiberLenMin;		// minimum fiber length
	int	nFiberLenMax;		// maximum fiber length

	// fiber selection by ROI
	char*	szBinRoiFile[MAX_ROI_FILE_Nr];	// binary ROI files
	// fiber selection by ROI, mgh format, add by sumiko 02-14-07 UCSD
	int	szRoiFileType[MAX_ROI_FILE_Nr];	// roi file format:    0: raw, 1:mgh
	int	nBinRoiOp[MAX_ROI_FILE_Nr];	// binary operation

	// Optional Output files name
	char	*szFiberAllFile;
	char	*szFiberSelFile;

	char	*szFiberAllTxtFile;
	char	*szFiberSelTxtFile;
	char	*szFiberSelMghFile;
	char	*szFiberSelVolFile;

} FACT_PARA;

// DTI maps
typedef struct tagFACT_Map {
	FACT_PARA*	pstFactPara;				// image parameters

	float*		pfAniImg;					// original FA map
	XYZ_TRIPLE*	pxyzVecImg;					// original Vector image

	FIBER*		pstFiberAry;				// fiber data array
	int**		pnFiberIdxAry;				// fiber data index of each voxel, 
	int*		pnFiberCntAry;				// number of fibers through each voxel

	// fiber data statistics
	long		lFiberNrTotal;				// number of fibers in the fiber data array (array size)
	long		lFiberLenMax;				// maximub length of all fibers
	float		fFiberLenMean;				// mean length of all fibers
	float		fFiberLenStdDev;			// standard derivation of fiber length
} FACT_MAP;

typedef struct optionStruct {
	int optionflag;
	int optiontype;
} optionStruct;

//==================================================================
// Function protocal
//==================================================================
// main entry
bool	DispErrorMsg(int nErrCode=DM_ERR_Unknow, void* pPara0=NULL, void* pPara1=NULL, void* pPara2=NULL, void* pPara3=NULL);
void	FreeMemFactPara(FACT_PARA*	pstFactPara);						// release memory for FACT_PARA
void	FreeMemFactMap(FACT_MAP*	pstFactMap);						// release memory for FACT_MAP

// read parameter file
int		ReadLn(FILE *inFile, char* pBuf, int nBufLen);						// read a line from text file
int		ReadFactParaFile(char* szParaFileName,	FACT_PARA*	pstFactPara);	// parsing the parameter file

// read original image
int		ReadImgFiles(FACT_MAP*	pstFactMap);			// read files and reconstruct the images
void	SwapFloatImgBytes(float* pfImg, long lImgSize); // swap float image bytes

// fiber tracking
int FiberTrackingFACT(FACT_MAP*	pstFactMap);
int RoiFiberTrackingFACT(FACT_MAP*	pstFactMap, roiInfo RoiInfo);
roiInfo FindRoiMask(FACT_MAP*      pstFactMap, unsigned char *ROI);
XYZ_TRIPLE* GetFiberChainFact(float* pfImgFA, XYZ_TRIPLE* pxyzImgVec, int nImgWidth, int nImgHeight, int nImgSlices,
							  XYZ_TRIPLE xyzSeed, float fStopFA, float fTurnDeg, int nFiberLenMin, int nFiberLenMax, int* pnFiberLen);
int	GetHalfChainFact(float* pfImgFA, XYZ_TRIPLE* pxyzImgVec, int nImgWidth, int nImgHeight, int nImgSlices,
							  XYZ_TRIPLE xyzSeedPt, float fStopFA, float fTurnDeg, int nFiberLenMin, int nFiberLenMax,
							  int	nTrackDir, XYZ_TRIPLE xyzFiberPt[]);
float CosAlphaVec(XYZ_TRIPLE vecA, XYZ_TRIPLE vecB);
XYZ_TRIPLE GetLineCubeInterPt(XYZ_TRIPLE xyzPt0, XYZ_TRIPLE xyzDir, float* pfStepLen);
int CreateFiberIndexAry(FACT_MAP*	pstFactMap);

// fiber selection
int FiberSelectByRoiMap(FACT_MAP*	pstFactMap);
void LabelFiberViaRoiImg(FIBER	stFiberAry[], int nFiberArySize, int* nFiberIdxAry[], int nFiberCntAry[],
						float* pByteImg,	int nImgWidth, int nImgHeight, int nImgSlices, int nRoiOp);

// fiber save
int FiberSave(FACT_MAP*	pstFactMap, bool bSaveSel=true);
int FiberSaveTxt(FACT_MAP*	pstFactMap, bool bSaveSel=true);

