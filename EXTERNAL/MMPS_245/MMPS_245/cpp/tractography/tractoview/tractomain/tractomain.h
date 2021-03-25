#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>
#include <qaction.h>
#include <qapplication.h>
#include <qcombobox.h>
#include <qfile.h>
#include <qfiledialog.h>
#include <qfont.h>
#include <qfontdialog.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qpixmap.h>
#include <qpopupmenu.h>
#include <qprinter.h>
#include <qradiobutton.h>
#include <qsettings.h>
#include <qspinbox.h>
#include <qstatusbar.h>
#include <qtoolbar.h>
#include <qtoolbutton.h>
#include <qlabel.h>
#include <qcheckbox.h>
#include <qlistbox.h>
#include <qtable.h>
#include <qlistview.h>
#include <qbuttongroup.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qtextedit.h>
#include <qnamespace.h>
#include <qsplitter.h>
#include <assert.h>
#include <stdio.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/mui.h>
#include <Vec3.h>
#include <Vec2.h>
#include <ColorRGBA.h>
#include <Matrix4f.h>
#include <Material.h>
#include <Light.h>
#include <Mouse.h>
#include <Keyboard.h>
#include <HyperBall.h>
#include <dirent.h>
#include <iostream>
#include <new>
#include <string>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cctype>
#include <tiffio.h>
#include <ContextData.h>
#include "mgh.h"

using namespace std;

/* generic constants and macros */
#define STRLEN           1000
#define MATCH(A,B)       (!strcmp(A,B))
#define RATIO(X,Y)        ((float)X/(float)Y)
#ifdef Darwin
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),INT_MAX)*((H)-(L)))
#else
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),MAXINT)*((H)-(L)))
#endif
#define SQR(X)            ((X)*(X))
#define DIST(X0,Y0,X1,Y1) (sqrt(SQR(X1-X0) + SQR(Y1-Y0)))
#define BIGF              10e9
#define MTEST(ptr) \
  if((ptr)==NULL) ( fprintf(stderr,"*** Memory allocation error\n"), exit(0))
#define FALSE 0
#define TRUE 1

/* constants */
#define DEFAULT_WIN_WIDTH     1024
#define DEFAULT_WIN_HEIGHT    1024
#define MAX_FIBER_FILES       1000
#define PI                    3.1415925 
#define STEPNUM               359 
#define AXIAL                 0
#define SAGITTAL              1
#define CORONAL               2
#define DEFAULT_WIDTH         1.0
#define DEFAULT_LEVEL         0.5
#define MAX_ROI_FILE_Nr       16

// todo: is a "base" application really necessary?

/*****************************************************************/
// from GlutKernelApplication.h

// vector data type
typedef struct tagXYZ_TRIPLE {
        float   x, y, z;
} XYZ_TRIPLE;

// R-G-B data type
typedef struct tagRGB_TRIPLE {
        unsigned char   r, g, b;
} RGB_TRIPLE;

typedef struct stFiberFileHeader{
       char    sFiberFileTag[8];       // file tag = FiberDat
       int     nFiberNr;               // total number of fibers
       int     nFiberLenMax;           // max-length of fibers
       float   fFiberLenMean;          // mean-length of fibers
       int     nImgWidth;              // image dimension
       int     nImgHeight;
       int     nImgSlices;
       float   fPixelSizeWidth;        // voxel size
       float   fPixelSizeHeight;
       float   fSliceThickness;
       int                     enumSliceOrientation;   // slice orientation
       int                     enumSliceSequencing;    // slice sequencing
} stFiberFileHeader;


// fiber data
typedef struct tagFIBER {
        int             nLength;                        // fiber length
        XYZ_TRIPLE*     pxyzChain;                      // pointer of fiber data chain
        unsigned char   nSelStatus;                     // fiber has been selected (>= 1) or deselected (== 0).
        int             nSelBeginIdx        ;           // start index of the selected fiber-segment
        int             nSelEndIdx;                     // end index of the selected fiber-segment
        int             nSelCnt;                        // number of being selected by ROI
        RGB_TRIPLE      rgbFiberClr;
        XYZ_TRIPLE      xyzPtSeed;                      // seed point
        XYZ_TRIPLE      xyzPtROI;                       // ROI point where the fiber been selected
} FIBER;

// parameter
typedef struct tagFACT_Para {
        // image diemnsion
        int             nImgWidth;                      // image diemnsion
        int             nImgHeight;
        int             nImgSlices;

        int             nImgOrientation;
        int             nImgSequencing;

        // FOV or pixel size
        float   fFovWidth;                      // FOV
        float   fFovHeight;
        float   fPixelSizeWidth;        // Pixel size
        float   fPixelSizeHeight;
        float   fSliceThickness;

        // input file name
        int     mghFileInput[2];
        char    *szImgVecFile;          // principal vector image 
        char    *szImgAniFile;          // anisotropy image 
        char    *szMghVecFile;          // MGH format principal vector image 
        char    *szMghAniFile;          // MGH format anisotropy image 

        bool    bSwapBytes;
        bool    bFlipVecX, bFlipVecY, bFlipVecZ;

        // tracking threshold values
        float   fStartFA;
        float   fStopFA;
        float   fTurnDeg;
        int     nFiberLenMin;           // minimum fiber length

        // fiber selection by ROI
        char*   szBinRoiFile[MAX_ROI_FILE_Nr];  // binary ROI files
        // fiber selection by ROI, mgh format, add by sumiko 02-14-07 UCSD
        int     szRoiFileType[MAX_ROI_FILE_Nr]; // roi file format:    0: raw, 1:mgh
        int     nBinRoiOp[MAX_ROI_FILE_Nr];     // binary operation

        // Optional Output files name
        char    *szFiberAllFile;
        char    *szFiberSelFile;

        char    *szFiberAllTxtFile;
        char    *szFiberSelTxtFile;
        char    *szFiberSelMghFile;
        char    *szFiberSelVolFile;

} FACT_PARA;


typedef struct imageStruct
{
  int width;
  int height;
  int depth;
  float voxel_w, voxel_h, voxel_t;
  int Orientation; // 0: axial, 1: coronal, 2: sagittal
  float *ImageArray;
  int AxialSliceLocation ;
  int CoronalSliceLocation ;
  int SagittalSliceLocation;
  int imageType;
  int AsliceLoc, CsliceLoc, SsliceLoc;
  int MaxAsliceLoc, MaxCsliceLoc, MaxSsliceLoc;
  float VoxelSize[3];
  int slideRedraw, lineFlag, mouseFlag;
  float AxialLineH ;
  float AxialLineL;
  float CoronalLineH;
  float CoronalLineL;
  float SagittalLineH;
  float SagittalLineL;

  char *FAmghFile;
  char *V0mghFile;
  char *FArawFile;
  char *V0rawFile;

  int fiberFileNumber;
  char **FiberFileNameArray;

  int fiberNum;
  int *totalFiberNum;
  RGB_TRIPLE *FiberColor;
  FIBER **pasFiberAry;
} imageStruct;

static imageStruct imageinfo;
static imageStruct temp_imageinfo;

//  make sure you remember to initialize the static App::_instance = new MyApp
class GlutKernelApplication
{
// Application attributes.
public:
	inline void              setWidth( const int& w ) { _width = w; }
	inline void              setHeight( const int& h ) { _height = h; }
	inline int				       width() const { return _width; }
	inline int               height() const { return _height; }
	inline const char* const name() const { return _name.data(); }
	inline void              setName( const std::string& name) { _name = name; }
	inline void              setName( const char* const name) { _name = name; }
	inline const Mouse&      mouse() const { return _mouse; }
	inline const Keyboard&   keyboard() const { return _keyboard; }
	inline Keyboard&         keyboard() { return _keyboard; }
	inline const char&       modifier() const { return _keyboardModifier; }

// 
public:
	//virtual void OnAppDraw() {}
	virtual void OnContextDraw() {}

  //: Do your frame calculations here.
	//  NOTE: Called every frame before draw
	virtual void OnAppIdle();

// Init methods
public:
	virtual void OnAppInit(int argc = 0, char *argv[] = NULL) {};
	virtual void AxialSliceInit(int argc, char *argv[] ) ;
	virtual void AxialDisplay( ) ;
	virtual void CoronalSliceInit(int argc, char *argv[] ) ;
	virtual void CoronalDisplay( ) ;
	virtual void SagittalSliceInit(int argc, char *argv[] ) ;
	virtual void SagittalDisplay( ) ;
	virtual void OperationPanelInit(int argc, char *argv[] ) ;
	virtual void OperationPanelEvent( ) ;
	virtual void OperationPanelOpenFilesInit(int argc, char *argv[] ) ;
	virtual void OperationPanelOpenFilesEvent( ) ;
	virtual void OnContextInit() {};

// Reshape methods
public:
	//: Window reshape
	//  Called on reshape of a window.
	//  Put application code in here.
	virtual void OnAppReshape(int width, int height);

	//: Window reshape for each GL context (window)
	//  Called on reshape of each gl window.
	//  Put GL code in here.
	virtual void OnContextReshape(int width, int height) {}

// Keyboard interaction methods
public:
	virtual void OnKeyboardDown(unsigned char k, int x, int y);
	virtual void OnKeyboardUp(unsigned char k, int x, int y);
	virtual void OnSpecialKeyboardDown(int k, int x, int y);
	virtual void OnSpecialKeyboardUp(int k, int x, int y);

// Mouse interaction methods
public:
	virtual void OnMouseEvent() {}
	virtual void OnMouseEvent4() {}
	virtual void OnMousePos( int x, int y );	// This is called when mouse changes position
	virtual void OnMouseClick( int a, int b, int c, int d );	// This is called when mouse clicks
	virtual void AxialMouseClick( int a, int b, int c, int d );
	virtual void CoronalMouseClick( int a, int b, int c, int d );
	virtual void SagittalMouseClick( int a, int b, int c, int d );
	virtual void OnMouseClick2( int a, int b, int c, int d );
	virtual void OnMouseClick3( int a, int b, int c, int d );
	virtual void OnMouseClick4( int a, int b, int c, int d );

	float AxialLineH ;
	float AxialLineL ;
	float CoronalLineH ;
	float CoronalLineL ;
	float SagittalLineH ;
	float SagittalLineL ;
	float AxialCurrentWidth, AxialCurrentHeight;
	float CoronalCurrentWidth, CoronalCurrentHeight;
	float SagittalCurrentWidth, SagittalCurrentHeight;

	int AxialMoustPositionX, AxialMoustPositionY;
	int CoronalMoustPositionX, CoronalMoustPositionY;
	int SagittalMoustPositionX, SagittalMoustPositionY;
	
	float AxialSliceLocation, CoronalSliceLocation, SagittalSliceLocation;

  imageStruct	imageinfo;

// Internal application data
private:
	int					_width, _height;
	Mouse				_mouse;
	Keyboard				_keyboard;
	std::string			_name;
	char				_keyboardModifier;

// registration
public:
	int 				imageType;
	int FileMenu, EditMenu, HelpMenu;
	muiObject *b1, *b2, *b3, *pd;

	float xx , yy , zz ;
	float cx , cy , cz ;
	float rr , gg , bb ;
	int sphereControl ;

	muiObject *Aslider, *Cslider, *Sslider;
	muiObject *Alabel, *Clabel, *Slabel;
	muiObject *noSphereB, *smallSphereB, *largeSphereB;

   // list of all applications in the system...
   static std::vector<GlutKernelApplication*>& applications()
   {
      static std::vector<GlutKernelApplication*> mRegisteredApplications;
      return mRegisteredApplications;
   }

   // create an instance of this type to register your application
   // delete it to unregister it.
   template< class applicationType >
   class Register
   {
   public:
      Register() : mApplication()
      {
         GlutKernelApplication::applications().push_back( &mApplication );
         int size = GlutKernelApplication::applications().size();
      }
      ~Register()
      {
         std::vector<GlutKernelApplication*>::iterator it;
         for (it = GlutKernelApplication::applications().begin(); 
              it != GlutKernelApplication::applications().end(); 
              ++it)
         {
            if (&mApplication == (*it))
            {
               GlutKernelApplication::applications().erase( it );
               return;
            }
            assert( false && "the application wasn't registered" );
         }
      }
      
      applicationType mApplication;
   };
};


/*****************************************************************/

class tractomain : public GlutKernelApplication
{
public:
   tractomain() : trackBall( &hyperBall )
   {
   }

public:
	virtual void OnContextDraw();

// Init methods
public:
	virtual void OnAppInit(int argc = 0, char *argv[] = NULL);
	virtual void OnContextInit();

// Reshape methods
public:
	//: Window reshape
	//  Called on reshape of a window.
	//  Put application code in here.
	virtual void OnAppReshape(int width, int height);

	//: Window reshape for each GL context (window)
	//  Called on reshape of each gl window.
	//  Put GL code in here.
	virtual void OnContextReshape(int width, int height);

// Keyboard interaction methods
public:
	virtual void OnKeyboardDown(unsigned char k, int x, int y);
	virtual void OnKeyboardUp(unsigned char k, int x, int y);
	virtual void OnSpecialKeyboardDown(int k, int x, int y);
	virtual void OnSpecialKeyboardUp(int k, int x, int y);
// Mouse interaction methods
public:
	virtual void OnMouseEvent();
	virtual void OnMouseEvent2();
	virtual void OnMouseEvent3();
	virtual void OnMouseEvent4();
	virtual void OnMousePos( int x, int y );
	virtual void OnMouseClick( int a, int b, int c, int d );
	virtual void OnMouseClick2( int a, int b, int c, int d );
	virtual void OnMouseClick3( int a, int b, int c, int d );
	virtual void OnMouseClick4( int a, int b, int c, int d );
	virtual void AxialMouseClick( int a, int b, int c, int d );
	virtual void CoronalMouseClick( int a, int b, int c, int d );
	virtual void SagittalMouseClick( int a, int b, int c, int d );
	int ReadMghFile(const char *mghFileName,  mgh_header *mghInfor, bool endianInfor);
	void bswap( char* p1, char* p2 );
	void bswap( short& n );
	void bswap( int& n );
	void bswap( long& n );
	void bswap( float& n );
	void bswap( double& n );
	FIBER* ReadFiber(string fname, int *numfibers);
	void AxCoSaRedraw( );
	void OperationPanelEvent();
	void AxialDisplay(void);
	void CoronalDisplay(void);
	void SagittalDisplay(void);
	imageStruct readfiles(int argc, char *argv[]);
	int ChangeFiberColors( );

private:
   void drawText();

   int columns, rows, slices;
   float xPixelSpacing, zSliceThickness;
   short *image;
   float *image_float;
   float **imageVector_float;
   int imageType;
   int ns, np, nv;
   float WW; 
   float WL;
   char **optionArray;
   char *colorPalletFile;

   static const GLfloat nearPlane =5.0f;
   static const GLfloat farPlane  = 2000.0f;
   static const GLfloat nearPlane2 = 60.0f;
   static const GLfloat farPlane2  = 100.0f;

   //interaction objects
   HyperBall hyperBall;
   TrackBall *trackBall;

   Light light, light2, light3;
   int fontBase;

   class TextBeacon
   {
   public:
       char text[256];
       int timeLeft;  
   };

   TextBeacon beacon;
};


/*****************************************************************/
/* function prototypes */
void parseopts(int argc, char **argv);
void usage();
int FileExists(const char *fname);
void getAxialSliceLocation( muiObject *obj, enum muiReturnValue rv );
void getCoronalSliceLocation( muiObject *obj, enum muiReturnValue rv );
void getSagittalSliceLocation( muiObject *obj, enum muiReturnValue rv );
void getXMovieMotion( muiObject *obj, enum muiReturnValue rv );
void getXReverseMovieMotion( muiObject *obj, enum muiReturnValue rv );
void getYMovieMotion( muiObject *obj, enum muiReturnValue rv );
void getYReverseMovieMotion( muiObject *obj, enum muiReturnValue rv );
void getZMovieMotion( muiObject *obj, enum muiReturnValue rv );
void getZReverseMovieMotion( muiObject *obj, enum muiReturnValue rv );
void getSDirection( muiObject *obj, enum muiReturnValue rv );
void getIDirection( muiObject *obj, enum muiReturnValue rv );
void getRDirection( muiObject *obj, enum muiReturnValue rv );
void getLDirection( muiObject *obj, enum muiReturnValue rv );
void getPDirection( muiObject *obj, enum muiReturnValue rv );
void getADirection( muiObject *obj, enum muiReturnValue rv );
void flipFiberXDir( muiObject *obj, enum muiReturnValue rv );
void flipFiberYDir( muiObject *obj, enum muiReturnValue rv );
void flipFiberZDir( muiObject *obj, enum muiReturnValue rv );
void changeContrast( muiObject *obj, enum muiReturnValue rv );
void changeBrightness( muiObject *obj, enum muiReturnValue rv );
void getxrot( muiObject *obj, enum muiReturnValue rv );
void getyrot( muiObject *obj, enum muiReturnValue rv );
void getzrot( muiObject *obj, enum muiReturnValue rv );
void optionsSetData2(muiObject *obj, enum muiReturnValue rv );
void glRenderImage(float* dataT2W, int originalMode, 
  float xResolution, float yResolution, float zResolution,
  int width, int height, int slices, float WL, float WW,
  float axialSliceLoc, float sagitalSliceLoc, float coronalSliceLoc,
  int axialflag, int sagitalflag, int coronalFlag, 
  int totalFiberSheafNum, int *totalFiberNum, RGB_TRIPLE *fiberColor, FIBER** Fibers, 
  float axialTranspa, float coronalTranspa, float sagittalTranspa, 
  int xflipflag, int yflipflag, int zflipflag, float ****fiberPointRGB, 
  int opacityFlag);
void glClearViewport( float c1R, float c1G, float c1B, 
                      float c2R, float c2G, float c2B, bool glclear = true );
void glRender( const Light& light );
void changeView(string view);
int save_tiff (const char* fname);

/*****************************************************************/
/* from glutMain: functions for slices viewing */
void AxialSliceInit(int argc, char *argv[]);
void CoronalSliceInit(int argc, char *argv[]);
void SagittalSliceInit(int argc, char *argv[]);
void AxialDisplay(void);
void CoronalDisplay(void);
void SagittalDisplay(void);

void bswap( char* p1, char* p2 );
void bswap( short& n );
void bswap( int& n );
void bswap( float& n );
void bswap( long& n );
void bswap( double& n );

