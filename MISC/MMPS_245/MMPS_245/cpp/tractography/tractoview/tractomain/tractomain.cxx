/* tractomain.cxx

fiber viewer

created:  05/20/07 Sumiko Abe
last mod: 02/23/08 Don Hagler
last mod: 10/27/08 Don Hagler

*/

// todo: option for offscreen render
// todo: options for window width and height
// todo: with savetiffflag=1, don't draw operation panel or slice views


#include "tractomain.h"

/* variables */
static char *progname = NULL;
int nfibers = 0;
muiObject *b1, *b2, *b3, *b4, *b5, *b6;
muiObject *Aslider, *Cslider, *Sslider;
muiObject *ContrastSlider, *BrightnessSlider;
muiObject *Alabel, *Clabel, *Slabel;
muiObject *SRadio, *IRadio, *RRadio, *LRadio, *ARadio, *PRadio;
muiObject *XFlipRadio, *YFlipRadio, *ZFlipRadio;
muiObject *FiberViewEffectRadio, *FiberViewOpacityRadio;
muiObject *xrotTB, *yrotTB, *zrotTB;
static float maxDist = 0;
int stepNum = 0;
float rotationStep[STEPNUM];
int InitialFlag = 0;
int contrastInitFlag = 0;
int brightnessInitFlag =0;
float *mm;
float **mmY;
int movieXFlag = 0;
int movieYFlag = 0;
int movieZFlag = 0;
float delta = 0;
int SFlag=0, IFlag=0, RFlag=0;
int LFlag=0, AFlag=0, PFlag=0;
int oldSFlag=0, oldIFlag=0, oldRFlag=0;
int oldLFlag=0, oldAFlag=0, oldPFlag=0;
int axialViewFlag=1, sagittalViewFlag=1, coronalViewFlag=1;
RGB_TRIPLE* rgb;
float ****fiberPointRGBA;
float slideraxial, slidercoronal, slidersagittal;
float oldslideraxial;
float oldslidercoronal;
float oldslidersagittal;
float brightnessPara = 0.50;
float contrastPara = 1.0;
float **floatFiberClr;
GlutKernelApplication::Register< tractomain > pva;
float oldxrot = 0, oldyrot = 0, oldzrot = 0;
int current_xflipflag = 0;
int current_yflipflag = 0;
int current_zflipflag = 0;
int ovly_flag = 0;
int alpha_flag = 0;

int drawcount = 0; // hack to get screen cap to work
int min_drawcount = 3;

/* parameters -- possibly set by user */
string fname_image;
string fname_ovly;
string fname_alpha;
vector<string> fname_fibers;
string fname_color;
float fiber_alpha = 0.4;
float min_fiber_alpha = 0.1;
float image_alpha = 0.7;
float ovly_offset = 0;
float ovly_thresh = 0;
float ovly_slope = 1;
float tract_xoffset = 0;
float tract_yoffset = 0;
float tract_zoffset = 0;
string initview = "S";
float xrot = 0;
float yrot = 0;
float zrot = 0;
int rendcorflag = 1;
int rendsagflag = 1;
int rendhorflag = 1;
int xflipflag = 0;
int yflipflag = 0;
int zflipflag = 1;
int savetifflag = 0;
string fname_tif = "./fibers.tif";

/*****************************************************************/
// from GlutKernelApplication.cxx

void GlutKernelApplication::OnAppReshape(int width, int height)
{
	_width = width;
	_height = height;
}


void GlutKernelApplication::OnMousePos(int xpos, int ypos)
{
    // Set the position, then update
    _mouse.setPosition( xpos, ypos );
    _mouse.updateEdgeStates();
    
    this->OnMouseEvent();
}

void GlutKernelApplication::OnMouseClick4(int button, int state, int xpos, int ypos)
{

}

void GlutKernelApplication::OnMouseClick3(int button, int state, int xpos, int ypos)
{
	cout << "click 3" << endl;

}
void GlutKernelApplication::OnMouseClick2(int button, int state, int xpos, int ypos)
{
	cout << "click 2" << endl;
}

void GlutKernelApplication::OnMouseClick(int button, int state, int xpos, int ypos)
{
  _keyboardModifier = glutGetModifiers();

  Mouse::Button b;
  Mouse::BinaryState binaryState;

  // set current wind 
  int winID = glutGetWindow();
  glutSetWindow(winID);

  switch(button) {
    case GLUT_LEFT_BUTTON: b = Mouse::LEFT; break;
    case GLUT_MIDDLE_BUTTON: b = Mouse::MIDDLE; break;
    case GLUT_RIGHT_BUTTON: b = Mouse::RIGHT; break;
    default: assert(false);
  }

  switch(state) {
    case GLUT_DOWN: binaryState = Mouse::ON; break;
    case GLUT_UP: binaryState = Mouse::OFF;  break;
    default: assert(false);
  }

  // Set the mousebutton state and the mouse position
  _mouse.setState(b, binaryState);
  _mouse.setPosition( xpos, ypos );
  _mouse.updateEdgeStates();

  this->OnMouseEvent();
}

void GlutKernelApplication::OnAppIdle()
{
   _keyboard.updateEdgeStates();
}


void GlutKernelApplication::OnKeyboardDown( unsigned char k, int x, int y )
{
   const Keyboard::BinaryState state = Keyboard::ON;
   cout << "key down!" << endl;
   getchar();
   if (k <= 127) // k >= 0 since k is unsigned
   {
      _keyboard.binaryState( k ) = state;
   }
   else
   {
      cout<<"unrecognized key = "<<k<<" = "<<(int)k<<"\n"<<flush;
      return;
   }
   switch(k)
   {
	case 'r':
		imageinfo.imageType = 1; 
		break;
	case 'g':
		imageinfo.imageType = 0; 
		break;
   }

   _keyboard.queue.enqueue( (Keyboard::Key)(int)k );
}

void GlutKernelApplication::OnKeyboardUp( unsigned char k, int x, int y )
{
   const Keyboard::BinaryState state = Keyboard::OFF;
   if (k <= 127) // k >= 0 since k is unsigned
   {
      _keyboard.binaryState( k ) = state;
   }
   else
   {
      cout<<"unrecognized key = "<<k<<" = "<<(int)k<<"\n"<<flush;
   }
}

void GlutKernelApplication::OnSpecialKeyboardDown( int k, int x, int y )
{
   const Keyboard::BinaryState state = Keyboard::ON;
   Keyboard::Key key;
   switch (k)
   {
      case GLUT_KEY_UP: _keyboard.binaryState( key = Keyboard::UPARROW ) = state; break;
      case GLUT_KEY_LEFT: _keyboard.binaryState( key = Keyboard::LEFTARROW ) = state; break;
      case GLUT_KEY_DOWN: _keyboard.binaryState( key = Keyboard::DOWNARROW ) = state; break;
      case GLUT_KEY_RIGHT: _keyboard.binaryState( key = Keyboard::RIGHTARROW ) = state; break;
      case GLUT_KEY_F1: _keyboard.binaryState( key = Keyboard::F1 ) = state; break;
      case GLUT_KEY_F2: _keyboard.binaryState( key = Keyboard::F2 ) = state; break;
      case GLUT_KEY_F3: _keyboard.binaryState( key = Keyboard::F3 ) = state; break;
      case GLUT_KEY_F4: _keyboard.binaryState( key = Keyboard::F4 ) = state; break;
      case GLUT_KEY_F5: _keyboard.binaryState( key = Keyboard::F5 ) = state; break;
      case GLUT_KEY_F6: _keyboard.binaryState( key = Keyboard::F6 ) = state; break;
      case GLUT_KEY_F7: _keyboard.binaryState( key = Keyboard::F7 ) = state; break;
      case GLUT_KEY_F8: _keyboard.binaryState( key = Keyboard::F8 ) = state; break;
      case GLUT_KEY_F9: _keyboard.binaryState( key = key = Keyboard::F9 ) = state; break;
      case GLUT_KEY_F10: _keyboard.binaryState( key = Keyboard::F10 ) = state; break;
      case GLUT_KEY_F11: _keyboard.binaryState( key = Keyboard::F11 ) = state; break;
      case GLUT_KEY_F12: _keyboard.binaryState( key = Keyboard::F12 ) = state; break;

      case GLUT_KEY_PAGE_UP: _keyboard.binaryState( key = Keyboard::PAGEUP ) = state; break;
      case GLUT_KEY_PAGE_DOWN: _keyboard.binaryState( key = Keyboard::PAGEDOWN ) = state; break;
      case GLUT_KEY_HOME: _keyboard.binaryState( key = Keyboard::HOME ) = state; break;
      case GLUT_KEY_END: _keyboard.binaryState( key = Keyboard::END ) = state; break;
      case GLUT_KEY_INSERT: _keyboard.binaryState( key = Keyboard::INSERT ) = state; break;
      default:
         cout<<"unrecognized key = "<<(int)k<<"\n"<<flush;
         return;
   }   
   //_keyboard.updateEdgeStates();
   _keyboard.queue.enqueue( key );
}

void GlutKernelApplication::OnSpecialKeyboardUp(int k, int x, int y)     
{
   const Keyboard::BinaryState state = Keyboard::OFF;
   Keyboard::Key key;
   switch (k)
   {
      case GLUT_KEY_UP: _keyboard.binaryState( key = Keyboard::UPARROW ) = state; break;
      case GLUT_KEY_LEFT: _keyboard.binaryState( key = Keyboard::LEFTARROW ) = state; break;
      case GLUT_KEY_DOWN: _keyboard.binaryState( key = Keyboard::DOWNARROW ) = state; break;
      case GLUT_KEY_RIGHT: _keyboard.binaryState( key = Keyboard::RIGHTARROW ) = state; break;
      case GLUT_KEY_F1: _keyboard.binaryState( key = Keyboard::F1 ) = state; break;
      case GLUT_KEY_F2: _keyboard.binaryState( key = Keyboard::F2 ) = state; break;
      case GLUT_KEY_F3: _keyboard.binaryState( key = Keyboard::F3 ) = state; break;
      case GLUT_KEY_F4: _keyboard.binaryState( key = Keyboard::F4 ) = state; break;
      case GLUT_KEY_F5: _keyboard.binaryState( key = Keyboard::F5 ) = state; break;
      case GLUT_KEY_F6: _keyboard.binaryState( key = Keyboard::F6 ) = state; break;
      case GLUT_KEY_F7: _keyboard.binaryState( key = Keyboard::F7 ) = state; break;
      case GLUT_KEY_F8: _keyboard.binaryState( key = Keyboard::F8 ) = state; break;
      case GLUT_KEY_F9: _keyboard.binaryState( key = key = Keyboard::F9 ) = state; break;
      case GLUT_KEY_F10: _keyboard.binaryState( key = Keyboard::F10 ) = state; break;
      case GLUT_KEY_F11: _keyboard.binaryState( key = Keyboard::F11 ) = state; break;
      case GLUT_KEY_F12: _keyboard.binaryState( key = Keyboard::F12 ) = state; break;

      case GLUT_KEY_PAGE_UP: _keyboard.binaryState( key = Keyboard::PAGEUP ) = state; break;
      case GLUT_KEY_PAGE_DOWN: _keyboard.binaryState( key = Keyboard::PAGEDOWN ) = state; break;
      case GLUT_KEY_HOME: _keyboard.binaryState( key = Keyboard::HOME ) = state; break;
      case GLUT_KEY_END: _keyboard.binaryState( key = Keyboard::END ) = state; break;
      case GLUT_KEY_INSERT: _keyboard.binaryState( key = Keyboard::INSERT ) = state; break;
      default:
         cout<<"unrecognized key = "<<(int)k<<"\n"<<flush;
         return;
   }
}

void GlutKernelApplication::AxialSliceInit(int argc, char *argv[])
{
  printf("%s: initializing axial slice\n",progname);
	imageinfo.imageType = 0;
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-imageinfo.width/2, imageinfo.width/2,  \
          -imageinfo.height/2, imageinfo.height/2,\
          -imageinfo.depth/2, imageinfo.depth/2);
  this->imageinfo.AxialLineH = 0.0;
  this->imageinfo.AxialLineL = 0.0;
  CoronalLineH = 0.0;
  CoronalLineL = 0.0;
  SagittalLineH = 0.0;
  SagittalLineL = 0.0;
	imageinfo.AsliceLoc = imageinfo.AxialSliceLocation;
	imageinfo.CsliceLoc = imageinfo.CoronalSliceLocation;
	imageinfo.SsliceLoc = imageinfo.SagittalSliceLocation;
}

void GlutKernelApplication::CoronalSliceInit(int argc, char *argv[])
{
  printf("%s: initializing coronal slice\n",progname);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-imageinfo.width/2, imageinfo.width/2,  \
          -imageinfo.height/2, imageinfo.height/2,\
          -imageinfo.depth/2, imageinfo.depth/2);
}
void GlutKernelApplication::SagittalSliceInit(int argc, char *argv[])
{
  printf("%s: initializing sagittal slice\n",progname);
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-imageinfo.width/2, imageinfo.width/2,  \
          -imageinfo.height/2, imageinfo.height/2,\
          -imageinfo.depth/2, imageinfo.depth/2);
}


void GlutKernelApplication::AxialDisplay(void)
{
	this->AxialDisplay();

}
void GlutKernelApplication::CoronalDisplay(void)
{
	this -> CoronalDisplay();
}

void GlutKernelApplication::SagittalDisplay(void)
{
	this -> SagittalDisplay();
}

void GlutKernelApplication::AxialMouseClick(int button, int state, int xpos, int ypos)
{                       
	this->AxialMouseClick(button, state, xpos, ypos);

}

void GlutKernelApplication::CoronalMouseClick(int button, int state, int xpos, int ypos)
{   
	this -> CoronalMouseClick(button, state, xpos, ypos);
}

void GlutKernelApplication::SagittalMouseClick(int button, int state, int xpos, int ypos)
{   
	this -> SagittalMouseClick(button, state, xpos, ypos);
}
void GlutKernelApplication::OperationPanelInit(int argc, char *argv[])
{
	char *CS[200];

	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE);
	muiInit();
	muiNewUIList(1);
	muiSetActiveUIList(1);

}
void GlutKernelApplication::OperationPanelOpenFilesInit(int argc, char *argv[])
{
	char *CS[200];

	//glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE);
	muiInit();
	muiNewUIList(2);
	muiSetActiveUIList(2);

}

void GlutKernelApplication::OperationPanelEvent( )
{

	this -> OperationPanelEvent();

}
void GlutKernelApplication::OperationPanelOpenFilesEvent( )
{

	this -> OperationPanelOpenFilesEvent();

}

/*****************************************************************/
// from glutMain.cxx

namespace GlutMain
{
	//ContextData<int> dispList;
	class DefaultFalseBool
	{
	private:
		bool _flag;
	public:
		DefaultFalseBool() : _flag(false) {}
		DefaultFalseBool(const bool& flag) : _flag(flag) {}
		inline bool& truth() { return _flag; }
		inline const bool& truth() const { return _flag; }
	};
	ContextData<GlutMain::DefaultFalseBool> oneTimeOnly;

	int currentContext = 0; //only one context for now...
	int mainWin_contextID = 0; // ID of the main (first) window
  static void OnRedisplay();

	static void OnIdle()
	{
    // According to the GLUT specification, the current window is 
    // undefined during an idle callback.  So we need to explicitly change
    // it if necessary
    if ( glutGetWindow() != mainWin_contextID ) {
       glutSetWindow( mainWin_contextID );  
    }
    // tell glut to call redisplay (which then calls OnRedisplay)
    glutPostRedisplay();
    int x;
    for (x = 0; x < GlutKernelApplication::applications().size(); ++x) {
       assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[x]->OnAppIdle();
    }
	}

	static void OnRedisplay() 
	{ 
		 // Initialize the context once and only once for every opengl window.
		 bool hasInitialized = GlutMain::oneTimeOnly( currentContext ).truth();
		 if (hasInitialized == false) {
			 GlutMain::oneTimeOnly( currentContext ).truth() = true;
       // init each application
       int x;
       for (x = 0; x < GlutKernelApplication::applications().size(); ++x) {
          assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
          GlutKernelApplication::applications()[x]->OnContextInit( );
       }
	   }
     // draw each application
     int x;
     for (x = 0; x < GlutKernelApplication::applications().size(); ++x) {
       assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[x]->OnContextDraw();
     }
     printf("%s: swapping buffers...\n",progname);
     glutSwapBuffers();
	}

	static void OnReshape(int width, int height) 
	{ 
		int x;
    for (x = 0; x < GlutKernelApplication::applications().size(); ++x) {
       assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[x]->OnAppReshape( width, height );
       GlutKernelApplication::applications()[x]->OnContextReshape(width, height);
    }
		OnIdle();
	}

	static void OnKeyboardDown(unsigned char k, int x, int y) 
  { 
    int a;
    for (a = 0; a < GlutKernelApplication::applications().size(); ++a) {
       assert( GlutKernelApplication::applications()[a] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[a]->OnKeyboardDown(k, x, y); 
    }
  }

	static void OnKeyboardUp(unsigned char k, int x, int y) 
  { 
    int a;
    for (a = 0; a < GlutKernelApplication::applications().size(); ++a) {
       assert( GlutKernelApplication::applications()[a] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[a]->OnKeyboardUp(k, x, y); 
    }
  }

	static void OnSpecialKeyboardDown(int k, int x, int y) 
  {
    int a;
    for (a = 0; a < GlutKernelApplication::applications().size(); ++a) {
       assert( GlutKernelApplication::applications()[a] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[a]->OnSpecialKeyboardDown(k, x, y); 
    }
  }

	static void OnSpecialKeyboardUp(int k, int x, int y) 
  { 
    int a;
    for (a = 0; a < GlutKernelApplication::applications().size(); ++a) {
       assert( GlutKernelApplication::applications()[a] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[a]->OnSpecialKeyboardUp(k, x, y); 
    }
  }

	static void OnMousePos( int x, int y )
  { 
     int a;
     for (a = 0; a < GlutKernelApplication::applications().size(); ++a) {
       assert( GlutKernelApplication::applications()[a] != NULL && "you registered a NULL application" );
       GlutKernelApplication::applications()[a]->OnMousePos( x, y ); 
     }
  }
	static void OnMouseClick( int a, int b, int c, int d ) 
  { 
     int x;
     for (x = 0; x < GlutKernelApplication::applications().size(); ++x) {
        assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
        GlutKernelApplication::applications()[x]->OnMouseClick(a, b, c, d); 
     }
  }
	static void AxialSliceInit( int argc, char *argv[]) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->AxialSliceInit(argc, argv); 
      }
   }
	static void AxialDisplay( void) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->AxialDisplay(); 
      }
   }
	static void CoronalSliceInit( int argc, char *argv[]) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->CoronalSliceInit(argc, argv); 
      }
   }
	static void CoronalDisplay( void) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->CoronalDisplay(); 
      }
   }
	static void SagittalSliceInit( int argc, char *argv[]) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->SagittalSliceInit(argc, argv); 
      }
   }
	static void SagittalDisplay( void) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->SagittalDisplay(); 
      }
   }
	static void AxialMouseClick( int a, int b, int c, int d ) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->AxialMouseClick(a, b, c, d); 
      }
   }
	static void CoronalMouseClick( int a, int b, int c, int d ) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->CoronalMouseClick(a, b, c, d); 
      }
   }
	static void SagittalMouseClick( int a, int b, int c, int d ) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->SagittalMouseClick(a, b, c, d); 
      }
   }
	static void OperationPanelInit( int argc, char *argv[]) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->OperationPanelInit(argc, argv); 
      }
   }
	static void OperationPanelEvent( ) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->OperationPanelEvent(); 
      }
   }
	static void OperationPanelOpenFilesInit( int argc, char *argv[]) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->OperationPanelOpenFilesInit(argc, argv); 
      }
   }
	static void OperationPanelOpenFilesEvent( ) 
   { 
      int x;
      for (x = 0; x < GlutKernelApplication::applications().size(); ++x)
      {
         assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
         GlutKernelApplication::applications()[x]->OperationPanelOpenFilesEvent(); 
      }
   }
};



int main( int argc, char *argv[] )
{	
  assert( GlutKernelApplication::applications().size() > 0 && "you must register at least one application" );

  //Initialize all registered applications, do this before initing glut, in case app
  // needs to set window position and name.
  int x;

  printf("glutmain: initializing application...\n");
  
  for (x = 0; x < GlutKernelApplication::applications().size(); ++x) {
    assert( GlutKernelApplication::applications()[x] != NULL && "you registered a NULL application" );
    GlutKernelApplication::applications()[x]->OnAppInit( argc, argv );
  }

  printf("glutmain: setting window position...\n");

  // Set the window position
  assert( GlutKernelApplication::applications()[0]->width() > 0 && 
  GlutKernelApplication::applications()[0]->height() > 0 &&
  "Width and height of application must be > 0" );
  ::glutInitWindowSize( GlutKernelApplication::applications()[0]->width(), GlutKernelApplication::applications()[0]->height() );
  ::glutInit( &argc, argv );
  ::glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );


  // Set the window title
  if ( GlutKernelApplication::applications()[0]->name() == "" ) {
    GlutKernelApplication::applications()[0]->setName( "OpenGL" );
  }

  ::glutInitWindowPosition( GlutKernelApplication::applications()[0]->width()*2/3+20, 0 );

  printf("glutmain: creating window...\n");

  GlutMain::mainWin_contextID = glutCreateWindow( GlutKernelApplication::applications()[0]->name() );
  glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );

  printf("glutmain: setting up callback functions...\n");

  // display callbacks.
  glutDisplayFunc( GlutMain::OnRedisplay );
  glutReshapeFunc( GlutMain::OnReshape );
  glutIdleFunc( GlutMain::OnIdle );

  // keyboard callback functions.
  glutKeyboardFunc( GlutMain::OnKeyboardDown );
  glutKeyboardUpFunc( GlutMain::OnKeyboardUp );
  glutSpecialFunc( GlutMain::OnSpecialKeyboardDown );
  glutSpecialUpFunc( GlutMain::OnSpecialKeyboardUp );

  // mouse callback functions...
  glutMouseFunc( GlutMain::OnMouseClick );
  glutMotionFunc( GlutMain::OnMousePos );
  glutPassiveMotionFunc( GlutMain::OnMousePos );
  glutPopWindow();

  printf("glutmain: creating additional windows...\n");

  ::glutInitWindowSize( GlutKernelApplication::applications()[0]->width()/3, GlutKernelApplication::applications()[0]->height()/3 );
  ::glutInitWindowPosition( GlutKernelApplication::applications()[0]->width()/3+10, 0 );
  int tempID1 = glutCreateWindow( "Axial" );
  GlutMain::AxialSliceInit(argc, argv);
  glutDisplayFunc(GlutMain::AxialDisplay);
  glutMouseFunc( GlutMain::AxialMouseClick );

  ::glutInitWindowSize( GlutKernelApplication::applications()[0]->width()/3, GlutKernelApplication::applications()[0]->height()/3 );
  ::glutInitWindowPosition( GlutKernelApplication::applications()[0]->width()/3+10,  GlutKernelApplication::applications()[0]->height()/3+25 );
  int tempID2 = glutCreateWindow( "Coronal" );
  GlutMain::CoronalSliceInit(argc, argv);
  glutDisplayFunc(GlutMain::CoronalDisplay);
  glutMouseFunc( GlutMain::CoronalMouseClick );

  ::glutInitWindowSize( GlutKernelApplication::applications()[0]->width()/3, GlutKernelApplication::applications()[0]->height()/3 );
  ::glutInitWindowPosition( GlutKernelApplication::applications()[0]->width()/3+10, 2*GlutKernelApplication::applications()[0]->height()/3 +50);
  int tempID3 = glutCreateWindow( "Sagittal" );
  GlutMain::SagittalSliceInit(argc, argv);
  glutDisplayFunc(GlutMain::SagittalDisplay);
  glutMouseFunc( GlutMain::SagittalMouseClick );

  ::glutInitWindowSize( GlutKernelApplication::applications()[0]->width()/3, GlutKernelApplication::applications()[0]->height()+50 );
  ::glutInitWindowPosition( 0, 0);
  int tempID4 = glutCreateWindow( "Operation Panel" );
  GlutMain::OperationPanelInit(argc, argv);
  GlutMain::OperationPanelEvent();

  printf("glutmain: starting main loop...\n");

  ::glutMainLoop();

  return (1);
}

/*****************************************************************/

void tractomain::OnAppInit( int argc, char *argv[] )
{
  parseopts(argc,argv);

  this->setName( "Neuro Fiber Tracking 3D Visualization" );
  this->setWidth( DEFAULT_WIN_WIDTH );
  this->setHeight( DEFAULT_WIN_HEIGHT );

  mmY = (float **)calloc(100, sizeof(float *));

  for(int i = 0; i<STEPNUM; i++){
    rotationStep[i] =  1.0/180.0 * PI;
  }

  // set up a light
  light.setNumber( 0 );
  light.setColor( Light::ambient, 0.5f, 0.5f, 0.5f );    
  light.setColor( Light::diffuse, 0.7f, 0.7f, 0.7f );    
  light.setColor( Light::specular, 0.9f, 0.9f, 0.9f );    
  light.setPos( 0, 0, 1, 0 );
  light.setAtten( 1, 0, 0 );

  // set up a second light
  light2.setNumber( 1 );
  light2.setColor( Light::ambient, 0.3f, 0.3f, 0.5f );    
  light2.setColor( Light::diffuse, 0.6f, 0.6f, 0.8f );    
  light2.setColor( Light::specular, 0.9f, 0.9f, 1.0f );    
  light2.setPos( 0, 0, -1, 0.0 );
  light2.setAtten( 1, 0, 0 );

  beacon.timeLeft = 180;

  this->imageinfo.fiberNum = 0;
  this->imageinfo = readfiles(argc, argv);

  WW = DEFAULT_WIDTH;
  WL = DEFAULT_LEVEL;
}

int FileExists(const char *fname)
{
  FILE *fp ;

  fp = fopen(fname, "r") ;
  if (fp) fclose(fp) ;
  return(fp != NULL) ;
}

void usage()
{
  printf("\nusage: %s [opts]\n\n", progname);
  printf(" Options:\n");
  printf("   -fname_image <string>        : background image filename\n");
  printf("   -fname_fibers <string list>  : list of fiber files\n");
  printf("   -fname_ovly <string>         : overlay image filename\n");
  printf("       displayed as colorscale\n");
  printf("   -ovly_offset <float>         : overlay offset value\n");
  printf("       overlay values at offset are purple\n");
  printf("       values above are red, below are blue\n");
  printf("                                  (default = 0)\n");
  printf("   -ovly_slope <float>          : overlay slope\n");
  printf("                                  (default = 1)\n");
  printf("   -ovly_thresh <float>          : overlay threshold\n");
  printf("                                  (default = 0)\n");
  printf("   -tract_xoffset <float>        : x shift of fibers relative to image\n");
  printf("                                  (default = 0)\n");
  printf("   -tract_yoffset <float>        : y shift of fibers relative to image\n");
  printf("                                  (default = 0)\n");
  printf("   -tract_zoffset <float>        : z shift of fibers relative to image\n");
  printf("                                  (default = 0)\n");
  printf("   -fname_alpha <string>        : fiber alpha filename\n");
  printf("       controls fiber transparency\n");
  printf("       values in this file will be normalized to max value\n");
  printf("   -fiber_alpha <float>         : fiber transparency max value (0-1)\n");
  printf("                                  (default = 0.4)\n");
  printf("   -min_fiber_alpha <float>     : minimum fiber transparency value (0-1)\n");
  printf("       ignored unless fname_alpha is supplied\n");
  printf("       must be less than fiber_alpha\n");
  printf("                                  (default = 0.1)\n");
  printf("   -fname_color <string>        : color pallet filename\n");
  printf("                                  (ignored if fname_overlay supplied)\n");
  printf("   -image_alpha <float>         : image transparency value (0-1)\n");
  printf("                                  (default = 0.4)\n");
  printf("   -initview <string>           : initial view (S, I, A, P, L, or R)\n");
  printf("                                  (default = \"S\")\n");
  printf("   -xrot <float>                : initial rotation about x axis (degrees)\n");
  printf("   -yrot <float>                : initial rotation about y axis (degrees)\n");
  printf("   -zrot <float>                : initial rotation about z axis (degrees)\n");
  printf("   -[no]rendcor                 : render coronal slice\n");
  printf("   -[no]rendsag                 : render sagittal slice\n");
  printf("   -[no]rendhor                 : render horizontal slice\n");
  printf("   -[no]xflip                   : flip x axis for fibers\n");
  printf("                                  (default = noxflip)\n");
  printf("   -[no]yflip                   : flip y axis for fibers\n");
  printf("                                  (default = noyflip)\n");
  printf("   -[no]zflip                   : flip z axis for fibers\n");
  printf("                                  (default = zflip)\n");
  printf("   -savetif                     : save image as tif and quit\n");
  printf("   -fname_tif <string>          : output tif filename\n");
  printf("                                  (default = ./fibers.tif)\n");
  printf("\n");
}

inline int checkopt(string *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt = argv[i];
  return(0);
}

inline int checkopt(char *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  strcpy(opt,argv[i]);
  return(0);
}

inline int checkopt(int *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt=atoi(argv[i]);
  return(0);
}

inline long checkopt(long *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt=atol(argv[i]);
  return(0);
}

inline int checkopt(float *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt=atof(argv[i]);
  return(0);
}

inline void nooptexit(const char *optname)
{
  printf("%s: ### missing parameter value \"%s\"\n",progname, optname);
  exit(1);
}

void parseopts(int argc, char *argv[])
{
  int i, k, j;
  char temp[STRLEN];
  char *tmpstr;
  char *pch;
  string tmp;

  progname = argv[0];
  /* strip off path */
  pch = strtok(progname,"/");
  while (pch!=NULL) {
    strcpy(temp,pch);
    pch = strtok (NULL,"/");
  }
  strcpy(progname,temp);

  if (argc==1) {usage(); exit(0);}

  printf("%s: parsing options...\n",progname);

  // options override defaults
  for (i=1;i<argc;i++) {
    strcpy(temp,argv[i]);
    if (MATCH(temp,"-h")) {usage(); exit(0);}
    else if ((MATCH(temp,"-fnameimage")) ||
            (MATCH(temp,"-fname_image")))
            {if(checkopt(&fname_image,i+1,argc,argv)) nooptexit("fname_image"); i++;}
    else if (MATCH(temp,"-fname_ovly"))
            {if(checkopt(&fname_ovly,i+1,argc,argv)) nooptexit("fname_ovly"); i++;}
    else if (MATCH(temp,"-ovly_offset"))       {if(checkopt(&ovly_offset,i+1,argc,argv)) nooptexit("ovly_offset"); i++;}
    else if (MATCH(temp,"-ovly_slope"))        {if(checkopt(&ovly_slope,i+1,argc,argv)) nooptexit("ovly_slope"); i++;}
    else if (MATCH(temp,"-ovly_thresh"))       {if(checkopt(&ovly_thresh,i+1,argc,argv)) nooptexit("ovly_thresh"); i++;}
    else if (MATCH(temp,"-tract_xoffset"))     {if(checkopt(&tract_xoffset,i+1,argc,argv)) nooptexit("tract_xoffset"); i++;}
    else if (MATCH(temp,"-tract_yoffset"))     {if(checkopt(&tract_yoffset,i+1,argc,argv)) nooptexit("tract_yoffset"); i++;}
    else if (MATCH(temp,"-tract_zoffset"))     {if(checkopt(&tract_zoffset,i+1,argc,argv)) nooptexit("tract_zoffset"); i++;}
    else if (MATCH(temp,"-fname_alpha"))
            {if(checkopt(&fname_alpha,i+1,argc,argv)) nooptexit("fname_alpha"); i++;}
    else if ((MATCH(temp,"-fnamecolor")) ||
            (MATCH(temp,"-fname_color")))
            {if(checkopt(&fname_color,i+1,argc,argv)) nooptexit("fname_color"); i++;}
    else if ((MATCH(temp,"-fnamefibers")) ||
            (MATCH(temp,"-fname_fibers"))) {
      for(j=0;j<MAX_FIBER_FILES;j++) {
        if(checkopt(&tmp,i+1,argc,argv)) {
          if(nfibers==0) nooptexit("fname_fibers"); else break;
        }
        if (tmp.c_str()[0]=='-') break;
        fname_fibers.push_back(tmp);
        nfibers++;
        i++;
      }
    }
    else if (MATCH(temp,"-fiber_alpha"))       {if(checkopt(&fiber_alpha,i+1,argc,argv)) nooptexit("fiber_alpha"); i++;}
    else if (MATCH(temp,"-min_fiber_alpha"))   {if(checkopt(&min_fiber_alpha,i+1,argc,argv)) nooptexit("min_fiber_alpha"); i++;}
    else if (MATCH(temp,"-image_alpha"))       {if(checkopt(&image_alpha,i+1,argc,argv)) nooptexit("image_alpha"); i++;}
    else if (MATCH(temp,"-initview"))          {if(checkopt(&initview,i+1,argc,argv)) nooptexit("initview"); i++;}
    else if (MATCH(temp,"-xrot"))              {if(checkopt(&xrot,i+1,argc,argv)) nooptexit("xrot"); i++;}
    else if (MATCH(temp,"-yrot"))              {if(checkopt(&yrot,i+1,argc,argv)) nooptexit("yrot"); i++;}
    else if (MATCH(temp,"-zrot"))              {if(checkopt(&zrot,i+1,argc,argv)) nooptexit("zrot"); i++;}
    else if (MATCH(temp,"-rendcor"))           {rendcorflag = 1;}
    else if (MATCH(temp,"-norendcor"))         {rendcorflag = 0;}
    else if (MATCH(temp,"-rendsag"))           {rendsagflag = 1;}
    else if (MATCH(temp,"-norendsag"))         {rendsagflag = 0;}
    else if (MATCH(temp,"-rendhor"))           {rendhorflag = 1;}
    else if (MATCH(temp,"-norendhor"))         {rendhorflag = 0;}
    else if (MATCH(temp,"-xflip"))             {xflipflag = 1;}
    else if (MATCH(temp,"-noxflip"))           {xflipflag = 0;}
    else if (MATCH(temp,"-yflip"))             {yflipflag = 1;}
    else if (MATCH(temp,"-noyflip"))           {yflipflag = 0;}
    else if (MATCH(temp,"-zflip"))             {zflipflag = 1;}
    else if (MATCH(temp,"-nozflip"))           {zflipflag = 0;}
    else if (MATCH(temp,"-savetif"))           {savetifflag = 1;}
    else if (MATCH(temp,"-nosavetif"))         {savetifflag = 0;}
    else if (MATCH(temp,"-fname_tif"))         {if(checkopt(&fname_tif,i+1,argc,argv)) nooptexit("fname_tif"); i++;}
    else {
      usage(); 
      printf("\n%s: ERROR: unrecognized option \"%s\"\n", progname, argv[i]);
      exit(1);
    }
  }

  // check options
  printf("%s: checking options...\n",progname);

  if (fname_image.size()>0) {
    if(!FileExists(fname_image.c_str())) {
      printf("%s: ERROR: fname_image %s not found\n",progname,fname_image.c_str());
      exit(1);    
    }    
  } else {
    printf("%s: ERROR: fname_image unspecified\n",progname);
    exit(1);    
  }

  if (fname_ovly.size()>0) {
    if(!FileExists(fname_ovly.c_str())) {
      printf("%s: ERROR: fname_ovly %s not found\n",progname,fname_ovly.c_str());
      exit(1);    
    }    
    ovly_flag = 1;
  } else if (fname_color.size()>0) {
    if(!FileExists(fname_color.c_str())) {
      printf("%s: ERROR: fname_color %s not found\n",progname,fname_color.c_str());
      exit(1);    
    }    
  }

  if (fiber_alpha<0 || fiber_alpha>1) {
    printf("%s: ERROR: invalid fiber_alpha (%0.2f) (should be positive, less than 1)\n",
      progname,fiber_alpha);
    exit(1);
  }
  if (image_alpha<0 || image_alpha>1) {
    printf("%s: ERROR: invalid image_alpha (%0.2f) (should be positive, less than 1)\n",
      progname,image_alpha);
    exit(1);
  }
  if (fname_alpha.size()>0) {
    if(!FileExists(fname_alpha.c_str())) {
      printf("%s: ERROR: fname_alpha %s not found\n",progname,fname_alpha.c_str());
      exit(1);
    }
    alpha_flag = 1;
    if (min_fiber_alpha<0 || min_fiber_alpha>1) {
      printf("%s: ERROR: invalid min_fiber_alpha (%0.2f) (should be positive, less than 1)\n",
        progname,min_fiber_alpha);
      exit(1);
    }
    if (min_fiber_alpha>fiber_alpha) {
      printf("%s: ERROR: min_fiber_alpha (%0.2f) must be less than fiber_alpha (%0.2f)\n",
        progname,min_fiber_alpha,fiber_alpha);
      exit(1);
    }
  }
  
  if(nfibers==0) {
    printf("%s: ERROR: fiber files unspecified\n",progname);
    exit(1);    
  }
  for(j=0;j<nfibers;j++) {
    if(!FileExists(fname_fibers[j].c_str())) {
      printf("%s: ERROR: fname_fiber %s not found\n",progname,fname_fibers[j].c_str());
      exit(1);    
    }
  }
  if(rendhorflag) {
    axialViewFlag = 1;
  } else {
    axialViewFlag = 0;
  }
  if(rendsagflag) {
    sagittalViewFlag = 1;
  } else {
    sagittalViewFlag = 0;
  }
  if(rendcorflag) {
    coronalViewFlag = 1;
  } else {
    coronalViewFlag = 0;
  }

  printf("%s: finished checking options\n",progname);
}

imageStruct tractomain::readfiles(int argc, char *argv[])
{
  FILE *fp;
  FIBER* fiber_array;
  mgh_header mgh_image, mgh_ovly, mgh_alpha;
  float *ImageArray;
  imageStruct imageInfo;
  int i, k, j, l;
  int fiberFilePosition;
  char temp[STRLEN];
  int index;
  float r, g, b, a, val, sval;
  int numfibers;
  float maxval;
  int shift_xc, shift_yc, shift_zc;

  // read image file
  printf("%s: reading image file %s...\n",progname,fname_image.c_str());
  ReadMghFile(fname_image.c_str(),&mgh_image,1);
  rows =  mgh_image.width;
  columns = mgh_image.height;
  slices = mgh_image.depth;
  xPixelSpacing = mgh_image.xyz_size[0] ;
  zSliceThickness =  mgh_image.xyz_size[2];
  ns = rows * columns;
  np = ns * slices;

  // allocate memory for image, rescale if necessary
  printf("%s: allocating memory for image data...\n",progname);
  ImageArray = (float *) calloc ( np, sizeof(float));
  printf("%s: copying image values...\n",progname);
  for(i = 0; i<np; i++) ImageArray[i] = mgh_image.floatImageArray[i];
  printf("%s: checking background image scaling...\n",progname);
  maxval = -BIGF;
  for(i = 0; i<np; i++) {
    if(ImageArray[i]>maxval) maxval=ImageArray[i];
  }
//  if(maxval>256) maxval=256;
  printf("%s: image max value = %0.2f...\n",progname,maxval);
  if(maxval>2) {
    printf("%s: rescaling image...\n",progname);
    for(i = 0; i<np; i++) {
      ImageArray[i] = ImageArray[i] / maxval;
    }
  }  

  if(ovly_flag) {
    // read overlay file
    printf("%s: reading overlay file %s...\n",progname,fname_ovly.c_str());
    ReadMghFile(fname_ovly.c_str(),&mgh_ovly,1);
    if ((rows !=  mgh_ovly.width) ||
        (columns != mgh_ovly.height) ||
        (slices != mgh_ovly.depth) ||
        (xPixelSpacing != mgh_ovly.xyz_size[0]) ||
        (zSliceThickness !=  mgh_ovly.xyz_size[2])) {
      printf("%s: ERROR: mismatch between image and overlay files\n",progname);
      exit(1);
    }
  }

  if(alpha_flag) {
    // read alpha file
    printf("%s: reading alpha file %s...\n",progname,fname_alpha.c_str());
    ReadMghFile(fname_alpha.c_str(),&mgh_alpha,1);
    if ((rows !=  mgh_alpha.width) ||
        (columns != mgh_alpha.height) ||
        (slices != mgh_alpha.depth) ||
        (xPixelSpacing != mgh_alpha.xyz_size[0]) ||
        (zSliceThickness !=  mgh_alpha.xyz_size[2])) {
      printf("%s: ERROR: mismatch between image and alpha files\n",progname);
      exit(1);
    }
    printf("%s: checking alpha image scaling...\n",progname,ns,np);
    maxval = -BIGF;
    for(i = 0; i<np; i++) {
      if(fabs(mgh_alpha.floatImageArray[i])>maxval)
        maxval=fabs(mgh_alpha.floatImageArray[i]);
    }
//    if(maxval>256) maxval=256;
    printf("%s: image max value = %0.2f...\n",progname,maxval);
    printf("%s: rescaling image...\n",progname);
    for(i = 0; i<np; i++) {
      mgh_alpha.floatImageArray[i] = fabs(mgh_alpha.floatImageArray[i]) / maxval;
    }
  }

  // allocate space for fibers
  printf("%s: allocating space for fibers...\n",progname);
  imageInfo.pasFiberAry = (FIBER **)calloc(nfibers, sizeof(FIBER *)); 
  imageInfo.totalFiberNum = (int *) calloc(nfibers, sizeof(int));

  // read fiber files
  printf("%s: reading fiber files...\n",progname);
  this->imageinfo.fiberNum = nfibers;
  for(j=0;j<nfibers;j++) {
    printf("%s: reading fiber file %s...\n",progname,fname_fibers[j].c_str());
    fiber_array = ReadFiber(fname_fibers[j],&numfibers);
    imageInfo.pasFiberAry[j] = fiber_array;
    imageInfo.totalFiberNum[j] = numfibers;
  }
  printf("%s: finished reading fiber files\n",progname);


  // add tract offset
  if (tract_xoffset!=0 || tract_yoffset!=0 || tract_zoffset!=0) {
    printf("%s: adding tract offsets...\n",progname);
    for(i=0; i<nfibers; i++) {
      for (j=0; j<imageInfo.totalFiberNum[i]; j++) {
        for (k=0; k<imageInfo.pasFiberAry[i][j].nLength; k++) {
          imageInfo.pasFiberAry[i][j].pxyzChain[k].x += tract_xoffset;
          imageInfo.pasFiberAry[i][j].pxyzChain[k].y += tract_yoffset;
          imageInfo.pasFiberAry[i][j].pxyzChain[k].z += tract_zoffset;
        }
      }
    }
    printf("%s: finished adding tract offsets\n",progname);
  }

  if (!ovly_flag) {
    // initialize colors for each fiber
    imageInfo.FiberColor = (RGB_TRIPLE *)calloc(this->imageinfo.fiberNum, sizeof(RGB_TRIPLE));
    rgb = (RGB_TRIPLE *)calloc(this->imageinfo.fiberNum, sizeof(RGB_TRIPLE));
    floatFiberClr = (float **)calloc(this->imageinfo.fiberNum, sizeof(float *));
    for (k=0; k<this->imageinfo.fiberNum; k++) {
      floatFiberClr[k] = (float *)calloc(3, sizeof(float));
    }
    if (fname_color.size()>0) {
      printf("%s: loading color file %s...\n",progname,fname_color.c_str());
      if ( (fp = fopen(fname_color.c_str(), "r")) == NULL) {
        printf("%s: failed to open color file %s\n", fname_color.c_str());
        exit(1);
      }
      unsigned char rgbtmp[this->imageinfo.fiberNum][3];
      char tempchar[STRLEN];
      for(k=0; k<this->imageinfo.fiberNum; k++) {
        char rgbtmp[8];
        int tempint[3];
        fgets(tempchar,STRLEN,fp);
        // todo: why is tempchar2 necessary?
        char *tempchar2 = strdup(tempchar);
        j = 0;
        l = 0;
        for(i=0; i<strlen(tempchar2); i++) {
          if(l>=3) {
            break;
          } else if(tempchar2[i] == ' ' || tempchar2[i] == '\t') {
            tempint[l] = atoi(rgbtmp);
            for(int kk = 0; kk < 8; kk++) {
              rgbtmp[kk] = ' ';
            }
            j=0;
            l++;
            continue;
          } else if(tempchar2[i] == '\0' || tempchar2[i] == '\n') {
            tempint[l] = atoi(rgbtmp);    
            for(int kk = 0; kk < 8; kk++) {
              rgbtmp[kk] = ' ';
            }
            j=0;
            l++;
            break;
          } else {
            rgbtmp[j] = tempchar2[i];
            j++;
          }
        }
        rgb[k].r = (unsigned char)tempint[0];
        rgb[k].g = (unsigned char)tempint[1];
        rgb[k].b = (unsigned char)tempint[2];
      }
      imageInfo.FiberColor = rgb; 
      fclose(fp);
    } else {
      printf("%s: generating random colors...\n",progname);
      for(i=0; i<this->imageinfo.fiberNum; i++) {
        rgb[i].r = uchar((255.0*rand() / RAND_MAX));
        rgb[i].g = uchar((255.0*rand() / RAND_MAX));
        rgb[i].b = uchar((255.0*rand() / RAND_MAX));
      }
      imageInfo.FiberColor = rgb; 
    }
    for(i =0; i<this->imageinfo.fiberNum; i++) {
      floatFiberClr[i][0] = (float ) imageInfo.FiberColor[i].r/255.0;
      floatFiberClr[i][1] = (float ) imageInfo.FiberColor[i].g/255.0;
      floatFiberClr[i][2] = (float ) imageInfo.FiberColor[i].b/255.0;
    }
  }
  
  printf("%s: setting imageInfo...\n",progname);
  imageInfo.width = columns;
  imageInfo.height = rows;
  imageInfo.depth = slices ;
  imageInfo.voxel_w = xPixelSpacing;
  imageInfo.voxel_h = xPixelSpacing;
  imageInfo.voxel_t = zSliceThickness;
  imageInfo.ImageArray = ImageArray;
  imageInfo.AxialSliceLocation = mgh_image.depth/2;
  imageInfo.CoronalSliceLocation = mgh_image.width/2;
  imageInfo.SagittalSliceLocation = mgh_image.height/2;
  imageInfo.AsliceLoc = imageInfo.AxialSliceLocation;
  imageInfo.CsliceLoc = imageInfo.CoronalSliceLocation;
  imageInfo.SsliceLoc = imageInfo.SagittalSliceLocation;
  imageInfo.AxialLineL = 0.0;
  imageInfo.AxialLineH = 0.0;
  imageInfo.CoronalLineL = 0.0;
  imageInfo.CoronalLineH = 0.0;
  imageInfo.SagittalLineL = 0.0;
  imageInfo.SagittalLineH = 0.0;
  imageInfo.fiberNum =  this->imageinfo.fiberNum;

  shift_xc = (int)((float)imageinfo.width/2+(float)0.0);
  shift_yc = (int)((float)imageinfo.height/2+(float)0.0);
  shift_zc = (int)((float)imageinfo.depth/2+(float)0.0);

  printf("%s: summarizing fiber info...\n",progname);
  printf("%s: number of fiber bundles = %d\n",progname,imageInfo.fiberNum);
  for(i=0; i<imageInfo.fiberNum; i++) {
    printf("%s: total fibers for fiber %d = %d\n",progname,i,imageInfo.totalFiberNum[i]);
  }
  int max_fiberNum = -100;
  for(i=0; i<imageInfo.fiberNum; i++) {
    if(imageInfo.totalFiberNum[i] > max_fiberNum) {
      max_fiberNum = imageInfo.totalFiberNum[i];
    }
  }
  int max_fiberLen= -100;
  for(i=0; i<imageInfo.fiberNum; i++) {
    for(j=0; j<imageInfo.totalFiberNum[i]; j++) {
      if(imageInfo.pasFiberAry[i][j].nLength  > max_fiberLen) {
        max_fiberLen = imageInfo.pasFiberAry[i][j].nLength;
      }
    }
  }
  printf("%s: max number of fiber bundles = %d\n",progname,max_fiberNum);
  printf("%s: max fiber length = %d\n",progname,max_fiberLen);

  // set fiber colors
  printf("%s: setting fiber colors...\n",progname);
  if ( (fiberPointRGBA = (float ****)calloc(imageInfo.fiberNum, sizeof(float ***))) == NULL) {
    cout << "failed to allocate memory for fiberPointRGBA****"<<endl;
    exit(1);
  }
  for(i=0; i<imageInfo.fiberNum; i++) {
    if( (fiberPointRGBA[i] = (float ***)calloc(max_fiberNum, sizeof(float **))) == NULL) {
      cout << "failed to allocate memory for fiberPointRGBA***"<<endl;
      exit(1);
    }
  }
  for(i=0; i<imageInfo.fiberNum; i++) {
    for (k=0; k<max_fiberNum; k++) {
      if( (fiberPointRGBA[i][k] = (float **)calloc(max_fiberLen, sizeof(float *))) == NULL) {
        cout << "failed to allocate memory for fiberPointRGBA**"<<endl;
        exit(1);
      }
    }
  }
  for(i=0; i<imageInfo.fiberNum; i++) {
    for (k=0; k<max_fiberNum; k++) {
      for (j=0; j<max_fiberLen; j++) {
        if( (fiberPointRGBA[i][k][j] = (float *)calloc(4, sizeof(float))) == NULL) {
          cout << "failed to allocate memory for fiberPointRGBA*"<<endl;
          exit(1);
        }
      }
    }
  }
  if(ovly_flag) {
    // set colors for each fiber point based on overlay data
    for(i=0; i<imageInfo.fiberNum; i++) {
      for (j=0; j<imageInfo.totalFiberNum[i]; j++) {
        for (k=0; k<imageInfo.pasFiberAry[i][j].nLength; k++) {
          index = (int)(imageInfo.pasFiberAry[i][j].pxyzChain[k].z+shift_zc) * ns +
                   (int)(imageInfo.pasFiberAry[i][j].pxyzChain[k].y+shift_yc) * rows +
                    (int)(imageInfo.pasFiberAry[i][j].pxyzChain[k].x+shift_xc);
          val = mgh_ovly.floatImageArray[index];
          if (fabs(val)<ovly_thresh)
            val = 0;
          else
            val = val-ovly_offset;
          
          sval = val * ovly_slope;

          // todo: add more colormaps
          
          // bpr: blue = -, purple = 0, red = +
          r = (1+sval)/2;
          g = 0;
          b = (1-sval)/2;

          if (r>1) r=1; if (r<0) r=0;
          if (b>1) b=1; if (b<0) b=0;

          if(alpha_flag) {
            a =  min_fiber_alpha + 
                 (fiber_alpha - min_fiber_alpha) * 
                 mgh_alpha.floatImageArray[index];
          } else {
            a =  fiber_alpha;
          }
          if (a>1) a=1;
          if (fabs(val)<ovly_thresh) a=min_fiber_alpha;
          if (a<min_fiber_alpha) a=min_fiber_alpha;

          fiberPointRGBA[i][j][k][0] =  r;
          fiberPointRGBA[i][j][k][1] =  g;
          fiberPointRGBA[i][j][k][2] =  b;
          fiberPointRGBA[i][j][k][3] =  a;
        }
      }
    }
  } else {
    // set colors for each fiber point (same color for entire fiber)
    for(i=0; i<imageInfo.fiberNum; i++) {
      for (j=0; j<imageInfo.totalFiberNum[i]; j++) {
        for (k=0; k<imageInfo.pasFiberAry[i][j].nLength; k++) {
          fiberPointRGBA[i][j][k][0] =  floatFiberClr[i][0];
          fiberPointRGBA[i][j][k][1] =  floatFiberClr[i][1];
          fiberPointRGBA[i][j][k][2] =  floatFiberClr[i][2];
          fiberPointRGBA[i][j][k][3] =  fiber_alpha;
        }
      }
    }
  }

  printf("%s: finished setting fiber colors\n",progname);

  return imageInfo;
}

int tractomain::ReadMghFile(const char *mghFileName,mgh_header *mghInfor,bool endianInfor)
{
  int i, j, k;
  int tmp_int[7];
  short tmp_short;
  int n=1;
  int endianFlag=-1;      // big endian: 1; little endian: 0; initial: -1
  FILE *fp;
  int unused_space_size = UNUSED_SPACE_SIZE-2;
  int nv;
  int fileLen;
  int headerSize;

  // check the endian type of local machine       
  if( *(char*)&n == 1 ) {
    endianFlag=0;
    //cout << "\t#####Little endian machine." << endl; 
  } else {
    endianFlag=1;
    //cout << "Big endian machine." << endl;
  }

  if ( (fp = fopen(mghFileName, "r")) == NULL) {
    printf("Could not open %s, exiting ...\n", mghFileName);
    return (-1);
  }

  fseek(fp, 0, SEEK_END);
  fileLen = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  //printf("%s: length of file is %d\n",progname,fileLen);

  fread( tmp_int, sizeof(int), 7, fp);
  if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0)) {
    for ( i=0; i<=6; i++) {
      if (endianInfor == 1) {
        bswap(tmp_int[i]);
      }
    }
  }
  fread(&tmp_short , sizeof(short), 1, fp);
  if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0)) {
    bswap(tmp_short);
  }
  mghInfor->v = tmp_int[0];
  mghInfor->width = tmp_int[1];
  mghInfor->height = tmp_int[2];
  mghInfor->depth = tmp_int[3];
  mghInfor->nframes = tmp_int[4];
  mghInfor->type = tmp_int[5];
  mghInfor->dof = tmp_int[6];
  mghInfor->ras_good_flag = tmp_short;
  if(mghInfor->ras_good_flag) {
    unused_space_size = unused_space_size - USED_SPACE_SIZE;
    fread( mghInfor->xyz_size, sizeof(float), 3, fp);
    fread( mghInfor->x_ras, sizeof(float), 3, fp);
    fread( mghInfor->y_ras, sizeof(float), 3, fp);
    fread( mghInfor->z_ras, sizeof(float), 3, fp);
    fread( mghInfor->c_ras, sizeof(float), 3, fp);
    if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0)) {
      for(i=0; i<3; i++) {
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
  nv = mghInfor->width * mghInfor->height * mghInfor->depth *mghInfor->nframes;
//  printf("%s: width = %d, height = %d, depth = %d, nframes = %d, num values = %d\n",
//    progname,mghInfor->width,mghInfor->height,mghInfor->depth,mghInfor->nframes,nv);

  if(mghInfor->type ==MRI_UCHAR) {
    //cout << "Type: unsigned char" << endl;
    fseek(fp, 0, SEEK_END);
    fileLen = ftell(fp);
    headerSize = fileLen - nv*sizeof(unsigned char);
    mghInfor->charImageArray = (unsigned char *)calloc(nv, sizeof(unsigned char));
    fseek(fp, headerSize-16, SEEK_SET) ;
    fread(mghInfor->charImageArray, sizeof(unsigned char), nv, fp);
  } else if(mghInfor->type ==MRI_INT) {
    //cout << "Type: int" << endl;
    fseek(fp, 0, SEEK_END);
    fileLen = ftell(fp);
    headerSize = fileLen - nv*sizeof(int);
    mghInfor->intImageArray = (int *)calloc(nv, sizeof(int));
    fseek(fp, headerSize-16, SEEK_SET) ;
    fread(mghInfor->intImageArray, sizeof(int), nv, fp);
    if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0)) {
      for(i=0; i<nv; i++) {
        bswap(mghInfor->intImageArray[i]);
      }
    }
  } else if(mghInfor->type ==MRI_LONG) {
    //cout << "Type: long" << endl;
    fseek(fp, 0, SEEK_END);
    fileLen = ftell(fp);
    headerSize = fileLen - nv*sizeof(MRI_LONG);
    mghInfor->longImageArray = (long *)calloc(nv, sizeof(long));
    fseek(fp, headerSize-16, SEEK_SET) ;
    fread(mghInfor->longImageArray, sizeof(long), nv, fp);
    if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0)) {
      for(i=0; i<nv; i++) {
        bswap(mghInfor->longImageArray[i]);
      }
    }
  } else if(mghInfor->type ==MRI_FLOAT) {
    //cout << "Type: float" << endl;
    mghInfor->floatImageArray = (float *)calloc(nv, sizeof(float));
    fseek(fp, 0, SEEK_END) ;
    fileLen = ftell(fp);
    headerSize = fileLen - nv*sizeof(MRI_FLOAT);
    fseek(fp, headerSize-16, SEEK_SET) ;
    fread(mghInfor->floatImageArray, sizeof(float), nv, fp);
    if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0)) {
      for(i=0; i<nv; i++) {
        bswap(mghInfor->floatImageArray[i]);
      }
    }
  } else if(mghInfor->type ==MRI_SHORT) {
    //cout << "Type: short" << endl;
    fseek(fp, 0, SEEK_END) ;
    fileLen = ftell(fp);
    headerSize = fileLen - nv*sizeof(short);
    mghInfor->shortImageArray = (short *)calloc(nv, sizeof(short));
    fseek(fp, headerSize-16, SEEK_SET) ;
    fread(mghInfor->shortImageArray, sizeof(short), nv, fp);
    if ( (endianFlag == 0 && endianInfor == 1) || (endianFlag == 1 && endianInfor == 0)) {
      for(i=0; i<nv; i++) {
        bswap(mghInfor->shortImageArray[i]);
      }
    }
  } else {
    printf("%s: ERROR: mgh type %d not defined\n",progname,mghInfor->type);
    exit(1);    
//    return (-2);
  }
  fclose(fp);
  return 1;
}

FIBER* tractomain::ReadFiber(string fname,int *numfibers)
{
  FIBER*            fiber_array;
  stFiberFileHeader fiber_hdr;
  int i,j;

  FILE* fp = fopen(fname.c_str(), "r" );
  if (!fp) {
    printf("%s: ERROR: failed to open file %s\n",progname,fname.c_str());
    exit(1);
  }
  fread(&fiber_hdr,sizeof(fiber_hdr),1,fp);      // read file header

  // save fibers one-by-one, starts from file_offset = 128;
  fseek(fp, 128, SEEK_SET);

  fiber_array = (FIBER *)calloc (fiber_hdr.nFiberNr, sizeof(FIBER));

  for(i=0; i<fiber_hdr.nFiberNr; i++) {
    fread(&fiber_array[i].nLength, sizeof(int), 1, fp);                            
    fread(&fiber_array[i].nSelStatus, sizeof(unsigned char), 1, fp);                          
    fread(&fiber_array[i].rgbFiberClr, sizeof(RGB_TRIPLE), 1, fp);                            
    fread(&fiber_array[i].nSelBeginIdx, sizeof(int), 1, fp);                            
    fread(&fiber_array[i].nSelEndIdx, sizeof(int), 1, fp);                           
    fiber_array[i].pxyzChain = (XYZ_TRIPLE *)calloc(fiber_array[i].nLength, sizeof(XYZ_TRIPLE));
    for(j = 0; j<fiber_array[i].nLength; j++) {
      fread(&fiber_array[i].pxyzChain[j], sizeof(XYZ_TRIPLE), 1, fp); 
    }
  }
  fclose(fp);
  *numfibers = fiber_hdr.nFiberNr;
  return fiber_array;
} 

/***************************************************************************/
// Open GL routines for rendering images

void glRenderImage(imageStruct imageinfo, float WL, float WW)
{
  float* bgimage = imageinfo.ImageArray;

  int width = imageinfo.width;
  int height = imageinfo.height;
  int slices = imageinfo.depth;

  float axialSliceLoc = imageinfo.AxialSliceLocation - imageinfo.depth/2; 
  float sagittalSliceLoc = -(imageinfo.width/2 - imageinfo.SagittalSliceLocation);
  float coronalSliceLoc =  -(imageinfo.height/2- imageinfo.CoronalSliceLocation);

  int totalFiberSheafNum = imageinfo.fiberNum;
  int *totalFiberNum = imageinfo.totalFiberNum;
  RGB_TRIPLE *fiberColor = imageinfo.FiberColor;
  FIBER **Fibers = imageinfo.pasFiberAry;

  float ratioXZ = imageinfo.voxel_t/imageinfo.voxel_w;
  float halfRatioXZ = ratioXZ/2;
  float ratioYZ = imageinfo.voxel_t/imageinfo.voxel_h;
  float halfRatioYZ = ratioYZ/2;
  float zLocation = axialSliceLoc*ratioXZ;

  float r, g, b, a;	
  int i, j, k, l, m;
  int ns, np;
  float grayScale;
  static int distFound = 0;
  float xc, yc, zc;
  ns = width * height;
  np = ns * slices;
  int shift_xc = (int)((float)width/2+(float)0.0);
  int shift_yc = (int)((float)height/2+(float)0.0);
  int shift_zc = (int)((float)slices/2+(float)0.0);
  float *linePoint;
  FILE *fp;
  GLfloat alphaAxial = image_alpha;
  GLfloat alphaCoronal = image_alpha;
  GLfloat alphaSagittal = image_alpha;
  int offset, index, index1, index2;
  GLfloat maxGrayLevel = 2.0f;
  GLfloat upper = WL + WW/2;
  GLfloat lower = WL - WW/2;
  if(lower<0) lower = 0;
  if(upper>2.0) upper = 2.0;

  // draw fibers
  glPushMatrix();

  /************************
  Flip the fiber direction as checking the x, y and z flipping flag.
  ************************/
  int ii;
  if(zflipflag != current_zflipflag) {
    for(ii=0; ii<totalFiberSheafNum; ii++) {
      for(i=0; i<totalFiberNum[ii]; i++) {
        for(j = 0; j<Fibers[ii][i].nLength; j++) {
          Fibers[ii][i].pxyzChain[j].z = slices - Fibers[ii][i].pxyzChain[j].z -1;
        }
      }
    }
    current_zflipflag = zflipflag;
  }
  if(xflipflag != current_xflipflag) {
    for(ii=0; ii<totalFiberSheafNum; ii++) {
      for(i=0; i<totalFiberNum[ii]; i++) {
        for(j = 0; j<Fibers[ii][i].nLength; j++) {
          Fibers[ii][i].pxyzChain[j].x =  height - Fibers[ii][i].pxyzChain[j].x -1;
        }
      }
    }
    current_xflipflag = xflipflag;
  }
  if(yflipflag != current_yflipflag) {
    for(ii=0; ii<totalFiberSheafNum; ii++) {
      for(i=0; i<totalFiberNum[ii]; i++) {
        for(j = 0; j<Fibers[ii][i].nLength; j++) {
          Fibers[ii][i].pxyzChain[j].y =  width - Fibers[ii][i].pxyzChain[j].y -1;
        }
      }
    }
    current_yflipflag = yflipflag;
  }

  for(ii=0; ii<totalFiberSheafNum; ii++) {
    for(i=0; i<totalFiberNum[ii]; i++) {
      for(j = Fibers[ii][i].nSelBeginIdx; j<Fibers[ii][i].nSelEndIdx-1; j++) {
        glBegin(GL_LINE_LOOP);
        
        // draw line connecting two points
        for(k=j; k<=j+1; k++) {
          r = fiberPointRGBA[ii][i][k][0];
          g = fiberPointRGBA[ii][i][k][1];
          b = fiberPointRGBA[ii][i][k][2];
          a = fiberPointRGBA[ii][i][k][3];
          glColor4f(r,g,b,a);
          glVertex4f(Fibers[ii][i].pxyzChain[k].y - shift_yc, 
            Fibers[ii][i].pxyzChain[k].x - shift_xc, 
            Fibers[ii][i].pxyzChain[k].z - slices + shift_zc, 
            1.0);
        }
        glEnd( );
      }
    }
  }
  glPopMatrix();

  /* Contrast control*/
  float *gray = (float*)calloc(np, sizeof(float));
  for ( i = 0; i<np; i++) {
    if(bgimage[i] >= upper) {
      gray[i] = maxGrayLevel;
    } else if (bgimage[i] <= lower) {
      gray[i] = 0.0f;
    } else {
      gray[i] = maxGrayLevel - maxGrayLevel * (upper-bgimage[i])/WW;
    }
  }

  /* Brightness control*/
  if (axialViewFlag) {
    for ( i = -shift_yc+1; i<shift_yc-1; i ++) {
      for ( j =-shift_xc+1; j<shift_xc-1; j++) {
        if(bgimage[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc)*width+j+shift_xc] == 0) {
          continue;
        }
        glBegin(GL_QUADS);	
        glColor4f(gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc)*width+j+shift_xc],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc)*width+j+shift_xc],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc)*width+j+shift_xc],
          alphaAxial);
        glVertex3f(i-0.5, j-0.5, zLocation);

        if(distFound == 0) {
          xc = i-0.5;  yc = j-0.5;  zc = zLocation;
          if ( maxDist < sqrt(xc*xc + yc*yc + zc*zc) ) {
            maxDist = sqrt(xc*xc + yc*yc + zc*zc);
          }
        }
        glColor4f(gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc+1)*width+j+shift_xc],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc+1)*width+j+shift_xc],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc+1)*width+j+shift_xc],
          alphaAxial);
        glVertex3f(i+0.5, j-0.5, zLocation);

        if(distFound == 0) {
          xc = i+0.5;  yc = j-0.5;  zc = zLocation;
          if ( maxDist < sqrt(xc*xc + yc*yc + zc*zc) ) {
            maxDist = sqrt(xc*xc + yc*yc + zc*zc);
          }
        }
        glColor4f(gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc+1)*width+j+shift_xc+1],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc+1)*width+j+shift_xc+1],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc+1)*width+j+shift_xc+1],
          alphaAxial);
        glVertex3f(i+0.5, j+0.5, zLocation);

        glColor4f(gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc)*width+j+shift_xc+1],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc)*width+j+shift_xc+1],
          gray[(shift_zc+(int)axialSliceLoc)*ns+(i+shift_yc)*width+j+shift_xc+1],
          alphaAxial);
        glVertex3f(i-0.5, j+0.5, zLocation);
        glEnd();
      }
    }
  }
  distFound = 1;

  // Draw sagittal slice
  if(sagittalViewFlag) {
    for ( i = -shift_zc; i<shift_zc-2; i ++) {
      for ( j =-shift_yc; j<shift_yc-1; j++) {
        if (bgimage[(i+shift_zc)*ns+(j+shift_yc)*width+shift_xc+(int)sagittalSliceLoc] == 0) {
          continue;
        }
        glBegin(GL_QUADS);
        glColor4f(gray[(i+shift_zc)*ns+(j+shift_yc)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc)*ns+(j+shift_yc)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc)*ns+(j+shift_yc)*width+shift_xc+(int)sagittalSliceLoc],
          alphaSagittal);
        glVertex3f(j-0.5, sagittalSliceLoc, i*ratioXZ-halfRatioXZ);

        glColor4f(gray[(i+shift_zc)*ns+(j+shift_yc+1)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc)*ns+(j+shift_yc+1)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc)*ns+(j+shift_yc+1)*width+shift_xc+(int)sagittalSliceLoc],
          alphaSagittal);
        glVertex3f(j+0.5, sagittalSliceLoc,  i*ratioXZ-halfRatioXZ);

        glColor4f(gray[(i+shift_zc+1)*ns+(j+shift_yc+1)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc+1)*ns+(j+shift_yc+1)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc+1)*ns+(j+shift_yc+1)*width+shift_xc+(int)sagittalSliceLoc],
          alphaSagittal);
        glVertex3f(j+0.5, sagittalSliceLoc, i*ratioXZ+halfRatioXZ);

        glColor4f(gray[(i+shift_zc+1)*ns+(j+shift_yc)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc+1)*ns+(j+shift_yc)*width+shift_xc+(int)sagittalSliceLoc],
          gray[(i+shift_zc+1)*ns+(j+shift_yc)*width+shift_xc+(int)sagittalSliceLoc],
          alphaSagittal);
        glVertex3f(j-0.5, sagittalSliceLoc, i*ratioXZ+halfRatioXZ);
        glEnd();
      }
    }
  }

  // Draw CORONAL slice
  if(coronalViewFlag) {
    for ( i = -shift_zc; i<shift_zc-2; i++) {
      for ( j = -shift_xc; j<shift_xc-2; j++) {
        if (bgimage[(i+shift_zc)*ns + (shift_yc+(int)coronalSliceLoc)*width + j+shift_xc] == 0) {
          continue;
        }
        glBegin(GL_POLYGON);
        glColor4f(gray[(i+shift_zc)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc],
          gray[(i+shift_zc)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc],
          gray[(i+shift_zc)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc],
          alphaCoronal);
        glVertex3f(coronalSliceLoc, j-0.5, i*ratioYZ-halfRatioYZ);

        glColor4f(gray[(i+shift_zc)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc+1],
          gray[(i+shift_zc)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc+1],
          gray[(i+shift_zc)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc+1],
          alphaCoronal);
        glVertex3f(coronalSliceLoc, j+0.5,  i*ratioYZ-halfRatioYZ);

        glColor4f(gray[(i+shift_zc+1)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc+1],
          gray[(i+shift_zc+1)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc+1],
          gray[(i+shift_zc+1)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc+1],
          alphaCoronal);
        glVertex3f(coronalSliceLoc, j+0.5, i*ratioYZ+halfRatioYZ);

        glColor4f(gray[(i+shift_zc+1)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc],
          gray[(i+shift_zc+1)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc],
          gray[(i+shift_zc+1)*ns+(shift_yc+(int)coronalSliceLoc)*width+j+shift_xc],
          alphaCoronal);
        glVertex3f(coronalSliceLoc, j-0.5, i*ratioYZ+halfRatioYZ);
        glEnd();
      }
    }
  }
  free(gray);
}
                            
void glClearViewport( float c1R, float c1G, float c1B, 
                      float c2R, float c2G, float c2B, bool glclear)
{
  // clear the viewport...
  if (glclear) {
     glClearColor( c1R, c1G, c1B, 1.0f );
     glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  }

  glPushAttrib( GL_ENABLE_BIT );
  glDisable( GL_DEPTH_TEST );
  glDisable( GL_LIGHTING );

  glMatrixMode( GL_PROJECTION );
  glPushMatrix();
  glLoadIdentity();                     
  glOrtho( 0, 1, 0, 1, 0, 1 );

  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  glLoadIdentity();      

  // draw over the entire screen
  const float depth = -0.9999999f;
  glBegin( GL_TRIANGLE_STRIP );
  glColor3f( c1R, c1G, c1B );
  glVertex3f( 0, 1, depth );

  glColor3f( c2R, c2G, c2B );
  glVertex3f( 0, 0, depth );

  glColor3f( c1R, c1G, c1B );
  glVertex3f( 1, 1, depth );

  glColor3f( c2R, c2G, c2B );
  glVertex3f( 1, 0, depth );
  glEnd();

  glMatrixMode( GL_MODELVIEW );
  glPopMatrix();
  glMatrixMode( GL_PROJECTION );
  glPopMatrix();
  glPopAttrib();
}

void glRender( const Light& light )
{
  int lightNumber = GL_LIGHT0;

  // figure out which GL light number the light is...
  switch (light.number()) {
    case 0: lightNumber = GL_LIGHT0; break;
    case 1: lightNumber = GL_LIGHT1; break;
    case 2: lightNumber = GL_LIGHT2; break;
    case 3: lightNumber = GL_LIGHT3; break;
    case 4: lightNumber = GL_LIGHT4; break;
    case 5: lightNumber = GL_LIGHT5; break;
    case 6: lightNumber = GL_LIGHT6; break;
    case 7: lightNumber = GL_LIGHT7; break;
    default:
    // light number has to be in the range [0,7] 
    // when rendering with OpenGL.
    cout<<"light number "<<light.number()<<" has to be in the range [0,7]\n"<<flush;
    assert(light.number()>=0 && light.number()<=7);
    lightNumber = GL_LIGHT0; break;
  }

  if (light.isOn()) {
    // gather the data to be used to setup the opengl light.
    float ambient[3], diffuse[3], specular[3], position[4];
    float constantAttenuation, linearAttenuation, quadraticAttenuation;
    float cutoff, direction[3], exponent;
    light.getColor( Light::ambient, ambient[0], ambient[1], ambient[2] );
    light.getColor( Light::diffuse, diffuse[0], diffuse[1], diffuse[2] );
    light.getColor( Light::specular, specular[0], specular[1], specular[2] );
    light.getPos( position[0], position[1], position[2], position[3] );
    light.getAtten( constantAttenuation, linearAttenuation, quadraticAttenuation );
    light.getSpotDir( direction[0], direction[1], direction[2] );
    light.getSpotCone( exponent, cutoff );

    // light color.
    glLightfv(lightNumber, GL_AMBIENT, ambient);
    glLightfv(lightNumber, GL_DIFFUSE, diffuse);
    glLightfv(lightNumber, GL_SPECULAR, specular);

    // position
    glLightfv(lightNumber, GL_POSITION, position);

    // attenuation
    glLightfv(lightNumber, GL_CONSTANT_ATTENUATION, &constantAttenuation);
    glLightfv(lightNumber, GL_LINEAR_ATTENUATION, &linearAttenuation);
    glLightfv(lightNumber, GL_QUADRATIC_ATTENUATION, &quadraticAttenuation);

    // spotlight
    glLightfv( lightNumber, GL_SPOT_CUTOFF, &cutoff );
    glLightfv( lightNumber, GL_SPOT_DIRECTION, direction );
    glLightfv( lightNumber, GL_SPOT_EXPONENT, &exponent );

    glEnable( lightNumber );

    const float thisIsFalse = GL_FALSE;
    const float thisIsTrue = GL_TRUE;
    glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, &thisIsFalse);
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, &thisIsFalse);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
  } else {
    glDisable( lightNumber );
  }
}



/*******************************************************************************/


void tractomain::drawText()
{
  glPushAttrib( GL_ENABLE_BIT );

  //: Draw the text overlay in ortho mode...
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode( GL_PROJECTION );
  glPushMatrix();
  glLoadIdentity();

  glOrtho( -50,50,-50,50,0,500 );
  char buff[256], buff2[256];
  glDisable(GL_LIGHTING);

  // Display status line
  if( beacon.timeLeft > 0 ) {
    glColor4f(1,1,1,1);
    glListBase(fontBase);
    glRasterPos3f( -48, 48, 0);
    glCallLists(strlen(beacon.text), GL_BYTE, beacon.text);
    --beacon.timeLeft;
  }

  // Display usage
  sprintf(buff2,"Keys:  T - mouse help / A - axial flag / C - coronal flag / S - sagittal flag / Q - Exit");
  glColor4f(.7, .7, 1, 1);
  glListBase(fontBase);
  glRasterPos3f(-48,-49,0);
  glCallLists(strlen(buff2), GL_BYTE, buff2);

  glPopMatrix();

  glMatrixMode( GL_MODELVIEW );
  glPopMatrix();

  glPopAttrib();
}

void tractomain::OnContextDraw()
{
  int xx, yy, zz, status;

  printf("%s: preparing to render image...\n",progname);

  if(stepNum >= STEPNUM) // check the angle of auto rotation over 360 degree
  stepNum = 0;
  //get original section view mode
  oldSFlag=SFlag; oldIFlag=IFlag; oldRFlag=RFlag;
  oldLFlag=LFlag; oldAFlag=AFlag; oldPFlag=PFlag;

  //update the viewport dimensions
  glViewport( 0, 0, this->width(), this->height() );

  glClearViewport( 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f );

  if(!savetifflag) drawText();

  //update the perspective matrix
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPerspective( 20.0f, ((float)this->width())/((float)this->height()), nearPlane, farPlane);

  // Affect the "modelview" matrix from now on.
  glMatrixMode( GL_MODELVIEW );
  // Clear the "modelview" matrix
  glLoadIdentity();

  //render the light
  glRender( light );
  glRender( light2 );

  // transform scene by the TrackBall
  if(InitialFlag == 0) {
    InitialFlag = 1;
    cout << progname << ": setting view direction..." << endl;
    mm = trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).data();
    mm[14] = mm[14]-350.0;
  } else {
    mm[14] = mm[14];
  }

  if (movieXFlag == 1) {
    if(SFlag || IFlag || RFlag || LFlag || AFlag || PFlag ) {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(stepNum*rotationStep[stepNum], 1.0f, 0.0f, 0.0f);
    } else {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(rotationStep[stepNum], 1.0f, 0.0f, 0.0f);
    }
    mm = trackBall->getMatrix().data();
    if(mm[14] == 0.0) mm[14] = mm[14]-350.0;
    stepNum++;
  }
  if (movieXFlag == -1) {
    if(SFlag || IFlag || RFlag || LFlag || AFlag || PFlag ) {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(stepNum*rotationStep[stepNum], -1.0f, 0.0f, 0.0f);
    } else {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(rotationStep[stepNum], -1.0f, 0.0f, 0.0f);
    }
    mm = trackBall->getMatrix().data();
    if(mm[14] == 0.0) mm[14] = mm[14]-350.0;
    stepNum++;
  }
  if (movieYFlag == 1) {
    if(SFlag || IFlag || RFlag || LFlag || AFlag || PFlag ) {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(stepNum*rotationStep[stepNum], 0.0f, 1.0f, 0.0f);
    } else {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(rotationStep[stepNum], 0.0f, 1.0f, 0.0f);
    }
    mm = trackBall->getMatrix().data();
    if(mm[14] == 0.0) mm[14] = mm[14]-350.0;
    stepNum++;
  }
  if (movieYFlag == -1) {
    if(SFlag || IFlag || RFlag || LFlag || AFlag || PFlag ) {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(stepNum*rotationStep[stepNum], 0.0f, -1.0f, 0.0f);
    } else {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(rotationStep[stepNum], 0.0f, -1.0f, 0.0f);
    }
    mm = trackBall->getMatrix().data();
    if(mm[14] == 0.0) mm[14] = mm[14]-350.0;
    stepNum++;
  }
  if (movieZFlag == -1) {
    if(SFlag || IFlag || RFlag || LFlag || AFlag || PFlag ) {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(stepNum*rotationStep[stepNum], 0.0f, 0.0f, -1.0f);
    } else {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(rotationStep[stepNum], 0.0f, 0.0f, -1.0f);
    }
    mm = trackBall->getMatrix().data();
    if(mm[14] == 0.0) mm[14] = mm[14]-350.0;
    stepNum++;    
  }
  if (movieZFlag == 1) {
    if(SFlag || IFlag || RFlag || LFlag || AFlag || PFlag ) {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(stepNum*rotationStep[stepNum], 0.0f, 0.0f, 1.0f);
    } else {
      trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).rotate(rotationStep[stepNum], 0.0f, 0.0f, 1.0f);
      mm = trackBall->getMatrix().data();
      if(mm[14] == 0.0) mm[14] = mm[14]-350.0;
      stepNum++;
    }
  }
  if( (SFlag || IFlag || RFlag || LFlag || AFlag || PFlag ) && movieXFlag == 0 && movieYFlag == 0 && movieZFlag == 0 ) {
      mm = trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).data();
      mm[14] = mm[14]-350.0;
  }
  if(oldxrot != xrot || oldyrot != yrot || oldzrot != zrot) {
    SFlag = 0; IFlag = 0; RFlag = 0;
    LFlag = 0; AFlag = 0; PFlag = 0;
    muiSetActive(SRadio, 0);
    muiSetActive(IRadio, 0);
    muiSetActive(RRadio, 0);
    muiSetActive(LRadio, 0);
    muiSetActive(ARadio, 0);
    muiSetActive(PRadio, 0);
  }
  if(oldxrot != xrot) {   
    if(xrot > STEPNUM) xrot = xrot - STEPNUM ;
    float tmpxrot =  xrot/180.0 * PI;
    trackBall->getMatrix().rotate(tmpxrot, 1.0f, 0.0f, 0.0f);
    oldxrot = xrot;
  }
  if(oldyrot != yrot) {
    if(yrot > STEPNUM) yrot = yrot - STEPNUM ;
    float tmpyrot =  yrot/180.0 * PI;
    trackBall->getMatrix().rotate(tmpyrot, 0.0f, 1.0f, 0.0f);
    oldyrot = yrot;
  }
  if(oldzrot != zrot) {
    if(zrot > STEPNUM) zrot = zrot - STEPNUM ;
    float tmpzrot =  zrot/180.0 * PI;
    trackBall->getMatrix().rotate(tmpzrot, 0.0f, 0.0f, 1.0f);
    oldzrot = zrot;
  }
  glMultMatrixf( mm );

  glRotatef(90.0f, 0.0f, 0.0f, -1.0f);

  printf("%s: rendering image...\n",progname);
  glRenderImage(this->imageinfo, WL, WW);

  if(contrastInitFlag == 1) WW = 2*(1-contrastPara);
  if(brightnessInitFlag == 1) WL = (1-brightnessPara);

  if (slideraxial != oldslideraxial) {
    glutSetWindow(2);
    this->imageinfo.AsliceLoc = (int)((float)this->imageinfo.depth * slideraxial );
    this->imageinfo.AxialSliceLocation = this->imageinfo.AsliceLoc;
    float brainArea = (this->imageinfo.voxel_t*this->imageinfo.depth)/2+this->imageinfo.voxel_t;
    float emptyArea =  (this->imageinfo.height - (int)brainArea)/2;
    CoronalLineH = (float)this->imageinfo.height/2-emptyArea -  \
      (this->imageinfo.AxialSliceLocation * this->imageinfo.voxel_t/2+this->imageinfo.voxel_t);
    CoronalLineL = CoronalLineL;
    SagittalLineH =  CoronalLineH;
    SagittalLineL = SagittalLineL;
    AxCoSaRedraw();
  }
  if (slidercoronal != oldslidercoronal) {
    glutSetWindow(3);
    this->imageinfo.CsliceLoc = (int)((float)this->imageinfo.height * slidercoronal );
    this->imageinfo.CoronalSliceLocation = this->imageinfo.CsliceLoc;
    AxialLineH=this->imageinfo.height/2 - this->imageinfo.CoronalSliceLocation  ;
    SagittalLineL = this->imageinfo.CoronalSliceLocation - this->imageinfo.height/2; 
    AxCoSaRedraw();
  }
  if (slidersagittal != oldslidersagittal) {
    glutSetWindow(4);
    this->imageinfo.SsliceLoc = (int)((float)this->imageinfo.width * slidersagittal );
    this->imageinfo.SagittalSliceLocation = this->imageinfo.SsliceLoc;
    AxialLineL=this->imageinfo.width/2 - this->imageinfo.SagittalSliceLocation;
    CoronalLineL = this->imageinfo.width/2- this->imageinfo.SagittalSliceLocation; 
    AxCoSaRedraw();
  }

  if (savetifflag && drawcount>min_drawcount) {
    printf("%s: saving image as file %s...\n",
      progname,fname_tif.c_str());
    printf("%s: swapping buffers...\n",progname);
    glutSwapBuffers();
    status = save_tiff(fname_tif.c_str());
    if(status) {
      printf("%s: failed to save as tif\n",progname);
    }
    exit(status);
  } else {
    printf("%s: not saving image\n",progname);
  }

  drawcount++;
}

void tractomain::OnContextInit()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glShadeModel(GL_SMOOTH);
    glEnable( GL_BLEND );

    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA ); 
    glOrtho(-this->imageinfo.width, this->imageinfo.width,  \
            -this->imageinfo.height, this->imageinfo.height, \
            -this->imageinfo.depth, this->imageinfo.depth);

    //set up font in display list
//    cout<<"Loading Helvetica Font\n\n"<<flush;
    fontBase = glGenLists(128);
    for (int i=0; i<128; i++)
    {
       glNewList(fontBase+i, GL_COMPILE);
       glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, i);
       glEndList();
    }
}

//: Window reshape
//  Called on reshape of a window.
//  Put application code in here.
void tractomain::OnAppReshape( int width, int height )
{
   this->setWidth( width );
   this->setHeight( height );

   //update the trackBall
   hyperBall.setViewport( TrackBall::top, 0 );
   hyperBall.setViewport( TrackBall::right, width );
   hyperBall.setViewport( TrackBall::left, 0 );
   hyperBall.setViewport( TrackBall::bottom, height );
}

//: Window reshape for each GL context (window)
//  Called on reshape of each gl window.
//  Put GL code in here.
void tractomain::OnContextReshape( int width, int height )
{
}

/*******************************************************************************/

void tractomain::OnKeyboardDown( unsigned char k, int x, int y )
{
  static bool fullscreen = false;
  static float windowedX = 0, windowedY, windowedWidth, windowedHeight;

  switch(k) {
    case 'f':
      if(fullscreen) {
        glutPositionWindow(int(windowedX), int(windowedY));
        glutReshapeWindow(int(windowedWidth), int(windowedHeight));
        fullscreen = false;
        sprintf(beacon.text,"Windowed Mode");
      } else {
        windowedX = glutGet(GLUT_WINDOW_X);
        windowedY = glutGet(GLUT_WINDOW_Y);
        windowedWidth = glutGet(GLUT_WINDOW_WIDTH);
        windowedHeight = glutGet(GLUT_WINDOW_HEIGHT);
        glutFullScreen();
        fullscreen = true;
        sprintf(beacon.text,"Fullscreen Mode");
      }
      break;

    case 't':
      sprintf(beacon.text,"Mouse: Right - Rotate / CTRL-Right - Translate XY / SHIFT-Right - Translate Z");
      break;

    case 'a':
      if(axialViewFlag == 0) {
        axialViewFlag = 1;
      } else {
        axialViewFlag = 0;
      }
      sprintf(beacon.text,"axialViewFlag = %d",axialViewFlag);
      break;
            
    case 's':
      if(sagittalViewFlag == 0) {
        sagittalViewFlag = 1;
      } else {
        sagittalViewFlag = 0;
      }
      sprintf(beacon.text,"sagittalViewFlag = %d",sagittalViewFlag);
      break;

    case 'c':
      if(coronalViewFlag == 0) {
        coronalViewFlag = 1;
      } else {
        coronalViewFlag = 0;
      }
      sprintf(beacon.text,"coronalViewFlag = %d",coronalViewFlag);
      break;
  }

  beacon.timeLeft = 60;
  if(k=='q' || k==27) 
  exit(0);
}

void tractomain::OnKeyboardUp(unsigned char k, int x, int y)
{
   beacon.timeLeft = 60;
}

void tractomain::OnSpecialKeyboardDown(int k, int x, int y)
{
}

void tractomain::OnSpecialKeyboardUp(int k, int x, int y)
{
}

/*******************************************************************************/

////////////////////////////////
//: This is called on any mouse event
//  override this, and use the _mouse 
//  member to get current mouse state.
////////////////////////////////
void tractomain::OnMouseEvent4()
{
    cout << "on mouse event4" << endl;
}
void tractomain::OnMouseEvent3()
{
    cout << "on mouse event3" << endl;
}
void tractomain::OnMouseEvent2()
{
    cout << "on mouse event2" << endl;
}
void tractomain::OnMouseEvent()
{
  GlutKernelApplication::OnMouseEvent();
  if(mouse().rightEdgeState() == Mouse::DOWN)
  {
    trackBall->setRotation(trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag));
    if( this->modifier() == GLUT_ACTIVE_CTRL ) {
      trackBall->applyXYtrans(mouse().x(), mouse().y(), mouse().xOld(), mouse().yOld());
      mm = trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).data();
      mm[14] = mm[14]-350.0;
    } else if( this->modifier() == GLUT_ACTIVE_SHIFT ) {
      trackBall->applyZtrans(mouse().y(), mouse().yOld());
      mm = trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).data();
      mm[14] = mm[14]-350.0;
    } else {
      trackBall->applyGeneralRot(mouse().x(), mouse().y(), mouse().xOld(), mouse().yOld());
      mm = trackBall->getMatrix(SFlag, IFlag, RFlag, LFlag, AFlag, PFlag).data();
      mm[14] = mm[14]-350.0;
    }
    beacon.timeLeft = 60;
  }
}

// called when mouse changes position
void tractomain::OnMousePos( int x, int y )
{
   GlutKernelApplication::OnMousePos( x, y );
   // Set the position, then update
    //mouse().setPosition( xpos, ypos );
    //mouse().updateEdgeStates();
}

// called when mouse button is clicked
void tractomain::OnMouseClick2( int a, int b, int c, int d )
{
   GlutKernelApplication::OnMouseClick2( a, b, c, d );
}

void tractomain::OnMouseClick3( int a, int b, int c, int d )
{
   GlutKernelApplication::OnMouseClick3( a, b, c, d );
}

void tractomain::OnMouseClick4( int a, int b, int c, int d )
{
   GlutKernelApplication::OnMouseClick4( a, b, c, d );
}

void tractomain::OnMouseClick( int a, int b, int c, int d )
{
   SFlag=0; IFlag=0; RFlag=0;
   LFlag=0; AFlag=0; PFlag=0;
   movieXFlag = 0;
   movieYFlag = 0;
   movieZFlag = 0;
   GlutKernelApplication::OnMouseClick( a, b, c, d );
}

/*******************************************************************************/


void flipFiberXDir( muiObject *obj, enum muiReturnValue rv )
{
  //cout << muiGetActive(obj) << endl;
  xflipflag = muiGetActive(obj);
}
void flipFiberYDir( muiObject *obj, enum muiReturnValue rv )
{
  //cout << muiGetActive(obj) << endl;
  yflipflag = muiGetActive(obj);
}
void flipFiberZDir( muiObject *obj, enum muiReturnValue rv )
{
  //cout << muiGetActive(obj) << endl;
  zflipflag = muiGetActive(obj);
}

void changeContrast( muiObject *obj, enum muiReturnValue rv )
{
  contrastInitFlag = 1;
  contrastPara  = muiGetHSVal( obj);
}
void changeBrightness( muiObject *obj, enum muiReturnValue rv )
{
  brightnessInitFlag = 1;
  brightnessPara =  muiGetHSVal( obj);
}

void getxrot( muiObject *obj, enum muiReturnValue rv )
{
    oldxrot = xrot;
    xrot = atof(muiGetTBString(obj));
    if(oldxrot == xrot)
        oldxrot = xrot - 1;
    SFlag=0; IFlag=0; RFlag=0;
    LFlag=0; AFlag=0; PFlag=0;
}
void getyrot( muiObject *obj, enum muiReturnValue rv )
{
    oldyrot = yrot;
    yrot = atof(muiGetTBString(obj));
    if(oldyrot == yrot)
        oldyrot = yrot - 1;
    SFlag=0; IFlag=0; RFlag=0;
    LFlag=0; AFlag=0; PFlag=0;
}
void getzrot( muiObject *obj, enum muiReturnValue rv )
{
    oldzrot = zrot;
    zrot = atof(muiGetTBString(obj));
    if(oldzrot == zrot)
        oldzrot = zrot - 1;
    SFlag=0; IFlag=0; RFlag=0;
    LFlag=0; AFlag=0; PFlag=0;
}

/*******************************************************************************/


void tractomain::bswap( char* p1, char* p2 )
{
  char t = *p1; *p1 = *p2; *p2 = t;
}

void tractomain::bswap( short& n ){
  bswap( 0+(char*)&n, 1+(char*)&n );
}

void tractomain::bswap( int& n ){
  bswap( 0+(char*)&n, 3+(char*)&n );
  bswap( 1+(char*)&n, 2+(char*)&n );
}

void tractomain::bswap( long& n ){
  bswap( 0+(char*)&n, 3+(char*)&n );
  bswap( 1+(char*)&n, 2+(char*)&n );
}

void tractomain::bswap( float& n ){
  bswap( 0+(char*)&n, 3+(char*)&n );
  bswap( 1+(char*)&n, 2+(char*)&n );
}

void tractomain::bswap( double& n ){
  bswap( 0+(char*)&n, 7+(char*)&n );
  bswap( 1+(char*)&n, 6+(char*)&n );
  bswap( 2+(char*)&n, 5+(char*)&n );
  bswap( 3+(char*)&n, 4+(char*)&n );
}

/*******************************************************************************/

void tractomain::AxCoSaRedraw( )
{
  glutSetWindow(2);
  this->AxialDisplay();

  glutSetWindow(3);
  this->CoronalDisplay();

  glutSetWindow(4);
  this->SagittalDisplay();
}

void tractomain::OperationPanelEvent( )
{
  char as[32], cs[32], ss[32];
  glutSetWindow(5);
  int PanelWidth = glutGet(GLUT_WINDOW_WIDTH);
  int PanelHeight = glutGet(GLUT_WINDOW_HEIGHT);
  int x, y;

  muiNewBoldLabel(5, 660, "Rotation Angle ( Degree )");
  muiNewLabel(10, 635, "X Dir:");
  muiNewLabel(115, 635, "Y Dir:");
  muiNewLabel(220, 635, "Z Dir:");
  xrotTB = muiNewTextbox(45, 100, 625);
  yrotTB = muiNewTextbox(150, 205, 625);
  zrotTB = muiNewTextbox(255, 310, 625);
  muiSetTBString(xrotTB, "0");
  muiSetTBString(yrotTB, "0");
  muiSetTBString(zrotTB, "0");
  muiSetCallback(xrotTB, getxrot);
  muiSetCallback(yrotTB, getyrot);
  muiSetCallback(zrotTB, getzrot);

/*
  muiNewBoldLabel(5, 600, "Fiber Visualization Effect");
  FiberViewEffectRadio = muiNewRadioButton( 10, 570 ) ;
  muiNewLabel(40, 575, "Glossy");
  muiSetActive(FiberViewEffectRadio, glossyFlag);
  muiSetCallback(FiberViewEffectRadio, fiberGlossyView);
  FiberViewOpacityRadio = muiNewRadioButton( 80, 570 ) ;
  muiNewLabel(110, 575, "Opacity");
  muiSetActive(FiberViewOpacityRadio, opacityFlag);
  muiSetCallback(FiberViewOpacityRadio, fiberOpacityView);
*/

  muiNewBoldLabel(5, 150, "Contrast controls");
  sprintf(as,"Contrast\0");
  muiNewLabel(10, 120, as);
  ContrastSlider = muiNewHSlider(75, 115, 280, 180, 5);
  muiSetCallback(ContrastSlider, changeContrast);
  sprintf(as,"Brightness\0");
  muiNewLabel(10, 90, as);
  BrightnessSlider = muiNewHSlider(75, 85, 280, 180, 5);
  muiSetCallback(BrightnessSlider, changeBrightness);

  muiNewBoldLabel(5, 545, "Tracto Fiber Direction Flip");
  XFlipRadio = muiNewRadioButton( 10, 505 ) ;
  muiNewLabel(40, 510, "X Dir");
  YFlipRadio = muiNewRadioButton ( 80, 505 ) ;
  muiNewLabel(110, 510, "Y Dir");
  ZFlipRadio = muiNewRadioButton ( 150, 505 ) ;
  muiNewLabel(180, 510, "Z Dir");

  muiSetCallback(XFlipRadio, flipFiberXDir);
  muiSetCallback(YFlipRadio, flipFiberYDir);
  muiSetCallback(ZFlipRadio, flipFiberZDir);

  muiSetActive(XFlipRadio, xflipflag);
  muiSetActive(YFlipRadio, yflipflag);
  muiSetActive(ZFlipRadio, zflipflag);

  muiNewBoldLabel(5, 295, "Slice controls");

  //get slice location from slider event
  sprintf(as,"Axial\0");
  muiNewLabel(10, 265, as);
  Aslider = muiNewHSlider(75, 260, 280, 180, 5);
  muiSetCallback(Aslider, getAxialSliceLocation);

  sprintf(cs,"Coronal\0");
  muiNewLabel(10, 235, cs);
  Cslider = muiNewHSlider(75, 230, 280, 180, 5);
  muiSetCallback(Cslider, getCoronalSliceLocation);

  sprintf(ss,"Sagittal\0");
  muiNewLabel(10, 205, ss);
  Sslider = muiNewHSlider(75, 200, 280, 180, 5);
  muiSetCallback(Sslider, getSagittalSliceLocation);

  muiNewBoldLabel(5, 365, "View");

  SRadio = muiNewRadioButton( 10, 330 ) ;
  muiNewLabel(40, 335, "S");
  IRadio = muiNewRadioButton ( 60, 330 ) ;
  muiNewLabel(90, 335, "I");
  RRadio = muiNewRadioButton ( 110, 330 ) ;
  muiNewLabel(140, 335, "R");
  LRadio = muiNewRadioButton ( 160, 330 ) ;
  muiNewLabel(190, 335, "L");
  ARadio = muiNewRadioButton ( 210, 330 ) ;
  muiNewLabel(240, 335, "A");
  PRadio = muiNewRadioButton ( 260, 330) ;
  muiNewLabel(290, 335, "P");

  if (initview == "S") {
    muiSetActive(SRadio, 1);
    SFlag = 1;
  } else if (initview == "I") {
    muiSetActive(IRadio, 1);
    IFlag = 1;
  } else if (initview == "R") {
    muiSetActive(RRadio, 1);
    RFlag = 1;
  } else if (initview == "L") {
    muiSetActive(LRadio, 1);
    LFlag = 1;
  } else if (initview == "A") {
    muiSetActive(ARadio, 1);
    AFlag = 1;
  } else if (initview == "P") {
    muiSetActive(PRadio, 1);
    PFlag = 1;
  }

  muiSetCallback(SRadio, getSDirection);
  muiSetCallback(IRadio, getIDirection);
  muiSetCallback(RRadio, getRDirection);
  muiSetCallback(LRadio, getLDirection);
  muiSetCallback(ARadio, getADirection);
  muiSetCallback(PRadio, getPDirection);

  muiNewBoldLabel(5, 475, "Rotation");

  b3 = muiNewButton(10, 100, 435, 460);
  b4 = muiNewButton(10, 100, 400, 425);
  b1 = muiNewButton(110, 200, 435, 460);
  b2 = muiNewButton(110, 200, 400, 425);
  b5 = muiNewButton(210, 300, 435, 460);
  b6 = muiNewButton(210, 300, 400, 425);

  muiLoadButton(b1, "Y Rotation");
  muiLoadButton(b2, "Y Rotation(-)");
  muiLoadButton(b3, "X Rotation");
  muiLoadButton(b4, "X Rotation(-)");
  muiLoadButton(b5, "Z Rotation");
  muiLoadButton(b6, "Z Rotation(-)");

  muiSetCallback(b1, getYMovieMotion);
  muiSetCallback(b2, getYReverseMovieMotion);
  muiSetCallback(b3, getXMovieMotion);
  muiSetCallback(b4, getXReverseMovieMotion);
  muiSetCallback(b5, getZMovieMotion);
  muiSetCallback(b6, getZReverseMovieMotion);

  muiAddToUIList( 1, Aslider );
  muiAddToUIList( 1, Cslider );
  muiAddToUIList( 1, Sslider );
  muiAddToUIList( 1, SRadio );
  muiAddToUIList( 1, IRadio );
  muiAddToUIList( 1, RRadio );
  muiAddToUIList( 1, LRadio );
  muiAddToUIList( 1, ARadio );
  muiAddToUIList( 1, PRadio );
  muiAddToUIList( 1, XFlipRadio);
  muiAddToUIList( 1, YFlipRadio);
  muiAddToUIList( 1, ZFlipRadio);
}

void tractomain::AxialDisplay( )
{
  imageStruct imageinfo;
  imageinfo = this ->imageinfo;

  int i, j, ns, np;
  ns = imageinfo.height * imageinfo.width;
  np = ns * imageinfo.depth;
  int shift_xc = imageinfo.width/2;
  int shift_yc = imageinfo.height/2;
  int shift_zc = imageinfo.depth/2;

  glutSetWindow(2);

  glClear(GL_COLOR_BUFFER_BIT);
  glPushMatrix();
  glRotatef(-90.0, 0, 0, 1);
  for ( i = -shift_yc; i < shift_yc-1 ; i++) {
    for ( j = -shift_xc; j < shift_xc-1 ; j++) {
      int index = (imageinfo.AxialSliceLocation)*ns + (j+shift_yc)*imageinfo.width + shift_xc +i;
      glBegin(GL_QUADS);
      glColor4f(imageinfo.ImageArray[index],
        imageinfo.ImageArray[index],
        imageinfo.ImageArray[index], 1);
      glVertex3f(j, i, 0);
      glColor4f(imageinfo.ImageArray[index+imageinfo.height],
        imageinfo.ImageArray[index+imageinfo.height],
        imageinfo.ImageArray[index+imageinfo.height], 1);
      glVertex3f(j+1, i, 0);
      glColor4f(imageinfo.ImageArray[index+imageinfo.height+1],
        imageinfo.ImageArray[index+imageinfo.height+1],
        imageinfo.ImageArray[index+imageinfo.height+1], 1);
      glVertex3f(j+1, i+1, 0);
      glColor4f(imageinfo.ImageArray[index+1],
        imageinfo.ImageArray[index+1],
        imageinfo.ImageArray[index+1], 1);
      glVertex3f(j, i+1, 0);
      glEnd();
    }
  }
  glPopMatrix();
  glLineWidth (1);
  glBegin (GL_LINES);
  glColor4f(1.0, 0.0, 0.0, 1.0);
  glVertex2f (-shift_xc, AxialLineH);
  glVertex2f (shift_xc-1, AxialLineH);
  glVertex2f (AxialLineL, -shift_yc);
  glVertex2f (AxialLineL, shift_yc-1);
  glEnd ();
  glutSwapBuffers();

  this->imageinfo = imageinfo;
}
void tractomain::CoronalDisplay(void)
{
  imageStruct imageinfo;
  imageinfo = this ->imageinfo;

  int i, j, ns, np;
  ns = imageinfo.height * imageinfo.width;
  np = ns * imageinfo.depth;
  int shift_xc = imageinfo.width/2;
  int shift_yc = imageinfo.height/2;
  int shift_zc = imageinfo.depth/2;

  float ratioXZ = imageinfo.voxel_t/imageinfo.voxel_w;
  float halfRatioXZ = ratioXZ/2;

  glClear(GL_COLOR_BUFFER_BIT); glPushMatrix();
  glRotatef(-180.0, 0, 0, 1);
  for ( i = -shift_zc; i < shift_zc-1 ; i++) {
    for ( j = -shift_xc+1; j < shift_xc-1 ; j++) {
      int index = (shift_zc+i)*ns + (imageinfo.CoronalSliceLocation)*imageinfo.width + shift_xc+j-1;
      glBegin(GL_QUADS);
      glColor4f(imageinfo.ImageArray[index],
                imageinfo.ImageArray[index],
                imageinfo.ImageArray[index], 0.5);
      glVertex3f((float)j-0.5, (float)i*ratioXZ-halfRatioXZ, 0);
      glColor4f(imageinfo.ImageArray[index+1],
                imageinfo.ImageArray[index+1],
                imageinfo.ImageArray[index+1], 0.5);
      glVertex3f((float)j+0.5, (float)i*ratioXZ-halfRatioXZ, 0);
      glColor4f(imageinfo.ImageArray[index+ns+1],
                imageinfo.ImageArray[index+ns+1],
                imageinfo.ImageArray[index+ns+1], 0.5);
      glVertex3f((float)j+0.5, (float)i*ratioXZ+halfRatioXZ, 0);
      glColor4f(imageinfo.ImageArray[index+ns],
                imageinfo.ImageArray[index+ns],
                imageinfo.ImageArray[index+ns], 0.5);
      glVertex3f((float)j-0.5, (float)i*ratioXZ+halfRatioXZ, 0);
      glEnd();
    }
  }
  glPopMatrix();
  glLineWidth (1);
  glLineWidth (1);
  glBegin (GL_LINES);
  glColor4f(0.0, 1.0, 0.0, 1.0);
  glVertex2f (-shift_xc, CoronalLineH);
  glVertex2f (shift_xc-1, CoronalLineH);
  glVertex2f (CoronalLineL, (-shift_zc+1)*ratioXZ);
  glVertex2f (CoronalLineL, (shift_zc-1)*ratioXZ);
  glEnd ();
  glutSwapBuffers();
  this->imageinfo = imageinfo;
}

void tractomain::SagittalDisplay(void)
{
  imageStruct imageinfo;
  imageinfo = this ->imageinfo;

  int i, j, ns, np;
  ns = imageinfo.height * imageinfo.width;
  np = ns * imageinfo.depth;
  int shift_xc = imageinfo.width/2;
  int shift_yc = imageinfo.height/2;
  int shift_zc = imageinfo.depth/2;

  //imageinfo.SagittalSliceLocation = SsliceLoc;
  float ratioYZ = imageinfo.voxel_t/imageinfo.voxel_h;
  float halfRatioYZ = ratioYZ/2;
  glClear(GL_COLOR_BUFFER_BIT);
  glPushMatrix();
  //glRotatef(-180.0, 0, 0, 1);
  for ( i = -shift_zc; i < shift_zc-1 ; i++) {
    for ( j = -shift_yc; j < shift_yc-1 ; j++) { 
      int index = (shift_zc+i)*ns + (shift_yc+j)*imageinfo.width + imageinfo.SagittalSliceLocation;
      glBegin(GL_QUADS);
      glColor4f(imageinfo.ImageArray[index],
                imageinfo.ImageArray[index],
                imageinfo.ImageArray[index], 0.5);
      glVertex3f((float)(j-0.5), (float)-i*ratioYZ+halfRatioYZ, 0);
      glColor4f(imageinfo.ImageArray[index+imageinfo.width],
                imageinfo.ImageArray[index+imageinfo.width],
                imageinfo.ImageArray[index+imageinfo.width], 0.5);
      glVertex3f((float)(j+0.5), (float)-i*ratioYZ+halfRatioYZ, 0);
      glColor4f(imageinfo.ImageArray[index+ns+imageinfo.width],
                imageinfo.ImageArray[index+ns+imageinfo.width],
                imageinfo.ImageArray[index+ns+imageinfo.width], 0.5);
      glVertex3f((float)(j+0.5), (float)-i*ratioYZ-halfRatioYZ, 0);
      glColor4f(imageinfo.ImageArray[index+ns],
                imageinfo.ImageArray[index+ns],
                imageinfo.ImageArray[index+ns], 0.5);
      glVertex3f((float)(j-0.5), (float)-i*ratioYZ-halfRatioYZ, 0);
      glEnd();
    }
  }
  glPopMatrix();
  glLineWidth (1);
  glBegin (GL_LINES);
  glColor4f(0.0, 0.0, 1.0, 1.0);
  glVertex2f (-shift_yc, SagittalLineH);
  glVertex2f (shift_yc-1, SagittalLineH);
  glVertex2f (SagittalLineL, (-shift_zc+1)*ratioYZ);
  glVertex2f (SagittalLineL, (shift_zc-1)*ratioYZ);
  glEnd ();
  glutSwapBuffers();
  this->imageinfo = imageinfo;
}

void tractomain::AxialMouseClick(int button, int state, int xpos, int ypos)
{                       
  // set current wind 
  int winID = glutGetWindow();
  glutSetWindow(winID);

  AxialMoustPositionX = xpos;
  AxialMoustPositionY = ypos;

  AxialCurrentWidth = glutGet(GLUT_WINDOW_WIDTH);
  AxialCurrentHeight = glutGet(GLUT_WINDOW_HEIGHT);

  float h_ratio = ypos/AxialCurrentWidth;
  float l_ratio = xpos/AxialCurrentHeight;

  AxialLineH = -(h_ratio * imageinfo.width - imageinfo.width/2);
  AxialLineL = l_ratio * imageinfo.height - imageinfo.height/2;

  CoronalSliceLocation = h_ratio * imageinfo.width;
  SagittalSliceLocation = l_ratio * imageinfo.height;
  imageinfo.CoronalSliceLocation = (int)(h_ratio * imageinfo.width);
  imageinfo.SagittalSliceLocation = (int)(l_ratio * imageinfo.height);
  CoronalLineL =  imageinfo.SagittalSliceLocation - imageinfo.width/2;
  SagittalLineL =  imageinfo.CoronalSliceLocation - imageinfo.width/2;

  glutSetWindow(3);
  CoronalDisplay();
  glutSetWindow(4);
  SagittalDisplay();
  glutSetWindow(2);
  AxialDisplay();
}
void tractomain::CoronalMouseClick(int button, int state, int xpos, int ypos)
{   
  // set current wind 
  int winID = glutGetWindow();
  glutSetWindow(winID);

  CoronalCurrentWidth = glutGet(GLUT_WINDOW_WIDTH);
  CoronalCurrentHeight = glutGet(GLUT_WINDOW_HEIGHT);

  CoronalMoustPositionX = xpos;
  CoronalMoustPositionY = ypos;

  float h_ratio = ypos/CoronalCurrentWidth;
  float l_ratio = xpos/CoronalCurrentHeight;

  CoronalLineH = -(h_ratio * imageinfo.height - imageinfo.height/2);
  CoronalLineL = l_ratio * imageinfo.width -  imageinfo.width/2;

  float brainArea = (imageinfo.voxel_t*imageinfo.depth)/2+imageinfo.voxel_t;
  float emptyArea =  (imageinfo.height - (int)brainArea)/2;
  imageinfo.AxialSliceLocation = (int)((h_ratio*imageinfo.height-emptyArea)/imageinfo.voxel_t*2-0.5);
  imageinfo.SagittalSliceLocation =  (int)(l_ratio * imageinfo.width);
  AxialSliceLocation =  imageinfo.AxialSliceLocation;
  SagittalSliceLocation =  imageinfo.SagittalSliceLocation;

  AxialLineL = imageinfo.SagittalSliceLocation - imageinfo.width/2;
  SagittalLineH = CoronalLineH;

  glutSetWindow(4);
  SagittalDisplay();
  glutSetWindow(3);
  CoronalDisplay();
  glutSetWindow(2);
  AxialDisplay();
}


void tractomain::SagittalMouseClick(int button, int state, int xpos, int ypos)
{   
  // set current wind 
  int winID = glutGetWindow();
  glutSetWindow(winID);

  SagittalCurrentWidth = glutGet(GLUT_WINDOW_WIDTH);
  SagittalCurrentHeight = glutGet(GLUT_WINDOW_HEIGHT);

  SagittalMoustPositionX = xpos;
  SagittalMoustPositionY = ypos;

  float h_ratio = ypos/SagittalCurrentWidth;
  float l_ratio = xpos/SagittalCurrentHeight;

  SagittalLineH = -(h_ratio * imageinfo.height - imageinfo.height/2);
  SagittalLineL = (l_ratio * imageinfo.width -  imageinfo.width/2);

  float brainArea = (imageinfo.voxel_t*imageinfo.depth)/2+imageinfo.voxel_t;
  float emptyArea =  (imageinfo.height - (int)brainArea)/2;
  imageinfo.AxialSliceLocation = (int)((h_ratio*imageinfo.height-emptyArea)/imageinfo.voxel_t*2-0.5);
  imageinfo.CoronalSliceLocation =  (int)(l_ratio * imageinfo.width);
  AxialSliceLocation =  imageinfo.AxialSliceLocation;
  CoronalSliceLocation =  imageinfo.CoronalSliceLocation;

  AxialLineH = -(imageinfo.CoronalSliceLocation - imageinfo.height/2);
  CoronalLineH = SagittalLineH;

  glutSetWindow(2);
  AxialDisplay();
  glutSetWindow(3);
  CoronalDisplay();
  glutSetWindow(4);
  SagittalDisplay();
}

void getAxialSliceLocation( muiObject *obj, enum muiReturnValue rv )
{
  glutPostRedisplay();
  oldslideraxial = slideraxial;
  slideraxial = muiGetHSVal( obj);
}

void getCoronalSliceLocation( muiObject *obj, enum muiReturnValue rv )
{
  glutPostRedisplay();
  oldslidercoronal = slidercoronal;
  slidercoronal = muiGetHSVal( obj);
}

void getSagittalSliceLocation( muiObject *obj, enum muiReturnValue rv )
{
  glutPostRedisplay();
  oldslidersagittal = slidersagittal;
  slidersagittal = muiGetHSVal( obj);
  slidersagittal = 1.0 - slidersagittal;
}

void getXMovieMotion( muiObject *obj, enum muiReturnValue rv )
{
  if(movieXFlag == 0 ||  movieXFlag == -1) {
    movieXFlag = 1;    
    movieYFlag = 0;    
    movieZFlag = 0;    
  } else if(movieXFlag == 1) {
    delta = 0;
  }
  SFlag = 0; IFlag = 0; RFlag = 0;
  LFlag = 0; AFlag = 0; PFlag = 0;
//  cout << stepNum << endl;
}

void getYMovieMotion( muiObject *obj, enum muiReturnValue rv )
{
  if(movieYFlag == 0 ||  movieYFlag == -1) {
    movieYFlag = 1;    
    movieXFlag = 0;    
    movieZFlag = 0;    
  }
  else if(movieYFlag == 1) {
    delta = 0;
  }
  SFlag = 0; IFlag = 0; RFlag = 0;
  LFlag = 0; AFlag = 0; PFlag = 0;
//  cout << stepNum << endl;
}

void getZMovieMotion( muiObject *obj, enum muiReturnValue rv )
{
  if(movieZFlag == 0 ||  movieZFlag == -1) {
    movieZFlag = 1;    
    movieYFlag = 0;    
    movieXFlag = 0;    
  } else if(movieZFlag == 1) {
    delta = 0;
  }
  SFlag = 0; IFlag = 0; RFlag = 0;
  LFlag = 0; AFlag = 0; PFlag = 0;
//  cout << stepNum << endl;
}

void getXReverseMovieMotion( muiObject *obj, enum muiReturnValue rv )
{
  if(movieXFlag == 0 || movieXFlag == 1) {
    movieXFlag = -1;    
    movieYFlag = 0;    
    movieZFlag = 0;    
  } else if(movieXFlag == -1) {
    delta = 0;
  }
  SFlag = 0; IFlag = 0; RFlag = 0;
  LFlag = 0; AFlag = 0; PFlag = 0;
//  cout << stepNum << endl;
}

void getYReverseMovieMotion( muiObject *obj, enum muiReturnValue rv )
{
  if(movieYFlag == 0 || movieYFlag == 1) {
    movieXFlag = 0;    
    movieYFlag = -1;    
    movieZFlag = 0;    
  } else if(movieYFlag == -1) {
    delta = 0;
  }
  SFlag = 0; IFlag = 0; RFlag = 0;
  LFlag = 0; AFlag = 0; PFlag = 0;
//  cout << stepNum << endl;
}

void getZReverseMovieMotion( muiObject *obj, enum muiReturnValue rv )
{
  if(movieZFlag == 0 || movieZFlag == 1) {
    movieZFlag = -1;    
    movieXFlag = 0;    
    movieYFlag = 0;    
  } else if(movieZFlag == -1) {
    delta = 0;
  }
  SFlag = 0; IFlag = 0; RFlag = 0;
  LFlag = 0; AFlag = 0; PFlag = 0;
//  cout << stepNum << endl;
}

void getSDirection( muiObject *obj, enum muiReturnValue rv )
{
  int i;
//  cout << muiGetActive(obj) << endl;
  SFlag = muiGetActive(obj);
  muiSetActive(IRadio, 0);
  muiSetActive(RRadio, 0);
  muiSetActive(LRadio, 0);
  muiSetActive(ARadio, 0);
  muiSetActive(PRadio, 0);
  IFlag = 0;
  RFlag = 0;
  LFlag = 0;
  AFlag = 0;
  PFlag = 0;
  movieYFlag = 0;
  initview = "S";
}

void getIDirection( muiObject *obj, enum muiReturnValue rv )
{
  int i;
//  cout << muiGetActive(obj) << endl;
  IFlag = muiGetActive(obj);
  muiSetActive(SRadio, 0);
  muiSetActive(RRadio, 0);
  muiSetActive(LRadio, 0);
  muiSetActive(ARadio, 0);
  muiSetActive(PRadio, 0);
  SFlag = 0;
  RFlag = 0;
  LFlag = 0;
  AFlag = 0;
  PFlag = 0;
  movieXFlag = 0;
  movieYFlag = 0;
  movieZFlag = 0;
  initview = "I";
}

void getLDirection( muiObject *obj, enum muiReturnValue rv )
{
  int i;
//  cout << muiGetActive(obj) << endl;
  LFlag = muiGetActive(obj);
  muiSetActive(SRadio, 0);
  muiSetActive(IRadio, 0);
  muiSetActive(RRadio, 0);
  muiSetActive(ARadio, 0);
  muiSetActive(PRadio, 0);
  SFlag = 0;
  IFlag = 0;
  RFlag = 0;
  AFlag = 0;
  PFlag = 0;
  movieXFlag = 0;
  movieYFlag = 0;
  movieZFlag = 0;
  initview = "L";
}

void getRDirection( muiObject *obj, enum muiReturnValue rv )
{
  int i;
//  cout << muiGetActive(obj) << endl;
  RFlag = muiGetActive(obj);
  muiSetActive(SRadio, 0);
  muiSetActive(IRadio, 0);
  muiSetActive(LRadio, 0);
  muiSetActive(ARadio, 0);
  muiSetActive(PRadio, 0);
  SFlag = 0;
  IFlag = 0;
  LFlag = 0;
  AFlag = 0;
  PFlag = 0;
  movieXFlag = 0;
  movieYFlag = 0;
  movieZFlag = 0;
  initview = "R";
}

void getADirection( muiObject *obj, enum muiReturnValue rv )
{
  int i;
//  cout << muiGetActive(obj) << endl;
  AFlag = muiGetActive(obj);
  muiSetActive(SRadio, 0);
  muiSetActive(IRadio, 0);
  muiSetActive(LRadio, 0);
  muiSetActive(RRadio, 0);
  muiSetActive(PRadio, 0);
  SFlag = 0;
  IFlag = 0;
  RFlag = 0;
  LFlag = 0;
  PFlag = 0;
  movieXFlag = 0;
  movieYFlag = 0;
  movieZFlag = 0;
  initview = "A";
}

void getPDirection( muiObject *obj, enum muiReturnValue rv )
{
  int i;
//  cout << muiGetActive(obj) << endl;
  PFlag = muiGetActive(obj);
  muiSetActive(SRadio, 0);
  muiSetActive(IRadio, 0);
  muiSetActive(LRadio, 0);
  muiSetActive(RRadio, 0);
  muiSetActive(ARadio, 0);
  SFlag = 0;
  IFlag = 0;
  RFlag = 0;
  LFlag = 0;
  AFlag = 0;
  movieXFlag = 0;
  movieYFlag = 0;
  movieZFlag = 0;
  initview = "P";
}


int save_tiff (const char* fname)
{
  GLint rowlength, skiprows, skippixels, alignment;
  GLboolean swapbytes, lsbfirst;
  GLubyte* pixel_data = NULL;
  GLenum gl_error;
  TIFF *tiff = NULL;
  tsize_t line_bytes;
  unsigned char* line_buffer = NULL;
  int scan_line_size;
  int strip_size;
  int row;
  int height, width;

  width = glutGet(GLUT_WINDOW_WIDTH);
  height = glutGet(GLUT_WINDOW_HEIGHT);  

  /* Allocate a buffer for pixels. */
  pixel_data = (GLubyte*) malloc (width * height * 3);
  if (NULL == pixel_data) {
    printf("%s: ERROR: failed to allocate pixel storage\n",progname);
    exit(1);    
  }

  /* Read from the front buffer. */
  glReadBuffer (GL_FRONT);
  
  /* Save our unpack attributes. */
  glGetBooleanv (GL_PACK_SWAP_BYTES, &swapbytes);
  glGetBooleanv (GL_PACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv (GL_PACK_ROW_LENGTH, &rowlength);
  glGetIntegerv (GL_PACK_SKIP_ROWS, &skiprows);
  glGetIntegerv (GL_PACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv (GL_PACK_ALIGNMENT, &alignment);
  
  /* Set them. */
  glPixelStorei (GL_PACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei (GL_PACK_ROW_LENGTH, 0);
  glPixelStorei (GL_PACK_SKIP_ROWS, 0);
  glPixelStorei (GL_PACK_SKIP_PIXELS, 0);
  glPixelStorei (GL_PACK_ALIGNMENT, 1);
  
  /* Read RGB pixel data. */
  glReadPixels (0, 0, width, height, GL_RGB,
		GL_UNSIGNED_BYTE, (GLvoid*)pixel_data);

  /* Check error at this point. */
  gl_error = glGetError ();

  /* Restore the attributes. */
  glPixelStorei (GL_PACK_SWAP_BYTES, swapbytes);
  glPixelStorei (GL_PACK_LSB_FIRST, lsbfirst);
  glPixelStorei (GL_PACK_ROW_LENGTH, rowlength);
  glPixelStorei (GL_PACK_SKIP_ROWS, skiprows);
  glPixelStorei (GL_PACK_SKIP_PIXELS, skippixels);
  glPixelStorei (GL_PACK_ALIGNMENT, alignment);

  /* Handle error now. We can bail safely as we've restored the GL
     context. */
  if (GL_NO_ERROR != gl_error) {
    free (pixel_data);
    printf("%s: ERROR: failed to read pixels\n",progname);
    exit(1);    
  }

  /* Open a TIFF. */
  tiff = TIFFOpen( fname, "w" );
  if (NULL == tiff) {
    free (pixel_data);
    printf("%s: ERROR: failed to create file %s\n",progname,fname);
    exit(1);    
  }

  /* Set the TIFF info. */
  TIFFSetField (tiff, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField (tiff, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField (tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField (tiff, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField (tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField (tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField (tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  
  /* Calculate some sizes and allocate a line buffer. */
  line_bytes = 3 * width;
  line_buffer = NULL;
  scan_line_size = TIFFScanlineSize (tiff);
  if (scan_line_size != line_bytes) {
    fprintf(stderr,"%s: scan_line_size %d, line_bytes %d\n",
	    progname,scan_line_size, (int)line_bytes);
  }
  
  line_buffer = (unsigned char*) _TIFFmalloc( scan_line_size  );
  if (NULL == line_buffer) {
    free (pixel_data);
    TIFFClose (tiff);
    printf("%s: ERROR: failed to create tiff storage\n",progname);
    exit(1);    
  }
  
  /* Set the strip size to default. */
  strip_size = TIFFDefaultStripSize (tiff, width * 3);
  TIFFSetField (tiff, TIFFTAG_ROWSPERSTRIP, strip_size);
  
  /* Write line by line (bottom to top). */
  for (row = 0; row < height; row++) {
    memcpy (line_buffer, &pixel_data[(height-row-1) * line_bytes], 
	    line_bytes);
    TIFFWriteScanline (tiff, line_buffer, row, 0);
  }

  /* Close the tiff file and free the line buffer. */
  TIFFClose (tiff);
  _TIFFfree (line_buffer);
    
  free (pixel_data);

  return (0);
}

