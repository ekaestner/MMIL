
#include <iostream> 
#include <string> 
#include <fstream> 

using namespace std; 

#include "mri.h"
#include "mesh.h"
#include "statistical.h" 

int main(int argc, char** argv)
{

  if(argc < 2) {
    cout << "SYNOPSIS: *** fileName" << endl;
    abort(); 
  }

  enum {NO=0, YES, BOTH}; 
  enum {LOCAL=0, GLOBAL};  
  int surfaceDecimation, rasShift; 

  /* 4 directories for data input to be specified 

  1) mriInputString (mgh-file) - compulsory input from the screen 

  2) triInputDirectory (ico-files) 

  3) projectDirectory (mean, template, eof and talairach) 

  4) projectDirectory/methodSubdirectory (surfaces) 

  */

  string mriInputString=argv[1]; 

  string triInputDirectory="/home/igor/TRIANGLE/NEW_TRI_FILES/"; 

  string projectDirectory="/space/monkeys/1/home/igor/PROJECTS/VETSA_MGH15T_ATL/"; 

  string methodSubdirectory="DETERMINISTIC/"; 

  int mriInputLength=mriInputString.length();  
  int mriInputPivot=mriInputString.rfind("/", mriInputLength); 

  string mriInputDirectory = mriInputString.substr(0, mriInputPivot+1);  
  string mriFileName = mriInputString.substr(mriInputPivot+1, mriInputLength);
  string outputDirectory=mriInputDirectory;

  int mriFileTypeLength = mriFileName.length()-4; 
  string subjectName = 
    mriInputString.substr(mriInputPivot+1,mriFileTypeLength); 

  string surfaceInputDirectory= 
    projectDirectory + methodSubdirectory + "TRI/";

  string meanSurfaceInputDirectory=projectDirectory + "MEAN/";

  string templateInputDirectory=projectDirectory + "MEAN/";

  string eofInputDirectory=projectDirectory + "MEAN/";

  string talairachInputDirectory=projectDirectory + "XFM/";

  string talairachFileName = 
    talairachInputDirectory + subjectName + ".xfm"; 

  string surfaceType, surfaceName; 

  /*  number of surface movements, i.e. call to move Vertices */
  int kMove, nMove; 

  /* number of local smoothing operations at the end of segmentation */ 
  int kRest, nRest; 

  /* number of global smoothing operations after each move */    
  int kRelax, nRelax;

  /* number of global smoothing operations before and after each dilation */  
  int kAdjust, nAdjust; 

  /* number of dilations */ 
  int kDilate, nDilate; 

  int nIntegral, nTemplate; 

  float templateSize; 

  float mghIntensityScaling=1.;

  float coefTangential, coefNormal, coefMRI, coefRepelling, coefTemplate; 
  float directionMRI; 
  float intensityThreshould; 

  float stepDilate; 

  float coefSmoothLocal, coefSmoothGlobal, coefSmoothGlobalDilate;

  float move, rest; 

  int meshNumber=5; 
  int talairachFlag=YES; 

  int nTemplateOut=41; 

  statistical *pStatistical; 

  pStatistical = new statistical(mriFileName,
				 mriInputDirectory,  
				 talairachFlag, 
				 talairachFileName,
				 triInputDirectory, 
				 meshNumber);

  /* INNER SKULL */

  surfaceType = "_inner_skull"; 
  surfaceName = subjectName + surfaceType; 

  pStatistical->initializeSurface(templateSize=1.);  

  pStatistical->importSurface(surfaceName, 
			      surfaceInputDirectory);  

  pStatistical->exportTemplate(surfaceName, outputDirectory, nTemplateOut);


  /* OUTER SKULL */

  surfaceType = "_outer_skull"; 
  surfaceName = subjectName + surfaceType; 

  pStatistical->importSurface(surfaceName, 
			      surfaceInputDirectory);  

  pStatistical->exportTemplate(surfaceName, outputDirectory, nTemplateOut);


  /* OUTER SCALP */


  surfaceType = "_outer_scalp"; 
  surfaceName = subjectName + surfaceType; 

  pStatistical->importSurface(surfaceName, 
			      surfaceInputDirectory);  

  pStatistical->exportTemplate(surfaceName, outputDirectory, nTemplateOut);


  /* WARNING - free memory in statisitical class */ 

  //delete pStatistical;

  return 0;
}
