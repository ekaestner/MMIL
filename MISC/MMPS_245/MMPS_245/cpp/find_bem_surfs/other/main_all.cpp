
#ifndef INNER_SKULL
#define INNER_SKULL
//#undef  INNER_SKULL 
#endif

#ifndef OUTER_SKULL
#define OUTER_SKULL
//#undef  OUTER_SKULL
#endif

#ifndef OUTER_SCALP
#define OUTER_SCALP
//#undef  OUTER_SCALP
#endif

#ifndef DETERMINISTIC
#define DETERMINISTIC
//#undef  DETERMINISTIC
#endif 

#ifndef TEMPLATE
#define TEMPLATE
#undef  TEMPLATE
#endif   

#ifndef VOLUME
#define VOLUME
#undef  VOLUME
#endif   

#ifndef STATISTICAL
#define STATISTICAL
#undef  STATISTICAL
#endif  

#ifndef SPECTRAL
#define SPECTRAL
#undef  SPECTRAL
#endif  

#include <iostream> 
#include <string> 
#include <fstream> 

using namespace std; 

#include "mri.h"
#include "mesh.h"
#include "deterministic.h"
#include "isotrack.h" 

void saveDeterministicParameters(string outputDirectory, 
				 string surfaceName, 
                                 float templateSize, 
				 int nMove, 
				 int nRest, 
				 int nRelax, 
				 int nAdjust, 
				 int nDilate,  
				 float stepDilate, 
				 float coefSmoothGlobalDilate, 
				 float coefTangential, 
				 float coefNormal, 
				 float coefMRI, 
				 float coefRepelling,  
				 float directionMRI, 
				 float intensityThreshould,  
				 float coefSmoothLocal, 
				 float coefSmoothGlobal)

{
  string parameterOutputFileName = outputDirectory +  
    surfaceName + ".par";

  ofstream parameterOutputFile(parameterOutputFileName.c_str());

  parameterOutputFile << "templateSize  " << templateSize << endl; 

  parameterOutputFile << "nMove  " << nMove << endl; 
  parameterOutputFile << "nRest  " << nRest << endl; 
  parameterOutputFile << "nRelax  " << nRelax << endl; 
  parameterOutputFile << "nAdjust  " << nAdjust << endl; 
  parameterOutputFile << "nDilate  " << nDilate << endl; 

  parameterOutputFile << "stepDilate  " << stepDilate << endl; 
  parameterOutputFile << "coefSmoothGlobalDilate  " 
		      << coefSmoothGlobalDilate  << endl;

  parameterOutputFile << "coefTangential  " << coefTangential << endl; 
  parameterOutputFile << "coefNormal  " << coefNormal << endl;
  parameterOutputFile << "coefMRI  " << coefMRI << endl;
  parameterOutputFile << "directionMRI  " << directionMRI << endl;

  parameterOutputFile << "intensityThreshould  " << intensityThreshould 
		      << endl;

  parameterOutputFile << "coefSmoothGlobal  " << coefSmoothGlobal << endl;
  parameterOutputFile << "coefSmoothLocal  " << coefSmoothLocal << endl;

  parameterOutputFile.close();
}


int main(int argc, char** argv)
{

  if(argc < 2) {
    cout << "SYNOPSIS: *** fileName" << endl;
    abort(); 
  }

  enum {NO=0, YES, BOTH}; 

  int surfaceDecimation, rasShift; 

  /* 4 directories for data input to be specified 

  1) mriInputString (mgh-file) - compulsory input from the screen 

  2) triInputDirectory (ico-files) 

  3) projectDirectory (mean, template, eof and talairach) 

  4) projectDirectory/methodSubdirectory (surfaces) 

  */

  string mriInputString=argv[1]; 

  string triInputDirectory="/home/igor/TRIANGLE/NEW_TRI_FILES/"; 

  string projectDirectory="/space/monkeys/1/home/PROJECTS/VETSA_UCSD15T/"; 

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

  float coefTangential, coefNormal, coefMRI, coefRepelling, coefTemplate; 
  float directionMRI; 
  float intensityThreshould; 

  float stepDilate; 


  /* WARNING - remove later */ 

  enum {LOCAL=0, GLOBAL}; 
  float coefSmoothLocal, coefSmoothGlobal, coefSmoothGlobalDilate;

  float move, rest; 

  //mesh *pMesh;

  deterministic *pMesh; 

  /* Mesh Number    File       Vertices     Faces  
   *  
   *     0         ic0.tri        12          20   
   *     1         ic1.tri        42          80  
   *     2         ic2.tri       162         320  
   *     3         ic3.tri       642        1280  
   *     4         ic4.tri      2562        5120  
   *     5         ic5.tri     10242       20480  
   *     6         ic6.tri     40962       81920  
   *     7         ic7.tri    163842      327680  
   *                                                    */ 
 
  /* initialize algorithm with the Mesh Number  
   *
   *  which should be > 0 to allow decimation  
   * 
   *  currently tri-files do not support decimation for 
   *  Mesh Number > 4 
   *                                                 */ 

  /* currently RAS coordinates are adjusted in getIntensity subroutine, 
   * set rasShift to NO everywhere */ 

  int meshNumber=5; 
  int talairachFlag=YES; 

  //  pMesh = new mesh(mriFileName,


  pMesh = new deterministic(mriFileName,
			    mriInputDirectory,  
			    talairachFlag, 
			    talairachFileName,
			    triInputDirectory, 
			    meshNumber);

  /* define common parameters */ 

  stepDilate=1.;

#ifdef INNER_SKULL /* beginning of the inner skull processing */ 

  /* set parameters for the inner skull segmentation */ 

 /* Optimal coefficients for the inner skull 
  * using "flash" files, i.e. without 1.e08 scaling
  * see old versions for the scaling case  
  *  
  *
  * tri  nStep  Templ. Tang. Normal   Smooth   MRI  Repel.  Threshould 
  *
  * ic5   300    30    0.2    0.15    0./0.1   0.3   0.        350. 
  *                                                                    */ 

  surfaceType = "_inner_skull"; 
  surfaceName = subjectName + surfaceType; 

  pMesh->initializeSurface(templateSize=30.);  

#ifdef DETERMINISTIC /* beginning of the inner skull processing 
			using deterministic deformable model */ 

  nAdjust=0; 
  nDilate=0; 

  pMesh->segmentSurface(nMove=200, // 300
                        nRest=30, // 30 
                        nRelax=1, // 1
			coefTangential=0.2, // 0.2
			coefNormal=0.15, // 0.15
			coefMRI=0.3, // 0.3
			directionMRI=1, 
			coefRepelling=0., 
			intensityThreshould=275.); // 325.

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       surfaceName, 
		       outputDirectory);

  saveDeterministicParameters(outputDirectory, surfaceName, 
                              templateSize, 
			      nMove, nRest, nRelax, nAdjust, 
			      nDilate, stepDilate, coefSmoothGlobalDilate, 
			      coefTangential, coefNormal, coefMRI, 
			      coefRepelling, directionMRI, 
			      intensityThreshould,  
			      coefSmoothLocal, coefSmoothGlobal); 


#endif /* DETERMINISTIC - end of the inner skull processing 
	  using deterministic deformable model */ 


#ifdef TEMPLATE /* beginning of computing templates for the pre-segmented 
		   inner skull surface */ 

  pMesh->importSurface(surfaceName, 
		       surfaceInputDirectory);  

  pMesh->exportTemplate(surfaceName, outputDirectory);


#endif /* TEMPLATE - end of computing templates for the pre-segmented 
	  inner skull surface */ 


#ifdef VOLUME /* beginning of computing volumes for the pre-segmented 
		   inner skull surface */ 

  pMesh->importSurface(surfaceName, 
		       surfaceInputDirectory);

  //pstxDeformableModel->inflateSurface(0.75);   

  //pstxDeformableModel->talairachForward(); 

  float volume = pMesh->fillContour()/1.e03; 

  cout <<  subjectName << "    " << volume << endl; 

  //pstxDeformableModel->exportFilling(surfaceName, outputDirectory);


#endif /* VOLUME - end of computing volumes for the pre-segmented 
	  inner skull surface */ 


#ifdef STATISTICAL /* beginning of the inner skull processing 
			using statistical deformable model */ 


  surfaceName="mean" + surfaceType; 

  pMesh->importSurface(surfaceName, 
		       meanSurfaceInputDirectory);
 
  pMesh->talairachBackward(); 
 
  pMesh->importMeanTemplate(surfaceType, 
			    templateInputDirectory); 

  for(kDilate=0; kDilate < 100; kDilate++) { //40 - 100

    pMesh->adjustSurface(coefTangential=0.2, // 0.05 - 0.2
			 coefNormal=0.1, // 0.05 - 0.1
			 coefTemplate=0.05, // 0.2 - 0.05
			 coefRepelling=0.,
			 nIntegral=41, // 21
			 nTemplate=41); // 15


    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);
    
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);

  }

  for(kDilate=0; kDilate < 40; kDilate++) { //40 - 100

    pMesh->adjustSurface(coefTangential=0.15, // 0.05 - 0.2
			 coefNormal=0.1, // 0.05 - 0.1
			 coefTemplate=0.15, // 0.2 - 0.05
			 coefRepelling=0.,
			 nIntegral=15, // 21
			 nTemplate=11); // 15

    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);

  }
  
  surfaceName=subjectName + surfaceType; 

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       surfaceName, 
		       outputDirectory);


#endif /* STATISTICAL - end of the inner skull processing 
	  using deterministic deformable model */ 



#endif /* INNER_SKULL - end of the inner skull processing */


  /******************************************************************/



#ifdef OUTER_SKULL /* beginning of the outer skull processing */ 

  /* set parameters for the outer skull segmentation 
  *
  * tri  nStep  Templ. Tang. Normal   Smooth   MRI  Repel.  Threshould 
  *
  * ic5   100    --    0.15    0.1    0./0.1   0.3   0.2       250. 
  *                                                                    */ 

  surfaceType="_outer_skull"; 
  surfaceName=subjectName + surfaceType; 

#ifdef DETERMINISTIC /* beginning of the outer skull processing 
			using deterministic deformable model */ 

  templateSize=0.; 
  nMove=100; 
  nRest=30; 
  nRelax=2; // 1
  nAdjust=1; // 1
  nDilate=3; // 4

  for(kDilate=0; kDilate < nDilate; kDilate++) { 
 
  for(kAdjust=0; kAdjust < nAdjust; kAdjust++) 
    rest = pMesh->smoothSurface(GLOBAL); 
     
    pMesh->dilateSurface(stepDilate); 

  for(kAdjust=0; kAdjust < nAdjust; kAdjust++) 
    rest = pMesh->smoothSurface(GLOBAL);  

  }

  pMesh->fillContour(); 

  pMesh->segmentSurface(nMove=100, // 100
                        nRest=30, // 30 
                        nRelax=2, // 1
			coefTangential=0.15, // 0.15
			coefNormal=0.1, // 0.1
			coefMRI=0.2, // 0.2
			directionMRI=-1, 
			coefRepelling=0.5, // 0.2 
			intensityThreshould=250.); // 250.

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       surfaceName, 
		       outputDirectory);

  saveDeterministicParameters(outputDirectory, surfaceName,
                              templateSize,  
			      nMove, nRest, nRelax, nAdjust, 
			      nDilate, stepDilate, coefSmoothGlobalDilate, 
			      coefTangential, coefNormal, coefMRI, 
			      coefRepelling, directionMRI, 
			      intensityThreshould,  
			      coefSmoothLocal, coefSmoothGlobal); 


#endif /* DETERMINISTIC - end of the outer skull processing 
	  using deterministic deformable model */ 


#ifdef TEMPLATE /* beginning of computing templates for the pre-segmented 
		   outer skull surface */ 

  pMesh->initializeSurface();  

  pMesh->importSurface(surfaceName, 
		       surfaceInputDirectory);  

  pMesh->exportTemplate(surfaceName, outputDirectory);

#endif /* TEMPLATE - end of computing templates for the pre-segmented 
	  outer skull surface */ 


#ifdef VOLUME /* beginning of computing volumes for the pre-segmented 
		   outer skull surface */ 

  pMesh->initializeSurface();  

  pMesh->importSurface(surfaceName, 
		       surfaceInputDirectory);

  pMesh->inflateSurface(0.75);   

  pMesh->talairachForward(); 

  pMesh->fillContour(); 

  pMesh->exportFilling(surfaceName, outputDirectory);

#endif /* VOLUME - end of computing volumes for the pre-segmented 
	  outer skull surface */ 


#ifdef STATISTICAL /* beginning of the outer skull processing 
			using statistical deformable model */ 


  nAdjust=1; // 1
  nDilate=4; // 4

  for(kDilate=0; kDilate < nDilate; kDilate++) { 
 
  for(kAdjust=0; kAdjust < nAdjust; kAdjust++) 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobalDilate); 
     
    pMesh->dilateSurface(stepDilate); 

  for(kAdjust=0; kAdjust < nAdjust; kAdjust++) 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobalDilate);  

  }

  coefSmoothGlobalDilate=coefSmoothGlobal; 

  pMesh->fillContour(); 

  /*
  surfaceName="mean" + surfaceType; 

  pstxDeformableModel->importSurface(surfaceName, 
				     meanSurfaceInputDirectory);
 
  pstxDeformableModel->talairachBackward(); 
  */
 
  pMesh->importMeanTemplate(surfaceType, 
			    templateInputDirectory); 

  for(kDilate=0; kDilate < 40; kDilate++) { // 40

    pMesh->adjustSurface(coefTangential=0.15, // 0.15
			 coefNormal=0.1, // 0.1
			 coefTemplate=0.1, // 0.1 
			 coefRepelling=0.8, // 0.8
			 nIntegral=21, // 21
			 nTemplate=15); // 15
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);
  }
  
  surfaceName=subjectName + surfaceType; 

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       surfaceName, 
		       outputDirectory);


#endif /* STATISTICAL - end of the outer skull processing 
	  using deterministic deformable model */ 


#endif /* OUTER_SKULL - end of the outer skull processing */

  /*******************************************************************/


#ifdef OUTER_SCALP /* beginning of the outer scalp processing */ 

  /* set parameters for the outer scalp segmentation */ 

 /* Optimal coefficients for the inner skull 
  * using *me* files, i.e. with 1.e08 scaling 
  *  
  *
  * tri  nStep  Templ. Tang. Normal   Smooth   MRI  Repel.  Threshould 
  *
  * ic5   200    150   0.03   0.02    0./0.1   0.7    0.       125. 
  *                                                                    */ 

  surfaceType="_outer_scalp"; 
  surfaceName=subjectName + surfaceType; 

  pMesh->initializeSurface(templateSize=150.); 

#ifdef DETERMINISTIC /* beginning of the outer scalp processing 
			using deterministic deformable model */ 


  nMove=300; 
  nRest=30; 
  nRelax=2; // 1  
  nAdjust=0; 
  nDilate=0; 

  pMesh->segmentSurface(nMove=300, // 300
                        nRest=30, // 30 
                        nRelax=1, // 1
			coefTangential=0.04, // 0.04
			coefNormal=0.02, 
			coefMRI=0.3,
			directionMRI=1., 
			coefRepelling=0., 
			intensityThreshould=125.);// 125.

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       surfaceName, 
		       outputDirectory);

  saveDeterministicParameters(outputDirectory, surfaceName, 
                              templateSize, 
			      nMove, nRest, nRelax, nAdjust, 
			      nDilate, stepDilate, coefSmoothGlobalDilate, 
			      coefTangential, coefNormal, coefMRI, 
			      coefRepelling, directionMRI, 
			      intensityThreshould,  
			      coefSmoothLocal, coefSmoothGlobal); 


#endif /* DETERMINISTIC - end of the outer scalp processing 
	  using deterministic deformable model */ 


#ifdef TEMPLATE /* beginning of computing templates for the pre-segmented 
		   outer scalp surface */ 

  pMesh->importSurface(surfaceName, 
		       surfaceInputDirectory);  

  pMesh->exportTemplate(surfaceName, outputDirectory);

#endif /* TEMPLATE - end of computing templates for the pre-segmented 
	  outer scalp surface */ 


#ifdef VOLUME /* beginning of computing volumes for the pre-segmented 
		   outer scalp surface */ 

  pMesh->initializeSurface();  

  pMesh->importSurface(surfaceName, 
		       surfaceInputDirectory);

  pMesh->inflateSurface(0.75);   

  pMesh->talairachForward(); 

  pMesh->fillContour(); 

  pMesh->exportFilling(surfaceName, outputDirectory);

#endif /* VOLUME - end of computing volumes for the pre-segmented 
	  outer scalp surface */ 


#ifdef STATISTICAL /* beginning of the outer scalp processing 
			using statistical deformable model */ 


  nAdjust=2; // 1
  nDilate=4; // 4

  for(kDilate=0; kDilate < nDilate; kDilate++) { 
 
  for(kAdjust=0; kAdjust < nAdjust; kAdjust++) 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobalDilate); 
     
    pMesh->dilateSurface(stepDilate); 

  for(kAdjust=0; kAdjust < nAdjust; kAdjust++) 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobalDilate);  

  }

  pMesh->fillContour(); 


  surfaceName="mean" + surfaceType; 

  pMesh->importSurface(surfaceName, 
		       meanSurfaceInputDirectory);
 
  pMesh->talairachBackward(); 

 
  pMesh->importMeanTemplate(surfaceType, 
			    templateInputDirectory); 

  for(kDilate=0; kDilate < 100; kDilate++) { // 40 - 100

    pMesh->adjustSurface(coefTangential=0.15, // 0.15
			 coefNormal=0.1, // 0.1
			 coefTemplate=0.3, // 0.3
			 coefRepelling=1., // 1.
			 nIntegral=31, // 21 - 31
			 nTemplate=31); // 15 - 31

    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 

    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 

  }


  for(kDilate=0; kDilate < 40; kDilate++) { // 40 - 100

    pMesh->adjustSurface(coefTangential=0.15, // 0.15
			 coefNormal=0.1, // 0.1
			 coefTemplate=0.3, // 0.3
			 coefRepelling=1., // 1.
			 nIntegral=21, // 21 - 31
			 nTemplate=15); // 15 - 31

    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal); 
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);
    rest = pMesh->smoothSurface(GLOBAL,coefSmoothGlobal);

  }

  
  surfaceName=subjectName + surfaceType; 

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       surfaceName, 
		       outputDirectory);


#endif /* STATISTICAL - end of the outer scalp processing 
	  using deterministic deformable model */ 


#endif /* OUTER_SCALP - end of the outer scalp processing */


  delete pMesh;

  return 0;
}
