
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


#include <iostream> 
#include <string> 
#include <fstream> 
#include <sstream> 

using namespace std; 

#include "mri.h"
#include "directory.h"
#include "mesh.h"
#include "deterministic.h"
#include "isotrack.h" 

int main()
{

  enum {NO=0, YES, BOTH}; 

  enum {innerSkull, outerSkull, outerScalp};  

  int surfaceDecimation, rasShift; 

  directory* pDirectory = new directory; 

  string talairachFileName = ".xfm"; // REMOVE 

   /* size of spherical template */ 
  float templateSize; 

  int nSmoothness; 

  /* the overall magnitude of surface dilation */  
  float movingDistance; 

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
 
  /* initialize algorithm with the Mesh Number > 0 */  

  /* currently RAS coordinates are adjusted in getIntensity subroutine, 
   * set rasShift to NO everywhere */ 

  int talairachFlag=YES; 

  deterministic* pMesh = new deterministic(pDirectory->mriInputString, 
					   talairachFlag, 
					   talairachFileName,
					   pDirectory->getTriInputDirectory());

#ifdef INNER_SKULL /* beginning of the inner skull processing */ 

  pMesh->initializeSurface(templateSize=30.);  

  pMesh->segmentSurface(pDirectory->getSurfaceName(innerSkull), 
			pDirectory->getMriInputDirectory());  

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       pDirectory->getSurfaceName(innerSkull), 
		       pDirectory->getMriInputDirectory());


#endif /* INNER_SKULL - end of the inner skull processing */


/******************************************************************/


#ifdef OUTER_SKULL /* beginning of the outer skull processing */ 

     
  pMesh->createRepellingVolume(movingDistance=3., // 3   
			       nSmoothness=1);  // 1  
  

  pMesh->segmentSurface(pDirectory->getSurfaceName(outerSkull), 
			pDirectory->getMriInputDirectory()); 

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       pDirectory->getSurfaceName(outerSkull), 
		       pDirectory->getMriInputDirectory());


#endif /* OUTER_SKULL - end of the outer skull processing */

  /*******************************************************************/


#ifdef OUTER_SCALP /* beginning of the outer scalp processing */ 


  pMesh->createRepellingVolume(movingDistance=3., // 3   
			       nSmoothness=1);  // 1  
 

  pMesh->initializeSurface(templateSize=150.);  

  pMesh->segmentSurface(pDirectory->getSurfaceName(outerScalp), 
			pDirectory->getMriInputDirectory()); 

  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       pDirectory->getSurfaceName(outerScalp), 
		       pDirectory->getMriInputDirectory());


#endif /* OUTER_SCALP - end of the outer scalp processing */

  delete pMesh;

  return 0;
}
