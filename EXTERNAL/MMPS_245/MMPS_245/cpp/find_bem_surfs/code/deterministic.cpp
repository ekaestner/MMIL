
#include <iostream> 
#include <sstream> 
#include <string> 
#include <fstream> 
#include <cmath> 

using namespace std;

#include "mri.h"
#include "mesh.h"
#include "deterministic.h"

/* deterministic constructor */

deterministic::deterministic(string mriFileName, 
			     string mriInputDirectory, 
			     int talairachFlag, 
			     string talairachFileName, 
			     string triInputDirectory_ini, 
			     int meshNumber_ini): 
  mesh(mriFileName, 
       mriInputDirectory, 
       talairachFlag, 
       talairachFileName, 
       triInputDirectory_ini, 
       meshNumber_ini) {

  /* add initialization of MRI force computation method here */ 

}

/* deterministic destructor */

deterministic::~deterministic(){

  // to be added 

}


/* compute and add MRI force if repelling force is zero */ 

void deterministic::addMriForce(float coefMRI,
				float directionMRI,  
				float intensityThreshold) {

  int i,j; 

  const int nCost=10; 
  const float stepCost=0.5;

  float alpha=log(sqrt(3.)/stepCost); 

  float cost, factor1, factor2;  
  float forceMRI[3]; 
  
  float dummy, intensity; 

  for(j=0; j < nVertex; j++) {

    dummy=0.; 
    for(i=0; i < 3; i++) 
      dummy += vertex[j].repellingForce[i]*vertex[j].repellingForce[i]; 

    if(!dummy) {


      /* Add this functionality later  

      cost=1.; 

      for(int l=0; l < nCost; l++) {
	for(int i=0; i < 3; i++) 
	  coor[i] = vertex[j].r[i]-(l+1)*stepCost*vertex[j].fNormal[i]; 

	factor1=getIntensity(coor[0], coor[1], coor[2])-
	  intensityThreshold+stepCost; 
        factor2= (factor1 < 0.0) ? 0. : tanh(alpha*factor1); 
        cost *= factor2; 
      }
      
      WARNING - implement local approximation below  */ 

      intensity=getIntensity(vertex[j].r[0],vertex[j].r[1],vertex[j].r[2]); 
      cost=tanh(directionMRI*(intensity-intensityThreshold)); 

      for(int i=0; i <3; i++) 
	vertex[j].force[i] += coefMRI*cost*vertex[j].fNormal[i];  

    }

  }
}


/* segment surface using deformable model and MRI force */ 

void deterministic::segmentSurface(string surfaceName, 
                                   string outputDirectory, 
				   int nMove, 
                                   int nRest, 
                                   int nRelax, 
                                   float coefTangential, 
				   float coefNormal, 
				   float coefMRI,
				   float directionMRI, 
				   float coefRepelling, 
				   float intensityThreshold) {


  float move, rest; 

  for(int kMove=0; kMove < nMove; kMove++) {

    move = moveVertices(coefTangential, 
			coefNormal, 
			coefMRI, 
			directionMRI, 
			coefRepelling, 
			intensityThreshold); 

//    cout << "iteration # " << (kMove+1) << " " << "(" << move << ")" << endl;
    cout << "move # " << (kMove+1) << " " << "(" << move << ")" << endl;

    if(nRelax) for(int kRelax=0; kRelax < nRelax; kRelax++) 
      rest = smoothSurface(GLOBAL); 
  }

  for(int kRest=0; kRest < nRest; kRest++) {

    rest = smoothSurface(LOCAL);  
//    cout << "iteration # " << (kRest+1) << " " << "(" << rest << ")" << endl;  
    cout << "rest # " << (kRest+1) << " " << "(" << rest << ")" << endl;  
  }

  string parameterOutputFileName = outputDirectory +  
    surfaceName + ".par";

  ofstream parameterOutputFile(parameterOutputFileName.c_str());

  parameterOutputFile << "nMove  " << nMove << endl; 
  parameterOutputFile << "nRest  " << nRest << endl; 
  parameterOutputFile << "nRelax  " << nRelax << endl; 

  parameterOutputFile << "coefTangential  " << coefTangential << endl; 
  parameterOutputFile << "coefNormal  " << coefNormal << endl;
  parameterOutputFile << "coefMRI  " << coefMRI << endl;

  parameterOutputFile << "coefRepelling  " << coefRepelling << endl;

  parameterOutputFile << "intensityThreshold  " << intensityThreshold 
		      << endl;

  parameterOutputFile.close();

}


/* move all vertices by one step following an algorithm by 
 * 
 *  Dale, AM, B. Fischl, and M.I. Sereno, Cortical surface-based 
 *  analysis  I: Segmentation and surface reconstruction, 
 *  Neuroimage, 9, 179-194, 1999. 
 *                                       
 *  and compute average displacement            */

float deterministic::moveVertices(float coefTangential, 
				  float coefNormal, 
				  float coefMRI,
				  float directionMRI, 
				  float coefRepelling, 
				  float intensityThreshold) {

  int i, j; 

  float dummy; 
  float globalDisplacement, vertexDisplacement; 

  localCoordinateSystem();

  computeRepellingForce(coefRepelling); 

  addInertialForce(coefTangential, coefNormal); 

  addMriForce(coefMRI, directionMRI, intensityThreshold);

  globalDisplacement=0.; 

  for(int j=0; j < nVertex; j++) { 
    
    vertexDisplacement=0.; 

    for(int i=0; i < 3; i++) { 
      dummy = vertex[j].force[i];
      vertex[j].r[i] += dummy; 
      vertexDisplacement += dummy*dummy; 
    }

    globalDisplacement += sqrt(vertexDisplacement); 

  }

  return globalDisplacement/nVertex; 

}
