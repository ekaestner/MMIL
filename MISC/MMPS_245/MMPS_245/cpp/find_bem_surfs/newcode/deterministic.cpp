
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

deterministic::deterministic(const string &mriInputFileName, 
			     int talairachFlag, 
			     string talairachFileName, 
			     string triInputDirectory): 
  mesh(mriInputFileName, 
       talairachFlag, 
       talairachFileName, 
       triInputDirectory) {

  /* add initialization of MRI force computation method here */ 

  surface = new structModelCoefficient; 

}

/* deterministic destructor */

deterministic::~deterministic(){

  delete surface; 

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
                                   string outputDirectory) { 

  float move, rest; 

  cout << "Reading model coefficients for " << surfaceName << endl; 

  cout << "nMove" << endl; cin >> surface->nMove; 
  cout << "nRest" << endl; cin >> surface->nRest; 
  cout << "nRelax" << endl; cin >> surface->nRelax; 

  cout << "coefTangential" << endl; cin >> surface->coefTangential;  
  cout << "coefNormal" << endl; cin >> surface->coefNormal;
  cout << "coefMRI" << endl; cin >> surface->coefMRI;

  cout << "directionMRI" << endl; cin >> surface->directionMRI;
  cout << "coefRepelling" << endl; cin >> surface->coefRepelling;
  cout << "intensityThreshould" << endl; cin >> surface->intensityThreshould;


  for(int kMove=0; kMove < surface->nMove; kMove++) {

    move = moveVertices(surface->coefTangential, 
			surface->coefNormal, 
			surface->coefMRI, 
			surface->directionMRI, 
			surface->coefRepelling, 
			surface->intensityThreshould); 

    cout << surfaceName << " "; 
    cout << "step # " << (kMove+1) << " " << "(" << move << ")" << endl;

    if(surface->nRelax) for(int kRelax=0; kRelax < surface->nRelax; kRelax++) 
      rest = smoothSurface(GLOBAL); 
  }

  for(int kRest=0; kRest < surface->nRest; kRest++) {

    rest = smoothSurface(LOCAL);

    cout << surfaceName << " (smoothing) ";   
    cout << "step # " << (kRest+1) << " " << "(" << rest << ")" << endl;  
  }

  /* WARNING - no parameter output 

  string parameterOutputFileName = outputDirectory +  
    surfaceName + ".par";

  ofstream parameterOutputFile(parameterOutputFileName.c_str());

  parameterOutputFile << "nMove  " << surface->nMove << endl; 
  parameterOutputFile << "nRest  " << surface->nRest << endl; 
  parameterOutputFile << "nRelax  " << surface->nRelax << endl; 

  parameterOutputFile << "coefTangential  " << surface->coefTangential 
		      << endl; 
  parameterOutputFile << "coefNormal  " << surface->coefNormal << endl;
  parameterOutputFile << "coefMRI  " << surface->coefMRI << endl;

  parameterOutputFile << "coefRepelling  " << surface->coefRepelling << endl;

  parameterOutputFile << "intensityThreshould  " 
		      << surface->intensityThreshould  << endl;

  parameterOutputFile.close();

  */

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
				  float intensityThreshould) {

  int i, j; 

  float dummy; 
  float globalDisplacement, vertexDisplacement; 

  localCoordinateSystem();

  computeRepellingForce(coefRepelling); 

  addInertialForce(coefTangential, coefNormal); 

  addMriForce(coefMRI, directionMRI, intensityThreshould);

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

/* read model coefficients - currently not used */

void deterministic::readModelCoefficients() {

  surface->nMove=int(getNumber()); 
  surface->nRest=int(getNumber()); 
  surface->nRelax=int(getNumber()); 
  surface->coefTangential=getNumber(); 
  surface->coefNormal=getNumber(); 
  surface->coefMRI=getNumber();
  surface->directionMRI=getNumber();
  surface->coefRepelling=getNumber();
  surface->intensityThreshould=getNumber();

} 

/* read one coeffcient - currently not used */  


float deterministic::getNumber() {

  string inputLine; 
  float inputNumber; 

  getline(inputFile,inputLine); 
  istringstream inputStream(inputLine.substr(0,inputLine.find(" ")));
  inputStream >> inputNumber;

  // cout << inputNumber << endl; 
 
  return inputNumber; 
}
