
#include <iostream> 
#include <sstream> 
#include <string> 
#include <fstream> 
#include <cmath> 

using namespace std;

#include "mri.h"
#include "mesh.h"
#include "isotrack.h"

/* isotrack constructor */

isotrack::isotrack(string mriFileName, 
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

  int i, j, k; 
  int nDummy; 

  nParameter=6; 

  /* WARNING */


  angularScale=0.5; 
  translationalScale=20.; 

  transformationMatrix = new float* [nRas];
  transformationMatrix0 = new float* [nRas];
  for (int i=0; i < nRas; i++) {
    transformationMatrix[i] = new float [nRas];
    transformationMatrix0[i] = new float [nRas];
  }

  for(j=0; j < nRas; j++) {
    for(i=0; i < nRas; i++) transformationMatrix[j][i]=0.;
    for(i=0; i < nRas; i++) transformationMatrix0[j][i]=0.;
    transformationMatrix[j][j]=1.; 
    transformationMatrix0[j][j]=1.;    
  }

  /* WARNING - temporary fix for transformationMatrix0 */ 


  transformationMatrix0[1][1]=0.;   
  transformationMatrix0[2][2]=0.; 
  transformationMatrix0[1][2]=-1.; 
  transformationMatrix0[2][1]=-1.; 

  transformationMatrix0[0][3]=0.; // 126.
  transformationMatrix0[1][3]=0.; // 126.
  transformationMatrix0[2][3]=0.; // 126.

  cout << "Type Isotrack Input File Name -> "; 
  cin >> isotrackInputFile; 

/* input data file should be similar to ic** RAS coordinate format, i.e. 
   nIsotrack
   # rasX, rasY, rasZ */ 


  ifstream isotrackFile(isotrackInputFile.c_str());

  if(!isotrackFile.is_open()) {
    cout << "Unable to open Isotrack file in function isotrack::isotrack" 
	 << endl; 
    abort(); 
  }


  isotrackFile >> nIsotrack; 

  if(!nIsotrack) {
    cout << "Zero number of isotrack points" << endl;
    abort(); 
  }

  isotrackCoordinate0 = new float* [nIsotrack];
  isotrackCoordinate = new float* [nIsotrack];
  for (j=0; j < nIsotrack; j++) {
    isotrackCoordinate0[j] = new float [3]; 
    isotrackCoordinate[j] = new float [3]; 
  }

  for (j=0; j < nIsotrack; j++) 
    for(i=0; i < 3; i++) isotrackFile >> isotrackCoordinate0[j][i];

  isotrackFile.close();


  /* WARNING - adjust all isotrack points */ 


  float v[4]; 
  v[3]=1.; 

  for (k=0; k < nIsotrack; k++) {

    for(i=0; i < 3; i++) v[i]=isotrackCoordinate0[k][i];

    for(j=0; j < 3; j++) {
      isotrackCoordinate[k][j]=0.; 
      for(i=0; i < 4; i++) 
	isotrackCoordinate[k][j]  += transformationMatrix0[j][i]*v[i]; 
    }

    for(i=0; i < 3; i++) isotrackCoordinate0[k][i]=isotrackCoordinate[k][i];
  }


}

/* isotrack destructor */

isotrack::~isotrack(){

  /* delete memory for transformationMatrix */ 

  for (int i=0; i < nRas; i++) 
    delete [] transformationMatrix[i];

  delete [] transformationMatrix;

  /* delete isotrack array */ 

  for (int i=0; i < 3; i++) {
    delete [] isotrackCoordinate0[i];
    delete [] isotrackCoordinate[i];
  }

  delete [] isotrackCoordinate0;
  delete [] isotrackCoordinate;

}

/* main solver for the matrix */ 

void isotrack::getIsotrackMatrix() {

  int i,j, k; 

  float **point; 
  float *value; 
  float *dummy; 

  dummy = new float [nParameter]; 
  dummy[0]=dummy[1]=dummy[2]=angularScale; 
  dummy[3]=dummy[4]=dummy[5]=translationalScale; 

  /* initialize vertices of starting simplex here */ 

  /*
  point = new float* [nParameter+1];
  for (j=0; j < (nParameter+1); j++) {
    point[j] = new float [nParameter];
    for(i=0; i < nParameter; i++) point[j][i]=0.; 
  }

  for(j=0; j < nParameter; j++) point[j][j]=dummy[j];  

  value = new float [nParameter+1]; 

  for(j=0; j < (nParameter+1); j++) {
    for(i=0; i < nParameter; i++) dummy[i]=point[j][i]; 
    value[j]=costFunction(dummy); 
  }

  simplexND(point,value,nParameter); 

  */

  for (k=0; k < nIsotrack; k++) { 
    cout << (k+1) << " "; 
    for(i=0; i < 3; i++) cout << isotrackCoordinate0[k][i] << " ";  
    cout << endl; 
  }

 
}

/* recompute matrix from 3 angular and 3 translational coordinates */ 

void isotrack::computeMatrix() {

  float cosX, sinX, cosY, sinY, cosZ, sinZ; 

  cosX=cos(angle[0]); 
  sinX=sin(angle[0]); 

  cosY=cos(angle[1]); 
  sinY=sin(angle[1]); 

  cosZ=cos(angle[2]); 
  sinZ=sin(angle[2]); 

  transformationMatrix[0][0]=cosY*cosZ;  
  transformationMatrix[0][1]=cosY*sinZ;  
  transformationMatrix[0][2]=sinY; 

  transformationMatrix[1][0]=-sinX*sinY*cosZ-cosX*sinZ;  
  transformationMatrix[1][1]=-sinX*sinY*sinZ+cosX*cosZ;  
  transformationMatrix[1][2]=sinX*cosY; 

  transformationMatrix[2][0]=-cosX*sinY*cosZ+sinX*sinZ;  
  transformationMatrix[2][1]=-cosX*sinY*sinZ-sinX*cosZ;  
  transformationMatrix[2][2]=cosX*cosY; 

  for (int i=0; i < 3; i++) transformationMatrix[i][nRas-1]=translation[i]; 

}

/* rotate and shift all vertices by 6 given parameters */ 

void isotrack::rotateMesh(float aX, float aY, float aZ, 
			  float tX, float tY, float tZ) {

  int i,j,k; 
  float v[4]; 
  v[3]=1.; 

  angle[0]=aX; 
  angle[1]=aY; 
  angle[2]=aZ; 

  translation[0]=tX; 
  translation[1]=tY; 
  translation[2]=tZ; 

  computeMatrix(); 

  for (j=0; j < nVertex; j++) {

    for(i=0; i < 3; i++) v[i]=vertex[j].r[i];

    for(k=0; k < 3; k++) {
      vertex[j].r[k]=0.; 
      for(i=0; i < 4; i++) 
	vertex[j].r[k]  += transformationMatrix[k][i]*v[i]; 
    }
  }

}


/* compute pseudo distance between isotrack points and scalp surface */ 

float isotrack::costFunction(float *dummy) {

  /* define weights here temporarily */

  int i, j, k;
  int volumeIndex;  
  float v[4]; 
  v[3]=1.; 

  float cost=0.; 

  for(i=0; i < 3; i++) {
    angle[i]=dummy[i]; 
    translation[i]=dummy[i+3]; 
  }

  cout << "coor: " << angle[0] << " " << angle[1] << " " << angle[2] << " " 
       << translation[0] << " " << translation[1] << " " << translation[2] 
       << endl; 

  computeMatrix(); 

  /*
  cout << "matrix" << endl; 
  for(j=0; j < 4; j++) {
    for(i=0; i < 4; i++) cout << transformationMatrix[j][i] << " "; 
    cout << endl; 
  }
  */


  /* begin cycle over isotrack points */

  for (k=0; k < nIsotrack; k++) {

    for(i=0; i < 3; i++) v[i]=isotrackCoordinate0[k][i];

    for(j=0; j < 3; j++) {
      isotrackCoordinate[k][j]=0.; 
      for(i=0; i < 4; i++) 
	isotrackCoordinate[k][j]  += transformationMatrix[j][i]*v[i]; 
    }

    /*

    cout << "is. coor " << isotrackCoordinate0[k][0] << " " 
	 << isotrackCoordinate0[k][1] << " " 
	 << isotrackCoordinate0[k][2] << endl; 
    */


    cost += getBorderMap(isotrackCoordinate[k][0], 
			 isotrackCoordinate[k][1], 
			 isotrackCoordinate[k][2]);  
	
  }

  //cout << endl; 

  /* end cycle over isotrack points */

  cout << "cost function: " << cost << endl; 

  return cost; 

}

