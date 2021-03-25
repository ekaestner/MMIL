
#include <iostream> 
#include <sstream> 
#include <string> 
#include <fstream> 
#include <cmath> 
#include <stdlib.h>

using namespace std;

#include "mri.h"
#include "mesh.h"

/* mesh constructor */

mesh::mesh(string mriFileName, 
	   string mriInputDirectory, 
	   int talairachFlag, 
	   string talairachFileName, 
	   string triInputDirectory_ini, 
	   int meshNumber_ini): 
  mri(mriFileName, 
      mriInputDirectory) {

  const int maxMeshNumber=8; 

  /* Initialize smoothing coeffcients here */

  LOCAL=0; 
  GLOBAL=1; 

  coefSmoothLocal=1.; 
  coefSmoothGlobal=1.; 

  weightInside=1; 
  weightBorder=2; 
  weightOutside=0; 

  int i, j, k;
  int nDummy;
  string sDummy; 

  triInputDirectory=triInputDirectory_ini; 
  meshNumber=meshNumber_ini; 

  triFileInfo = new structTriFileInfoType[maxMeshNumber];  

  triFileInfo[0].nFace=20; 
  for(i=1; i < maxMeshNumber; i++) 
    triFileInfo[i].nFace = 4*triFileInfo[i-1].nFace; 
  for(i=0; i < maxMeshNumber; i++) 
    triFileInfo[i].nVertex = triFileInfo[i].nFace/2 + 2;

  nVertex=triFileInfo[meshNumber].nVertex; 
  nFace=triFileInfo[meshNumber].nFace; 

  /* Read tesselated icosahedron vertices and faces */ 
  
  icoCoordinate = new float* [nVertex];
  icoFaceElement = new int* [nFace];

  for (j=0; j < nVertex; j++) icoCoordinate[j] = new float [3]; 
  for (k=0; k < nFace; k++) icoFaceElement[k] = new int [3]; 

  ifstream triFile;
  stringstream converter; 
  converter << meshNumber; 
  sDummy = triInputDirectory + "ic" + converter.str() + ".tri"; 
  triFile.open(sDummy.c_str());

  triFile >> nDummy;
  for (j=0; j < nVertex; j++) {
    triFile >> nDummy; 
    for(i=0; i < 3; i++) triFile >> icoCoordinate[j][i];
  }  

  triFile >> nDummy;
  for(k=0; k < nFace; k++) { 
    triFile >> nDummy; 
    for(i=0; i < 3; i++) triFile >> icoFaceElement[k][i];  
  }

  triFile.close();

  vertex = new structVertexType[nVertex]; 
  face = new structFaceType[nFace]; 

  /* allocate memory for Talairach matrices, read the Talairach file 
     if necessary and compute inverse Talairach matrix */

  talairach = new float* [4];
  inverseTalairach = new float* [4];

  for (int i=0; i < 4; i++) {
    talairach[i] = new float [4];
    inverseTalairach[i] = new float [4]; 
  }

  ifstream talairachFileInput(talairachFileName.c_str());

  if(talairachFileInput.is_open()) {

    for(j=0; j < 3; j++) 
      for(i=0; i < 4; i++) talairachFileInput >> talairach[j][i]; 
 
    talairach[3][0]=0.;
    talairach[3][1]=0.;
    talairach[3][2]=0.;

    talairach[3][3]=1.;

    talairachFileInput.close();
  }

  else { 

    cout << "Unable to open Talairach file in  mesh::mesh -"; 
    cout << " using unity matrix instead" << endl; 

    for(j=0; j < 4; j++) 
      for(i=0; i < 4; i++) talairach[j][i]=0.; 

    for(j=0; j < 4; j++) talairach[j][j]=1.; 

  }

  computeInverseMatrix(talairach, 4, inverseTalairach); 

  /* Initialize the filling arrays */

  int jByte, jRow, jSlice; 

  contourFilling = new int** [nByte];
  for (jByte = 0; jByte < nByte; jByte++) {
    contourFilling[jByte] = new int* [nRow];
    for (jRow = 0; jRow < nRow; jRow++) {
      contourFilling[jByte][jRow] = new int [nSlice];
    }
  }
 
  for (jByte = 0; jByte < nByte; jByte++)
    for (jRow = 0; jRow < nRow; jRow++)
      for (jSlice = 0; jSlice < nSlice; jSlice++)
	contourFilling[jByte][jRow][jSlice]=0; 


  borderMap = new float** [nByte];
  for (jByte = 0; jByte < nByte; jByte++) {
    borderMap[jByte] = new float* [nRow];
    for (jRow = 0; jRow < nRow; jRow++) {
      borderMap[jByte][jRow] = new float [nSlice];
    }
  }
 
  for (jByte = 0; jByte < nByte; jByte++)
    for (jRow = 0; jRow < nRow; jRow++)
      for (jSlice = 0; jSlice < nSlice; jSlice++)
	borderMap[jByte][jRow][jSlice]=0.; 

}




/* mesh destructor */

mesh::~mesh(){

  delete [] triFileInfo; 
  delete [] vertex; 
  delete [] face; 

  /* Destroy both Talairach matrices */ 

  for (int i=0; i < 4; i++) {
    delete [] talairach[i];
    delete [] inverseTalairach[i];
  }

  delete [] talairach;
  delete [] inverseTalairach; 


  /* Destroy the filling arrays */ 

  for(int jByte = 0; jByte < nByte; jByte++)
  {
    for(int jRow = 0; jRow < nRow; jRow++) 
      delete [] contourFilling[jByte][jRow];
    delete [] contourFilling[jByte];
  }
  delete [] contourFilling;

  for(int jByte = 0; jByte < nByte; jByte++)
  {
    for(int jRow = 0; jRow < nRow; jRow++) 
      delete [] borderMap[jByte][jRow];
    delete [] borderMap[jByte];
  }
  delete [] borderMap;



}

  /* Initialize deformable mesh and create a sphere of a given radious */ 

void mesh::initializeSurface(float templateSize) {

  const bool FALSE=0, TRUE=1;

  int i, j, k, m;
  int last, next, skiplast, skipnext; 
  int array1[7], array2[36]; 
  int n1, n2;
  int jFlag;  

  for(j=0; j < nVertex; j++) {
 
    for(i=0; i < 3; i++) {
      vertex[j].r[i]=icoCoordinate[j][i];  
    }
    vertex[j].vnum = 0;
    vertex[j].vnum2 = 0;   
    vertex[j].fnum = 0;
  }

  for(k=0; k < nFace; k++) { 
    for (i=0; i < 3; i++) { 
      face[k].nTri[i]=icoFaceElement[k][i];  
      face[k].nTri[i]--;
    }
  }

  for(k=0; k < nFace; k++) {

    /* Determine first order neighbor vertices */ 
                                                                               
    for (i=0; i < 3; i++) {
                                                                               
      vertex[face[k].nTri[i]].f[vertex[face[k].nTri[i]].fnum++] = k;
      last = (i>0)?i-1:2;
      next = (i<2)?i+1:0;
      skiplast = skipnext = FALSE;
                                                                               
      for (j=0; j < vertex[face[k].nTri[i]].vnum; j++) {
                                                                               
        if (vertex[face[k].nTri[i]].v[j]==face[k].nTri[last]) skiplast = TRUE;
                                                                               
        if (vertex[face[k].nTri[i]].v[j]==face[k].nTri[next]) skipnext = TRUE;
      }
                                                                               
      if (!skiplast)
        vertex[face[k].nTri[i]].v[vertex[face[k].nTri[i]].vnum++]=
	  face[k].nTri[last];
                                                                               
      if (!skipnext)
        vertex[face[k].nTri[i]].v[vertex[face[k].nTri[i]].vnum++]=
	  face[k].nTri[next];
    }

    /* Determine second order neighbor vertices */ 

  }

  for(j=0; j < nVertex; j++) {

    n1=vertex[j].vnum+1;
    array1[0]=j;  
    n2=0; 

    for (m=0; m < vertex[j].vnum; m++) { 
      array1[m+1]=vertex[j].v[m];
 
      for (int k=0; k < vertex[vertex[j].v[m]].vnum; k++) 
	array2[n2+k] =vertex[vertex[j].v[m]].v[k]; 

      n2 += vertex[vertex[j].v[m]].vnum;  
    }

    vertex[j].vnum2=0; 

    for(m=0; m < n2; m++) { 
      jFlag=0; 

      for(k=0; k < n1; k++) if(array2[m]==array1[k]) jFlag++; 

      for(k=0; k < vertex[j].vnum2; k++) 
	if(array2[m]==vertex[j].v2[k]) jFlag++; 

      if(!jFlag) { 
	vertex[j].v2[vertex[j].vnum2]=array2[m]; 
	vertex[j].vnum2++; 
      }
    }
  }

  for (int j=0; j < nVertex; j++) 
    for(int i=0; i < 3; i++) vertex[j].r[i] *= templateSize; 

}

/* perform Talairach transform over the mesh */

void mesh::talairachForward() {

  int i, j, k; 
  float v[4]; 

  for (k=0; k < nVertex; k++) {

    for(i=0; i < 3; i++) v[i]=vertex[k].r[i];
    v[3]=1.; 

    for(j=0; j < 3; j++) {
      vertex[k].r[j]=0.; 
      for(i=0; i < 4; i++) vertex[k].r[j] += talairach[j][i]*v[i]; 
    }	
  }
}

/* perform inverse Talairach transform over the mesh */

void mesh::talairachBackward() {

  int i, j, k; 
  float v[4]; 
  v[3]=1.; 

  for (k=0; k < nVertex; k++) {

    for(i=0; i < 3; i++) v[i]=vertex[k].r[i];

    for(j=0; j < 3; j++) {
      vertex[k].r[j]=0.; 
      for(i=0; i < 4; i++) vertex[k].r[j] += inverseTalairach[j][i]*v[i]; 
    }	
  }
}

/* Establish local coordinate systems */

void mesh::localCoordinateSystem() {

  const float quasiUnity=0.999; 

  float v1[3], v2[3], v[3], product;

  float dummy; 

  /* Compute normals to all triangles */ 


  for(int k=0; k < nFace; k++) { 

    for(int i=0; i < 3; i++) {  
      v1[i] = vertex[face[k].nTri[0]].r[i]-vertex[face[k].nTri[1]].r[i];
      v2[i] = vertex[face[k].nTri[2]].r[i]-vertex[face[k].nTri[1]].r[i];
    }

    face[k].fNormal[0] = v1[1]*v2[2]-v1[2]*v2[1];
    face[k].fNormal[1] = v1[2]*v2[0]-v1[0]*v2[2];
    face[k].fNormal[2] = v1[0]*v2[1]-v1[1]*v2[0];

    product=sqrt(face[k].fNormal[0]*face[k].fNormal[0]+
                 face[k].fNormal[1]*face[k].fNormal[1]+
                 face[k].fNormal[2]*face[k].fNormal[2]); 

    if (product == 0.0) {
      cout << "normal to face is zero" << endl;
      return;
    }
    for(int i=0; i < 3; i++) face[k].fNormal[i] /= product;
  }

  for(int j=0; j < nVertex; j++) {

    /* Compute normals to the vertices */ 

    for(int i=0; i < 3; i++) v[i]=0.; 

    for(int m=0; m < vertex[j].fnum; m++) 
      for(int i=0; i < 3; i++) v[i] += face[vertex[j].f[m]].fNormal[i]; 

    for(int i=0; i < 3; i++) vertex[j].fNormal[i]=v[i]/vertex[j].fnum; 

    product=sqrt(vertex[j].fNormal[0]*vertex[j].fNormal[0]+
                 vertex[j].fNormal[1]*vertex[j].fNormal[1]+
                 vertex[j].fNormal[2]*vertex[j].fNormal[2]); 

    if (product == 0.0) {
      cout << "normal to vertex is zero" << endl;
      return;
    }
    for(int i=0; i < 3; i++) vertex[j].fNormal[i] /= product;


    /* Compute the first tangential vector */ 

    dummy = vertex[j].fNormal[0]* vertex[j].fNormal[0]; 

    /* check if normal is not aligned with (1,0,0) */ 

    if(dummy < quasiUnity) { 
      vertex[j].xTangent[0]=1.-dummy; 
      vertex[j].xTangent[1]=-vertex[j].fNormal[0]* vertex[j].fNormal[1]; 
      vertex[j].xTangent[2]=-vertex[j].fNormal[0]* vertex[j].fNormal[2]; 
    }

    else {
      dummy = vertex[j].fNormal[1]* vertex[j].fNormal[1];
      vertex[j].xTangent[0]=-vertex[j].fNormal[0]* vertex[j].fNormal[1];
      vertex[j].xTangent[1]=1-dummy;  
      vertex[j].xTangent[2]=-vertex[j].fNormal[1]* vertex[j].fNormal[2]; 
    }

    product=sqrt(vertex[j].xTangent[0]*vertex[j].xTangent[0]+
                 vertex[j].xTangent[1]*vertex[j].xTangent[1]+
                 vertex[j].xTangent[2]*vertex[j].xTangent[2]); 

    if (product == 0.0) {
      cout << "tangential to vertex is zero" << endl;
      return;
    }
    for(int i=0; i < 3; i++) vertex[j].xTangent[i] /= product;


    /* Compute the second tangential vector as a vector product */ 

      vertex[j].yTangent[0]=vertex[j].fNormal[1]*vertex[j].xTangent[2]-
	vertex[j].fNormal[2]*vertex[j].xTangent[1]; 
      vertex[j].yTangent[1]=vertex[j].fNormal[2]*vertex[j].xTangent[0]-
	vertex[j].fNormal[0]*vertex[j].xTangent[2]; 
      vertex[j].yTangent[2]=vertex[j].fNormal[0]*vertex[j].xTangent[1]-
	vertex[j].fNormal[1]*vertex[j].xTangent[0]; 

  }

}

/* compute local curvature using quadrics */ 

void mesh::computeCurvature() {

  float a, b, c; 
  float a2, b2, c2; 

  float xLocal[6], yLocal[6], zLocal[6];
  float xLocal2[12], yLocal2[12], zLocal2[12];
  float x[6], y[6], z[6];
  float x2[12], y2[12], z2[12];   
  float x4, x3y, x2y2, xy3, y4, x2z, xyz, y2z;  
  float dummy_x2, dummy_y2, det; 

  for(int j=0; j < nVertex; j++) {

    x4=x3y=x2y2=xy3=y4=x2z=xyz=y2z=0.;

    // Compute curvature 

    for(int k=0; k < vertex[j].vnum; k++) { 
      x[k]=vertex[vertex[j].v[k]].r[0]-vertex[j].r[0]; 
      y[k]=vertex[vertex[j].v[k]].r[1]-vertex[j].r[1]; 
      z[k]=vertex[vertex[j].v[k]].r[2]-vertex[j].r[2]; 

      xLocal[k]=vertex[j].xTangent[0]*x[k]+ 
	        vertex[j].xTangent[1]*y[k]+ 
	        vertex[j].xTangent[2]*z[k];  

      yLocal[k]=vertex[j].yTangent[0]*x[k]+ 
	        vertex[j].yTangent[1]*y[k]+ 
	        vertex[j].yTangent[2]*z[k];  

      zLocal[k]=vertex[j].fNormal[0]*x[k]+ 
	        vertex[j].fNormal[1]*y[k]+ 
	        vertex[j].fNormal[2]*z[k];  


      dummy_x2=xLocal[k]*xLocal[k]; 
      dummy_y2=yLocal[k]*yLocal[k]; 

      x4 += dummy_x2*dummy_x2; 
      x3y +=  dummy_x2*xLocal[k]*yLocal[k]; 
      x2y2 += dummy_x2*dummy_y2; 
      xy3 += xLocal[k]*yLocal[k]*dummy_y2; 
      y4 += dummy_y2*dummy_y2; 

      x2z += dummy_x2*zLocal[k];
      xyz += xLocal[k]*yLocal[k]*zLocal[k];
      y2z += dummy_y2*zLocal[k];
    }

    det = x4*(x2y2*y4-xy3*xy3)-x3y*(x3y*y4-x2y2*xy3)+x2y2*(x3y*xy3-x2y2*x2y2); 
    a = x2z*(x2y2*y4-xy3*xy3)-xyz*(x3y*y4-x2y2*xy3)+y2z*(x3y*xy3-x2y2*x2y2); 
    b = x4*(xyz*y4-x2z*xy3)-x3y*(x2z*y4-y2z*x2y2)+x2y2*(x2z*xy3-xyz*x2y2); 
    c = x4*(x2y2*y2z-xy3*xyz)-x3y*(x3y*y2z-xy3*x2z)+x2y2*(x3y*xyz-x2y2*x2z); 
    
    if(det) { a /= det; b /= det; c /= det; }

    else cout << "Problem with curvature calculation" << endl; 

    vertex[j].curvature=a+c; 

    // Compute curvature2 

    for(int k=0; k < vertex[j].vnum2; k++) { 
      x2[k]=vertex[vertex[j].v2[k]].r[0]-vertex[j].r[0]; 
      y2[k]=vertex[vertex[j].v2[k]].r[1]-vertex[j].r[1]; 
      z2[k]=vertex[vertex[j].v2[k]].r[2]-vertex[j].r[2]; 

      xLocal2[k]=vertex[j].xTangent[0]*x2[k]+ 
	         vertex[j].xTangent[1]*y2[k]+ 
	         vertex[j].xTangent[2]*z2[k];  

      yLocal2[k]=vertex[j].yTangent[0]*x2[k]+ 
	         vertex[j].yTangent[1]*y2[k]+ 
	         vertex[j].yTangent[2]*z2[k];  

      zLocal2[k]=vertex[j].fNormal[0]*x2[k]+ 
	         vertex[j].fNormal[1]*y2[k]+ 
	         vertex[j].fNormal[2]*z2[k];  


      dummy_x2=xLocal2[k]*xLocal2[k]; 
      dummy_y2=yLocal2[k]*yLocal2[k]; 

      x4 += dummy_x2*dummy_x2; 
      x3y +=  dummy_x2*xLocal2[k]*yLocal2[k]; 
      x2y2 += dummy_x2*dummy_y2; 
      xy3 += xLocal2[k]*yLocal2[k]*dummy_y2; 
      y4 += dummy_y2*dummy_y2; 

      x2z += dummy_x2*zLocal2[k];
      xyz += xLocal2[k]*yLocal2[k]*zLocal2[k];
      y2z += dummy_y2*zLocal2[k];
    }

    det = x4*(x2y2*y4-xy3*xy3)-x3y*(x3y*y4-x2y2*xy3)+x2y2*(x3y*xy3-x2y2*x2y2); 
    a2 = x2z*(x2y2*y4-xy3*xy3)-xyz*(x3y*y4-x2y2*xy3)+y2z*(x3y*xy3-x2y2*x2y2); 
    b2 = x4*(xyz*y4-x2z*xy3)-x3y*(x2z*y4-y2z*x2y2)+x2y2*(x2z*xy3-xyz*x2y2); 
    c2 = x4*(x2y2*y2z-xy3*xyz)-x3y*(x3y*y2z-xy3*x2z)+x2y2*(x3y*xyz-x2y2*x2z); 
    
    if(det) { a2 /= det; b2 /= det; c2 /= det; }

    else cout << "Problem with curvature2 calculation" << endl; 

    vertex[j].curvature2=a2+c2; 

  }
} 

/* dilate the surface by a given distance and fill the volume */ 

void mesh::createRepellingVolume(float movingDistance, int nSmoothness) {

  int i, j, k, n; 

  float stepDilate=1., dummy; 

  int nDilate=int(movingDistance/stepDilate); 

  for(n=0; n < nDilate; n++) { 

    for(k=0; k < nSmoothness; k++) dummy = smoothSurface(GLOBAL); // 1st 

    localCoordinateSystem();
    for(j=0; j < nVertex; j++) 
      for(i=0; i < 3; i++)  
	vertex[j].r[i] += stepDilate*vertex[j].fNormal[i]; 

    for(k=0; k < nSmoothness; k++) dummy = smoothSurface(GLOBAL); // 2nd 

  }

  fillContour(); 

} 


void mesh::createBorderMap(float movingDistance, int nSmoothness) {

  int i, j, k, n; 

  int jByte, jRow, jSlice; 

  float stepDilate=1., dummy; 

  int nDilate=int(movingDistance/stepDilate); 

  if(nDilate < 0) {
    nDilate = -nDilate; 
    stepDilate=-1.; 
  }

  for (jByte = 0; jByte < nByte; jByte++)
    for (jRow = 0; jRow < nRow; jRow++)
      for (jSlice = 0; jSlice < nSlice; jSlice++)
	borderMap[jByte][jRow][jSlice]=0.; 



  for(n=0; n < nDilate; n++) { 

    for(k=0; k < nSmoothness; k++) dummy = smoothSurface(GLOBAL); // 1st 

    localCoordinateSystem();
    for(j=0; j < nVertex; j++) 
      for(i=0; i < 3; i++)  
	vertex[j].r[i] += stepDilate*vertex[j].fNormal[i]; 

    for(k=0; k < nSmoothness; k++) dummy = smoothSurface(GLOBAL); // 2nd 

    fillContour(); 

    for (jByte = 0; jByte < nByte; jByte++)
      for (jRow = 0; jRow < nRow; jRow++)
	for (jSlice = 0; jSlice < nSlice; jSlice++)
	  if(contourFilling[jByte][jRow][jSlice] != 0) 
	    borderMap[jByte][jRow][jSlice] += 1.; 
  }

  for (jByte = 0; jByte < nByte; jByte++)
    for (jRow = 0; jRow < nRow; jRow++)
      for (jSlice = 0; jSlice < nSlice; jSlice++)
	borderMap[jByte][jRow][jSlice] /= nDilate; 


  /* WARNING - temporary output */


  ofstream file1("slice1");
  ofstream file2("slice2");

  jByte=128; 

  for (jRow = 0; jRow < nRow; jRow++) {
    for (jSlice = 0; jSlice < nSlice; jSlice++) 
      file1 << borderMap[jByte][jRow][jSlice] << " "; 
    file1 << endl; 
  }	

  jSlice=64; 
  for (jRow = 0; jRow < nRow; jRow++) {
    for (jByte = 0; jByte < nByte; jByte++) 
      file2 << borderMap[jByte][jRow][jSlice] << " "; 
    file2 << endl; 
  }

  file1.close();
  file2.close();
  
  /* end of temporary output */ 

} 


/* Implement global smoothing here 

GLOBAL SMOOTHING     
Option 1: coefSmooth*(vertex[j].curvature - meanCurvature)
Option 2: coefSmooth*(vertex[j].curvature2 - meanCurvature2)

LOCAL SMOOTHING 
Option 1: coefSmooth*(vertex[j].curvature1 - vertex[j].curvature2) */

float mesh::smoothSurface(int typeSmooth) {

  int k, j, i;

  float meanCurvature, meanCurvature2;  
  float dummy; 
  float globalDisplacement=0., vertexDisplacement; 

  /* Global smoothing */ 

  localCoordinateSystem();
  computeCurvature();

  switch(typeSmooth) {

  /* Local smoothing */ 

  case 0: 

    for(j=0; j < nVertex; j++) { 
    
      vertexDisplacement=0.; 

      for(i=0; i < 3; i++) { 
	dummy = coefSmoothLocal*(vertex[j].curvature - vertex[j].curvature2)*
	  vertex[j].fNormal[i]; 
	vertex[j].r[i] += dummy; 
	vertexDisplacement += dummy*dummy; 
      }

      globalDisplacement += sqrt(vertexDisplacement); 

    }

    break; 

  /* Global smoothing */ 

  case 1: 

    meanCurvature=0.; 
    meanCurvature2=0.;

    for(j=0; j < nVertex; j++) {
      meanCurvature += vertex[j].curvature;
      meanCurvature2 += vertex[j].curvature2; 
    }

    meanCurvature /= nVertex;
    meanCurvature2 /= nVertex; 

    for(j=0; j < nVertex; j++) { 
    
      vertexDisplacement=0.; 

      for(i=0; i < 3; i++) { 
	dummy = coefSmoothGlobal*(vertex[j].curvature2 - meanCurvature2)*
	  vertex[j].fNormal[i]; 
	vertex[j].r[i] += dummy; 
	vertexDisplacement += dummy*dummy; 
      }

      globalDisplacement += sqrt(vertexDisplacement); 

    }


    break; 


  default: 

    cout << "smoothSurface: smoothing info missing" << endl; 

  }

  return globalDisplacement/nVertex; 

} 


/* WARNING - check for the optimal value of step and if 
   the use of dmax is necessary */ 


/* mark inside voxels by 0, and outside - by 1 and rearrange the weights */ 


float mesh::fillContour() {

  const float step=0.25; 

  int jFlag; 

  int jByte, jRow, jSlice; 
  int n1, n2;
  int jVoxel[nRas]; 

  float fVoxel[nRas]; 
  float ras[nRas], ras0[nRas];   
  float dist1, dist2; 
  float vect1[3], vect2[3]; 
  float unit1[3], unit2[3]; 

  ras[3]=ras0[3]=1.; 

  for (jByte = 0; jByte < nByte; jByte++)
    for (jRow = 0; jRow < nRow; jRow++)
      for (jSlice = 0; jSlice < nSlice; jSlice++)
	contourFilling[jByte][jRow][jSlice]=0; 

  /* Mark all voxels transected by the surface */ 

  for(int kFace=0; kFace < nFace; kFace++) { 

    dist1=dist2=0.; 

    for(int i=0; i < 3; i++) {

      ras0[i]= vertex[face[kFace].nTri[0]].r[i]; 
      
      vect1[i] = vertex[face[kFace].nTri[1]].r[i]- 
	vertex[face[kFace].nTri[0]].r[i]; 
      vect2[i] = vertex[face[kFace].nTri[2]].r[i]- 
	vertex[face[kFace].nTri[0]].r[i]; 

      dist1 += vect1[i]*vect1[i]; 
      dist2 += vect2[i]*vect2[i]; 
    }

    dist1 = sqrt(dist1); 
    dist2 = sqrt(dist2); 

    n1=int(dist1/step)+1; 
    for(int i=0; i < 3; i++) unit1[i] = vect1[i]/n1; 

    for(int j1=0; j1 <= n1; j1++) { 
      n2=int(dist2*(n1-j1)/step/n1)+1; 
      for(int i=0; i < 3; i++) unit2[i] = vect2[i]/n2;       

      for(int j2=0; j2 < n2; j2++) { 
	for(int i=0; i < 3; i++) ras[i] = ras0[i] + j1*unit1[i] + j2*unit2[i];

  /* WARNING - check if need to adjust RAS coordinates here to comply 
     with FreeSurfer */ 
 
	for(int i=0; i < 3; i++) ras[i] += rasCenter[i]; 

	for(int k=0; k < nRas; k++) {
	  fVoxel[k]=0.;
	  for(int i=0; i < nRas; i++) {
	    fVoxel[k] += fMRas2Vox[k][i]*ras[i];
	  }
	  jVoxel[k]=int(fVoxel[k]);
	}
  
	if( (jVoxel[0] < nByte) && (jVoxel[0] >= 0) &&
	    (jVoxel[1] < nRow) && (jVoxel[1] >= 0) &&
	    (jVoxel[2] < nSlice) && (jVoxel[2] >= 0) ) {

	  contourFilling[jVoxel[0]][jVoxel[1]][jVoxel[2]]=2;
	}

	else {

	  cout << "Surface is outside of volume " << endl; 

        cout <<  "ras: " << ras[0] << " " << ras[1] << " " << ras[2] << endl; 

        cout << "jVox: " << jVoxel[0] << " " << jVoxel[1] << " " << 
	  jVoxel[2] << endl; 

	}
      } 
    }
  }


  /* mark all voxels inside the surface */ 

  int totalfilled, newfilled; 

  contourFilling[1][1][1] = 1;

  totalfilled = newfilled = 1;


  while(newfilled > 0)

    {
    newfilled = 0;

    for (jByte=1; jByte < (nByte-1); jByte++)
      for (jRow=1; jRow < (nRow-1); jRow++)
	for (jSlice=1; jSlice < (nSlice-1); jSlice++)
	  if (contourFilling[jByte][jRow][jSlice]==0)
	    if (contourFilling[jByte-1][jRow][jSlice]==1 ||  
		contourFilling[jByte][jRow-1][jSlice]==1 ||
		contourFilling[jByte][jRow][jSlice-1]==1) {

		contourFilling[jByte][jRow][jSlice] = 1;
		newfilled++;
	      }


    for (jByte=(nByte-2); jByte >= 1 ; jByte--)
      for (jRow=(nRow-2); jRow >= 1; jRow--)
	for (jSlice=(nSlice-2); jSlice >= 1; jSlice--)
	  if (contourFilling[jByte][jRow][jSlice]==0)
	    if (contourFilling[jByte+1][jRow][jSlice]==1 ||
                contourFilling[jByte][jRow+1][jSlice]==1 ||
                contourFilling[jByte][jRow][jSlice+1]==1) {

		contourFilling[jByte][jRow][jSlice] = 1;
		newfilled++;
	      }

	      totalfilled += newfilled; 

    }

  for (jByte = 0; jByte < nByte; jByte++)
    for (jRow = 0; jRow < nRow; jRow++)
      for (jSlice = 0; jSlice < nSlice; jSlice++) {

	switch(contourFilling[jByte][jRow][jSlice]) {

	case 0: 
	  contourFilling[jByte][jRow][jSlice]=weightInside; 
          break;   

	case 1: 
	  contourFilling[jByte][jRow][jSlice]=weightOutside; 
          break;   

	case 2: 
	  contourFilling[jByte][jRow][jSlice]=weightBorder; 
          break;   

	default: 
	  cout << "Problem is contour filling" << endl; 
          abort(); 

	}
	
	if(jByte==0 || jRow==0 || jSlice==0 || 
	   jByte==(nByte-1) || jRow==(nRow-1) || jSlice==(nSlice-1) )   
	  contourFilling[jByte][jRow][jSlice]=weightOutside; 

      }

  /* WARNING - temporary output 


  ofstream file1("slice1");
  ofstream file2("slice2");

  jByte=128; 

  for (jRow = 0; jRow < nRow; jRow++) {
    for (jSlice = 0; jSlice < nSlice; jSlice++) 
      file1 << contourFilling[jByte][jRow][jSlice] << " "; 
    file1 << endl; 
  }	

  jSlice=64; 
  for (jRow = 0; jRow < nRow; jRow++) {
    for (jByte = 0; jByte < nByte; jByte++) 
      file2 << contourFilling[jByte][jRow][jSlice] << " "; 
    file2 << endl; 
  }

  file1.close();
  file2.close();
  
  end of temporary output */ 



  /* compute volume of the filled area */ 

  /* WARNING - do not use the sides of the cube */ 

  float volume=0; 

  for (jByte = 1; jByte < (nByte-1); jByte++)
    for (jRow = 1; jRow < (nRow-1); jRow++)
      for (jSlice = 1; jSlice < (nSlice-1); jSlice++) {
	if(contourFilling[jByte][jRow][jSlice]==weightInside) volume += 1.; 
	if(contourFilling[jByte][jRow][jSlice]==weightBorder) volume += 0.5;
      }

  return voxelVolume*volume; 
}


/* export fillied volume */

void mesh::exportFilling(string surfaceName, 
			 string outputDirectory) {

  int jByte, jRow, jSlice; 

  stringstream converter; 
  converter << meshNumber; 

  string volumeOutputFileName=outputDirectory +  
    surfaceName + converter.str()+ "_volume.dat";

  ofstream volumeOutputFile(volumeOutputFileName.c_str());


  for(jSlice=0; jSlice < nSlice; jSlice++)  
    for (jRow = 0; jRow < nRow; jRow++) {
      for (jByte = 0; jByte < nByte; jByte++) 
	volumeOutputFile << 
	  bool(contourFilling[jByte][jRow][jSlice]) << " "; 
      volumeOutputFile << endl; 
    }

  volumeOutputFile.close();

}

/* determine if voxel is marked in contourFilling array */ 

int mesh::getFilling(float r, float a, float s) {

  int jVoxel[4]; 

  float ras[4], fVoxel[4];  
  ras[0]=r; ras[1]=a; ras[2]=s; ras[3]=1.; 


  /* WARNING - adjust RAS coordinates here to comply with FreeSurfer */ 

  for(int i=0; i < 3; i++) ras[i] += rasCenter[i]; 

  for(int k=0; k < nRas; k++) {
    fVoxel[k]=0.;
    for(int i=0; i < nRas; i++) {
      fVoxel[k] += fMRas2Vox[k][i]*ras[i];
    }
    jVoxel[k]=int(fVoxel[k])-1;
  }

  
  if( (jVoxel[0] < nRow)   && (jVoxel[0] >= 0) &&
      (jVoxel[1] < nByte)  && (jVoxel[1] >= 0) &&
      (jVoxel[2] < nSlice) && (jVoxel[2] >= 0) ) {

      return contourFilling[jVoxel[0]][jVoxel[1]][jVoxel[2]];
  }

  /* image intensity is zero outside of the volume */ 
                                                                               
  else {

    /* WARNING */

                                                                               
    return weightOutside;
  }

  cout << "problem in getFilling subroutine" << endl;

}

/* similar to getFilling */


float mesh::getBorderMap(float r, float a, float s) {

  int jVoxel[4]; 

  float ras[4], fVoxel[4];  
  ras[0]=r; ras[1]=a; ras[2]=s; ras[3]=1.; 


  /* WARNING - adjust RAS coordinates here to comply with FreeSurfer */ 

  for(int i=0; i < 3; i++) ras[i] += rasCenter[i]; 

  for(int k=0; k < nRas; k++) {
    fVoxel[k]=0.;
    for(int i=0; i < nRas; i++) {
      fVoxel[k] += fMRas2Vox[k][i]*ras[i];
    }
    jVoxel[k]=int(fVoxel[k])-1;
  }

  
  if( (jVoxel[0] < nRow)   && (jVoxel[0] >= 0) &&
      (jVoxel[1] < nByte)  && (jVoxel[1] >= 0) &&
      (jVoxel[2] < nSlice) && (jVoxel[2] >= 0) ) {

      return borderMap[jVoxel[0]][jVoxel[1]][jVoxel[2]];
  }

  /* image intensity is zero outside of the volume */ 
                                                                               
  else {

    /* WARNING */

                                                                               
    return 0.;
  }

  cout << "problem in getMap subroutine" << endl;

}


/* export vertex coordinates and normal vectors */

void mesh::exportSurface(int rasShift, 
			 int surfaceDecimation, 
			 string surfaceName, 
			 string outputDirectory) {

  enum {NO=0, YES}; 
  string sDummy;  

  /* write out bem_surface file */ 

  if(surfaceDecimation != YES) { 

    stringstream converter; 
    converter << meshNumber; 

    string bemOutputFileName=outputDirectory +  
      surfaceName + converter.str()+ ".tri";

    ofstream bemOutputFile(bemOutputFileName.c_str());

#ifdef NORMAL_DISK_COPY
    string normalOutputFileName=outputDirectory +
      surfaceName + converter.str()+ ".nrm"; 
    ofstream normalOutputFile(normalOutputFileName.c_str());
#endif /* NORMAL_DISK_COPY */

    bemOutputFile << nVertex << endl; 

    for(int j=0; j < nVertex; j++) {

      bemOutputFile << (j+1) << " ";

#ifdef NORMAL_DISK_COPY  
      normalOutputFile << (j+1) << " "; 
#endif /* NORMAL_DISK_COPY */

      /* decide on RAS shift */ 

      for(int i=0; i <3; i++) {

	switch (rasShift) {

	case YES: 
	  bemOutputFile << (vertex[j].r[i]-rasCenter[i]) << " "; 
	  break; 

	case NO:
	  bemOutputFile << vertex[j].r[i] << " "; 
	  break; 

	default: 
	  cout << "Wrong rasShift value" << endl; 
	  break; 

	}

#ifdef NORMAL_DISK_COPY 
	normalOutputFile << vertex[j].fNormal[i] << " "; 
#endif /* NORMAL_DISK_COPY */

      }
      bemOutputFile << endl; 

#ifdef NORMAL_DISK_COPY 
      normalOutputFile << endl; 
#endif /* NORMAL_DISK_COPY */
    }

    bemOutputFile << nFace << endl; 

    for(int k=0; k < nFace; k++) {
      bemOutputFile << (k+1) << " "; 
      for(int i=0; i < 3; i++) {
	bemOutputFile << (face[k].nTri[i]+1) << " "; 
      }
      bemOutputFile << endl; 
    }
    bemOutputFile.close();

#ifdef NORMAL_DISK_COPY  
    normalOutputFile.close(); 
#endif /* NORMAL_DISK_COPY */
  }

  /* write out bem_decim_surface file */

  if(surfaceDecimation != NO) { 

    stringstream converterDecim; 
    converterDecim << (meshNumber-1); 

    string bemDecimOutputFileName=outputDirectory +  
      surfaceName + converterDecim.str()+ ".tri"; 
    ofstream bemDecimOutputFile(bemDecimOutputFileName.c_str());

#ifdef NORMAL_DISK_COPY 
    string normalDecimOutputFileName=outputDirectory +
      surfaceName + converterDecim.str()+ ".nrm"; 
    ofstream normalDecimOutputFile(normalDecimOutputFileName.c_str());
#endif /* NORMAL_DISK_COPY */

  /* input file name for decimated deformable model */ 

    ifstream triFile;
    sDummy = triInputDirectory + "ic" + converterDecim.str() + ".tri"; 
    triFile.open(sDummy.c_str());

    int nVertexDecim=triFileInfo[meshNumber-1].nVertex; 
    int nFaceDecim=triFileInfo[meshNumber-1].nFace;  

    getline(triFile, sDummy);
    bemDecimOutputFile << sDummy << endl; 

    for(int j=0; j < nVertexDecim; j++) {

      getline(triFile, sDummy);
      bemDecimOutputFile << (j+1) << " "; 

#ifdef NORMAL_DISK_COPY
      normalDecimOutputFile << (j+1) << " "; 
#endif /* NORMAL_DISK_COPY */

      /* decide on RAS shift */ 

      for(int i=0; i <3; i++) {

	switch (rasShift) {

	case YES: 
	  bemDecimOutputFile << (vertex[j].r[i]-rasCenter[i]) << " "; 
	  break; 

	case NO:
	  bemDecimOutputFile << vertex[j].r[i] << " "; 
	  break; 

	default: 
	  cout << "Wrong rasShift value" << endl; 
	  break; 

	}
#ifdef NORMAL_DISK_COPY
	normalDecimOutputFile << vertex[j].fNormal[i] << " "; 
#endif /* NORMAL_DISK_COPY */

      }
      bemDecimOutputFile << endl; 

#ifdef NORMAL_DISK_COPY
      normalDecimOutputFile << endl; 
#endif /* NORMAL_DISK_COPY */
    }

    for(int k=0; k < (nFaceDecim+1); k++) {
      getline(triFile, sDummy);
      bemDecimOutputFile << sDummy << endl; 
    }
    bemDecimOutputFile.close(); 
 
#ifdef NORMAL_DISK_COPY 
    normalDecimOutputFile.close(); 
#endif /* NORMAL_DISK_COPY */

    triFile.close();
  }

}


/* Read a pre-computed surface */ 

void mesh::importSurface(string surfaceName, 
			 string surfaceInputDirectory) {

  int i, j, k;
  int nDummy; 
  string sDummy; 

  stringstream converter; 
  converter << meshNumber; 
  sDummy = surfaceInputDirectory + surfaceName + converter.str() + ".tri"; 

  ifstream surfaceInputFile(sDummy.c_str());

  if(!surfaceInputFile.is_open()) {
    cout << "Unable to open surface input file in function "  
	 <<  "mesh::importSurface" << endl; 
    abort(); 
  }

  surfaceInputFile >> nDummy;
  for (j=0; j < nVertex; j++) {
    surfaceInputFile >> nDummy; 
    for(i=0; i < 3; i++) surfaceInputFile >> vertex[j].r[i];
  }  

  surfaceInputFile.close();

}

/* compute repelling force */ 

void mesh::computeRepellingForce(float coefRepelling) {

  int i,j; 

  const float stepBack=1.;   

  float coorBack[3]; 
  float forceRepelling[3]; 

  if(coefRepelling != 0.) {

    for(j=0; j < nVertex; j++) {

      for(i=0; i <3; i++)
	coorBack[i]=vertex[j].r[i]-stepBack*vertex[j].fNormal[i]; 

      if(getFilling(coorBack[0],coorBack[1],coorBack[2])==weightInside)
	for(i=0; i <3; i++) {
	  vertex[j].repellingForce[i]=coefRepelling*vertex[j].fNormal[i];
	  vertex[j].force[i]=vertex[j].repellingForce[i];  
	} 
   
      else
	for(i=0; i <3; i++) {
	  vertex[j].repellingForce[i]=0.;
	  vertex[j].force[i]=0.;  
	} 
    } 
  }

  else {

    for(j=0; j < nVertex; j++) 
      for(i=0; i < 3; i++) {
	vertex[j].repellingForce[i]=0.; 
	vertex[j].force[i]=0.; 
      }

  }

}

/* add normal and tangential components to the force */ 

void mesh::addInertialForce(float coefTangential,
			    float coefNormal) {

  int i,j,m; 

  float delta[6][3]; 
  float product[6];
  float forceTangential[3], forceNormal[3]; 
  
  for(j=0; j < nVertex; j++) {

     for(i=0; i < 3; i++) { 
      forceTangential[i]=0.; 
      forceNormal[i]=0.; 
    }

     /* compute normal and tangential forces */

    for(m=0; m < vertex[j].fnum; m++) { 

      product[m]=0.; 
      for(i=0; i <3; i++) {
	delta[m][i]=vertex[vertex[j].v[m]].r[i]-vertex[j].r[i];
        product[m] += delta[m][i]*vertex[j].fNormal[i];  
      }

      for(i=0; i <3; i++) {
	forceTangential[i] += delta[m][i] - product[m]*vertex[j].fNormal[i]; 
	forceNormal[i] += product[m]*vertex[j].fNormal[i];
      }
    }

    for(i=0; i < 3; i++)
      vertex[j].force[i] += 
	(coefTangential*forceTangential[i] + coefNormal*forceNormal[i]); 

  }
}

/* to be called by simplexND to compute cost function for a test vector */

float mesh::trial(float **point, float *value, float *psum, 
		  int indexHigh, float factor, int nParameter)
{

  int j; 
  float factor1, factor2, testValue; 

  float* testPoint = new float [nParameter]; 
 
  factor1=(1.-factor)/nParameter; 
  factor2=factor1-factor; 

  for (j=0; j < nParameter; j++) 
    testPoint[j]=psum[j]*factor1-point[indexHigh][j]*factor2;

  testValue=costFunction(testPoint);

  if (testValue < value[indexHigh]) {

    value[indexHigh]=testValue; 

    for (j=0; j < nParameter; j++) { 
      psum[j] += testPoint[j]-point[indexHigh][j]; 
      point[indexHigh][j]=testPoint[j]; 
    } 
  } 

  delete [] testPoint; 

  return testValue; 

}

/* downhill simplex method of Nelder and Mead */ 

void mesh::simplexND(float **point, float *value, int nParameter)

{
  /* WARNING */

  const int nMax=200; 

  //const int nMax=400; 

  int i, j, ihi, ilo, inhi; 

  float dummy, savedValue, testValue; 

  float* psum = new float [nParameter]; 
   
  for (j=0; j < nParameter; j++) {
    psum[j]=0.; 
    for (i=0; i < (nParameter+1); i++) psum[j] += point[i][j];
  }

  for (int n=0; ; n++) { 
    ilo=0;
    ihi = value[0] > value[1] ? (inhi=1,0) : (inhi=0,1);

    for (i=0; i < (nParameter+1); i++) {
      if (value[i] <= value[ilo]) ilo=i; 
      if (value[i] > value[ihi]) { inhi=ihi; ihi=i; } 
      else if (value[i] > value[inhi] && i != ihi) inhi=i; 
    }

    if (n > nMax ) {

      dummy=value[0]; 
      value[0]=value[ilo]; 
      value[ilo]=dummy;       

      for (i=0; i < nParameter; i++) { 
	dummy=point[0][i]; 
	point[0][i]=point[ilo][i]; 
	point[ilo][i]=dummy;   
      }
	break; 
    } 

    testValue=trial(point,value,psum,ihi,-1.,nParameter); 

    if (testValue <= value[ilo]) 
      testValue=trial(point,value,psum,ihi,2.,nParameter);

    else if(testValue >= value[inhi]) {

      savedValue=value[ihi]; 
      testValue=trial(point,value,psum,ihi,0.5,nParameter);

      if (testValue >= savedValue) {

	for (i=0; i < (nParameter+1); i++) {
	  if (i != ilo) { 
	    for (j=0; j < nParameter; j++) point[i][j]=psum[j]
					     =0.5*(point[i][j]+point[ilo][j]);
	    value[i]=costFunction(psum);
	  } 
	} 

	for (j=0; j< nParameter;j++) {
	  psum[j]=0.; 
	  for (i=0; i < (nParameter+1); i++) psum[j] += point[i][j];
	}
      } 
    } 
  }

  delete [] psum; 

}

/* abstract function to be never called */ 

float mesh::costFunction(float *dummy) {

  cout << "mesh::costFunction called by mistake" << endl; 
  abort(); 

}
