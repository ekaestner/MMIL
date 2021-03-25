
#include <iostream> 
#include <sstream> 
#include <string> 
#include <fstream> 
#include <cmath> 

using namespace std;

#include "mri.h"
#include "mesh.h"
#include "statistical.h" 

/* statistical constructor */

statistical::statistical(string mriFileName, 
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

  stepTemplate=1.; 

    meanSurface = new float* [nVertex];
    for (int j=0; j < nVertex; j++) meanSurface[j] = new float [3]; 
   
}

/* statistical destructor */

statistical::~statistical(){

  // WARNING - delete mean surface here  

}


/* move all vertices by a variable distance determined by templates */ 

float statistical::adjustSurface(float coefTangential, 
				 float coefNormal, 
				 float coefTemplate, 
				 float coefRepelling, 
				 int nIntegral, 
				 int nTemplate) {
  int j, i; 

  float dummy; 
  float globalDisplacement, vertexDisplacement; 

  localCoordinateSystem();

  computeShift(nIntegral, nTemplate);

  computeRepellingForce(coefRepelling); 

  addInertialForce(coefTangential, coefNormal);  

  globalDisplacement=0.;
 
  for(int j=0; j < nVertex; j++) { 
    
    vertexDisplacement=0.; 

    for(int i=0; i < 3; i++) { 
      dummy = vertex[j].force[i] + 
	coefTemplate*vertex[j].displacement*vertex[j].fNormal[i];
      vertex[j].r[i] += dummy; 
      vertexDisplacement += dummy*dummy; 
    }

    globalDisplacement += sqrt(vertexDisplacement); 

  }

  return globalDisplacement/nVertex; 

} 


/* Use a template to the location of the boundary between tissues in the 
   local coordinate system. For every vertex, move a template along the 
   normal back and forth by 20 mm to find the point where the integral of  
   [(f(s)-template(s))/sigma(s)]**2 is minimum. */ 


void statistical::computeShift(int nIntegral, int nTemplate) {

  int nFrame=(nTemplateMax-nTemplate)/2; 

  int nMargin=(nTemplate-1)/2; 
  int nHalfIntegral=(nIntegral-1)/2; 

  int nData=nIntegral+2*nMargin; 

  int kValue=0;  

  int limit=nMargin; // define maximum displacement in either direction 

  int interpolation=limit-1; // define interpolation area 

  int i,j,k; 
   
  float* fData; 
  float* fIntegral; 
  float* tVal; 
  float* tDev; 

  float shift, coor[3]; 
  float start=-stepTemplate*(nData-1)/2.;  
  float f, f1, f2, fNorm, value; 

  float fLimit=limit*stepTemplate;
  float dummy, fInt_1, fInt0, fInt1, coef; 

  fData = new float [nData];
  fIntegral = new float [nIntegral];
  tVal = new float [nTemplate]; 
  tDev = new float [nTemplate]; 

  localCoordinateSystem(); 

  /* WARNING */ 
  //ofstream file5("displacement.dat");

  for(j=0; j < nVertex; j++) {

    /* select template values */ 

    for(k=0; k < nTemplate; k++) {
      tVal[k]=templateVal[j][k+nFrame];    
      tDev[k]=templateDev[j][k+nFrame]; 
    }

    /* pre-compute coordinates and date along the normal to the vertex */ 

    for(k=0; k < nData; k++) { 
      shift = start + k*stepTemplate; 

      for(i=0; i <3; i++)
	coor[i]=vertex[j].r[i]+shift*vertex[j].fNormal[i]; 
      fData[k]=getIntensity(coor[0], coor[1], coor[2]); 

      /* WARNING - PADDING */ 

      //if(fData[k] == 0.) fData[k]=vertex[j].templateVal[nTemplate-1]; 


    }

    /* find convolution of the template and data */
    /* find the minimum, keep fIntegral array for future modifications */  

    for(k=0; k < nIntegral; k++) { 

      fNorm=0.; 
      f1=0.;
      f2=0.; 

      if(1) {
	f1=(fData[k]-tVal[0]);
	//f1=(fData[k]-templateVal[0])/templateDev[0];
	fNorm += 0.5; 
      }

      if(1) { 
	f2=(fData[k+nTemplate-1]-tVal[nTemplate-1]);
	//f2=(fData[k+nTemplate-1]-templateVal[nTemplate-1])/
	//templateDev[nTemplate-1];
	fNorm += 0.5; 
      }

      fIntegral[k]=(f1*f1+f2*f2)/2.; 
      
      for(i=1; i < (nTemplate-1); i++) {

	f=0.; 
	if(1) {
	  f=(fData[k+i]-tVal[i]);
	  //f=(fData[k+i]-tVal[i])/tDev[i];
	  fNorm += 1.; 
	}

	fIntegral[k] += f*f; 
      }

      fIntegral[k] /= fNorm; 
	
      if(k==0) {
	value=fIntegral[k]; 
	kValue=0; 
      }

      if(fIntegral[k] < value) { 
	value=fIntegral[k]; 
	kValue=k; 
      }
    }

  
    /* find the first guess for the displacement 
       adjust it if outside of [-fLimit, fLimit] or 
       interpolate if inside of [-fLimit, fLimit] */ 

    dummy=(kValue-nHalfIntegral)*stepTemplate; 
    vertex[j].displacement=dummy;

    /*


    if(dummy >= fLimit) vertex[j].displacement=fLimit; 

    else if(dummy <= (-fLimit)) vertex[j].displacement=-fLimit; 

    else if( ((kValue-nHalfIntegral) <= interpolation) && 
             ((kValue-nHalfIntegral) >= -interpolation) ) {

      fInt_1=fIntegral[kValue-1]; 
      fInt0=fIntegral[kValue]; 
      fInt1=fIntegral[kValue+1]; 

      if(fInt_1 != fInt1) { 

	coef=(fInt_1-fInt1)/(fInt_1-2.*fInt0+fInt1); 
	vertex[j].displacement += stepTemplate*coef/2.;      

      }

    }

    */

    /*

    if(j==9213) { 
      cout << "j= " << j << endl; 
      cout << "kValue " << kValue << endl; 
      cout << "fIntegral[kValue] " << fIntegral[kValue] << endl; 
      cout << "Displacement " << vertex[j].displacement << endl; 

      cout << endl; 

      cout << "Integral" << endl; 
      for(k=0; k < nIntegral; k++) cout << (k-20) << " "; 
      cout << endl; 
      for(k=0; k < nIntegral; k++) cout << fIntegral[k] << " "; 
      cout << endl; 

      cout << "Template" << endl; 
      for(k=0; k < nTemplate; k++) cout << (k-20) << " "; 
      cout << endl; 
      for(k=0; k < nTemplate; k++) cout << vertex[j].templateVal[k] << " "; 
      cout << endl; 

      cout << endl; 

      cout << "Data" << endl; 
      for(k=0; k < nData; k++) cout << (k-40) << " "; 
      cout << endl; 
      for(k=0; k < nData; k++) cout << fData[k] << " "; 
      cout << endl; 

    }

    file5 << j  << " " << vertex[j].displacement << endl; 
    */
                                     
  }


  /* WARNING */  
  //file5.close(); 

  delete [] fData; 
  delete [] fIntegral; 
  delete [] tVal; 
  delete [] tDev; 

}


/* import pre-computed template values for all verteces */

void statistical::importMeanTemplate(string surfaceType, 
				     string templateInputDirectory) {

  int nDummy; 
  int j, k;  

  ifstream templateInputFile;

  stringstream converter; 
  converter << meshNumber; 

  string templateInputFileName = templateInputDirectory +  
    "mean" + surfaceType + converter.str()+ "_template.dat";

  templateInputFile.open(templateInputFileName.c_str());

  /* read the mean template file and allocate memory */  

  nTemplateMax=0; 
  templateInputFile >> nTemplateMax; 

  if(nTemplateMax != 0) { 

    templateVal = new float* [nVertex];
    templateDev = new float* [nVertex];

    for (j=0; j < nVertex; j++) {
      templateVal[j] = new float [nTemplateMax]; 
      templateDev[j] = new float [nTemplateMax]; 
    }
    
    for(j=0; j < nVertex; j++) {
      
      templateInputFile >> nDummy;    
      for(k=0; k < nTemplateMax; k++) 
	templateInputFile >>  templateVal[j][k];       
      
      templateInputFile >> nDummy;    
      for(k=0; k < nTemplateMax; k++) 
	templateInputFile >>  templateDev[j][k]; 
    }
  }
   
  else { 
    cout << "Problem in reading mean template file" << endl;
    abort(); 
  }

  templateInputFile.close(); 

}

/* import pre-computed mean surface */

void statistical::importMeanSurface(string surfaceType, 
				     string meanSurfaceInputDirectory) {

  int i,j; 

  string surfaceName="mean" + surfaceType; 

  importSurface(surfaceName, meanSurfaceInputDirectory);

  for(j=0; j < nVertex; j++) 
    for(i=0; i < 3; i++) meanSurface[j][i]=vertex[j].r[i];

}


/* import pre-computed EOFs */

void statistical::importEof(string surfaceType, 
			    string eofInputDirectory) {

  int nDummy; 
  int i,j,k,n; 
  
  float v1[4], v0[4], fnorm[3]; 
  v1[3]=1.;
  v0[3]=1.;  

  ifstream eofInputFile;

  stringstream converter; 
  converter << meshNumber; 

  string eofInputFileName = eofInputDirectory +  
    "eof" + surfaceType + converter.str()+ ".dat";

  eofInputFile.open(eofInputFileName.c_str());

  /* read the eof file and allocate memory */  

  nEof=0; 
  eofInputFile >> nEof; 

  if(nEof != 0) { 

    coefEof = new float [nEof]; 

    eof = new float** [nEof];

    for (n=0; n < nEof; n++) {

      eof[n] = new float* [nVertex];
      for (j=0; j < nVertex; j++) eof[n][j] = new float [3];
    }

    for(n=0; n < nEof; n++) {

      for(j=0; j < nVertex; j++) {
	eofInputFile >> nDummy; 
	for(i=0; i < 3; i++) eofInputFile >> eof[n][j][i]; 

	/* apply backward Talairach transform to the EOF */

	for(i=0; i < 3; i++) {
	  v1[i]=eof[n][j][i];
	  v0[i]=0.;
	}

	for(k=0; k < 3; k++) {
	  eof[n][j][k]=0.; 
	  fnorm[k]=0.;
	  for(i=0; i < 4; i++) { 
	    eof[n][j][k] += inverseTalairach[k][i]*v1[i]; 
	    fnorm[k] += inverseTalairach[k][i]*v0[i]; 
	  }	
	} 

	for(i=0; i < 3; i++) eof[n][j][i] -= fnorm[i]; 

      }
    }
  }
   
  else { 
    cout << "Problem in reading eof file" << endl;
    abort(); 
  }

  eofInputFile.close(); 

}



/* compute and export template values for all verteces */

void statistical::exportTemplate(string surfaceName, 
				 string outputDirectory, 
				 int nTemplateOut) {

  string sDummy;  
  float shift, coor[3];

  float start=-stepTemplate*(nTemplateOut-1)/2.;  

  stringstream converter; 
  converter << meshNumber; 

  string templateOutputFileName=outputDirectory +  
    surfaceName + converter.str()+ "_template.dat";

  ofstream templateOutputFile(templateOutputFileName.c_str());

  templateOutputFile << nTemplateOut << endl; 

  localCoordinateSystem(); 

  for(int j=0; j < nVertex; j++) {

    templateOutputFile << (j+1) << " "; 

    for(int k=0; k < nTemplateOut; k++) { 
      shift = start + k*stepTemplate; 

      for(int i=0; i <3; i++)
	coor[i]=vertex[j].r[i]+shift*vertex[j].fNormal[i]; 

      templateOutputFile << getIntensity(coor[0], coor[1], coor[2]) << " ";  
    }

    templateOutputFile << endl; 

  }

  templateOutputFile.close();
}

/* optimize coefficients for the EOFs */

void statistical::optimizeCoefficient() {

  int i,j; 

  float **point; 
  float *value; 
  float *dummy; 

  point = new float* [nEof+1];
  for (j=0; j < (nEof+1); j++) {
    point[j] = new float [nEof];
    for(i=0; i < nEof; i++) point[j][i]=0.;     
  }


  /* WARNING */


  for(j=0; j < nEof; j++) point[j][j]=0.1; 

  value = new float [nEof+1]; 
  dummy = new float [nEof]; 

  for(j=0; j < (nEof+1); j++) {
    for(i=0; i < nEof; i++) dummy[i]=point[j][i]; 
    value[j]=costFunction(dummy); 
  }

  simplexND(point, value, nEof); 

}

/* move the surface using EOFs*/

void statistical::moveEof() {

  int i,j,n; 

  for(j=0; j < nVertex; j++) 
    for(i=0; i < 3; i++) {
      vertex[j].r[i]=meanSurface[j][i]; 
      for(n=0; n < nEof; n++) vertex[j].r[i] += coefEof[n]*eof[n][j][i];
    }

}

/* compute cost function of the surface against the edge */ 

float statistical::costFunction(float *dummy) {

  int nIntegral, nTemplate; 

  int j, n; 

  const float coefBorder=0.; // 200

  for(int n=0; n < nEof; n++) coefEof[n]=dummy[n]; 

  moveEof(); 

  float cost=0., costBorder=0.; 

  nIntegral=nTemplateMax; 

  nTemplate=nTemplateMax; 

  /* WARNING */ 

  nIntegral=21; 
  nTemplate=21; 

  localCoordinateSystem();

  computeShift(nIntegral, nTemplate);

  for(j=0; j < nVertex; j++) { 
    cost += vertex[j].displacement*vertex[j].displacement; 
    costBorder += 
      getBorderMap(vertex[j].r[0],vertex[j].r[1],  vertex[j].r[2]); 
  }

  cost += coefBorder*costBorder; 

  for(n=0; n < nEof; n++) cout << coefEof[n] << " "; 
  cout << endl; 
  cout << sqrt(cost/nVertex) << endl; 

  return sqrt(cost/nVertex); 

}

