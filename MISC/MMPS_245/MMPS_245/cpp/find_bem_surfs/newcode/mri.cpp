 
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <unistd.h> // check for an alternative name 

using namespace std;

#include "mri.h"

/* mri constructor */ 

mri::mri(const string &mriInputFileName){

  if(checkIfLittleEndian() == 0) {

    cout << "Problem with byte swapping in mri::mri" << endl; 
    abort(); 

  }  

  /* WARNING - remove writing */ 

  mriOutputFileName="/home/igor/slice";

  mghIntensityScaling=1.;
  nRas=4;  

  readMghFile(mriInputFileName);   

}

/* mri destructor */

mri::~mri(){

  int jByte, jRow; 

  /*Delete memory for Vox_to_Ras and Ras_to_Vox matrices */ 

  for (int i=0; i < nRas; i++) {
    delete [] fMVox2Ras[i];
    delete [] fMRas2Vox[i];
  }

  delete [] fMVox2Ras;
  delete [] fMRas2Vox; 

  /* Delete memory for MGH arrays */

  for (jByte = 0; jByte < nByte; jByte++)
  {
    for (jRow = 0; jRow < nRow; jRow++) delete [] dataMghVolume[jByte][jRow];
    delete [] dataMghVolume[jByte];
  }
  delete [] dataMghVolume;

}

/* Read MRI PD data and re-scale data by a predefined factor  
 *
 * Future modification - implement automatic treatment of byte-ordering 
 * currently the code assumes little-endian byte ordering for input */  

void mri::readMghFile(const string &mriInputFileName) {

  const int nHeader1=7, nHeader2=15, nHeader3=48;

  enum {MRI_UCHAR=0, MRI_INT, MRI_LONG, MRI_FLOAT, MRI_SHORT, MRI_BITMAP};  
  
  int shift=2;
  int mriType; 
  int header1[nHeader1];

  int jByte, jRow, jSlice;
  int k, i; 

  float header2[nHeader2];

  float ci, cj, ck;

  float xsize, ysize, zsize; 
  float x_r, x_a, x_s;
  float y_r, y_a, y_s;
  float z_r, z_a, z_s;
  float c_r, c_a, c_s;

  unsigned char* mghSliceUchar; 
  short* mghSliceShort; 
  float* mghSliceFloat;

  ifstream inputFile(mriInputFileName.c_str());

  if(!inputFile.is_open()) {
    cout << "Unable to open MGH file in function mri::readMghFile" << endl; 
    abort(); 
  }


  inputFile.read ( (char *) header1, 4*nHeader1);

  nByte=intSwap(header1[1]);
  nRow=intSwap(header1[2]);
  nSlice=intSwap(header1[3]);

  mriType=intSwap(header1[5]);  

  switch(mriType) {

  case MRI_UCHAR: 
    mghSliceUchar = new unsigned char[nByte*nRow];       
    break; 

  case MRI_FLOAT: 
    mghSliceFloat = new float[nByte*nRow];       
    break; 

  case MRI_SHORT: 
    mghSliceShort = new short[nByte*nRow]; 
    break; 

  default:
    cout << "Wrong MRI type in MGH file in function mri::readMghFile" << endl;
    abort(); 

  }

 /* Allocate memory for RAS-to-VOX and VOX-to-RAS matrices */ 

 
  fMVox2Ras = new float* [nRas];
  fMRas2Vox = new float* [nRas];

  for (int i=0; i < nRas; i++) {
    fMVox2Ras[i] = new float [nRas];
    fMRas2Vox[i] = new float [nRas]; 
  } 

    /* Allocate memory for 3D arrays */

  dataMghVolume = new float** [nByte];
  for (int jByte = 0; jByte < nByte; jByte++) {

    dataMghVolume[jByte] = new float* [nRow];

    for (int jRow = 0; jRow < nRow; jRow++) {
      dataMghVolume[jByte][jRow] = new float [nSlice];
    }
  }

  
  inputFile.seekg(inputFile.tellg() + long(shift));
  inputFile.read ( (char *) header2, 4*nHeader2);

  xsize=floatSwap(header2[0]);
  ysize=floatSwap(header2[1]);
  zsize=floatSwap(header2[2]);

  voxelVolume=xsize*ysize*zsize; 

  x_r=floatSwap(header2[3]);
  x_a=floatSwap(header2[4]);
  x_s=floatSwap(header2[5]);

  y_r=floatSwap(header2[6]);
  y_a=floatSwap(header2[7]);
  y_s=floatSwap(header2[8]);

  z_r=floatSwap(header2[9]);
  z_a=floatSwap(header2[10]);
  z_s=floatSwap(header2[11]);

  c_r=floatSwap(header2[12]);
  c_a=floatSwap(header2[13]);
  c_s=floatSwap(header2[14]);

  /* Compute RAS_to_voxel matrix  */ 

  fMVox2Ras[0][0] = xsize * x_r ;
  fMVox2Ras[0][1] = ysize * y_r ;
  fMVox2Ras[0][2] = zsize * z_r ;

  fMVox2Ras[1][0] = xsize * x_a ;
  fMVox2Ras[1][1] = ysize * y_a ;
  fMVox2Ras[1][2] = zsize * z_a ;
  
  fMVox2Ras[2][0] = xsize * x_s ;
  fMVox2Ras[2][1] = ysize * y_s ;
  fMVox2Ras[2][2] = zsize * z_s ;

  fMVox2Ras[3][0] = 0.;
  fMVox2Ras[3][1] = 0.;
  fMVox2Ras[3][2] = 0.;

  ci = (float(nByte)-1.)/2.;
  cj = (float(nRow)-1.)/2.;
  ck = (float(nSlice)-1.)/2.;

  fMVox2Ras[0][3] = c_r - 
    (fMVox2Ras[0][0]*ci + fMVox2Ras[0][1]*cj + fMVox2Ras[0][2]*ck);
  fMVox2Ras[1][3] = c_a - 
    (fMVox2Ras[1][0]*ci + fMVox2Ras[1][1]*cj + fMVox2Ras[1][2]*ck);
  fMVox2Ras[2][3] = c_s - 
    (fMVox2Ras[2][0]*ci + fMVox2Ras[2][1]*cj + fMVox2Ras[2][2]*ck);
  fMVox2Ras[3][3] = 1.;

  /* Compute central Index and RAS coordinates */
  /* Future modification - address half-a-voxel issue */ 

  fIndexCenter[0]=float(nByte)/2.+1.; 
  fIndexCenter[1]=float(nRow)/2.+1.; 
  fIndexCenter[2]=float(nSlice)/2.+1.; 
  fIndexCenter[3]=1.; 

  for(k=0; k < nRas; k++) {
    rasCenter[k]=0.;
    for(i=0; i < nRas; i++) {
      rasCenter[k] += fMVox2Ras[k][i]*fIndexCenter[i];
    }
  }

  /* Compute voxel_to_RAS matrix  */ 

  computeInverseMatrix(fMVox2Ras, nRas, fMRas2Vox); 

  inputFile.seekg(inputFile.tellg() + long(shift+4*nHeader3)); 


  for(jSlice=0; jSlice < nSlice; jSlice++) { 


    switch(mriType) {

    case MRI_UCHAR: 
      inputFile.read ( (char*) mghSliceUchar, nByte*nRow);
      break; 

    case MRI_FLOAT: 
      inputFile.read ( (char *) mghSliceFloat, 4*nByte*nRow);
      break; 

    case MRI_SHORT: 
      inputFile.read ( (char *) mghSliceShort, 2*nByte*nRow);
      break; 

    default:
      cout << "Wrong MRI type in mgh-file" << endl; 

    }

    for(jByte=0; jByte < nByte; jByte++) {
      for(jRow=0; jRow < nRow; jRow++) {
	k=nRow*jByte+jRow;

	switch(mriType) {

	case MRI_UCHAR: 
	  dataMghVolume[jRow][jByte][jSlice]=float(int(mghSliceUchar [k])); 
	  break; 

	case MRI_FLOAT: 
	  dataMghVolume[jRow][jByte][jSlice]=floatSwap(mghSliceFloat [k]);
	  break; 

	case MRI_SHORT: 
	  dataMghVolume[jRow][jByte][jSlice]=shortSwap(mghSliceShort[k]);
	  break; 

	default:
	  cout << "Wrong MRI type mgh-file" << endl; 

	}
      }
    }
  }

  for(jSlice=0; jSlice < nSlice; jSlice++) 
    for(jRow=0; jRow < nRow; jRow++)
      for(jByte=0; jByte < nByte; jByte++) 
      	dataMghVolume[jRow][jByte][jSlice] *= mghIntensityScaling; 


  /* WARNING - temporary output

  print2dSlice(dataMghVolume, "slice", 1); 
  print2dSlice(dataMghVolume, "slice", 2); 
  print2dSlice(dataMghVolume, "slice", 3); 

  */

#ifdef PD_DISK_COPY

  ofstream outputFile(mriOutputFileName.c_str());

  for(jSlice=0; jSlice < nSlice; jSlice++)  
    for(jRow=0; jRow < nRow; jRow++) 
      for(jByte=0; jByte < nByte; jByte++) 
	outputFile << dataMghVolume[jRow][jByte][jSlice] << " "; 
      
  outputFile << endl; 
    
  outputFile.close();

#endif /* PD_DISK_COPY */

  inputFile.close(); 

  switch(mriType) {

  case MRI_UCHAR: 
    delete [] mghSliceUchar; 
    break; 

  case MRI_FLOAT: 
    delete [] mghSliceFloat; 
    break; 

  case MRI_SHORT: 
    delete [] mghSliceShort; 
    break; 

  default:
    cout << "Wrong MRI type in mgh-file" << endl; 

  }

}

/* Get image intensity at a given point with RAS coordinates 
 * 
 * Future modification - interpolate intensity using neighbor voxels */    

float mri::getIntensity(float r, float a, float s) {

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
      (jVoxel[2] < nSlice) && (jVoxel[2] >= 0) ) 

    return dataMghVolume[jVoxel[0]][jVoxel[1]][jVoxel[2]];
 

  /* image intensity is zero outside of the volume */ 
                                                                               
  else {
                                                                               
    return 0.;
  }

  cout << "problem in getIntensity subroutine" << endl;

}

/* WARNING - merge both functions together */ 

/* Set image intensity at a given point with RAS coordinates 
 * 
 * Future modification - interpolate intensity using neighbor voxels */    

int mri::setIntensity(float fIntensity, float r, float a, float s) {

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

    dataMghVolume[jVoxel[0]][jVoxel[1]][jVoxel[2]] = fIntensity;

    cout << jVoxel[0] << " " << jVoxel[1] << " " << jVoxel[2] << " " 
	 << dataMghVolume[jVoxel[0]][jVoxel[1]][jVoxel[2]] << endl; 

    return 1; 
  }

  /* no voxel for a given point */ 
                                                                               
  else {
                                                                               
    return 0;
  }

  cout << "problem in setIntensity subroutine" << endl;

}


/* Check if the architecture uses Little Endian */

int mri::checkIfLittleEndian()
{

  int test=1; 
  bool *p=  (bool*) &test; 

  return int(p[0]);  

}

/* Swap short integer bytes */ 

short mri::shortSwap(short s)
{
  unsigned char b1, b2;
  
  b1 = s & 255;
  b2 = (s >> 8) & 255;

  return (b1 << 8) + b2;
}


/* Swap integer bytes */ 

int mri::intSwap (int i)
{
  unsigned char b1, b2, b3, b4;
                                                                               
  b1 = i & 255;
  b2 = ( i >> 8 ) & 255;
  b3 = ( i >> 16 ) & 255;
  b4 = ( i >> 24 ) & 255;
                                                                               
  return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
}

/* Swap float bytes */ 

float mri::floatSwap(float f)
{
  union { float f; unsigned char b[4]; } dat1, dat2;
                                                                               
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
                                                                               
  return dat2.f;
}

/* Recursive definition of determinate using expansion by minors 
 *
 * Future modification - implement linear algebra routines using 
 * Inter Math Kernel Library  
 * see http://www.intel.com/software/products/mkl/ */ 


float mri::computeDeterminant(float **a,int n)
{
   int i,j,j1,j2;
   float det = 0;
   float **m = NULL;
                                                                               
   if (n < 1) { /* Error */
                                                                               
   } else if (n == 1) { 
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = new float* [4*(n-1)];
         for (i=0;i<n-1;i++)
            m[i] = new float [4*(n-1)];
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * computeDeterminant(m,n-1);
         for (i=0;i<n-1;i++)
            delete [] m[i];
         delete [] m;
      }
   }
   return(det);
}

/* Compute inverse of a square matrix */

void mri::computeInverseMatrix(float **a,int n, float **b)
{
   int i,j,ii,jj,i1,j1;
   float detC, detA;
   float **c;

   detA=computeDeterminant(a,n);
                                                                               
   c = new float* [4*(n-1)];
   for (i=0;i<n-1;i++)
     c[i] = new float [4*(n-1)];
                                                                               
   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {
                                                                               
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }
                                                                               
         detC = computeDeterminant(c,n-1);
                                                                               
	 b[j][i] = pow(-1.0,i+j+2.0)*detC/detA;
      }
   }
   for (i=0;i<n-1;i++)
      delete [] c[i];
   delete [] c;
}

/* Auxilary function - get path to working directory */

string mri::getWorkingDirectory(void) 
{
  const int maxSize=255; // WARNING - may not be the case
  char dummy[maxSize];

  getcwd(dummy, maxSize);

  return dummy;
}
                                                                               
/* Auxilary function - print a 2D slice in the middle of 3D volume */

void mri::print2dSlice(float*** array3d, string sliceName, int sliceNumber)
{

  int jByte, jRow, jSlice; 

  string workingDirectory = getWorkingDirectory();

  string outputFileName; 

  switch(sliceNumber) {

  case 1: 
    { 
      outputFileName = workingDirectory + "/" + sliceName + "1"; 
      ofstream file1(outputFileName.c_str());

      jByte = nByte/2; 
      for (jRow = 0; jRow < nRow; jRow++) {
	for (jSlice = 0; jSlice < nSlice; jSlice++) 
	  file1 << array3d[jByte][jRow][jSlice] << " "; 
	file1 << endl; 
      }
      file1.close();
    }
    break; 

  case 2: 
    { 
      outputFileName = workingDirectory + "/" + sliceName + "2"; 
      ofstream file2(outputFileName.c_str());

      jRow = nRow/2; 
      for (jByte = 0; jByte < nByte; jByte++) {
	for (jSlice = 0; jSlice < nSlice; jSlice++) 
	  file2 << array3d[jByte][jRow][jSlice] << " "; 
	file2 << endl; 
      }
      file2.close();
    }
    break; 

  case 3: 
    { 
      outputFileName = workingDirectory + "/" + sliceName + "3"; 
      ofstream file3(outputFileName.c_str());

      jSlice = nSlice/2; 
      for (jByte = 0; jByte < nByte; jByte++) {
	for (jRow = 0; jRow < nRow; jRow++) 
	  file3 << array3d[jByte][jRow][jSlice] << " "; 
	file3 << endl; 
      }
      file3.close();
    }
    break; 

  default: 
    {
      cout << "Problem in mri::print2dSlice" << endl; 
      abort(); 
    }

  }

}
