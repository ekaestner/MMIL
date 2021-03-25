
#ifndef PD_DISK_COPY
#define PD_DISK_COPY
#undef  PD_DISK_COPY
#endif

class mri {
  
public:

  mri(string, string);
  
  virtual ~mri();

protected: 

  /* PD refers to proton density MRI imaging */ 

  string mriInputFileName, mriOutputFileName;

  float ***dataMghVolume;

  float mghIntensityScaling; 

  int nByte, nRow, nSlice;

  float voxelVolume; 

  /* Matrix dimension */

  int nRas; 

  /* fMVox2Ras - refers to voxel_to_RAS transform and vice versa */ 

  float **fMVox2Ras, **fMRas2Vox; 

  float fIndexCenter[4], rasCenter[4]; 

  void readMghFile();  

  float getIntensity(float, float, float); 

  int checkIfLittleEndian(void); 

  short shortSwap(short); 

  int intSwap(int);

  float floatSwap(float);

  float computeDeterminant(float**, int); 

  void computeInverseMatrix(float**, int, float**);

};


