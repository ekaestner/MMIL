
#ifndef PD_DISK_COPY
#define PD_DISK_COPY
#undef  PD_DISK_COPY
#endif

class mri {
  
 public:

  mri(const string&);
  
  virtual ~mri();

 protected: 

  /* PD refers to proton density MRI imaging */ 

  string mriOutputFileName; // REMOVE

  float ***dataMghVolume;

  float mghIntensityScaling; 

  int nByte, nRow, nSlice;

  float voxelVolume; 

  /* Matrix dimension */

  int nRas; 

  /* fMVox2Ras - refers to voxel_to_RAS transform and vice versa */ 

  float **fMVox2Ras, **fMRas2Vox; 

  float fIndexCenter[4], rasCenter[4]; 

  void readMghFile(const string&); 

  string getWorkingDirectory(void); 

  void print2dSlice(float***, string, int); 

  float getIntensity(float, float, float); 

  int setIntensity(float, float, float, float); 

  int checkIfLittleEndian(void); 

  short shortSwap(short); 

  int intSwap(int);

  float floatSwap(float);

  float computeDeterminant(float**, int); 

  void computeInverseMatrix(float**, int, float**);

};


