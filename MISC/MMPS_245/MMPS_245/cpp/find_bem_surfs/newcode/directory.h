
class directory {

 public: 

  directory(); 

  ~directory(); 

  int surfaceNumber; 

  string *surfaceName; 

  string mriInputString; // ??  

  string getTriInputDirectory(void); 

  string getMriInputDirectory(void); 

  string getSurfaceName(int); 

  string mriFileName; // ?? 

 protected: 

  string mriInputDirectory; 

  string triInputDirectory;  

  string subjectName;

  string projectDirectory;

  string methodSubdirectory;

  string surfaceInputDirectory;

  string meanSurfaceInputDirectory;

  string templateInputDirectory; 

  string eofInputDirectory; 

  string talairachInputDirectory; 

  void getMriFileInfo(void); 

};
