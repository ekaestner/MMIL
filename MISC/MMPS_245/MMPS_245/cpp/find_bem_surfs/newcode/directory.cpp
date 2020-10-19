
#include <string>
#include <iostream> 
#include <string> 

using namespace std;

#include "directory.h"

/* directory constructor */ 

directory::directory() {

  /* list all pre-defined directories here */ 

  // triInputDirectory = "/home/igor/TRIANGLE/NEW_TRI_FILES/"; 

  projectDirectory = "/space/monkeys/1/home/PROJECTS/VETSA/";

  methodSubdirectory = "DETERMINISTIC/"; 

  surfaceInputDirectory = projectDirectory + methodSubdirectory + "TRI/";

  meanSurfaceInputDirectory = projectDirectory + "MEAN/";

  templateInputDirectory = projectDirectory + "MEAN/";

  eofInputDirectory = projectDirectory + "MEAN/";

  talairachInputDirectory = projectDirectory + "XFM/";


  /* process input info here */ 

  surfaceNumber=3; 

  surfaceName = new string [surfaceNumber]; 

  surfaceName[0]="_inner_skull"; 
  surfaceName[1]="_outer_skull"; 
  surfaceName[2]="_outer_scalp"; 

  cout << "Reading MRI file name" << endl; 
  cin >> mriInputString; 

  getMriFileInfo(); 

  cout << "Getting TRI file directory" << endl; 
  cin >> triInputDirectory; 

}

/* directory destructor */ 

directory::~directory() {


  delete [] surfaceName;  


}

/* grt input MRI directory */ 


string directory::getMriInputDirectory() {

  return mriInputDirectory; 

}


/* WARNING - check of the tri-file here */ 


string directory::getTriInputDirectory() {

  return triInputDirectory; 

}

/* get MRI file name, subject name and extension */ 

void directory::getMriFileInfo() {

  int length=mriInputString.length(); 
  int pivot=mriInputString.rfind("/", length); 

  mriInputDirectory = mriInputString.substr(0, pivot+1);

  mriFileName = mriInputString.substr(pivot+1, length);

  subjectName = mriInputString.substr(pivot+1, mriFileName.length()-4 ); 
}

/* get surface name */ 

string directory::getSurfaceName(int index) {

  if( (index==0) || (index==1) || (index==2) ) 
    return ( subjectName + surfaceName[index]); 

  else {

    cout << "Problem in directory::getSurfaceName" << endl; 
    abort(); 

  }

}
