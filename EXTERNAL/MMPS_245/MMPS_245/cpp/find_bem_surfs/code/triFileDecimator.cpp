
#include <iostream> 
#include <string> 
#include <sstream> 
#include <fstream> 
#include <stdlib.h>

using namespace std; 

int main(int argc, char** argv)
{

  if(argc < 2) {
    cout << "SYNOPSIS: *** fileName" << endl;
    abort(); 
  }

  const int maxMeshNumber=8; 

  struct structTriFileInfoType { 
    int nVertex; 
    int nFace; 
  } *triFileInfo; 

  triFileInfo = new structTriFileInfoType[maxMeshNumber];  

  int nVertex, nFace; 
  int i, j, k, nDummy; 

  triFileInfo[0].nFace=20; 
  for(i=1; i < maxMeshNumber; i++) 
    triFileInfo[i].nFace = 4*triFileInfo[i-1].nFace; 
  for(i=0; i < maxMeshNumber; i++) 
    triFileInfo[i].nVertex = triFileInfo[i].nFace/2 + 2;

  float **vdata; 
  int **tindices; 

//  string triInputDirectory="/home/igor/TRIANGLE/NEW_TRI_FILES/"; 
  string triInputDirectory="/home/halgdev/misc/BEM_TRI_FILES/"; 

  string fileInputString=argv[1]; 

  int fileInputLength=fileInputString.length();  
  int fileInputPivot=fileInputString.rfind("/", fileInputLength); 

  string fileInputDirectory = fileInputString.substr(0, fileInputPivot+1);  
  string fileName = fileInputString.substr(fileInputPivot+1, fileInputLength);
  string fileOutputDirectory=fileInputDirectory;

  int fileTypeLength = fileName.length()-4; 
  string fileType = fileInputString.substr(fileInputPivot+1,fileTypeLength); 

  string meshNumberString 
    = fileInputString.substr(fileInputLength-5,1); 

  istringstream meshStream(meshNumberString);

  int meshNumber; 

  meshStream >> meshNumber; 

  meshNumber--; 

  ifstream triFile;
  stringstream converter; 
  converter << meshNumber;
  meshNumberString=converter.str(); 
  string sDummy = triInputDirectory + "ic" + meshNumberString + ".tri"; 


  cout << endl; 
  cout << "ICOSAHEDRON FILE INFO" << endl;   
  cout << "Name: " << sDummy << endl; 
  cout << "Order: " << meshNumber << endl; 


  triFile.open(sDummy.c_str());

  /* Abort the code if ic-file is unaivailable for reading */ 

  if(!triFile.is_open()) {
    cout << "Unable to open icosahedron tri-file" << endl; 
    abort(); 
  }

  triFile >> nVertex;

  cout << "Number of vertices: " << nVertex << endl;

  if(nVertex != triFileInfo[meshNumber].nVertex) {
    cout << "Wrong number of vertices" << endl; 
    abort(); 
  }


  vdata = new float* [nVertex]; 
  for(i=0; i < nVertex; i++) vdata[i] = new float [3]; 

  for (j=0; j < nVertex; j++) {
    triFile >> nDummy; 
    for(i=0; i < 3; i++) triFile >> vdata[j][i];
  }  

  triFile >> nFace;

  cout << "Number of triangles: " << nFace << endl; 

  if(nFace != triFileInfo[meshNumber].nFace) {
    cout << "Wrong number of triangles" << endl; 
    abort(); 
  }

  tindices = new int* [nFace]; 
  for(j=0; j < nFace; j++) tindices[j] = new int [3];  

  for(k=0; k < nFace; k++) { 
    triFile >> nDummy; 
    for(i=0; i < 3; i++) triFile >> tindices[k][i];  
  }

  triFile.close();


  ifstream inputFile;
  inputFile.open(fileInputString.c_str()); 

  /* Abort the code if ic-file is unaivailable for reading */ 

  if(!inputFile.is_open()) {
    cout << "Unable to open input tri-file" << endl; 
    abort(); 
  }

  cout << endl; 
  cout << "INPUT FILE INFO" << endl; 
  cout << "Name: " << fileInputString << endl; 
  cout << "Order: " << (meshNumber+1) << endl; 

  inputFile >> nDummy;

  cout << "Number of vertices: " << nDummy << endl;

  if(nDummy != triFileInfo[meshNumber+1].nVertex) {
    cout << "Wrong number of vertices" << endl; 
    abort(); 
  }

  cout << "Number of triangles: N/A" << endl; 

  for (j=0; j < nVertex; j++) {
    inputFile >> nDummy; 
    for(i=0; i < 3; i++) inputFile >> vdata[j][i];
  }  

  inputFile.close();


  string outputFileName = fileInputString.substr(0,fileInputLength-5) + 
                          meshNumberString + ".tri"; 


  cout << endl; 
  cout << "OUTPUT FILE INFO" << endl; 
  cout << "Name: " << outputFileName << endl; 
  cout << "Order: " << meshNumber << endl; 
  cout << "Number of vertices: " << nVertex << endl;

  if(nVertex != triFileInfo[meshNumber].nVertex) {
    cout << "Wrong number of vertices" << endl; 
    abort(); 
  }

  cout << "Number of triangles: " << nFace << endl; 

  if(nFace != triFileInfo[meshNumber].nFace) {
    cout << "Wrong number of triangles" << endl; 
    abort(); 
  }

  ofstream outputFile(outputFileName.c_str()); 

  outputFile << nVertex << endl; 

    for(int j=0; j < nVertex; j++) {
      outputFile << (j+1) << " ";
      for(int i=0; i <3; i++) outputFile << vdata[j][i] << " "; 
      outputFile << endl; 
    }

  outputFile << nFace << endl; 

  for(int k=0; k < nFace; k++) {
      outputFile << (k+1) << " "; 
      for(int i=0; i < 3; i++) outputFile << tindices[k][i] << " "; 
      outputFile << endl; 
  }

  outputFile.close(); 

}
