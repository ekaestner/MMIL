
/* Author: Igor A. Podgorny (e-mail: igor@fiji.ucsd.edu)  
   
   Project - "Generating tri-files for FreeSurfer" 
   Method of Solution - Subdividing polygons approximating a sphere  

   The code 
   - reads a polygonal approximation of n-th order 
   - generates and writes a polygonal approximation of (n+1)-th order 
   - displays the results using OpenGL 

     Order        File       Vertices    Triangles  
    
       0         ic0.tri        12          20   
       1         ic1.tri        42          80  
       2         ic2.tri       162         320  
       3         ic3.tri       642        1280  
       4         ic4.tri      2562        5120  
       5         ic5.tri     10242       20480  
       6         ic6.tri     40962       81920  
       7         ic7.tri    163842      327680

   Zero-order tri-file: 

   12
    1   .0000   .0000  1.0000
    2   .8944   .0000   .4472
    3   .2764   .8507   .4472
    4  -.7236   .5257   .4472
    5  -.7236  -.5257   .4472
    6   .2764  -.8507   .4472
    7   .7236  -.5257  -.4472
    8   .7236   .5257  -.4472
    9  -.2764   .8507  -.4472
   10  -.8944   .0000  -.4472
   11  -.2764  -.8507  -.4472
   12   .0000   .0000 -1.0000
   20
    1   1   5   4
    2   1   6   5
    3   1   2   6
    4   1   3   2
    5   1   4   3
    6   4   9   3
    7   4  10   9
    8   4   5  10
    9   5  11  10
   10   5   6  11
   11   6   7  11
   12   6   2   7
   13   2   8   7
   14   2   3   8
   15   3   9   8
   16   9  10  12
   17  10  11  12
   18  11   7  12
   19   7   8  12
   20   8   9  12

   Compilation: 

   g++ -O3 triFileGenerator.cpp -lglut -lGL -lGLU -lX11 -lXmu -lXi -lm 
       -L /usr/X11R6/lib/                        */ 


#include <GL/glut.h>
#include <cstdlib>
#include <iostream> 
#include <cmath> 
#include <fstream> 
#include <sstream> 
#include <string> 

using namespace std; 

void normalize(float v[3]) 
{    
   GLfloat d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 

   if (d == 0.0) {
     cout << "Error: zero length vector" << endl;    
     return;
   }

   v[0] /= d; v[1] /= d; v[2] /= d; 
}

void normcrossprod(float v1[3], float v2[3], float out[3]) 
{ 
   GLint i, j; 
   GLfloat length;

   out[0] = v1[1]*v2[2] - v1[2]*v2[1]; 
   out[1] = v1[2]*v2[0] - v1[0]*v2[2]; 
   out[2] = v1[0]*v2[1] - v1[1]*v2[0]; 
   normalize(out); 
}

void init(void) 
{
   GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
   GLfloat mat_shininess[] = { 50.0 };
   GLfloat light_position[] = { 0.0, 0.5, -1.0, 0.0 };

   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel (GL_SMOOTH);

   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);
}
 
void subdivide(GLfloat **vdata0, GLuint **tindices0, 
	       GLfloat **vdata, GLuint **tindices, 
               int nVertex0, int nFace0) 
{
  GLuint jCount, jMatch, jFace, i, j, k;
  GLuint j0, j1, j2; 
  GLuint numberNew[3]; 

  GLfloat coorNew[3][3], coorV[3];

  jCount=nVertex0-1; 

  for (j = 0; j < nVertex0; j++) 
    for (i = 0; i < 3; i++) vdata[j][i]=vdata0[j][i]; 
  for (j=0; j < nFace0; j++) 
    for (i=0; i < 3; i++) tindices[j][i]=tindices0[j][i]; 

  for (jFace=0; jFace < nFace0; jFace++) {

    /* step 1: compute new vertices on the sphere */ 

	j0=tindices0[jFace][0]; 
	j1=tindices0[jFace][1]; 
	j2=tindices0[jFace][2]; 

      for(i=0; i < 3; i++) {
	coorNew[0][i]=(vdata0[j0][i]+vdata0[j1][i])/2.; 
	coorNew[1][i]=(vdata0[j1][i]+vdata0[j2][i])/2.; 
	coorNew[2][i]=(vdata0[j0][i]+vdata0[j2][i])/2.; 
      }

      for(i=0; i < 3; i++) {
	for(k=0; k < 3; k++) coorV[k]=coorNew[i][k];
	normalize(coorV);
	for(k=0; k < 3; k++) coorNew[i][k]=coorV[k];	
      }

      /* step 2: check if the vertices are listed and add them if not */ 

      for(i=0; i < 3; i++) {

	jMatch=0; 
	for(j=0; j <= jCount; j++) {
	  if(coorNew[i][0]==vdata[j][0] && 
	     coorNew[i][1]==vdata[j][1] && 
	     coorNew[i][2]==vdata[j][2] ) {
	    jMatch=1; 
	    numberNew[i]=j; 
	  } 
	}

	if(jMatch==0) {
	  jCount++; 
	  numberNew[i]=jCount; 
	  for(k=0; k < 3; k++) vdata[jCount][k]=coorNew[i][k]; 	    
	}	            
      }      

      /* step 3: add new triangles and indices to the list */ 

      tindices[4*jFace][0] = j0; 
      tindices[4*jFace][1] = numberNew[0]; 
      tindices[4*jFace][2] = numberNew[2];

      tindices[4*jFace+1][0] = numberNew[0]; 
      tindices[4*jFace+1][1] = j1; 
      tindices[4*jFace+1][2] = numberNew[1];

      tindices[4*jFace+2][0] = numberNew[2]; 
      tindices[4*jFace+2][1] = numberNew[1];
      tindices[4*jFace+2][2] = j2;

      tindices[4*jFace+3][0] = numberNew[0]; 
      tindices[4*jFace+3][1] = numberNew[1];
      tindices[4*jFace+3][2] = numberNew[2]; 

  }

} 

void display(void)
{
  string sDummy; 
  GLuint faceNumber[8] = {20, 80, 320, 1280, 5120, 20480, 81920, 327680}; 

  GLuint inputNumber=4, outputNumber=inputNumber+1;  

  GLuint nFace0=faceNumber[inputNumber], nFace=4*nFace0; 
  GLuint nVertex0=nFace0/2+2, nVertex=nFace/2+2; 

  GLuint nDummy, i, j, k; 

  GLfloat **vdata0, **vdata; 
  GLuint **tindices0, **tindices; 

  /* allocate memory and read the data */

  vdata0 = new GLfloat* [nVertex0]; 
  for(i=0; i < nVertex0; i++) vdata0[i] = new GLfloat [3]; 
  tindices0 = new GLuint* [nFace0]; 
  for(j=0; j < nFace0; j++) tindices0[j] = new GLuint [3];  

  vdata = new GLfloat* [nVertex]; 
  for(i=0; i < nVertex; i++) vdata[i] = new GLfloat [3]; 
  tindices = new GLuint* [nFace]; 
  for(j=0; j < nFace; j++) tindices[j] = new GLuint [3];  

  stringstream converterInput; 
  converterInput << inputNumber; 
  sDummy = "ic" + converterInput.str() + "_new.tri"; 
  ifstream triFileInput(sDummy.c_str());

  triFileInput >> nDummy;
  for (j=0; j < nVertex0; j++) {
    triFileInput >> nDummy; 
    for(i=0; i < 3; i++) triFileInput >> vdata0[j][i];
  }  

  triFileInput >> nDummy;
  for(k=0; k < nFace0; k++) { 
    triFileInput >> nDummy; 
    for(i=0; i < 3; i++) {
      triFileInput >> tindices0[k][i];  
      tindices0[k][i]--;
    }
  }

  triFileInput.close();


  GLfloat d1[3], d2[3], norm[3]; 

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBegin(GL_TRIANGLES);    

  /* generate polygonal approximation of (n+1)-th order */ 

  subdivide(vdata0, tindices0, vdata, tindices, nVertex0, nFace0);  

  stringstream converterOutput; 
  converterOutput << outputNumber; 
  sDummy = "ic" + converterOutput.str() + "_new.tri"; 
  ofstream triFileOutput(sDummy.c_str());

  triFileOutput << nVertex << endl;
  for (j=0; j < nVertex; j++) {
    triFileOutput << (j+1) << " "; 
    for(i=0; i < 3; i++) triFileOutput << vdata[j][i] << " ";
    triFileOutput << endl;   
  }  

  triFileOutput << nFace << endl; 
  for(k=0; k < nFace; k++) { 
    triFileOutput << (k+1) << " "; 
    for(i=0; i < 3; i++) triFileOutput << (tindices[k][i]+1) << " ";  
    triFileOutput << endl; 
  }

  triFileOutput.close();

  for (int i = 0; i < nFace; i++) {    

    for (int j = 0; j < 3; j++) {  
      d1[j] = vdata[tindices[i][0]][j] - vdata[tindices[i][1]][j];    
      d2[j] = vdata[tindices[i][1]][j] - vdata[tindices[i][2]][j];    
    }

    normcrossprod(d1, d2, norm); 

    glNormal3fv(norm);

    glVertex3fv(&vdata[tindices[i][0]][0]); 
    glVertex3fv(&vdata[tindices[i][1]][0]); 
    glVertex3fv(&vdata[tindices[i][2]][0]); 
  }

  glEnd();

  glFlush ();
}

void reshape (int w, int h)
{
   glViewport (0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode (GL_PROJECTION);
   glLoadIdentity();
   if (w <= h)
      glOrtho (-1.5, 1.5, -1.5*(GLfloat)h/(GLfloat)w,
         1.5*(GLfloat)h/(GLfloat)w, -10.0, 10.0);
   else
      glOrtho (-1.5*(GLfloat)w/(GLfloat)h,
         1.5*(GLfloat)w/(GLfloat)h, -1.5, 1.5, -10.0, 10.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
}

int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
   glutInitWindowSize (500, 500); 
   glutInitWindowPosition (100, 100);
   glutCreateWindow (argv[0]);
   init ();
   glutDisplayFunc(display); 
   glutReshapeFunc(reshape);
   glutMainLoop();
   return 0;
}
