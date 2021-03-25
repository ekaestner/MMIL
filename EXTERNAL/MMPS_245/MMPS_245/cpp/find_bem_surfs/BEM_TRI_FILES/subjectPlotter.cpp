
/* same as triFilePlotter, but now includes Talairach transform and 
   uses VETSA subject number as input */ 

#include <GL/glut.h>
#include <cstdlib>
#include <iostream> 
#include <cmath> 
#include <fstream> 
#include <sstream> 
#include <string> 

using namespace std; 

/* define global variables here */ 

string surfacePath="/space/monkeys/1/home/igor/PROJECTS/VETSA_MGH15T/EOF/TRI/"; 
string talairachPath="/space/monkeys/1/home/igor/PROJECTS/VETSA_MGH15T/XFM/"; 
string subjectName;  
int nView;
int nSurface;  

float factor=85.; 
float shift[3][3] = { {0., 0.4, 0.}, {0., 0.4, 0.}, {0., 0., 0.} };  

/* no more global variables below */


void talairachMatrix(float** talairach)
{

  /* open and read Talairach file, but keep only rotational 
     components of the matrix */ 

  string sDummy = talairachPath + subjectName + ".xfm";

  ifstream talairachFileInput(sDummy.c_str());

  for(int j=0; j < 3; j++) 
    for(int i=0; i < 4; i++) talairachFileInput >> talairach[j][i];  

  talairachFileInput.close();

  talairach[0][3]=0.;
  talairach[1][3]=0.;
  talairach[2][3]=0.;

  talairach[3][0]=0.;
  talairach[3][1]=0.;
  talairach[3][2]=0.;

  talairach[3][3]=1.;

}

void rotationMatrix(float** r) 
{
  float r01[3][3]= { {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.} };
  /* top, nose up */ 

  float r02[3][3]= { {1., 0., 0.}, {0., -1., 0.}, {0., 0., -1.} };
  /* bottom, nose down */ 

  float r03[3][3]= { {1., 0., 0.}, {0., 0., 1.}, {0., -1., 0.} };
  /* back, legs down */ 

  float r04[3][3]= { {1., 0., 0.}, {0., 0., -1.}, {0., 1., 0.} };
  /* front, legs down */ 

  float r05[3][3] = { {-1., 0., 0.}, {0., 1., 0.}, {0., 0., -1.} }; 
  /* bottom, nose up */ 

  float r06[3][3] = { {-1., 0., 0.}, {0., -1., 0.}, {0., 0., 1.} }; // #3  
  /* top, nose down */ 

  float r07[3][3] = { {-1., 0., 0.}, {0., 0., 1.}, {0., 1., 0.} }; // #2  
  /* front, legs down */ 

  float r08[3][3] = { {-1., 0., 0.}, {0., 0., -1.}, {0., -1., 0.} };   
  /* back, legs up */ 

  float r09[3][3] = { {0., 1., 0.}, {1., 0., 0.}, {0., 0., -1.} };   
  /* bottom, nose right */ 

  float r10[3][3] = { {0., -1., 0.}, {1., 0., 0.}, {0., 0., 1.} };   
  /* top, nose left */ 

  float r11[3][3] = { {0., 0., 1.}, {1., 0., 0.}, {0., 0., 1.} };   
  /* front, legs left */ 

  float r12[3][3] = { {0., 0., -1.}, {1., 0., 0.}, {0., -1., 0.} };   
  /* back, legs right */ 

  float r13[3][3] = { {0., 1., 0.}, {-1., 0., 0.}, {0., 0., 1.} };   
  /* top, nose right */ 

  float r14[3][3] = { {0., -1., 0.}, {-1., 0., 0.}, {0., 0., -1.} };   
  /* bottom, nose left */ 

  float r15[3][3] = { {0., 0., 1.}, {-1., 0., 0.}, {0., -1., 0.} };   
  /* back, legs left */ 

  float r16[3][3] = { {0., 0., -1.}, {-1., 0., 0.}, {0., 1., 0.} };   
  /* face, legs right */ 

  float r17[3][3] = { {0., 0., -1.}, {0., 1., 0.}, {1., 0., 0.} };   
  /* right, nose up */ 

  float r18[3][3] = { {0., 0., 1.}, {0., -1., 0.}, {1., 0., 0.} };   
  /* right, nose down */ 

  float r19[3][3] = { {0., 1., 0.}, {0., 0., 1.}, {1., 0., 0.} };   
  /* right, nose right */ 

  float r20[3][3] = { {0., -1., 0.}, {0., 0., -1.}, {1., 0., 0.} };   
  /* right, legs up */ 

  float r21[3][3] = { {0., 0., 1.}, {0., -1., 0.}, {-1., 0., 0.} };   
  /* left, nose up */ 

  float r22[3][3] = { {0., 0., -1.}, {0., -1., 0.}, {-1., 0., 0.} };   
  /* left, nose down */ 

  float r23[3][3] = { {0., -1., 0.}, {0., 0., 1.}, {-1., 0., 0.} }; // #1    
  /* left, nose left */ 

  float r24[3][3] = { {0., 1., 0.}, {0., 0., -1.}, {-1., 0., 0.} };     
  /* left, legs up */ 

  switch(nView) {

  case 1: 

    for(int j=0; j < 3; j++) 
      for(int i=0; i < 3; i++) r[j][i]=r23[j][i];     
    break; 

  case 2: 

    for(int j=0; j < 3; j++) 
      for(int i=0; i < 3; i++) r[j][i]=r07[j][i];     
    break; 

  case 3: 

    for(int j=0; j < 3; j++) 
      for(int i=0; i < 3; i++) r[j][i]=r06[j][i];     
    break; 

  default: 
    cout << "Wrong projection number" << endl; 

  }
}


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
 
void display(void)
{
  string surfaceName; 

  switch(nSurface) {

  case 1: 

    surfaceName="_inner_skull5"; 
    break; 

  case 2: 

    surfaceName="_outer_skull5"; 
    break; 

  case 3: 

    surfaceName="_outer_scalp5"; 
    break; 

  default: 
    cout << "Wrong surface number" << endl; 

  }

  string sDummy = surfacePath + subjectName + surfaceName + ".tri"; 

  GLuint nVertex, nFace; 
  GLuint i, j, k, nDummy; 

  GLfloat **vdata; 
  GLuint **tindices; 

  GLfloat v[4]; 

  ifstream triFileInput(sDummy.c_str());

  triFileInput >> nVertex;

  vdata = new GLfloat* [nVertex]; 
  for(i=0; i < nVertex; i++) vdata[i] = new GLfloat [3]; 

  for (j=0; j < nVertex; j++) {
    triFileInput >> nDummy; 
    for(i=0; i < 3; i++) {
      triFileInput >> vdata[j][i];
      vdata[j][i] /= factor;      
    }
  }  

  triFileInput >> nFace;

  tindices = new GLuint* [nFace]; 
  for(j=0; j < nFace; j++) tindices[j] = new GLuint [3];  

  for(k=0; k < nFace; k++) { 
    triFileInput >> nDummy; 
    for(i=0; i < 3; i++) {
      triFileInput >> tindices[k][i];  
      tindices[k][i]--;
    }
  }

  triFileInput.close();


  /* align the object with the Talairach system */ 

  float** talairach; 
  talairach = new float*[4]; 
  for(i=0; i < 4; i++) talairach[i] = new float [4];   
  talairachMatrix(talairach);  

  v[3]=1.; 

  for (k=0; k < nVertex; k++) {

    for(i=0; i < 3; i++) v[i]=vdata[k][i];

    for(j=0; j < 3; j++) {
      vdata[k][j]=0.; 
      for(i=0; i < 4; i++) vdata[k][j] += talairach[j][i]*v[i]; 
    }	
  }

  /* rotate the object */ 

  float** r; 
  r = new float*[3]; 
  for(i=0; i < 3; i++) r[i] = new float [3];   
  rotationMatrix(r);   

  for (j=0; j < nVertex; j++) {

    for(i=0; i < 3; i++) v[i]=vdata[j][i]; 

    for(k=0; k < 3; k++) {
      vdata[j][k]=0.; 
      for(i=0; i < 3; i++) vdata[j][k] += r[k][i]*v[i]; 
    }
  }

  /* move the object upwards */


    for (j=0; j < nVertex; j++) 
      for(i=0; i < 3; i++) vdata[j][i] += shift[nView-1][i]; 


  GLfloat d1[3], d2[3], norm[3]; 

  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glBegin(GL_TRIANGLES);    

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

  if(argc < 2) {
    cout << "SYNOPSIS: *** subjectName" << endl;
    abort(); 
  }

  subjectName=argv[1]; 

  cout << "Projection number? "; 
  cin >> nView; 

  cout << "Inner skull (1), outer skull (2) or scalp (3)? "; 
  cin >> nSurface; 
  
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
