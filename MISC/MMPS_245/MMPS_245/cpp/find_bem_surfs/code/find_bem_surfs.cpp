#include <iostream> 
#include <string> 
#include <fstream> 
#include <sstream> 
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>

using namespace std; 

#include "mri.h"
#include "mesh.h"
#include "deterministic.h"
#include "isotrack.h" 

#define MINARGC 2
#define UNDEFSTR "unknown"
#define STRLEN      1000
#define MATCH(A,B)   (!strcmp(A,B))
#define DEFAULT_TRI_FILE_DIR "/home/mmildev/bem_surfs/BEM_TRI_FILES/"

// variables
static char *progname = NULL;

// parameter defaults
string infile = UNDEFSTR;
string tri_file_dir = DEFAULT_TRI_FILE_DIR;
string outdir = UNDEFSTR;
string outstem = UNDEFSTR;
int ico = 5;  // icosahedral mesh order
int outer_skull_flag = 1;
int outer_scalp_flag = 1;

struct structModelCoefficient {
  int nMove; 
  int nRest; 
  int nRelax;
  float coefTangential; 
  float coefNormal; 
  float coefMRI; 
  float coefRepelling; 
  float intensityThreshold; 

} inner_skull, outer_skull, outer_scalp; 

void defaultCoefs(void) {
  inner_skull.nMove = 300;
  inner_skull.nRest = 30;
  inner_skull.nRelax = 2;
  inner_skull.coefTangential = 0.1;
  inner_skull.coefNormal = 0.1;
  inner_skull.coefMRI = 0.3;
  inner_skull.coefRepelling = 0.0;
  inner_skull.intensityThreshold = 80.0;

  outer_skull.nMove = 100;
  outer_skull.nRest = 30;
  outer_skull.nRelax = 3;
  outer_skull.coefTangential = 0.2;
  outer_skull.coefNormal = 0.2;
  outer_skull.coefMRI = 0.2;
  outer_skull.coefRepelling = 1.0;
  outer_skull.intensityThreshold = 50.0;

  outer_scalp.nMove = 200;
  outer_scalp.nRest = 30;
  outer_scalp.nRelax = 1;
  outer_scalp.coefTangential = 0.04;
  outer_scalp.coefNormal = 0.02;
  outer_scalp.coefMRI = 0.5;
  outer_scalp.coefRepelling = 0.5;
  outer_scalp.intensityThreshold = 35.0;
}

int
FileExists(const char *fname)
{
  FILE *fp ;

  fp = fopen(fname, "r") ;
  if (fp) fclose(fp) ;
  return(fp != NULL) ;
}

int isadir(const char *path)
{
  struct stat stbuf;

  stat(path,&stbuf);
  if ((stbuf.st_mode & S_IFMT) == S_IFDIR) return 1; else return 0;
}

void usage()
{
  printf("Usage: %s infile [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    infile               full path of input MRI volume (mgh format)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -ico                 [5]    icosahedral mesh order\n");
  printf("    -tri_file_dir        [ ]    directory with input ico tri files\n");
  printf("      default: %s\n",tri_file_dir.c_str());
  printf("    -outdir              [ ]    output directory\n");
  printf("      default: input directory\n");
  printf("    -outstem             [ ]    output file stem\n");
  printf("      default: (none)\n");
  printf("    -[no]inner_skull            create inner skull surface\n");
  printf("    -[no]outer_skull            create outer skull surface\n");
  printf("    -[no]outer_scalp            create outer scalp surface\n");
  printf("      default: create all three surfaces\n");
  printf("    -inner_skull_nMove   [300]  movement steps\n");
  printf("    -inner_skull_nRest   [30]   local smoothing steps\n");
  printf("    -inner_skull_nRelax  [2]    global smoothing steps per search\n");
  printf("    -inner_skull_cTang   [0.1]  tangential component coefficient\n");
  printf("    -inner_skull_cNorm   [0.1]  normal component coefficient\n");
  printf("    -inner_skull_cMRI    [0.3]  mri force coefficient\n");
  printf("    -inner_skull_cRepell [0.0]  repelling force coefficient\n");
  printf("    -inner_skull_thresh  [80.0] intensity threshold\n");
  printf("    -outer_skull_nMove   [100]  movement steps\n");
  printf("    -outer_skull_nRest   [30]   local smoothing steps\n");
  printf("    -outer_skull_nRelax  [3]    global smoothing steps per search\n");
  printf("    -outer_skull_cTang   [0.2]  tangential component coefficient\n");
  printf("    -outer_skull_cNorm   [0.2]  normal component coefficient\n");
  printf("    -outer_skull_cMRI    [0.2]  mri force coefficient\n");
  printf("    -outer_skull_cRepell [1.0]  repelling force coefficient\n");
  printf("    -outer_skull_thresh  [50.0] intensity threshold\n");
  printf("    -outer_scalp_nMove   [200]  movement steps\n");
  printf("    -outer_scalp_nRest   [30]   local smoothing steps\n");
  printf("    -outer_scalp_nRelax  [1]    global smoothing steps per search\n");
  printf("    -outer_scalp_cTang   [0.04] tangential component coefficient\n");
  printf("    -outer_scalp_cNorm   [0.02] normal component coefficient\n");
  printf("    -outer_scalp_cMRI    [0.5]  mri force coefficient\n");
  printf("    -outer_scalp_cRepell [0.5]  repelling force coefficient\n");
  printf("    -outer_scalp_thresh  [35.0] intensity threshold\n");
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i;
  char tempstr[STRLEN];
  char *pch;

  progname = argv[0];
  /* strip off path */
  pch = strtok(progname,"/");
  while (pch!=NULL) {
    strcpy(tempstr,pch);
    pch = strtok (NULL,"/");
  }
  strcpy(progname,tempstr);

  if (argc<MINARGC) {usage(); exit(0);}

  /* parse arguments */
  infile = argv[1];
  defaultCoefs();
  for (i=2;i<argc;i++) {
    if (MATCH(argv[i],"-ico") && i+1<argc) {
      ico = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outdir") && i+1<argc) {
      outdir = argv[i+1]; i++;
    } else 
    if (MATCH(argv[i],"-outstem") && i+1<argc) {
      outstem = argv[i+1]; i++;
    } else 
    if (MATCH(argv[i],"-tri_file_dir") && i+1<argc) {
      tri_file_dir = argv[i+1]; i++;
    } else 
    if (MATCH(argv[i],"-noouter_skull")) {
      outer_skull_flag = 0;
    } else 
    if (MATCH(argv[i],"-noouter_scalp")) {
      outer_scalp_flag = 0;
    } else 
    if (MATCH(argv[i],"-inner_skull_nMove") && i+1<argc) {
      inner_skull.nMove = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-inner_skull_nRest") && i+1<argc) {
      inner_skull.nRest = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-inner_skull_nRelax") && i+1<argc) {
      inner_skull.nRelax = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-inner_skull_cTang") && i+1<argc) {
      inner_skull.coefTangential = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-inner_skull_cNorm") && i+1<argc) {
      inner_skull.coefNormal = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-inner_skull_cMRI") && i+1<argc) {
      inner_skull.coefMRI = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-inner_skull_cRepell") && i+1<argc) {
      inner_skull.coefRepelling = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-inner_skull_thresh") && i+1<argc) {
      inner_skull.intensityThreshold = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_nMove") && i+1<argc) {
      outer_skull.nMove = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_nRest") && i+1<argc) {
      outer_skull.nRest = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_nRelax") && i+1<argc) {
      outer_skull.nRelax = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_cTang") && i+1<argc) {
      outer_skull.coefTangential = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_cNorm") && i+1<argc) {
      outer_skull.coefNormal = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_cMRI") && i+1<argc) {
      outer_skull.coefMRI = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_cRepell") && i+1<argc) {
      outer_skull.coefRepelling = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_skull_thresh") && i+1<argc) {
      outer_skull.intensityThreshold = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_nMove") && i+1<argc) {
      outer_scalp.nMove = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_nRest") && i+1<argc) {
      outer_scalp.nRest = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_nRelax") && i+1<argc) {
      outer_scalp.nRelax = atoi(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_cTang") && i+1<argc) {
      outer_scalp.coefTangential = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_cNorm") && i+1<argc) {
      outer_scalp.coefNormal = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_cMRI") && i+1<argc) {
      outer_scalp.coefMRI = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_cRepell") && i+1<argc) {
      outer_scalp.coefRepelling = atof(argv[i+1]); i++;
    } else 
    if (MATCH(argv[i],"-outer_scalp_thresh") && i+1<argc) {
      outer_scalp.intensityThreshold = atof(argv[i+1]); i++;
    } else 
    {
      printf("%s: ### unrecognized option: %s ...quitting\n",progname,argv[i]);
      usage();
      exit(0);
    }
  }

  /* check arguments */
  if (!FileExists(infile.c_str())) {
    printf("%s: ### input file %s not found ...quitting\n",
      progname,infile.c_str());
    exit(0);
  }
  if (!FileExists(tri_file_dir.c_str())) {
    printf("%s: ### tri_file_dir %s not found ...quitting\n",
      progname,tri_file_dir.c_str());
    exit(0);
  }
  if(!isadir(tri_file_dir.c_str())) {
    printf("%s: ### %s exists, but not a dir ...quitting\n",
      progname,tri_file_dir.c_str());
    exit(0);
  }
  if (!MATCH(outdir.c_str(),UNDEFSTR)) {
    if (!FileExists(outdir.c_str())) {
      printf("%s: ### outdir %s not found ...quitting\n",
        progname,outdir.c_str());
      exit(0);
    }
    if(!isadir(outdir.c_str())) {
      printf("%s: ### %s exists, but not a dir ...quitting\n",
        progname,outdir.c_str());
      exit(0);
    }
  }
  if (ico < 0 || ico > 7) {
    printf("%s: ### ico must be and integer between 0 and 7 ...quitting\n",
      progname);
    exit(0);
  };

  if (outer_scalp_flag)
    outer_skull_flag = 1;
}

int main(int argc, char** argv)
{

  parse_args(argc,argv);

  enum {NO=0, YES, BOTH}; 

  int surfaceDecimation, rasShift; 

  /* 4 directories for data input to be specified 

  1) mriInputString (mgh-file)

  2) triInputDirectory (ico-files) 

  3) projectDirectory (mean, template, eof and talairach) 

  4) projectDirectory/methodSubdirectory (surfaces) 

  */

  string mriInputString=infile;

  string triInputDirectory=tri_file_dir; 

  string projectDirectory="/space/monkeys/1/home/PROJECTS/VETSA_UCSD15T/"; 

  string methodSubdirectory="DETERMINISTIC/"; 

  int mriInputLength=mriInputString.length();  
  int mriInputPivot=mriInputString.rfind("/", mriInputLength); 

  string mriInputDirectory = mriInputString.substr(0, mriInputPivot+1);  
  string mriFileName = mriInputString.substr(mriInputPivot+1, mriInputLength);
  string outputDirectory;
  if MATCH(outdir.c_str(),UNDEFSTR)
    outputDirectory=mriInputDirectory;
  else
    outputDirectory=outdir + "/";

  int mriFileTypeLength = mriFileName.length()-4; 
  string subjectName = 
    mriInputString.substr(mriInputPivot+1,mriFileTypeLength); 

  string surfaceInputDirectory= 
    projectDirectory + methodSubdirectory + "TRI/";

  string meanSurfaceInputDirectory=projectDirectory + "MEAN/";

  string templateInputDirectory=projectDirectory + "MEAN/";

  string eofInputDirectory=projectDirectory + "MEAN/";

  string talairachInputDirectory=projectDirectory + "XFM/";

  string talairachFileName = 
    talairachInputDirectory + subjectName + ".xfm"; 

  string surfaceType, surfaceName; 

  /* size of spherical template */ 
  float templateSize; 

  /*  number of surface movements for surface segmentation */
  int kMove, nMove; 

  /* number of local smoothing operations at the end of segmentation */ 
  int kRest, nRest; 

  /* number of global smoothing operations after each move */    
  int kRelax, nRelax;

  /* number of smoothing operation before and after 
     each surface dilation */
  int nSmoothness; 

  /* the overall magnitude of surface dilation */  
  float movingDistance; 

  float coefTangential, coefNormal, coefMRI, coefRepelling; 
  float directionMRI; 
  float intensityThreshold; 

  deterministic *pMesh; 

  /* Mesh Number    File       Vertices     Faces  
   *  
   *     0         ic0.tri        12          20   
   *     1         ic1.tri        42          80  
   *     2         ic2.tri       162         320  
   *     3         ic3.tri       642        1280  
   *     4         ic4.tri      2562        5120  
   *     5         ic5.tri     10242       20480  
   *     6         ic6.tri     40962       81920  
   *     7         ic7.tri    163842      327680  
   *                                                    */ 
 
  /* initialize algorithm with the Mesh Number > 0 */  

  /* currently RAS coordinates are adjusted in getIntensity subroutine, 
   * set rasShift to NO everywhere */ 

  int talairachFlag=YES; 

  pMesh = new deterministic(mriFileName,
			    mriInputDirectory,  
			    talairachFlag, 
			    talairachFileName,
			    triInputDirectory, 
			    ico);

  /* beginning of the inner skull processing */ 

  /* set parameters for the inner skull segmentation */ 

 /* Optimal coefficients for the inner skull 
  * using "flash" files, i.e. without 1.e08 scaling
  * see old versions for the scaling case  
  *  
  *
  * tri  nStep  Templ. Tang. Normal   Smooth   MRI  Repel.  Threshold 
  *
  * ic5   300    30    0.2    0.15    0./0.1   0.3   0.        350. 
  *                                                                    */ 
  cout << "### inner skull ###" << endl;
  surfaceType = "inner_skull"; 
  if MATCH(outstem.c_str(),UNDEFSTR)
    surfaceName = surfaceType; 
  else
    surfaceName = outstem + "_" + surfaceType; 

  pMesh->initializeSurface(templateSize=30.);  

  pMesh->segmentSurface(surfaceName, 
			outputDirectory, 
			inner_skull.nMove, // 300
                        inner_skull.nRest, // 30 
                        inner_skull.nRelax, // 1
			inner_skull.coefTangential, // 0.2
			inner_skull.coefNormal, // 0.15
			inner_skull.coefMRI, // 0.3
			directionMRI=1, 
			inner_skull.coefRepelling, // 0.
			inner_skull.intensityThreshold); // 325.

  printf("%s: writing %s%s.tri\n",progname,outputDirectory.c_str(),surfaceName.c_str());
  pMesh->exportSurface(rasShift=NO, 
		       surfaceDecimation=NO, 
		       surfaceName, 
		       outputDirectory);

/******************************************************************/


  if(outer_skull_flag) {
  /* beginning of the outer skull processing */ 
  /* set parameters for the outer skull segmentation 
  *
  * tri  nStep  Templ. Tang. Normal   Smooth   MRI  Repel.  Threshold 
  *
  * ic5   100    --    0.15    0.1    0./0.1   0.3   0.2       250. 
  *                                                                    */ 
    cout << "### outer skull ###" << endl;
    surfaceType="outer_skull"; 
    if MATCH(outstem.c_str(),UNDEFSTR)
      surfaceName = surfaceType; 
    else
      surfaceName = outstem + "_" + surfaceType; 

    pMesh->createRepellingVolume(movingDistance=3., // 3   
			         nSmoothness=1);  // 1  


    pMesh->segmentSurface(surfaceName, 
			  outputDirectory, 
			  outer_skull.nMove, // 100
                          outer_skull.nRest, // 30 
                          outer_skull.nRelax, // 1
			  outer_skull.coefTangential, // 0.15
			  outer_skull.coefNormal, // 0.1
			  outer_skull.coefMRI, // 0.2
			  directionMRI=-1, 
			  outer_skull.coefRepelling, // 0.2 
			  outer_skull.intensityThreshold); // 250.

    printf("%s: writing %s%s.tri\n",progname,outputDirectory.c_str(),surfaceName.c_str());
    pMesh->exportSurface(rasShift=NO, 
		         surfaceDecimation=NO, 
		         surfaceName, 
		         outputDirectory);
  }

  /*******************************************************************/


  if(outer_scalp_flag) {
  /* beginning of the outer scalp processing */ 
  /* set parameters for the outer scalp segmentation */ 
 /* Optimal coefficients for the outer scalp segmentation  
  *  
  *
  * tri  nStep  Templ. Tang. Normal   Smooth   MRI  Repel.  Threshold 
  *
  * ic5   200    150   0.03   0.02    0./0.1   0.7    0.       125. 
  *                                                                    */ 

    cout << "### outer scalp ###" << endl;
    surfaceType="outer_scalp"; 
    if MATCH(outstem.c_str(),UNDEFSTR)
      surfaceName = surfaceType; 
    else
      surfaceName = outstem + "_" + surfaceType; 

    pMesh->createRepellingVolume(movingDistance=3., // 3   
			         nSmoothness=1);  // 1  

    pMesh->initializeSurface(templateSize=150.);  

    pMesh->segmentSurface(surfaceName, 
			  outputDirectory, 
			  outer_scalp.nMove, // 300
                          outer_scalp.nRest, // 30 
                          outer_scalp.nRelax, // 1
			  outer_scalp.coefTangential, // 0.04
			  outer_scalp.coefNormal, 
			  outer_scalp.coefMRI,
			  directionMRI=1., 
			  outer_scalp.coefRepelling, 
			  outer_scalp.intensityThreshold);// 125.

    printf("%s: writing %s%s.tri\n",progname,outputDirectory.c_str(),surfaceName.c_str());
    pMesh->exportSurface(rasShift=NO, 
		         surfaceDecimation=NO, 
		         surfaceName, 
		         outputDirectory);
  }

  delete pMesh;

  return 0;
}
