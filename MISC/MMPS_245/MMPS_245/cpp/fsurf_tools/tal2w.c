/* tal2w.c: convert list of talairach points to w file
      created: 09/17/03 DH
     last mod: 09/23/03 DH

   purpose:
     painting lists of talairach points

   input:
     text file with list of xyz coordinates

   output:
     surface value file (.w extension)
*/
#include "surflib.h"

#define MINARGC 3

/* global variables */
char progname[20]="tal2w";
MRI_SURFACE *mris = NULL ;

/* parameter defaults */
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char pointsfile[STRLEN]=UNDEFSTR;
char outdir[STRLEN]=".";
char outstem[STRLEN]="output";
float paintvalue = 100.0;
int smoothsteps = 0;

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);
int readPoint(FILE *fp, double *x, double *y, double *z);

int main(int argc, char **argv)
{
  Transform mult_transform;
  int vindex;
  char fname[STRLEN];
  FILE *fp;
  double x,y,z;

  if (argc<MINARGC) {usage(); exit(0);}
  parse_args(argc,argv);

  mris = openSurface(subj,hemi,"orig");

  if(mris->transform_loaded) {
    printTransform(mris->linear_transform,"linear_transform");
  } else {
    ErrorExit(ERROR_BADFILE, "%s: transform NOT loaded succesfully\n", progname);
  }
  
  if(mris->inverse_transform_loaded) {
    printTransform(mris->inverse_linear_transform,"inverse_linear_transform");
  } else {
    ErrorExit(ERROR_BADFILE, "%s: inverse transform NOT loaded succesfully\n",
              progname);
  }

  /* multiply matrices to check that inverse is correct */
  multTransform(&mult_transform,mris->linear_transform,
                                mris->inverse_linear_transform);
  printTransform(&mult_transform,"mult_transform");

  fp = fopen(pointsfile,"r");
  if(fp==NULL) {
    ErrorExit(ERROR_NOFILE,"%s: ### Points file %s not found\n",
              progname,pointsfile);
  }

  while(readPoint(fp,&x,&y,&z)!=-1) {
    selectTalairachPoint(mris,&vindex,x,y,z);
    mris->vertices[vindex].val = paintvalue;
  }
  fclose(fp);

  MRISsmoothValues(mris,smoothsteps);
  sprintf(fname,"%s/%s-%s.w",outdir,outstem,hemi);
  MsgPrintf("%s: painting vertices to file: %s\n",progname,fname);
  MRISwriteValues(mris,fname);

  exit(0);
}

void usage()
{
  printf("\n");
  printf("Usage: %s pointsfile subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    pointsfile   full path name of file\n");
  printf("                   containing list of x y z talairach coordinates\n");
  printf("    subjname     name of subject in freesurfer subject dir\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -hemi    <str>  [rh]       hemisphere (rh or lh)\n");
  printf("    -outdir  <str>  [.]        output dir (must exist)\n");
  printf("    -outstem <str>  [output]   output file stem\n");
  printf("    -value  <float> [100.0]    value to assign to selected vertices\n");
  printf("    -smooth  <int>  [0]        smoothing steps (on surface)\n");
  printf("    -quiet                     suppress messages\n");
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i,a=1;
  FILE *fp;
  
  /* parse arguments */
  strcpy(pointsfile,argv[a++]);
  strcpy(subj,argv[a++]);
  for (i=a;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outstem") && i+1<argc) {
        strcpy(outstem,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-value") && i+1<argc) {
        paintvalue = atof(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-smooth") && i+1<argc) {
        smoothsteps = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      }
      else {
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
      }
    } else {
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    }
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (MATCH(pointsfile,"pointsfile")) {
    ErrorExit(ERROR_BADPARM,"%s: ### points file name not specified ...quitting\n",
              progname);
  }
  fp = fopen(outdir,"r");
  if(fp==NULL) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  } else fclose(fp);
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  if(smoothsteps<0) {
    ErrorExit(ERROR_BADPARM,"%s: smooth steps (%d) must be > 0\n\n",
              progname,smoothsteps);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

int readPoint(FILE *fp, double *x, double *y, double *z)
{
  char funcname[STRLEN]="readPoint";
  char line[MAX_LINE_LENGTH];

  if(fp==NULL) {
    printf("%s: ### Point file not open!\n",funcname);
    return(-1);
  }
  if(fgets(line,MAX_LINE_LENGTH-1,fp)==NULL) return(-1);
  sscanf(line,"%lf %lf %lf",x,y,z);

  return(0);
}
