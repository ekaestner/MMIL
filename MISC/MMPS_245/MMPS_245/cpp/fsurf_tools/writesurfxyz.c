/* writesurfxyz.c: write surface xyz coordinates as three w files
      created: 03/13/05 DH
     last mod: 03/13/05 DH

   purpose:
     getting easy access to xyz coordinates

   input:
     none

   output:
     w value files
*/

#include "surflib.h"

#define MINARGC 3

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char outstem[STRLEN]=UNDEFSTR;
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char surf[STRLEN]="orig";
char hemi[STRLEN]="rh";

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);

int main(int argc, char **argv)
{
  int k, nverts;
  float *x, *y, *z;
  char fname[STRLEN];
  MRIS *mris;

  parse_args(argc,argv);

  mris = openSurface(subj,hemi,surf);
  nverts = mris->nvertices;

  x = (float *)calloc(nverts,sizeof(float));  MTEST(x);
  y = (float *)calloc(nverts,sizeof(float));  MTEST(y);
  z = (float *)calloc(nverts,sizeof(float));  MTEST(z);
  
  for(k=0;k<nverts;k++) {
    x[k] = mris->vertices[k].x;
    y[k] = mris->vertices[k].y;
    z[k] = mris->vertices[k].z;
  }

  sprintf(fname,"%s-%s-x-%s.w",outstem,surf,hemi);
  writeSurfVals(fname,x,nverts);
  sprintf(fname,"%s-%s-y-%s.w",outstem,surf,hemi);
  writeSurfVals(fname,y,nverts);
  sprintf(fname,"%s-%s-z-%s.w",outstem,surf,hemi);
  writeSurfVals(fname,z,nverts);

  free(x);
  free(y);
  free(z);

  exit(0);
}

void usage()
{
  printf("\n");
  printf("Usage: %s -name subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -name subjname             subject name\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outdir  <str>  [.]        output dir (must exist)\n");
  printf("    -outstem <str>  [subjname] outstem\n");
  printf("    -surf    <str>  [orig]     surface (orig, smoothwm, etc)\n");
  printf("    -hemi    <str>  [rh]       hemisphere (rh or lh)\n");
  printf("    -quiet                     suppress messages\n");
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
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-name") && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-oustem") && i+1<argc) {
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      } else 
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    } else
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### must supply a subject name ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    strcpy(outstem,subj);
  }
  if (!FileExists(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  }
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }

  MsgPrintf("%s: finished parsing arguments\n",progname);
}
