/* randsurfvals.c: generate random values for a surface
      created: 12/22/04 DH
     last mod: 12/23/04 DH

   purpose:
     generating random values at each vertex of a surface

   input:
     none

   output:
     random stats file (w file)

    acknowlegdgements:
      portions of code adapted from AFNI's AlphaSim
*/

#include "surflib.h"

#define MINARGC 3
#define BIGFLOAT 1e10
#define DEFOUTSTEM "randsurfvals-output"
#define RATIO(X,Y)        ((float)X/(float)Y)
#ifdef Darwin
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),INT_MAX)*((H)-(L)))
#else
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),MAXINT)*((H)-(L)))
#endif

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char outstem[STRLEN]=UNDEFSTR;
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char outdir[STRLEN]=".";
float stdev = 1.0;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -subj   subjname   subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem outstem    omit extension, hemi: <outstem>-rh.w\n");
  printf("    -outdir   [.]       output dir\n");
  printf("    -hemi     [rh]      hemisphere (rh or lh)\n");
  printf("    -stdev    [1.0]     standard deviation of random values\n");
  printf("    -quiet              suppress messages\n");
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
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-stdev") && i+1<argc) {
        stdev = atof(argv[i+1]); i++;
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
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    strcpy(outstem,DEFOUTSTEM);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  if (!FileExists(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  }
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

float uniform ()
{
  return ( (float)drand48() );
}

void normal (float * n1, float * n2)
{
  float u1, u2;
  float r;

  u1 = 0.0;
  while (u1 <= 0.0) {u1 = uniform();}
  u2 = uniform();

#if 0
  r = sqrt(-2.0*log(u1));
#else
  r = stdev * sqrt(-2.0*log(u1));
#endif
  *n1 = r * cos(2.0*M_PI*u2);
  *n2 = r * sin(2.0*M_PI*u2);
}

void generate_randsurf(MRIS *mris)
{
  int k, nverts, nvertsdiv2;
  float n1, n2;
  
  nverts = mris->nvertices;
  nvertsdiv2 = nverts/2;
  
  MsgPrintf("%s: generating random values\n",progname);
  for (k=0;k<nvertsdiv2;k++) {
    normal(&n1,&n2);
    mris->vertices[k].val = n1;
    mris->vertices[k+nvertsdiv2].val = n2;
  }
  normal(&n1,&n2);
  mris->vertices[nverts-1].val = n1;
}

int main(int argc, char **argv)
{
  int ecode;
  char valfile[STRLEN];
  MRIS *mris;

  randseed();
  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,"orig");

  generate_randsurf(mris);

  /* write to file */
  MsgPrintf("%s: writing output\n",progname);
  sprintf(valfile,"%s/%s-%s.w",outdir,outstem,hemi);
  ecode = MRISwriteValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error writing value file %s\n",
              progname,valfile);
  }
}

