/* smoothsurfvals.c: smooth functional data across cortical surface
      created: 09/08/03 DH
     last mod: 06/24/06 DH

   purpose:
     2D surface smoothing of functional data

   input:
     w value files

   output:
     w value files
*/

#include "surflib.h"

#define MINARGC 3
#define SMALLFLOAT 1e-4

/* global variables */
char progname[20]="smoothsurfvals";
MRIS *mris;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char surf[STRLEN]="smoothwm";
int niter = 10;
int maxiter = 1000;
int sparseflag = 0;
int heatflag = 0;
float sigma = 0.2;
int nsoap = 1; /* number of soap bubble iterations (1=normal smoothing) */
int sparse_pre_niter = 0; /* number of sparse pre-smoothing iterations */

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);

int main(int argc, char **argv)
{
  int ecode;
  char fname[STRLEN];
  int i,k;
  float *orig_vals;


  if (argc<MINARGC) {usage(); exit(0);}
  parse_args(argc,argv);

  mris = openSurface(subj,hemi,surf);

  /* read data */
  sprintf(fname,"%s/%s-%s.w",indir,instem,hemi);
  MsgPrintf("%s: value fname = %s\n",progname,fname);
  ecode = MRISreadValues(mris,fname);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,fname);
  }

  /* copy original values */
  orig_vals = (float *)calloc(mris->nvertices,sizeof(float));    MTEST(orig_vals);
  for(k=0;k<mris->nvertices;k++) {
    orig_vals[k] = mris->vertices[k].val;
  }

  for(i=0;i<nsoap;i++) {
    if(i>0) {
      /* refresh values */
      for(k=0;k<mris->nvertices;k++) {
        mris->vertices[k].val = orig_vals[k];
      }
    }

    /* smooth values */
    if (heatflag) {
      MsgPrintf("%s: heat kernel smoothing values %d times\n",progname,niter);
      MRISsmoothValuesHeat(mris,niter,sigma);
    } else {
      if(sparse_pre_niter) {
        MsgPrintf("%s: sparse pre-smoothing values %d times\n",
          progname,sparse_pre_niter);
        MRISsmoothValuesSparse(mris,sparse_pre_niter);
      }
      MsgPrintf("%s: smoothing values %d times\n",progname,niter);
      if(sparseflag) {
        MRISsmoothValuesSparse(mris,niter);
      } else {
        MRISsmoothValues(mris,niter);
      }
    }
  }

  /* write values to file */
  sprintf(fname,"%s/%s-%s.w",outdir,outstem,hemi);
  MsgPrintf("%s: value file name = %s\n",progname,fname);
  ecode = MRISwriteValues(mris,fname);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode, "%s: ### error writing value file %s\n",
              progname,fname);
  }
  exit(0);
}

void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  instem              omit infix, extension, hemi: <instem>-rh.w\n");
  printf("    -subj   <str>                subject name\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem outstem  [instem-sm] omit extension, hemi: <outstem>-rh.w\n");
  printf("    -indir   <str>    [.]         input dir\n");
  printf("    -outdir  <str>    [.]         output dir\n");
  printf("    -hemi    <str>    [rh]        hemisphere (rh or lh)\n");
  printf("    -sparse_pre_niter [0]         number of sparse smoothing iterations\n");
  printf("                                   applied before regular smoothing\n");
  printf("    -niter    <d>     [10]        number of smoothing iterations\n");
  printf("    -sparse                       sparse smoothing (fill in undefined vals)\n");
  printf("    -heat                         use heat kernel smoothing\n");
  printf("    -sigma    <f>     [0.2]       sigma for heat kernel smoothing\n");
  printf("    -nsoap    <i>     [1]         number of soap-bubble smoothing steps\n");
  printf("                                  nsoap=1 is normal smoothing\n");
  printf("    -quiet                        suppress messages\n");
  printf("\n");
}



void parse_args(int argc, char **argv)
{
  int i;
  
  /* parse arguments */
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-sparse_pre_niter") && i+1<argc){
        sparse_pre_niter = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-niter") && i+1<argc){
        niter = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-sparse")){
        sparseflag = 1;
      }
      else if (MATCH(argv[i],"-heat")){
        heatflag = 1;
      }
      else if (MATCH(argv[i],"-sigma") && i+1<argc){
        sigma = atof(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-nsoap") && i+1<argc) {
        nsoap = atoi(argv[i+1]); i++;
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
  if (MATCH(instem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### must supply an instem ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    sprintf(outstem,"%s-sm%d",instem,niter);
  }
  if (!FileExists(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",
              progname,indir);
  }
  if(!isadir(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
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
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",
              progname,hemi);
  }
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (sparse_pre_niter < 0) {
    ErrorExit(ERROR_BADPARM,"%s: ### number of sparse pre-smoothing iterations must be greater than zero ...quitting\n",
              progname);
  }
  if (niter < 0) {
    ErrorExit(ERROR_BADPARM,"%s: ### number of iterations must be greater than zero ...quitting\n",
              progname);
  }
  if (sparse_pre_niter > maxiter) {
    ErrorExit(ERROR_BADPARM,"%s: ### maximum number of iterations allowed is %d...quitting\n",
              progname,maxiter);
  }
  if (niter > maxiter) {
    ErrorExit(ERROR_BADPARM,"%s: ### maximum number of iterations allowed is %d...quitting\n",
              progname,maxiter);
  }
  if (sigma <= 0) {
    ErrorExit(ERROR_BADPARM,"%s: ### sigma must be greater than zero ...quitting\n",
              progname);
  }
  if (nsoap < 1) nsoap=1;

  MsgPrintf("%s: finished parsing arguments\n",progname);
}
