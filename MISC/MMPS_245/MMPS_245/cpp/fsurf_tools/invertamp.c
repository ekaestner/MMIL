/* invertamp.c: read surface stats file and invert amplitudes
      created: 12/10/04 DH
     last mod: 12/27/04 DH

   purpose:
     making minima into maxima

   input:
     stats file (w file) -- real or complex

   output:
     stats file (w file) -- real or complex
*/

#include "surflib.h"

#define MINARGC 3
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
int complex = 0;
int subtractflag = 0;
char *real_infix=NULL;
char *imag_infix=NULL;

/* functions */

void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -name subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  instem           omit extension, infixes, hemi: <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -name    <str>            subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem  <str> [instem-grad] output file stem\n");
  printf("    -complex                      for complex data sets (real + imaginary)\n");
  printf("    -subtract                     invert by subtraction (rather than division)\n");
  printf("    -indir   <str>  [.]           input dir\n");
  printf("    -hemi    <str>  [rh]          hemisphere (rh or lh)\n");
  printf("    -infixes <r> <i>              real,imaginary infixes--default: %s %s\n",
         REAL_INFIX,IMAG_INFIX);
  printf("    -quiet                        suppress messages\n");
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

  real_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(real_infix);
  imag_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(imag_infix);
  strcpy(real_infix,REAL_INFIX);
  strcpy(imag_infix,IMAG_INFIX);

  /* parse arguments */
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-name") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-complex")){
        complex = 1;
      } else
      if (MATCH(argv[i],"-subtract")){
        subtractflag = 1;
      } else
      if ((MATCH(argv[i],"-infixes")) && i+2<argc){
        strcpy(real_infix,argv[i+1]); strcpy(imag_infix,argv[i+2]); i+=2;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
      {
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
      }
    } else {
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    }
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(instem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### instem not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    sprintf(outstem,"%s-grad",instem);
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
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}


void complex2polar(MRIS *mris)
{
  int k;
  float phase,amplitude;
  
  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].ripflag) {
      phase = atan2(mris->vertices[k].imag_val,mris->vertices[k].val);
      amplitude = hypot(mris->vertices[k].val,mris->vertices[k].imag_val);
      mris->vertices[k].imag_val = phase;
      mris->vertices[k].val = amplitude;
    }
  }
}

void polar2complex(MRIS *mris)
{
  int k;
  float phase,amplitude;
  
  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].ripflag) {
      phase = mris->vertices[k].imag_val;
      amplitude = mris->vertices[k].val;
      mris->vertices[k].val = amplitude*cos(phase);
      mris->vertices[k].imag_val = amplitude*sin(phase);
    }
  }
}


int main(int argc, char **argv)
{
  int k, ecode=NO_ERROR;
  float val,maxval=0;
  char tempstr[STRLEN];
  MRIS *mris;

  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,"orig");

  /* read input files */
  if(complex) {
    sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,real_infix,hemi);
    ecode = MRISreadValues(mris,tempstr);

    sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,imag_infix,hemi);
    ecode = MRISreadImagValues(mris,tempstr);
  } else {
    sprintf(tempstr,"%s/%s-%s.w",indir,instem,hemi);
    ecode = MRISreadValues(mris,tempstr);
  }

  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  }

  if(complex) complex2polar(mris);

  if(!subtractflag){
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].ripflag) {
        val = mris->vertices[k].val;
        if(val!=0) val = 1/val;
        mris->vertices[k].val = val;
      }
    }
  } else {
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].ripflag) {
        if(mris->vertices[k].val>maxval) maxval=mris->vertices[k].val;
      }
    }
    if(maxval==0)
      ErrorExit(ERROR_BADFILE,"%s: max value = 0 ... quitting\n",progname);
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].ripflag) {
        val = mris->vertices[k].val;
        val = 10.0*(1.0-val/maxval);
        mris->vertices[k].val = val;
      }
    }
  }  

  if(complex) polar2complex(mris);
  
  /* write output to file */
  if(complex) {
    sprintf(tempstr,"%s/%s%s-%s.w",outdir,outstem,real_infix,hemi);
    ecode = MRISwriteValues(mris,tempstr);

    sprintf(tempstr,"%s/%s%s-%s.w",outdir,outstem,imag_infix,hemi);
    ecode = MRISwriteImagValues(mris,tempstr);
  } else {
    sprintf(tempstr,"%s/%s-%s.w",outdir,outstem,hemi);
    ecode = MRISwriteValues(mris,tempstr);
  }

  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error writing output files %s\n",
              progname,tempstr);
  }

  exit(0);
}

