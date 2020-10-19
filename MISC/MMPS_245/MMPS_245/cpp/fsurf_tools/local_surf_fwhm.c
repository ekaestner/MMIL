/* local_surf_fwhm.c: read surface stats file and calculate 2d gaussian smoothness
      created: 08/04/04 DH
     last mod: 01/30/06 DH

   purpose:
     calculating 2d gaussian smoothness to input to alphasim
       alphasim is used to simulate clustering, to generate
       cluster-size thresholds for multiple comparison correction

   input:
     stats file (w file)

   output:
     stdout only

   acknowledgements:
     this code is adapted from AFNI's 3dFWHM
     (B.D. Ward, Medical College of Wisconsin 1997)
    and uses equations and resels concept from Worsley et al., 1999
      "Detecting Changes in Nonisotropic Images", Human Brain Mapping 8:98-101
    and Hayasaka et al., 2004
      "Nonstationary cluster-size inference with random field and permutation
       methods", NeuroImage 22:676-687
*/

#include "surflib.h"

#define MINARGC 3

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char surf[STRLEN]="smoothwm";
int reselflag = 0;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -name subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  instem              omit extension, hemi: <instem>-rh.w\n");
  printf("    -subj    <str>               subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -indir   <str> [.]           input dir\n");
  printf("    -outstem <str> [instem-fwhm] omit extension, hemi: <outstem>-rh.w\n");
  printf("    -outdir  <str> [.]           output dir\n");
  printf("    -surf          [smoothwm]    surface used for area calculations\n");
  printf("    -hemi    <str> [rh]          hemisphere (rh or lh)\n");
  printf("    -resel                       output resel's instead of fwhm\n");
  printf("    -quiet                       suppress messages\n");
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
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc) {
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-resel")){
        reselflag=1;
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
    sprintf(outstem,"%s-fwhm",instem);
  }
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
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
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

int main(int argc, char **argv)
{
  int i, j, k, dfcount=0, ecode, zerocount=0;
  char valfile[STRLEN];
  float sum=0, sumsq=0, var=0, df=0, dfsum=0, dfsumsq=0, dist=0, fwhm;
  MRIS *mris;

  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,surf);

  /* read values */
  sprintf(valfile,"%s/%s-%s.w",indir,instem,hemi);
  ecode = MRISreadValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,valfile);
  }

  /* calculate sums and sums of squares */
  for (i=0;i<mris->nvertices;i++) {
    sum   += mris->vertices[i].val;
    sumsq += mris->vertices[i].val * mris->vertices[i].val;

    /* estimage partial derivatives between neighboring vertices */
    dfsum=dfsumsq=0;
    dfcount=0;
    for (j=0;j<mris->vertices[i].vnum;j++) {
      k = mris->vertices[i].v[j];
      dist = mris->vertices[i].dist[j];
      if(dist<=0) continue;
      df = (mris->vertices[k].val - mris->vertices[i].val)/dist;
      dfsum   += df;
      dfsumsq += df * df;
      dfcount++;
    }
    /* estimate variance of partial derivatives */
    if(dfcount>1) {
      mris->vertices[i].val2 = (dfsumsq - (dfsum*dfsum)/dfcount)/(dfcount-1);
      if(mris->vertices[i].val2<0) {
        mris->vertices[i].val2 = 0;
      }
    } else {
      mris->vertices[i].val2 = 0;
    }
  }

  /* estimate variance of data across all vertices */
  var =  (sumsq - (sum*sum)/mris->nvertices)/(mris->nvertices-1);
  if(var==0)
    ErrorExit(ERROR_BADFILE,"%s: ### error: overal var is zero\n",
              progname);

  /* calculate local smoothness for each vertex */
  for (i=0;i<mris->nvertices;i++) {
    df = mris->vertices[i].val2/var;
//    df = 1.0 - 0.5*(mris->vertices[i].val2/var);
    if(df>0) {
//      fwhm = sqrt(-1.0 / (4.0 * log(df))) * 2.0*sqrt(2.0*log(2.0));
//      fwhm = sqrt(-2.0*log(2.0)/log(df));
      fwhm = sqrt(4.0*log(2.0)/df);
    } else {
      fwhm=0;
    }
    if(!finite(fwhm)) fwhm=0;
    if(reselflag) {
      if(!finite(mris->vertices[i].area))
        mris->vertices[i].val = 0;
      else if(fwhm!=0)
        mris->vertices[i].val = mris->vertices[i].area/(fwhm*fwhm);
      else
        mris->vertices[i].val = mris->vertices[i].area;
    } else
      mris->vertices[i].val = fwhm;
    if(fwhm==0) zerocount++;
  }
  
  if(zerocount)
    printf("%s: %d vertices with fwhm = 0\n",progname,zerocount);

  /* write output file */
  sprintf(valfile,"%s/%s-%s.w",outdir,outstem,hemi);
  ecode = MRISwriteValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error writing value file %s\n",
              progname,valfile);
  }

  exit(0);
}

