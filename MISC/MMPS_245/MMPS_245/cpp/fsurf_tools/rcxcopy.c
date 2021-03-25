/* rcxcopy.c: copies a real or complex bfloat or w data set
  for complex data: optionally reverses phase or truncates values outside
                    a selected range of phases
                    or convert to polar coordinates (amplitude and phase)
  for real data: optionally truncates values above or below a selected value
      created: 09/20/03 DH
     last mod: 02/24/06 DH

   purpose:
     reversing phase of polar angle mapping data before averaging with
       opposite direction data
     removing ipsilateral noise from polar angle mapping data
     thresholding stats

   input:
     bfloat volumes (not timeseries), 
     w value files

   output:
     bfloat volumes
     w value files
*/

#include "surflib.h"

#define MINARGC 2
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"
#define AMPL_INFIX "_a"
#define PHAS_INFIX "_p"
#define SAMPL_INFIX "_sa"
#define BIGFLOAT 1e10
#define DEFOUTSTEM "rcxcopy-output"

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char outstem[STRLEN]=UNDEFSTR;
char maskstem[STRLEN]=UNDEFSTR;
char maskdir[STRLEN]=".";
char ext[STRLEN]="bfloat";
char out_ext[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char subj[STRLEN]=UNDEFSTR;
char *real_infix=NULL;
char *imag_infix=NULL;
char *ampl_infix=NULL;
char *phas_infix=NULL;
char *sampl_infix=NULL;
int nslices=0;
int xnum=0,ynum=0,depth=0;
int outtext = 0;
int halfmaxflag = 0;
int complex = 0;
int cmplx2polar = 0;
int polar2cmplx = 0;
int cmplx2samp = 0; /* convert complex to signed amplitude */
int normcxamp = 0; /* set complex amplitude to 1 */
int truncflag = 0;
int revflag = 0;
int maskflag = 0;
int threshflag = 0;
int threshabsflag = 0;
float threshmin = -BIGFLOAT;
float threshmax = BIGFLOAT;
float offset = 0.0;
float phase0 = 0.0;
float phase1 = 0.5;
float slope = 10.0;

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);
void read_stats(float ***stats, char *indir, char *instem, char *infix);
void write_stats(float ***stats, char *outdir, char *outstem, \
                                               char *prefix, char *infix);
void write_stats_text(float ***stats, char *outdir, char *outstem, \
                                                    char *prefix, char *infix);
void subphase_stats(float ***stats_r, float ***stats_i);
void revphase_stats(float ***stats_r, float ***stats_i);
void truncphase_stats(float ***stats_r, float ***stats_i);
void cmplx2polar_stats(float ***stats_r, float ***stats_i);
void polar2cmplx_stats(float ***stats_r, float ***stats_i);
void cmplx2samp_stats(float ***stats_r, float ***stats_i);
void normcxamp_stats(float ***stats_r, float ***stats_i);
void thresh_stats(float ***stats, float min, float max, int absflag);
void thresh_cx_stats(float ***stats_r, float ***stats_i, float min, float max);
void mask_stats(float ***stats_r, float ***stats_i);
float getMax(float ***stats);
float getcxMax(float ***stats_r, float ***stats_i);

int main(int argc, char **argv)
{
  int i,n;
  float ***dat_r=NULL;
  float ***dat_i=NULL;
  char fname[STRLEN];

  parse_args(argc,argv);

  if (MATCH(ext,"bfloat")) {
    if (MATCH(out_ext, "w")) {
      if(complex) {
        sprintf(fname,"%s/%s%s-%s_000.hdr",indir,instem,real_infix,hemi);
      } else {
        sprintf(fname,"%s/%s-%s_000.hdr",indir,instem,hemi);
      }
      if(readFloatHeader(fname,&xnum,&ynum,&depth)==-1) {
        ErrorExit(ERROR_BADFILE,"%s: ### error reading header file...quitting\n",
               progname);
      }
      if (xnum!=1) {
        ErrorExit(ERROR_BADFILE,"%s: ### bfloat to w conversion only for 1D bfloats (xsize=1)...quitting\n",
                  progname);
      }
    } else {
      if(complex) {
        sprintf(fname,"%s/%s%s_000.hdr",indir,instem,real_infix);
      } else {
        sprintf(fname,"%s/%s_000.hdr",indir,instem);
      }
      if(readFloatHeader(fname,&xnum,&ynum,&depth)==-1) {
        ErrorExit(ERROR_BADFILE,"%s: ### error reading header file...quitting\n",
               progname);
      }
    }
  }
  
  if (MATCH(ext,"w") || MATCH(out_ext,"w")) {
    ynum = nverticesSurf(subj,hemi);
    if(ynum==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error reading number of vertices...quitting\n",
             progname);
    }
    xnum = 1;
  }

  /* allocate memory for data and calculations */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  dat_r = (float ***)calloc(nslices,sizeof(float **));       MTEST(dat_r);
  for (n=0;n<nslices;n++) {
    dat_r[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_r);
    for (i=0;i<ynum;i++) {
      dat_r[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_r);
    }
  }
  if(complex) {
    dat_i = (float ***)calloc(nslices,sizeof(float **));       MTEST(dat_i);
    for (n=0;n<nslices;n++) {
      dat_i[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_i);
      for (i=0;i<ynum;i++) {
        dat_i[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_i);
      }
    }
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* read data images */
  if(polar2cmplx) {
    read_stats(dat_r,indir,instem,ampl_infix);
    read_stats(dat_i,indir,instem,phas_infix);
  } else {
    read_stats(dat_r,indir,instem,real_infix);
    if(complex)
      read_stats(dat_i,indir,instem,imag_infix);
  }

  /* optional manipulations */
  if(complex) {
    if(threshflag) {
      if(halfmaxflag)
        threshmin = getcxMax(dat_r,dat_i)/2.0;
      thresh_cx_stats(dat_r,dat_i,threshmin,threshmax);
    }
    if(offset!=0.0) subphase_stats(dat_r,dat_i);
    if(revflag)     revphase_stats(dat_r,dat_i);
    if(cmplx2samp)  cmplx2samp_stats(dat_r,dat_i);
    else {
      if(polar2cmplx) polar2cmplx_stats(dat_r,dat_i);
      if(truncflag)   truncphase_stats(dat_r,dat_i);
      if(cmplx2polar) cmplx2polar_stats(dat_r,dat_i);
      if(normcxamp)   normcxamp_stats(dat_r,dat_i);
    }
  } else {
    if(threshflag) {
      if(halfmaxflag)
        threshmin = getMax(dat_r)/2.0;
      thresh_stats(dat_r,threshmin,threshmax,threshabsflag);
    }
  }
  if(maskflag) mask_stats(dat_r,dat_i);

  /* write output to file */
  if(cmplx2polar) {
    strcpy(real_infix,ampl_infix);
    strcpy(imag_infix,phas_infix);
  }
  if(cmplx2samp) {
    strcpy(real_infix,sampl_infix);
  }
  write_stats(dat_r,outdir,outstem,NULL,real_infix);
  if(complex && !cmplx2samp)
    write_stats(dat_i,outdir,outstem,NULL,imag_infix);
  if (outtext) {
    write_stats_text(dat_r,outdir,outstem,NULL,real_infix);
    if(complex && !cmplx2samp)
      write_stats_text(dat_i,outdir,outstem,NULL,imag_infix);
  }

  /* free memory */
  MsgPrintf("%s: starting to free memory\n",progname);
  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++) {
      if(dat_r[n][i]!=NULL) free(dat_r[n][i]);
    }
    if(dat_r[n]!=NULL) free(dat_r[n]);
  }
  if(dat_r!=NULL) free(dat_r);

  if(complex) {
    for (n=0;n<nslices;n++) {
      for (i=0;i<ynum;i++) {
        if(dat_i[n][i]!=NULL) free(dat_i[n][i]);
      }
      if(dat_i[n]!=NULL) free(dat_i[n]);
    }
    if(dat_i!=NULL) free(dat_i);
  }
  MsgPrintf("%s: finished freeing memory\n",progname);

  exit(0);
}

void usage()
{
  printf("\n");
  printf("Usage: %s instem {-subj subjname} [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    instem   omit infixes, suffix: <instem>_{r,i}_000.bfloat\n");
  printf("                                or <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -subj   <str>  subject name (required only if ext = w)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -indir    <str>  [.]       input dir\n");
  printf("    -outdir   <str>  [.]       output dir (must exist)\n");
  printf("    -outstem  <str>  [output]  output file stem\n");
  printf("    -maskdir  <str>  [.]       dir containing maskfile\n");
  printf("    -maskstem <str>  [NONE]    mask file stem\n");
  printf("    -nslices  <int>  [0]       slice count (bfloat ext only)\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -ext      <str>  [bfloat]  input extension (bfloat or w)\n");
  printf("    -out_ext  <str>  [ext]     output extension (bfloat or w)\n");
  printf("    -hemi     <str>  [rh]      hemisphere (rh or lh) (w or 1D-bfloat only)\n");
  printf("    -outtext                   save stats to additional text files\n");
  printf("    -quiet                     suppress messages\n");
  printf("\n");
  printf("  Optional parameters for real data only:\n");
  printf("    -threshmin <f>  [-10^10]   set to zero if values are less\n");
  printf("    -threshmax <f>  [+10^10]   set to zero if values are greater\n");
  printf("    -threshmin_halfmax         use half of global max as minimum threshold\n");
  printf("    -threshabs                 when thresholding, use absolute values\n");
  printf("\n");
  printf("  Optional parameters for complex data only:\n");
  printf("    -complex                   for complex data sets (real + imaginary)\n");
  printf("    -cmplx2polar               convert complex to polar (amplitude + phase)\n");
  printf("    -polar2cmplx               convert polar to complex (real + imaginary)\n");
  printf("    -cmplx2samp                convert complex to signed amplitude\n");
  printf("    -normcxamp                 set amplitude to 1 (but keep as complex)\n");
  printf("    -offset    <f>  [0.0]      phase offset value (fraction of 2*Pi)\n");
  printf("       subtraction of offset is done prior to reversal and truncation\n");
  printf("    -revphase                  takes negative of phase values\n");
  printf("       phase reversal is done prior to truncation\n");
  printf("    -truncphase                sets values to zero if phase lies outside range\n");
  printf("    -phase0    <f>  [0.0]      start of preserved phases\n");
  printf("    -phase1    <f>  [0.5]      end of preserved phases\n");
  printf("    -slope     <f>  [10.0]     slope of fall-off\n");
  printf("    -cx_infixes <r> <i>        real,imaginary infixes--default: %s %s\n",
         REAL_INFIX,IMAG_INFIX);
  printf("    -polar_infixes <r> <i>     amplitude,phase infixes--default: %s %s\n",
         AMPL_INFIX,PHAS_INFIX);
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i;
  char *tempfix=NULL;
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

  tempfix     = (char *)malloc(STRLEN*sizeof(char));      MTEST(tempfix);
  real_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(real_infix);
  imag_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(imag_infix);
  ampl_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(ampl_infix);
  phas_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(phas_infix);
  sampl_infix = (char *)malloc(STRLEN*sizeof(char));      MTEST(sampl_infix);
  strcpy(real_infix,REAL_INFIX);
  strcpy(imag_infix,IMAG_INFIX);
  strcpy(ampl_infix,AMPL_INFIX);
  strcpy(phas_infix,PHAS_INFIX);
  strcpy(sampl_infix,SAMPL_INFIX);

  /* parse arguments */
  strcpy(instem,argv[1]);  
  for (i=2;i<argc;i++) {
    if (argv[i][0]=='-') {
      if ((MATCH(argv[i],"-nslices") || MATCH(argv[i],"-n")) && i+1<argc) {
        nslices = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-maskdir") && i+1<argc) {
        strcpy(maskdir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-maskstem") && i+1<argc){
        strcpy(maskstem,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-ext") && i+1<argc){
        strcpy(ext,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-out_ext") && i+1<argc){
        strcpy(out_ext,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outtext")){
        outtext = 1;
      }
      else if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      }
      else if (MATCH(argv[i],"-threshmin") && i+1<argc){
        threshmin = atof(argv[i+1]); i++;
        threshflag = 1;
      }
      else if (MATCH(argv[i],"-threshmax") && i+1<argc){
        threshmax = atof(argv[i+1]); i++;
        threshflag = 1;
      }
      else if (MATCH(argv[i],"-threshmin_halfmax")){
        halfmaxflag = 1;
        threshflag = 1;
      }
      else if (MATCH(argv[i],"-threshabs")){
        threshabsflag = 1;
      }
      else if (MATCH(argv[i],"-complex")){
        complex = 1;
      }
      else if (MATCH(argv[i],"-cmplx2polar")){
        complex = 1;
        cmplx2polar = 1;
      }
      else if (MATCH(argv[i],"-polar2cmplx")){
        complex = 1;
        polar2cmplx = 1;
      }
      else if (MATCH(argv[i],"-cmplx2samp")){
        complex = 1;
        cmplx2samp = 1;
      }
      else if (MATCH(argv[i],"-normcxamp")){
        complex = 1;
        normcxamp = 1;
      }
      else if (MATCH(argv[i],"-offset") && i+1<argc){
        offset = atof(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-revphase")){
        revflag = 1;
      }
      else if (MATCH(argv[i],"-truncphase")){
        truncflag = 1;
      }
      else if ((MATCH(argv[i],"-phase0")) && i+1<argc) {
        phase0 = atof(argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-phase1")) && i+1<argc) {
        phase1 = atof(argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-slope")) && i+1<argc){
        slope = atof(argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-cx_infixes")) && i+2<argc){
        strcpy(real_infix,argv[i+1]); strcpy(imag_infix,argv[i+2]); i+=2;
      }
      else if ((MATCH(argv[i],"-polar_infixes")) && i+2<argc){
        strcpy(ampl_infix,argv[i+1]); strcpy(phas_infix,argv[i+2]); i+=2;
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
    strcpy(outstem,DEFOUTSTEM);
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
  if (!FileExists(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### maskdir %s not found ...quitting\n",
              progname,maskdir);
  }
  if(!isadir(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,maskdir);
  }
  if (!MATCH(maskstem,UNDEFSTR)) {
    maskflag=1;
  }
  if (!complex) {
    free(real_infix);
    real_infix = NULL;
    free(imag_infix);
    imag_infix = NULL;
  }
  if (!polar2cmplx && !cmplx2polar) {
    free(ampl_infix);
    ampl_infix = NULL;
    free(phas_infix);
    phas_infix = NULL;
  }
  if (!cmplx2samp) {
    free(sampl_infix);
    sampl_infix = NULL;
  }
  if (polar2cmplx && cmplx2polar) {
    ErrorExit(ERROR_BADPARM,"%s: ### cannot do both polar2cmplx and cmplx2polar ...quitting\n",
              progname,ext);
  }

  if (MATCH(out_ext,UNDEFSTR)) strcpy(out_ext,ext);
  if (!MATCH(ext,"bfloat") && !MATCH(ext,"w")) {
    ErrorExit(ERROR_BADPARM,"%s: ### extension %s not supported ...quitting\n",
              progname,ext);
  }
  if (!MATCH(out_ext,"bfloat") && !MATCH(out_ext,"w")) {
    ErrorExit(ERROR_BADPARM,"%s: ### extension %s not supported ...quitting\n",
              progname,out_ext);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",
              progname,hemi);
  }
  if (MATCH(ext,"w")) {
    nslices = 1;
  } else if (MATCH(ext, "bfloat")) {
    /* if converting from 1D-bfloat (surface) to w, then use hemi in file name */
    if (MATCH(out_ext, "w")) {
      if(real_infix==NULL) {
        sprintf(tempfix,"-%s",hemi);
      } else {
        sprintf(tempfix,"%s-%s",real_infix,hemi);
      }
    } else {
      if(real_infix==NULL) {
        tempfix = NULL;
      } else {
        strcpy(tempfix,real_infix);
      }
    }
    i = getNumBfloats(indir,instem,tempfix);
    if(i<=0) {
      if(tempfix==NULL) {
        ErrorExit(ERROR_BADPARM,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
               progname,indir,instem);
      } else {
        ErrorExit(ERROR_BADPARM,"%s: ### no bfloats found in %s/ with stem %s%s ...quitting\n",
               progname,indir,instem,tempfix);
      }

    }
    if(nslices==0 || nslices>i) {
      nslices = i;
      if(tempfix==NULL) {
        MsgPrintf("%s: %d slice(s) found in %s/ with stem %s\n",
                  progname,i,indir,instem);
      } else {
        MsgPrintf("%s: %d slice(s) found in %s/ with stem %s%s\n",
                  progname,i,indir,instem,tempfix);
      }
    }
    if (MATCH(out_ext, "w") && nslices > 1) {
      ErrorExit(ERROR_BADPARM,"%s: ### bfloat to w conversion only for 1D blfoats (xsize=1) ...quitting\n",
                progname);
    }
  }
  if (MATCH(subj,UNDEFSTR) && (MATCH(ext,"w") || MATCH(out_ext, "w"))) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (threshmin > threshmax) {
    ErrorExit(ERROR_BADPARM,"%s: ### threshmin (%f) > threshmax (%f) ...quitting\n",
              progname,phase0);
  }
  if (threshmin > BIGFLOAT && halfmaxflag) {
    MsgPrintf("%s: threshmin ignored since threshmin_halfmax is selected\n",progname);
  }
  if (phase0 < 0 || phase0 > 1) {
    ErrorExit(ERROR_BADPARM,"%s: ### bad phase0: %f => must be between 0.0 and 1.0 ...quitting\n",
              progname,phase0);
  }
  if (phase1 < 0 || phase1 > 1) {
    ErrorExit(ERROR_BADPARM,"%s: ### bad phase1: %f => must be between 0.0 and 1.0 ...quitting\n",
              progname,phase1);
  }
  if (phase0 == phase1) {
    ErrorExit(ERROR_BADPARM,"%s: ### bad phase1: %f => phase1 must not equal phase0 ...quitting\n",
              progname,phase1);
  }
  if (slope < 0) slope = -slope;
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

void read_stats(float ***stats, char *indir, char *instem, char *infix)
{
  int i,n;
  char fname[STRLEN];
  float *vals;

  if (MATCH(ext,"bfloat")){
    for (n=0;n<nslices;n++) {
      if (MATCH(out_ext, "w")) {
        if (infix==NULL) {
          sprintf(fname,"%s/%s-%s_%03d.bfloat",indir,instem,hemi,n);
        } else {
          sprintf(fname,"%s/%s%s-%s_%03d.bfloat",indir,instem,infix,hemi,n);
        }
      } else {
        if (infix==NULL) {
          sprintf(fname,"%s/%s_%03d.bfloat",indir,instem,n);
        } else {
          sprintf(fname,"%s/%s%s_%03d.bfloat",indir,instem,infix,n);
        }
      }
      if(readFloatImage(fname,stats[n],xnum,ynum)==-1)
        ErrorExit(ERROR_BADPARM,"%s: ### error reading float image...quitting\n",
                  progname);
    }
  } else if (MATCH(ext, "w")) {
    if (infix==NULL)
      sprintf(fname,"%s/%s-%s.w",indir,instem,hemi);
    else
      sprintf(fname,"%s/%s%s-%s.w",indir,instem,infix,hemi);
    vals = (float *)calloc(ynum,sizeof(float));             MTEST(vals);
    if(readSurfVals(fname,vals,ynum)==-1)
      ErrorExit(ERROR_BADPARM,"%s: ### error reading surface values...quitting\n",
                progname);
    for(i=0;i<ynum;i++) stats[0][i][0]=vals[i];
  }    
}

void write_stats(float ***stats, char *outdir, char *outstem, \
                                               char *prefix, char *infix)
{
  int i,n;
  char imgname[STRLEN];
  char hdrname[STRLEN];
  float *vals;

  if (MATCH(out_ext,"bfloat")) {
    for (n=0;n<nslices;n++) {
      if (MATCH(ext, "w")) {
      /* include hemi if data came from surface */
        if (prefix == NULL && infix == NULL) {
          sprintf(imgname,"%s/%s-%s_%03d.bfloat",outdir,outstem,hemi,n);
          sprintf(hdrname,"%s/%s-%s_%03d.hdr",outdir,outstem,hemi,n);
        } else if (prefix == NULL) {
          sprintf(imgname,"%s/%s%s-%s_%03d.bfloat",outdir,outstem,infix,hemi,n);
          sprintf(hdrname,"%s/%s%s-%s_%03d.hdr",outdir,outstem,infix,hemi,n);
        } else if (infix == NULL) {
          sprintf(imgname,"%s/%s-%s-%s_%03d.bfloat",outdir,prefix,outstem,hemi,n);
          sprintf(hdrname,"%s/%s-%s-%s_%03d.hdr",outdir,prefix,outstem,hemi,n);
        } else {
          sprintf(imgname,"%s/%s-%s%s-%s_%03d.bfloat",outdir,prefix,outstem,infix,hemi,n);
          sprintf(hdrname,"%s/%s-%s%s-%s_%03d.hdr",outdir,prefix,outstem,infix,hemi,n);
        }
      } else {
        if (prefix == NULL && infix == NULL) {
          sprintf(imgname,"%s/%s_%03d.bfloat",outdir,outstem,n);
          sprintf(hdrname,"%s/%s_%03d.hdr",outdir,outstem,n);
        } else if (prefix == NULL) {
          sprintf(imgname,"%s/%s%s_%03d.bfloat",outdir,outstem,infix,n);
          sprintf(hdrname,"%s/%s%s_%03d.hdr",outdir,outstem,infix,n);
        } else if (infix == NULL) {
          sprintf(imgname,"%s/%s-%s_%03d.bfloat",outdir,prefix,outstem,n);
          sprintf(hdrname,"%s/%s-%s_%03d.hdr",outdir,prefix,outstem,n);
        } else {
          sprintf(imgname,"%s/%s-%s%s_%03d.bfloat",outdir,prefix,outstem,infix,n);
          sprintf(hdrname,"%s/%s-%s%s_%03d.hdr",outdir,prefix,outstem,infix,n);
        }
      }
      if(writeFloatImage(imgname,stats[n],xnum,ynum)==-1) {
        ErrorExit(ERROR_BADFILE,"%s: ### error writing float image...quitting\n",
                  progname);
      }
      if(writeFloatHeader(hdrname,xnum,ynum,1,0)==-1) {
        ErrorExit(ERROR_BADFILE,"%s: ### error writing float header...quitting\n",
                  progname);
      }
    }
  } else if (MATCH(out_ext, "w")) {
    if (prefix == NULL && infix == NULL) {
      sprintf(imgname,"%s/%s-%s.w",outdir,outstem,hemi);
    } else if (prefix == NULL) {
      sprintf(imgname,"%s/%s%s-%s.w",outdir,outstem,infix,hemi);
    } else if (infix == NULL) {
      sprintf(imgname,"%s/%s-%s-%s.w",outdir,prefix,outstem,hemi);
    } else {
      sprintf(imgname,"%s/%s-%s%s-%s.w",outdir,prefix,outstem,infix,hemi);
    }
    vals = (float *)calloc(ynum,sizeof(float));             MTEST(vals);
    for(i=0;i<ynum;i++) vals[i]=stats[0][i][0];
    if(writeSurfVals(imgname,vals,ynum)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error writing surface values...quitting\n",
           progname);
    }
  }  
}

void write_stats_text(float ***stats, char *outdir, char *outstem, \
                                                    char *prefix, char *infix)
{
  int i,n;
  char fname[STRLEN];
  float *vals;

  if (MATCH(out_ext,"bfloat")) {
    for (n=0;n<nslices;n++) {
      if (prefix == NULL && infix == NULL) {
        sprintf(fname,"%s/%s_%03d.txt",outdir,outstem,n);
      } else if (prefix == NULL) {
        sprintf(fname,"%s/%s%s_%03d.txt",outdir,outstem,infix,n);
      } else if (infix == NULL) {
        sprintf(fname,"%s/%s-%s_%03d.txt",outdir,prefix,outstem,n);
      } else {
        sprintf(fname,"%s/%s-%s%s_%03d.txt",outdir,prefix,outstem,infix,n);
      }
      if(writeTextImage(fname,stats[n],xnum,ynum)==-1) {
        ErrorExit(ERROR_BADFILE,"%s: ### error writing float image...quitting\n",
                  progname);
      }
    }
  } else if (MATCH(out_ext, "w")) {
    if (prefix == NULL && infix == NULL) {
      sprintf(fname,"%s/%s-%s.txt",outdir,outstem,hemi);
    } else if (prefix == NULL) {
      sprintf(fname,"%s/%s%s-%s.txt",outdir,outstem,infix,hemi);
    } else if (infix == NULL) {
      sprintf(fname,"%s/%s-%s-%s.txt",outdir,prefix,outstem,hemi);
    } else {
      sprintf(fname,"%s/%s-%s%s-%s.txt",outdir,prefix,outstem,infix,hemi);
    }
    vals = (float *)calloc(ynum,sizeof(float));             MTEST(vals);
    for(i=0;i<ynum;i++) vals[i]=stats[0][i][0];
    if(writeTextVals(fname,vals,ynum)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error writing surface values...quitting\n",
                progname);
    }
  }  
}

void subphase_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float ampl, phas, real, imag;

  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++) {
      for (j=0;j<xnum;j++) {
        real = stats_r[n][i][j];
        imag = stats_i[n][i][j];
        ampl = hypot(real,imag);
        phas = atan2(imag,real) - offset*2.0*M_PI;
        real = ampl*cos(phas);
        imag = ampl*sin(phas);
        stats_r[n][i][j] = real;
        stats_i[n][i][j] = imag;
      }
    }
  }    
}

/* these two versions are equivalent */
#if 0
void revphase_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float ampl, phas, real, imag;

  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++) {
      for (j=0;j<xnum;j++) {
        real = stats_r[n][i][j];
        imag = stats_i[n][i][j];
        ampl = hypot(real,imag);
        phas = -atan2(imag,real);
        real = ampl*cos(phas);
        imag = ampl*sin(phas);
        stats_r[n][i][j] = real;
        stats_i[n][i][j] = imag;
      }
    }
  }    
}
#else
void revphase_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;

  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++) {
      for (j=0;j<xnum;j++) {
        stats_i[n][i][j] = -stats_i[n][i][j];
      }
    }
  }    
}
#endif

void truncphase_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float x,y;
  float a,a0,a1; /* amplitude */
  float p,p0,p1,corner0,corner1; /* phase */

  p0 = phase0*2.0*M_PI;
  p1 = phase1*2.0*M_PI;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    x = stats_r[n][i][j];
    y = stats_i[n][i][j];
    if (x==0 && y==0) { continue;}

    p = atan2(y,x); a = hypot(x,y); 
    if (p < 0) p+= 2.0*M_PI;
    if (p > p0) corner0 = p0 + 2.0*M_PI;
    else        corner0 = p0;
    if (p < p1) corner1 = p1 - 2.0*M_PI;
    else        corner1 = p1;

    if (((p0 < p1) && (p < p0 || p > p1)) ||
        ((p0 > p1) && (p < p0 && p > p1))) {
      a0 = a*(1.0-slope*(corner0-p)/M_PI);
      a1 = a*(1.0-slope*(p-corner1)/M_PI);
      /* a0 and a1 are amplitudes for two tails of trunc dropoff */
      if (a0 > a1) a = a0; else a = a1;
      if (a < 0) a = 0;
    }

    stats_r[n][i][j] = a*cos(p);
    stats_i[n][i][j] = a*sin(p);
  }
}

void cmplx2polar_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float x,y;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    x = stats_r[n][i][j];
    y = stats_i[n][i][j];
    stats_r[n][i][j] = hypot(x,y);
    stats_i[n][i][j] = atan2(y,x);
  }
}

void cmplx2samp_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float x,y;
  float a,p,p0,p1; /* amplitude and phase*/

  p0 = phase0*2.0*M_PI;
  p1 = phase1*2.0*M_PI;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    x = stats_r[n][i][j];
    y = stats_i[n][i][j];

    if (x==0 && y==0) {
      p = 0;
      a = 0;
    } else {
      p = atan2(y,x);
      a = hypot(x,y); 
      if (p < 0) p+= 2.0*M_PI;
      if (((p0 < p1) && (p < p0 || p > p1)) ||
          ((p0 > p1) && (p < p0 && p > p1))) {
        a = -a;
      }
    }
    stats_r[n][i][j] = a;
    stats_i[n][i][j] = p;
  }
}

void polar2cmplx_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float a,p;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    a = stats_r[n][i][j];
    p = stats_i[n][i][j];
    stats_r[n][i][j] = a*cos(p);
    stats_i[n][i][j] = a*sin(p);
  }
}

void normcxamp_stats(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float x,y,a,p;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    x = stats_r[n][i][j];
    y = stats_i[n][i][j];
    if (x==0 && y==0) {
      p = 0;
      a = 0;
    } else {
      p = atan2(y,x);
      a = 1;
    }
    stats_r[n][i][j] = a*cos(p);
    stats_i[n][i][j] = a*sin(p);
  }
}

void thresh_stats(float ***stats, float min, float max, int absflag)
{
  int i,j,n;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
    if (absflag) {
      if (fabs(stats[n][i][j]) < fabs(min) || fabs(stats[n][i][j]) > fabs(max))
        stats[n][i][j] = 0;
    } else {
      if (stats[n][i][j] < min || stats[n][i][j] > max)
        stats[n][i][j] = 0;
    }
}

void thresh_cx_stats(float ***stats_r, float ***stats_i, float min, float max)
{
  int i,j,n;
  float x,y,a;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    x = stats_r[n][i][j];
    y = stats_i[n][i][j];
    a = hypot(x,y);
    if (a < min || a > max) {
      stats_r[n][i][j] = 0;
      stats_i[n][i][j] = 0;
    }
  }
}

void mask_stats(float ***stats_r, float ***stats_i)
{
  float ***maskvals;
  int n,i,j;

  /* allocate memory for data and calculations */
  MsgPrintf("%s: starting to allocate memory for mask vals\n",progname);
  maskvals = (float ***)calloc(nslices,sizeof(float **));       MTEST(maskvals);
  for (n=0;n<nslices;n++) {
    maskvals[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*maskvals);
    for (i=0;i<ynum;i++) {
      maskvals[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**maskvals);
    }
  }
  MsgPrintf("%s: finished allocating memory for mask vals\n",progname);

  read_stats(maskvals,maskdir,maskstem,NULL);

  /* apply mask to stat values */
  for (n=0;n<nslices;n++)
    for (i=0;i<ynum;i++)
      for (j=0;j<xnum;j++) {
        if (maskvals[n][i][j]==0) {
          stats_r[n][i][j] = 0;
          if (complex) stats_i[n][i][j] = 0;
        }
      }

  MsgPrintf("%s: freeing memory for mask vals\n",progname);
  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++) {
      free(maskvals[n][i]);
    }
    free(maskvals[n]);
  }
  free(maskvals);
}

float getMax(float ***stats)
{
  int i,j,n;
  float max=-BIGFLOAT;

  for (n=0;n<nslices;n++)
    for (i=0;i<ynum;i++)
      for (j=0;j<xnum;j++)
        if (stats[n][i][j] > max)
          max = stats[n][i][j];
  return(max);
}

float getcxMax(float ***stats_r, float ***stats_i)
{
  int i,j,n;
  float x,y,a;
  float max=-BIGFLOAT;

  for (n=0;n<nslices;n++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    x = stats_r[n][i][j];
    y = stats_i[n][i][j];
    a = hypot(x,y);
    if (a > max)
      max = a;
  }
  return(max);
}

