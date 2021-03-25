/* cxfstat.c: calculate vector average and f-stat from complex datasets
      created: 03/16/06 DH
     last mod: 03/16/06 DH

   purpose:
     calculating group averages and f-stats for complex (phase-encoded) data

   input:
     multi-subject complex (real and imaginary) data
       these should be the raw fourier components for the stimulus frequency
       sampled onto the average ico sphere

     these can be either w files or bfloats (including 1D surface bfloats)
     
   output:
     complex (real and imaginary) average
     complex (real and imaginary) stdev
     f-stat
     p-value
     
     output is in same format as input (w files or bfloats)
*/

#include "surflib.h"
#include "dhstatlib.h"

#define MINARGC 6

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char **input;
char outdir[STRLEN]=".";
char outstem[STRLEN]=UNDEFSTR;
char ext[STRLEN]="w";
char hemi[STRLEN]="rh";
char subj[STRLEN]="ico";
char real_infix[STRLEN]="_r";
char imag_infix[STRLEN]="_i";
int nslices=0;
int xnum=0,ynum=0,depth=0;
int N=0;
float offset=0.0;
int zeromeanflag=0;
int ignoreampflag=0;

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);
void read_stats(float ***stats, char *indir, char *instem, char *infix);
void write_stats(float ***stats, char *outdir, char *outstem, \
                                               char *prefix, char *infix);
void write_stats_text(float ***stats, char *outdir, char *outstem, \
                                                    char *prefix, char *infix);

void usage() {
  printf("\n");
  printf("Usage: %s -input indir1 instem1 indir2 instem2 [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -input     subsequent dirs and stems used for average or one-sample ttest\n");
  printf("      indir      each contains a data set\n");
  printf("      instem     omit infixes, suffix: <instem>_{r,i}_000.bfloat\n");
  printf("                                    or <instem>_{r,i}-{rh,lh}.w\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -ext     <str>  [w]        input extension (w or bfloat)\n");
  printf("    -subj    <str>  [ico]      subject name (w ext only)\n");
  printf("    -hemi    <str>  [rh]       hemisphere (rh or lh) (w ext only)\n");
  printf("    -nslices <int>  [0]        input slice count (bfloat ext only)\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -outdir  <str>  [.]        output dir (must exist)\n");
  printf("    -outstem <str>  [ ]        output file stem\n");
  printf("       will output outstem-avg, outstem-stdev, outstem-fratio,\n");
  printf("         and outstem-pval\n");
  printf("       avg and stdev are complex, fratio and pval are real\n");
  printf("       pval is -log10(significance p-value)\n");
  printf("    -zeromean                  assume zero mean\n");
  printf("    -ignoreamp                 ignore amplitude (norm to 1)\n");
  printf("    -infixes <r> <i>           real,imaginary infixes--default: %s %s\n",
         real_infix,imag_infix);
  printf("    -quiet                     suppress messages\n");
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i,j;
  char tempstr[STRLEN];
  char *pch;
  FILE *fp;

  progname = argv[0];
  /* strip off path */
  pch = strtok(progname,"/");
  while (pch!=NULL) {
    strcpy(tempstr,pch);
    pch = strtok (NULL,"/");
  }
  strcpy(progname,tempstr);

  if (argc<MINARGC) {usage(); exit(0);}
  
  /* allocate memory */
  input = (char **)malloc(argc*sizeof(char *));          MTEST(input);
  for(i=0;i<argc;i++)
    input[i] = (char *)malloc(STRLEN*sizeof(char));      MTEST(*input);

  /* parse arguments */
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if ((MATCH(argv[i],"-nslices") || MATCH(argv[i],"-n")) && i+1<argc) {
        nslices = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-ext") && i+1<argc){
        strcpy(ext,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      }
      else if (MATCH(argv[i],"-zeromean")){
        zeromeanflag=1;
      }
      else if (MATCH(argv[i],"-ignoreamp")){
        ignoreampflag=1;
      }
      else if (MATCH(argv[i],"-offset") && i+1<argc){
        offset = atof(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-input")){
        if (++i>argc) {
          ErrorExit(ERROR_BADPARM,"%s: ### no datasets following -input ...quitting\n",
                    progname);
        }
        for (j=i;j<argc;j++) {
          if (argv[j][0]=='-') break;
          strcpy(input[N],argv[j]);
          N++;
        }
        i=j-1;
      }
      else if ((MATCH(argv[i],"-infixes")) && i+2<argc){
        strcpy(real_infix,argv[i+1]); strcpy(imag_infix,argv[i+2]); i+=2;
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
  if (N==0) {
    ErrorExit(ERROR_BADPARM,"%s: ### no input datasets specified ...quitting\n",
              progname);
  } 
  if (N%2) {
    ErrorExit(ERROR_BADPARM,"%s: ### odd input dirs+formats count ...quitting\n",
              progname);
  } 
  N /= 2;
  fp = fopen(outdir,"r");
  if (fp==NULL) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  } else fclose(fp);
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  if (!MATCH(ext,"bfloat") && !MATCH(ext,"w")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: extension not supported\n",
              progname,ext);
  }


  if (MATCH(ext,"w")) {
    nslices = 1;
  } else {
    i = getNumBfloats(input[0],input[1],real_infix);
    if(i<=0) {
      if(real_infix==NULL)
        ErrorExit(ERROR_NOFILE,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
               progname,input[0],input[1]);
      else
        ErrorExit(ERROR_NOFILE,"%s: ### no bfloats found in %s/ with stem %s%s ...quitting\n",
               progname,input[0],input[1],real_infix);
    }
    if(nslices==0 || nslices>i) {
      nslices = i;
      if(real_infix==NULL)
        MsgPrintf("%s: %d slices found in %s/ with stem %s\n",progname,i,input[0],input[1]);
      else
        MsgPrintf("%s: %d slices found in %s/ with stem %s%s\n",progname,i,input[0],input[1],real_infix);
    }
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",
              progname,hemi);
  }
  if (MATCH(subj,UNDEFSTR) && MATCH(ext,"w")) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

void read_stats(float ***stats, char *indir, char *instem, char *infix)
{
  int i,n;
  char fname[STRLEN];
  float *vals;

  if (MATCH(ext,"bfloat")){
    for (n=0;n<nslices;n++) {
      if (infix==NULL)
        sprintf(fname,"%s/%s_%03d.bfloat",indir,instem,n);
      else
        sprintf(fname,"%s/%s%s_%03d.bfloat",indir,instem,infix,n);
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

  if (MATCH(ext,"bfloat")) {
    for (n=0;n<nslices;n++) {
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
      if(writeFloatImage(imgname,stats[n],xnum,ynum)==-1) {
        ErrorExit(ERROR_BADPARM,"%s: ### error writing float image...quitting\n",
                  progname);
      }
      if(writeFloatHeader(hdrname,xnum,ynum,1,0)==-1) {
        ErrorExit(ERROR_BADPARM,"%s: ### error writing float header...quitting\n",
                  progname);
      }
    }
  } else if (MATCH(ext, "w")) {
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
      ErrorExit(ERROR_BADPARM,"%s: ### error writing surface values...quitting\n",
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

  if (MATCH(ext,"bfloat")) {
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
        ErrorExit(ERROR_BADPARM,"%s: ### error writing float image...quitting\n",
                  progname);
      }
    }
  } else if (MATCH(ext, "w")) {
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
      ErrorExit(ERROR_BADPARM,"%s: ### error writing surface values...quitting\n",
                progname);
    }
  }  
}


int main(int argc, char **argv)
{
  int i,j,k,n;
  float ***dat_r=NULL;
  float ***dat_i=NULL;
  float ***sum_r=NULL;
  float ***sum_i=NULL;
  float ***sumsq_r=NULL;
  float ***sumsq_i=NULL;
  float ***fratio=NULL;
  float ***pval=NULL;
  float N_inv, Nm1_inv,phase;
  int dof1,dof2;
  double numer,denom;
  char tempstr[STRLEN];

  parse_args(argc,argv);

  if (MATCH(ext,"bfloat")) {
    if(readFloatHeaders(N,input,real_infix,&xnum,&ynum,&depth)==-1) {
      ErrorExit(ERROR_BADFILE, "%s: ### error reading header file...quitting\n",
             progname);
    }
  } else if (MATCH(ext,"w")) {
    ynum = nverticesSurf(subj,hemi);
    if(ynum==-1) {
      ErrorExit(ERROR_BADFILE, "%s: ### error reading number of vertices...quitting\n",
             progname);
    }
    xnum = 1;
  }

  /* allocate memory for data and calculations */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  dat_r = (float ***)calloc(nslices,sizeof(float **));       MTEST(dat_r);
  sum_r = (float ***)calloc(nslices,sizeof(float **));       MTEST(sum_r);
  sumsq_r = (float ***)calloc(nslices,sizeof(float **));     MTEST(sumsq_r);
  fratio = (float ***)calloc(nslices,sizeof(float **));      MTEST(fratio);
  pval = (float ***)calloc(nslices,sizeof(float **));        MTEST(pval);
  for (n=0;n<nslices;n++) {
    dat_r[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_r);
    sum_r[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*sum_r);
    sumsq_r[n] = (float **)calloc(ynum,sizeof(float *));     MTEST(*sumsq_r);
    fratio[n] = (float **)calloc(ynum,sizeof(float *));      MTEST(*fratio);
    pval[n] = (float **)calloc(ynum,sizeof(float *));        MTEST(*pval);
    for (i=0;i<ynum;i++) {
      dat_r[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_r);
      sum_r[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**sum_r);
      sumsq_r[n][i] = (float *)calloc(xnum,sizeof(float));   MTEST(**sumsq_r);
      fratio[n][i] = (float *)calloc(xnum,sizeof(float));    MTEST(**fratio);
      pval[n][i] = (float *)calloc(xnum,sizeof(float));      MTEST(**pval);
    }
  }
  dat_i = (float ***)calloc(nslices,sizeof(float **));       MTEST(dat_i);
  sum_i = (float ***)calloc(nslices,sizeof(float **));       MTEST(sum_i);
  sumsq_i = (float ***)calloc(nslices,sizeof(float **));     MTEST(sumsq_i);
  for (n=0;n<nslices;n++) {
    dat_i[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_i);
    sum_i[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*sum_i);
    sumsq_i[n] = (float **)calloc(ynum,sizeof(float *));     MTEST(*sumsq_i);
    for (i=0;i<ynum;i++) {
      dat_i[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_i);
      sum_i[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**sum_i);
      sumsq_i[n][i] = (float *)calloc(xnum,sizeof(float));   MTEST(**sumsq_i);
    }
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* calculate sum and sum of squares for input datasets */
  MsgPrintf("%s: calculating sum and sum of squares...\n",progname);
  for (k=0;k<N;k++) {
    read_stats(dat_r,input[2*k],input[2*k+1],real_infix);
    read_stats(dat_i,input[2*k],input[2*k+1],imag_infix);
    for (n=0;n<nslices;n++) {
      for (i=0;i<ynum;i++) {
        for (j=0;j<xnum;j++) {
          /* normalize amplitudes to 1 */
          if(ignoreampflag) {
            phase = atan2(dat_i[n][i][j],dat_r[n][i][j]);
            dat_r[n][i][j] = cos(phase); 
            dat_i[n][i][j] = sin(phase);
          }
          sum_r[n][i][j] += dat_r[n][i][j];
          sumsq_r[n][i][j] += dat_r[n][i][j]*dat_r[n][i][j];
          sum_i[n][i][j] += dat_i[n][i][j];
          sumsq_i[n][i][j] += dat_i[n][i][j]*dat_i[n][i][j];
        } /* for j */
      } /* for i */
    } /* for n */
  } /* for k */

  /* calculate f-ratio */
  MsgPrintf("%s: calculating average and f-ratio...\n",progname);
  N_inv = 1.0/N;
  Nm1_inv = 1.0/(N-1);
  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++) {
      for (j=0;j<xnum;j++) {
        /* calcuate average */
        sum_r[n][i][j] *= N_inv;
        sum_i[n][i][j] *= N_inv;
        
        /* calcuate variance */
        if(zeromeanflag) {
          /* null hypothesis is zero mean, so could calculate 
             variance assuming mean is zero */
          sumsq_r[n][i][j] *= N_inv;
          sumsq_i[n][i][j] *= N_inv;
          dof1 = 2;
          dof2 = 2*N;
        } else {
          sumsq_r[n][i][j] = sumsq_r[n][i][j] - (double)N*sum_r[n][i][j]*sum_r[n][i][j];
          sumsq_i[n][i][j] = sumsq_i[n][i][j] - (double)N*sum_i[n][i][j]*sum_i[n][i][j];
          sumsq_r[n][i][j] *= Nm1_inv;
          sumsq_i[n][i][j] *= Nm1_inv;
          dof1 = 2;
          dof2 = 2*N-2;
        }

        /* numer equals sum of squared averages */
        numer = sum_r[n][i][j]*sum_r[n][i][j] +
                sum_i[n][i][j]*sum_i[n][i][j];

        /* denom equals sum of squared variance */
        denom = (sumsq_r[n][i][j] + sumsq_i[n][i][j]);

        /* fratio is (X2/dof1)/(X2/dof2) with X2 a chi-squared stat
           sum of squares is a chi-squared stat
           dof is the number of elements in the sum
           dof1 = 2 for real and imag
           dof2 = 2*N
           fratio = (numer/2)/(denom/2N) */
        numer /= dof1;
        denom /= dof2;
        fratio[n][i][j] = numer/denom;

        /* calcuate pval */
        pval[n][i][j] = sigf(fratio[n][i][j],dof1,dof2);
        if (pval[n][i][j]>1) pval[n][i][j]=1;
        if (pval[n][i][j]<=0) pval[n][i][j]=1e-40;
        pval[n][i][j] = -log10(pval[n][i][j]);

        /* calcuate stdev */
        sumsq_r[n][i][j] = (sumsq_r[n][i][j] > 0.0) ? sqrt(sumsq_r[n][i][j]) : 0.0;
        sumsq_i[n][i][j] = (sumsq_r[n][i][j] > 0.0) ? sqrt(sumsq_r[n][i][j]) : 0.0;
      }
    }
  }

  /* write output to file */
  MsgPrintf("%s: writing output to file...\n",progname);
  if(MATCH(outstem,UNDEFSTR))
    sprintf(tempstr,"avg");
  else
    sprintf(tempstr,"%s-avg",outstem);
  write_stats(sum_r,outdir,tempstr,NULL,real_infix);
  write_stats(sum_i,outdir,tempstr,NULL,imag_infix);
  if(MATCH(outstem,UNDEFSTR))
    sprintf(tempstr,"stdev");
  else
    sprintf(tempstr,"%s-stdev",outstem);
  write_stats(sumsq_r,outdir,tempstr,NULL,real_infix);
  write_stats(sumsq_i,outdir,tempstr,NULL,imag_infix);
  if(MATCH(outstem,UNDEFSTR))
    sprintf(tempstr,"fratio");
  else
    sprintf(tempstr,"%s-fratio",outstem);
  write_stats(fratio,outdir,tempstr,NULL,NULL);
  if(MATCH(outstem,UNDEFSTR))
    sprintf(tempstr,"pval");
  else
    sprintf(tempstr,"%s-pval",outstem);
  write_stats(pval,outdir,tempstr,NULL,NULL);

  exit(0);
}
