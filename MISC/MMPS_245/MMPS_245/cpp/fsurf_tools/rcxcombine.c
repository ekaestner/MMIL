/* rcxcombine.c: calculates average from real or complex
    bfloat or w (painted surface) data sets
      created: 08/11/03 DH
     last mod: 03/17/06 DH

   purpose:
     calculating multi-scan averages for complex or real data
     calculating group averages complex or real data

   input:
     bfloat volumes or single "slice" output from mri_surf2surf (not timeseries)
     w value files

   output:
     bfloat volumes (or single slice)
     w value files
*/

#include "surflib.h"

#define MINARGC 8
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"

/* global variables */
char progname[20]="rcxcombine";

/* parameter defaults */
char **input;
char outdir[STRLEN]=".";
char outstem[STRLEN]=UNDEFSTR;
char ext[STRLEN]="bfloat";
char hemi[STRLEN]="rh";
char subj[STRLEN]=UNDEFSTR;
char *real_infix=NULL;
char *imag_infix=NULL;
int *revphase, *invphase, *negvals;
int nslices=0;
int xnum=0,ynum=0,depth=0;
int N=0;
int Nrev=0;
int Ninv=0;
int Nneg=0;
int outtext = 0;
int complex = 0;
float offset=0.0;

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);
void read_stats(float ***stats, char *indir, char *instem, char *infix);
void write_stats(float ***stats, char *outdir, char *outstem, \
                                               char *prefix, char *infix);
void write_stats_text(float ***stats, char *outdir, char *outstem, \
                                                    char *prefix, char *infix);

int main(int argc, char **argv)
{
  int i,j,k,n;
  float ***dat_r=NULL;
  float ***dat_i=NULL;
  float ***avg_r=NULL;
  float ***avg_i=NULL;
  float real,imag;
  float ampl,phas;
  float N_inv;

  if (argc<MINARGC) {usage(); exit(0);}
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
  avg_r = (float ***)calloc(nslices,sizeof(float **));       MTEST(avg_r);
  for (n=0;n<nslices;n++) {
    dat_r[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_r);
    avg_r[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*avg_r);
    for (i=0;i<ynum;i++) {
      dat_r[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_r);
      avg_r[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**avg_r);
    }
  }
  if(complex) {
    dat_i = (float ***)calloc(nslices,sizeof(float **));       MTEST(dat_i);
    avg_i = (float ***)calloc(nslices,sizeof(float **));       MTEST(avg_i);
    for (n=0;n<nslices;n++) {
      dat_i[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_i);
      avg_i[n] = (float **)calloc(ynum,sizeof(float *));       MTEST(*avg_i);
      for (i=0;i<ynum;i++) {
        dat_i[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_i);
        avg_i[n][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**avg_i);
      }
    }
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* calculate sum for input datasets */
  MsgPrintf("%s: calculating sum and sum of squares...\n",progname);
  for (k=0;k<N;k++) {
    read_stats(dat_r,input[2*k],input[2*k+1],real_infix);
    if(complex)
      read_stats(dat_i,input[2*k],input[2*k+1],imag_infix);
    for (n=0;n<nslices;n++) {
      for (i=0;i<ynum;i++) {
        for (j=0;j<xnum;j++) {
          real = dat_r[n][i][j];
          if (negvals[k]) real = -real;
          if(complex) {
            imag = dat_i[n][i][j];
            if (offset!=0.0) {
              ampl = hypot(real,imag);
              phas = atan2(imag,real) - offset*2.0*M_PI;
              real = ampl*cos(phas);
              imag = ampl*sin(phas);
            }
            if (Nrev!=0) if (revphase[k]) imag = -imag;
            if (Ninv!=0) if (invphase[k]) {
              imag = -imag;
              real = -real;
            }
            avg_i[n][i][j] += imag;
          } /* if complex */
          avg_r[n][i][j] += real;
        } /* for j */
      } /* for i */
    } /* for n */
  } /* for k */

  /* calculate averages */
  MsgPrintf("%s: calculating average...\n",progname);
  N_inv = 1.0/N;
  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++)
    for (j=0;j<xnum;j++) {
      avg_r[n][i][j] *= N_inv;
      if(complex) {
        avg_i[n][i][j] *= N_inv;
      }
    }
  }

  /* write output to file */
  MsgPrintf("%s: writing output to file...\n",progname);
  write_stats(avg_r,outdir,outstem,NULL,real_infix);
  if(complex)
    write_stats(avg_i,outdir,outstem,NULL,imag_infix);
  if (outtext) {
    write_stats_text(avg_r,outdir,outstem,NULL,real_infix);
    if(complex)
      write_stats_text(avg_i,outdir,outstem,NULL,imag_infix);
  }

  exit(0);
}

void usage() {
  printf("\n");
  printf("Usage: %s -input indir1 instem1 indir2 instem2 ... \\\n",progname);
  printf("          -subj subjname\n");
  printf("         [-revphase 0 1 ... -options]\n");
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -input     subsequent dirs and stems used for average\n");
  printf("      indir      each contains a data set\n");
  printf("      instem     omit infixes, suffix: <instem>_{r,i}_000.bfloat\n");
  printf("                                    or <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -subj    <str> subject name (required only if ext = w)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -ext     <str> [bfloat]     input extension (bfloat or w)\n");
  printf("    -nslices <int> [0]          input slice count (bfloat ext only)\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -hemi    <str> [rh]         hemisphere (rh or lh) (w ext only)\n");
  printf("    -outdir  <str> [.]          output dir (must exist)\n");
  printf("    -outstem <str> [instem-avg] output file stem\n");
  printf("    -outtext                    save stats to additional text files\n");
  printf("    -negvals subsequent flags (0 or 1) indicate whether to make\n");
  printf("       values of input data negative (for non-complex data)\n");
  printf("    -quiet                      suppress messages\n");
  printf("\n");
  printf("  Optional parameters for complex data only:\n");
  printf("    -complex                   for complex data sets (real + imaginary)\n");
  printf("    -revphase  subsequent flags (0 or 1) indicate whether to reverse\n");
  printf("       phase of input data (set imag component negative)\n");
  printf("    -invphase  subsequent flags (0 or 1) indicate whether to invert\n");
  printf("       phase of input data (set real and imag components negative)\n");
  printf("       if invphase is supplied, revphase is ignored\n");
  printf("    -offset <float> [0.0]      phase offset value\n");
  printf("       subtraction of phase offset is done prior to phase reversal\n");
  printf("       use phase offset to correct for estimated hemodynamic delay\n");
  printf("    -infixes <r> <i>           real,imaginary infixes--default: %s %s\n",
         REAL_INFIX,IMAG_INFIX);
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i,j;
  FILE *fp;
  
  /* allocate memory */
  revphase = (int *)malloc(argc*sizeof(int));            MTEST(revphase);
  invphase = (int *)malloc(argc*sizeof(int));            MTEST(invphase);
  negvals  = (int *)malloc(argc*sizeof(int));            MTEST(negvals);
  input    = (char **)malloc(argc*sizeof(char *));       MTEST(input);
  for(i=0;i<argc;i++)
    input[i] = (char *)malloc(STRLEN*sizeof(char));      MTEST(*input);
  real_infix = (char *)malloc(STRLEN*sizeof(char));      MTEST(real_infix);
  imag_infix = (char *)malloc(STRLEN*sizeof(char));      MTEST(imag_infix);

  strcpy(real_infix,REAL_INFIX);
  strcpy(imag_infix,IMAG_INFIX);

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
      else if (MATCH(argv[i],"-outtext")){
        outtext = 1;
      }
      else if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      }
      else if (MATCH(argv[i],"-offset") && i+1<argc){
        offset = atof(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outtype") && i+1<argc) {
        i++;
        /* no longer used */
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
      else if (MATCH(argv[i],"-complex")){
        complex = 1;
      }
      else if (MATCH(argv[i],"-revphase")){
        MsgPrintf("%s: searching for revphase flags...\n",progname);
        if (++i>argc) {
          ErrorExit(ERROR_BADPARM,"%s: ### no flags following -revphase ...quitting\n",
                    progname);
        }
        for (j=i;j<argc;j++) {
          if (argv[j][0]=='-') break;
          revphase[Nrev]=atoi(argv[j]);
          Nrev++;
        }
        i=j-1;
      }
      else if (MATCH(argv[i],"-invphase")){
        MsgPrintf("%s: searching for invphase flags...\n",progname);
        if (++i>argc) {
          ErrorExit(ERROR_BADPARM,"%s: ### no flags following -invphase ...quitting\n",
                    progname);
        }
        for (j=i;j<argc;j++) {
          if (argv[j][0]=='-') break;
          invphase[Ninv]=atoi(argv[j]);
          Ninv++;
        }
        i=j-1;
      }
      else if (MATCH(argv[i],"-negvals")){
        MsgPrintf("%s: searching for negvals flags...\n",progname);
        if (++i>argc) {
          ErrorExit(ERROR_BADPARM,"%s: ### no flags following -negvals ...quitting\n",
                    progname);
        }
        for (j=i;j<argc;j++) {
          if (argv[j][0]=='-') break;
          negvals[Ninv]=atoi(argv[j]);
          Ninv++;
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
  if (MATCH(outstem,UNDEFSTR)) sprintf(outstem,"avg");
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
  if (Ninv) Nrev=0;
  if (Nrev!=0 && Nrev!=N) {
    ErrorExit(ERROR_BADPARM,"%s: ### number of revphase flags (%d) must equal number of datasets (%d) ...quitting\n",
              progname,Nrev,N);
  }
  for (i=0;i<Nrev;i++) {
    if (revphase[i] != 0 && revphase[i] != 1) {
      ErrorExit(ERROR_BADPARM,"%s: ### revphase flags (%d) must be 0 or 1 ...quitting\n",
                progname, revphase[i]);
    }
  }
  if (Ninv!=0 && Ninv!=N) {
    ErrorExit(ERROR_BADPARM,"%s: ### number of invphase flags (%d) must equal number of datasets (%d) ...quitting\n",
              progname,Ninv,N);
  }
  for (i=0;i<Ninv;i++) {
    if (invphase[i] != 0 && invphase[i] != 1) {
      ErrorExit(ERROR_BADPARM,"%s: ### invphase flags (%d) must be 0 or 1 ...quitting\n",
                progname, invphase[i]);
    }
  }
  if (Nneg!=0 && Nneg!=N) {
    ErrorExit(ERROR_BADPARM,"%s: ### number of negvals flags (%d) must equal number of datasets (%d) ...quitting\n",
              progname,Nneg,N);
  }
  for (i=0;i<Nneg;i++) {
    if (negvals[i] != 0 && negvals[i] != 1) {
      ErrorExit(ERROR_BADPARM,"%s: ### negvals flags (%d) must be 0 or 1 ...quitting\n",
                progname, negvals[i]);
    }
  }
  if (!complex) {
    free(real_infix);
    real_infix = NULL;
    free(imag_infix);
    imag_infix = NULL;
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

