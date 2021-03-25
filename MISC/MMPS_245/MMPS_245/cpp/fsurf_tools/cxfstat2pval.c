/* cxfstat2pval.c: calculate pvals from complex f-stats
      created: 05/01/06 DH
     last mod: 05/01/06 DH

   purpose:
     calculating pvals from complex f-stats (see cxfstat)

   input:
     either w file or bfloat (including 1D surface bfloats)
     
   output:
     p-value
     
     output is in same format as input (w files or bfloats)
*/

#include "surflib.h"
#include "dhstatlib.h"

#define MINARGC 6

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
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

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);
void read_stats(float ***stats, char *indir, char *instem, char *infix);
void write_stats(float ***stats, char *outdir, char *outstem, \
                                               char *prefix, char *infix);

void usage() {
  printf("\n");
  printf("Usage: %s -instem instem [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    instem   omit infixes, suffix: <instem>_000.bfloat\n");
  printf("                                or <instem>-{rh,lh}.w\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -ext     <str>  [w]        input extension (w or bfloat)\n");
  printf("    -subj    <str>  [ico]      subject name (w ext only)\n");
  printf("    -hemi    <str>  [rh]       hemisphere (rh or lh) (w ext only)\n");
  printf("    -nslices <int>  [0]        input slice count (bfloat ext only)\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -N       <int>  [10]       number of datasets used to calculate fratio\n");
  printf("    -outdir  <str>  [.]        output dir (must exist)\n");
  printf("    -outstem <str>  [ ]        output file stem\n");
  printf("       will output outstem-pval\n");
  printf("       pval is -log10(significance p-value)\n");
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
      if ((MATCH(argv[i],"-nslices") || MATCH(argv[i],"-n")) && i+1<argc) {
        nslices = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-instem") && i+1<argc){
        strcpy(instem,argv[i+1]); i++;
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
      else if (MATCH(argv[i],"-ext") && i+1<argc){
        strcpy(ext,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-N") && i+1<argc){
        N = atoi(argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
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
  if (!MATCH(ext,"bfloat") && !MATCH(ext,"w")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: extension not supported\n",
              progname,ext);
  }

  if (MATCH(ext,"w")) {
    nslices = 1;
  } else {
    i = getNumBfloats(indir,instem,NULL);
    if(i<=0) {
      ErrorExit(ERROR_NOFILE,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
               progname,indir,instem);
    }
    if(nslices==0 || nslices>i) {
      nslices = i;
      MsgPrintf("%s: %d slices found in %s/ with stem %s\n",progname,i,indir,instem);
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

int main(int argc, char **argv)
{
  int i,j,n;
  float ***fratio=NULL;
  float ***pval=NULL;
  int dof1,dof2;
  char tempstr[STRLEN];

  parse_args(argc,argv);

  if (MATCH(ext,"bfloat")) {
    sprintf(tempstr,"%s/%s-%s_000.hdr",indir,instem,hemi);
    if(readFloatHeader(tempstr,&xnum,&ynum,&depth)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error reading header file...quitting\n",
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
  fratio = (float ***)calloc(nslices,sizeof(float **));      MTEST(fratio);
  pval = (float ***)calloc(nslices,sizeof(float **));        MTEST(pval);
  for (n=0;n<nslices;n++) {
    fratio[n] = (float **)calloc(ynum,sizeof(float *));      MTEST(*fratio);
    pval[n] = (float **)calloc(ynum,sizeof(float *));        MTEST(*pval);
    for (i=0;i<ynum;i++) {
      fratio[n][i] = (float *)calloc(xnum,sizeof(float));    MTEST(**fratio);
      pval[n][i] = (float *)calloc(xnum,sizeof(float));      MTEST(**pval);
    }
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* load data */
  MsgPrintf("%s: loading f-ratios...\n",progname);
  read_stats(fratio,indir,instem,NULL);

  /* calculate f-ratio */
  MsgPrintf("%s: calculating pvals...\n",progname);
  dof1 = 2;
  dof2 = 2*N;
  for (n=0;n<nslices;n++) {
    for (i=0;i<ynum;i++) {
      for (j=0;j<xnum;j++) {
        /* calcuate pval */
        if(fratio[n][i][j]<=0) {
          pval[n][i][j]=0;
        } else {
          pval[n][i][j] = sigf(fratio[n][i][j],dof1,dof2);
          if (pval[n][i][j]>1) pval[n][i][j]=1;
          else if (pval[n][i][j]<=0) pval[n][i][j]=1e-40;
          pval[n][i][j] = -log10(pval[n][i][j]);
        }
      }
    }
  }

  /* write output to file */
  MsgPrintf("%s: writing output to file...\n",progname);
  if(MATCH(outstem,UNDEFSTR))
    sprintf(tempstr,"pval");
  else
    sprintf(tempstr,"%s-pval",outstem);
  write_stats(pval,outdir,tempstr,NULL,NULL);

  exit(0);
}
