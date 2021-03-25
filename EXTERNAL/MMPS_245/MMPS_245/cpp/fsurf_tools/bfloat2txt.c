/* bfloat2txt.c: convert from bfloat to txt
      created: 09/02/03 DH
     last mod: 09/22/03 DH

   purpose:
     converting binary bfloat images into human readable text files

   input:
     bfloat volumes (not timeseries)

   output:
     text files with .txt extension

   credits:
     this program was created with code stolen from:
       surfer's  phasecombine.c  (arguement parsing, file i/o, byte swapping)
*/

#include "surflib.h"

#define MINARGC 2

/* global variables */
char progname[20]="bfloat2txt";

/* parameter defaults */
int nslices=0;
int xnum=0,ynum=0,depth=0;
char instem[STRLEN]="unknown";
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
int outtext = 0;

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);

int main(int argc, char **argv)
{
  int i,n;
  float **stats;
  char fname[STRLEN];

  if (argc<MINARGC) {usage(); exit(0);}
  parse_args(argc,argv);

  sprintf(fname,"%s/%s_000.hdr",indir,instem);
  if(readFloatHeader(fname,&xnum,&ynum,&depth)==-1) {
    ErrorExit(ERROR_BADPARM,"%s: ### error reading header file...quitting\n",
           progname);
  }

  /* allocate memory for data and calculations */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  stats = (float **)calloc(ynum,sizeof(float *));         MTEST(stats);
  for (i=0;i<ynum;i++) {
    stats[i] = (float *)calloc(xnum,sizeof(float));       MTEST(*stats);
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* convert images */
  for (n=0;n<nslices;n++) {
    sprintf(fname,"%s/%s_%03d.bfloat",indir,instem,n);
    if(readFloatImage(fname,stats,xnum,ynum)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error reading float image...quitting\n",
           progname);
    }

    sprintf(fname,"%s/%s_%03d.txt",outdir,instem,n);
    if(writeTextImage(fname,stats,xnum,ynum)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error writing text image...quitting\n",
           progname);
    }
  }

  exit(0);
}

void usage()
{
  printf("\n");
  printf("Usage: %s instem [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    instem     omit infixes, suffix: <instem>_000.bfloat\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -indir   <str>  [.]        input dir\n");
  printf("    -outdir  <str>  [.]        output dir (must exist)\n");
  printf("    -nslices <int>  [0]        input slice count\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -quiet                     suppress messages\n");
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i;
  FILE *fp;
  
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
  if (MATCH(instem,"unknown")) {
    ErrorExit(ERROR_BADPARM,"%s: ### must supply an instem ...quitting\n",
              progname);
  }
  fp = fopen(indir,"r");
  if (fp==NULL) {
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",
              progname,indir);
  } else fclose(fp);
  if(!isadir(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
  }
  fp = fopen(outdir,"r");
  if (fp==NULL) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  } else fclose(fp);
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  i = getNumBfloats(indir,instem,NULL);
  if(i<=0) {
    ErrorExit(ERROR_NOFILE,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
           progname,indir,instem);
  }
  if(nslices==0 || nslices>i) {
    nslices = i;
    MsgPrintf("%s: %d slices found in %s/ with stem %s\n",progname,i,indir,instem);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

