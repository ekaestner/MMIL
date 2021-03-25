/* wmask.c: apply mask to surface stats
      created: 05/18/04 DH
     last mod: 03/09/05 DH

   purpose:
     thresholding surface stats using a mask file

   input:
     stats file and mask file (both w files)

   output:
     masked stats file (w file)
*/

#include "surflib.h"

#define MINARGC 7

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char subj[STRLEN]=UNDEFSTR;
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char maskstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char hemi[STRLEN]="rh";
int overwrite = 0;


void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -maskstem maskstem -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem   instem                omit extension, hemi: <instem>-rh.w\n");
  printf("    -maskstem maskstem              omit extension, hemi: <instem>-rh.w\n");
  printf("    -subj     subjname              subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem  outstem [instem-mask] omit extension, hemi: <outstem>-rh.w\n");
  printf("    -indir    indir   [.]           input dir\n");
  printf("    -outdir   outdir  [.]           output dir\n");
  printf("    -hemi     hemi    [rh]          hemisphere (rh or lh)\n");
  printf("    -overwrite                      force overwrite of existing output files\n");
  printf("    -quiet                          suppress messages\n");
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
      if (MATCH(argv[i],"-maskstem") && i+1<argc) {
        strcpy(maskstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-overwrite")){
        overwrite = 1;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    } else {
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    }
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### must supply a subject name ...quitting\n",
              progname);
  }
  if (MATCH(instem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### must supply an instem ...quitting\n",
              progname);
  }
  if (MATCH(maskstem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### must supply an maskstem ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) sprintf(outstem,"%s-mask",instem);
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
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  sprintf(tempstr,"%s/%s-%s.w",outdir,outstem,hemi);
  if(FileExists(tempstr)) {
    if(overwrite) {
      MsgPrintf("%s: output file %s exists, will overwrite\n",progname,tempstr);
    } else {
      ErrorExit(ERROR_BADPARM,"%s: output file %s exists\n      (use -overwrite or supply a different outstem)  ...quitting\n",
                progname,tempstr);
    }
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}


int main(int argc, char **argv)
{
  int i;
  float *statvals;
  float *maskvals;
  char statname[STRLEN];
  char maskname[STRLEN];
  char outname[STRLEN];
  int nverts;

  parse_args(argc,argv);

  /* set file names */
  sprintf(statname,"%s/%s-%s.w",indir,instem,hemi);
  sprintf(maskname,"%s/%s-%s.w",indir,maskstem,hemi);
  sprintf(outname,"%s/%s-%s.w",outdir,outstem,hemi);

  nverts = nverticesSurf(subj,hemi);
  if(nverts==-1) {
    ErrorExit(ERROR_BADFILE,"%s: ### error reading number of vertices...quitting\n",
           progname);
  }

  /* allocate memory */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  statvals = (float *)calloc(nverts,sizeof(float)); MTEST(statvals);
  maskvals = (float *)calloc(nverts,sizeof(float)); MTEST(maskvals);
  MsgPrintf("%s: finished allocating memory\n",progname);


  /* load values */
  if(readSurfVals(statname,statvals,nverts)==-1) {
    ErrorExit(ERROR_BADFILE,"%s: ### error reading surface stat values...quitting\n",
              progname);
  }
  if(readSurfVals(maskname,maskvals,nverts)==-1) {
    ErrorExit(ERROR_BADFILE,"%s: ### error reading surface mask values...quitting\n",
              progname);
  }

  /* apply mask to stat values */
  for(i=0;i<nverts;i++)
   if (maskvals[i]==0)
     statvals[i]=0;

  /* write out file */
  if(writeSurfVals(outname,statvals,nverts)==-1) {
    ErrorExit(ERROR_BADFILE,"%s: ### error writing surface values...quitting\n",
         progname);
  }

  exit(0);
}
