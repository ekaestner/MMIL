/* wt2w.c: extract time points from wt file (surface time series)
      created: 06/01/05 DH
     last mod: 12/31/05 DH

   purpose:
     extracting single or multiple time points from a wt file
     (w file with multiple time points)

   input:
     wt file

   output:
     w file
*/

#include "surflib.h"

#define MINARGC 2

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char infix[STRLEN]=UNDEFSTR;
int tfirst = 0;
int tlast = 0;

void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -subj subj [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  <instem>      omit extension: <instem>-{rh,lh}.wt\n");
  printf("    -subj    <subj>        subject name\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem  [instem]   output file stem\n");
  printf("    -infix    [null]     infix added before hemi (e.g. _r)\n");
  printf("    -indir    [.]        input dir\n");
  printf("    -outdir   [.]        output dir (must exist)\n");
  printf("    -hemi     [rh]       hemisphere (rh or lh)\n");
  printf("    -tfirst   [0]        index of first time point to extract\n");
  printf("    -tlast    [0]        index of last time point to extract\n");
  printf("    -quiet               suppress messages\n");
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
      if (MATCH(argv[i],"-instem") && i+1<argc){
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
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
      if (MATCH(argv[i],"-infix") && i+1<argc) {
        strcpy(infix,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-tfirst") && i+1<argc) {
        tfirst = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-tlast") && i+1<argc) {
        tlast = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    } else
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
  }

  /* check arguments */
  MsgPrintf("%s: starting to check arguments\n",progname);
  if (MATCH(instem,UNDEFSTR))
    ErrorExit(ERROR_BADPARM,"%s: ### must supply an instem ...quitting\n",
              progname);
  if (MATCH(outstem,UNDEFSTR))
    strcpy(outstem,instem);
  if (!FileExists(indir))
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",
              progname,indir);
  if(!isadir(indir))
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
  if (!FileExists(outdir))
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  if(!isadir(outdir))
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  if (MATCH(subj,UNDEFSTR))
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh"))
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  MsgPrintf("%s: finished checking arguments\n",progname);
}

int main(int argc, char **argv)
{
  FILE  *in_fp=NULL, *out_fp=NULL;
  char tempstr[STRLEN];
  int i,t,nverts,numsurfvals,tpoints,vnum;
  float val;

  parse_args(argc,argv);

  /* get number of vertices for subject's surface */
  nverts = nverticesSurf(subj,hemi);

  /* open input file */
  if(MATCH(infix,UNDEFSTR)) {
    sprintf(tempstr,"%s/%s-%s.wt",indir,instem,hemi);
  } else {
    sprintf(tempstr,"%s/%s%s-%s.wt",indir,instem,infix,hemi);
  }  
  MsgPrintf("%s: opening input wt file %s\n",progname,tempstr);
  in_fp = fopen(tempstr,"r");
  if (in_fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,tempstr);
  /* read number of tpoints */
  fread2(&tpoints,in_fp);
  MsgPrintf("%s: %d tpoints found in %s\n",progname,tpoints,tempstr);

  /* check bounds */
  if(tfirst < 0) tfirst = 0;
  if(tfirst >= tpoints) tfirst = tpoints-1;
  if(tlast < tfirst) tlast = tfirst;
  if(tlast >= tpoints) tlast = tpoints-1;

  for(t=0;t<=tlast;t++) {
    fread3(&numsurfvals,in_fp);
/*
    MsgPrintf("%s: %d vals to read at time point %d\n",
              progname,numsurfvals,t);
*/
    if (t >= tfirst) {
      /* create new w file */
      if (MATCH(infix,UNDEFSTR)) {
        if (tpoints >= 10000) {
          sprintf(tempstr,"%s/%s-tpoint%05d-%s.w",outdir,outstem,t,hemi);
        } else if (tpoints >= 1000) {
          sprintf(tempstr,"%s/%s-tpoint%04d-%s.w",outdir,outstem,t,hemi);
        } else if (tpoints >= 100) {
          sprintf(tempstr,"%s/%s-tpoint%03d-%s.w",outdir,outstem,t,hemi);
        } else if (tpoints >= 10) {
          sprintf(tempstr,"%s/%s-tpoint%02d-%s.w",outdir,outstem,t,hemi);
        } else {
          sprintf(tempstr,"%s/%s-tpoint%d-%s.w",outdir,outstem,t,hemi);
        }
      } else {
        if (tpoints >= 10000) {
          sprintf(tempstr,"%s/%s-tpoint%05d%s-%s.w",outdir,outstem,t,infix,hemi);
        } else if (tpoints >= 1000) {
          sprintf(tempstr,"%s/%s-tpoint%04d%s-%s.w",outdir,outstem,t,infix,hemi);
        } else if (tpoints >= 100) {
          sprintf(tempstr,"%s/%s-tpoint%03d%s-%s.w",outdir,outstem,t,infix,hemi);
        } else if (tpoints >= 10) {
          sprintf(tempstr,"%s/%s-tpoint%02d%s-%s.w",outdir,outstem,t,infix,hemi);
        } else {
          sprintf(tempstr,"%s/%s-tpoint%d%s-%s.w",outdir,outstem,t,infix,hemi);
        }
      }
      MsgPrintf("%s: creating output w file %s\n",progname,tempstr);
      out_fp = fopen(tempstr,"w");
      if (out_fp==NULL)
        ErrorExit(ERROR_BADFILE,"%s: ### unable to create file %s\n",progname,tempstr);
      /* write number of tpoints at start of output file */
      fwrite2(0,out_fp);
      fwrite3(numsurfvals, out_fp);
    }
    for(i=0;i<numsurfvals;i++) {
      fread3(&vnum, in_fp);
      fread4(&val, in_fp);
/*      
      MsgPrintf("%s: vnum=%d, val=%0.4f\n",
                progname,vnum,val);
*/
      if (t >= tfirst) {
        fwrite3(vnum, out_fp);
        fwrite4(val, out_fp);
      }
    }
    if (t >= tfirst) fclose(out_fp);  
  }
  fclose(in_fp);
  
  MsgPrintf("%s: finished.\n",progname);
  exit(0);
}

