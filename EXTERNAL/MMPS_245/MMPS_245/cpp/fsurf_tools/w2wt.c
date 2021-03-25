/* w2wt.c: combine multiple w time points into a wt file (surface time series)
      created: 06/07/05 DH
     last mod: 06/07/05 DH

   purpose:
     combining multiple time points into a wt file

   input:
     w file

   output:
     wt file
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
int tlast=0;

void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -subj subj [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  <instem>      omit extension: <instem>-{rh,lh}.w\n");
  printf("    -subj    <subj>        subject name\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -tlast    [0]        index of last time point to find\n");
  printf("    -outstem  [instem]   output file stem\n");
  printf("    -indir    [.]        input dir\n");
  printf("    -outdir   [.]        output dir (must exist)\n");
  printf("    -hemi     [rh]       hemisphere (rh or lh)\n");
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
  if (tlast < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### %s: tlast must be >= 0\n",progname,hemi);
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

  /* todo: use dirent to search for files matching instem */

  sprintf(tempstr,"%s/%s-%s.wt",outdir,outstem,hemi);
  MsgPrintf("%s: creating output wt file %s\n",progname,tempstr);
  out_fp = fopen(tempstr,"w");
  if (out_fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### unable to create file %s\n",progname,tempstr);
  fwrite2(tlast,out_fp);

  /* open input files and write to wt file */
  for(t=0;t<=tlast;t++) {
    if (tlast >= 10000) {
      sprintf(tempstr,"%s/%s-tpoint%05d-%s.w",indir,instem,t,hemi);
    } else if (tlast >= 1000) {
      sprintf(tempstr,"%s/%s-tpoint%04d-%s.w",indir,instem,t,hemi);
    } else if (tlast >= 100) {
      sprintf(tempstr,"%s/%s-tpoint%03d-%s.w",indir,instem,t,hemi);
    } else if (tlast >= 10) {
      sprintf(tempstr,"%s/%s-tpoint%02d-%s.w",indir,instem,t,hemi);
    } else {
      sprintf(tempstr,"%s/%s-tpoint%d-%s.w",indir,instem,t,hemi);
    }
    MsgPrintf("%s: opening input wt file %s\n",progname,tempstr);
    in_fp = fopen(tempstr,"r");
    if (in_fp==NULL)
      ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,tempstr);
    /* read number of tpoints */
    fread2(&tpoints,in_fp);
    /* read number of non-zero values */
    fread3(&numsurfvals,in_fp);
    fwrite3(numsurfvals,out_fp);

    for(i=0;i<numsurfvals;i++) {
      fread3(&vnum, in_fp);
      fread4(&val, in_fp);
/*      
      MsgPrintf("%s: vnum=%d, val=%0.4f\n",
                progname,vnum,val);
*/
      fwrite3(vnum, out_fp);
      fwrite4(val, out_fp);
    }
    fclose(in_fp);
  }
  fclose(out_fp);
  
  MsgPrintf("%s: finished.\n",progname);
  exit(0);
}

