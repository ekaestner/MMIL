/* stc2wt.c: convert from stc file (source time course) to wt file (surface time series)
      created: 06/01/05 DH
     last mod: 06/01/05 DH

   purpose:
     converting stc files (containing mne estimated MEG/EEG source time series) to
     wt (w file with multiple time points)

   input:
     stc file

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
int swapflag = 0;

void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -subj subj [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  <instem>      omit extension: <instem>-{rh,lh}.stc\n");
  printf("    -subj    <subj>        subject name\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem  [instem]   output file stem\n");
  printf("    -indir    [.]        input dir\n");
  printf("    -outdir   [.]        output dir (must exist)\n");
  printf("    -hemi     [rh]       hemisphere (rh or lh)\n");
  printf("    -swap                swap bytes after reading from file\n");
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
      if (MATCH(argv[i],"-swap")){
        swapflag=1;
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
  int i,vnum,t,nverts,num_source_vertices,tpoints;
  char tempstr[STRLEN];
  FILE  *in_fp, *out_fp, *info_fp;
  float val, starttime, sample_period;
  int *source_vertices=NULL;
  unsigned int temp_ui;

  parse_args(argc,argv);

  /* get number of vertices for subject's surface */
  nverts = nverticesSurf(subj,hemi);

  /* open input file */
  sprintf(tempstr,"%s/%s-%s.stc",indir,instem,hemi);
  MsgPrintf("%s: opening input stc file %s\n",progname,tempstr);
  in_fp = fopen(tempstr,"r");
  if (in_fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,tempstr);

  /* read info from start of input file */
  starttime=freadFloat(in_fp);
  sample_period=freadFloat(in_fp);
  if (swapflag) {
    starttime = swapFloat(starttime);
    sample_period = swapFloat(sample_period);
  }
  MsgPrintf("%s: start time (ms) = %0.5f\n",progname,starttime);
  MsgPrintf("%s: sample period (ms) = %0.5f\n",progname,sample_period);

  temp_ui = freadUInt(in_fp);
  if (swapflag) temp_ui = swapUInt(temp_ui);
  num_source_vertices = (int)temp_ui;;
  MsgPrintf("%s: # of sources (vertices) = %d\n",progname,num_source_vertices);

  if (num_source_vertices <= 0)
    ErrorExit(ERROR_BADPARM,"%s: ### bad num_source_vertices (%d)\n",
              progname,num_source_vertices);
  else if (num_source_vertices > nverts)
    ErrorExit(ERROR_BADPARM,"%s: ### num_source_vertices (%d) > nvertices (%d)\n",
              progname,num_source_vertices,nverts);
  else {
    /* allocate space for source_vertces */
    source_vertices = (int *)calloc(num_source_vertices,sizeof(int)); MTEST(source_vertices);
  }
  for (i=0;i<num_source_vertices;i++) {
    temp_ui = freadUInt(in_fp);
    if (swapflag) temp_ui = swapUInt(temp_ui);
    source_vertices[i]=(int)temp_ui;;
  }
  temp_ui = freadUInt(in_fp);
  if (swapflag) temp_ui = swapUInt(temp_ui);
  tpoints = (int)temp_ui;;
  MsgPrintf("%s: # of time points = %d\n",progname,tpoints);

  /* write extra info to info file */
  sprintf(tempstr,"%s/%s-%s.info",outdir,outstem,hemi);
  MsgPrintf("%s: opening output info file %s\n",progname,tempstr);
  info_fp = fopen(tempstr,"w");
  if (info_fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### unable to create file %s\n",progname,tempstr);
  MsgPrintf("%s: writing info file\n",progname);
  fprintf(info_fp,"start time (ms) = %0.2f\n",starttime);
  fprintf(info_fp,"sample period (ms) = %0.2f\n",sample_period);
  fprintf(info_fp,"# of sources (vertices) = %d\n",num_source_vertices);
  fprintf(info_fp,"number of time points = %d\n",tpoints);
  fprintf(info_fp,"vertex list:\n");
  for (i=0;i<num_source_vertices;i++) {
    fprintf(info_fp,"   %d\n",source_vertices[i]);
  }
  fclose(info_fp);

  /* create output file */
  sprintf(tempstr,"%s/%s-%s.wt",outdir,outstem,hemi);
  MsgPrintf("%s: opening output wt file %s\n",progname,tempstr);
  out_fp = fopen(tempstr,"w");
  if (out_fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### unable to create file %s\n",progname,tempstr);
  /* write number of tpoints at start of output file */
  fwrite2(tpoints,out_fp);

  /* read each time point from stc file, write to wt file */
  MsgPrintf("%s: writing output file\n",progname);
  for (t=0;t<tpoints;t++) {
/*
    MsgPrintf("%s: writing %d values at timepoint %d\n",
              progname,num_source_vertices,t);
*/
    /* write number of non zero vertices */
    fwrite3(num_source_vertices,out_fp);
    for (i=0;i<num_source_vertices;i++) {
      vnum = source_vertices[i];
      /* write vertex number */
      fwrite3(vnum,out_fp);
      val = freadFloat(in_fp);
      if (swapflag) val = swapFloat(val);
      /* write value */
      fwriteFloat(val,out_fp);
/*
      MsgPrintf("%s: vnum=%d, val=%0.4f\n",
                progname,vnum,val);
*/
    }
  }

  free(source_vertices);
  fclose(in_fp);
  fclose(out_fp);
  
  MsgPrintf("%s: finished.\n",progname);
  exit(0);
}

