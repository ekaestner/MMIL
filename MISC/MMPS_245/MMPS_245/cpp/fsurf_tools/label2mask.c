/* label2mask.c: convert label file to value (.w) file
      created: 09/29/04 DH
     last mod: 09/29/04 DH

   purpose:
     creating a mask file from thresholding surface stats using a mask file

   input:
     label file (ascii .label file)

   output:
     masked file (.w file)
*/

#include "surflib.h"

#define MINARGC 2

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char labelstem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
int overwrite = 0;
float maskval = 1;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s labelstem [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    labelstem                        omit extension, hemi: <labelstem>-rh.label\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -hemi    <str>  [rh]             hemisphere (rh or lh)\n");
  printf("    -indir   <str>  [.]              input dir\n");
  printf("    -outdir  <str>  [.]              output dir (must exist)\n");
  printf("    -outstem <str>  [labelstem-mask] output file stem (for w file)\n");
  printf("    -maskval  <f>   [1]              value to assign to vertices in label\n");
  printf("    -overwrite                       force overwrite of existing output files\n");
  printf("    -quiet                           suppress messages\n");
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
  strcpy(labelstem,argv[1]);
  for (i=2;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-overwrite")){
        overwrite = 1;
      }
      else if (MATCH(argv[i],"-maskval") && i+1<argc){
        maskval = atof(argv[i+1]); i++;
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
  if (MATCH(labelstem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### labelstem not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) sprintf(outstem,"%s-mask",labelstem);
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
  char labelname[STRLEN];
  char outname[STRLEN];
  MRIS *mris;
  FILE *fp;
  int vnum;
  float x, y, z, rip;
  char line[STRLEN];
  char subj[STRLEN];
  int i;
  unsigned char c;

  parse_args(argc,argv);

  /* set file names */
  sprintf(labelname,"%s/%s-%s.label",indir,labelstem,hemi);
  sprintf(outname,"%s/%s-%s",outdir,outstem,hemi);

  /* first open label file as binary and get subject name*/
  fp = fopen(labelname, "rb") ;
  if (!fp)
    ErrorExit(ERROR_BADFILE,"%s: ### File %s not found\n",progname,labelname);
  for(i=0;i<STRLEN;i++) {fread(&c,1,1,fp); line[i]=c;}    
  fclose(fp);
  sscanf(line, "%*s %*s %*s %*s %*s %s", subj);
  printf("%s: subject = %s\n",progname,subj);
  
  /* load surface */
  mris = openSurface(subj,hemi,"orig");

  /* clear values from surface */
  for(vnum=0;vnum<mris->nvertices;vnum++) mris->vertices[vnum].val = 0;

  /* open label file again as text */
  fp = fopen(labelname, "r") ;
  if (!fp)
    ErrorExit(ERROR_BADFILE,"%s: ### File %s not found\n",progname,labelname);

  /* skip first line (contains number of vertices in label) */
  fgetl(line, STRLEN, fp);

  /* now read vertex numbers from label file */
  while(!feof(fp)){
    fgetl(line, STRLEN, fp);
    if(feof(fp)) break;
    sscanf(line, "%d %f %f %f %f\n", &vnum, &x, &y, &z, &rip);
    if(vnum>=0 && vnum<mris->nvertices) {
      mris->vertices[vnum].val = maskval;
    } else
      ErrorExit(ERROR_BADFILE,"%s: ### bad vertex number %d\n",progname,vnum);
  }
  fclose(fp);

  /* write out file */
  MRISwriteValues(mris,outname);

  exit(0);
}

