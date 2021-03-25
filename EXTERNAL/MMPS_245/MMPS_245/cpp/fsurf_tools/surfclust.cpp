/* surfclust.cpp: read surface stats file and gets stats on clusters
      created: 12/22/04 DH
     last mod: 08/18/04 DH

   purpose:
     finding the size of supra-threshold clusters for cluster-size exclusion

   input:
     stats file (w file)

   output:
     masked stats file (w file)

   acknowledgements:
     Doug Greve's mri_surfcluster
*/

#include "surflib.h"
#include "clustlib.h"

#define MINARGC 3
#define SMALLFLOAT 1e-10
#define DEFOUTSTEM "surfclust-output"

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
int threshabsflag = 0;
float thresh = SMALLFLOAT;
int minarea = 1;
char surf[STRLEN]="smoothwm";
int maxvmaskflag = 0;
int summaryflag = 0;
int talflag = 0;
int fixmniflag = 1;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  instem            omit extension, hemi: <instem>-rh.w\n");
  printf("    -subj    <str>             subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem outstem           omit extension, hemi: <outstem>-rh.w\n");
  printf("    -indir   <str>  [.]        input dir\n");
  printf("    -outdir  <str>  [.]        output dir\n");
  printf("    -hemi    <str>  [rh]       hemisphere (rh or lh)\n");
  printf("    -thresh   <f>   [10^-10]   minimimum value for clusters\n");
  printf("    -threshabs                 when thresholding, use absolute values\n");
  printf("    -minarea  <f>   [1]        minimum cluster size (mm^2)\n");
  printf("    -surf     <str> [smoothwm] surface used for area calculations\n");
  printf("    -maxvmask                  output w file with only maximum vertices\n");
  printf("                                 for each cluster set to 1\n");
  printf("    -summary                   output text summary of each cluster\n");
  printf("    -talcoords                 transform coords to talairach space\n");
  printf("    -nofixmni                  do not fix mni coordinates (Matthew Brett's txfm)\n");
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
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
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
      if (MATCH(argv[i],"-thresh") && i+1<argc){
        thresh = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-threshabs")){
        threshabsflag = 1;
      } else
      if (MATCH(argv[i],"-minarea") && i+1<argc){
        minarea = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-maxvmask")){
        maxvmaskflag = 1;
      } else
      if (MATCH(argv[i],"-summary")){
        summaryflag = 1;
      } else
      if (MATCH(argv[i],"-talcoords")){
        talflag = 1;
      } else
      if (MATCH(argv[i],"-nofixmni")){
        fixmniflag = 0;
      } else
      if (MATCH(argv[i],"-fixmni")){
        fixmniflag = 1;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    } else
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(instem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### instem not specified ...quitting\n",
              progname);
  }
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    strcpy(outstem,DEFOUTSTEM);
  }
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
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

int main(int argc, char **argv)
{
  int i,k,m,ci;
  int ecode;
  char valfile[STRLEN];
  MRIS *mris;
  ClusterList clusters;
  Cluster clust1, clust2;

  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,surf);

  /* load talairach transform */
  if (talflag) {
    if(mris->transform_loaded) {
      if(!getQuiet()) printTransform(mris->linear_transform,"linear_transform");
    } else {
      MsgPrintf("%s: transform NOT loaded succesfully -- cannot do talairach transform\n", progname);
      talflag = 0;
    }
  }
  
  /* read values */
  sprintf(valfile,"%s/%s-%s.w",indir,instem,hemi);
  ecode = MRISreadValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,valfile);
  }
  MsgPrintf("%s: finding clusters\n",progname);
  clusters.init(mris->nvertices);
  clusters.findClusters(mris,thresh,threshabsflag);
  clusters.calcStats(mris);
  if(talflag) {
    clusters.transformCoords(mris);
    if(fixmniflag) clusters.fixMNICoords();
  }
  clusters.sortByArea(mris);

  k=m=0;
  for (i=0;i<clusters.size();i++) {
    if (clusters[i].size())
      k++;
    if (clusters[i].Area()>=minarea) m++;
  }
  MsgPrintf("%s: number of clusters = %d\n",progname,k);
  MsgPrintf("%s: number of clusters as large as %d mm^2 = %d\n",progname,minarea,m);

  /* applying cluster exclusion */
  if (maxvmaskflag) {
    MsgPrintf("%s: setting maxverts to 1, all others to 0\n",progname);
    for (k=0;k<mris->nvertices;k++) mris->vertices[k].val = 0;
    for (i=0;i<clusters.size();i++) {
      k = clusters[i].Max_v();
      mris->vertices[k].val = 1;
    }
  } else {
    MsgPrintf("%s: setting excluded vals to zero\n",progname);
    for (k=0;k<mris->nvertices;k++) {
      ci = clusters.getClusterIndex(k);
      if (ci<0)
        mris->vertices[k].val = 0;
      else if (clusters[ci].Area() < minarea)
        mris->vertices[k].val = 0;
    }
  
  }

  /* write to file */
  MsgPrintf("%s: writing output\n",progname);
  sprintf(valfile,"%s/%s-%s.w",outdir,outstem,hemi);
  ecode = MRISwriteValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error writing value file %s\n",
              progname,valfile);
  }

  /* output summary */
  if(summaryflag) {
    printf("\nsurfclust output summary:\n");
    printf("ClustNo  #Vertices  Area(mm^2)  MinVal   MaxVal   AvgVal   StDev   CentX   CentY   CentZ    MaxX    MaxY    MaxZ\n");
    for (i=0;i<clusters.size();i++) {
      printf(" %3d %10d %10.2f %9.2f %8.2f %9.2f %8.2f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n",
        clusters[i].Index(),clusters[i].size(),clusters[i].Area(),
        clusters[i].Minval(),clusters[i].Maxval(),
        clusters[i].Avgval(),clusters[i].Stdev(),
        clusters[i].Cent_x(),clusters[i].Cent_y(),clusters[i].Cent_z(),
        clusters[i].Max_x(),clusters[i].Max_y(),clusters[i].Max_z());
    }
  }  
}

