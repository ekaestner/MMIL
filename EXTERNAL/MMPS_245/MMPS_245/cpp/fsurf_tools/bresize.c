/* bresize.c: calculates average, t-stats, and p-values from bfloat data sets
      created: 09/03/03 DH
     last mod: 09/22/03 DH

   purpose:
     resizing bfloats, especially for converting 1D spherical surface bfloats
     into 3D bfloats which can be converted into afni briks

   input:
     bfloat volumes or single "slice" output from mri_vol2surf (not timeseries)

   output:
     bfloat volumes or single slices
*/

#include "surflib.h"

#define MINARGC 2
#define DEFNEWX 64
#define DEFNEWY 64

/* global variables */
char progname[20]="bresize";

/* parameter defaults */
int nslices=0;
int xnum=0,ynum=0,depth=0;
int newxnum=DEFNEWX,newynum=DEFNEWY;
int newnslices=0;
int outvox=0;
int truncvox=0;
char instem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char outstem[STRLEN]="output";
char subjname[STRLEN]=UNDEFSTR;
char hemi[STRLEN]=UNDEFSTR;
int unfold = 0;

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);
void read_stats(float ***stats, char *indir, char *instem);
void write_stats(float ***stats, char *outdir, char *outstem);

int main(int argc, char **argv)
{
  int i,j,n,k,nvoxels,newnvoxels;
  float f;
  float ***instats;
  float ***outstats;
  float *unfolded = NULL;
  char fname[STRLEN];

  sprintf(progname,"bresize");

  if (argc<MINARGC) {usage(); exit(0);}
  parse_args(argc,argv);

  sprintf(fname,"%s/%s_000.hdr",indir,instem);
  if(readFloatHeader(fname,&xnum,&ynum,&depth)==-1) {
    ErrorExit(ERROR_BADFILE,"%s: ### error reading header %s...quitting\n",
           progname,fname);
  }

  nvoxels = nslices*ynum*xnum;
  MsgPrintf("%s: total number of voxels in input volume: %d\n",progname,nvoxels);

  if(unfold) {
    newnslices=1;
    newxnum=1;
    if(outvox > 0) {
      newnvoxels=outvox;
    } else {
      newnvoxels=nvoxels-truncvox;
    }
    newynum=newnvoxels;
  } else {
    f = (float)nvoxels/(float)(newxnum*newynum);
    newnslices = (f>(int)f) ? (int)f + 1 : (int)f;
    newnvoxels = newnslices*newxnum*newynum;
  }

  MsgPrintf("%s: number of slices in output volume: %d\n",progname,newnslices);
  MsgPrintf("%s: total number of voxels in output volume: %d\n",progname,newnvoxels);
  if(newnvoxels>nvoxels) MsgPrintf("         (extras will be set to zero)\n");

  /* allocate memory */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  instats = (float ***)calloc(nslices,sizeof(float **));        MTEST(instats);
  for (n=0;n<nslices;n++) {
    instats[n] = (float **)calloc(ynum,sizeof(float *));        MTEST(*instats);
    for (i=0;i<ynum;i++) {
      instats[n][i] = (float *)calloc(xnum,sizeof(float));      MTEST(**instats);
    }
  }
  outstats = (float ***)calloc(newnslices,sizeof(float **));    MTEST(outstats);
  for (n=0;n<newnslices;n++) {
    outstats[n] = (float **)calloc(newynum,sizeof(float *));    MTEST(*outstats);
    for (i=0;i<newynum;i++) {
      outstats[n][i] = (float *)calloc(newxnum,sizeof(float));  MTEST(**outstats);
    }
  }
  if(newnvoxels>nvoxels) {
    unfolded = (float *)calloc(newnvoxels,sizeof(float));       MTEST(unfolded);
  } else {
    unfolded = (float *)calloc(nvoxels,sizeof(float));          MTEST(unfolded);
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* read data images */
  MsgPrintf("%s: reading stats from file...\n",progname);
  read_stats(instats,indir,instem);

  /* unfold into linear array */
  MsgPrintf("%s: unfolding into linear array...\n",progname);
  k=0;
  for (n=0;n<nslices;n++)
    for (i=0;i<ynum;i++)
      for (j=0;j<xnum;j++)
        unfolded[k++]=instats[n][i][j];

  /* fill in remainder of array with zeros */
  if(newnvoxels>nvoxels) {
    MsgPrintf("%s: filling in remainder of array with zeros...\n",progname);
    for (k=nvoxels;k<newnvoxels;k++) unfolded[k]=0;
  }

  /* refold into output array */
  MsgPrintf("%s: refolding into output array...\n",progname);
  k=0;
  for (n=0;n<newnslices;n++)
    for (i=0;i<newynum;i++)
      for (j=0;j<newxnum;j++)
        outstats[n][i][j] = unfolded[k++];

  /* write output to file */
  MsgPrintf("%s: writing output to file...\n",progname);
  write_stats(outstats,outdir,outstem);

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
  printf("    -nslices   <int> [0]        input slice count\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -indir     <str> [.]        input dir\n");
  printf("    -outdir    <str> [.]        output dir (must exist)\n");
  printf("    -outstem   <str> [output]   output file stem\n");
  printf("    -newx      <int> [%d]       new number of elements in x-dimension\n",
                                newxnum);
  printf("    -newy      <int> [%d]       new number of elements in y-dimension\n",
                                newynum);
  printf("    -unfold                     save as one-dimensional array\n");
  printf("    -surfverts <str> [unknown]  subject name\n");
  printf("       output as many voxels as vertices in subject\'s orig surface (unfold only)\n");
  printf("    -hemi      <str> [unknown]  hemisphere (rh or lh)\n");
  printf("       required when surfverts option is used\n");
  printf("    -outvox    <int> [-invox-]  output file with this many voxels (unfold only)\n");
  printf("    -trunc     <int> [0]        truncate this many voxels from end (unfold only)\n");
  printf("       note: surfverts, outvox, and trunc are mutually exclusive\n");
  printf("             order of precedence is same as order above\n");
  printf("    -quiet                      suppress messages\n");
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
      else if (MATCH(argv[i],"-newx") && i+1<argc) {
        newxnum = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-newy") && i+1<argc) {
        newynum = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-trunc") && i+1<argc) {
        truncvox = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outvox") && i+1<argc) {
        outvox = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-surfverts") && i+1<argc) {
        strcpy(subjname,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-outstem") && i+1<argc) {
        strcpy(outstem,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-unfold")){
        unfold = 1;
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
  if(MATCH(instem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### must supply an instem ...quitting\n",progname);
  }
  fp = fopen(indir,"r");
  if(fp==NULL) {
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",progname,indir);
  } else fclose(fp);
  if(!isadir(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",progname,indir);
  }
  fp = fopen(outdir,"r");
  if(fp==NULL) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",progname,outdir);
  } else fclose(fp);
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",progname,outdir);
  }
  if(newxnum<0 || newynum<0) {
    ErrorExit(ERROR_BADPARM,"%s: ### new x (%d) and y (%d) must be greater than 0 ...quitting\n",
           progname,newxnum,newynum);
  }
  if(!unfold && (outvox || truncvox || !MATCH(subjname,UNDEFSTR))) {
    MsgPrintf("%s: ## no voxels will be truncated (only when -unfold flag is set)\n",
             progname);
  } else if(!MATCH(subjname,UNDEFSTR)) {
    if(MATCH(hemi,UNDEFSTR)) {
      ErrorExit(ERROR_BADPARM,"%s: ### must supply hemi if using surfverts options ...quitting\n",
             progname);
    }
    if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
      ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
    }
    if(outvox>0) {
      MsgPrintf("%s: ## ignoring outvox because surfverts takes precedence\n",
               progname);
    }
    if(truncvox>0) {
      truncvox = 0;
      MsgPrintf("%s: ## ignoring truncvox because surfverts takes precedence\n",
                progname);
    }
    outvox = nverticesSurf(subjname, hemi);
    if (outvox == -1) {    
      ErrorExit(ERROR_NOFILE,"%s: ### error reading number of %s surface vertices for %s ...quitting\n",
             progname,hemi,subjname);
    }
  } else {
    if(outvox && truncvox) {
      truncvox = 0;
      MsgPrintf("%s: ## ignoring trunc because outvox takes precedence\n",
               progname);
    }
    if(truncvox<0) {
      ErrorExit(ERROR_BADPARM,"%s: ### number of voxels to truncate (%d) must be greater than 0 ...quitting\n",
             progname,truncvox);
    }
    if(outvox<0) {
      ErrorExit(ERROR_BADPARM,"%s: ### number of voxels to output (%d) must be greater than 0 ...quitting\n",
             progname,outvox);
    }
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

void read_stats(float ***stats, char *indir, char *instem)
{
  int n;
  char fname[STRLEN];

  for (n=0;n<nslices;n++) {
    sprintf(fname,"%s/%s_%03d.bfloat",indir,instem,n);
    if(readFloatImage(fname,stats[n],xnum,ynum)==-1) {
      ErrorExit(ERROR_NOFILE,"%s: ### error reading float image...quitting\n",
           progname);
    }
  }
}

void write_stats(float ***stats, char *outdir, char *outstem)
{
  int n;
  char fname[STRLEN];

  for (n=0;n<newnslices;n++) {
    sprintf(fname,"%s/%s_%03d.bfloat",outdir,outstem,n);
    if(writeFloatImage(fname,stats[n],newxnum,newynum)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error writing float image...quitting\n",
           progname);
    }
    sprintf(fname,"%s/%s_%03d.hdr",outdir,outstem,n);
    if(writeFloatHeader(fname,newxnum,newynum,1,0)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error writing float header...quitting\n",
           progname);
    }
  }
}

