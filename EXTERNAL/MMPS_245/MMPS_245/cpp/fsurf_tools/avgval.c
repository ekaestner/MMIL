/* avgvals.c: open a value file (.w) and average all values
      created: 09/29/04 DH
     last mod: 09/29/04 DH

   purpose:
     

   input:
     w value files

   output:
     stdout
*/
#include "surflib.h"

#define MINARGC 2
#define MATCH(A,B)   (!strcmp(A,B))
#define SQR(x) ((x)*(x))

static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char hemi[STRLEN]="rh";
char subj[STRLEN]=UNDEFSTR;
char ext[STRLEN]="bfloat";
int nslices=0;
int xnum=0,ynum=0,depth=0;
int zeroflag=0;

void usage()
{
  printf("\n");
  printf("Usage: %s instem {-name subjname} [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    instem   omit infixes, suffix: <instem>_000.bfloat\n");
  printf("                                or <instem>-{rh,lh}.w\n");
  printf("    -name   <str>  subject name (required only if ext = w)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -indir    <str>  [.]       input dir\n");
  printf("    -nslices  <int>  [0]       slice count (bfloat ext only)\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -ext      <str>  [bfloat]  input extension (bfloat or w)\n");
  printf("    -hemi     <str>  [rh]      hemisphere (rh or lh) (w or 1D-bfloat only)\n");
  printf("    -countzeros                include zero values in average\n");
  printf("    -debug                     print debug statements\n");
  printf("    -h                         print this usage statement\n");
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i;
  char *tempfix=NULL;
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

  tempfix    = (char *)malloc(STRLEN*sizeof(char));      MTEST(tempfix);

  setQuiet(1);

  /* parse arguments */
  strcpy(instem,argv[1]);  
  for (i=2;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-name") && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-nslices") || MATCH(argv[i],"-n")) && i+1<argc) {
        nslices = atoi(argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-ext") && i+1<argc){
        strcpy(ext,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      }
      else if (MATCH(argv[i],"-debug")){
        setQuiet(0);
      }
      else if (MATCH(argv[i],"-countzero")){
        zeroflag=1;
      }
      else if (MATCH(argv[i],"-h")){
        usage();
        exit(0);
      }
      else {printf("-1\n");exit(1);}
    } else {printf("-1\n");exit(1);}
  }

  /* check arguments */
  if (MATCH(instem,UNDEFSTR)) {
    if(getQuiet()) printf("-1\n");
    MsgPrintf("%s: ### must supply an instem ...quitting\n",
              progname);
    exit(1);
  }

  if (!MATCH(ext,"bfloat") && !MATCH(ext,"w")) {
    if(getQuiet()) printf("-1\n");
    MsgPrintf("%s: ### extension %s not supported ...quitting\n",
              progname,ext);
    exit(1);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    if(getQuiet()) printf("-1\n");
    MsgPrintf("%s: ### %s: hemi must be rh or lh\n",
              progname,hemi);
    exit(1);
  }
  if (MATCH(ext,"w")) {
    nslices = 1;
  } else if (MATCH(ext, "bfloat")) {
    i = getNumBfloats(indir,instem,NULL);
    if(i<=0) {
      if(getQuiet()) printf("-1\n");
      MsgPrintf("%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
                progname,indir,instem);
      exit(1);
    }
    if(nslices==0 || nslices>i) {
      nslices = i;
      if(tempfix==NULL) {
        MsgPrintf("%s: %d slice(s) found in %s/ with stem %s\n",
                  progname,i,indir,instem);
      } else {
        MsgPrintf("%s: %d slice(s) found in %s/ with stem %s%s\n",
                  progname,i,indir,instem,tempfix);
      }
    }
  }
  if (MATCH(subj,UNDEFSTR) && (MATCH(ext,"w"))) {
    if(getQuiet()) printf("-1\n");
    MsgPrintf("%s: ### subject name not specified ...quitting\n",
              progname);
    exit(1);
  }
}

void read_stats(float ***stats, char *indir, char *instem)
{
  int i,n;
  char fname[STRLEN];
  float *vals;

  if (MATCH(ext,"bfloat")){
    for (n=0;n<nslices;n++) {
      sprintf(fname,"%s/%s_%03d.bfloat",indir,instem,n);
      if(readFloatImage(fname,stats[n],xnum,ynum)==-1) {
        if(getQuiet()) printf("-1\n");
        MsgPrintf("%s: ### error reading float image...quitting\n",
                  progname);
        exit(1);
      }
    }
  } else if (MATCH(ext, "w")) {
    sprintf(fname,"%s/%s-%s.w",indir,instem,hemi);
    vals = (float *)calloc(ynum,sizeof(float));             MTEST(vals);
    if(readSurfVals(fname,vals,ynum)==-1) {
      if(getQuiet()) printf("-1\n");
      MsgPrintf("%s: ### error reading surface values...quitting\n",
                progname);
      exit(1);
    }
    for(i=0;i<ynum;i++) stats[0][i][0]=vals[i];
  }    
}

int main(int argc, char **argv)
{
  int n,i,j,k;
  float ***dat=NULL;
  char fname[STRLEN];
  float avg;

  parse_args(argc,argv);

  /* get dimensions of input value file */
  if (MATCH(ext,"bfloat")) {
    sprintf(fname,"%s/%s_000.hdr",indir,instem);
    if(readFloatHeader(fname,&xnum,&ynum,&depth)==-1) {
      if(getQuiet()) printf("-1\n");
      MsgPrintf("%s: ### error reading header file...quitting\n",
                progname);
      exit(1);
    }
  } else if (MATCH(ext,"w")) {
    ynum = nverticesSurf(subj,hemi);
    if(ynum==-1) {
      if(getQuiet()) printf("-1\n");
      MsgPrintf("%s: ### error reading number of vertices...quitting\n",
                progname);
      exit(1);
    }
    xnum = 1;
  }
  
  /* allocate memory for data and calculations */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  dat = (float ***)calloc(nslices,sizeof(float **));        MTEST(dat);
  for (i=0;i<nslices;i++) {
    dat[i] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat);
    for (j=0;j<ynum;j++) {
      dat[i][j] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat);
    }
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* read data */
  read_stats(dat,indir,instem);

  MsgPrintf("%s: finished reading stats\n",progname);

  /* calculate average */
  n=0;
  avg=0;
  for (i=0;i<nslices;i++)
    for (j=0;j<ynum;j++)
      for (k=0;k<xnum;k++) {
        if(zeroflag || dat[i][j][k]!=0) {
          avg += dat[i][j][k];
          n++;
        }
      }

  MsgPrintf("%s: number of values = %d\n",progname,n);
  MsgPrintf("%s: sum of values = %f\n",progname,avg);
  MsgPrintf("%s: average value = ",progname);

  if (n==0) avg = 0;
  else      avg /= n;

  printf("%f\n",avg);
  exit(0);
}

