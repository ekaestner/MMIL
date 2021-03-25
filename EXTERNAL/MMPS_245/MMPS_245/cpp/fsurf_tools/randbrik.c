/* randbrik.c: generate a brik with random values
      created: 08/05/05 DH
     last mod: 08/05/05 DH

   purpose:
     generating noise brik with dimensions of input brik

   input:
     afni BRIK + HEAD file

   optional input:
     mask BRIK + HEAD file

   output:
     BRIK + HEAD file

   acknowlegdgements:
     portions of code adapted from AFNI's AlphaSim
*/

#include "surflib.h"

#define MINARGC 3

#define BRIK_BYTE    0
#define BRIK_SHORT   1
#define BRIK_LONG    2
#define BRIK_FLOAT   3
#define BRIK_DOUBLE  4
#define BRIK_COMPLEX 5


/* parameter defaults */
char  instem[STRLEN]=UNDEFSTR;
char  indir[STRLEN]=".";
char  outstem[STRLEN]=UNDEFSTR;
char  outdir[STRLEN]=".";
int   imgindex=0;
float stdev = 1.0;

/* other global variables */
static char *progname = NULL;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem   <instem>         omit extension: <instem>.BRIK\n");
  printf("\n");
  printf("  Optional parameters:  (defaults in []'s)\n");
  printf("    -indir    [.]              input dir\n");
  printf("    -outstem  [randvals+orig]  output file stem: <outstem>.BRIK\n");
  printf("    -outdir   [.]              output dir\n");
  printf("    -stdev    [1.0]            standard deviation of random values\n");
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
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-stdev") && i+1<argc) {
        stdev = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
        ErrorExit(ERROR_BADPARM,"%s: ### error parsing option: %s\n",progname,argv[i]);
    } else
      ErrorExit(ERROR_BADPARM,"%s: ### error parsing option: %s\n",progname,argv[i]);
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments before proceeding */
  if (MATCH(instem,UNDEFSTR))
    ErrorExit(ERROR_BADPARM,"%s: ### instem not specified ...quitting\n",
              progname);
  if (MATCH(outstem,UNDEFSTR))
    sprintf(outstem,"randvals+orig");
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
  sprintf(tempstr,"%s/%s.BRIK",indir,instem);
  if (!FileExists(tempstr))
    ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s.HEAD",indir,instem);
  if (!FileExists(tempstr))
    ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
              progname,tempstr);
  if (imgindex < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### bad imgindex: %d => must be >= 0 ...quitting\n",
              progname,imgindex);
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

float uniform ()
{
  return ( (float)drand48() );
}

void normal (float * n1, float * n2)
{
  float u1, u2;
  float r;

  u1 = 0.0;
  while (u1 <= 0.0) {u1 = uniform();}
  u2 = uniform();

#if 0
  r = sqrt(-2.0*log(u1));
#else
  r = stdev * sqrt(-2.0*log(u1));
#endif
  *n1 = r * cos(2.0*M_PI*u2);
  *n2 = r * sin(2.0*M_PI*u2);
}

void generate_randvals(float *vals, int nvals)
{
  int k, nvalsdiv2;
  float n1, n2;
  
  nvalsdiv2 = nvals/2;
  
  MsgPrintf("%s: generating random values\n",progname);
  for (k=0;k<nvalsdiv2;k++) {
    normal(&n1,&n2);
    vals[k] = n1;
    vals[k+nvalsdiv2] = n2;
  }
  normal(&n1,&n2);
  vals[nvals-1] = n1;
}

void read_HEAD(char *fname, 
               int *xsize, int *ysize, int *zsize, int *tsize,
               float **ffacts, int **btypes, char *byteorder)
{
  FILE *fp;
  char line[STRLEN];
  char entry[STRLEN];
  char value[STRLEN];
  int i;

  float ffact;
  int btype;

  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,fname);
  while (fgets(line,99,fp) != NULL) {
    sscanf(line,"%s %*s %s",entry,value);
    if (MATCH(value,"DATASET_DIMENSIONS")) {
      fgets(line,99,fp);
      sscanf(line,"%s %*s %s",entry,value);
      fgets(line,99,fp);
      sscanf(line,"%d %d %d",xsize,ysize,zsize);
    }
    if (MATCH(value,"BRICK_FLOAT_FACS")) {
      fgets(line,99,fp);
      sscanf(line,"%s %*s %s",entry,value);
      *tsize = atoi(value);
      if (*ffacts) {
        MsgPrintf("%s: freeing old float facs...\n",progname);
        free(*ffacts);
      }
      *ffacts = (float *)calloc(*tsize,sizeof(float)); MTEST(*ffacts);
      MsgPrintf("%s: %d float_facs in %s:  ",progname,*tsize,fname);
      for(i=0;i<*tsize;i++) {
        fscanf(fp,"%f",&ffact);
        MsgPrintf("%1.4f ",ffact);
        (*ffacts)[i] = ffact;
      }
      MsgPrintf("\n");
    }
    if (MATCH(value,"BRICK_TYPES")) {
      fgets(line,99,fp);
      sscanf(line,"%s %*s %s",entry,value);
      *tsize = atoi(value);
      if (*btypes) {
        MsgPrintf("%s: freeing old brik types...\n",progname);
        free(*btypes);
      }
      *btypes = (int *)calloc(*tsize,sizeof(int)); MTEST(*btypes);

      MsgPrintf("%s: %d brick_types in %s:  ",progname,*tsize,fname);
      for(i=0;i<*tsize;i++) {
        fscanf(fp,"%d",&btype);
        MsgPrintf("%d ",btype);
        (*btypes)[i] = btype;
      }
      MsgPrintf("\n");
    }
    if (MATCH(value,"BYTEORDER_STRING")) {
      fgets(line,99,fp);
      sscanf(line,"%s %*s %s",entry,value);
      fgets(line,99,fp);
      sscanf(line,"%s",byteorder);
      if(MATCH(byteorder,"'LSB_FIRST~"))
        strcpy(byteorder,"LSB");
      else
        strcpy(byteorder,"MSB");
      MsgPrintf("%s: byteorder of %s: %s\n",progname,fname,byteorder);
    }
  }
  fclose(fp);
}

/* write subbrik */
void write_BRIK(float *values, char *outdir, char *outstem,
               int xnum, int ynum, int nslices,
               float ffact, int btype, char *byteorder)
{
  FILE *fp;
  char tempstr[STRLEN];
  int i, nvox;
  double sum=0,sum2=0,avg,sd;
  float max=-BIGFLOAT,min=BIGFLOAT,f;
  short s;
  int c;

  if (ffact==0) ffact = 1;
  sprintf(tempstr,"%s/%s.BRIK",outdir,outstem);
  MsgPrintf("%s: writing AFNI BRIK %s\n",
            progname,tempstr);
  fp = fopen(tempstr,"wb");
  if (fp==NULL)
    ErrorExit(ERROR_BADPARM,"%s: ### unable to create file %s ...quitting\n",
            progname,tempstr);

  if (btype==BRIK_BYTE)
    MsgPrintf("%s: BRIK data type is BYTE\n",progname);
  else if (btype==BRIK_SHORT)
    MsgPrintf("%s: BRIK data type is SHORT\n",progname);
  else if (btype==BRIK_FLOAT)
    MsgPrintf("%s: BRIK data type is FLOAT\n",progname);
  else
    MsgPrintf("%s: BRIK data type is unsupported\n",progname);

  nvox = nslices*ynum*xnum;
  for (i=0;i<nvox;i++) {
    f = values[i]/ffact;
    if (f>max) max = f;
    if (f<min) min = f;
    sum += f;
    sum2 += f*f;
    if (btype==BRIK_BYTE) {
      c = (int)(f);
      fwrite1(c,fp);
    } else if (btype==BRIK_SHORT) {
      s = (short)(f);
      if(MATCH(byteorder,"LSB")) s = swapShort(s);
      fwrite2(s,fp);
    } else if (btype==BRIK_FLOAT) {
      if(MATCH(byteorder,"LSB")) f = swapFloat(f);
      fwrite4(f,fp);
    } else
      f = 0;
  }
  avg = sum/nvox;
  sd = sqrt(sum2/nvox-avg*avg);
  MsgPrintf("%s: voxels=%d, avg=%6.2lf, stdev=%6.2lf, min=%6.2f, max=%6.2f\n",
            progname,nvox,avg,sd,min,max);
}

int main(int argc, char **argv)
{
  char infile[STRLEN], outfile[STRLEN], tempstr[STRLEN];
  int nvox;
  float *scalefacts=NULL;
  int *briktypes=NULL;
  char byteorder[STRLEN]=UNDEFSTR;
  float *dat=NULL;
  int xnum,ynum,nslices,tpoints;

  randseed();

  parse_args(argc,argv);

  /* get dimensions of dataset */
  sprintf(infile,"%s/%s.HEAD",indir,instem);
  MsgPrintf("%s: reading AFNI header %s\n",progname,infile);
  read_HEAD(infile,&xnum,&ynum,&nslices,&tpoints,
            &scalefacts,&briktypes,byteorder);
  if (imgindex>=tpoints)
    ErrorExit(ERROR_BADPARM,"%s: ### image offset (%d) >= # time points (%d) ...quitting\n",
              progname,imgindex,tpoints);
  nvox = xnum*ynum*nslices;

  /* allocate memory for data */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  dat = (float *)calloc(nvox,sizeof(float));           MTEST(dat);
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* generate random values for each voxel */
  generate_randvals(dat,nvox);

  /* write brik */
  sprintf(infile,"%s/%s.HEAD",indir,instem);
  sprintf(outfile,"%s/%s.HEAD",outdir,outstem);
  sprintf(tempstr,"cp %s %s", infile,outfile);
  system(tempstr);
  write_BRIK(dat,outdir,outstem,
             xnum,ynum,nslices,
             scalefacts[imgindex],briktypes[imgindex],byteorder);

  exit(0);
}


