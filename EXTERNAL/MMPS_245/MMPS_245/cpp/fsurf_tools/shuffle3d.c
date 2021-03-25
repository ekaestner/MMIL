/* shufflebrik.c: read brik and scramble voxels
      created: 07/18/05 DH
     last mod: 07/18/05 DH

   purpose:
     generating noise brik with identical distribution as input brik

   input:
     afni BRIK + HEAD file

   optional input:
     mask BRIK + HEAD file

   output:
     BRIK + HEAD file
*/

// todo: either require briks with only 1 subbrik or modify output HEAD file (how?)

// todo: accept bfloat as input

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
char  maskstem[STRLEN]=UNDEFSTR;
char  indir[STRLEN]=".";
char  outstem[STRLEN]=UNDEFSTR;
char  outdir[STRLEN]=".";
int   imgindex=0;

/* other global variables */
static char *progname = NULL;
int maskflag=0;

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
  printf("    -maskstem <maskstem>       omit extension: <maskstem>.BRIK\n");
  printf("    -imgindex [0]              sub-BRIK index to shuffle\n");
  printf("    -indir    [.]              input dir\n");
  printf("    -outstem  [shuffle-instem] output file stem: <outstem>.BRIK\n");
  printf("    -outdir   [.]              output dir\n");
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
      if (MATCH(argv[i],"-maskstem") && i+1<argc) {
        strcpy(maskstem,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-imgindex")) && i+1<argc) {
        imgindex = atoi(argv[i+1]); i++;
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
  if (!MATCH(maskstem,UNDEFSTR))
    maskflag=1;
  if (MATCH(outstem,UNDEFSTR))
    sprintf(outstem,"shuffle-%s",instem);
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

/* read subbrik for a single time point */
void read_BRIK(float *values, char *indir, char *instem, int subbrik,
               int xnum, int ynum, int nslices,
               float ffact, int btype, char *byteorder)
{
  FILE *fp;
  char tempstr[STRLEN];
  long offset=0;
  int i, nvox;
  double sum=0,sum2=0,avg,stdev;
  float max=-BIGFLOAT,min=BIGFLOAT,f;
  short s;
  int c;

  if (ffact==0) ffact = 1;
  switch (btype) {
    case BRIK_BYTE:
      offset = subbrik*nslices*ynum*xnum*sizeof(char);
      break;
    case BRIK_SHORT:
      offset = subbrik*nslices*ynum*xnum*sizeof(short);
      break;
    case BRIK_FLOAT:
      offset = subbrik*nslices*ynum*xnum*sizeof(float);
      break;
    default:
      ErrorExit(ERROR_BADFILE,"%s: ### brik type %d is currently unsupported ...quitting\n",
              progname,btype);
  }  
  sprintf(tempstr,"%s/%s.BRIK",indir,instem);
  MsgPrintf("%s: reading subbrik %d of AFNI BRIK %s\n",
            progname,subbrik,tempstr);
  fp = fopen(tempstr,"rb");
  if (fp==NULL)
    ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
            progname,tempstr);
  fseek(fp,offset,SEEK_SET);  /* TODO: check file size */
  
  nvox = nslices*ynum*xnum;
  for (i=0;i<nvox;i++) {
    if (btype==BRIK_BYTE) {
      fread1(&c,fp);
      f = values[i] = ffact*(float)c;
    } else if (btype==BRIK_SHORT) {
      s = freadShort(fp);
      if(MATCH(byteorder,"LSB")) s = swapShort(s);
      f = values[i] = ffact*(float)s;
    } else if (btype==BRIK_FLOAT) {
      f = freadFloat(fp);
      if(MATCH(byteorder,"LSB")) f = swapFloat(f);
      f = values[i] = ffact*f;
    } else
      f = 0;
    if (f>max) max = f;
    if (f<min) min = f;
    sum += f;
    sum2 += f*f;
  }
  avg = sum/nvox;
  stdev = sqrt(sum2/nvox-avg*avg);
  MsgPrintf("%s: voxels=%d, avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            progname,nvox,avg,stdev,min,max);
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
  int i,j,k,n,nvox,nvox_m;
  float f;
  float *scalefacts=NULL, *scalefacts_m=NULL;
  int *briktypes=NULL, *briktypes_m=NULL;
  char byteorder[STRLEN]=UNDEFSTR, byteorder_m[STRLEN]=UNDEFSTR;
  float *dat=NULL, *maskdat=NULL;
  int *maskvox=NULL;
  int xnum,ynum,nslices,tpoints;
  int xnum_m,ynum_m,nslices_m,tpoints_m;

  randseed();  // different random series each time program is run

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

  /* get dimensions of mask, require match */
  if (maskflag) {
    sprintf(infile,"%s/%s.HEAD",indir,maskstem);
    MsgPrintf("%s: reading AFNI header %s\n",progname,infile);
    read_HEAD(infile,&xnum_m,&ynum_m,&nslices_m,&tpoints_m,
              &scalefacts_m,&briktypes_m,byteorder_m);
    if (tpoints_m>1)
      MsgPrintf("%s: mask BRIK has %d tpoints... using 0th\n",
                progname,tpoints_m);
    if(xnum_m!=xnum || ynum_m!=ynum || nslices_m!=nslices)
      ErrorExit(ERROR_BADFILE,
        "%s: ### image dimensions for mask do not match dataset...quitting\n",
        progname);
  }

  /* allocate memory for data */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  dat = (float *)calloc(nvox,sizeof(float));           MTEST(dat);
  if (maskflag) {
    maskdat = (float *)calloc(nvox,sizeof(float));     MTEST(maskdat);
  }
  maskvox = (int *)calloc(nvox,sizeof(int));           MTEST(maskvox);
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* read brik(s) */
  MsgPrintf("%s: reading input data set\n",progname);
  read_BRIK(dat,indir,instem,imgindex,
            xnum,ynum,nslices,
            scalefacts[imgindex],briktypes[imgindex],byteorder);
  if (maskflag) {
    read_BRIK(maskdat,indir,maskstem,0,
              xnum,ynum,nslices,
              scalefacts_m[0],briktypes_m[0],byteorder_m);
    for(i=0,nvox_m=0;i<nvox;i++) {
      if(maskdat[i]!=0) {
        maskvox[nvox_m]=i;
        nvox_m++;
      }
    }
  } else {
    for(i=0;i<nvox;i++) maskvox[i]=i;
    nvox_m = nvox;
  }

  /* shuffle voxels inside the mask */
  for(k=0;k<nvox_m;k++) {
    n=k+random()%(nvox_m-k);
    i = maskvox[k];
    j = maskvox[n];
    f = dat[j];
    dat[j]=dat[i];
    dat[i]=f;
  }

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

