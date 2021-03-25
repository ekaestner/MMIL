/* avgtseries.c: read 3d+t BRIK or wt file, extract average time series for surface roi
     created:  10/01/04 DH
     last mod: 05/05/04 DH

   purpose:
     averaging 3d or surface timeseries for surface roi

   input:
     afni BRIK + HEAD files or wt file
     register.dat file
     surface mask file

   output:
     ascii file containing average timeseries


   todo:
     allow register.dat to be full path name

*/

#include "surflib.h"

#define MINARGC 5

#define BRIK_BYTE    0
#define BRIK_SHORT   1
#define BRIK_LONG    2
#define BRIK_FLOAT   3
#define BRIK_DOUBLE  4
#define BRIK_COMPLEX 5

/* parameter defaults */
char datatype[STRLEN]="vol";
int voldataflag = 1;
char  instem[STRLEN]=UNDEFSTR;
char  maskstem[STRLEN]=UNDEFSTR;
char  subj[STRLEN]=UNDEFSTR;
char  indir[STRLEN]=".";
char  maskdir[STRLEN]=".";
char  outdir[STRLEN]=".";
char  outfile[STRLEN]=UNDEFSTR;
char  hemi[STRLEN]="rh";
char  regdat[STRLEN]="register.dat";
float dpaint=0;
int   sparseavgflag=0;

/* other global variables */
static char *progname = NULL;
MRIS *mris;
float ***dat=NULL;
float *maskvals=NULL;
int xnum,ynum,nslices=0,tpoints;
float inplane;
float thick;
float regmatrix[4][4];
float x_min,x_max,y_min,y_max,z_min,z_max;
float *scalefacts;
int *briktypes;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -maskstem maskstem [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem     <instem>        input time series dataset (omit extension)\n");
  printf("                    vol (3d+t BRIK): <instem>.BRIK \n");
  printf("                  surf (surface ts): <instem>-{rh,lh}.wt\n");
  printf("    -maskstem   <maskstem>      w file containing non-zero values in ROI\n");
  printf("                   omit extension: <maskstem>-{rh,lh}.w\n");
  printf("\n");
  printf("  Optional parameters:  (defaults in []'s)\n");
  printf("    -datatype [vol]            input data type: vol or surf\n");
  printf("    -indir    [.]              dir containing BRIK file\n");
  printf("    -maskdir  [.]              dir containing mask w file\n");
  printf("    -outfile  [instem-avg.txt] output file name\n");
  printf("    -outdir   [.]              output dir\n");
  printf("    -hemi     [rh]             hemisphere (rh or lh)\n");
  printf("    -quiet                     suppress messages\n");
  printf("\n");
  printf("  Parameters specific to vol data type:\n");
  printf("    -regdat   [register.dat]   registration file\n");
  printf("                                 (full path name or relative to indir)\n");
  printf("    -dpaint   [0.0]            dist (mm) from surface to project along normal\n");
  printf("\n");
  printf("  Parameters specific to surf data type:\n");
  printf("    -name     <subjname>       subject name (required)\n");
  printf("    -sparseavg                 average values from non-zero vertices only\n");
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
      if (MATCH(argv[i],"-datatype") && i+1<argc) {
        strcpy(datatype,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-name") && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-maskstem") && i+1<argc) {
        strcpy(maskstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-maskdir") && i+1<argc) {
        strcpy(maskdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outfile") && i+1<argc){
        strcpy(outfile,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-regdat") && i+1<argc) {
        strcpy(regdat,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-dpaint")) && i+1<argc) {
        dpaint = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-sparseavg")){
        sparseavgflag = 1;
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
  if (MATCH(datatype,"vol"))
    voldataflag = 1;
  else if (MATCH(datatype,"surf")) {
    voldataflag = 0;
    if (MATCH(subj,UNDEFSTR))
      ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
                progname);
  } else
    ErrorExit(ERROR_BADPARM,"%s: ### datatype %s not supported ...quitting\n",
                progname,datatype);
  if (MATCH(instem,UNDEFSTR))
    ErrorExit(ERROR_BADPARM,"%s: ### instem not specified ...quitting\n",
              progname);
  if (MATCH(maskstem,UNDEFSTR))
    ErrorExit(ERROR_BADPARM,"%s: ### maskstem not specified ...quitting\n",
              progname);
  if (MATCH(outfile,UNDEFSTR))
    sprintf(outfile,"%s-avg.txt",instem);
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh"))
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  if (!FileExists(indir))
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",
              progname,indir);
  if(!isadir(indir))
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
  if (!FileExists(maskdir))
    ErrorExit(ERROR_BADPARM,"%s: ### maskdir %s not found ...quitting\n",
              progname,maskdir);
  if(!isadir(indir))
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
  if (!FileExists(outdir))
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  if(!isadir(outdir))
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
  if (!FileExists(tempstr))
    ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
              progname,tempstr);
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

void read_HEAD(char *fname, int *xsize, int *ysize, int *zsize, int *tsize)
{
  FILE *fp;
  char line[STRLEN];
  char entry[STRLEN];
  char value[STRLEN];
  int i;

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
      MsgPrintf("%s: %d float_facs in %s\n",progname,*tsize,fname);
      scalefacts = (float *)calloc(*tsize,sizeof(float)); MTEST(scalefacts);
      for(i=0;i<*tsize;i++) {
        fscanf(fp,"%f",&scalefacts[i]);
        MsgPrintf("%s: float scale factor for subbrik %i: %f\n",progname,i,scalefacts[i]);
      }
    }
    if (MATCH(value,"BRICK_TYPES")) {
      fgets(line,99,fp);
      sscanf(line,"%s %*s %s",entry,value);
      *tsize = atoi(value);
      MsgPrintf("%s: %d brick_types in %s\n",progname,*tsize,fname);
      briktypes = (int *)calloc(*tsize,sizeof(int)); MTEST(briktypes);
      for(i=0;i<*tsize;i++) {
        fscanf(fp,"%d",&briktypes[i]);
        MsgPrintf("%s: brick data type for subbrik %i: %d\n",progname,i,briktypes[i]);
      }
    }
  }
  fclose(fp);
}

/* read x,y,z subbrik for a single time point */
void read_BRIK(float ***xyz, char *indir, char *instem, int subbrik,
               float scalefact, int briktype)
{
  FILE *fp;
  char tempstr[STRLEN];
  long offset=0;
  int z,y,x;
  float max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0,avg,stdev,f;
  int npix=0;

  if (scalefact==0) scalefact = 1;
  switch (briktype) {
    case BRIK_SHORT:
      offset = subbrik*nslices*ynum*xnum*sizeof(short);
      break;
    case BRIK_FLOAT:
      offset = subbrik*nslices*ynum*xnum*sizeof(float);
      break;
    default:
      ErrorExit(ERROR_BADFILE,"%s: ### briktype %d is currently unsupported ...quitting\n",
              progname,briktype);
  }  
  sprintf(tempstr,"%s/%s.BRIK",indir,instem);
  MsgPrintf("%s: reading subbrik %d of AFNI BRIK %s\n",
            progname,subbrik,tempstr);
  fp = fopen(tempstr,"rb");
  if (fp==NULL)
    ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
            progname,tempstr);
  fseek(fp,offset,SEEK_SET);  /* TODO: check file size */
  for (z=0;z<nslices;z++)
  for (y=0;y<ynum;y++)
  for (x=0;x<xnum;x++) {
    if (briktype==BRIK_SHORT) {
      f = xyz[z][y][x] = scalefact*(float)freadShort(fp);
    } else if (briktype==BRIK_FLOAT) {
      f = xyz[z][y][x] = scalefact*freadFloat(fp);
    } else
      f = 0;
    if (f>max) max = f;
    if (f<min) min = f;
    sum += f;
    sum2 += f*f;
    npix++;
  }
  avg = sum/(xnum*ynum);
  stdev = sqrt(sum2/(xnum*ynum)-avg*avg);

  MsgPrintf("%s: pixels(x*y*z)=%d, avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            progname,npix,avg,stdev,min,max);
}

void tranform(float x1, float y1, float z1, float *x2, float *y2,
                     float *z2, float M[4][4])
{
  *x2 = x1*M[0][0]+y1*M[0][1]+z1*M[0][2]+M[0][3];
  *y2 = x1*M[1][0]+y1*M[1][1]+z1*M[1][2]+M[1][3];
  *z2 = x1*M[2][0]+y1*M[2][1]+z1*M[2][2]+M[2][3];
}

float avg_surf_roi()
{
  static float x_orig, y_orig, z_orig;
  static float x_tran, y_tran, z_tran;
  static VERTEX *v;
  static int k,i,j,n;
  int nverts=0;
  float avg=0;

  for (k=0;k<mris->nvertices;k++) {
    if (maskvals[k]==0) continue;
    v = &mris->vertices[k];
    x_orig = v->x + dpaint*v->nx;
    y_orig = v->y + dpaint*v->ny;
    z_orig = v->z + dpaint*v->nz;
    tranform(x_orig,y_orig,z_orig,&x_tran,&z_tran,&y_tran,regmatrix);
    /* swapping y and z for some reason */

    /* should do trilinear interpolation here */
    n = (int)(floor((z_tran-z_min)/thick));
    i = (int)(ynum-1-(y_max-y_tran)/inplane);
    j = (int)((x_max-x_tran)/inplane);
    if (n >= 0 && n < nslices &&
        i >= 0 && i < ynum  && 
        j >= 0 && j < xnum) {
      avg += dat[n][ynum-1-i][j];
      nverts++;
    }
  }
  if (nverts==0) return 0;
  else return(avg/nverts);
}

int main(int argc, char **argv)
{
  char tempstr[STRLEN];
  FILE  *fp=NULL, *fp_wt=NULL;
  int i,j,k,t; /*y,x,z,time*/
  float f,avg;
  int numsurfvals,n;

  parse_args(argc,argv);

  if (voldataflag) {
    /* read register.dat */
    if (regdat[0]=='/')
      sprintf(tempstr,"%s",regdat);
    else
      sprintf(tempstr,"%s/%s",indir,regdat);
    fp = fopen(tempstr,"r");
    if (fp==NULL)
      ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                progname,tempstr);
    fscanf(fp,"%s",subj);
    fscanf(fp,"%f",&inplane);
    fscanf(fp,"%f",&thick);
    fscanf(fp,"%*f");
    for (i=0;i<4;i++)
    for (j=0;j<4;j++) {
      fscanf(fp,"%f",&regmatrix[i][j]);
    }
    fclose(fp);
  }

  /* load surface */
  mris = openSurface(subj,hemi,"orig");

  /* get dimensions of dataset */
  if (voldataflag) {
    sprintf(tempstr,"%s/%s.HEAD",indir,instem);
    MsgPrintf("%s: reading AFNI header %s\n",progname,tempstr);
    read_HEAD(tempstr,&xnum,&ynum,&nslices,&tpoints);
  } else {
    sprintf(tempstr,"%s/%s-%s.wt",indir,instem,hemi);
    MsgPrintf("%s: opening time series file %s\n",progname,tempstr);
    fp_wt = fopen(tempstr,"r");
    if (fp_wt==NULL)
      ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,tempstr);
    fread2(&tpoints,fp_wt);
    MsgPrintf("%s: number of surface time points: %d\n",progname,tpoints);
  }

  /* check mask file */
  sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
  i = maxVert(tempstr)+1;
  if (i==0)
    ErrorExit(ERROR_BADFILE,"%s: ### error reading surface mask file %s ...quitting\n",
              progname,tempstr);
  if (i > mris->nvertices)
    ErrorExit(ERROR_BADFILE,"%s: ### mask file %s has vertex number %d > nvertices %d...quitting\n",
              progname,tempstr,i,mris->nvertices);

  /* allocate memory for data */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  dat = (float ***)calloc(nslices,sizeof(float **));         MTEST(dat);
  for (k=0;k<nslices;k++) {
    dat[k] = (float **)calloc(ynum,sizeof(float *));         MTEST(*dat);
    for (i=0;i<ynum;i++) {
      dat[k][i] = (float *)calloc(xnum,sizeof(float));       MTEST(**dat);
    }
  }
  maskvals = (float *)calloc(mris->nvertices,sizeof(float)); MTEST(maskvals);
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* read mask */
  sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
  if(readSurfVals(tempstr,maskvals,mris->nvertices)==-1)
    ErrorExit(ERROR_BADFILE,"%s: ### error reading surface threshold values...quitting\n",
              progname);

  sprintf(tempstr,"%s/%s",outdir,outfile);
  fp = fopen(tempstr,"w");
  if (fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### unable to create output file %s ...quitting\n",
              progname,tempstr);
  if (voldataflag) {
    /* define bounds for painting */
    z_min = -thick*nslices/2.0;
    z_max = thick*nslices/2.0;
    x_min = -inplane*xnum/2.0;
    x_max = inplane*xnum/2.0;
    y_min = -inplane*ynum/2.0;
    y_max = inplane*ynum/2.0;

    MsgPrintf("%s: reading brik and calculating averages\n",progname);
    for(t=0;t<tpoints;t++) {
      /* read one time point from brik */
      read_BRIK(dat,indir,instem,t,scalefacts[t],briktypes[t]);
      avg=avg_surf_roi();
/*      MsgPrintf("%s: avg(%d) = %f\n",progname,t,avg);*/
      fprintf(fp,"%f\n",avg);
    }
  } else {
    /* read values from wt file for each vertex at each timepoint */
    MsgPrintf("%s: reading wt file and calculating averages\n",progname);
    for(t=0;t<tpoints;t++) {
      fread3(&numsurfvals,fp_wt);
/*      MsgPrintf("%s: number of surface values for time point %d: %d\n",
                progname,t,numsurfvals);*/
      for (i=0,n=0;i<numsurfvals;i++) {
        fread3(&k,fp_wt);
        fread4(&f,fp_wt);
        if (maskvals[k]!=0) {
          avg += f;
          n++;
        }
      }
      if (sparseavgflag)
        avg /= n;
      else {
        for(k=0,n=0;k<mris->nvertices;k++)
          if (maskvals[k]!=0) n++;
        avg /= n;
      }
/*      MsgPrintf("%s: avg(%d) = %f\n",progname,t,avg);*/
      fprintf(fp,"%f\n",avg);
    }
    fclose(fp_wt);
  }
  fclose(fp);

  exit(0);
}

