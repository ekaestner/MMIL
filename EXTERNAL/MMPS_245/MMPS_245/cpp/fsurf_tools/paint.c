/* paint.c: read brik or bfloats and paint to surface
      created: ? MS
     last mod: 10/01/04 DH

   purpose:
     sampling 3d stats onto 3d surface

   input:
     afni BRIK + HEAD files or freesurfer bfloat and hdr files
     register.dat file

   optional input:
     threshold dataset

   output:
     stats file (w file)
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
char  subj[STRLEN]=UNDEFSTR;
char  type[STRLEN]="BRIK";
int   imgindex=0;
int   complexflag=0;
char  instem_i[STRLEN]=UNDEFSTR;
int   imgindex_i=0;
int   tseriesflag=0;
int   threshflag=0;
int   threshabsflag=0;
int   threshsurfflag=0;
char  threshstem[STRLEN]=UNDEFSTR;
char  threshtype[STRLEN]="BRIK";
int   threshindex=0;
float threshval=0;
char  indir[STRLEN]=".";
char  outstem[STRLEN]=UNDEFSTR;
char  outstem_i[STRLEN]=UNDEFSTR;
char  outdir[STRLEN]=".";
char  hemi[STRLEN]="rh";
char  regdat[STRLEN]="register.dat";
char  paintsurf[STRLEN]="orig";
float dmax=0;
float dstep=0.25;
int   truncflag=0;
float minphase=0.25;
float maxphase=0.75;
float slope=10.0;

/* other global variables */
static char *progname = NULL;
MRIS *mris;
float ***dat_r=NULL;
float ***dat_i=NULL;
float ***threshdat=NULL;
float *threshvals=NULL;
int xnum,ynum,nslices,tpoints;
float inplane;
float thick;
float regmatrix[4][4];

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem      <instem>       omit extension: <instem>_000.bfloat\n");
  printf("                                            or: <instem>.BRIK\n");
  printf("\n");
  printf("  Optional parameters:  (defaults in []'s)\n");
  printf("    -subj                       freesurfer subject name (overrides register.dat)\n");
  printf("    -type        [BRIK]         input type: BRIK or bfloat\n");
  printf("    -imgindex    [0]            sub-BRIK index or imageoffset to paint\n");
  printf("    -complex                    complex data (real and imaginary)\n");
  printf("    -instem_i    <instem_i>     instem for imaginary dataset (complex only)\n");
  printf("    -imgindex_i  [0]            index or offset of imaginary dataset\n");
  printf("    -timeseries                 ignore imgindex and paint all to wt file\n");
  printf("    -thresh                     threshold instem data using threshstem data\n");
  printf("    -threshabs                  threshold using absolute values (+/-)\n");
  printf("    -threshstem  [instem]       stem for threshold dataset\n");
  printf("    -threshtype  [BRIK]         type of threshold dataset (BRIK, bfloat, or w)\n");
  printf("    -threshindex [0]            index or offset for threshold dataset\n");
  printf("    -threshval   [0.0]          minimum value for threshold data\n");
  printf("    -indir       [.]            input dir\n");
  printf("    -outstem     [instem]       output file stem: <outstem>-rh.w\n");
  printf("    -outstem_i   [instem_i]     outstem for imaginary data (complex)\n");
  printf("    -outdir      [.]            output dir\n");
  printf("    -hemi        [rh]           hemisphere (rh or lh)\n");
  printf("    -regdat      [register.dat] name of register file\n");
  printf("    -paintsurf   [orig]         surface to paint onto\n");
  printf("    -dmax        [0.0]          dist (mm) to project along normal\n");
  printf("    -dstep       [0.25]         size (mm) of projection step\n");
  printf("    -truncphase                 zeroes values with phase out of range (complex)\n");
  printf("    -minphase    [0.25]         start of preserved phases (for truncphase)\n");
  printf("    -maxphase    [0.75]         end of preserved phases (for truncphase)\n");
  printf("    -slope       [10.0]         slope of fall-off (for truncphase)\n");
  printf("    -quiet                      suppress messages\n");
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
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-type") && i+1<argc) {
        strcpy(type,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-imgindex")) && i+1<argc) {
        imgindex = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-complex")){
        complexflag = 1;
      } else
      if (MATCH(argv[i],"-instem_i") && i+1<argc) {
        strcpy(instem_i,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-imgindex_i")) && i+1<argc) {
        imgindex = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-timeseries")){
        tseriesflag = 1;
      } else
      if (MATCH(argv[i],"-thresh")){
        threshflag = 1;
      } else
      if (MATCH(argv[i],"-threshabs")){
        threshflag = 1;
        threshabsflag = 1;
      } else
      if (MATCH(argv[i],"-threshstem") && i+1<argc) {
        strcpy(threshstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-threshtype") && i+1<argc) {
        strcpy(threshtype,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-threshindex")) && i+1<argc) {
        threshindex = atoi(argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-threshval")) && i+1<argc) {
        threshval = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem_i") && i+1<argc){
        strcpy(outstem_i,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-regdat") && i+1<argc) {
        strcpy(regdat,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-paintsurf") && i+1<argc) {
        strcpy(paintsurf,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-dmax")) && i+1<argc) {
        dmax = atof(argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-dstep")) && i+1<argc) {
        dstep = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-truncphase")){
        truncflag = 1;
      } else
      if ((MATCH(argv[i],"-minphase")) && i+1<argc) {
        minphase = atof(argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-maxphase")) && i+1<argc) {
        maxphase = atof(argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-slope")) && i+1<argc){
        slope = atof(argv[i+1]); i++;
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
  if (MATCH(instem_i,UNDEFSTR) && complexflag) {
    sprintf(instem_i,"%s_i",instem);
    sprintf(instem,"%s_r",instem);
  }
  if (MATCH(threshstem,UNDEFSTR))
    strcpy(threshstem,instem);
  if (!complexflag) {
    if (MATCH(outstem,UNDEFSTR))
      strcpy(outstem,instem);
  } else {
    if (MATCH(outstem_i,UNDEFSTR)) {
      if (MATCH(outstem,UNDEFSTR)) {
        strcpy(outstem,instem);
        strcpy(outstem_i,instem_i);
      } else {
        sprintf(outstem_i,"%s_i",outstem);
        sprintf(outstem,"%s_r",outstem);
      }
    } else {
      if (MATCH(outstem,UNDEFSTR))
        strcpy(outstem,instem);
      /* if outstem and outstem_i are both defined, do nothing */
    }
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh"))
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
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
  if (MATCH(type,"BRIK")) {
    sprintf(tempstr,"%s/%s.BRIK",indir,instem);
    if (!FileExists(tempstr))
      ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                progname,tempstr);
    sprintf(tempstr,"%s/%s.HEAD",indir,instem);
    if (!FileExists(tempstr))
      ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                progname,tempstr);
    if(complexflag) {
      sprintf(tempstr,"%s/%s.BRIK",indir,instem_i);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
      sprintf(tempstr,"%s/%s.HEAD",indir,instem_i);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
    }
  } else if (MATCH(type,"bfloat")) {
    nslices = getNumBfloats(indir,instem,NULL);
    if(nslices<=0)
      ErrorExit(ERROR_BADPARM,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
             progname,indir,instem);
    if(complexflag) {
      i = getNumBfloats(indir,instem_i,NULL);
      if(i<=0)
        ErrorExit(ERROR_BADPARM,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
               progname,indir,instem_i);
      if(i!=nslices)
        ErrorExit(ERROR_BADPARM,"%s: ### number of real and imaginary bfloats (%d and %d) do not match ...quitting\n",
               progname,nslices,i);
    }
  } else
    ErrorExit(ERROR_BADPARM,"%s: ### type %s not supported ...quitting\n",
              progname,type);
  
  if (threshflag) {
    if (MATCH(threshtype,"BRIK")) {
      sprintf(tempstr,"%s/%s.BRIK",indir,threshstem);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
      sprintf(tempstr,"%s/%s.HEAD",indir,threshstem);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
    } else if (MATCH(threshtype, "bfloat")) {
      i = getNumBfloats(indir,threshstem,NULL);
      if(i<=0)
        ErrorExit(ERROR_BADPARM,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
               progname,indir,threshstem);
      if(i!=nslices)
        ErrorExit(ERROR_BADPARM,"%s: ### number of instem and threshstem bfloats (%d and %d) do not match ...quitting\n",
               progname,nslices,i);
    } else if (MATCH(threshtype, "w")) {
      sprintf(tempstr,"%s/%s-%s.w",indir,threshstem,hemi);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
      threshsurfflag=1; /* apply threshold after painting  */
    } else
      ErrorExit(ERROR_BADPARM,"%s: ### threshtype %s not supported ...quitting\n",
                progname,threshtype);
  }
  if (imgindex < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### bad imgindex: %d => must be >= 0 ...quitting\n",
              progname,imgindex);
  if (imgindex_i < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### bad imgindex_i: %d => must be >= 0 ...quitting\n",
              progname,imgindex_i);
  if (threshindex < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### bad threshindex: %d => must be >= 0 ...quitting\n",
              progname,threshindex);
  if (dmax < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### bad dmax: %f => must be >= 0 ...quitting\n",
              progname,dmax);
  if (dstep <= 0)
    ErrorExit(ERROR_BADPARM,"%s: ### bad dstep: %f => must be > 0 ...quitting\n",
              progname,dstep);
  if (minphase < 0 || minphase > 1)
    ErrorExit(ERROR_BADPARM,"%s: ### bad minphase: %f => must be between 0.0 and 1.0 ...quitting\n",
              progname,minphase);
  if (maxphase < 0 || maxphase > 1)
    ErrorExit(ERROR_BADPARM,"%s: ### bad maxphase: %f => must be between 0.0 and 1.0 ...quitting\n",
              progname,maxphase);
  if (minphase == maxphase)
    ErrorExit(ERROR_BADPARM,"%s: ### bad maxphase: %f => maxphase must not equal minphase ...quitting\n",
              progname,maxphase);
  if (slope < 0) slope = -slope;
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

/* read x,y,z subbrik for a single time point */
void read_BRIK(float ***xyz, char *indir, char *instem, int subbrik,
               float ffact, int btype, char *byteorder)
{
  FILE *fp;
  char tempstr[STRLEN];
  long offset=0;
  int z,y,x;
  double sum=0,sum2=0,avg,stdev;
  float max=-BIGFLOAT,min=BIGFLOAT,f;
  int npix=0;
  short s;

  if (ffact==0) ffact = 1;
  switch (btype) {
    case BRIK_SHORT:
      offset = subbrik*nslices*ynum*xnum*sizeof(short);
      break;
    case BRIK_FLOAT:
      offset = subbrik*nslices*ynum*xnum*sizeof(float);
      break;
    default:
      ErrorExit(ERROR_BADFILE,"%s: ### btype %d is currently unsupported ...quitting\n",
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
  for (z=0;z<nslices;z++)
  for (y=0;y<ynum;y++)
  for (x=0;x<xnum;x++) {
    if (btype==BRIK_SHORT) {
      s = freadShort(fp);
      if(MATCH(byteorder,"LSB")) s = swapShort(s);
      f = xyz[z][y][x] = ffact*(float)s;
    } else if (btype==BRIK_FLOAT) {
      f = freadFloat(fp);
      if(MATCH(byteorder,"LSB")) f = swapFloat(f);
      f = xyz[z][y][x] = ffact*f;
    } else
      f = 0;
    if (f>max) max = f;
    if (f<min) min = f;
    sum += f;
    sum2 += f*f;
    npix++;
  }
  avg = sum/npix;
  stdev = sqrt(sum2/npix-avg*avg);
  MsgPrintf("%s: pixels(x*y*z)=%d, avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            progname,npix,avg,stdev,min,max);
}

void read_bfloats(float ***stats, char *indir, char *instem, int imageoffset)
{
  FILE *fp;
  char tempstr[STRLEN];
  long offset;
  int i,j,k,num;
  float f;
  double sum=0,sum2=0,max=-1000,min=1000;

  for (k=0;k<nslices;k++) {
    sprintf(tempstr,"%s/%s_%03d.bfloat",indir,instem,k);
    MsgPrintf("%s: reading sub-image %d of bfloat %s\n",
              progname,imageoffset,tempstr);
    fp = fopen(tempstr,"rb");
    if (fp==NULL)
      ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
              progname,tempstr);
    offset = imageoffset*ynum*xnum*sizeof(float);
    fseek(fp,offset,SEEK_SET);  /* TODO: check file size */
    for (i=0;i<ynum;i++)
    for (j=0;j<xnum;j++) {
      f = stats[k][i][j] = freadFloat(fp);
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
    fclose(fp);
    num=xnum*ynum;
    sum /= num;
    sum2 = sqrt(sum2/num-sum*sum);
    MsgPrintf("%s: file %s read (%dx%d)\n",progname,tempstr,xnum,ynum);
    MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);
  }
}

void tranform(float x1, float y1, float z1, float *x2, float *y2,
                     float *z2, float M[4][4])
{
  *x2 = x1*M[0][0]+y1*M[0][1]+z1*M[0][2]+M[0][3];
  *y2 = x1*M[1][0]+y1*M[1][1]+z1*M[1][2]+M[1][3];
  *z2 = x1*M[2][0]+y1*M[2][1]+z1*M[2][2]+M[2][3];
}

void paint_surface()
{
  float x_min,x_max,y_min,y_max,z_min,z_max;
  float x_orig, y_orig, z_orig;
  float x_tran, y_tran, z_tran;
  float fa, fr, fi;
  VERTEX *v;
  float dist, famax, frmax, fimax;
  int k,i,j,n;

  /* define bounds */
  z_min = -thick*nslices/2.0;
  z_max = thick*nslices/2.0;
  x_min = -inplane*xnum/2.0;
  x_max = inplane*xnum/2.0;
  y_min = -inplane*ynum/2.0;
  y_max = inplane*ynum/2.0;

  for (k=0;k<mris->nvertices;k++) {
    v = &mris->vertices[k];
    famax = frmax = fimax = 0;
    for (dist=0; dist<=dmax; dist+=dstep) {
      x_orig = v->x + dist*v->nx;
      y_orig = v->y + dist*v->ny;
      z_orig = v->z + dist*v->nz;
      tranform(x_orig,y_orig,z_orig,&x_tran,&z_tran,&y_tran,regmatrix);
      /* swapping y and z for some reason */

      /* should do trilinear interpolation here */
      n = (int)(floor((z_tran-z_min)/thick));
      i = (int)(ynum-1-(y_max-y_tran)/inplane);
      j = (int)((x_max-x_tran)/inplane);
      if (n >= 0 && n < nslices &&
          i >= 0 && i < ynum  && 
          j >= 0 && j < xnum) {
        fr = dat_r[n][ynum-1-i][j];
        if (complexflag) fi = dat_i[n][ynum-1-i][j];
      }
      else fr=fi=0;
      if (complexflag) {
        fa = hypot(fr,fi);
        if (fabs(fa) > fabs(famax)) {famax=fa; frmax=fr; fimax=fi;}
      } else {
        if(fabs(fr) > fabs(frmax)) frmax=fr;
      }
    } /*dist*/
    mris->vertices[k].val = frmax;
    if (complexflag) mris->vertices[k].imag_val = fimax;
  }
}

int main(int argc, char **argv)
{
  char tempstr[STRLEN];
  int xnum_tmp, ynum_tmp, nslices_tmp, tpoints_tmp;
  FILE  *fp, *fp_i=NULL;
  int i,j,k,t; /*y,x,z,t*/
  int nonzero_verts;
  float max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0,avg,stdev,f;
  float max_i=-BIGFLOAT,min_i=BIGFLOAT,sum_i=0,sum2_i=0;
  float *scalefacts=NULL, *scalefacts_i=NULL, *scalefacts_t=NULL;
  int *briktypes=NULL, *briktypes_i=NULL, *briktypes_t=NULL;
  char byteorder[STRLEN]=UNDEFSTR;
  char byteorder_i[STRLEN]=UNDEFSTR;
  char byteorder_t[STRLEN]=UNDEFSTR;
  char tmp_subj[STRLEN]=UNDEFSTR;

  parse_args(argc,argv);

  /* read register.dat */
  if(regdat[0]=='/')
    sprintf(tempstr,"%s",regdat);
  else if(regdat[0]=='.' && regdat[1]=='/')
    sprintf(tempstr,"%s",regdat);
  else
    sprintf(tempstr,"%s/%s",indir,regdat);
  MsgPrintf("%s: reading %s\n",progname,tempstr);
  fp = fopen(tempstr,"r");
  if (fp==NULL)
    ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
              progname,tempstr);
  fscanf(fp,"%s",tmp_subj);
  if (MATCH(subj,UNDEFSTR)) strcpy(subj,tmp_subj);
  fscanf(fp,"%f",&inplane);
  fscanf(fp,"%f",&thick);
  fscanf(fp,"%*f");
  for (i=0;i<4;i++)
  for (j=0;j<4;j++) {
    fscanf(fp,"%f",&regmatrix[i][j]);
  }
  fclose(fp);

  /* load surface */
  mris = openSurface(subj,hemi,paintsurf);

  /* get dimensions of dataset */
  if (MATCH(type,"BRIK")) {
    sprintf(tempstr,"%s/%s.HEAD",indir,instem);
    MsgPrintf("%s: reading AFNI header %s\n",progname,tempstr);
    read_HEAD(tempstr,&xnum,&ynum,&nslices,&tpoints,
              &scalefacts,&briktypes,byteorder);
  } else if (MATCH(type,"bfloat")) {
    sprintf(tempstr,"%s/%s_000.hdr",indir,instem);
    MsgPrintf("%s: reading bfloat header %s\n",progname,tempstr);
    readFloatHeader(tempstr,&xnum,&ynum,&tpoints);
  }
  if (!tseriesflag && imgindex>=tpoints)
    ErrorExit(ERROR_BADPARM,"%s: ### image offset (%d) >= # time points (%d) ...quitting\n",
              progname,imgindex,tpoints);

  /* make sure dimensions match for imaginary */
  if (complexflag) {
    if (MATCH(type,"BRIK")) {
      sprintf(tempstr,"%s/%s.HEAD",indir,instem_i);
      MsgPrintf("%s: reading AFNI header %s\n",progname,tempstr);
      read_HEAD(tempstr,&xnum_tmp,&ynum_tmp,&nslices_tmp,&tpoints_tmp,
                &scalefacts_i,&briktypes_i,byteorder_i);
      if(xnum_tmp!=xnum || ynum_tmp!=ynum || nslices_tmp!=nslices)
        ErrorExit(ERROR_BADFILE,
          "%s: ### image dimensions for imaginary do not match real...quitting\n",
          progname);
    } else if (MATCH(type,"bfloat")) {
      sprintf(tempstr,"%s/%s_000.hdr",indir,instem_i);
      MsgPrintf("%s: reading bfloat header %s\n",progname,tempstr);
      readFloatHeader(tempstr,&xnum_tmp,&ynum_tmp,&tpoints_tmp);
      if(xnum_tmp!=xnum || ynum_tmp!=ynum)
        ErrorExit(ERROR_BADFILE,
          "%s: ### image dimensions for imaginary do not match real...quitting\n",
          progname);
    }
    if (!tseriesflag && imgindex_i>=tpoints_tmp) {
      ErrorExit(ERROR_BADPARM,"%s: ### image offset (%d) >= # time points (%d) ...quitting\n",
                progname,imgindex_i,tpoints_tmp);
    } else if (tseriesflag && tpoints!=tpoints_tmp) {
      ErrorExit(ERROR_BADFILE,
        "%s: ### number of time points for imaginary do not match real...quitting\n",
        progname);
    }
  }
  /* make sure dimensions match for thresh dataset */
  if (threshflag) {
    if (MATCH(threshtype,"w")) {
      sprintf(tempstr,"%s/%s-%s.w",indir,threshstem,hemi);
      i = maxVert(tempstr)+1;
      if (i==0)
        ErrorExit(ERROR_BADFILE,"%s: ### error reading surface threshold file %s ...quitting\n",
                  progname,tempstr);
      if (i > mris->nvertices)
        ErrorExit(ERROR_BADFILE,"%s: ### threshold file %s has vertex number %d > nvertices %d...quitting\n",
                  progname,tempstr,i,mris->nvertices);
    } else {
      if (MATCH(threshtype,"BRIK")) {
        sprintf(tempstr,"%s/%s.HEAD",indir,threshstem);
        MsgPrintf("%s: reading AFNI header %s\n",progname,tempstr);
        read_HEAD(tempstr,&xnum_tmp,&ynum_tmp,&nslices_tmp,&tpoints_tmp,
                  &scalefacts_t,&briktypes_t,byteorder_t);
        if(xnum_tmp!=xnum || ynum_tmp!=ynum || nslices_tmp!=nslices)
          ErrorExit(ERROR_BADFILE,
            "%s: ### image dimensions for threshstem do not match instem...quitting\n",
            progname);
      } else if (MATCH(threshtype,"bfloat")) {
        sprintf(tempstr,"%s/%s_000.hdr",indir,threshstem);
        MsgPrintf("%s: reading bfloat header %s\n",progname,tempstr);
        readFloatHeader(tempstr,&xnum_tmp,&ynum_tmp,&tpoints_tmp);
        if(xnum_tmp!=xnum || ynum_tmp!=ynum)
          ErrorExit(ERROR_BADFILE,
            "%s: ### image dimensions for threshstem do not match instem...quitting\n",
            progname);
      }
      if (threshindex>=tpoints_tmp)
        ErrorExit(ERROR_BADPARM,"%s: ### image offset (%d) >= # time points (%d) ...quitting\n",
                progname,threshindex,tpoints_tmp);
    }
  }

  /* allocate memory for data */
  MsgPrintf("%s: starting to allocate memory\n",progname);
  dat_r = (float ***)calloc(nslices,sizeof(float **));         MTEST(dat_r);
  for (k=0;k<nslices;k++) {
    dat_r[k] = (float **)calloc(ynum,sizeof(float *));         MTEST(*dat_r);
    for (i=0;i<ynum;i++) {
      dat_r[k][i] = (float *)calloc(xnum,sizeof(float));       MTEST(**dat_r);
    }
  }
  if (complexflag) {
    dat_i = (float ***)calloc(nslices,sizeof(float **));       MTEST(dat_i);
    for (k=0;k<nslices;k++) {
      dat_i[k] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_i);
      for (i=0;i<ynum;i++) {
        dat_i[k][i] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_i);
      }
    }
  }
  if (threshflag) {
    if (threshsurfflag) {
      threshvals = (float *)calloc(mris->nvertices,sizeof(float)); MTEST(threshvals);
    } else {
      threshdat = (float ***)calloc(nslices,sizeof(float **));     MTEST(threshdat);
      for (k=0;k<nslices;k++) {
        threshdat[k] = (float **)calloc(ynum,sizeof(float *));     MTEST(*threshdat);
        for (i=0;i<ynum;i++) {
          threshdat[k][i] = (float *)calloc(xnum,sizeof(float));   MTEST(**threshdat);
        }
      }
    }
  }
  MsgPrintf("%s: finished allocating memory\n",progname);

  /* read threshold data */
  if (threshflag) {
    if(MATCH(threshtype,"BRIK")) {
      read_BRIK(threshdat,indir,threshstem,threshindex,
                scalefacts_t[threshindex],briktypes_t[threshindex],byteorder_t);
    } else if (MATCH(threshtype,"bfloat")) {
      read_bfloats(threshdat,indir,threshstem,threshindex);
    } else if (MATCH(threshtype,"w")) {
      sprintf(tempstr,"%s/%s-%s.w",indir,threshstem,hemi);
      if(readSurfVals(tempstr,threshvals,mris->nvertices)==-1)
        ErrorExit(ERROR_BADFILE,"%s: ### error reading surface threshold values...quitting\n",
                  progname);
    }
  }

  if (tseriesflag) {
    /* open output file(s) */
    sprintf(tempstr,"%s/%s-%s.wt",outdir,outstem,hemi);
    fp = fopen(tempstr,"w");
    if(fp==NULL){
      MsgPrintf("%s: ### can't create file %s\n",progname,tempstr);
      return(-1);
    }
    /* write number of tpoints at start of file */
    fwrite2(tpoints,fp);

    if (complexflag) {
      sprintf(tempstr,"%s/%s-%s.wt",outdir,outstem_i,hemi);
      fp_i = fopen(tempstr,"w");
      if(fp==NULL){
        MsgPrintf("%s: ### can't create file %s\n",progname,tempstr);
        return(-1);
      }
      fwrite2(tpoints,fp_i);
    }

    /* loop through each time point */
    for (t=0;t<tpoints;t++) {
      /* read brik or bfloat */
      MsgPrintf("%s: reading data set\n",progname);
      if(MATCH(type,"BRIK")) {
        read_BRIK(dat_r,indir,instem,t,
                  scalefacts[t],briktypes[t],byteorder);
        if(complexflag)
          read_BRIK(dat_i,indir,instem_i,t,
                  scalefacts_i[t],briktypes_i[t],byteorder_i);
      } else if (MATCH(type,"bfloat")) {
        read_bfloats(dat_r,indir,instem,t);
        if(complexflag)
          read_bfloats(dat_i,indir,instem_i,t);
      }

      /* apply 3d threshold */
      if (threshflag && !threshsurfflag) {
        MsgPrintf("%s: thresholding 3d data\n",progname);
        for (k=0;k<nslices;k++)
        for (i=0;i<ynum;i++)
        for (j=0;j<xnum;j++)
          if (threshabsflag) {
            if (fabs(threshdat[k][i][j]) < fabs(threshval)) {
              dat_r[k][i][j]=0;
              if (complexflag)
                dat_i[k][i][j]=0;
            }
          } else {
            if (threshdat[k][i][j] < threshval) {
              dat_r[k][i][j]=0;
              if (complexflag)
                dat_i[k][i][j]=0;
            }
          }
      }

      /* paint to surface */
      MsgPrintf("%s: painting time point %d onto %s surface from subject %s\n",
                progname,t,paintsurf,subj);
      paint_surface();

      /* apply surface threshold */
      if (threshsurfflag) {
        MsgPrintf("%s: thresholding surface data\n",progname);
        for (k=0;k<mris->nvertices;k++) {
          if (threshabsflag) {
            if (fabs(threshvals[k]) < fabs(threshval)) {
              mris->vertices[k].val = 0;
              if (complexflag) mris->vertices[k].imag_val = 0;
            }
          } else {
            if (threshvals[k] < threshval) {
              mris->vertices[k].val = 0;
              if (complexflag) mris->vertices[k].imag_val = 0;
            }
          }
        }
      }

      /* write to wt file(s) */
      MsgPrintf("%s: writing values for time point %d\n",progname,t);
      sum = sum2 = 0;
      max = -BIGFLOAT;
      min = BIGFLOAT;
      /* write number of non-zero vertices at start of time-point */
      for (k=0,nonzero_verts=0;k<mris->nvertices;k++) {
        if (mris->vertices[k].val!=0 || 
           (complexflag && mris->vertices[k].imag_val!=0))
          nonzero_verts++;
      }
      fwrite3(nonzero_verts,fp);
      if (complexflag) fwrite3(nonzero_verts,fp_i);
      for (k=0;k<mris->nvertices;k++) {
        if (mris->vertices[k].val!=0 || 
           (complexflag && mris->vertices[k].imag_val!=0)) {
          fwrite3(k,fp);
          f = mris->vertices[k].val;
          fwriteFloat(f, fp);
          sum += f;
          sum2 += f*f;
          if (f>max) max=f;
          if (f<min) min=f;
          if (complexflag) {
            fwrite3(k,fp_i);
            f = mris->vertices[k].imag_val;
            fwriteFloat(f, fp_i);
            sum_i += f;
            sum2_i += f*f;
            if (f>max_i) max_i=f;
            if (f<min_i) min_i=f;
          }
        }
      }
      avg = sum/mris->nvertices;
      stdev = sqrt(sum2/mris->nvertices-avg*avg);
      MsgPrintf("%s: %d of %d vertices, avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
                progname,nonzero_verts,mris->nvertices,avg,stdev,min,max);
      if (complexflag) {
        avg = sum_i/mris->nvertices;
        stdev = sqrt(sum2_i/mris->nvertices-avg*avg);
        MsgPrintf("%s: imag: %d of %d vertices, avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
                  progname,nonzero_verts,mris->nvertices,avg,stdev,min_i,max_i);
      }
    }
    fclose(fp);
    if (complexflag) fclose(fp_i);
  } else {
    /* read brik or bfloat */
    MsgPrintf("%s: reading data set\n",progname);

    if(MATCH(type,"BRIK")) {
      read_BRIK(dat_r,indir,instem,imgindex,
                scalefacts[imgindex],briktypes[imgindex],byteorder);
      if(complexflag)
        read_BRIK(dat_i,indir,instem_i,imgindex_i,
                scalefacts_i[imgindex],briktypes_i[imgindex],byteorder_i);
    } else if (MATCH(type,"bfloat")) {
      read_bfloats(dat_r,indir,instem,imgindex);
      if(complexflag)
        read_bfloats(dat_i,indir,instem_i,imgindex_i);
    }

    /* apply 3d threshold */
    if (threshflag && !threshsurfflag) {
      MsgPrintf("%s: thresholding 3d data\n",progname);
      for (k=0;k<nslices;k++)
      for (i=0;i<ynum;i++)
      for (j=0;j<xnum;j++) {
        if (threshabsflag) {
          if (fabs(threshdat[k][i][j]) < fabs(threshval)) {
            dat_r[k][i][j]=0;
            if (complexflag)
              dat_i[k][i][j]=0;
          }
        } else {
          if (threshdat[k][i][j] < threshval) {
            dat_r[k][i][j]=0;
            if (complexflag)
              dat_i[k][i][j]=0;
          }
        }
      }
    }

    /* paint to surface */
    MsgPrintf("%s: painting onto %s surface from subject %s\n",
      progname,paintsurf,subj);
    paint_surface();

    /* apply surface threshold */
    if (threshsurfflag) {
      MsgPrintf("%s: thresholding surface data\n",progname);
      for (k=0;k<mris->nvertices;k++) {
        if (threshabsflag) {
          if (fabs(threshvals[k]) < fabs(threshval)) {
            mris->vertices[k].val = 0;
            if (complexflag) mris->vertices[k].imag_val = 0;
          }
        } else {
          if (threshvals[k] < threshval) {
            mris->vertices[k].val = 0;
            if (complexflag) mris->vertices[k].imag_val = 0;
          }
        }
      }
    }

    /* write value file(s) */
    sprintf(tempstr,"%s/%s-%s",outdir,outstem,hemi);
    MsgPrintf("%s: writing values to file %s\n",progname,tempstr);
    MRISwriteValues(mris,tempstr);
    if (complexflag) {
      sprintf(tempstr,"%s/%s-%s",outdir,outstem_i,hemi);
      MsgPrintf("%s: writing values to file %s\n",progname,tempstr);
      MRISwriteImagValues(mris,tempstr);
    }
  }

  exit(0);
}

