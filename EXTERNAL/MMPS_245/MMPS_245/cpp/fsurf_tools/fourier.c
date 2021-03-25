/* fourier.c: read BRIKs or wt files containing time series or
              bfloats or wt files containing complex fourier series,
              calculate significance
     created:  11/18/04 DH (based on fourier by AD, MS)
     last mod: 05/11/05 DH

   purpose:
     generating fourier stats or fourier series

   input:
      afni BRIK + HEAD files containing time series
          or
      freesurfer bshort files containing time series
          or          
      painted time series file (wt)

          or

      complex (real and imaginary) freesurfer bfloat and hdr files
        with depth = number of frequencies
          or
      complex (real and imaginary) painted time series file (wt)
        with time points = number of frequencies

   output:
      bfloat and hdr files containing fourier stats
                                (r+i with amp = -log10(pval))
        or
      painted stats file (w) containing fourier stats

        and optionally

      complex (real and imaginary) freesurfer bfloat and hdr files
          with depth = number of frequencies
        or
      complex (real and imaginary) painted time series file (wt)
          with time points = number of frequencies
   
*/

/*
   todo: allow multiple time series datasets as input and sum fourier series
         across scans before calculating f-ratio
*/


#include "surflib.h"

#define MINARGC 3

#define BRIK_BYTE    0
#define BRIK_SHORT   1
#define BRIK_LONG    2
#define BRIK_FLOAT   3
#define BRIK_DOUBLE  4
#define BRIK_COMPLEX 5

#define MAXIM 2048
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* parameter defaults */
char datatype[STRLEN]="BRIK";
int voldataflag = 1;
char subj[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char outfsstem[STRLEN]=UNDEFSTR;
char real_infix[20]="_r";
char imag_infix[20]="_i";
char hemi[STRLEN]="rh";
float fstimfreq=8;
int stimfreq=8;
int lofreq = 3;
int hifreq = -1;
int omitfreq = -1;
int trendorder = 1;
int inputfsflag = 0;
int outputfsflag = 0;
int fratioflag = 0; /* output fratio, not significance */
int rawflag = 0; /* output raw real and imaginary for stim freq */
int numscans = 1;

/* other global variables */
static char *progname = NULL;
float ***tseries=NULL;      /* raw time series for a single slice */
float ***fseries_r=NULL;    /* real component of fourier series for a single slice */
float ***fseries_i=NULL;    /* imaginary component of fourier series for a single slice */
float **fsig_r=NULL;        /* real component of ftest sig vals + phase for a single slice */
float **fsig_i=NULL;        /* imaginary component of ftest sig vals + phase for a single slice */

int *surf_vnums=NULL;
float **surf_tseries=NULL;
float **surf_fseries_r=NULL;
float **surf_fseries_i=NULL;

int xnum,ynum,nslices,tpoints,numfreqs;
float ts[MAXIM];
float cfs[MAXIM];
float fs_r[MAXIM];
float fs_i[MAXIM];
float amp[MAXIM];
float phase[MAXIM];
float fsig_real,fsig_imag;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem   <instem>          time or fourier series dataset (omit extension)\n");
  printf("      time series (ts) types: bshort, BRIK, wt\n");
  printf("      fourier series (fs) types: bfloat, wt\n");
  printf("               BRIK ts: <instem>.BRIK \n");
  printf("             bshort ts: <instem>_%%03d.bshort\n");
  printf("             bfloat fs: <instem>_{r,i}_%%03d.bfloat\n");
  printf("      surface ts or fs: <instem>-{rh,lh}.wt\n");
  printf("\n");
  printf("  Optional parameters:  (defaults in []'s)\n");
  printf("    -datatype    [BRIK]         input data type\n");
  printf("       datatypes allowed: BRIK, bshort, bfloat, wtts, wtfs\n");
  printf("    -subj      <subjname>       required if data type = surf\n");
  printf("    -hemi        [rh]           hemisphere (rh or lh) (surf data types only)\n");
  printf("    -indir       [.]            dir containing input file\n");
  printf("    -outdir      [.]            output dir\n");
  printf("    -outstem     [instem-fourier]   file stem for complex stats output\n");
  printf("        volume: <outstem>_{r,i}_%%03d.bfloat\n");
  printf("       surface: <outstem>_{r,i}-{rh,lh}.w (if input is wt)\n");
  printf("    -outfsstem   <outfsstem>    file stem for output fourier series (optional)\n");
  printf("        volume: <outfsstem>_{r,i}_%%03d.bfloat\n");
  printf("       surface: <outfsstem>_{r,i}-{rh,lh}.wt\n");
  printf("    -stimfreq    [%d]            number of stim cycles per scan\n",stimfreq);
  printf("    -freqrange   [%d] [%d]       ignore freqs outside this range\n",lofreq,hifreq);
  printf("    -omitfreq    [%d]           ignore this freq for sig calc\n",omitfreq);
/*  printf("    -trendorder  [%d]           remove linear (1) or quadratic (2) trend\n",trendorder);*/
  printf("    -fratio                     output complex stats with old style amplitudes\n");
  printf("                                  (f-ratio not f-sig)\n");
  printf("    -raw                        output raw fourier components for stim freq\n");
  printf("                                  (neither f-ratio nor f-sig)\n");
  printf("                                  \"-raw\" takes precedence over \"-fratio\"\n");
  printf("    -numscans    [1]            number of dof's are multiplied by this number\n");
  printf("    -infixes     [%s] [%s]      real and imaginary infixes\n",real_infix,imag_infix);
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
      if (MATCH(argv[i],"-datatype") && i+1<argc) {
        strcpy(datatype,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-inputfs")) {
        /* for backwards compatibility */
        inputfsflag = 1;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outfsstem") && i+1<argc){
        strcpy(outfsstem,argv[i+1]); i++;
        outputfsflag = 1;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-stimfreq") && i+1<argc){
        fstimfreq = atof(argv[i+1]); i++;
        stimfreq = (int)rint(fstimfreq);
      } else
      if (MATCH(argv[i],"-freqrange") && i+2<argc){
        lofreq = (int)rint(atof(argv[i+1])); i++;
        hifreq = (int)rint(atof(argv[i+1])); i++;
      } else
      if (MATCH(argv[i],"-omitfreq") && i+1<argc){
        omitfreq = (int)rint(atof(argv[i+1])); i++;
      } else
      if (MATCH(argv[i],"-trendorder") && i+1<argc){
        trendorder = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-fratio")) {
        /* to revert to old behavior of setting amplitudes based on raw fratio */
        fratioflag = 1;
      } else
      if (MATCH(argv[i],"-raw")) {
        /* to output raw, unweighted, real and imaginary components for stim freq */
        rawflag = 1;
      } else
      if (MATCH(argv[i],"-numscans") && i+1<argc){
        numscans = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-infixes") && i+2<argc) {
        strcpy(real_infix,argv[i+1]); strcpy(imag_infix,argv[i+2]); i+=2;
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
  /* vol included for backward compatibility */
  if (MATCH(datatype,"BRIK") ||
      MATCH(datatype,"bshort") ||
      MATCH(datatype,"bfloat") ||
      MATCH(datatype,"vol")) {
    voldataflag = 1;
  } else if (MATCH(datatype,"wtts") ||
             MATCH(datatype,"wtfs") ||
             MATCH(datatype,"surf")) {
    voldataflag = 0;
  } else {
    ErrorExit(ERROR_BADPARM,"%s: ### datatype %s not supported ...quitting\n",
                progname,datatype);
  }

  if(inputfsflag) {
    if (MATCH(datatype,"vol"))
      sprintf(datatype,"bfloat");
    else if (MATCH(datatype,"surf"))
      sprintf(datatype,"wtfs");
  } else {
    if (MATCH(datatype,"vol"))
      sprintf(datatype,"BRIK");
    else if (MATCH(datatype,"surf"))
      sprintf(datatype,"wtts");
    else if (MATCH(datatype,"bfloat") || MATCH(datatype,"wtfs"))
      inputfsflag = 1;
  }

  if (MATCH(instem,UNDEFSTR))
    ErrorExit(ERROR_BADPARM,"%s: ### instem not specified ...quitting\n",
              progname);
  if (MATCH(outstem,UNDEFSTR))
    sprintf(outstem, "%s-fourier", instem);
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
  if(inputfsflag && outputfsflag)
    ErrorExit(ERROR_BADPARM,"%s: ### fourier series output not allowed when input is also fourier series ...quitting\n",
              progname,outdir);
  if (!voldataflag) {
    if (MATCH(subj,UNDEFSTR))
      ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
                progname);
    if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh"))
      ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",
                progname,hemi);
  }
  if(trendorder < 0 || trendorder > 2)
    ErrorExit(ERROR_BADPARM,"%s: ### trend order must be 0 (none), 1 (linear), or 2 (quadratic) ...quitting\n",
              progname);
/* todo: actually implement quadratic detrend */
  if(numscans < 1)
    ErrorExit(ERROR_BADPARM,"%s: ### numscans must be >= 1 ...quitting\n",
              progname);

  /* check that input files exist */
  if (inputfsflag) {
    if (voldataflag) {
      nslices = getNumBfloats(indir,instem,real_infix);
      if(nslices<=0)
        ErrorExit(ERROR_BADPARM,"%s: ### no bfloats found in %s/ with stem %s%s ...quitting\n",
               progname,indir,instem,real_infix);
      i = getNumBfloats(indir,instem,imag_infix);
      if(i<=0)
        ErrorExit(ERROR_BADPARM,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
               progname,indir,instem,imag_infix);
      if(i!=nslices)
        ErrorExit(ERROR_BADPARM,"%s: ### number of real and imaginary bfloats (%d and %d) do not match ...quitting\n",
               progname,nslices,i);
    } else { /* surface data (wt file) */
      sprintf(tempstr,"%s/%s%s-%s.wt",indir,instem,real_infix,hemi);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
      sprintf(tempstr,"%s/%s%s-%s.wt",indir,instem,imag_infix,hemi);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
    }
  } else {
    if (voldataflag) {
      if (MATCH(datatype,"BRIK")) {
        sprintf(tempstr,"%s/%s.BRIK",indir,instem);
        if (!FileExists(tempstr))
          ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                    progname,tempstr);
        sprintf(tempstr,"%s/%s.HEAD",indir,instem);
        if (!FileExists(tempstr))
          ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                    progname,tempstr);
      } else if (MATCH(datatype,"bshort")) {
        nslices = getNumBshorts(indir,instem);
        if(nslices<=0)
          ErrorExit(ERROR_BADPARM,"%s: ### no bshorts found in %s/ with stem %s ...quitting\n",
                 progname,indir,instem);
      } else {
        ErrorExit(ERROR_BADPARM,"%s: wrong datatype %s for volume data timeseries ...quitting\n",
                progname,datatype);
      }
    } else {
      sprintf(tempstr,"%s/%s-%s.wt",indir,instem,hemi);
      if (!FileExists(tempstr))
        ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                  progname,tempstr);
    }
  }

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

/* read x,y,t slice for a single slice */
void read_BRIK_xyt(float ***xyt, char *indir, char *instem, int slicenum,
                   float *ffacts, int *btypes, char *byteorder)
{
  FILE *fp;
  char tempstr[STRLEN];
  long offset=0;
  int t,y,x;
  float ffact;
  int btype;
  float max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0,avg,stdev,f;
  int npix=0;
  short s;
  
  if(slicenum < 0 || slicenum > nslices)
    ErrorExit(ERROR_BADPARM,"%s: ### trying to read brik with bad slice index ...quitting\n",
            progname,tempstr);

  sprintf(tempstr,"%s/%s.BRIK",indir,instem);
  MsgPrintf("%s: reading slice %d of AFNI BRIK %s\n",
            progname,slicenum,tempstr);
  fp = fopen(tempstr,"rb");
  if (fp==NULL)
    ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
            progname,tempstr);

  for (t=0;t<tpoints;t++) {
    ffact = ffacts[t];
    if (ffact==0) ffact = 1;
    btype = btypes[t];
    switch (btype) {
      case BRIK_SHORT:
        offset = (t*nslices+slicenum)*ynum*xnum*sizeof(short);
        break;
      case BRIK_FLOAT:
        offset = (t*nslices+slicenum)*ynum*xnum*sizeof(float);
        break;
      default:
        ErrorExit(ERROR_BADFILE,"%s: ### btype %d is currently unsupported ...quitting\n",
                progname,btype);
    }  

    fseek(fp,offset,SEEK_SET);  /* TODO: check file size */
    max = -BIGFLOAT;
    min = BIGFLOAT;
    npix = 0;
    sum = sum2 = 0;
    for (y=0;y<ynum;y++)
    for (x=0;x<xnum;x++) {
      if (btype==BRIK_SHORT) {
        s = freadShort(fp);
        if(MATCH(byteorder,"LSB")) s = swapShort(s);
        f = xyt[t][y][x] = ffact*(float)s;
      } else if (btype==BRIK_FLOAT) {
        f = freadFloat(fp);
        if(MATCH(byteorder,"LSB")) f = swapFloat(f);
        f = xyt[t][y][x] = ffact*f;
      } else
        f = 0;
      if (f>max) max = f;
      if (f<min) min = f;
      sum += f;
      sum2 += f*f;
      npix++;
    }
  }
  avg = sum/npix;
  stdev = sqrt(sum2/npix-avg*avg);
  
  MsgPrintf("%s: last time step = %d:\n",progname,t);
  MsgPrintf("%s: pixels(x*y)=%d, avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            progname,npix,avg,stdev,min,max);
  fclose(fp);
}

/* todo: remove quadratic trend */
void remove_trend(float *v,int n)
{
  int k;
  double a,b,c1,c2,c11,c12,c21,c22,f;

  c1 = c2 = c11 = c12 = c21 = c22 = 0;
  for (k=0;k<n;k++) {
    f = v[k];
    c1 += f;
    c2 += k*f;
    c11 += k;
    c12 += 1;
    c21 += k*k;
    c22 += k;
  }
  a = (c1*c22-c2*c12)/(c11*c22-c12*c21);
  b = (c2*c11-c1*c21)/(c11*c22-c12*c21);
  for (k=0;k<n;k++)
    v[k] -= a*k+b;
}

void fourier(float *data,int nn,int isign)
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  n=nn << 1; /* multiply by 2 -- numerical recipes bullshit presumably */
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}

void ftrans(float *data,float *fdata,int nsamp)
{
  float tmpvec[MAXIM]; /* todo: dynamically allocate this array */
  int j;

  for (j=0;j<nsamp;j++) {
    tmpvec[2*j] = data[j];
    tmpvec[2*j+1] = 0;
  }
  fourier(tmpvec-1,nsamp,1); /* tmpvec-1 is presumably passed because fourier uses stupid 1-based indexing */
  for (j=0;j<nsamp*2;j++)
    fdata[j] = tmpvec[j];
}

void compute_fourier(float ***xyt, float ***xyf_r, float ***xyf_i)
{
  int x,y,t;
  
  MsgPrintf("%s: running fft's...\n",progname);

  for (y=0;y<ynum;y++)
  for (x=0;x<xnum;x++) {
    for (t=0;t<tpoints;t++) ts[t] = xyt[t][y][x];
    if (trendorder) remove_trend(ts,tpoints);
    ftrans(ts,cfs,tpoints);
    for (t=0;t<numfreqs;t++) {
      xyf_r[t][y][x] = cfs[2*t];
      xyf_i[t][y][x] = cfs[2*t+1];
    }
  }
}

void write_fourier(float ***xyf_r, float ***xyf_i, int slicenum, char *outdir, char *outstem)
{
  char tempstr[STRLEN];

  sprintf(tempstr,"%s/%s%s_%03d.hdr",outdir,outstem,real_infix,slicenum);
  writeFloatHeader(tempstr,xnum,ynum,numfreqs,0);
  sprintf(tempstr,"%s/%s%s_%03d.bfloat",outdir,outstem,real_infix,slicenum);
  writeFloatImageTS(tempstr,fseries_r,xnum,ynum,numfreqs);

  sprintf(tempstr,"%s/%s%s_%03d.hdr",outdir,outstem,imag_infix,slicenum);
  writeFloatHeader(tempstr,xnum,ynum,numfreqs,0);
  sprintf(tempstr,"%s/%s%s_%03d.bfloat",outdir,outstem,imag_infix,slicenum);
  writeFloatImageTS(tempstr,fseries_i,xnum,ynum,numfreqs);
}

/*====================
 Begin from numrec_c
=====================*/
#define ITMAX 100
#define EPS 3.0e-7
float betacf(float a, float b, float x)
{
  float qap,qam,qab,em,tem,d;
  float bz,bm=1.0,bp,bpp;
  float az=1.0,am=1.0,ap,app,aold;
  int m;

  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  bz=1.0-qab*x/qap;
  for (m=1;m<=ITMAX;m++) {
    em=(float) m;
    tem=em+em;
    d=em*(b-em)*x/((qam+tem)*(a+tem));
    ap=az+d*am;
    bp=bz+d*bm;
    d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem));
    app=ap+d*az;
    bpp=bp+d*bz;
    aold=az;
    am=ap/bpp;
    bm=bp/bpp;
    az=app/bpp;
    bz=1.0;
    if (fabs(az-aold) < (EPS*fabs(az))) return az;
  }
  MsgPrintf("%s: bad input to betacf: a=%2.2f, b=%2.2f, x=%2.2f\n",
            progname,a,b,x);
  return 0;
}
#undef ITMAX
#undef EPS

float gammln(float xx)
{
  double x,tmp,ser,res;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
    -1.231739516,0.120858003e-2,-0.536382e-5};
  int j;

  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  res = -tmp+log(2.50662827465*ser);
  return res;
}

float betai(float a, float b, float x)
{
  float bt;

  if (x < 0.0 || x > 1.0) 
    MsgPrintf("%s: betai error: bad x\n",progname);
  if (x == 0.0 || x == 1.0) bt=0.0;
  else
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}

/*====================
 End from numrec_c
=====================*/

float sigf(float F,int df1,int df2)
{
  float sig;

  sig = betai(0.5*df2,0.5*df1,df2/(df2+df1*F));
  return sig;
}

void compute_fsig(float *fsig_real, float *fsig_imag, float *fs_r, float *fs_i)
{
  int k;
  int n1=0,n2=0,dof1,dof2;
  double sum1=0,sum2=0,numer,denom;
  float F,fsig;
  float real,imag;

  if(rawflag) {
    *fsig_real = fs_r[stimfreq];
    *fsig_imag = fs_i[stimfreq];
    return;
  }

  for (k=0;k<numfreqs;k++) {
    real = fs_r[k];
    imag = fs_i[k];
    amp[k]=hypot(real,imag);
    phase[k]=atan2(imag,real);
  }
  sum1 = amp[stimfreq]*amp[stimfreq];
  n1 = 1;
  sum2 = n2 = 0;
  for (k=lofreq;k<=hifreq;k++) {
    if (((k<fstimfreq-1)||(k>fstimfreq+1))&&
        ((k<2*fstimfreq-1)||(k>2*fstimfreq+1))&&
        ((k<3*fstimfreq-1)||(k>3*fstimfreq+1))&&
         (k!=omitfreq)) {
      sum2 += amp[k]*amp[k];
      n2++;
    }
  }
  dof1 = n1*2; /* why * 2? */
  dof2 = n2*2; /* is it because amplitude comes from two components? */
  if(sum2==0) sum2=1;
  numer = sum1/dof1;
  denom = sum2/dof2;
  F = numer/denom;
  fsig = sigf(F,dof1*numscans,dof2*numscans);
  if (fsig>1) fsig=1;
  if (fsig<=0) fsig=1e-40;
  fsig = -log10(fsig);
  if(fratioflag) {
    *fsig_real = sqrt(F)*cos(phase[stimfreq]);
    *fsig_imag = sqrt(F)*sin(phase[stimfreq]);
  } else {
    *fsig_real = fsig*cos(phase[stimfreq]);
    *fsig_imag = fsig*sin(phase[stimfreq]);
  }
}

void read_cx_bfloat_ts(float ***xyt_r, float ***xyt_i, int slicenum, char *indir, char *instem)
{
  char tempstr[STRLEN];

  sprintf(tempstr,"%s/%s%s_%03d.bfloat",indir,instem,real_infix,slicenum);
  readFloatImageTS(tempstr,xyt_r,xnum,ynum,numfreqs);
  sprintf(tempstr,"%s/%s%s_%03d.bfloat",indir,instem,imag_infix,slicenum);
  readFloatImageTS(tempstr,xyt_i,xnum,ynum,numfreqs);
}

void write_cx_bfloat(float **xy_r, float **xy_i, int slicenum, char *outdir, char *outstem)
{
  char tempstr[STRLEN];

  sprintf(tempstr,"%s/%s%s_%03d.hdr",outdir,outstem,real_infix,slicenum);
  writeFloatHeader(tempstr,xnum,ynum,1,0);
  sprintf(tempstr,"%s/%s%s_%03d.bfloat",outdir,outstem,real_infix,slicenum);
  writeFloatImage(tempstr,xy_r,xnum,ynum);

  sprintf(tempstr,"%s/%s%s_%03d.hdr",outdir,outstem,imag_infix,slicenum);
  writeFloatHeader(tempstr,xnum,ynum,1,0);
  sprintf(tempstr,"%s/%s%s_%03d.bfloat",outdir,outstem,imag_infix,slicenum);
  writeFloatImage(tempstr,xy_i,xnum,ynum);
}

int main(int argc, char **argv)
{
  int t,x,y,z,i,k;
  char tempstr[STRLEN];
  int xnum_tmp,ynum_tmp,numfreqs_tmp,numsurfvals,numsurfvals_tmp,nverts;
  FILE *fp_ts=NULL,*fp_fs_r=NULL,*fp_fs_i=NULL,*fp_fsig_r=NULL,*fp_fsig_i=NULL;
  int vnum, vnum_tmp;
  float fsig_real,fsig_imag,real,imag;
  float *scalefacts=NULL;
  int *briktypes=NULL;
  char byteorder[STRLEN]=UNDEFSTR;
  
  parse_args(argc,argv);

  if (voldataflag) {
    /* get dimensions of dataset */
    if (inputfsflag) {
      sprintf(tempstr,"%s/%s%s_000.hdr",indir,instem,real_infix);
      MsgPrintf("%s: reading real bfloat header\n  %s\n",progname,tempstr);
      readFloatHeader(tempstr,&xnum,&ynum,&numfreqs);
      /* make sure dimensions match for imaginary */
      sprintf(tempstr,"%s/%s%s_000.hdr",indir,instem,imag_infix);
      MsgPrintf("%s: reading imag bfloat header\n  %s\n",progname,tempstr);
      readFloatHeader(tempstr,&xnum_tmp,&ynum_tmp,&numfreqs_tmp);
      if(xnum_tmp!=xnum || ynum_tmp!=ynum || numfreqs_tmp!=numfreqs)
        ErrorExit(ERROR_BADFILE,
        "%s: ### image dimensions for imaginary do not match real...quitting\n",
        progname);
      if (numfreqs<=1) {
        ErrorExit(ERROR_BADFILE,
        "%s: ### bfloat images have depth <= 1: not a fourier series? ...quitting\n",
        progname);
      }
    } else {
      if (MATCH(datatype,"BRIK")) {
        sprintf(tempstr,"%s/%s.HEAD",indir,instem);
        MsgPrintf("%s: reading AFNI header %s\n",progname,tempstr);
        read_HEAD(tempstr,&xnum,&ynum,&nslices,&tpoints,
                  &scalefacts,&briktypes,byteorder);
      } else if (MATCH(datatype,"bshort")) {
        sprintf(tempstr,"%s/%s_000.hdr",indir,instem);
        MsgPrintf("%s: reading bshort header\n  %s\n",progname,tempstr);
        readBshortHeader(tempstr,&xnum,&ynum,&tpoints);
      }
      /* todo: should use dynamic allocation instead */
      if (tpoints>MAXIM/2)
        ErrorExit(ERROR_BADFILE,"%s: ### reps > %d ...quitting\n",
                  progname, MAXIM/2);
      numfreqs = tpoints/2;
    }
    if(hifreq<lofreq || hifreq>=numfreqs) hifreq = numfreqs-1;

    /* allocate memory for data */
    MsgPrintf("%s: starting to allocate memory\n",progname);
    if (!inputfsflag) {
      tseries = (float ***)calloc(tpoints,sizeof(float **));         MTEST(tseries);
      for (t=0;t<tpoints;t++) {
        tseries[t] = (float **)calloc(ynum,sizeof(float *));         MTEST(*tseries);
        for (y=0;y<ynum;y++) {
          tseries[t][y] = (float *)calloc(xnum,sizeof(float));       MTEST(**tseries);
        }
      }
    }
    fseries_r = (float ***)calloc(numfreqs,sizeof(float **));      MTEST(fseries_r);
    fseries_i = (float ***)calloc(numfreqs,sizeof(float **));      MTEST(fseries_i);
    for (t=0;t<numfreqs;t++) {
      fseries_r[t] = (float **)calloc(ynum,sizeof(float *));       MTEST(*fseries_r);
      fseries_i[t] = (float **)calloc(ynum,sizeof(float *));       MTEST(*fseries_i);
      for (y=0;y<ynum;y++) {
        fseries_r[t][y] = (float *)calloc(xnum,sizeof(float));     MTEST(**fseries_r);
        fseries_i[t][y] = (float *)calloc(xnum,sizeof(float));     MTEST(**fseries_i);
      }
    }
    fsig_r = (float **)calloc(ynum,sizeof(float *));               MTEST(fsig_r);
    fsig_i = (float **)calloc(ynum,sizeof(float *));               MTEST(fsig_i);
    for (y=0;y<ynum;y++) {
      fsig_r[y] = (float *)calloc(xnum,sizeof(float));             MTEST(*fsig_r);
      fsig_i[y] = (float *)calloc(xnum,sizeof(float));             MTEST(*fsig_i);
    }
    MsgPrintf("%s: finished allocating memory\n",progname);

    for (z=0;z<nslices;z++) {
      if (inputfsflag) {
        read_cx_bfloat_ts(fseries_r,fseries_i,z,indir,instem);
      } else {
        if (MATCH(datatype,"BRIK")) {
          read_BRIK_xyt(tseries,indir,instem,z,
                        scalefacts,briktypes,byteorder);
        } else if (MATCH(datatype,"bshort")) {
          sprintf(tempstr,"%s/%s_%03d.bshort",indir,instem,z);
          readBshortImageTS(tempstr,tseries,xnum,ynum,tpoints);
        }
        compute_fourier(tseries,fseries_r,fseries_i);
        if (outputfsflag)
          write_fourier(fseries_r,fseries_i,z,outdir,outfsstem);    
      }
      /* compute fsig for each pixel in slice */
      for (y=0;y<ynum;y++)
      for (x=0;x<xnum;x++) {
        for (k=0;k<numfreqs;k++) {
          fs_r[k]=fseries_r[k][y][x];
          fs_i[k]=fseries_i[k][y][x];
        }
        compute_fsig(&(fsig_r[y][x]),&(fsig_i[y][x]),fs_r,fs_i);
      }
      write_cx_bfloat(fsig_r,fsig_i,z,outdir,outstem);
    }
  } else {
    nverts = nverticesSurf(subj,hemi);
    if (inputfsflag) {
      /* open fourier series surface data file */
      sprintf(tempstr,"%s/%s%s-%s.wt",indir,instem,real_infix,hemi);
      MsgPrintf("%s: opening fourier series file %s\n",progname,tempstr);
      fp_fs_r = fopen(tempstr,"r");
      if (fp_fs_r==NULL)
        ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,tempstr);
      fread2(&numfreqs,fp_fs_r);
      MsgPrintf("%s: number of surface values: %d\n",progname,numsurfvals);

      /* open imaginary data file -- make sure dimensions match */
      sprintf(tempstr,"%s/%s%s-%s.wt",indir,instem,imag_infix,hemi);
      MsgPrintf("%s: opening fourier series file %s\n",progname,tempstr);
      fp_fs_i = fopen(tempstr,"r");
      if (fp_fs_i==NULL)
        ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,tempstr);
      fread2(&numfreqs_tmp,fp_fs_i);
      if(numfreqs_tmp!=numfreqs) {
        MsgPrintf("%s: number of imag surface frequencies: %d\n",progname,numfreqs_tmp);
        ErrorExit(ERROR_BADFILE,
        "%s: ### image dimensions for imaginary do not match real...quitting\n",
        progname);
      }
    } else {
      /* open time series surface data file */
      sprintf(tempstr,"%s/%s-%s.wt",indir,instem,hemi);
      MsgPrintf("%s: opening time series file %s\n",progname,tempstr);
      fp_ts = fopen(tempstr,"r");
      if (fp_ts==NULL)
        ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,tempstr);
      fread2(&tpoints,fp_ts);
      MsgPrintf("%s: number of surface time points: %d\n",progname,tpoints);
      numfreqs = tpoints/2;
    }
    if (2*numfreqs>MAXIM/2)
      ErrorExit(ERROR_BADFILE,"%s: ### reps > %d ...quitting\n",
                progname, MAXIM/2);
    if(hifreq<lofreq || hifreq>=numfreqs) hifreq = numfreqs-1;

    /* allocate memory for data */
    MsgPrintf("%s: starting to allocate memory\n",progname);
    surf_vnums = (int *)calloc(nverts,sizeof(int));              MTEST(surf_vnums);
    if (!inputfsflag) {
      surf_tseries = (float **)calloc(tpoints,sizeof(float *));  MTEST(surf_tseries);
      for (t=0;t<tpoints;t++) {
        surf_tseries[t] = (float *)calloc(nverts,sizeof(float)); MTEST(*surf_tseries);
      }
    }
    surf_fseries_r = (float **)calloc(numfreqs,sizeof(float *)); MTEST(surf_fseries_r);
    surf_fseries_i = (float **)calloc(numfreqs,sizeof(float *)); MTEST(surf_fseries_i);
    for (k=0;k<numfreqs;k++) {
      surf_fseries_r[k] = (float *)calloc(nverts,sizeof(float)); MTEST(*surf_fseries_r);
      surf_fseries_i[k] = (float *)calloc(nverts,sizeof(float)); MTEST(*surf_fseries_i);
    }
    MsgPrintf("%s: finished allocating memory\n",progname);

    if (inputfsflag) {
      /* read fourier series files */
      MsgPrintf("%s: reading fourier series files\n",progname);
      for (k=0;k<numfreqs;k++) {
        fread3(&numsurfvals,fp_fs_r);
        MsgPrintf("%s: number of surface values for freq %d: %d\n",
                  progname,k,numsurfvals);
        fread3(&numsurfvals_tmp,fp_fs_i);
        if(numsurfvals_tmp!=numsurfvals) {
          ErrorExit(ERROR_BADFILE,
          "%s: ### image dimensions for imaginary do not match real...quitting\n",
          progname);
        }
        for (i=0;i<numsurfvals;i++) {
          fread3(&vnum,fp_fs_r);
          fread3(&vnum_tmp,fp_fs_i);
          if (vnum!=vnum_tmp)
            ErrorExit(ERROR_BADFILE,
            "%s: ### vertex number for imaginary does not match real...quitting\n",
            progname);
          fread4(&real,fp_fs_r);
          fread4(&imag,fp_fs_i);
          surf_fseries_r[k][i]=real;
          surf_fseries_i[k][i]=imag;
          surf_vnums[i] = vnum;
        }
      }
    } else {
      /* read time series file */
      MsgPrintf("%s: reading time series file\n",progname);
      for (t=0;t<tpoints;t++) {
        fread3(&numsurfvals,fp_ts);
        MsgPrintf("%s: number of surface values for time point %d: %d\n",
                  progname,t,numsurfvals);
        for (i=0;i<numsurfvals;i++) {
          fread3(&vnum,fp_ts);
          fread4(&real,fp_ts);
          surf_tseries[t][i]=real;
          surf_vnums[i] = vnum;
        }
      }

      /* compute fourier series foreach vertex */
      for (i=0;i<numsurfvals;i++) {
        for (t=0;t<tpoints;t++) {
          ts[t] = surf_tseries[t][i];
        }        

        if (trendorder) remove_trend(ts,tpoints);
        ftrans(ts,cfs,tpoints);
        for (k=0;k<numfreqs;k++) {
          surf_fseries_r[k][i] = cfs[2*k];
          surf_fseries_i[k][i] = cfs[2*k+1];
        }
      }

      /* write real and imag to output file */
      if (outputfsflag) {
        sprintf(tempstr,"%s/%s%s-%s.wt",outdir,outfsstem,real_infix,hemi);
        MsgPrintf("%s: initializing fourier series output file %s...\n", progname,tempstr);
        fp_fs_r = fopen(tempstr,"w");
        if (fp_fs_r==NULL)
          ErrorExit(ERROR_BADFILE,"%s: ### cannot create file %s\n",progname,tempstr);

        sprintf(tempstr,"%s/%s%s-%s.wt",outdir,outfsstem,imag_infix,hemi);
        MsgPrintf("%s: initializing fourier series output file %s...\n", progname,tempstr);
        fp_fs_i = fopen(tempstr,"w");
        if (fp_fs_i==NULL)
          ErrorExit(ERROR_BADFILE,"%s: ### cannot create file %s\n",progname,tempstr);

        MsgPrintf("%s: number of surface values: %d\n",progname,numsurfvals);
        MsgPrintf("%s: number of surface frequency points: %d\n",progname,numfreqs);
        fwrite2(numfreqs,fp_fs_r);
        fwrite2(numfreqs,fp_fs_i);
        for (k=0;k<numfreqs;k++) {
          for (i=0,numsurfvals=0;i<nverts;i++)
            if (surf_fseries_r[k][i]!=0 || surf_fseries_i[k][i]!=0) numsurfvals++;
          fwrite3(numsurfvals,fp_fs_r);
          fwrite3(numsurfvals,fp_fs_i);
          for (i=0;i<nverts;i++) {
            real = surf_fseries_r[k][i];
            imag = surf_fseries_i[k][i];
            if (real!=0 || imag!=0) {
              vnum = surf_vnums[i];
              fwrite3(vnum,fp_fs_r);
              fwrite3(vnum,fp_fs_i);
              fwrite4(real,fp_fs_r);
              fwrite4(imag,fp_fs_i);
            }
          }
        }
      }
    }

    /* open output surface data files */
    sprintf(tempstr,"%s/%s%s-%s.w",outdir,outstem,real_infix,hemi);
    MsgPrintf("%s: creating output file %s\n",progname,tempstr);
    fp_fsig_r = fopen(tempstr,"w");
    if (fp_fsig_r==NULL)
      ErrorExit(ERROR_BADFILE,"%s: ### unable to create file %s\n",progname,tempstr);

    sprintf(tempstr,"%s/%s%s-%s.w",outdir,outstem,imag_infix,hemi);
    MsgPrintf("%s: creating output file %s\n",progname,tempstr);
    fp_fsig_i = fopen(tempstr,"w");
    if (fp_fsig_i==NULL)
      ErrorExit(ERROR_BADFILE,"%s: ### unable to create file %s\n",progname,tempstr);

    /* write "header" info for fsig output file */
    fwrite2(0,fp_fsig_r);
    fwrite2(0,fp_fsig_i);
    fwrite3(nverts,fp_fsig_r);
    fwrite3(nverts,fp_fsig_i);

    /* compute fsig for each vertex and write to file */
    for (i=0;i<nverts;i++) {
      vnum = surf_vnums[i];
      for (k=0;k<numfreqs;k++) {
        fs_r[k]=surf_fseries_r[k][i];
        fs_i[k]=surf_fseries_i[k][i];
      }
      compute_fsig(&fsig_real,&fsig_imag,fs_r,fs_i);
      
      /* write to output file */
      fwrite3(vnum,fp_fsig_r);
      fwrite3(vnum,fp_fsig_i);
      fwrite4(fsig_real,fp_fsig_r);
      fwrite4(fsig_imag,fp_fsig_i);
    }

    if (fp_ts) fclose(fp_ts);
    if (fp_fs_r) fclose(fp_fs_r);
    if (fp_fs_i) fclose(fp_fs_i);
    fclose(fp_fsig_r);
    fclose(fp_fsig_i);
  }

  MsgPrintf("%s: finished.\n",progname,tempstr);
  exit(0);
}

