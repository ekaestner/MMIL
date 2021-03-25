#include "surflib.h"

/* sections in this file: */
/* file io */
/* error */
/* random */
/* bfloats */
/* matrix */
/* surface */
/* transform */
/* numerical recipes */

/*##########################################################################*/
/*---------------------------------------------------------------------------*/
/* file io */
/*---------------------------------------------------------------------------*/
short
swapShort(short s)
{
  SWAP_SHORT ss;
  char       c;

  /* first swap bytes in word */
  ss.s = s;
  c = ss.buf[0];
  ss.buf[0] = ss.buf[1] ;
  ss.buf[1] = c;

  return(ss.s);
}

/*---------------------------------------------------------------------------*/
int
swapInt(int i)
{
  SWAP_LONG  sl ;
  short      s ;

  sl.i = i ;  /* now swap words */
  sl.s[0] = swapShort(sl.s[0]) ;
  sl.s[1] = swapShort(sl.s[1]) ;
  s = sl.s[0] ;  /* now swap words */
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;

  return(sl.i) ;
}

/*---------------------------------------------------------------------------*/
unsigned int
swapUInt(unsigned int ui)
{
  SWAP_LONG  sl ;
  short      s ;

  sl.ui = ui ;  /* now swap words */
  sl.s[0] = swapShort(sl.s[0]) ;
  sl.s[1] = swapShort(sl.s[1]) ;
  s = sl.s[0] ;  /* now swap words */
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;

  return(sl.ui) ;
}

/*---------------------------------------------------------------------------*/
float
swapFloat(float f)
{
  SWAP_LONG  sl;
  short      s;

  sl.f = f;  /* first swap bytes in each word */
  sl.s[0] = swapShort(sl.s[0]);
  sl.s[1] = swapShort(sl.s[1]);
  s = sl.s[0];  /* now swap words */
  sl.s[0] = sl.s[1];
  sl.s[1] = s;

  return(sl.f);
}

/*---------------------------------------------------------------------------*/
int
fread1(int *v, FILE *fp)
{
  unsigned char c;
  int  ret ;

  ret = fread(&c,1,1,fp);
  *v = c;
  return(ret) ;
}

/*---------------------------------------------------------------------------*/
int
fread2(int *v, FILE *fp)
{
  short s;
  int   ret ;

  ret = fread(&s,2,1,fp);
#ifdef Linux
  s = swapShort(s) ;
#endif
  *v = s;
  return(ret) ;
}

/*---------------------------------------------------------------------------*/
short
freadShort(FILE *fp)
{
  int   nread;
  short s;

  nread = fread(&s,sizeof(short),1,fp);
#ifdef Linux
  s = swapShort(s);
#endif
  return(s);
}

/*---------------------------------------------------------------------------*/
int
fread3(int *v, FILE *fp)
{
  unsigned int i = 0;
  int  ret ;

  ret = fread(&i,3,1,fp);
#ifdef Linux
  i = (unsigned int)swapInt(i) ;
#endif
  *v = ((i>>8) & 0xffffff);
  return(ret) ;
}

/*---------------------------------------------------------------------------*/
int
fread4(float *v, FILE *fp)
{
  float f;
  int ret;

  ret = fread(&f,4,1,fp);
#ifdef Linux
  f = swapFloat(f);
#endif
  *v = f;
  return(ret);
}

/*---------------------------------------------------------------------------*/
int
freadInt(FILE *fp)
{
  int  i, nread ;

  nread = fread(&i,sizeof(int),1,fp);
#ifdef Linux
  i = swapInt(i) ;
#endif
  return(i) ;
}

/*---------------------------------------------------------------------------*/
unsigned int
freadUInt(FILE *fp)
{
  unsigned int i;
  int nread ;

  nread = fread(&i,sizeof(int),1,fp);
#ifdef Linux
  i = swapUInt(i) ;
#endif
  return(i) ;
}

/*---------------------------------------------------------------------------*/
float
freadFloat(FILE *fp)
{
  float f;
  int ret;

  ret = fread(&f,4,1,fp);
#ifdef Linux
  f = swapFloat(f);
#endif
  return(f);
}

/*---------------------------------------------------------------------------*/
int
fwrite1(int v, FILE *fp)
{
  unsigned char c = (unsigned char)v;

  return(fwrite(&c,1,1,fp));
}

/*---------------------------------------------------------------------------*/
int
fwrite2(int v, FILE *fp)
{
  short s ;

  if (v > 0x7fff)    /* don't let it overflow */
    v = 0x7fff ;
  else if (v < -0x7fff)
    v = -0x7fff ;
  s = (short)v;
#ifdef Linux
  s = swapShort(s) ;
#endif
  return(fwrite(&s,2,1,fp));
}

/*---------------------------------------------------------------------------*/
int
fwrite3(int v, FILE *fp)
{
  unsigned int i = (unsigned int)(v<<8);

#ifdef Linux
  i = (unsigned int)swapInt(i) ;
#endif
  return(fwrite(&i,3,1,fp));
}

/*---------------------------------------------------------------------------*/
int
fwrite4(float v, FILE *fp)
{
#ifdef Linux
  v = swapFloat(v);
#endif
  return(fwrite(&v,4,1,fp));
}

/*---------------------------------------------------------------------------*/
int
fwriteInt(int v, FILE *fp)
{
#ifdef Linux
  v = swapInt(v) ;
#endif
  return(fwrite(&v,sizeof(int),1,fp));
}

/*---------------------------------------------------------------------------*/
int
fwriteFloat(float f, FILE *fp)
{
#ifdef Linux
  f = swapFloat(f) ;
#endif
  return(fwrite(&f,sizeof(float),1,fp));
}

/*---------------------------------------------------------------------------*/
char *
fgetl(char *s, int n, FILE *fp)
{
  char *cp, *cp2 ;
  int  len ;

  do {
    cp = fgets(s, n, fp) ;
    if (!cp)
      return(NULL) ;

    while (isspace((int)*cp))
      cp++ ;

  } while (((*cp) == '#') || ((*cp) == '\n') || ((*cp) == 0)) ;

  for (cp2 = cp ; *cp2 ; cp2++)
    if (*cp2 == '#')
      *cp2 = 0 ;

  len = strlen(cp) ;
  if (cp[len-1] == '\n')  /* strip newline */
    cp[len-1] = 0 ;
  return(cp) ;
}

/*---------------------------------------------------------------------------*/
int
FileExists(char *fname)
{
  FILE *fp ;

  fp = fopen(fname, "r") ;
  if (fp) fclose(fp) ;
  return(fp != NULL) ;
}

/*---------------------------------------------------------------------------*/
char *
FileNamePath(char *fname, char *pathName)
{
  char *slash ;

  strcpy(pathName, fname) ;
  slash = strrchr(pathName, '/') ;
  if (slash)
    *slash = 0 ;          /* remove file name */
  else
#ifndef Linux
    getwd(pathName)  ;    /* no path at all, must be cwd */
#else
    sprintf(pathName, ".") ;
#endif
  return(pathName) ;
}

/*---------------------------------------------------------------------------*/
int isadir(char *path)
{
  struct stat stbuf;

  stat(path,&stbuf);
  if ((stbuf.st_mode & S_IFMT) == S_IFDIR) return 1; else return 0;
}

/* ##########################################################################*/
/*---------------------------------------------------------------------------*/
/* error */
/*---------------------------------------------------------------------------*/
unsigned long  Gdiag      = 0 ;
int            Gdiag_no   = -1 ;
int Gerror = NO_ERROR ;
static int (*error_vfprintf)(FILE *fp,const char *fmt,va_list args) = vfprintf;
#if 0
static void (*error_exit)(int ecode) = (void *)(int)exit ;
#endif
static int quiet = 0;

void
ErrorExit(int ecode, char *fmt, ...)
{
  va_list  args ;

  Gerror = ecode ;
  va_start(args, fmt) ;
  vfprintf(stderr, fmt, args) ;
  fprintf(stderr, "\n") ;
#if 0
  (*error_exit)(ecode) ;
#endif
  exit(ecode);
}

/*---------------------------------------------------------------------------*/
int
ErrorPrintf(int ecode, char *fmt, ...)
{
  va_list  args ;
  FILE     *fp ;

  Gerror = ecode ;
  va_start(args, fmt) ;
  (*error_vfprintf)(stderr, fmt, args) ;
  (*error_vfprintf)(stderr, "\n", NULL) ;

  fp = fopen(ERROR_FNAME, "a") ;
  if (fp) {
    (*error_vfprintf)(fp, fmt, args) ;
    (*error_vfprintf)(fp, "\n", NULL) ;
    fclose(fp) ;     /* close file to flush changes */
  }
  return(ecode) ;
}

/*---------------------------------------------------------------------------*/
void
DiagBreak(void)
{
}

/*---------------------------------------------------------------------------*/
void
setQuiet(int qval)
{
  quiet = qval;
}

/*---------------------------------------------------------------------------*/
int
getQuiet(void)
{
  return(quiet);
}

/*---------------------------------------------------------------------------*/
int
MsgPrintf(char *fmt, ...)
{
  va_list  args ;

  if(!quiet) {
    va_start(args, fmt) ;
    vfprintf(stdout, fmt, args) ;
    fflush(stdout);
    return(0);
  } else {
    return(1);
  }
}

/*##########################################################################*/
/*---------------------------------------------------------------------------*/
/* random */
/*---------------------------------------------------------------------------*/
static long idum = 0L ;
static float ran1(long *idum);

int
setRandomSeed(long seed)
{
  idum = seed ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static float
ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[32];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=32+7;j>=0;j--) {
      k=(*idum)/127773;
      *idum=16807*(*idum-k*127773)-2836*k;
      if (*idum < 0) *idum += 2147483647;
      if (j < 32) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/127773;
  *idum=16807*(*idum-k*127773)-2836*k;
  if (*idum < 0) *idum += 2147483647;
  j= (int)(iy/(1+(2147483647-1)/32));
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=((float)(1.0/2147483647))*(float)iy) > ((float)(1.0-1.2e-7)))
    return ((float)(1.0-1.2e-7));
  else return temp;
}

/*---------------------------------------------------------------------------*/
double
randomNumber(double low, double hi)
{
  double val, range ;

  if (low > hi) {
    val = low ;
    low = hi ;
    hi = val ;
  }
  if (idum == 0L) /* change seed from run to run */
    idum = -1L * (long)(abs((int)time(NULL))) ;
  range = hi - low ;
  val = ran1(&idum) * range + low ;
  if ((val < low) || (val > hi))
    ErrorPrintf(ERROR_BADPARM, "randomNumber(%2.1f, %2.1f) - %2.1f\n",
               (float)low, (float)hi, (float)val) ;
   return(val) ;
}

/*---------------------------------------------------------------------------*/
void
randseed()
{
  unsigned int seed;
  struct timeval tv;
  struct timezone tz;

  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;
  gettimeofday(&tv,&tz);
  seed = (unsigned int)(tv.tv_usec);

  srandom(seed);
  srand48(seed);
}

/*---------------------------------------------------------------------------*/
/* bfloats */
/*---------------------------------------------------------------------------*/
int
getNumBfloats(char *indir, char *instem, char *infix)
{
  FILE *fp;
  char fname[STRLEN];
  int n = 0;
 
  while(1) {
    if (infix==NULL) {
      sprintf(fname,"%s/%s_%03d.bfloat",indir,instem,n);
    } else {
      sprintf(fname,"%s/%s%s_%03d.bfloat",indir,instem,infix,n);
    }
    fp = fopen(fname,"r");
    if (fp==NULL) break;
    n++;
    fclose(fp);
  }
  return(n);
}

int
readFloatHeader(char *fname, int *xsize, int *ysize, int *depth)
{
  FILE *fp;
  char funcname[STRLEN]="readFloatHeader";

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  fscanf(fp,"%d %d %d %*d",xsize,ysize,depth);
  fclose(fp);
  MsgPrintf("%s: stat image size: %dx%d with depth=%d\n",
            funcname,*xsize,*ysize,*depth);

  return(0);
}

/*---------------------------------------------------------------------------*/
int
readFloatHeaders(int N, char **input, char *infix, int *xsize, int *ysize, int *depth)
{
  /*
    input should contain a list of pairs of indirs and instems 
    e.g.: indir1 instem1 indir2 instem2 indir3 instem3 ...     

    this function will read the xsize and ysize for each and make sure
    they are all the same
  */


  int n,x,y,d;
  char fname[STRLEN];
  char funcname[STRLEN]="readFloatHeaders";

  for (n=0;n<N;n++) {
    if(infix==NULL)
      sprintf(fname,"%s/%s_%03d.hdr",input[2*n],input[2*n+1],0);
    else
      sprintf(fname,"%s/%s%s_%03d.hdr",input[2*n],input[2*n+1],infix,0);

    if(n==0) {
      if(readFloatHeader(fname,xsize,ysize,depth)==-1) return(-1);
    } else {
      if(readFloatHeader(fname,&x,&y,&d)==-1) return(-1);
      if (*xsize!=x || *ysize!=y) {
        MsgPrintf("%s: ### diff imsize (%s,%s)\n",funcname,input[0],input[2*n]);
        return(-1);
      }
      if (*depth!=d) {
        MsgPrintf("%s: ### non-matching depth (%s,%s)\n",funcname,input[0],input[2*n]);
        return(-1);
      }
    }
  }
  MsgPrintf("%s: stat image size: %dx%d with depth=%d\n",
            funcname,*xsize,*ysize,depth);

  return(0);
}

/*---------------------------------------------------------------------------*/
int
readFloatImage(char *fname, float **fim, int xsize, int ysize)
{
  FILE *fp;
  int x,y,npix=0;
  double max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0;
  float f;
  char funcname[STRLEN]="readFloatImage";

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  max = -BIGFLOAT;
  min = BIGFLOAT;
  sum = sum2 = npix = 0;
  for (y=0;y<ysize;y++)
  for (x=0;x<xsize;x++) {
    f = fim[y][x] = freadFloat(fp);
    sum += f;
    sum2 += f*f;
    if (f>max) max=f;
    if (f<min) min=f;
    npix++;
  }
  fclose(fp);
  sum /= npix;
  sum2 = sqrt(sum2/npix-sum*sum);
  MsgPrintf("%s: file %s read (%dx%d)\n",funcname,fname,xsize,ysize);
  MsgPrintf("%s:   avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            funcname,sum,sum2,min,max);
  return(0);
}

/*---------------------------------------------------------------------------*/
int
readFloatImageTS(char *fname, float ***fimts, int xsize, int ysize, int tsize)
{
  FILE *fp;
  int x,y,t,npix=0;
  double max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0;
  float f;
  char funcname[STRLEN]="readFloatImageTS";

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  for (t=0;t<tsize;t++) {
    max = -BIGFLOAT;
    min = BIGFLOAT;
    sum = sum2 = npix = 0;
    for (y=0;y<ysize;y++)
    for (x=0;x<xsize;x++) {
      f = fimts[t][y][x] = freadFloat(fp);
      if (f>max) max = f;
      if (f<min) min = f;
      sum += f;
      sum2 += f*f;
      npix++;
    }
  }
  fclose(fp);
  sum /= npix;
  sum2 = sqrt(sum2/npix-sum*sum);
  MsgPrintf("%s: file %s read (%dx%d)\n",funcname,fname,xsize,ysize);
  MsgPrintf("%s: last time step = %d:\n",funcname,t);
  MsgPrintf("%s:   avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            funcname,sum,sum2,min,max);
  return(0);
}

/*---------------------------------------------------------------------------*/
int
writeFloatImage(char *fname, float **fim, int xsize, int ysize)
{
  FILE *fp;
  int x,y,npix=0;
  float f;
  double max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0;
  char funcname[STRLEN]="writeFloatImage";

  fp = fopen(fname,"w");
  if(fp==NULL){
    MsgPrintf("%s: ### can't create file %s\n",funcname,fname);
    return(-1);
  }
  for (y=0;y<ysize;y++) {
    for (x=0;x<xsize;x++)
    {
      f = fim[y][x];
      if (f>max) max=f;
      if (f<min) min=f;
      sum += f;
      sum2 += f*f;
      npix++;
      
      fwriteFloat(f, fp);
    }
  }
  fclose(fp);  
  sum /= npix;
  sum2 = sqrt(sum2/npix-sum*sum);
  MsgPrintf("%s: file %s written (%dx%d)\n",funcname,fname,xsize,ysize);
  MsgPrintf("%s:   avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            funcname,sum,sum2,min,max);
  return(0);
}

/*---------------------------------------------------------------------------*/
int
writeFloatImageTS(char *fname, float ***fimts, int xsize, int ysize, int tsize)
{
  FILE *fp;
  int x,y,t,npix=0;
  float f;
  double max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0;
  char funcname[STRLEN]="writeFloatImage";

  fp = fopen(fname,"w");
  if(fp==NULL){
    MsgPrintf("%s: ### can't create file %s\n",funcname,fname);
    return(-1);
  }
  for (t=0;t<tsize;t++) {
    max = -BIGFLOAT;
    min = BIGFLOAT;
    npix = 0;
    sum = sum2 = 0;
    for (y=0;y<ysize;y++) {
      for (x=0;x<xsize;x++)
      {
        f = fimts[t][y][x];
        if (f>max) max = f;
        if (f<min) min = f;
        sum += f;
        sum2 += f*f;
        npix++;

        fwriteFloat(f, fp);
      }
    }
  }
  fclose(fp);
  sum /= npix;
  sum2 = sqrt(sum2/npix-sum*sum);
  MsgPrintf("%s: file %s written (%dx%d)\n",funcname,fname,xsize,ysize);
  MsgPrintf("%s: last time step = %d:\n",funcname,t);
  MsgPrintf("%s:   avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            funcname,sum,sum2,min,max);
  return(0);
}

/*---------------------------------------------------------------------------*/
int
writeFloatHeader(char *fname, int xsize, int ysize, int t, int num)
{
  FILE *fp;
  char funcname[STRLEN]="writeFloatHeader";

  fp = fopen(fname,"w");
  if(fp==NULL){
    MsgPrintf("%s: ### can't create file %s\n",funcname,fname);
    return(-1);
  }
  fprintf(fp,"%d %d %d %d\n",xsize,ysize,t,num);
  fclose(fp);
  MsgPrintf("%s: file %s written\n",funcname,fname);

  return(0);
}

/*---------------------------------------------------------------------------*/
int
writeTextImage(char *fname, float **fim, int xsize, int ysize)
{
  FILE *fp;
  int i,j,num;
  float f;
  double sum=0,sum2=0,max= -1000,min=1000;
  char funcname[STRLEN]="writeTextImage";

  fp = fopen(fname,"w");
  if(fp==NULL){
    MsgPrintf("%s: ### can't create file %s\n",funcname,fname);
    return(-1);
  }
  for (i=0;i<ysize;i++) {
    for (j=0;j<xsize;j++) {
      f = fim[i][j];
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
      fprintf(fp,"%f  ",f);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  num=xsize*ysize;
  sum /= num;
  sum2 = sqrt(sum2/num-sum*sum);
  MsgPrintf("%s: file %s written (%dx%d)\n",funcname,fname,xsize,ysize);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(0);
}

/*---------------------------------------------------------------------------*/
/* bshorts */
/*---------------------------------------------------------------------------*/
int
getNumBshorts(char *indir, char *instem)
{
  FILE *fp;
  char fname[STRLEN];
  int n = 0;
 
  while(1) {
    sprintf(fname,"%s/%s_%03d.bshort",indir,instem,n);
    fp = fopen(fname,"r");
    if (fp==NULL) break;
    n++;
    fclose(fp);
  }
  return(n);
}

/*---------------------------------------------------------------------------*/
int
readBshortHeader(char *fname, int *xsize, int *ysize, int *tsize)
{
  FILE *fp;
  char funcname[STRLEN]="readFloatHeader";

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  fscanf(fp,"%d %d %d %*d",xsize,ysize,tsize);
  fclose(fp);
  MsgPrintf("%s: stat image size: %dx%d with tsize=%d\n",
            funcname,*xsize,*ysize,*tsize);

  return(0);
}

/*---------------------------------------------------------------------------*/
int
readBshortImageTS(char *fname, float ***fimts, int xsize, int ysize, int tsize)
{
  FILE *fp;
  int x,y,t,npix=0;
  double max=-BIGFLOAT,min=BIGFLOAT,sum=0,sum2=0;
  float f;
  short s;
  char funcname[STRLEN]="readBshortImageTS";

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  for (t=0;t<tsize;t++) {
    max = -BIGFLOAT;
    min = BIGFLOAT;
    sum = sum2 = npix = 0;
    for (y=0;y<ysize;y++)
    for (x=0;x<xsize;x++) {
      s = freadShort(fp);
      f = fimts[t][y][x] = (float)s;
      if (f>max) max = f;
      if (f<min) min = f;
      sum += f;
      sum2 += f*f;
      npix++;
    }
  }
  fclose(fp);
  sum /= npix;
  sum2 = sqrt(sum2/npix-sum*sum);
  MsgPrintf("%s: file %s read (%dx%d)\n",funcname,fname,xsize,ysize);
  MsgPrintf("%s: last time step = %d:\n",funcname,t);
  MsgPrintf("%s:   avg=%6.2f, stdev=%6.2f, min=%6.2f, max=%6.2f\n",
            funcname,sum,sum2,min,max);
  return(0);
}

/*---------------------------------------------------------------------------*/
/*##########################################################################*/
/*---------------------------------------------------------------------------*/
/* matrix */
/*---------------------------------------------------------------------------*/
static int EigenSystem(float *data, int n, float *evalues, float *evectors);
static int compare_evalues(const void *l1, const void *l2);

typedef struct
{
  int   eno ;
  float evalue ;
} EIGEN_VALUE, EVALUE ;

MATRIX *
MatrixAlloc(int rows, int cols, int type)
{
  MATRIX *mat ;
  int    row, nelts ;

  mat = (MATRIX *)calloc(1, sizeof(MATRIX)) ;
  if (!mat)
    ErrorExit(ERROR_NOMEMORY,
              "MatrixAlloc(%d, %d, %d): could not allocate mat",
              rows, cols, type) ;
  mat->rows = rows ;
  mat->cols = cols ;
  mat->type = type ;

  nelts = rows*cols ;
  if (type == MATRIX_COMPLEX)
    nelts *= 2 ;

  mat->data = (float *)calloc(nelts+2, sizeof(float)) ;
  if (!mat->data) {
    fprintf(stderr, "MatrixAlloc(%d, %d): allocation failed\n",
      rows, cols) ;
    exit(1) ;
  }
  mat->data += 2 ;
  mat->rptr = (float **)calloc(rows+1, sizeof(float *)) ;
  if (!mat->rptr) {
    free(mat->data) ;
    free(mat) ;
    ErrorExit(ERROR_NOMEMORY, "MatrixAlloc(%d, %d): could not allocate rptr",
              rows, cols) ;
  }
  for (row = 1 ; row <= rows ; row++) {
    switch (type) {
      case MATRIX_REAL:
        mat->rptr[row] = mat->data + (row-1)*cols - 1 ;
        break ;
      case MATRIX_COMPLEX:
        mat->rptr[row] = (float *)(((CPTR)mat->data) +
            (row-1)*cols - 1) ;
        break ;
      default:
        ErrorReturn(NULL,
                    (ERROR_BADPARM, "MatrixAlloc: unknown type %d\n",type)) ;
    }
  }
  return(mat) ;
}

/*---------------------------------------------------------------------------*/
int
MatrixFree(MATRIX **pmat)
{
  MATRIX *mat ;

  mat = *pmat ;
  *pmat = NULL;

  if (!mat)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixFree: NULL POINTER!\n"));
  mat->data -= 2 ;
  free(mat->data) ;
  free(mat->rptr) ;
  free(mat) ;

  return(0) ;
}

/*---------------------------------------------------------------------------*/

MATRIX *
MatrixTranspose(MATRIX *mIn, MATRIX *mOut)
{
  int  row, col, rows, cols ;

  if (!mOut)
  {
    mOut = MatrixAlloc(mIn->cols, mIn->rows, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
      mOut->rptr[col][row] = mIn->rptr[row][col] ;
  }

  return(mOut) ;
}

/*---------------------------------------------------------------------------*/

MATRIX *
MatrixMultiply(MATRIX *m1, MATRIX *m2, MATRIX *m3)
{
  int   col, row, i, rows, cols, m1_cols ;
  float *r3 ;
  register float val, *r1, *r2 ;
  MATRIX   *m_tmp1 = NULL, *m_tmp2 = NULL ;

  if (m1->cols != m2->rows)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                "MatrixMultiply: m1 cols %d does not match m2 rows %d\n",
                m1->cols, m2->rows)) ;

  if (!m3)
  {
    /* twitzel also did something here */
    if((m1->type == MATRIX_COMPLEX) || (m2->type == MATRIX_COMPLEX)) 
      m3 = MatrixAlloc(m1->rows, m2->cols, MATRIX_COMPLEX);
    else
      m3 = MatrixAlloc(m1->rows, m2->cols, m1->type) ;
    if (!m3)
      return(NULL) ;
  }
  else if ((m3->rows != m1->rows) || (m3->cols != m2->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM, 
                 "MatrixMultiply: (%d x %d) * (%d x %d) != (%d x %d)\n",
                 m1->rows, m1->cols, m2->rows, m2->cols, m3->rows, m3->cols)) ;

  if (m3 == m2)
  {
    m_tmp1 = MatrixCopy(m2, NULL) ;
    m2 = m_tmp1 ;
  }
  if (m3 == m1)
  {
    m_tmp2 = MatrixCopy(m1, NULL) ;
    m1 = m_tmp2 ;
  }
  /*  MatrixClear(m3) ;*/
  cols = m3->cols ;
  rows = m3->rows ;
  m1_cols = m1->cols ;

  /* twitzel modified here */
  if((m1->type == MATRIX_REAL) && (m2->type == MATRIX_REAL))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      r3 = &m3->rptr[row][1] ;
      for (col = 1 ; col <= cols ; col++)
      {
        val = 0.0 ;
        r1 = &m1->rptr[row][1] ;
        r2 = &m2->rptr[1][col] ;
        for (i = 1 ; i <= m1_cols ; i++, r2 += cols)
        {
#if 0
          m3->rptr[row][col] +=
            m1->rptr[row][i] * m2->rptr[i][col] ;
#else
          val += *r1++ * *r2 ;
#endif
        }
        *r3++ = val ;
      }
    }
  } else if((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_COMPLEX)) {
 
    for (row = 1 ; row <= rows ; row++)
      {
  for (col = 1 ; col <= cols ; col++)
    {
      for (i = 1 ; i <= m1->cols ; i++)
        {
    float a, b, c, d ;  /* a + ib and c + id */
    
    a = MATRIX_CELT_REAL(m1,row,i) ;
    b = MATRIX_CELT_IMAG(m1,row,i) ;
    c = MATRIX_CELT_REAL(m2,i,col) ;
    d = MATRIX_CELT_IMAG(m2,i,col) ;
    MATRIX_CELT_REAL(m3,row,col) += a*c - b*d ;
    MATRIX_CELT_IMAG(m3,row,col) += a*d + b*c ;
        }
    }
      }
  } else if((m1->type == MATRIX_REAL) && (m2->type == MATRIX_COMPLEX)) {
    for (row = 1 ; row <= rows ; row++)
      {
  for (col = 1 ; col <= cols ; col++)
    {
      for (i = 1 ; i <= m1->cols ; i++)
        {
    float a, c, d ;  /* a + ib and c + id and b=0 here*/
    
    a = *MATRIX_RELT(m1,row,i);
    c = MATRIX_CELT_REAL(m2,i,col);
    d = MATRIX_CELT_IMAG(m2,i,col);
    MATRIX_CELT_REAL(m3,row,col) += a*c;
    MATRIX_CELT_IMAG(m3,row,col) += a*d;
        }
    }
      }
  }
  if (m_tmp1)
    MatrixFree(&m_tmp1) ;
  if (m_tmp2)
    MatrixFree(&m_tmp2) ;
  return(m3) ;
}

/*---------------------------------------------------------------------------*/

/* calcluate the condition # of a matrix using svd */
float
MatrixConditionNumber(MATRIX *m)
{
  float cond ;
  VECTOR  *v_w ;
  MATRIX  *m_U, *m_V ;
  int     row, rows, cols ;
  float   wmax, wmin, wi ;

  cols = m->cols ;
  rows = m->rows ;
  m_U = MatrixCopy(m, NULL) ;
  v_w = RVectorAlloc(cols, MATRIX_REAL) ;
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;

  /* calculate condition # of matrix */
  svdcmp(m_U->rptr, rows, cols, v_w->rptr[1], m_V->rptr) ;
  wmax = 0.0f ;
  wmin = wmax = RVECTOR_ELT(v_w,1) ;
  for (row = 2 ; row <= rows ; row++)
  {
    wi = fabs(RVECTOR_ELT(v_w,row)) ;
    if (wi > wmax)
      wmax = wi ;
    if (wi < wmin)
      wmin = wi ;
  }

  if (FZERO(wmin))
    cond = 1e8 ;   /* something big */
  else
    cond = wmax / wmin ;

  MatrixFree(&m_U) ;
  VectorFree(&v_w) ;
  MatrixFree(&m_V) ;
  return(cond) ;
}

/*---------------------------------------------------------------------------*/

#define TOO_SMALL   1e-4

MATRIX *
MatrixSVDInverse(MATRIX *m, MATRIX *m_inverse)
{
  VECTOR  *v_w ;
  MATRIX  *m_U, *m_V, *m_w, *m_Ut, *m_tmp ;
  int     row, rows, cols ;
  float   wmax, wmin ;

  cols = m->cols ;
  rows = m->rows ;
  m_U = MatrixCopy(m, NULL) ;
  v_w = RVectorAlloc(cols, MATRIX_REAL) ;
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;
  m_w = MatrixAlloc(cols, cols, MATRIX_REAL) ;

  if (svdcmp(m_U->rptr, rows, cols, v_w->rptr[1], m_V->rptr) != NO_ERROR)
  {
    MatrixFree(&m_U) ;
    VectorFree(&v_w) ;
    MatrixFree(&m_V) ;
    MatrixFree(&m_w) ;
    return(NULL) ;
  }

  wmax = 0.0f ;
  for (row = 1 ; row <= rows ; row++)
    if (fabs(RVECTOR_ELT(v_w,row)) > wmax)
      wmax = fabs(RVECTOR_ELT(v_w,row)) ;
  wmin = TOO_SMALL * wmax ;
  for (row = 1 ; row <= rows ; row++)
  {
    if (fabs(RVECTOR_ELT(v_w, row)) < wmin)
      m_w->rptr[row][row] = 0.0f ;
    else
      m_w->rptr[row][row] = 1.0f / RVECTOR_ELT(v_w,row) ;
  }

  m_Ut = MatrixTranspose(m_U, NULL) ;
  m_tmp = MatrixMultiply(m_w, m_Ut, NULL) ;
  m_inverse = MatrixMultiply(m_V, m_tmp, m_inverse) ;

  MatrixFree(&m_U) ;
  VectorFree(&v_w) ;
  MatrixFree(&m_V) ;
  MatrixFree(&m_w) ;

  MatrixFree(&m_Ut) ;
  MatrixFree(&m_tmp) ;
  return(m_inverse) ;
}

/*---------------------------------------------------------------------------*/

MATRIX *
MatrixCopy(MATRIX *mIn, MATRIX *mOut)
{
  int row, rows, cols ;

  if(mIn == NULL)
    return(NULL);

  rows = mIn->rows ;
  cols = mIn->cols ;

  if (!mOut)
    mOut = MatrixAlloc(rows, cols, mIn->type) ;

  for (row = 1 ; row <= rows ; row++)
    memcpy((char *)(mOut->rptr[row]), (char *)mIn->rptr[row], 
           (cols+1)*sizeof(float)) ;

  return(mOut) ;
}

/*---------------------------------------------------------------------------*/
MATRIX *
MatrixEigenSystem(MATRIX *m, float *evalues, MATRIX *m_evectors)
{
  int     col, i, nevalues, row ;
  EVALUE  *eigen_values ;
  MATRIX  *mTmp ;

  nevalues = m->rows ;
  eigen_values = (EVALUE *)calloc((UINT)nevalues, sizeof(EIGEN_VALUE));
  if (!m_evectors)
    m_evectors = MatrixAlloc(m->rows, m->cols, MATRIX_REAL) ;

  mTmp = MatrixAlloc(m->rows, m->cols, MATRIX_REAL) ;
  if (EigenSystem(m->data, m->rows, evalues, mTmp->data) != NO_ERROR)
    return(NULL) ;

/* 
  sort eigenvalues in order of decreasing absolute value. The
  EIGEN_VALUE structure is needed to record final order of eigenvalue so 
  we can also sort the eigenvectors.
*/
  for (i = 0 ; i < nevalues ; i++)
  {
    eigen_values[i].eno = i ;
    eigen_values[i].evalue = evalues[i] ;
  }
  qsort((char *)eigen_values, nevalues, sizeof(EVALUE), compare_evalues) ;
  for (i = 0 ; i < nevalues ; i++)
    evalues[i] = eigen_values[i].evalue ;

  /* now sort the eigenvectors */
  for (col = 0 ; col < mTmp->cols ; col++)
  {
    for (row = 1 ; row <= mTmp->rows ; row++)
      m_evectors->rptr[row][col+1] = mTmp->rptr[row][eigen_values[col].eno+1] ;
  }

  free(eigen_values) ;
  MatrixFree(&mTmp) ;
  return(m_evectors) ;
}

/*---------------------------------------------------------------------------*/

static int
EigenSystem(float *data, int n, float *evalues, float *evectors)
{
  float   *fptr, *e, **mat, *nr_evalues, norm ;
  int     row, col ;
  
  mat = matrix(1, n, 1, n) ;
  e = fvector(1, n) ;
  nr_evalues = fvector(1, n) ;

  /* convert to numerical recipes format */
  for (fptr = data, row = 1 ; row <= n ; row++)
  {                                    
    for (col = 1 ; col <= n ; col++)
      mat[row][col] = *fptr++ ;
  }
 
/*
  note that tred2 modifies mat so that tqli can use it to
  calculate the eigenvectors of the original matrix.
*/
  tred2(mat, n, nr_evalues, e) ;   /* tridiagonalize it */

  /* find eigenvectors and eigenvalues */
  if (tqli(nr_evalues, e, n, mat) != NO_ERROR)
    return(Gerror) ;

  /* convert back to array format. Columns of evectors are eigenvectors */
  for (fptr = evectors, row = 1 ; row <= n ; row++)
  {
    for (col = 1 ; col <= n ; col++)
      *fptr++ = mat[row][col] ;
  }

  for (row = 1 ; row <= n ; row++)
    evalues[row-1] = nr_evalues[row] ;  /* I hate this 1-based stuff */

  /* normalize eigenvectors */
  for (col = 0 ; col < n ; col++)
  {
    fptr = evectors + col ;
    for (norm = 0.0f, row = 0 ; row < n ; row++)
    {
      norm += (float)SQR(*fptr) ;
      fptr += n ;
    }
    fptr = evectors + col ;
    norm = (float)sqrt(norm) ;
    for (row = 0 ; row < n ; row++)
    {
      *fptr /= norm ;
      fptr += n ;
    }
  }

  free_vector(e, 1, n) ;
  free_vector(nr_evalues, 1, n) ;
  free_matrix(mat, 1, n, 1, n) ; 
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/

static int
compare_evalues(const void *l1, const void *l2)
{
  EVALUE *e1, *e2 ;

  e1 = (EVALUE *)l1 ;
  e2 = (EVALUE *)l2 ;
  return(fabs(e1->evalue) < fabs(e2->evalue) ? 1 : -1) ;
}

float
MatrixSVDEigenValues(MATRIX *m, float *evalues)
{
  float cond ;
  VECTOR  *v_w ;
  MATRIX  *m_U, *m_V ;
  int     row, rows, cols, nevalues, i ;
  float   wmax, wmin, wi ;
  EVALUE  *eigen_values ;

  cols = m->cols ;
  rows = m->rows ;
  m_U = MatrixCopy(m, NULL) ;
  v_w = RVectorAlloc(cols, MATRIX_REAL) ;
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;
  nevalues = m->rows ;
  memset(evalues, 0, nevalues*sizeof(evalues[0])) ;
  
  /* calculate condition # of matrix */
  if (svdcmp(m_U->rptr, rows, cols, v_w->rptr[1], m_V->rptr) != NO_ERROR)
    return(Gerror) ;

  eigen_values = (EVALUE *)calloc((UINT)nevalues, sizeof(EIGEN_VALUE));
  for (i = 0 ; i < nevalues ; i++)
  {
    eigen_values[i].eno = i ;
    eigen_values[i].evalue = RVECTOR_ELT(v_w, i+1) ;
  }
  qsort((char *)eigen_values, nevalues, sizeof(EVALUE), compare_evalues) ;
  for (i = 0 ; i < nevalues ; i++)
    evalues[i] = eigen_values[i].evalue ;

  wmax = 0.0f ;
  wmin = wmax = RVECTOR_ELT(v_w,1) ;
  for (row = 2 ; row <= rows ; row++)
  {
    wi = fabs(RVECTOR_ELT(v_w,row)) ;
    if (wi > wmax)
      wmax = wi ;
    if (wi < wmin)
      wmin = wi ;
  }

  if (FZERO(wmin))
    cond = 1e8 ;  /* something big */
  else
    cond = wmax / wmin ;

  free(eigen_values) ;
  MatrixFree(&m_U) ;
  VectorFree(&v_w) ;
  MatrixFree(&m_V) ;
  return(cond) ;
}

/*---------------------------------------------------------------------------*/
double
Vector3Angle(VECTOR *v1, VECTOR *v2)
{
  double  angle, l1, l2, dot, norm, x, y, z ;

  x = V3_X(v1) ; y = V3_Y(v1) ; z = V3_Z(v1) ; l1 = sqrt(x*x+y*y+z*z) ;
  x = V3_X(v2) ; y = V3_Y(v2) ; z = V3_Z(v2) ; l2 = sqrt(x*x+y*y+z*z) ;
  norm = l1*l2 ;
  if (FZERO(norm))
    return(0.0f) ;
  dot = V3_DOT(v1, v2) ;
  if (dot > norm)
    angle = acos(1.0) ;
  else
    angle = acos(dot / norm) ;
  return(angle) ;
}

/*---------------------------------------------------------------------------*/
float
VectorTripleProduct(VECTOR *v1, VECTOR *v2, VECTOR *v3)
{
  float   x1, x2, y1, y2, z1, z2, x3, y3, z3, total ;

  if (v1->rows != 3 && v1->cols != 1)
    ErrorReturn(0.0f,(ERROR_BADPARM, "VectorCrossProduct: must be 3-vectors"));

  x1 = VECTOR_ELT(v1,1) ; y1 = VECTOR_ELT(v1,2) ; z1 = VECTOR_ELT(v1,3) ;
  x2 = VECTOR_ELT(v2,1) ; y2 = VECTOR_ELT(v2,2) ; z2 = VECTOR_ELT(v2,3) ;
  x3 = VECTOR_ELT(v3,1) ; y3 = VECTOR_ELT(v3,2) ; z3 = VECTOR_ELT(v3,3) ;

  total =  x3 * (y1*z2 - z1*y2) ;
  total += y3 * (z1*x2 - x1*z2) ;
  total += z3 * (x1*y2 - y1*x2) ;
  return(total) ;
}

/*###########################################################################*/
/*---------------------------------------------------------------------------*/
/* surface */
/*---------------------------------------------------------------------------*/
/* static = private */
static int mrisFileNameType(char *fname);
static int mrisFindNeighbors(MRI_SURFACE *mris);
static int mrisComputeVertexDistances(MRI_SURFACE *mris);
static int mrisOrientSurface(MRI_SURFACE *mris);
static int mrisComputeSurfaceDimensions(MRI_SURFACE *mris);
static void mrisNormalize(float v[3]);
static int mrisNormalFace(MRIS *mris, int fac,int n,float norm[]);
static float mrisTriangleArea(MRIS *mris, int fac, int n);
static int mrisOrientPlane(MRI_SURFACE *mris);
static int mrisOrientEllipsoid(MRI_SURFACE *mris);
static MRI_SURFACE *mrisReadAsciiFile(char *fname);
static MRI_SURFACE *mrisReadTriangleFile(char *fname, double pct_over);
static MRI_SURFACE *mrisReadGeoFile(char *fname);
static int mrisRemoveVertexLink(MRI_SURFACE *mris, int vno1, int vno2);
static int mrisComputeBoundaryNormals(MRI_SURFACE *mris);
static int mrisComputeTangentPlanes(MRI_SURFACE *mris);


int
MRISbuildFileName(MRI_SURFACE *mris, char *sname, char *fname)
{
  char   path[STRLEN], *slash, *dot ;

  slash = strchr(sname, '/') ;
  if (!slash) {             /* no path - use same one as mris was read from */
    dot = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (dot && (*(dot-1) == 'h') && (*(dot-2) == 'l' || *(dot-2) == 'r'))
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path,
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else strcpy(fname, sname) ;  /* path specified explicitly */
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisFileNameType(char *fname)
{
  int   type ;
  char  *dot, ext[STRLEN] ;

  ext[0] = 0 ;
  dot = strrchr(fname, '@') ;   /* forces the type of the file */
  if (dot) {
    *dot = 0 ;   /* remove it from file name so that fopen will work */
    strcpy(ext, dot+1) ;
  }
  else {    /* no explicit type - check extension */
    dot = strrchr(fname, '.') ;
    if (dot)
      strcpy(ext, dot+1) ;
  }
  StrUpper(ext) ;
  if (!strcmp(ext, "ASC"))                             type = MRIS_ASCII_FILE;
  else if (!strcmp(ext, "GEO"))                        type = MRIS_GEO_FILE ;
  else if (!strcmp(ext, "TRI") || !strcmp(ext, "ICO")) type = MRIS_ICO_FILE ;
  else if (!strcmp(ext, "VTK"))                        type = MRIS_VTK_FILE ;
  else                                                 type = MRIS_BINARY_FILE;

  return(type) ;
}

/*---------------------------------------------------------------------------*/
char *
StrUpper(char *str)
{
  char *cp ;

  for (cp = str ; *cp ; cp++)
    *cp = (char)toupper(*cp) ;
  return(str) ;
}

/*---------------------------------------------------------------------------*/
int
MRISwriteValues(MRI_SURFACE *mris, char *sname)
{
  char funcname[STRLEN]="MRISwriteValues";
  int k,num;                   /* loop counters */
  float f;
  char  fname[STRLEN], *cp ;
  FILE *fp;
  int maxk=-1;
  double sum=0,sum2=0,max= -1000,min=1000;

  MRISbuildFileName(mris, sname, fname) ;
  cp = strrchr(fname, '.') ;
  if (!cp || *(cp+1) != 'w') strcat(fname, ".w") ;
  MsgPrintf("%s: writing surf values to %s\n", funcname,fname);

  fp = fopen(fname,"wb");
  if (fp==NULL) ErrorExit(ERROR_NOFILE, "Can't create file %s\n",fname) ;

  for (k=0,num=0;k<mris->nvertices;k++)
    if (mris->vertices[k].val!=0) num++;
  fwrite2(0,fp);
  fwrite3(num,fp);
  for (k=0;k<mris->nvertices;k++) {
    if (mris->vertices[k].val!=0) {
      if (k>maxk) maxk=k;
      fwrite3(k,fp);
      f = mris->vertices[k].val;
      if (!finite(f))
        ErrorPrintf(ERROR_BADPARM,
                    "%s(%s): val at vertex %d is not finite",
                    funcname,fname, k) ;
      fwriteFloat(f, fp) ;
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);

  if(num>0) {
    sum /= num;
    sum2 = sqrt(sum2/num-sum*sum);
  } else {
    sum = sum2 = 0;
  }
  MsgPrintf("%s: file %s read \n",funcname,fname);
  MsgPrintf("      with %d vals and max vertex index = %d\n",num,maxk);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISwriteImagValues(MRI_SURFACE *mris, char *sname)
{
  char funcname[STRLEN]="MRISwriteImagValues";
  int k,num;                   /* loop counters */
  float f;
  char  fname[STRLEN], *cp ;
  FILE *fp;
  int maxk=-1;
  double sum=0,sum2=0,max= -1000,min=1000;

  MRISbuildFileName(mris, sname, fname) ;
  cp = strrchr(fname, '.') ;
  if (!cp || *(cp+1) != 'w') strcat(fname, ".w") ;
  MsgPrintf("%s: writing surf values to %s\n", funcname,fname);

  fp = fopen(fname,"wb");
  if (fp==NULL) ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",funcname,fname) ;

  for (k=0,num=0;k<mris->nvertices;k++)
    if (mris->vertices[k].imag_val!=0) num++;
  fwrite2(0,fp);
  fwrite3(num,fp);
  for (k=0;k<mris->nvertices;k++) {
    if (mris->vertices[k].imag_val!=0) {
      if (k>maxk) maxk=k;
      fwrite3(k,fp);
      f = mris->vertices[k].imag_val;
      if (!finite(f))
        ErrorPrintf(ERROR_BADPARM,
                    "%s(%s): imag_val at vertex %d is not finite",
                    funcname,fname, k) ;
      fwriteFloat(f, fp) ;
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);

  if(num>0) {
    sum /= num;
    sum2 = sqrt(sum2/num-sum*sum);
  } else {
    sum = sum2 = 0;
  }
  MsgPrintf("%s: file %s written \n",funcname,fname);
  MsgPrintf("      with %d vals and max vertex index = %d\n",num,maxk);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisFindNeighbors(MRI_SURFACE *mris)
{
  int          n0,n1,i,k,m,n, vno, vtotal, ntotal, vtmp[MAX_NEIGHBORS] ;
  FACE         *f;
  VERTEX       *v ;

  if (Gdiag&DIAG_SHOW&&DIAG_VERBOSE_ON) fprintf(stdout,"finding surf neigh") ;

  for (k=0;k<mris->nvertices;k++) {
    if (k == Gdiag_no) DiagBreak() ;
    v = &mris->vertices[k];
    v->vnum = 0;
    for (m=0;m<v->num;m++) {
      n = v->n[m];     /* # of this vertex in the mth face that it is in */
      f = &mris->faces[v->f[m]];  /* ptr to the mth face */
      /* index of vertex we are connected to */
      n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
      n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
      for (i=0;i<v->vnum && vtmp[i]!=f->v[n0];i++);
      if (i==v->vnum)
        vtmp[(int)v->vnum++] = f->v[n0];
      for (i=0;i<v->vnum && vtmp[i]!=f->v[n1];i++);
      if (i==v->vnum)
        vtmp[(int)v->vnum++] = f->v[n1];
    }
    if (mris->vertices[k].v)
      free(mris->vertices[k].v) ;
    mris->vertices[k].v = (int *)calloc(mris->vertices[k].vnum,sizeof(int));
    if (!mris->vertices[k].v)
      ErrorExit(ERROR_NOMEMORY,
                "mrisFindNeighbors: could not allocate nbr array") ;

    v->vtotal = v->vnum ;
    v->nsize = 1 ;
    for (i=0;i<v->vnum;i++)
      v->v[i] = vtmp[i];

    if (v->dist)      free(v->dist) ;
    if (v->dist_orig) free(v->dist_orig) ;

    v->dist = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist)
      ErrorExit(ERROR_NOMEMORY,
                "mrisFindNeighbors: could not allocate list of %d "
                "dists at v=%d", v->vnum, k) ;
    v->dist_orig = (float *)calloc(v->vnum, sizeof(float)) ;
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "mrisFindNeighbors: could not allocate list of %d "
                  "dists at v=%d", v->vnum, k) ;
  }

  for (k=0;k<mris->nfaces;k++) {
    f = &mris->faces[k];
    for (m=0;m<VERTICES_PER_FACE;m++) {
      v = &mris->vertices[f->v[m]];
      for (i=0;i<v->num && k!=v->f[i];i++);
      if (i==v->num)   /* face has vertex, but vertex doesn't have face */
        ErrorExit(ERROR_BADPARM,
                  "%s: face[%d].v[%d] = %d, but face %d not in vertex %d "
                  "face list\n", mris->fname,k,m,f->v[m], k, f->v[m]);
    }
  }

  for (vno = ntotal = vtotal = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag) continue ;
    vtotal += v->vtotal ;
    ntotal++ ;
  }
  mris->avg_nbrs = (float)vtotal / (float)ntotal ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRIScomputeNormals(MRI_SURFACE *mris)
{
  int       k,n, num, i ;
  VERTEX    *v ;
  FACE      *f;
  float     norm[3],snorm[3], len ;

  i = 0 ;
  for (k=0;k<mris->nfaces;k++) if (mris->faces[k].ripflag) {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
      mris->vertices[f->v[n]].border = TRUE;
  }

  for (k=0;k<mris->nvertices;k++) if (!mris->vertices[k].ripflag) {
    v = &mris->vertices[k];
    snorm[0]=snorm[1]=snorm[2]=0;
    v->area = 0;
    for (num = n=0;n<v->num;n++) if (!mris->faces[v->f[n]].ripflag) {
      num++ ;
      mrisNormalFace(mris, v->f[n], (int)v->n[n],norm);
      snorm[0] += norm[0];
      snorm[1] += norm[1];
      snorm[2] += norm[2];
      /* Note: overestimates area by *2 !! */
      /* I think this is wrong; should be 3; see below -- dhagler */
      v->area += mrisTriangleArea(mris, v->f[n], (int)v->n[n]);
    }
    if (!num)
      continue ;
    mrisNormalize(snorm);

    #if 0
        v->area /= 2.0 ;
    #else
        /* each vertex shares a face with 2 other vertices, so must divide
           area by 3, not 2! */
        v->area /= 3.0 ;
    #endif

    if (v->origarea<0)        /* has never been set */
      v->origarea = v->area;

    len = sqrt(snorm[0]*snorm[0] + snorm[1]*snorm[1] + snorm[2]*snorm[2]) ;
    if (!FZERO(len)) {
      v->nx = snorm[0];
      v->ny = snorm[1];
      v->nz = snorm[2];
      i = 0 ;
    }
    else {
      if (i++ > 5) continue ;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "vertex %d: degenerate normal\n", k) ;
      v->x += (float)randomNumber(-RAN, RAN) ;
      v->y += (float)randomNumber(-RAN, RAN) ;
      /* normal is always (0,0,+-1) anyway */
      if (mris->status == MRIS_PLANE || mris->status == MRIS_CUT) {
        v->nx = v->ny = 0.0f ; v->nz = 1.0f ;
        continue ;
      }

      v->z += (float)randomNumber(-RAN, RAN) ;
      for (n=0;n<v->vnum;n++) { /*if (!mris->faces[v->f[n]].ripflag)*/
        mris->vertices[v->v[n]].x += (float)randomNumber(-RAN, RAN) ;
        mris->vertices[v->v[n]].y += (float)randomNumber(-RAN, RAN) ;
        mris->vertices[v->v[n]].z += (float)randomNumber(-RAN, RAN) ;
      }
      k-- ;   /* recalculate the normal for this vertex */
    }
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisComputeVertexDistances(MRI_SURFACE *mris)
{
  int     vno, n, vtotal, *pv ;
  VERTEX  *v, *vn ;
  float   d, xd, yd, zd, circumference = 0.0f, angle ;
  VECTOR  *v1, *v2 ;

  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;

  for (vno=0;vno<mris->nvertices;vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->dist == NULL)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vtotal = v->vtotal ;
    switch (mris->status) {
    default:   /* don't really know what to do in other cases */
    case MRIS_PLANE:
      for (pv = v->v, n = 0 ; n < vtotal ; n++) {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        xd = v->x - vn->x ; yd = v->y - vn->y ; zd = v->z - vn->z ;
        d = xd*xd + yd*yd + zd*zd ;
        v->dist[n] = sqrt(d) ;
      }
      break ;
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
      VECTOR_LOAD(v1, v->x, v->y, v->z) ;  /* radius vector */
      if (FZERO(circumference))   /* only calculate once */
        circumference = M_PI * 2.0 * V3_LEN(v1) ;
      for (pv = v->v, n = 0 ; n < vtotal ; n++) {
        vn = &mris->vertices[*pv++] ;
        if (vn->ripflag)
          continue ;
        VECTOR_LOAD(v2, vn->x, vn->y, vn->z) ;  /* radius vector */
        angle = fabs(Vector3Angle(v1, v2)) ;
        d = circumference * angle / (2.0 * M_PI) ;
        v->dist[n] = d ;
      }
      break ;
    }
  }
  VectorFree(&v1) ; VectorFree(&v2) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
double
MRISaverageRadius(MRI_SURFACE *mris)
{
  double  radius ;
  int    vno, n ;
  VERTEX *vertex ;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0 ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    x = (double)vertex->x ; y = (double)vertex->y ; z = (double)vertex->z ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi)
      zhi=z;
    if (z<zlo)
      zlo=z;
    if (zlo < -1000)
      DiagBreak() ;
  }
  x0 = (xlo+xhi)/2.0f ; y0 = (ylo+yhi)/2.0f ; z0 = (zlo+zhi)/2.0f ;
  for (radius = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    n++ ;
    x = (double)vertex->x-x0 ;
    y = (double)vertex->y-y0 ;
    z = (double)vertex->z-z0 ;
    radius += sqrt(x*x + y*y + z*z) ;
  }

  radius /= (double)n ;
  return(radius) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisOrientSurface(MRI_SURFACE *mris)
{
  switch (mris->status) {
  case MRIS_RIGID_BODY:
  case MRIS_PARAMETERIZED_SPHERE:
  case MRIS_SPHERE:
  case MRIS_ELLIPSOID:
    MRISupdateEllipsoidSurface(mris) ;
    break ;
  case MRIS_PLANE:
    mrisOrientPlane(mris) ;
    break ;
  default:
    break ;
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRIScomputeTriangleProperties(MRI_SURFACE *mris)
{
  VECTOR  *v_a, *v_b, *v_n ;
  VERTEX  *v0, *v1, *v2, *va, *vb, *vo, *v ;
  FACE    *face ;
  int     fno, ano, vno  ;
  float   area, angle, dot, cross, dz ;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */

  mris->total_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;
    VERTEX_EDGE(v_a, v0, v1) ;
    VERTEX_EDGE(v_b, v0, v2) ;

    /* compute metric properties of first triangle */
    V3_CROSS_PRODUCT(v_a, v_b, v_n) ;
    area = V3_LEN(v_n) * 0.5f ;
    dot = V3_DOT(v_a, v_b) ;
    face->area = area ;
    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */
    face->nx = V3_X(v_n); face->ny = V3_Y(v_n); face->nz = V3_Z(v_n);
    mris->total_area += area ;

    /* now compute angles */
    VECTOR_LOAD(v_n, face->nx, face->ny, face->nz) ;
    if ((V3_X(v_n) < V3_Y(v_n)) && (V3_X(v_n) < V3_Z(v_n)))
      dz = fabs(V3_X(v_n)) ;
    else if (V3_Y(v_n) < V3_Z(v_n))
      dz = fabs(V3_Y(v_n)) ;
    else
      dz = fabs(V3_Z(v_n)) ;
    for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)  {
      switch (ano) {  /* vertices for triangle 1 */
      default:
      case 0: vo = v0 ; va = v2 ; vb = v1 ; break ;
      case 1: vo = v1 ; va = v0 ; vb = v2 ; break ;
      case 2: vo = v2 ; va = v1 ; vb = v0 ; break ;
      }

      VERTEX_EDGE(v_a, vo, va) ;VERTEX_EDGE(v_b, vo, vb) ;
      cross = VectorTripleProduct(v_b, v_a, v_n) ;
      dot = V3_DOT(v_a, v_b) ;
      angle = atan2(cross, dot) ;
      face->angle[ano] = angle ;
    }
  }

  /* calculate the "area" of the vertices */
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->area = 0.0 ;
    for (fno = 0 ; fno < v->num ; fno++)
      v->area += mris->faces[v->f[fno]].area ;
  #if 0
      v->area /= 2.0 ;
  #else
      /* each vertex shares a face with 2 other vertices, so must divide
         area by 3, not 2! */
      v->area /= 3.0 ;
  #endif
  }

  VectorFree(&v_a) ;
  VectorFree(&v_b) ;
  VectorFree(&v_n) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisComputeSurfaceDimensions(MRI_SURFACE *mris)
{
  int    vno ;
  VERTEX *vertex ;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    x = (double)vertex->x ; y = (double)vertex->y ; z = (double)vertex->z ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  mris->xlo = xlo ; mris->xhi = xhi ;
  mris->ylo = ylo ; mris->yhi = yhi ;
  mris->zlo = zlo ; mris->zhi = zhi ;
  mris->xctr = (xlo+xhi)/2.0f ;
  mris->yctr = (ylo+yhi)/2.0f ;
  mris->zctr = (zlo+zhi)/2.0f ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRIScomputeMetricProperties(MRI_SURFACE *mris)
{
  MRIScomputeNormals(mris) ;
  mrisComputeVertexDistances(mris) ;
  mrisComputeSurfaceDimensions(mris) ;
  MRIScomputeTriangleProperties(mris) ;  /* compute areas and normals */
  mrisOrientSurface(mris) ;
  if (mris->status == MRIS_PARAMETERIZED_SPHERE ||
      mris->status == MRIS_RIGID_BODY) {
    double old_area ;
    old_area = mris->total_area ;
    mris->total_area = M_PI * mris->radius * mris->radius * 4.0 ;
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISsaveVertexPositions(MRI_SURFACE *mris, int which)
{
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    switch (which) {
    case PIAL_VERTICES:
      v->pialx = v->x ; v->pialy = v->y ; v->pialz = v->z ;
      break ;
    case INFLATED_VERTICES:
      v->infx = v->x ; v->infy = v->y ; v->infz = v->z ;
      break ;
    case FLATTENED_VERTICES:
      v->fx = v->x ; v->fy = v->y ; v->fz = v->z ;
      break ;
    case CANONICAL_VERTICES:
      v->cx = v->x ; v->cy = v->y ; v->cz = v->z ;
      break ;
    case ORIGINAL_VERTICES:
      v->origx = v->x ; v->origy = v->y ; v->origz = v->z ;
      break ;
    default:
    case TMP_VERTICES:
      v->tx = v->x ; v->ty = v->y ; v->tz = v->z ;
      break ;
    }
  }
  if (which == CANONICAL_VERTICES) MRIScomputeCanonicalCoordinates(mris) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISstoreCurrentPositions(MRI_SURFACE *mris)
{
  return(MRISsaveVertexPositions(mris, TMP_VERTICES)) ;
}

/*---------------------------------------------------------------------------*/
static void
mrisNormalize(float v[3])
{
  float d;

  d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d>0) {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }
}

/*---------------------------------------------------------------------------*/
static int
mrisNormalFace(MRIS *mris, int fac,int n,float norm[])
{
  int     n0,n1, *pv ;
  FACE    *f;
  float   v0[3],v1[3];
  register VERTEX  *v, *vn0, *vn1 ;

  n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  pv = f->v ;
  vn0 = &mris->vertices[pv[n0]] ;
  vn1 = &mris->vertices[pv[n1]] ;
  v =  &mris->vertices[pv[n]] ;
  v0[0] = v->x - vn0->x; v0[1] = v->y - vn0->y; v0[2] = v->z - vn0->z;
  v1[0] = vn1->x - v->x; v1[1] = vn1->y - v->y; v1[2] = vn1->z - v->z;
  mrisNormalize(v0);
  mrisNormalize(v1);
  norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
  norm[1] = v1[0]*v0[2] - v0[0]*v1[2];
  norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static float
mrisTriangleArea(MRIS *mris, int fac, int n)
{
  int n0,n1;
  face_type *f;
  float v0[3],v1[3],d1,d2,d3;

  n0 = (n == 0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n == VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
  d1 = -v1[1]*v0[2] + v0[1]*v1[2];
  d2 = v1[0]*v0[2] - v0[0]*v1[2];
  d3 = -v1[0]*v0[1] + v0[0]*v1[1];
  return sqrt(d1*d1+d2*d2+d3*d3)/2;
}

/*---------------------------------------------------------------------------*/
static int
mrisOrientPlane(MRI_SURFACE *mris)
{
  int     fno, ano, vno ;
  FACE    *face ;
  VERTEX  *v ;

  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    if (face->ripflag) continue ;
    /* now give the area an orientation: if the unit normal is pointing
       downwards in the plane then the area should be negative.  */
    if (face->nz < 0.0f) { /* not in same dir, area < 0, reverse n */
      face->area *= -1.0f ;
      face->nx *= -1.0f; face->ny *= -1.0f; face->nz *= -1.0f;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        face->angle[ano] *= -1.0f ;
    }
  }
  /* now recompute the total surface area, ignoring negative areas */
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    if (face->ripflag) continue ;
    if (face->area >= 0.0f)
      mris->total_area += face->area ;
    else {
      mris->neg_area += -face->area ;
      mris->neg_orig_area += face->orig_area ;
    }
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)  continue ;
    if (v->nz < 0) {
      v->nz *= -1 ;
      v->neg = 1 ;
    }
    else v->neg = 0 ;
    v->area = 0 ;
    for (fno = 0 ; fno < v->num ; fno++) {
      face = &mris->faces[v->f[fno]] ;
      v->area += face->area ;
    }
#if 0
    v->area /= 2.0 ;
#else
    /* each vertex shares a face with 2 other vertices, so must divide
       area by 3, not 2! */
    v->area /= 3.0 ;
#endif
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisOrientEllipsoid(MRI_SURFACE *mris)
{
  int     fno, ano ;
  VERTEX  *v ;
  FACE    *face ;
  float   dot ;

  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    if (face->ripflag) continue ;
    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative. */
    v = &mris->vertices[face->v[0]] ;
    dot = v->x * face->nx + v->y * face->ny + v->z * face->nz;
    if (dot < 0.0f) {  /* not in same direction, area < 0 and reverse n */
      face->area *= -1.0f ;
      face->nx *= -1.0f; face->ny *= -1.0f; face->nz *= -1.0f;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        face->angle[ano] *= -1.0f ;
    }
  }
  /* now recompute the total surface area, ignoring negative areas */
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    if (face->ripflag) continue ;
    if (face->area >= 0.0f)
      mris->total_area += face->area ;
    else {
      mris->neg_area += -face->area ;
      mris->neg_orig_area += face->orig_area ;
    }
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRIScomputeCanonicalCoordinates(MRI_SURFACE *mris)
{
  float   theta, phi, r, d, x, y, z ;
  VERTEX  *v ;
  int     vno ;

  r = mris->radius = MRISaverageCanonicalRadius(mris) ;
  r = mris->radius = (float)nint(mris->radius) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    x = v->cx ; y = v->cy ; z = v->cz ;
    theta = atan2(y/r, x/r) ;
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
    d = r*r-z*z ; if (d < 0.0) d = 0.0 ;
    phi = atan2(sqrt(d), z) ;
    v->theta = theta ; v->phi = phi ;
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISupdateEllipsoidSurface(MRI_SURFACE *mris)
{
  mrisOrientEllipsoid(mris) ;      /* orient the normals and angles */
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
double
MRISaverageCanonicalRadius(MRI_SURFACE *mris)
{
  double  radius ;
  int    vno, n ;
  VERTEX *vertex ;
  double x, y, z, xlo, ylo, zlo, xhi, yhi, zhi, x0, y0, z0 ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    x = (double)vertex->cx ; y = (double)vertex->cy ; z = (double)vertex->cz ;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  x0 = (xlo+xhi)/2.0f ; y0 = (ylo+yhi)/2.0f ; z0 = (zlo+zhi)/2.0f ;
  for (radius = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    n++ ;
    x = (double)vertex->cx-x0 ;
    y = (double)vertex->cy-y0 ;
    z = (double)vertex->cz-z0 ;
    radius += sqrt(x*x + y*y + z*z) ;
  }
  radius /= (double)n ;
  return(radius) ;
}

/*---------------------------------------------------------------------------*/
MRI_SURFACE *
MRISalloc(int nvertices, int nfaces)
{
  MRI_SURFACE   *mris ;

  mris = (MRI_SURFACE *)calloc(1, sizeof(MRI_SURFACE)) ;
  if (!mris)
    ErrorExit(ERROR_NOMEMORY,
              "MRISalloc(%d, %d): could not allocate mris structure",
              nvertices, nfaces);

  mris->nsize = 1 ;  /* only 1-connected neighbors initially */
  mris->nvertices = nvertices ;
  mris->nfaces = nfaces ;
  mris->vertices = (VERTEX *)calloc(nvertices, sizeof(VERTEX)) ;
  if (!mris->vertices)
    ErrorExit(ERROR_NOMEMORY,
              "MRISalloc(%d, %d): could not allocate vertices",
              nvertices, sizeof(VERTEX));
  mris->faces = (FACE *)calloc(nfaces, sizeof(FACE)) ;
  if (!mris->faces)
    ErrorExit(ERROR_NOMEMORY,
              "MRISalloc(%d, %d): could not allocate faces",
              nfaces, sizeof(FACE));
  return(mris) ;
}

/*---------------------------------------------------------------------------*/
MRI_SURFACE  *
MRISoverAlloc(int max_vertices, int max_faces, int nvertices, int nfaces)
{
  MRI_SURFACE  *mris ;

  if (max_vertices <= 0)
    max_vertices = nvertices ;
  if (max_faces <= 0)
    max_faces = nfaces ;
  mris = MRISalloc(max_vertices, max_faces) ;
  mris->nvertices = nvertices ;
  mris->nfaces = nfaces ;
  mris->max_vertices = max_vertices ;
  mris->max_faces = max_faces ;
  return(mris) ;
}

/*---------------------------------------------------------------------------*/
MRI_SURFACE *
MRISreadOverAlloc(char *fname, double pct_over)
{
  MRI_SURFACE *mris = NULL ;
  int         nquads, nvertices, magic, version, ix, iy, iz, vno, fno, n, m,
              imnr, imnr0, imnr1, type, vertices[VERTICES_PER_FACE+1], num ;
  float       x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  FILE        *fp = NULL ;
  VERTEX      *vertex ;
  FACE        *face ;

int which ;
char *surf_name ;

  /*chklc();*/    /* check to make sure license.dat is present */
  type = mrisFileNameType(fname);
  if (type == MRIS_ASCII_TRIANGLE_FILE) {
    mris = mrisReadAsciiFile(fname);
    if (!mris) return(NULL);
    version = -3;
  }
  else if (type == MRIS_ICO_FILE) {
    mris = ICOreadOverAlloc(fname, pct_over) ;
    if (!mris) return(NULL);
    return(mris);
    /*version = -2;*/
  }
  else if (type == MRIS_GEO_TRIANGLE_FILE) {
    mris = mrisReadGeoFile(fname);
    if (!mris) return(NULL);
    version = -4;
  }
  else {
    fp = fopen(fname, "rb") ;
    if (!fp)
     ErrorReturn(NULL,(ERROR_NOFILE,"MRISread(%s):could not open file",fname));

    fread3(&magic, fp);
    if (magic == QUAD_FILE_MAGIC_NUMBER) {
      version = -1;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stdout, "new surface file format\n");
    }
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) {
      version = -2 ;
    }
    else if (magic == TRIANGLE_FILE_MAGIC_NUMBER) {
      fclose(fp) ;
      mris = mrisReadTriangleFile(fname, pct_over) ;
      if (!mris) ErrorReturn(NULL,(Gerror, "mrisReadTriangleFile failed.\n"));
      version = -3 ;
    }
    else  {
      rewind(fp);
      version = 0;
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        printf("surfer: old surface file format\n");
    }
  }
  if (version >= -2) { /* some type of quadrangle file */
    fread3(&nvertices, fp);
    fread3(&nquads, fp);   /* # of qaudrangles - not triangles */

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,"reading %d vertices and %d faces.\n",
      nvertices,2*nquads);

#if 0
    mris = MRISoverAlloc(pct_over*nvertices,pct_over*2*nquads,nvertices,
       2*nquads);
#endif
    mris = MRISoverAlloc((int)(pct_over*nvertices),(int)(pct_over*2*nquads),nvertices,
       2*nquads);
    mris->type = MRIS_BINARY_QUADRANGLE_FILE ;

    imnr0 = 1000 ;
    imnr1 = 0 ;
    for (vno = 0 ; vno < nvertices ; vno++) {
      vertex = &mris->vertices[vno];
      if (version == -1) {
        fread2(&ix,fp);
        fread2(&iy,fp);
        fread2(&iz,fp);
        vertex->x = ix/100.0;
        vertex->y = iy/100.0;
        vertex->z = iz/100.0;
      }
      else { /* version == -2 */
        vertex->x = freadFloat(fp);
        vertex->y = freadFloat(fp);
        vertex->z = freadFloat(fp);
      }
      imnr = (int)((vertex->y-START_Y)/SLICE_THICKNESS+0.5);
      if (imnr > imnr1) imnr1 = imnr ;
      if (imnr < imnr0) imnr0 = imnr ;
      if (version == 0) { /* old surface format */
        fread1(&num,fp);   /* # of faces we are part of */
        vertex->num = num ;
        vertex->f = (int *)calloc(vertex->num,sizeof(int));
        if (!vertex->f)
          ErrorExit(ERROR_NOMEMORY, "MRISread: could not allocate %d faces",
                    vertex->num) ;
        vertex->n = (uchar *)calloc(vertex->num,sizeof(uchar));
        if (!vertex->n)
          ErrorExit(ERROR_NOMEMORY, "MRISread: could not allocate %d nbrs",
                    vertex->n) ;
        for (n=0;n<vertex->num;n++)
          fread3(&vertex->f[n],fp);
      }
      else vertex->num = 0;   /* will figure it out */
    }

    for (fno = 0 ; fno < mris->nfaces ; fno += 2) {
      if (fno == 86) DiagBreak() ;
      for (n = 0 ; n < 4 ; n++) {  /* read quandrangular face */
        fread3(&vertices[n],fp);
        if (vertices[n] == 22) DiagBreak() ;
      }
/* if we're going to be arbitrary, we might as well be really arbitrary */
#define WHICH_FACE_SPLIT(vno0, vno1) (1*nint(sqrt(1.9*vno0) + sqrt(3.5*vno1)))
      /* NOTE: for this to work properly in the write, the first two
         vertices in the first face (EVEN and ODD) must be 0 and 1. */
      which = WHICH_FACE_SPLIT(vertices[0], vertices[1]) ;

      /* 1st triangle */
      if (EVEN(which)) {
        mris->faces[fno].v[0] = vertices[0] ;
        mris->faces[fno].v[1] = vertices[1] ;
        mris->faces[fno].v[2] = vertices[3] ;
        /* 2nd triangle */
        mris->faces[fno+1].v[0] = vertices[2] ;
        mris->faces[fno+1].v[1] = vertices[3] ;
        mris->faces[fno+1].v[2] = vertices[1] ;
      }
      else {
        mris->faces[fno].v[0] = vertices[0] ;
        mris->faces[fno].v[1] = vertices[1] ;
        mris->faces[fno].v[2] = vertices[2] ;
        /* 2nd triangle */
        mris->faces[fno+1].v[0] = vertices[0] ;
        mris->faces[fno+1].v[1] = vertices[2] ;
        mris->faces[fno+1].v[2] = vertices[3] ;
      }
      for (n = 0 ; n < VERTICES_PER_FACE ; n++) {
        mris->vertices[mris->faces[fno].v[n]].num++;
        mris->vertices[mris->faces[fno+1].v[n]].num++;
      }
    }
    fclose(fp);
  }
  strcpy(mris->fname, fname) ;

  surf_name = strrchr(fname, '/') ;
  if (surf_name == NULL) surf_name = fname ;
  else                   surf_name++ ;  /* past the last slash */
  if (toupper(*surf_name) == 'R') mris->hemisphere = RIGHT_HEMISPHERE ;
  else                            mris->hemisphere = LEFT_HEMISPHERE ;

  if ((version<0) || type == MRIS_ASCII_TRIANGLE_FILE) {
    for (vno = 0 ; vno< mris->nvertices ; vno++) {
      vertex = &mris->vertices[vno] ;
      mris->vertices[vno].f =
        (int *)calloc(mris->vertices[vno].num,sizeof(int));
      if (!mris->vertices[vno].f)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d faces at %dth vertex",
                  fname, vno, mris->vertices[vno].num) ;

      mris->vertices[vno].n =
        (uchar *)calloc(mris->vertices[vno].num,sizeof(uchar));
      if (!mris->vertices[vno].n)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISread(%s): could not allocate %d indices at %dth vertex",
                  fname, vno, mris->vertices[vno].num) ;
      mris->vertices[vno].num = 0 ;
    }
    for (fno = 0 ; fno < mris->nfaces ; fno++) {
      face = &mris->faces[fno] ;
      for (n = 0 ; n < VERTICES_PER_FACE ; n++)
        mris->vertices[face->v[n]].f[mris->vertices[face->v[n]].num++] = fno;
    }
  }

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    mris->vertices[vno].curv = 0;
    mris->vertices[vno].origarea = -1;
    mris->vertices[vno].border = 0;
    for (n=0;n<mris->vertices[vno].num;n++) {
      for (m=0;m<VERTICES_PER_FACE;m++) {
        if (mris->faces[mris->vertices[vno].f[n]].v[m] == vno)
          mris->vertices[vno].n[n] = m;
      }
    }
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }
  mris->xlo = xlo ; mris->ylo = ylo ; mris->zlo = zlo ;
  mris->xhi = xhi ; mris->yhi = yhi ; mris->zhi = zhi ;
  mris->xctr = (xhi+xlo)/2;
  mris->yctr = (yhi+ylo)/2;
  mris->zctr = (zhi+zlo)/2;

  mrisFindNeighbors(mris);
  MRIScomputeNormals(mris);
  mrisComputeVertexDistances(mris) ;
  mrisReadTransform(mris, fname) ;
  mris->radius = MRISaverageRadius(mris) ;
  MRIScomputeMetricProperties(mris) ;
  MRISstoreCurrentPositions(mris) ;

  return(mris) ;
}

/*---------------------------------------------------------------------------*/
MRI_SURFACE *
MRISread(char *fname)
{
  MRI_SURFACE  *mris;

  mris = MRISreadOverAlloc(fname, 0.0);
  return(mris);
}

/*---------------------------------------------------------------------------*/
static MRI_SURFACE *
mrisReadAsciiFile(char *fname)
{
  MRI_SURFACE   *mris ;
  char    line[STRLEN], *cp ;
  int     vno, fno, n, nvertices, nfaces, patch, rip ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
              (ERROR_NOFILE,
               "MRISreadAsciiFile: could not open file %s",fname));
  patch = 0 ;
  cp = fgetl(line, STRLEN, fp) ;
  sscanf(cp, "%d %d\n", &nvertices, &nfaces) ;
  mris = MRISalloc(nvertices, nfaces) ;
  mris->type = MRIS_TRIANGULAR_SURFACE ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    fscanf(fp, "%f  %f  %f  %d\n", &v->x, &v->y, &v->z, &rip) ;
    v->ripflag = rip ;
    if (v->ripflag)
      patch = 1 ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++) {
      fscanf(fp, "%d ", &face->v[n]) ;
      mris->vertices[face->v[n]].num++;
    }
    fscanf(fp, "%d\n", &rip) ;
    face->ripflag = rip ;
  }

  mris->patch = patch ;
  if (patch)
    mris->status = MRIS_PLANE ;
  fclose(fp) ;
  return(mris) ;
}

/*---------------------------------------------------------------------------*/
static MRI_SURFACE *
mrisReadTriangleFile(char *fname, double pct_over)
{
  VERTEX      *v ;
  FACE        *f ;
  int         nvertices, nfaces, magic, vno, fno, n ;
  char        line[STRLEN] ;
  FILE        *fp ;
  MRI_SURFACE *mris ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL,(ERROR_NOFILE,
                      "mrisReadTriangleFile(%s): could not open file",fname));

  fread3(&magic, fp) ;
  fgets(line, 200, fp) ;
  fscanf(fp, "\n") ;
  nvertices = freadInt(fp);
  nfaces = freadInt(fp);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout,"surface %s: %d vertices and %d faces.\n",
            fname, nvertices,nfaces);

#if 0
  mris = MRISoverAlloc((pct_over*nvertices, pct_over*nfaces,nvertices,nfaces) ;
#endif
  mris = MRISoverAlloc((int)(pct_over*nvertices), (int)(pct_over*nfaces),nvertices,nfaces) ;
  mris->type = MRIS_TRIANGULAR_SURFACE ;

  for (vno = 0 ; vno < nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->x = freadFloat(fp);
    v->y = freadFloat(fp);
    v->z = freadFloat(fp);
    v->num = 0;   /* will figure it out */
    if (fabs(v->x) > 10000 || !finite(v->x))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d x coordinate %f!",
                "mrisReadTriangleFile", vno, v->x) ;
    if (fabs(v->y) > 10000 || !finite(v->y))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d y coordinate %f!",
                "mrisReadTriangleFile", vno, v->y) ;
    if (fabs(v->z) > 10000 || !finite(v->z))
      ErrorExit(ERROR_BADFILE, "%s: vertex %d z coordinate %f!",
                "mrisReadTriangleFile", vno, v->z) ;
  }

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++) {
      f->v[n] = freadInt(fp);
      if (f->v[n] >= mris->nvertices || f->v[n] < 0)
        ErrorExit(ERROR_BADFILE, "f[%d]->v[%d] = %d - out of range!\n",
                  fno, n, f->v[n]) ;
    }

    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      mris->vertices[mris->faces[fno].v[n]].num++;
  }
  fclose(fp);
  return(mris) ;
}

/*---------------------------------------------------------------------------*/
static MRI_SURFACE *
mrisReadGeoFile(char *fname)
{
  MRI_SURFACE   *mris ;
  char    line[202], *cp ;
  int     vno, fno, n, nvertices, nfaces, patch, vertices_per_face, nedges ;
  VERTEX  *v ;
  FACE    *face ;
  FILE    *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
              (ERROR_NOFILE,
               "mrisReadGeoFile: could not open file %s",fname));

  patch = 0 ;
  cp = fgetl(line, 100, fp) ;
  if (sscanf(cp, "%*d %d %d %d\n", &nvertices, &nfaces, &nedges) != 3) {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "mrisReadGeoFile(%s): could not scan "
                       "dimensions from '%s'", fname, cp)) ;
  }
  vertices_per_face = nedges / nfaces ;
  if (vertices_per_face != VERTICES_PER_FACE) {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE,
                       "mrisReadGeoFile(%s): unsupported vertices/face %d.",
                       fname, vertices_per_face)) ;
  }

  cp = fgetl(line, 200, fp) ;   /* nfaces again */
  mris = MRISalloc(nvertices, nfaces) ;
  mris->type = MRIS_GEO_TRIANGLE_FILE ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    fscanf(fp, "%e %e %e", &v->x, &v->y, &v->z) ;
    if (ISODD(vno))
      fscanf(fp, "\n") ;
  }
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    int tmp ;
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE-1 ; n++) {
      fscanf(fp, "%d ", &face->v[n]) ; face->v[n]-- ;   /* make it 0-based */
      if (face->v[n] < 0)
        DiagBreak() ;
      mris->vertices[face->v[n]].num++;
    }
    n = VERTICES_PER_FACE-1 ;   /* already true - but make it explicit */
    fscanf(fp, "-%d\n", &face->v[n]) ; face->v[n]-- ;   /* make it 0-based */
    mris->vertices[face->v[n]].num++;

    /* swap positions so normal (via cross-product) will point outwards */
    tmp = face->v[1] ;
    face->v[1] = face->v[2] ; face->v[2] = tmp ;
  }

  fclose(fp) ;
  return(mris) ;
}


/*---------------------------------------------------------------------------*/
int
MRISreadValues(MRI_SURFACE *mris, char *sname)
{
  char funcname[STRLEN]="MRISreadValues";
  int   i,k,num,ilat, vno ;
  float f;
  float lat, *cvec;
  FILE  *fp;
  char  *cp, fname[STRLEN] ;
  int maxk=-1;
  double sum=0,sum2=0,max= -1000,min=1000;

  cvec = MRISreadCurvatureVector(mris, sname) ;
  if (cvec)
  {
    printf("reading values from curvature-format file...\n") ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      if (vno == Gdiag_no)
        DiagBreak() ;
      mris->vertices[vno].val = cvec[vno];
    }
    free(cvec) ;
  }
  else
  {
    strcpy(fname, sname) ;
    cp = strrchr(fname, '.') ;
    if (!cp || *(cp+1) != 'w')
      strcat(fname, ".w") ;
    fp = fopen(fname,"rb");
    if (fp==NULL) 
      ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE,
                                 "%s: File %s not found\n",funcname,fname));
    fread2(&ilat,fp);
    lat = ilat/10.0;
    
    for (k=0;k<mris->nvertices;k++) {
      mris->vertices[k].val=0;
      mris->vertices[k].undefval = 1;
      mris->vertices[k].fixedval = 0;
    }
    if (fread3(&num,fp) < 1)
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE, 
                   "%s(%s): could not read # of vertices",
                   funcname,fname)) ;
    for (i=0;i<num;i++)
    {
      if (fread3(&k,fp) < 1)
        ErrorReturn(ERROR_BADFILE,
                    (ERROR_BADFILE, 
                     "%s(%s): could not read %dth vno",
                     funcname,fname, i)) ;
      f = freadFloat(fp) ;
      if (k>maxk) maxk=k;
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;

      if (k>=mris->nvertices || k<0)
        printf("%s: vertex index out of range: %d f=%f\n",funcname,k,f);
      else
      {
        if (k == Gdiag_no)
          DiagBreak() ;
        mris->vertices[k].val = f;
        mris->vertices[k].undefval = 0;
        mris->vertices[k].fixedval = 1;
      }
    }
    fclose(fp);
  }
  
  if(num>0) {
    sum /= num;
    sum2 = sqrt(sum2/num-sum*sum);
  } else {
    sum = sum2 = 0;
  }
  MsgPrintf("%s: file %s read \n",funcname,fname);
  MsgPrintf("      with %d vals and max vertex index = %d\n",num,maxk);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(NO_ERROR) ;
}


/*---------------------------------------------------------------------------*/
int
MRISreadImagValues(MRI_SURFACE *mris, char *sname)
{
  char funcname[STRLEN]="MRISreadImagValues";
  int i,k,num,ilat;
  float f;
  float lat;
  FILE *fp;
  int maxk=-1;
  double sum=0,sum2=0,max= -1000,min=1000;
  char  *cp, fname[STRLEN] ;

  strcpy(fname, sname) ;
  cp = strrchr(fname, '.') ;
  if (!cp || *(cp+1) != 'w')
    strcat(fname, ".w") ;
  fp = fopen(fname,"r");
  if (fp==NULL) 
    ErrorReturn(ERROR_NOFILE, 
              (ERROR_NOFILE,"%s: File %s not found\n",funcname,fname));
  fread2(&ilat,fp);
  lat = ilat/10.0;

  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].imag_val=0;
  fread3(&num,fp);
  for (i=0;i<num;i++)
  {
    fread3(&k,fp);
    f = freadFloat(fp);

    if (k>maxk) maxk=k;
    sum += f;
    sum2 += f*f;
    if (f>max) max=f;
    if (f<min) min=f;
    
    if (k>=mris->nvertices||k<0)
      printf("%s: vertex index out of range: %d f=%f\n",funcname,k,f);
    else
      mris->vertices[k].imag_val = f;
  }
  fclose(fp);
  
  if(num>0) {
    sum /= num;
    sum2 = sqrt(sum2/num-sum*sum);
  } else {
    sum = sum2 = 0;
  }
  MsgPrintf("%s: file %s read \n",funcname,fname);
  MsgPrintf("      with %d vals and max vertex index = %d\n",num,maxk);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
float *
MRISreadCurvatureVector(MRI_SURFACE *mris, char *sname)
{
  int    k,i,vnum,fnum;
  float  *cvec ;
  FILE   *fp;
  char   *cp, path[STRLEN], fname[STRLEN] ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (cp)
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
    fprintf(stdout, "reading curvature file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    return(NULL) ;

  fread3(&vnum,fp);
  if (vnum == NEW_VERSION_MAGIC_NUMBER)
  {
    fclose(fp) ;
    return(MRISreadNewCurvatureVector(mris, fname)) ;
  }
  
  fread3(&fnum,fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    return(NULL) ;
  }
  cvec = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (!cvec)
    ErrorExit(ERROR_NOMEMORY, "MRISreadCurvatureVector(%s): calloc failed",
              fname) ;

  for (k=0;k<vnum;k++)
  {
    fread2(&i,fp);
    cvec[k] = i/100.0 ;
  }
  fclose(fp);
  return(cvec) ;
}

/*---------------------------------------------------------------------------*/
float *
MRISreadNewCurvatureVector(MRI_SURFACE *mris, char *sname)
{
  int    k,vnum,fnum, vals_per_vertex ;
  float  *cvec ;
  FILE   *fp;
  char   *cp, path[STRLEN], fname[STRLEN] ;
  
  cp = strchr(sname, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    cp = strchr(sname, '.') ;
    FileNamePath(mris->fname, path) ;
    if (cp)
      sprintf(fname, "%s/%s", path, sname) ;
    else   /* no hemisphere specified */
      sprintf(fname, "%s/%s.%s", path, 
              mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", sname) ;
  }
  else   
    strcpy(fname, sname) ;  /* path specified explcitly */

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
    fprintf(stdout, "reading curvature file...") ;

  fp = fopen(fname,"r");
  if (fp==NULL) 
    return(NULL) ;

  fread3(&vnum,fp);
  if (vnum != NEW_VERSION_MAGIC_NUMBER)
  {
    fclose(fp) ;
    return(MRISreadCurvatureVector(mris, fname)) ;
  }
  
  vnum = freadInt(fp);
  fnum = freadInt(fp);
  if (vnum!= mris->nvertices)
  {
    fclose(fp) ;
    return(NULL) ;
  }
  vals_per_vertex = freadInt(fp) ;
  if (vals_per_vertex != 1)
  {
    fclose(fp) ;
    return(NULL) ;
  }

  cvec = (float *)calloc(mris->nvertices, sizeof(float)) ;
  if (!cvec)
    ErrorExit(ERROR_NOMEMORY, "MRISreadNewCurvatureVector(%s): calloc failed",
              fname) ;
  for (k=0;k<vnum;k++)
  {
    cvec[k] = freadFloat(fp) ;
  }
  fclose(fp);
  return(cvec) ;
}

/*---------------------------------------------------------------------------*/
int
MRISremoveTriangleLinks(MRI_SURFACE *mris)
{
  int    fno, which ;
  FACE   *f ;

  if (!IS_QUADRANGULAR(mris))
    return(NO_ERROR) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "removing non-quadrangular links.\n") ;

  for (fno = 0 ; fno < mris->nfaces ; fno += 2)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    which = WHICH_FACE_SPLIT(f->v[0], f->v[1]) ;
    if (EVEN(which))
    {
      mrisRemoveVertexLink(mris, f->v[1], f->v[2]) ;
      mrisRemoveVertexLink(mris, f->v[2], f->v[1]) ;
    }
    else
    {
      mrisRemoveVertexLink(mris, f->v[0], f->v[2]) ;
      mrisRemoveVertexLink(mris, f->v[2], f->v[0]) ;
    }
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisRemoveVertexLink(MRI_SURFACE *mris, int vno1, int vno2)
{
  int    n ;
  VERTEX *v ;

  v = &mris->vertices[vno1] ;
  for (n = 0 ; n < v->vnum ; n++)
    if (v->v[n] == vno2)
      break ;

  if (n < v->vnum)
  {
    memmove(v->v+n, v->v+n+1, (v->vtotal-(n+1))*sizeof(int)) ;
    v->vnum-- ; v->vtotal-- ;
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
/* won't undo MRISremoveRipped */
int
MRISunrip(MRI_SURFACE *mris)
{
  int    vno, fno ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].ripflag = 0 ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
    mris->faces[fno].ripflag = 0 ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISfree(MRI_SURFACE **pmris)
{
  MRI_SURFACE  *mris ;
  int          vno ;

  mris = *pmris ;
  *pmris = NULL ;

  if (mris->dx2)
    free(mris->dx2) ;
  if (mris->dy2)
    free(mris->dy2) ;
  if (mris->dz2)
    free(mris->dz2) ;
  if (mris->labels)
    free(mris->labels) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (mris->vertices[vno].f)
      free(mris->vertices[vno].f) ;
    if (mris->vertices[vno].n)
      free(mris->vertices[vno].n) ;
    if (mris->vertices[vno].dist)
      free(mris->vertices[vno].dist) ;
    if (mris->vertices[vno].dist_orig)
      free(mris->vertices[vno].dist_orig) ;
    if (mris->vertices[vno].v)
      free(mris->vertices[vno].v) ;
  }

  if (mris->vertices)
    free(mris->vertices) ;
  if (mris->faces)
    free(mris->faces) ;
  free(mris) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISreadPatchNoRemove(MRI_SURFACE *mris, char *pname)
{

  int         ix, iy, iz, k, i, j, npts ;
  FILE        *fp ;
  char        fname[STRLEN] ;

#if 0
  char        path[STRLEN], *cp ;
  cp = strchr(pname, '/') ;
  if (cp)
    strcpy(fname, pname) ;    /* path already specified */
  else                        /* no path - use same as was used in MRISread */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/%s", path, pname) ;
  }
#else
  MRISbuildFileName(mris, pname, fname) ;
#endif
  fp = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,
                      "MRISreadPatch(%s): could not open file", fname));


  npts = freadInt(fp) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "reading patch %s with %d vertices (%2.1f%% of total)\n",
            pname, npts, 100.0f*(float)npts/(float)mris->nvertices) ;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].ripflag = TRUE;
  for (j=0;j<npts;j++)
  {
    i = freadInt(fp) ;
    if (i<0)
    {
      k = -i-1;
      if (k < 0 || k >= mris->nvertices)
        ErrorExit(ERROR_BADFILE, 
                  "MRISreadPatch: bad vertex # (%d) found in patch file", k) ;
      mris->vertices[k].border = TRUE;
    } 
    else
    {
      k = i-1;
      if (k < 0 || k >= mris->nvertices)
        ErrorExit(ERROR_BADFILE, 
                  "MRISreadPatch: bad vertex # (%d) found in patch file", k) ;
      mris->vertices[k].border = FALSE;
    }
    mris->vertices[k].ripflag = FALSE;
    fread2(&ix,fp);
    fread2(&iy,fp);
    fread2(&iz,fp);
    mris->vertices[k].x = ix/100.0;
    mris->vertices[k].y = iy/100.0;
    mris->vertices[k].z = iz/100.0;
    if (mris->vertices[k].x > mris->xhi) mris->xhi = mris->vertices[k].x;
    if (mris->vertices[k].x < mris->xlo) mris->xlo = mris->vertices[k].x;
    if (mris->vertices[k].y > mris->yhi) mris->yhi = mris->vertices[k].y;
    if (mris->vertices[k].y < mris->ylo) mris->ylo = mris->vertices[k].y;
    if (mris->vertices[k].z > mris->zhi) mris->zhi = mris->vertices[k].z;
    if (mris->vertices[k].z < mris->zlo) mris->zlo = mris->vertices[k].z;
    if (k == Gdiag_no && Gdiag & DIAG_SHOW)
      fprintf(stdout, "vertex %d read @ (%2.2f, %2.2f, %2.2f)\n",k,
              mris->vertices[k].x,mris->vertices[k].y,mris->vertices[k].z) ;
  }
  fclose(fp);
  MRISripFaces(mris);
  mris->patch = 1 ;
  mris->status = MRIS_CUT ;

  mrisComputeBoundaryNormals(mris);
  MRISupdateSurface(mris) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/

int
MRISreadPatch(MRI_SURFACE *mris, char *pname)
{
  int ret  ;

  ret = MRISreadPatchNoRemove(mris, pname) ;
  if (ret != NO_ERROR)
    return(ret) ;
  MRISremoveRipped(mris) ;
  MRISupdateSurface(mris) ;

  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISripFaces(MRI_SURFACE *mris)
{
  int n,k;
  face_type *f;

  for (k=0;k<mris->nfaces;k++)
    mris->faces[k].ripflag = FALSE;
  for (k=0;k<mris->nfaces;k++)
  {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
      if (mris->vertices[f->v[n]].ripflag)
        f->ripflag = TRUE;
  }
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].border = FALSE;
  for (k=0;k<mris->nfaces;k++)
  if (mris->faces[k].ripflag)
  {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
        mris->vertices[f->v[n]].border = TRUE;
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
static int
mrisComputeBoundaryNormals(MRI_SURFACE *mris)
{
#if 0
  int      k,m,n;
  VERTEX   *v;
  float    sumx,sumy,r,nx,ny,f;

  for (k=0;k<mris->nvertices;k++)
  if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
  {
    v = &mris->vertices[k];
    n = 0;
    sumx = 0;
    sumy = 0;
    for (m=0;m<v->vnum;m++)
    if (!mris->vertices[v->v[m]].ripflag)
    {
      sumx += v->x-mris->vertices[v->v[m]].x;
      sumy += v->y-mris->vertices[v->v[m]].y;
      n++;
    }
    v->bnx = (n>0)?sumx/n:0;
    v->bny = (n>0)?sumy/n:0;
  }
  for (k=0;k<mris->nvertices;k++)
  if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
  {
    v = &mris->vertices[k];
    n = 0;
    sumx = 0;
    sumy = 0;
    for (m=0;m<v->vnum;m++)
    if ((!mris->vertices[v->v[m]].ripflag)&&mris->vertices[v->v[m]].border)
    {
      nx = -(v->y-mris->vertices[v->v[m]].y); 
      ny = v->x-mris->vertices[v->v[m]].x; 
      f = nx*v->bnx+ny*v->bny;
/*
      f = (f<0)?-1.0:(f>0)?1.0:0.0;
*/
      sumx += f*nx;
      sumy += f*ny;
      n++;
    }
    v->bnx = (n>0)?sumx/n:0;
    v->bny = (n>0)?sumy/n:0;
  }
  for (k=0;k<mris->nvertices;k++)
  if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
  {
    r = sqrt(SQR(mris->vertices[k].bnx)+SQR(mris->vertices[k].bny));
    if (r>0)
    {
      mris->vertices[k].bnx /= r;
      mris->vertices[k].bny /= r;
    }
  }
#endif
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISupdateSurface(MRI_SURFACE *mris)
{
  MRIScomputeMetricProperties(mris) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/

int
MRISremoveRipped(MRI_SURFACE *mris)
{
  int     vno, n, fno, nripped ;
  VERTEX  *v ;
  FACE    *face ;

  if (Gdiag & DIAG_SHOW)
    fprintf(stdout, "removing ripped vertices and faces...\n") ;
  do
  {
    nripped = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        if (v->dist)
          free(v->dist) ;
        if (v->dist_orig)
          free(v->dist_orig) ;
        v->dist = v->dist_orig = NULL ;
        continue ;
      }

      for (n = 0 ; n < v->vtotal ; n++)
      {
        /* remove this vertex from neighbor list if it is ripped */
        if (mris->vertices[v->v[n]].ripflag)
        {
          if (n < v->vtotal-1)  /* not the last one in the list */
          {
            memmove(v->v+n, v->v+n+1, (v->vtotal-n-1)*sizeof(int)) ;
            memmove(v->dist+n, v->dist+n+1, (v->vtotal-n-1)*sizeof(float)) ;
            memmove(v->dist_orig+n, v->dist_orig+n+1, 
                    (v->vtotal-n-1)*sizeof(float)) ;
          }
          if (n < v->vnum)      /* it was a 1-neighbor */
            v->vnum-- ;
          if (n < v->v2num)     /* it was a 2-neighbor */
            v->v2num-- ;
          if (n < v->v3num)     /* it was a 3-neighbor */
            v->v3num-- ;
          n-- ; v->vtotal-- ;
        }
      }
      for (fno = 0 ; fno < v->num ; fno++)
      {
        /* remove this face from face list if it is ripped */
        if (mris->faces[v->f[fno]].ripflag)
        {
          if (fno < v->num-1)  /* not the last one in the list */
          {
            memmove(v->f+fno, v->f+fno+1, (v->num-fno-1)*sizeof(int)) ;
            memmove(v->n+fno, v->n+fno+1, (v->num-fno-1)*sizeof(uchar)) ;
          }
          v->num-- ; fno-- ;
        }
      }
      if (v->num <= 0 || v->vnum <= 0)  /* degenerate vertex */
      {
        v->ripflag = 1 ;
        nripped++ ;
      }
    }
  } while (nripped > 0) ;

  /* now recompute total original area for scaling */
  mris->orig_area = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;
    mris->orig_area += face->orig_area ;
  }


  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/

static int
mrisComputeTangentPlanes(MRI_SURFACE *mris)
{
  VECTOR  *v_n, *v_e1, *v_e2, *v ;
  int     vno ;
  VERTEX  *vertex ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    /* now find some other non-parallel vector */
#if 0
    if (!FZERO(vertex->nx) || !FZERO(vertex->ny))
    {VECTOR_LOAD(v, 0.0, 0.0, 1.0) ; }
    else
    {VECTOR_LOAD(v, 0.0, 1.0, 0.0) ; }
#else
    VECTOR_LOAD(v, vertex->ny, vertex->nz, vertex->nx) ;
#endif
    V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    if (FZERO(V3_LEN(v_e1)))  /* happened to pick a parallel vector */
    {
      VECTOR_LOAD(v, vertex->ny, -vertex->nz, vertex->nx) ;
      V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    }

    if (FZERO(V3_LEN(v_e1)) && DIAG_VERBOSE_ON)  /* happened to pick a parallel vector */
      fprintf(stderr, "vertex %d: degenerate tangent plane\n", vno) ;
    V3_CROSS_PRODUCT(v_n, v_e1, v_e2) ;
    V3_NORMALIZE(v_e1, v_e1) ;
    V3_NORMALIZE(v_e2, v_e2) ;
    vertex->e1x = V3_X(v_e1) ; vertex->e2x = V3_X(v_e2) ;
    vertex->e1y = V3_Y(v_e1) ; vertex->e2y = V3_Y(v_e2) ;
    vertex->e1z = V3_Z(v_e1) ; vertex->e2z = V3_Z(v_e2) ;
  }

  VectorFree(&v) ;
  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/

#define ILL_CONDITIONED   500000.0
int
MRIScomputeSecondFundamentalForm(MRI_SURFACE *mris)
{
  int    vno, i, n, vmax, nbad = 0 ;
  VERTEX *vertex, *vnb ;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse, *m_eigen, *m_Q ;
  VECTOR *v_c, *v_z, *v_n, *v_e1, *v_e2, *v_yi ;
  float  k1, k2, evalues[3], a11, a12, a21, a22, cond_no, kmax, kmin, rsq, k ;
  double ui, vi, total_area = 0.0, max_error ;

  if (mris->status == MRIS_PLANE)
    return(NO_ERROR) ;

  mrisComputeTangentPlanes(mris) ;

  v_c = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_yi = VectorAlloc(3, MATRIX_REAL) ;
  m_Q = MatrixAlloc(2, 2, MATRIX_REAL) ;   /* the quadratic form */
  m_eigen = MatrixAlloc(2, 2, MATRIX_REAL) ;

  mris->Kmin = mris->Hmin = 10000.0f ;
  mris->Kmax = mris->Hmax = -10000.0f ;
  mris->Ktotal = 0.0f ;
  vmax = -1 ; max_error = -1.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    vertex = &mris->vertices[vno] ;
    if (vertex->ripflag)
      continue ;

    if (vno == 142915)
      DiagBreak() ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z) ;
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z) ;

    if (vertex->vtotal <= 0)
      continue ;
    
    m_U = MatrixAlloc(vertex->vtotal, 3, MATRIX_REAL) ;
    v_z = VectorAlloc(vertex->vtotal, MATRIX_REAL) ;

    if (vno == Gdiag_no)
      DiagBreak() ;

    /* fit a quadratic form to the surface at this vertex */
    kmin = 10000.0f ; kmax = -kmin ;
    for (n = i = 0 ; i < vertex->vtotal ; i++)
    {
      vnb = &mris->vertices[vertex->v[i]] ;
      if (vnb->ripflag)
        continue ;
/* 
   calculate the projection of this vertex onto the local tangent plane 
*/
      VECTOR_LOAD(v_yi, vnb->x-vertex->x, vnb->y-vertex->y,vnb->z-vertex->z);
      ui = V3_DOT(v_yi, v_e1) ; vi = V3_DOT(v_yi, v_e2) ;
      *MATRIX_RELT(m_U, n+1, 1) = ui*ui ;
      *MATRIX_RELT(m_U, n+1, 2) = 2*ui*vi ;
      *MATRIX_RELT(m_U, n+1, 3) = vi*vi ;
      VECTOR_ELT(v_z, n+1) = V3_DOT(v_n, v_yi) ;  /* height above TpS */
      rsq = ui*ui + vi*vi ;
      if (!FZERO(rsq))
      {
        k = VECTOR_ELT(v_z, n+1) / rsq ;
        if (k > kmax)
          kmax = k ;
        if (k < kmin)
          kmin = k ;
      }
      n++ ;
    }

    m_Ut = MatrixTranspose(m_U, NULL) ;          /* Ut */
    m_tmp2 = MatrixMultiply(m_Ut, m_U, NULL) ;   /* Ut U */
    cond_no = MatrixConditionNumber(m_tmp2) ;
#if 0
    m_inverse = MatrixInverse(m_tmp2, NULL) ;    /* (Ut U)^-1 */
#else
    m_inverse = MatrixSVDInverse(m_tmp2, NULL) ;    /* (Ut U)^-1 */
#endif
    if (!m_inverse)   /* singular matrix - must be planar?? */
    {
      nbad++ ;
      evalues[0] = evalues[1] = 0.0 ;
    }
    else
    {
      m_tmp1 = MatrixMultiply(m_Ut, v_z, NULL) ;   /* Ut z */
      MatrixMultiply(m_inverse, m_tmp1, v_c) ;     /* (Ut U)^-1 Ut z */

      /* now build Hessian matrix */
      *MATRIX_RELT(m_Q,1,1) = 2*VECTOR_ELT(v_c, 1) ;
      *MATRIX_RELT(m_Q,1,2) = *MATRIX_RELT(m_Q,2,1) = 2*VECTOR_ELT(v_c, 2) ;
      *MATRIX_RELT(m_Q,2,2) = 2*VECTOR_ELT(v_c, 3) ;

      if (cond_no >= ILL_CONDITIONED)
      {
#if 0
        MatrixSVDEigenValues(m_Q, evalues) ;
        vertex->k1 = k1 = evalues[0] ;
        vertex->k2 = k2 = evalues[1] ;
#else
        vertex->k1 = k1 = kmax ;
        vertex->k2 = k2 = kmin ;
#endif
        vertex->K = k1*k2 ; vertex->H = (k1+k2)/2 ;
        MatrixFree(&m_Ut) ;
        MatrixFree(&m_tmp2) ;
        MatrixFree(&m_U) ;
        VectorFree(&v_z) ;
        MatrixFree(&m_tmp1) ;
        MatrixFree(&m_inverse) ;
        continue ;
      }

      /* the columns of m_eigen will be the eigenvectors of m_Q */
      if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL)
      {
        nbad++ ;
        MatrixSVDEigenValues(m_Q, evalues) ;
        vertex->k1 = k1 = evalues[0] ;
        vertex->k2 = k2 = evalues[1] ;
        vertex->K = k1*k2 ; vertex->H = (k1+k2)/2 ;
        MatrixFree(&m_Ut) ;
        MatrixFree(&m_tmp2) ;
        MatrixFree(&m_U) ;
        VectorFree(&v_z) ;
        MatrixFree(&m_tmp1) ;
        MatrixFree(&m_inverse) ;
        continue ;
      }

      MatrixFree(&m_tmp1) ;
      MatrixFree(&m_inverse) ;
    }
    k1 = evalues[0] ; k2 = evalues[1] ;
    vertex->k1 = k1 ; vertex->k2 = k2 ;
    vertex->K = k1 * k2 ;
    vertex->H = (k1 + k2) / 2 ;
    if (vno == Gdiag_no && (Gdiag & DIAG_SHOW))
      fprintf(stdout, "v %d: k1=%2.3f, k2=%2.3f, K=%2.3f, H=%2.3f\n",
              vno, vertex->k1, vertex->k2, vertex->K, vertex->H) ;
    if (vertex->K < mris->Kmin)
      mris->Kmin = vertex->K ;
    if (vertex->H < mris->Hmin)
      mris->Hmin = vertex->H ;
    if (vertex->K > mris->Kmax)
      mris->Kmax = vertex->K ;
    if (vertex->H > mris->Hmax)
      mris->Hmax = vertex->H ;
    mris->Ktotal += (double)k1 * (double)k2 * (double)vertex->area ;
    total_area += (double)vertex->area ;

    /* now update the basis vectors to be the principal directions */
    a11 = *MATRIX_RELT(m_eigen,1,1) ; a12 = *MATRIX_RELT(m_eigen,1,2) ;
    a21 = *MATRIX_RELT(m_eigen,2,1) ; a22 = *MATRIX_RELT(m_eigen,2,2) ;
    vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a21 ;
    vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a21 ;
    vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a21 ;
    vertex->e2x = V3_X(v_e1) * a12 + V3_X(v_e2) * a22 ;
    vertex->e2y = V3_Y(v_e1) * a12 + V3_Y(v_e2) * a22 ;
    vertex->e2z = V3_Z(v_e1) * a12 + V3_Z(v_e2) * a22 ;

    MatrixFree(&m_Ut) ;
    MatrixFree(&m_tmp2) ;
    MatrixFree(&m_U) ;
    VectorFree(&v_z) ;

  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "max H error=%2.3f at %d\n", max_error, vmax) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout, "total area = %2.3f\n", total_area);

  if (Gdiag & DIAG_SHOW && (nbad > 0))
    fprintf(stdout, "%d ill-conditioned points\n", nbad) ;
  MatrixFree(&m_eigen) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_c) ;
  VectorFree(&v_n) ;
  VectorFree(&v_yi) ;
  MatrixFree(&m_Q) ;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
int
MRISsmoothValues(MRIS *mris, int niter)
{
  char funcname[STRLEN]="smoothValues";
  int iter,k,m,n;
  VERTEX *v;
  float sum;

  if(niter<=0) return(-1);

  MsgPrintf("%s(%d): smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++)
  {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].tdx = mris->vertices[k].val;
    for (k=0;k<mris->nvertices;k++) {
      v = &mris->vertices[k];
      sum=v->tdx;
      if (k == Gdiag_no) DiagBreak();
      n = 1;
      for (m=0;m<v->vnum;m++) {
        sum += mris->vertices[v->v[m]].tdx;
        n++;
      }
      if (!finite(sum)) DiagBreak();
      if (n>0) 
        v->val = sum/n;
      if (!finite(v->val))
        DiagBreak() ;
    }
  }
  MsgPrintf("\n");

  return(0);
}

/*---------------------------------------------------------------------------*/
int
MRISsmoothValuesHeat(MRIS *mris, int niter, float sigma)
{
  char funcname[STRLEN]="smoothValues";
  int iter,k,m,n;
  VERTEX *v;
  float denomW,dist;

  if(niter<=0 || sigma<=0) return(-1);

  /* calculate weights */
  MsgPrintf("%s: calculating weights for heat kernel smoothing\n",funcname);
  for (k=0;k<mris->nvertices;k++) {
    v = &mris->vertices[k];
    n = v->vnum;
    denomW = 1;
    for (m=0;m<n;m++) {
      dist = v->dist[m];
      v->dist_orig[m] = exp(-0.5*dist*dist/(sigma*sigma));
      denomW += v->dist_orig[m];
    }
    v->d = 1/denomW;
    for (m=0;m<n;m++)
      v->dist_orig[m] /= denomW;
  }

  MsgPrintf("%s(%d): heat kernel smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++)
  {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].tdx = mris->vertices[k].val;
    for (k=0;k<mris->nvertices;k++) {
      v = &mris->vertices[k];
      n = v->vnum;
      v->val = v->d*v->tdx;
      for (m=0;m<n;m++)
        v->val += v->dist_orig[m]*mris->vertices[v->v[m]].tdx;
    }
  }
  MsgPrintf("\n");

  return(0);
}

/*---------------------------------------------------------------------------*/
int
MRISsmoothComplexValues(MRIS *mris, int niter)
{
  char funcname[STRLEN]="smoothComplexValues";
  int iter,k,m,n;
  VERTEX *v;
  float sum,imag_sum;

  if(niter<=0) return(-1);

  MsgPrintf("%s(%d): smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++)
  {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++) {
      mris->vertices[k].tdx = mris->vertices[k].val;
      mris->vertices[k].tdy = mris->vertices[k].imag_val;
    }
    for (k=0;k<mris->nvertices;k++) {
      v = &mris->vertices[k];
      sum=v->tdx;
      imag_sum=v->tdy;
      if (k == Gdiag_no) DiagBreak();
      n = 1;
      for (m=0;m<v->vnum;m++) {
        sum += mris->vertices[v->v[m]].tdx;
        imag_sum += mris->vertices[v->v[m]].tdy;
        n++;
      }
      if (!finite(sum)) DiagBreak();
      if (n>0) {
        v->val = sum/n;
        v->imag_val = imag_sum/n;
      }
      if (!finite(v->val) || !finite(v->imag_val))
        DiagBreak() ;
    }
  }
  MsgPrintf("\n");

  return(0);
}

/*---------------------------------------------------------------------------*/
int
MRISsmoothValuesSparse(MRIS *mris, int niter)
{
  char funcname[STRLEN]="smoothValuesSparse";
  int    iter,k,m,n;
  VERTEX *v;
  float  sum;
  
  if(niter<=0) return(-1);

  for (k=0;k<mris->nvertices;k++) {
    mris->vertices[k].fixedval=!mris->vertices[k].undefval;
  }

  MsgPrintf("%s(%d): smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++) {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++) {
      mris->vertices[k].tdx = mris->vertices[k].val;
      mris->vertices[k].old_undefval = mris->vertices[k].undefval;
    }
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].fixedval) {
        v = &mris->vertices[k];
        if(v->undefval) n=0; else n=1;
        sum = v->tdx;
        for (m=0;m<v->vnum;m++) {
          if (!mris->vertices[v->v[m]].old_undefval) {
            sum += mris->vertices[v->v[m]].tdx;
            n++;
          }
        }
        if (n>0) {
          v->val = sum/n;
          v->undefval = FALSE;
        }
        if (!finite(v->val))
          DiagBreak() ;
      } /* if !fixedval */
    } /* for k */
  } /* for iter */
  MsgPrintf("\n");

  return(0);
}

/*---------------------------------------------------------------------------*/
int
MRISsmoothComplexValuesSparse(MRIS *mris, int niter)
{
  char funcname[STRLEN]="smoothComplexValuesSparse";
  int iter,k,m,n;
  VERTEX *v;
  float sum,imag_sum;

  if(niter<=0) return(-1);

  MsgPrintf("%s(%d): smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++)
  {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].fixedval &&
          !mris->vertices[k].undefval) {
        mris->vertices[k].tdx = mris->vertices[k].val;
        mris->vertices[k].tdy = mris->vertices[k].imag_val;
      }
    }
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].fixedval) {
        v = &mris->vertices[k];
        if(v->undefval) n=0; else n=1;
        sum = v->tdx;
        imag_sum=v->tdy;
        for (m=0;m<v->vnum;m++) {
          if (!mris->vertices[v->v[m]].undefval) {
            sum += mris->vertices[v->v[m]].tdx;
            imag_sum += mris->vertices[v->v[m]].tdy;
            n++;
          }
        }
        if (!finite(sum)) DiagBreak();
        if (n>0) {
          v->val = sum/n;
          v->imag_val = imag_sum/n;
          v->undefval = FALSE;
        }
        if (!finite(v->val) || !finite(v->imag_val))
          DiagBreak() ;
      }
    }
  }
  MsgPrintf("\n");

  return(0);
}

/*---------------------------------------------------------------------------*/
int
MRISsmoothValuesROI(MRIS *mris, int niter)
{
  char funcname[STRLEN]="smoothValuesROI";
  int    iter,k,m,n;
  VERTEX *v;
  float  sum;
  
  if(niter<=0) return(-1);

  MsgPrintf("%s(%d): smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++) {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].undefval)
        mris->vertices[k].tdx = mris->vertices[k].val;
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].undefval) {
        v = &mris->vertices[k];
        n = 1;
        sum = v->tdx;
        for (m=0;m<v->vnum;m++) {
          if (!mris->vertices[v->v[m]].undefval) {
            sum += mris->vertices[v->v[m]].tdx;
            n++;
          }
        }
        if (n>0) {
          v->val = sum/n;
          v->undefval = FALSE;
        }
        if (!finite(v->val))
          DiagBreak() ;
      } /* if ripflag, etc. */
    } /* for k */
  } /* for iter */
  MsgPrintf("\n");

  return(0);
}

/*---------------------------------------------------------------------------*/
int
MRISsmoothComplexValuesROI(MRIS *mris, int niter)
{
  char funcname[STRLEN]="smoothComplexValuesROI";
  int iter,k,m,n;
  VERTEX *v;
  float sum,imag_sum;

  if(niter<=0) return(-1);

  MsgPrintf("%s(%d): smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++)
  {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].undefval) {
        mris->vertices[k].tdx = mris->vertices[k].val;
        mris->vertices[k].tdy = mris->vertices[k].imag_val;
      }
    }
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].undefval) {
        v = &mris->vertices[k];
        sum = v->tdx;
        imag_sum=v->tdy;
        if (k == Gdiag_no) DiagBreak();
        n = 1;
        for (m=0;m<v->vnum;m++) {
          if (!mris->vertices[v->v[m]].undefval) {
            sum += mris->vertices[v->v[m]].tdx;
            imag_sum += mris->vertices[v->v[m]].tdy;
            n++;
          }
        }
        if (!finite(sum)) DiagBreak();
        if (n>0) {
          v->val = sum/n;
          v->imag_val = imag_sum/n;
          v->undefval = FALSE;
        }
        if (!finite(v->val) || !finite(v->imag_val))
          DiagBreak() ;
      }
    }
  }
  MsgPrintf("\n");

  return(0);
}

/*---------------------------------------------------------------------------*/
int
MRISsmoothValuesHeatROI(MRIS *mris, int niter, float sigma)
{
  char funcname[STRLEN]="smoothValuesHeatROI";
  int iter,k,m,n;
  VERTEX *v;
  float denomW,dist;

  if(niter<=0 || sigma<=0) return(-1);

  /* calculate weights */
  MsgPrintf("%s: calculating weights for heat kernel smoothing\n",funcname);
  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].undefval) {
      v = &mris->vertices[k];
      n = v->vnum;
      denomW = 1;
      for (m=0;m<n;m++) {
        dist = v->dist[m];
        v->dist_orig[m] = exp(-0.5*dist*dist/(sigma*sigma));
        denomW += v->dist_orig[m];
      }
      v->d = 1/denomW;
      for (m=0;m<n;m++)
        v->dist_orig[m] /= denomW;
    }
  }

  MsgPrintf("%s(%d): heat kernel smoothing",funcname,niter);
  for (iter=0;iter<niter;iter++)
  {
    MsgPrintf(".");
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].undefval)
      mris->vertices[k].tdx = mris->vertices[k].val;
    for (k=0;k<mris->nvertices;k++) {
      if (!mris->vertices[k].undefval) {
        v = &mris->vertices[k];
        n = v->vnum;
        v->val = v->d*v->tdx;
        for (m=0;m<n;m++)
          v->val += v->dist_orig[m]*mris->vertices[v->v[m]].tdx;
      }
    }
  }
  MsgPrintf("\n");

  return(0);
}


/*---------------------------------------------------------------------------*/
int
nverticesSurf(char *subject, char *hemi)
{
  /* 
     this function uses older style surface reading and maybe not be
     compatible with some surfaces - DH 

     this function doesn't actually read the whole surface, it just
     reads the number of vertices
  */
  char funcname[STRLEN]="nverticesSurf";
  int version;
  FILE *fp;
  int magic,nvertices,nfaces;
  char line[512];
  char *SUBJECTS_DIR = NULL;
  char *MRI_DIR = NULL;
  char fname[STRLEN];
  int triangle_file = 0;
  int vperface = 4;

  if(MATCH(subject,"ico")) {
    MRI_DIR = getenv("MRI_DIR");
    if (MRI_DIR==NULL) {
      MsgPrintf("%s: environment variable MRI_DIR undefined (use setenv)\n",\
             funcname);
      return(-1);
    }
    sprintf(fname,"%s/lib/bem/ic%d.tri",MRI_DIR,7);
    MsgPrintf("%s: reading from %s...\n",funcname,fname);

    fp = fopen(fname,"r");
    if(fp==NULL){
      MsgPrintf("%s: ### File %s not found\n",funcname,fname);
      return(-1);
    }
    fgetl(line, 150, fp);
    sscanf(line, "%d", &nvertices);
    MsgPrintf("%s: vertices=%d\n",funcname,nvertices);
  } else {
    MsgPrintf("%s(%s,%s)\n",funcname,subject,hemi);
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL) {
      MsgPrintf("%s: environment variable SUBJECTS_DIR undefined (use setenv)\n",\
             funcname);
      return(-1);
    }
    sprintf(fname,"%s/%s/surf/%s.orig",SUBJECTS_DIR,subject,hemi);
    MsgPrintf("%s: reading from %s...\n",funcname,fname);

    fp = fopen(fname,"rb");
    if(fp==NULL){
      MsgPrintf("%s: ### File %s not found\n",funcname,fname);
      return(-1);
    }
    fread3(&magic,fp);
    switch(magic) {
      case QUAD_FILE_MAGIC_NUMBER:
        MsgPrintf("%s: QUAD_FILE_MAGIC_NUMBER\n", funcname);
        version = -1;
        break;
      case TRIANGLE_FILE_MAGIC_NUMBER:
        triangle_file = 1;
        version = -3;
        vperface = 3;
        MsgPrintf("%s: TRIANGLE_FILE_MAGIC_NUMBER\n", funcname);
        break;
      case NEW_QUAD_FILE_MAGIC_NUMBER:
        MsgPrintf("%s: NEW_QUAD_FILE_MAGIC_NUMBER\n", funcname);
        version = -2;
        break;
      default:
        MsgPrintf("%s: old surface format [%d]\n", funcname,magic);
        version = 0;
        rewind(fp);
    }
    if(triangle_file) {
      fgets(line, 200, fp);
      fscanf(fp, "\n");
      nvertices = freadInt(fp);
      nfaces = freadInt(fp);
    } else {
      fread3(&nvertices,fp);
      fread3(&nfaces,fp);
    }
    MsgPrintf("%s: vertices=%d, faces=%d\n",funcname,nvertices,nfaces);
  }
  fclose(fp) ;
  return(nvertices);
}


#if 0
static ICOSOHEDRON *
read_icosahedron(char *fname)
{
  FILE        *fp ;
  char        line[200], *cp ;
  int         vno, fno, vno1, vno2, vno3, n, nvertices, nfaces ;
  float       x, y, z ;
  ICOSOHEDRON *ico ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NOFILE, "read_icosahedron: could not open %s", fname));

  ico = (ICOSOHEDRON *)calloc(1, sizeof(ICOSOHEDRON)) ;

  fgetl(line, 150, fp) ;   /* discard # of vertices */
  sscanf(line, "%d", &nvertices) ;
  ico->nvertices = nvertices ;
  ico->vertices = 
    (IC_VERTEX *)calloc(nvertices, sizeof(IC_VERTEX)) ;
  if (!ico->vertices)
    ErrorExit(ERROR_NOMEMORY, "read_ico: could not allocate vertex list") ;


  /* first read vertices */
  n = 0 ;
  while ((cp = fgetl(line, 150, fp)) != NULL)
  {
    if (sscanf(cp, "%d %f %f %f\n", &vno, &x, &y, &z) < 4)
      break ;
    ico->vertices[vno-1].x = x ;
    ico->vertices[vno-1].y = y ;
    ico->vertices[vno-1].z = z ;
    if (++n >= ico->nvertices)
      break ;
  }
  n = 0 ; fgetl(line, 150, fp) ;   /* discard # of faces */
  sscanf(line, "%d", &nfaces) ;
  ico->nfaces = nfaces ;
  ico->faces = 
    (IC_FACE *)calloc(ico->nfaces, sizeof(IC_FACE)) ;
  if (!ico->faces)
    ErrorExit(ERROR_NOMEMORY, "read_ico: could not allocate vertex list") ;
  while ((cp = fgetl(line, 150, fp)) != NULL)
  {
    if (sscanf(cp, "%d %d %d %d\n", &fno, &vno1, &vno2, &vno3) < 4)
      break ;
    ico->faces[fno-1].vno[0] = vno1 ;
    ico->faces[fno-1].vno[1] = vno2 ;
    ico->faces[fno-1].vno[2] = vno3 ;
    if (++n >= ico->nfaces)
      break ;
  }
  fclose(fp) ;
  return(ico) ;
}
#endif

/*---------------------------------------------------------------------------*/
int
maxVert(char *fname)
{
  char funcname[STRLEN]="maxVert";
  int i,k,lat,num,maxk=0;
  float f;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  fread2(&lat,fp);
  fread3(&num,fp);
  MsgPrintf("%s: number of values: %d\n",funcname,num);
  for (i=0;i<num;i++) {
    fread3(&k,fp);
    fread4(&f,fp);
    if (k>maxk) maxk=k;
  }
  fclose(fp);
  MsgPrintf("%s: max vertex index = %d\n",funcname,maxk);
  return(maxk);
}

/*---------------------------------------------------------------------------*/
int
readSurfVals(char *fname, float *data, int nvertices)
{
  /*
     this function is a method for reading a value file independent of 
     the MRIS structure 

     nvertices must be supplied to make sure there are no accidental
     memory overruns
  */
  
  char funcname[STRLEN]="readSurfVals";
  int i,k,lat,num,maxk=-1;
  float f;
  FILE *fp;
  double sum=0,sum2=0,max= -1000,min=1000;

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  fread2(&lat,fp);
  fread3(&num,fp);
  if (num > nvertices) {
    MsgPrintf("%s: number of values (%d) > number of vertices (%d)\n",
           funcname,num,nvertices);
    return(-1);
  }
  for (i=0;i<nvertices;i++) data[i]=0.0;
  for (i=0;i<num;i++) {
    fread3(&k,fp);
    fread4(&f,fp);
    if (k>maxk) maxk=k;
    if (k>=nvertices||k<0) {
      MsgPrintf("%s: vertex index (%d) out of range; val=%f\n",funcname,k,f);
      return(-1);
    } else {
      data[k] = f;
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);
  sum /= num;
  sum2 = sqrt(sum2/num-sum*sum);
  MsgPrintf("%s: file %s read \n",funcname,fname);
  MsgPrintf("      with %d vals and max vertex index = %d\n",num,maxk);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(0);
}

/*---------------------------------------------------------------------------*/
int
writeSurfVals(char *fname, float *data, int nvertices)
{
  char funcname[STRLEN]="writeSurfVals";
  FILE *fp;
  float f;
  int k,num,maxk=-1;
  double sum=0,sum2=0,max= -1000,min=1000;

  fp = fopen(fname,"w");
  if(fp==NULL){
    MsgPrintf("%s: ### can't create file %s\n",funcname,fname);
    return(-1);
  }
  for (k=0,num=0;k<nvertices;k++) 
    if (data[k]!=0) num++;
  fwrite2(0,fp);  /* latency (not used) */
  fwrite3(num,fp);
  for (k=0;k<nvertices;k++) {
    if (data[k]!=0) {
      if (k>maxk) maxk=k;
      fwrite3(k,fp);
      f = data[k];
      fwrite4(f,fp);
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);
  sum /= num;
  sum2 = sqrt(sum2/num-sum*sum);
  MsgPrintf("%s: file %s written \n",funcname,fname);
  MsgPrintf("      with %d vals and max vertex index = %d\n",num,maxk);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(0);
}

/*---------------------------------------------------------------------------*/
int
writeTextVals(char *fname, float *data, int nvertices)
{
  char funcname[STRLEN]="writeTextVals";
  FILE *fp;
  float f;
  int k,num,maxk=-1;
  double sum=0,sum2=0,max= -1000,min=1000;

  fp = fopen(fname,"w");
  if(fp==NULL){
    MsgPrintf("%s: ### can't create file %s\n",funcname,fname);
    return(-1);
  }
  for (k=0,num=0;k<nvertices;k++) 
    if (data[k]!=0) num++;
  fprintf(fp,"0\n");
  fprintf(fp,"%d\n",num);
  for (k=0;k<nvertices;k++) {
    if (data[k]!=0) {
      if (k>maxk) maxk=k;
      fprintf(fp,"%d\n",k);
      f = data[k];
      fprintf(fp,"%f\n",f);
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);
  sum /= num;
  sum2 = sqrt(sum2/num-sum*sum);
  MsgPrintf("%s: file %s written \n",funcname,fname);
  MsgPrintf("      with %d vals and max vertex index = %d\n",num,maxk);
  MsgPrintf("      avg=%f, stdev=%f, min=%f, max=%f\n",sum,sum2,min,max);

  return(0);
}

/*---------------------------------------------------------------------------*/
MRI_SURFACE *
openSurface(char *subject, char *hemi, char *surface)
{
  char *funcname="open_surface";
  MRI_SURFACE  *mris;
  char fname[STRLEN];
  char *SUBJECTS_DIR = NULL;
  char *MRI_DIR = NULL;
  int IcoOrder = 7;
  float IcoRadius = 100.0;

  if(MATCH(subject,"ico")) {
    MRI_DIR = getenv("MRI_DIR");
    if (MRI_DIR==NULL)
      ErrorExit(ERROR_BADPARM,"%s: ### environment variable MRI_DIR undefined (use setenv)\n",
        funcname);
    sprintf(fname,"%s/lib/bem/ic%d.tri",MRI_DIR,IcoOrder);
    MsgPrintf("%s: %s...\n",funcname,fname);
    mris = ReadIcoByOrder(IcoOrder, IcoRadius);
  } else {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR==NULL)
      ErrorExit(ERROR_BADPARM,"%s: ### environment variable SUBJECTS_DIR undefined (use setenv)\n",
        funcname);
    sprintf(fname,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surface);
    MsgPrintf("%s: %s...\n",funcname,fname);
    mris = MRISread(fname);
  }
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", funcname, fname);
  return(mris);
}

/*---------------------------------------------------------------*/
int GetNVtxsFromWFile(char *wfile)
{
  FILE *fp;
  int i,ilat, num, nvertices;
  int *vtxnum;
  float *wval;
  char funcname[STRLEN]="GetNVtxsFromWFile";

  fp = fopen(wfile,"r");
  if (fp==NULL) {
    ErrorExit(ERROR_BADFILE,"%s: ### GetNVtxsFromWFile(): Could not open %s\n",
      funcname,wfile);
  }
  
  fread2(&ilat,fp);
  fread3(&num,fp);
  vtxnum = (int *)   calloc(sizeof(int),   num);
  wval   = (float *) calloc(sizeof(float), num);

  for (i=0;i<num;i++){
    fread3(&vtxnum[i],fp);
    wval[i] = freadFloat(fp) ;
  }
  fclose(fp);

  nvertices = vtxnum[num-1] + 1;

  free(vtxnum);
  free(wval);

  return(nvertices);
}

/*---------------------------------------------------------------------------*/
int
MRISclearValues(MRI_SURFACE *mris)
{
  int k;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].val = 0;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISclearAllValues(MRI_SURFACE *mris)
{
  int k;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].val = 0;
    mris->vertices[k].imag_val = 0;
    mris->vertices[k].valbak = 0;
    mris->vertices[k].val2 = 0;
    mris->vertices[k].val2bak = 0;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
MRISclearMarks(MRI_SURFACE *mris)
{
  int k;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].marked = 0;
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
float
MRISgaussFWHM(MRIS *mris, int sparseflag) {
  int i, j, k, dfcount=0, vcount=0;
  double sum=0, sumsq=0;
  double df=0, dfsum=0, dfsumsq=0;
  double var=0, dfvar=0, varratio;
  float avgdist=0, fwhm, dist=0;

  char funcname[STRLEN]="MRISgaussFWHM";

  /* calculate sums and sums of squares */
  for (i=0;i<mris->nvertices;i++) {
    if(sparseflag && mris->vertices[i].undefval) continue;
    vcount++;
    sum   += mris->vertices[i].val;
    sumsq += mris->vertices[i].val * mris->vertices[i].val;

    /* estimage partial derivatives between neighboring vertices */
    for (j=0;j<mris->vertices[i].vnum;j++) {
      k = mris->vertices[i].v[j];
      if(sparseflag && mris->vertices[k].undefval) continue;
      if(k<=i) continue; /* only count neighbor pair once */
      dist = mris->vertices[i].dist[j];
      if(dist<=0) continue;
      avgdist += dist;
      df = mris->vertices[k].val - mris->vertices[i].val;
      dfsum   += df;
      dfsumsq += df * df;
      dfcount++;
    }
  }
  avgdist /= dfcount;

  /* estimate variance of data across all vertices */
  var =  (sumsq - (sum*sum)/vcount)/(vcount-1);
  
  /* estimate variance of partial derivatives */
  dfvar = (dfsumsq - (dfsum*dfsum)/dfcount)/(dfcount-1);
  
  /* estimate equivalent Gaussian filter width */
  varratio = -log(1.0 - 0.5*(dfvar/var));
  if (varratio<=0)
    varratio = 0.5*(dfvar/var);
  if (varratio<=0)
    fwhm=0;
  else {
    fwhm = sqrt(2.0*log(2.0)/varratio)*avgdist;
  }

  MsgPrintf("%s: overall variance = %0.5lf\n",funcname,var);
  MsgPrintf("%s: neighbor variance = %0.5lf\n",funcname,dfvar);

  return fwhm;
}

/*---------------------------------------------------------------------------*/
float
MRISgaussCxFWHM(MRIS *mris, int sparseflag) {
  int i, j, k, dfcount=0, vcount=0;
  double r_sum=0, r_sumsq=0;
  double i_sum=0, i_sumsq=0;
  double df_r=0, df_r_sum=0, df_r_sumsq=0;
  double df_i=0, df_i_sum=0, df_i_sumsq=0;
  double var=0, dfvar=0, varratio;
  float avgdist=0, fwhm, dist=0;

  char funcname[STRLEN]="MRISgaussCxFWHM";

  /* calculate sums and sums of squares */
  for (i=0;i<mris->nvertices;i++) {
    if(sparseflag && mris->vertices[i].undefval) continue;
    vcount++;
    r_sum   += mris->vertices[i].val;
    r_sumsq += mris->vertices[i].val * mris->vertices[i].val;
    i_sum   += mris->vertices[i].imag_val;
    i_sumsq += mris->vertices[i].imag_val * mris->vertices[i].imag_val;

    /* estimage partial derivatives between neighboring vertices */
    for (j=0;j<mris->vertices[i].vnum;j++) {
      k = mris->vertices[i].v[j];
      if(sparseflag && mris->vertices[k].undefval) continue;
      if(k<=i) continue; /* only count neighbor pair once */
      dist = mris->vertices[i].dist[j];
      if(dist<=0) continue;
      avgdist += dist;
      df_r = (mris->vertices[k].val - mris->vertices[i].val)/dist;
      df_i = (mris->vertices[k].imag_val - mris->vertices[i].imag_val)/dist;
      df_r_sum   += df_r;
      df_r_sumsq += df_r * df_r;
      df_i_sum   += df_i;
      df_i_sumsq += df_i * df_i;
      dfcount++;
    }
  }
  avgdist /= dfcount;

  /* estimate variance of data across all vertices */
  var =  (r_sumsq - (r_sum*r_sum)/vcount)/(vcount-1) +
         (i_sumsq - (i_sum*i_sum)/vcount)/(vcount-1);
  
  /* estimate variance of partial derivatives */
  dfvar = (df_r_sumsq - (df_r_sum*df_r_sum)/dfcount)/(dfcount-1) +
          (df_i_sumsq - (df_i_sum*df_i_sum)/dfcount)/(dfcount-1);
  
  /* estimate equivalent Gaussian filter width */
  varratio = 1.0 - 0.5*(dfvar/var);
  if ((varratio<=0) || (dfvar<=0))
    fwhm=0;
  else {
    fwhm = sqrt(-2.0*log(2.0)/log(varratio))*avgdist;
  }

  MsgPrintf("%s: overal variance = %0.5lf\n",funcname,var);
  MsgPrintf("%s: neighbor variance = %0.5lf\n",funcname,dfvar);

  return fwhm;
}

/*---------------------------------------------------------------------------*/
/* transform */
/*---------------------------------------------------------------------------*/
int
mrisReadTransform(MRIS *mris, char *mris_fname)
{
  char funcname[STRLEN]="mrisReadTransform";
  char transform_fname[STRLEN], fpref[300] ;

  mris->transform_loaded = 0;
  mris->inverse_transform_loaded = 0;
   if(mris->linear_transform != NULL)
    free(mris->linear_transform);
  if(mris->inverse_linear_transform != NULL)
    free(mris->inverse_linear_transform);
  mris->linear_transform = (Transform *)calloc(1,sizeof(Transform));
  mris->inverse_linear_transform = (Transform *)calloc(1,sizeof(Transform));
  if(mris->linear_transform==NULL ||
     mris->inverse_linear_transform==NULL) {
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "%s: could not allocate memory for xform",funcname)) ;
  }

  FileNamePath(mris_fname, fpref) ;
  sprintf(transform_fname, "%s/../mri/transforms/talairach.xfm", fpref) ;
  if (!FileExists(transform_fname))
    return(ERROR_NOFILE) ;
  if (readTransform(transform_fname, mris->linear_transform)==-1) {
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "%s: could not read xform file '%s'",
                 funcname, transform_fname)) ;
  } else {
    mris->transform_loaded = 1;
    if (invertTransform(mris->linear_transform,
                        mris->inverse_linear_transform)==-1) {
      ErrorReturn(ERROR_BADFILE,
                  (ERROR_BADFILE,
                   "%s: problem inverting transform",
                   funcname, transform_fname)) ;
    } else {
      mris->inverse_transform_loaded = 1;
      mris->free_transform = 1 ;
    }
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
int
readTransform(char *fname, Transform *transform)
{
  FILE *fp;
  char line[MAX_LINE_LENGTH];
  char str[STRLEN];
  char funcname[STRLEN]="readTransform";
  int i=-1,j;
  double row[4];

  fp = fopen(fname,"r");
  if (fp==NULL) {
    MsgPrintf("%s: ### File %s not found\n",funcname,fname);
    return(-1);
  }
  if(transform==NULL) {
    printf("%s: ### No memory allocated for transform matrix!\n",funcname);
    return(-1);
  }

  do {
    if(fgets(line,MAX_LINE_LENGTH-1,fp)==NULL) return(-1);
    sscanf(line,"%s",str);
    if (MATCH(str,TRANSFORM_STRING)) i++;
    else if(i>-1) {
      j=0;
      sscanf(line,"%lf %lf %lf %lf",
        &row[0],&row[1],&row[2],&row[3]);
      for(j=0;j<4;j++)
        transform->m[i][j]=row[j];
      i++;
    }
    if(i>2) break;
  } while(!feof(fp));
  fclose(fp);
  transform->m[3][3]=1.0;

  return(0);
}

/*---------------------------------------------------------------------------*/
int
invertTransform(Transform *transform, Transform *inverse)
{
  char funcname[STRLEN]="invertTransform";
  int i,j,k,p,v,tmp;
  int row[4];
  double d,val,best_val,m,scale_factor;
  double a[4][4],solution[4][4], s[4];

  if(transform==NULL) {
    printf("%s: Transform matrix is NULL!\n",funcname);
    return(-1);
  } else if(inverse==NULL) {
    printf("%s: Inverse transform matrix is NULL!\n",funcname);
    return(-1);
  }

  for(i=0;i<4;i++) {
    for(j=0;j<4;j++)
      inverse->m[i][j] = 0.0;
    inverse->m[i][i] = 1.0;
  }

  for(i=0;i<4;i++) {
    for(j=0;j<4;j++)
      a[i][j] = transform->m[i][j];
    for(v=0;v<4;v++)
      solution[i][v] = inverse->m[v][i];
  }

  for(i=0;i<4;i++)
    row[i] = i;

  for(i=0;i<4;i++) {
    s[i] = FABS(a[i][0]);
    for(j=1;j<4;j++) {
      if(FABS(a[i][j])>s[i])
        s[i] = FABS(a[i][j]);
    }
    if( s[i] == 0.0 ) return(-1);
  }

  for(i=0;i<3;i++) {
    p = i;
    best_val = a[row[i]][i] / s[row[i]];
    best_val = FABS(best_val);
    for(j=i+1;j<4;j++)
    {
      val = a[row[j]][i] / s[row[j]];
      val = FABS(val);
      if(val>best_val) {
        best_val = val;
        p = j;
      }
    }

    if(a[row[p]][i]==0.0) return(-1);

    if(i!=p) {
      tmp = row[i];
      row[i] = row[p];
      row[p] = tmp;
    }

    for(j=i+1;j<4;j++) {
      if(a[row[i]][i]==0.0) return(-1);
      m = a[row[j]][i] / a[row[i]][i];
      for(k=i+1;k<4;k++)
        a[row[j]][k] -= m * a[row[i]][k];
      for(v=0;v<4;v++)
        solution[row[j]][v] -= m * solution[row[i]][v];
    }
  }

  if(a[row[3]][3] == 0.0) return(-1);

  for(i=3;i>=0;--i) {
    for(j=i+1;j<4;j++) {
      scale_factor = a[row[i]][j];
      for(v=0;v<4;v++)
        solution[row[i]][v] -= scale_factor * solution[row[j]][v];
    }
    for(v=0;v<4;v++)
      solution[row[i]][v] /= a[row[i]][i];
  }


  for(i=0;i<4;i++) {
    for(v=0;v<4;v++)
      inverse->m[v][i] = solution[row[i]][v];
  }

  for(i=0;i<3;i++) {
    for(j=i+1;j<4;j++) {
      d = inverse->m[i][j];
      inverse->m[i][j] = inverse->m[j][i];
      inverse->m[j][i] = d;
    }
  }

  return(0);
}

/*---------------------------------------------------------------------------*/
int
transformPoint(Transform *transform, double x, double y, double z,
  double *x_trans, double *y_trans, double *z_trans)
{
  double w=1.0,w_trans=1.0;

  if(transform==NULL) return(-1);

  *x_trans = transform->m[0][0] * x +
             transform->m[0][1] * y +
             transform->m[0][2] * z +
             transform->m[0][3] * w ; 

  *y_trans = transform->m[1][0] * x +
             transform->m[1][1] * y +
             transform->m[1][2] * z +
             transform->m[1][3] * w ;

  *z_trans = transform->m[2][0] * x +
             transform->m[2][1] * y +
             transform->m[2][2] * z +
             transform->m[2][3] * w ; 

   w_trans = transform->m[3][0] * x +
             transform->m[3][1] * y +
             transform->m[3][2] * z +
             transform->m[3][3] * w ; 

  if(w_trans!=0.0 && w_trans!=1.0) {
    *x_trans /= w_trans;
    *y_trans /= w_trans;
    *z_trans /= w_trans;
  }

  return(0);
}

/*---------------------------------------------------------------------------*/
int
printTransform(Transform *transform, char *transname)
{
  char funcname[STRLEN]="printTransform";
  int i,j;

  if(quiet) return(0);

  if(transform==NULL) {
    MsgPrintf("%s: Transform matrix %s is NULL!\n",funcname, transname);
    return(-1);
  }

  printf("%s(%s):\n",funcname,transname);
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      printf("   %f",transform->m[i][j]);
    }
    printf("\n");
  }

  return(0);
}

/*---------------------------------------------------------------------------*/
void
multTransform(Transform *t, Transform *t1, Transform *t2)
{
  int i,j,k;
  double sum;

  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      sum = 0.0;
      for(k=0;k<4;k++) {
        sum += t2->m[i][k] * t1->m[k][j];
      }
      t->m[i][j] = sum;
    }
  }
}

/*---------------------------------------------------------------------------*/
int
selectTalairachPoint(MRIS *mris, int *vindex,
  double x_tal, double y_tal, double z_tal)
{
  char funcname[STRLEN]="selectTalairachPoint";
  double x,y,z;

  if(!mris->inverse_transform_loaded) {
    printf("%s: ### inverse transform not loaded\n",funcname);
    return(-1);
  }

  if(transformPoint(mris->inverse_linear_transform,
                    x_tal, y_tal, z_tal, &x, &y, &z)==-1) {
    printf("%s: ### error transforming point\n",funcname);
    return(-1);
  }

  MsgPrintf("%s: talairach coords:\n",funcname);
  MsgPrintf("x=%f y=%f z=%f\n",x_tal,y_tal,z_tal);
  MsgPrintf("%s: surface coords:\n",funcname);
  MsgPrintf("x=%f y=%f z=%f\n",x,y,z);
  
  selectVertexCoord(mris,vindex,x,y,z);

  return(0) ;
}

/*---------------------------------------------------------------------------*/
int
selectVertexCoord(MRIS *mris, int *vindex,
  double x, double y, double z)
{
  char funcname[STRLEN]="selectVertexCoord";
  double d,dmin,vx=0,vy=0,vz=0,tx=0,ty=0,tz=0;
  int k;

  dmin = 1e10;
  for (k=0;k<mris->nvertices;k++) {
    tx = mris->vertices[k].x ;
    ty = mris->vertices[k].y ;
    tz = mris->vertices[k].z ;
    d = SQR(tx-x)+SQR(ty-y)+SQR(tz-z);
    if (d<dmin) {
      vx=tx;
      vy=ty;
      vz=tz;
      dmin=d;
      *vindex=k;
    }
  }

  MsgPrintf("%s: vertex %d: dist = %f\n",funcname,*vindex,sqrt(dmin));
  MsgPrintf("%s: vertex coords:\n",funcname);
  MsgPrintf("x=%f y=%f z=%f\n",vx,vy,vz);

  if(dmin==1e10)
    return(-1);
  else
    return(0);
}

/*from mri.c */
int stuff_four_by_four(MATRIX *m, float m11, float m12, float m13, float m14, 
                                  float m21, float m22, float m23, float m24, 
                                  float m31, float m32, float m33, float m34, 
                                  float m41, float m42, float m43, float m44)
{

  if(m == NULL)
  {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "stuff_four_by_four(): matrix is NULL"));
  }

  if(m->rows != 4 || m->cols != 4)
  {
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "stuff_four_by_four(): matrix is not four-by-four"));
  }

  *MATRIX_RELT(m, 1, 1) = m11;  *MATRIX_RELT(m, 1, 2) = m12;  *MATRIX_RELT(m, 1, 3) = m13;  *MATRIX_RELT(m, 1, 4) = m14;
  *MATRIX_RELT(m, 2, 1) = m21;  *MATRIX_RELT(m, 2, 2) = m22;  *MATRIX_RELT(m, 2, 3) = m23;  *MATRIX_RELT(m, 2, 4) = m24;
  *MATRIX_RELT(m, 3, 1) = m31;  *MATRIX_RELT(m, 3, 2) = m32;  *MATRIX_RELT(m, 3, 3) = m33;  *MATRIX_RELT(m, 3, 4) = m34;
  *MATRIX_RELT(m, 4, 1) = m41;  *MATRIX_RELT(m, 4, 2) = m42;  *MATRIX_RELT(m, 4, 3) = m43;  *MATRIX_RELT(m, 4, 4) = m44;

  return(NO_ERROR);

} /* end stuff_four_by_four() */

/*-----------------------------------------------------------
  FixMNITal() - function to compute the "real" talairach coordinates
  from the MNI talaiarch coordinates. Can be done in-place. This is
  based on Matthew Brett's 10/8/98 transform. See
  http://www.mrc-cbu.cam.ac.uk/Imaging/mnispace.html
  -----------------------------------------------------------*/
int fixMNITal(double  xmni, double  ymni, double  zmni,
        double *xtal, double *ytal, double *ztal)
{
  MATRIX *T, *xyzMNI, *xyzTal;


  T = MatrixAlloc(4, 4, MATRIX_REAL);
  if(zmni >= 0.0){
    stuff_four_by_four(T, 
           .9900,  .0000, .0000, 0,
           .0000,  .9688, .0460, 0,
           .0000, -.0485, .9189, 0,
           .0000,  .0000, .0000, 1);
  }
  else {
    stuff_four_by_four(T, 
           .9900,  .0000, .0000, 0,
           .0000,  .9688, .0420, 0,
           .0000, -.0485, .8390, 0,
           .0000,  .0000, .0000, 1);
  }

  xyzMNI = MatrixAlloc(4, 1, MATRIX_REAL);
  xyzMNI->rptr[1][1] = xmni;
  xyzMNI->rptr[2][1] = ymni;
  xyzMNI->rptr[3][1] = zmni;
  xyzMNI->rptr[4][1] = 1.0;

  xyzTal = MatrixMultiply(T,xyzMNI,NULL);

  *xtal = xyzTal->rptr[1][1];
  *ytal = xyzTal->rptr[2][1];
  *ztal = xyzTal->rptr[3][1];

  MatrixFree(&T);
  MatrixFree(&xyzMNI);
  MatrixFree(&xyzTal);

  return(0);
}

/*---------------------------------------------------------------------------*/
/* icosahedron */
/*---------------------------------------------------------------------------*/
static ICOSOHEDRON *read_icosahedron(char *fname) ;

/*---------------------------------------------------------------------------*/
IC_VERTEX ic4_vertices[2562] =
{
  {   .0000,   .0000,  1.0000 },
  {   .8944,   .0000,   .4472 },
  {   .2764,   .8507,   .4472 },
  {  -.7236,   .5257,   .4472 },
  {  -.7236,  -.5257,   .4472 },
  {   .2764,  -.8507,   .4472 },
  {   .7236,  -.5257,  -.4472 },
  {   .7236,   .5257,  -.4472 },
  {  -.2764,   .8507,  -.4472 },
  {  -.8944,   .0000,  -.4472 },
  {  -.2764,  -.8507,  -.4472 },
  {   .0000,   .0000, -1.0000 },
  {  -.4253,  -.3090,   .8507 },
  {  -.8507,   .0000,   .5257 },
  {  -.4253,   .3090,   .8507 },
  {   .1625,  -.5000,   .8507 },
  {  -.2629,  -.8090,   .5257 },
  {   .5257,   .0000,   .8507 },
  {   .6882,  -.5000,   .5257 },
  {   .1625,   .5000,   .8507 },
  {   .6882,   .5000,   .5257 },
  {  -.2629,   .8090,   .5257 },
  {  -.5878,   .8090,   .0000 },
  {   .0000,  1.0000,   .0000 },
  {  -.9511,   .3090,   .0000 },
  {  -.6882,   .5000,  -.5257 },
  {  -.9511,  -.3090,   .0000 },
  {  -.5878,  -.8090,   .0000 },
  {  -.6882,  -.5000,  -.5257 },
  {   .0000,  -1.0000,  .0000 },
  {   .5878,  -.8090,   .0000 },
  {   .2629,  -.8090,  -.5257 },
  {   .9511,  -.3090,   .0000 },
  {   .9511,   .3090,   .0000 },
  {   .8507,   .0000,  -.5257 },
  {   .5878,   .8090,   .0000 },
  {   .2629,   .8090,  -.5257 },
  {  -.5257,   .0000,  -.8507 },
  {  -.1625,   .5000,  -.8507 },
  {  -.1625,  -.5000,  -.8507 },
  {   .4253,  -.3090,  -.8507 },
  {   .4253,   .3090,  -.8507 },
  {  -.2211,  -.1606,   .9619 },
  {  -.4472,   .0000,   .8944 },
  {  -.2211,   .1606,   .9619 },
  {  -.5972,  -.4339,   .6746 },
  {  -.8183,  -.2733,   .5057 },
  {  -.6708,  -.1625,   .7236 },
  {  -.6708,   .1625,   .7236 },
  {  -.8183,   .2733,   .5057 },
  {  -.5972,   .4339,   .6746 },
  {   .0844,  -.2599,   .9619 },
  {  -.1382,  -.4253,   .8944 },
  {   .2281,  -.7020,   .6746 },
  {   .0070,  -.8627,   .5057 },
  {  -.0528,  -.6882,   .7236 },
  {  -.3618,  -.5878,   .7236 },
  {  -.5128,  -.6938,   .5057 },
  {   .2733,   .0000,   .9619 },
  {   .3618,  -.2629,   .8944 },
  {   .7382,   .0000,   .6746 },
  {   .8226,  -.2599,   .5057 },
  {   .6382,  -.2629,   .7236 },
  {  .4472,  -.5257,   .7236 },
  {  .5014,  -.7020,   .5057 },
  {  .0844,   .2599,   .9619 },
  {  .3618,   .2629,   .8944 },
  {  .2281,   .7020,   .6746 },
  {  .5014,   .7020,   .5057 },
  {  .4472,   .5257,   .7236 },
  {  .6382,   .2629,   .7236 },
  {  .8226,   .2599,   .5057 },
  { -.1382,   .4253,   .8944 },
  { -.5128,   .6938,   .5057 },
  { -.3618,   .5878,   .7236 },
  { -.0528,   .6882,   .7236 },
  {  .0070,   .8627,   .5057 },
  { -.6816,   .6938,   .2325 },
  { -.4472,   .8507,   .2764 },
  { -.4492,   .8627,  -.2325 },
  { -.1437,   .9619,  -.2325 },
  { -.3090,   .9511,   .0000 },
  { -.1382,   .9511,   .2764 },
  {  .1437,   .9619,   .2325 },
  { -.8705,   .4339,   .2325 },
  { -.8090,   .5878,   .0000 },
  { -.9593,   .1606,  -.2325 },
  { -.8226,   .2599,  -.5057 },
  { -.8618,   .4253,  -.2764 },
  { -.6708,   .6882,  -.2764 },
  { -.5014,   .7020,  -.5057 },
  { -.9472,   .1625,   .2764 },
  { -.8705,  -.4339,   .2325 },
  { -.9472,  -.1625,   .2764 },
  { -1.0000,  .0000,   .0000 },
  { -.9593,  -.1606,  -.2325 },
  { -.6816,  -.6938,   .2325 },
  { -.8090,  -.5878,   .0000 },
  { -.4492,  -.8627,  -.2325 },
  { -.5014,  -.7020,  -.5057 },
  { -.6708,  -.6882,  -.2764 },
  { -.8618,  -.4253,  -.2764 },
  { -.8226,  -.2599,  -.5057 },
  { -.4472,  -.8507,   .2764 },
  {  .1437,  -.9619,   .2325 },
  { -.1382,  -.9511,   .2764 },
  { -.3090,  -.9511,   .0000 },
  { -.1437,  -.9619,  -.2325 },
  {  .4492,  -.8627,   .2325 },
  {  .3090,  -.9511,   .0000 },
  {  .6816,  -.6938,  -.2325 },
  {  .5128,  -.6938,  -.5057 },
  {  .4472,  -.8507,  -.2764 },
  {  .1382,  -.9511,  -.2764 },
  { -.0070,  -.8627,  -.5057 },
  {  .6708,  -.6882,   .2764 },
  {  .9593,  -.1606,   .2325 },
  {  .8618,  -.4253,   .2764 },
  {  .8090,  -.5878,   .0000 },
  {  .8705,  -.4339,  -.2325 },
  {  .9593,   .1606,   .2325 },
  { 1.0000,   .0000,   .0000 },
  {  .8705,   .4339,  -.2325 },
  {  .8183,   .2733,  -.5057 },
  {  .9472,   .1625,  -.2764 },
  {  .9472,  -.1625,  -.2764 },
  {  .8183,  -.2733,  -.5057 },
  {  .8618,   .4253,   .2764 },
  {  .4492,   .8627,   .2325 },
  {  .6708,   .6882,   .2764 },
  {  .8090,   .5878,   .0000 },
  {  .6816,   .6938,  -.2325 },
  {  .3090,   .9511,   .0000 },
  { -.0070,   .8627,  -.5057 },
  {  .1382,   .9511,  -.2764 },
  {  .4472,   .8507,  -.2764 },
  {  .5128,   .6938,  -.5057 },
  { -.4472,   .5257,  -.7236 },
  { -.2281,   .7020,  -.6746 },
  { -.7382,   .0000,  -.6746 },
  { -.6382,   .2629,  -.7236 },
  { -.3618,   .2629,  -.8944 },
  { -.2733,   .0000,  -.9619 },
  { -.0844,   .2599,  -.9619 },
  { -.6382,  -.2629,  -.7236 },
  { -.2281,  -.7020,  -.6746 },
  { -.4472,  -.5257,  -.7236 },
  { -.3618,  -.2629,  -.8944 },
  { -.0844,  -.2599,  -.9619 },
  {  .0528,  -.6882,  -.7236 },
  {  .5972,  -.4339,  -.6746 },
  {  .3618,  -.5878,  -.7236 },
  {  .1382,  -.4253,  -.8944 },
  {  .2211,  -.1606,  -.9619 },
  {  .6708,  -.1625,  -.7236 },
  {  .5972,   .4339,  -.6746 },
  {  .6708,   .1625,  -.7236 },
  {  .4472,   .0000,  -.8944 },
  {  .2211,   .1606,  -.9619 },
  {  .3618,   .5878,  -.7236 },
  {  .0528,   .6882,  -.7236 },
  {  .1382,   .4253,  -.8944 },
  { -.1116,  -.0811,   .9904 },
  { -.2240,   .0000,   .9746 },
  { -.1116,   .0811,   .9904 },
  { -.3263,  -.2371,   .9150 },
  { -.4417,  -.1564,   .8834 },
  { -.3376,  -.0811,   .9378 },
  { -.3376,   .0811,   .9378 },
  { -.4417,   .1564,   .8834 },
  { -.3263,   .2371,   .9150 },
  { -.5162,  -.3750,   .7700 },
  { -.6406,  -.3013,   .7063 },
  { -.5549,  -.2387,   .7969 },
  { -.6668,  -.4844,   .5663 },
  { -.7784,  -.4034,   .4811 },
  { -.7170,  -.3582,   .5979 },
  { -.7522,  -.2201,   .6210 },
  { -.8425,  -.1380,   .5207 },
  { -.7702,  -.0822,   .6325 },
  { -.5665,   .0823,   .8199 },
  { -.5549,   .2387,   .7969 },
  { -.5665,  -.0823,   .8199 },
  { -.6799,   .0000,   .7333 },
  { -.7702,   .0822,   .6325 },
  { -.6406,   .3013,   .7063 },
  { -.5162,   .3750,   .7700 },
  { -.8425,   .1380,   .5207 },
  { -.7522,   .2201,   .6210 },
  { -.7170,   .3582,   .5979 },
  { -.7784,   .4034,   .4811 },
  { -.6668,   .4844,   .5663 },
  {  .0426,  -.1312,   .9904 },
  { -.0692,  -.2130,   .9746 },
  {  .1246,  -.3836,   .9150 },
  {  .0123,  -.4684,   .8834 },
  { -.0272,  -.3462,   .9378 },
  { -.1815,  -.2960,   .9378 },
  { -.2853,  -.3717,   .8834 },
  {  .1972,  -.6068,   .7700 },
  {  .0886,  -.7023,   .7063 },
  {  .0555,  -.6015,   .7969 },
  {  .2547,  -.7838,   .5663 },
  {  .1431,  -.8649,   .4811 },
  {  .1191,  -.7926,   .5979 },
  { -.0231,  -.7835,   .6210 },
  { -.1292,  -.8439,   .5207 },
  { -.1598,  -.7579,   .6325 },
  { -.2534,  -.5134,   .8199 },
  { -.3985,  -.4540,   .7969 },
  { -.0968,  -.5643,   .8199 },
  { -.2101,  -.6466,   .7333 },
  { -.3162,  -.7071,   .6325 },
  { -.4845,  -.5161,   .7063 },
  { -.3916,  -.7587,   .5207 },
  { -.4418,  -.6474,   .6210 },
  { -.5623,  -.5713,   .5979 },
  { -.6241,  -.6156,   .4811 },
  {  .1380,   .0000,   .9904 },
  {  .1812,  -.1317,   .9746 },
  {  .4034,   .0000,   .9150 },
  {  .4493,  -.1331,   .8834 },
  {  .3208,  -.1328,   .9378 },
  {  .2254,  -.2641,   .9378 },
  {  .2654,  -.3862,   .8834 },
  {  .6381,   .0000,   .7700 },
  {  .6953,  -.1328,   .7063 },
  {  .5892,  -.1331,   .7969 },
  {  .8242,   .0000,   .5663 },
  {  .8668,  -.1312,   .4811 },
  {  .7907,  -.1317,   .5979 },
  {  .7380,  -.2641,   .6210 },
  {  .7627,  -.3836,   .5207 },
  {  .6715,  -.3862,   .6325 },
  {  .4100,  -.3996,   .8199 },
  {  .3086,  -.5193,   .7969 },
  {  .5067,  -.2664,   .8199 },
  {  .5500,  -.3996,   .7333 },
  {  .5748,  -.5193,   .6325 },
  {  .3412,  -.6202,   .7063 },
  {  .6005,  -.6068,   .5207 },
  {  .4792,  -.6202,   .6210 },
  {  .3695,  -.7113,   .5979 },
  {  .3926,  -.7838,   .4811 },
  {  .0426,   .1312,   .9904 },
  {  .1812,   .1317,   .9746 },
  {  .1246,   .3836,   .9150 },
  {  .2654,   .3862,   .8834 },
  {  .2254,   .2641,   .9378 },
  {  .3208,   .1328,   .9378 },
  {  .4493,   .1331,   .8834 },
  {  .1972,   .6068,   .7700 },
  {  .3412,   .6202,   .7063 },
  {  .3086,   .5193,   .7969 },
  {  .2547,   .7838,   .5663 },
  {  .3926,   .7838,   .4811 },
  {  .3695,   .7113,   .5979 },
  {  .4792,   .6202,   .6210 },
  {  .6005,   .6068,   .5207 },
  {  .5748,   .5193,   .6325 },
  {  .5067,   .2664,   .8199 },
  {  .5892,   .1331,   .7969 },
  {  .4100,   .3996,   .8199 },
  {  .5500,   .3996,   .7333 },
  {  .6715,   .3862,   .6325 },
  {  .6953,   .1328,   .7063 },
  {  .7627,   .3836,   .5207 },
  {  .7380,   .2641,   .6210 },
  {  .7907,   .1317,   .5979 },
  {  .8668,   .1312,   .4811 },
  { -.0692,   .2130,   .9746 },
  { -.2853,   .3717,   .8834 },
  { -.1815,   .2960,   .9378 },
  { -.0272,   .3462,   .9378 },
  {  .0123,   .4684,   .8834 },
  { -.4845,   .5161,   .7063 },
  { -.3985,   .4540,   .7969 },
  { -.6241,   .6156,   .4811 },
  { -.5623,   .5713,   .5979 },
  { -.4418,   .6474,   .6210 },
  { -.3916,   .7587,   .5207 },
  { -.3162,   .7071,   .6325 },
  { -.0968,   .5643,   .8199 },
  {  .0555,   .6015,   .7969 },
  { -.2534,   .5134,   .8199 },
  { -.2101,   .6466,   .7333 },
  { -.1598,   .7579,   .6325 },
  {  .0886,   .7023,   .7063 },
  { -.1292,   .8439,   .5207 },
  { -.0231,   .7835,   .6210 },
  {  .1191,   .7926,   .5979 },
  {  .1431,   .8649,   .4811 },
  { -.7094,   .6156,   .3431 },
  { -.6051,   .7029,   .3739 },
  { -.6408,   .7587,   .1173 },
  { -.5240,   .8402,   .1399 },
  { -.5703,   .7802,   .2571 },
  { -.4849,   .7802,   .3951 },
  { -.3595,   .8402,   .4061 },
  { -.5235,   .8439,  -.1173 },
  { -.3830,   .9162,  -.1174 },
  { -.4540,   .8910,   .0000 },
  { -.3663,   .8649,  -.3431 },
  { -.2121,   .9150,  -.3431 },
  { -.3003,   .9243,  -.2355 },
  { -.2287,   .9664,  -.1174 },
  { -.0725,   .9904,  -.1173 },
  { -.1564,   .9877,   .0000 },
  { -.2966,   .9130,   .2801 },
  { -.2030,   .8910,   .4061 },
  { -.3832,   .9130,   .1401 },
  { -.2266,   .9639,   .1401 },
  { -.0700,   .9877,   .1399 },
  { -.0663,   .9162,   .3951 },
  {  .0725,   .9904,   .1173 },
  {  .0028,   .9664,   .2571 },
  {  .0763,   .9243,   .3739 },
  {  .2121,   .9150,   .3431 },
  { -.8047,   .4844,   .3431 },
  { -.7863,   .5713,   .2355 },
  { -.9196,   .3750,   .1173 },
  { -.8910,   .4540,   .0000 },
  { -.8484,   .5161,   .1174 },
  { -.7530,   .6474,   .1174 },
  { -.7071,   .7071,   .0000 },
  { -.9644,   .2371,  -.1173 },
  { -.9199,   .2960,  -.2571 },
  { -.9177,   .3717,  -.1399 },
  { -.9358,   .0811,  -.3431 },
  { -.8668,   .1312,  -.4811 },
  { -.9027,   .2130,  -.3739 },
  { -.8509,   .3462,  -.3951 },
  { -.7627,   .3836,  -.5207 },
  { -.7847,   .4684,  -.4061 },
  { -.7499,   .6466,  -.1401 },
  { -.6371,   .7579,  -.1399 },
  { -.8467,   .5134,  -.1401 },
  { -.7766,   .5643,  -.2801 },
  { -.6880,   .6015,  -.4061 },
  { -.5658,   .7835,  -.2571 },
  { -.6005,   .6068,  -.5207 },
  { -.5922,   .7023,  -.3951 },
  { -.4815,   .7926,  -.3739 },
  { -.3926,   .7838,  -.4811 },
  { -.8555,   .3582,   .3739 },
  { -.9101,   .0822,   .4061 },
  { -.8919,   .2201,   .3951 },
  { -.9182,   .3013,   .2571 },
  { -.9610,   .2387,   .1399 },
  { -.8919,  -.2201,   .3951 },
  { -.9101,  -.0822,   .4061 },
  { -.8047,  -.4844,   .3431 },
  { -.8555,  -.3582,   .3739 },
  { -.9182,  -.3013,   .2571 },
  { -.9196,  -.3750,   .1173 },
  { -.9610,  -.2387,   .1399 },
  { -.9867,   .0823,   .1401 },
  { -.9877,   .1564,   .0000 },
  { -.9600,   .0000,   .2801 },
  { -.9867,  -.0823,   .1401 },
  { -.9877,  -.1564,   .0000 },
  { -.9898,   .0811,  -.1174 },
  { -.9644,  -.2371,  -.1173 },
  { -.9898,  -.0811,  -.1174 },
  { -.9719,   .0000,  -.2355 },
  { -.9358,  -.0811,  -.3431 },
  { -.7094,  -.6156,   .3431 },
  { -.7863,  -.5713,   .2355 },
  { -.6408,  -.7587,   .1173 },
  { -.7071,  -.7071,   .0000 },
  { -.7530,  -.6474,   .1174 },
  { -.8484,  -.5161,   .1174 },
  { -.8910,  -.4540,   .0000 },
  { -.5235,  -.8439,  -.1173 },
  { -.5658,  -.7835,  -.2571 },
  { -.6371,  -.7579,  -.1399 },
  { -.3663,  -.8649,  -.3431 },
  { -.3926,  -.7838,  -.4811 },
  { -.4815,  -.7926,  -.3739 },
  { -.5922,  -.7023,  -.3951 },
  { -.6005,  -.6068,  -.5207 },
  { -.6880,  -.6015,  -.4061 },
  { -.8467,  -.5134,  -.1401 },
  { -.9177,  -.3717,  -.1399 },
  { -.7499,  -.6466,  -.1401 },
  { -.7766,  -.5643,  -.2801 },
  { -.7847,  -.4684,  -.4061 },
  { -.9199,  -.2960,  -.2571 },
  { -.7627,  -.3836,  -.5207 },
  { -.8509,  -.3462,  -.3951 },
  { -.9027,  -.2130,  -.3739 },
  { -.8668,  -.1312,  -.4811 },
  { -.6051,  -.7029,   .3739 },
  { -.3595,  -.8402,   .4061 },
  { -.4849,  -.7802,   .3951 },
  { -.5703,  -.7802,   .2571 },
  { -.5239,  -.8402,   .1399 },
  { -.0663,  -.9162,   .3951 },
  { -.2030,  -.8910,   .4061 },
  {  .2121,  -.9150,   .3431 },
  {  .0763,  -.9243,   .3739 },
  {  .0028,  -.9664,   .2571 },
  {  .0725,  -.9904,   .1173 },
  { -.0700,  -.9877,   .1399 },
  { -.3832,  -.9130,   .1401 },
  { -.4540,  -.8910,   .0000 },
  { -.2966,  -.9130,   .2801 },
  { -.2266,  -.9639,   .1401 },
  { -.1564,  -.9877,   .0000 },
  { -.3830,  -.9162,  -.1174 },
  { -.0725,  -.9904,  -.1173 },
  { -.2287,  -.9664,  -.1174 },
  { -.3003,  -.9243,  -.2355 },
  { -.2121,  -.9150,  -.3431 },
  {  .3663,  -.8649,   .3431 },
  {  .3003,  -.9243,   .2355 },
  {  .5235,  -.8439,   .1173 },
  {  .4540,  -.8910,   .0000 },
  {  .3830,  -.9162,   .1174 },
  {  .2287,  -.9664,   .1174 },
  {  .1564,  -.9877,   .0000 },
  {  .6408,  -.7587,  -.1173 },
  {  .5703,  -.7802,  -.2571 },
  {  .5240,  -.8402,  -.1399 },
  {  .7094,  -.6156,  -.3431 },
  {  .6241,  -.6156,  -.4811 },
  {  .6051,  -.7029,  -.3739 },
  {  .4849,  -.7802,  -.3951 },
  {  .3916,  -.7587,  -.5207 },
  {  .3595,  -.8402,  -.4061 },
  {  .2266,  -.9639,  -.1401 },
  {  .0700,  -.9877,  -.1399 },
  {  .3832,  -.9130,  -.1401 },
  {  .2966,  -.9130,  -.2801 },
  {  .2030,  -.8910,  -.4061 },
  { -.0028,  -.9664,  -.2571 },
  {  .1292,  -.8439,  -.5207 },
  {  .0663,  -.9162,  -.3951 },
  { -.0763,  -.9243,  -.3739 },
  { -.1431,  -.8649,  -.4811 },
  {  .4815,  -.7926,   .3739 },
  {  .6880,  -.6015,   .4061 },
  {  .5922,  -.7023,   .3951 },
  {  .5658,  -.7835,   .2571 },
  {  .6371,  -.7579,   .1399 },
  {  .8509,  -.3462,   .3951 },
  {  .7847,  -.4684,   .4061 },
  {  .9358,  -.0811,   .3431 },
  {  .9027,  -.2130,   .3739 },
  {  .9199,  -.2960,   .2571 },
  {  .9644,  -.2371,   .1173 },
  {  .9177,  -.3717,   .1399 },
  {  .7499,  -.6466,   .1401 },
  {  .7071,  -.7071,   .0000 },
  {  .7766,  -.5643,   .2801 },
  {  .8467,  -.5134,   .1401 },
  {  .8910,  -.4540,   .0000 },
  {  .7530,  -.6474,  -.1174 },
  {  .9196,  -.3750,  -.1173 },
  {  .8484,  -.5161,  -.1174 },
  {  .7863,  -.5713,  -.2355 },
  {  .8047,  -.4844,  -.3431 },
  {  .9358,   .0811,   .3431 },
  {  .9719,   .0000,   .2355 },
  {  .9644,   .2371,   .1173 },
  {  .9877,   .1564,   .0000 },
  {  .9898,   .0811,   .1174 },
  {  .9898,  -.0811,   .1174 },
  {  .9877,  -.1564,   .0000 },
  {  .9196,   .3750,  -.1173 },
  {  .9182,   .3013,  -.2571 },
  {  .9610,   .2387,  -.1399 },
  {  .8047,   .4844,  -.3431 },
  {  .7784,   .4034,  -.4811 },
  {  .8555,   .3582,  -.3739 },
  {  .8919,   .2201,  -.3951 },
  {  .8425,   .1380,  -.5207 },
  {  .9101,   .0822,  -.4061 },
  {  .9867,  -.0823,  -.1401 },
  {  .9610,  -.2387,  -.1399 },
  {  .9867,   .0823,  -.1401 },
  {  .9600,   .0000,  -.2801 },
  {  .9101,  -.0822,  -.4061 },
  {  .9182,  -.3013,  -.2571 },
  {  .8425,  -.1380,  -.5207 },
  {  .8919,  -.2201,  -.3951 },
  {  .8555,  -.3582,  -.3739 },
  {  .7784,  -.4034,  -.4811 },
  {  .9027,   .2130,   .3739 },
  {  .7847,   .4684,   .4061 },
  {  .8509,   .3462,   .3951 },
  {  .9199,   .2960,   .2571 },
  {  .9177,   .3717,   .1399 },
  {  .5922,   .7023,   .3951 },
  {  .6880,   .6015,   .4061 },
  {  .3663,   .8649,   .3431 },
  {  .4815,   .7926,   .3739 },
  {  .5658,   .7835,   .2571 },
  {  .5235,   .8439,   .1173 },
  {  .6371,   .7579,   .1399 },
  {  .8467,   .5134,   .1401 },
  {  .8910,   .4540,   .0000 },
  {  .7766,   .5643,   .2801 },
  {  .7499,   .6466,   .1401 },
  {  .7071,   .7071,   .0000 },
  {  .8484,   .5161,  -.1174 },
  {  .6408,   .7587,  -.1173 },
  {  .7530,   .6474,  -.1174 },
  {  .7863,   .5713,  -.2355 },
  {  .7094,   .6156,  -.3431 },
  {  .3003,   .9243,   .2355 },
  {  .1564,   .9877,   .0000 },
  {  .2287,   .9664,   .1174 },
  {  .3830,   .9162,   .1174 },
  {  .4540,   .8910,   .0000 },
  { -.0028,   .9664,  -.2571 },
  {  .0700,   .9877,  -.1399 },
  { -.1431,   .8649,  -.4811 },
  { -.0763,   .9243,  -.3739 },
  {  .0663,   .9162,  -.3951 },
  {  .1292,   .8439,  -.5207 },
  {  .2030,   .8910,  -.4061 },
  {  .3832,   .9130,  -.1401 },
  {  .5240,   .8402,  -.1399 },
  {  .2266,   .9639,  -.1401 },
  {  .2966,   .9130,  -.2801 },
  {  .3595,   .8402,  -.4061 },
  {  .5703,   .7802,  -.2571 },
  {  .3916,   .7587,  -.5207 },
  {  .4849,   .7802,  -.3951 },
  {  .6051,   .7029,  -.3739 },
  {  .6241,   .6156,  -.4811 },
  { -.3695,   .7113,  -.5979 },
  { -.2547,   .7838,  -.5663 },
  { -.5748,   .5193,  -.6325 },
  { -.4792,   .6202,  -.6210 },
  { -.3412,   .6202,  -.7063 },
  { -.3086,   .5193,  -.7969 },
  { -.1972,   .6068,  -.7700 },
  { -.7380,   .2641,  -.6210 },
  { -.6715,   .3862,  -.6325 },
  { -.8242,   .0000,  -.5663 },
  { -.7907,   .1317,  -.5979 },
  { -.6953,   .1328,  -.7063 },
  { -.6381,   .0000,  -.7700 },
  { -.5892,   .1331,  -.7969 },
  { -.4100,   .3996,  -.8199 },
  { -.2654,   .3862,  -.8834 },
  { -.5500,   .3996,  -.7333 },
  { -.5067,   .2664,  -.8199 },
  { -.4493,   .1331,  -.8834 },
  { -.2254,   .2641,  -.9378 },
  { -.1246,   .3836,  -.9150 },
  { -.4034,   .0000,  -.9150 },
  { -.3208,   .1328,  -.9378 },
  { -.1812,   .1317,  -.9746 },
  { -.1380,   .0000,  -.9904 },
  { -.0426,   .1312,  -.9904 },
  { -.7907,  -.1317,  -.5979 },
  { -.6715,  -.3862,  -.6325 },
  { -.7380,  -.2641,  -.6210 },
  { -.6953,  -.1328,  -.7063 },
  { -.5892,  -.1331,  -.7969 },
  { -.4792,  -.6202,  -.6210 },
  { -.5748,  -.5193,  -.6325 },
  { -.2547,  -.7838,  -.5663 },
  { -.3695,  -.7113,  -.5979 },
  { -.3412,  -.6202,  -.7063 },
  { -.1972,  -.6068,  -.7700 },
  { -.3086,  -.5193,  -.7969 },
  { -.5067,  -.2664,  -.8199 },
  { -.4493,  -.1331,  -.8834 },
  { -.5500,  -.3996,  -.7333 },
  { -.4100,  -.3996,  -.8199 },
  { -.2654,  -.3862,  -.8834 },
  { -.3208,  -.1328,  -.9378 },
  { -.1246,  -.3836,  -.9150 },
  { -.2254,  -.2641,  -.9378 },
  { -.1812,  -.1317,  -.9746 },
  { -.0426,  -.1312,  -.9904 },
  { -.1191,  -.7926,  -.5979 },
  {  .1598,  -.7579,  -.6325 },
  {  .0231,  -.7835,  -.6210 },
  { -.0886,  -.7023,  -.7063 },
  { -.0555,  -.6015,  -.7969 },
  {  .4418,  -.6474,  -.6210 },
  {  .3162,  -.7071,  -.6325 },
  {  .6668,  -.4844,  -.5663 },
  {  .5623,  -.5713,  -.5979 },
  {  .4845,  -.5161,  -.7063 },
  {  .5162,  -.3750,  -.7700 },
  {  .3985,  -.4540,  -.7969 },
  {  .0968,  -.5643,  -.8199 },
  { -.0123,  -.4684,  -.8834 },
  {  .2101,  -.6466,  -.7333 },
  {  .2534,  -.5134,  -.8199 },
  {  .2853,  -.3717,  -.8834 },
  {  .0272,  -.3462,  -.9378 },
  {  .3263,  -.2371,  -.9150 },
  {  .1815,  -.2960,  -.9378 },
  {  .0692,  -.2130,  -.9746 },
  {  .1116,  -.0811,  -.9904 },
  {  .7170,  -.3582,  -.5979 },
  {  .7702,  -.0822,  -.6325 },
  {  .7522,  -.2201,  -.6210 },
  {  .6406,  -.3013,  -.7063 },
  {  .5549,  -.2387,  -.7969 },
  {  .7522,   .2201,  -.6210 },
  {  .7702,   .0822,  -.6325 },
  {  .6668,   .4844,  -.5663 },
  {  .7170,   .3582,  -.5979 },
  {  .6406,   .3013,  -.7063 },
  {  .5162,   .3750,  -.7700 },
  {  .5549,   .2387,  -.7969 },
  {  .5665,  -.0823,  -.8199 },
  {  .4417,  -.1564,  -.8834 },
  {  .6799,   .0000,  -.7333 },
  {  .5665,   .0823,  -.8199 },
  {  .4417,   .1564,  -.8834 },
  {  .3376,  -.0811,  -.9378 },
  {  .3263,   .2371,  -.9150 },
  {  .3376,   .0811,  -.9378 },
  {  .2240,   .0000,  -.9746 },
  {  .1116,   .0811,  -.9904 },
  {  .5623,   .5713,  -.5979 },
  {  .3162,   .7071,  -.6325 },
  {  .4418,   .6474,  -.6210 },
  {  .4845,   .5161,  -.7063 },
  {  .3985,   .4540,  -.7969 },
  {  .0231,   .7835,  -.6210 },
  {  .1598,   .7579,  -.6325 },
  { -.1191,   .7926,  -.5979 },
  { -.0886,   .7023,  -.7063 },
  { -.0555,   .6015,  -.7969 },
  {  .2534,   .5134,  -.8199 },
  {  .2853,   .3717,  -.8834 },
  {  .2101,   .6466,  -.7333 },
  {  .0968,   .5643,  -.8199 },
  { -.0123,   .4684,  -.8834 },
  {  .1815,   .2960,  -.9378 },
  {  .0272,   .3462,  -.9378 },
  {  .0692,   .2130,  -.9746 },
  { -.0559,  -.0406,   .9976 },
  { -.1120,   .0000,   .9937 },
  { -.0559,   .0406,   .9976 },
  { -.1667,  -.1211,   .9785 },
  { -.2233,  -.0806,   .9714 },
  { -.1682,  -.0406,   .9849 },
  { -.1682,   .0406,   .9849 },
  { -.2233,   .0806,   .9714 },
  { -.1667,   .1211,   .9785 },
  { -.2744,  -.1993,   .9407 },
  { -.3330,  -.1596,   .9293 },
  { -.2801,  -.1212,   .9523 },
  { -.3767,  -.2737,   .8850 },
  { -.4349,  -.2334,   .8697 },
  { -.3850,  -.1973,   .9016 },
  { -.3906,  -.1191,   .9128 },
  { -.4458,  -.0785,   .8917 },
  { -.3934,  -.0407,   .9185 },
  { -.2815,   .0407,   .9587 },
  { -.2801,   .1212,   .9523 },
  { -.2815,  -.0407,   .9587 },
  { -.3387,   .0000,   .9409 },
  { -.3934,   .0407,   .9185 },
  { -.3330,   .1596,   .9293 },
  { -.2744,   .1993,   .9407 },
  { -.4458,   .0785,   .8917 },
  { -.3906,   .1191,   .9128 },
  { -.3850,   .1973,   .9016 },
  { -.4349,   .2334,   .8697 },
  { -.3767,   .2737,   .8850 },
  { -.4719,  -.3428,   .8123 },
  { -.5369,  -.3077,   .7855 },
  { -.4916,  -.2747,   .8263 },
  { -.5580,  -.4054,   .7240 },
  { -.6205,  -.3685,   .6922 },
  { -.5802,  -.3392,   .7405 },
  { -.5992,  -.2706,   .7535 },
  { -.6574,  -.2325,   .7168 },
  { -.6148,  -.2012,   .7626 },
  { -.6335,  -.4603,   .6220 },
  { -.6936,  -.4224,   .5835 },
  { -.6593,  -.3974,   .6383 },
  { -.6969,  -.5063,   .5080 },
  { -.7528,  -.4657,   .4653 },
  { -.7250,  -.4454,   .5254 },
  { -.7495,  -.3817,   .5408 },
  { -.8002,  -.3391,   .4946 },
  { -.7702,  -.3168,   .5536 },
  { -.6987,  -.2616,   .6659 },
  { -.7134,  -.1918,   .6740 },
  { -.6806,  -.3306,   .6538 },
  { -.7366,  -.2899,   .6111 },
  { -.7873,  -.2473,   .5648 },
  { -.7631,  -.1515,   .6283 },
  { -.7227,  -.1227,   .6801 },
  { -.8324,  -.2061,   .5144 },
  { -.7999,  -.1796,   .5727 },
  { -.8085,  -.1104,   .5781 },
  { -.8486,  -.0691,   .5245 },
  { -.8129,  -.0412,   .5809 },
  { -.5000,   .1982,   .8430 },
  { -.4916,   .2747,   .8263 },
  { -.5086,   .0413,   .8600 },
  { -.5057,   .1198,   .8543 },
  { -.5625,   .1610,   .8110 },
  { -.6207,   .1228,   .7743 },
  { -.6148,   .2012,   .7626 },
  { -.5057,  -.1198,   .8543 },
  { -.5086,  -.0413,   .8600 },
  { -.5000,  -.1982,   .8430 },
  { -.5625,  -.1610,   .8110 },
  { -.6207,  -.1228,   .7743 },
  { -.6253,   .0413,   .7793 },
  { -.6776,   .0815,   .7309 },
  { -.5685,   .0000,   .8227 },
  { -.6253,  -.0413,   .7793 },
  { -.6776,  -.0815,   .7309 },
  { -.7273,   .0413,   .6851 },
  { -.7227,   .1227,   .6801 },
  { -.7273,  -.0413,   .6851 },
  { -.7728,   .0000,   .6346 },
  { -.8129,   .0412,   .5809 },
  { -.5369,   .3077,   .7855 },
  { -.4719,   .3428,   .8123 },
  { -.6574,   .2325,   .7168 },
  { -.5992,   .2706,   .7535 },
  { -.5802,   .3392,   .7405 },
  { -.6205,   .3685,   .6922 },
  { -.5580,   .4054,   .7240 },
  { -.7631,   .1515,   .6283 },
  { -.7134,   .1918,   .6740 },
  { -.8486,   .0691,   .5245 },
  { -.8085,   .1104,   .5781 },
  { -.7999,   .1796,   .5727 },
  { -.8324,   .2061,   .5144 },
  { -.7873,   .2473,   .5648 },
  { -.6806,   .3306,   .6538 },
  { -.6593,   .3974,   .6383 },
  { -.6987,   .2616,   .6659 },
  { -.7366,   .2899,   .6111 },
  { -.7702,   .3168,   .5536 },
  { -.6936,   .4224,   .5835 },
  { -.6335,   .4603,   .6220 },
  { -.8002,   .3391,   .4946 },
  { -.7495,   .3817,   .5408 },
  { -.7250,   .4454,   .5254 },
  { -.7528,   .4657,   .4653 },
  { -.6969,   .5063,   .5080 },
  {  .0214,  -.0658,   .9976 },
  { -.0346,  -.1065,   .9937 },
  {  .0637,  -.1960,   .9785 },
  {  .0076,  -.2372,   .9714 },
  { -.0133,  -.1725,   .9849 },
  { -.0906,  -.1474,   .9849 },
  { -.1456,  -.1874,   .9714 },
  {  .1048,  -.3225,   .9407 },
  {  .0489,  -.3660,   .9293 },
  {  .0287,  -.3038,   .9523 },
  {  .1439,  -.4429,   .8850 },
  {  .0876,  -.4857,   .8697 },
  {  .0686,  -.4271,   .9016 },
  { -.0075,  -.4083,   .9128 },
  { -.0632,  -.4483,   .8917 },
  { -.0829,  -.3867,   .9185 },
  { -.1257,  -.2552,   .9587 },
  { -.2018,  -.2289,   .9523 },
  { -.0483,  -.2803,   .9587 },
  { -.1047,  -.3221,   .9409 },
  { -.1603,  -.3616,   .9185 },
  { -.2547,  -.2674,   .9293 },
  { -.2124,  -.3998,   .8917 },
  { -.2340,  -.3347,   .9128 },
  { -.3066,  -.3052,   .9016 },
  { -.3564,  -.3414,   .8697 },
  {  .1802,  -.5547,   .8123 },
  {  .1267,  -.6057,   .7855 },
  {  .1093,  -.5525,   .8263 },
  {  .2131,  -.6560,   .7240 },
  {  .1587,  -.7040,   .6922 },
  {  .1433,  -.6566,   .7405 },
  {  .0722,  -.6535,   .7535 },
  {  .0179,  -.6970,   .7168 },
  {  .0014,  -.6468,   .7626 },
  {  .2420,  -.7447,   .6220 },
  {  .1874,  -.7902,   .5835 },
  {  .1742,  -.7498,   .6383 },
  {  .2662,  -.8192,   .5080 },
  {  .2102,  -.8598,   .4653 },
  {  .1995,  -.8271,   .5254 },
  {  .1314,  -.8308,   .5408 },
  {  .0752,  -.8659,   .4946 },
  {  .0633,  -.8304,   .5536 },
  {  .0328,  -.7453,   .6659 },
  { -.0380,  -.7377,   .6740 },
  {  .1041,  -.7494,   .6538 },
  {  .0481,  -.7901,   .6111 },
  { -.0081,  -.8252,   .5648 },
  { -.0917,  -.7726,   .6283 },
  { -.1066,  -.7253,   .6801 },
  { -.0612,  -.8553,   .5144 },
  { -.0764,  -.8162,   .5727 },
  { -.1448,  -.8030,   .5781 },
  { -.1965,  -.8284,   .5245 },
  { -.2120,  -.7859,   .5809 },
  { -.3430,  -.4143,   .8430 },
  { -.4132,  -.3827,   .8263 },
  { -.1964,  -.4709,   .8600 },
  { -.2702,  -.4440,   .8543 },
  { -.3269,  -.4852,   .8110 },
  { -.3086,  -.5524,   .7743 },
  { -.3813,  -.5225,   .7626 },
  { -.0424,  -.5180,   .8543 },
  { -.1179,  -.4964,   .8600 },
  {  .0340,  -.5368,   .8430 },
  { -.0207,  -.5847,   .8110 },
  { -.0750,  -.6283,   .7743 },
  { -.2325,  -.5820,   .7793 },
  { -.2869,  -.6192,   .7309 },
  { -.1757,  -.5406,   .8227 },
  { -.1540,  -.6075,   .7793 },
  { -.1319,  -.6696,   .7309 },
  { -.2640,  -.6790,   .6851 },
  { -.3401,  -.6494,   .6801 },
  { -.1855,  -.7045,   .6851 },
  { -.2388,  -.7350,   .6346 },
  { -.2904,  -.7604,   .5809 },
  { -.4585,  -.4156,   .7855 },
  { -.4242,  -.5534,   .7168 },
  { -.4425,  -.4862,   .7535 },
  { -.5019,  -.4470,   .7405 },
  { -.5422,  -.4762,   .6922 },
  { -.3799,  -.6789,   .6283 },
  { -.4028,  -.6192,   .6740 },
  { -.3280,  -.7857,   .5245 },
  { -.3548,  -.7348,   .5781 },
  { -.4180,  -.7052,   .5727 },
  { -.4532,  -.7280,   .5144 },
  { -.4785,  -.6723,   .5648 },
  { -.5247,  -.5451,   .6538 },
  { -.5816,  -.5042,   .6383 },
  { -.4647,  -.5837,   .6659 },
  { -.5034,  -.6109,   .6111 },
  { -.5393,  -.6346,   .5536 },
  { -.6160,  -.5291,   .5835 },
  { -.5698,  -.6563,   .4946 },
  { -.5947,  -.5949,   .5408 },
  { -.6476,  -.5519,   .5254 },
  { -.6755,  -.5721,   .4653 },
  {  .0691,   .0000,   .9976 },
  {  .0906,  -.0658,   .9937 },
  {  .2061,   .0000,   .9785 },
  {  .2280,  -.0660,   .9714 },
  {  .1600,  -.0660,   .9849 },
  {  .1122,  -.1317,   .9849 },
  {  .1333,  -.1964,   .9714 },
  {  .3391,   .0000,   .9407 },
  {  .3632,  -.0666,   .9293 },
  {  .2978,  -.0666,   .9523 },
  {  .4657,   .0000,   .8850 },
  {  .4890,  -.0667,   .8697 },
  {  .4274,  -.0667,   .9016 },
  {  .3860,  -.1333,   .9128 },
  {  .4068,  -.1986,   .8917 },
  {  .3422,  -.1983,   .9185 },
  {  .2039,  -.1984,   .9587 },
  {  .1553,  -.2627,   .9523 },
  {  .2517,  -.1326,   .9587 },
  {  .2740,  -.1991,   .9409 },
  {  .2944,  -.2642,   .9185 },
  {  .1756,  -.3249,   .9293 },
  {  .3146,  -.3255,   .8917 },
  {  .2460,  -.3259,   .9128 },
  {  .1955,  -.3859,   .9016 },
  {  .2146,  -.4445,   .8697 },
  {  .5833,   .0000,   .8123 },
  {  .6152,  -.0667,   .7855 },
  {  .5592,  -.0667,   .8263 },
  {  .6898,   .0000,   .7240 },
  {  .7186,  -.0666,   .6922 },
  {  .6688,  -.0666,   .7405 },
  {  .6438,  -.1333,   .7535 },
  {  .6685,  -.1983,   .7168 },
  {  .6156,  -.1986,   .7626 },
  {  .7831,   .0000,   .6220 },
  {  .8094,  -.0660,   .5835 },
  {  .7669,  -.0660,   .6383 },
  {  .8614,   .0000,   .5080 },
  {  .8827,  -.0658,   .4653 },
  {  .8483,  -.0658,   .5254 },
  {  .8308,  -.1317,   .5408 },
  {  .8467,  -.1960,   .4946 },
  {  .8093,  -.1964,   .5536 },
  {  .7190,  -.1991,   .6659 },
  {  .6898,  -.2642,   .6740 },
  {  .7449,  -.1326,   .6538 },
  {  .7663,  -.1984,   .6111 },
  {  .7823,  -.2627,   .5648 },
  {  .7064,  -.3259,   .6283 },
  {  .6569,  -.3255,   .6801 },
  {  .7946,  -.3225,   .5144 },
  {  .7527,  -.3249,   .5727 },
  {  .7190,  -.3859,   .5781 },
  {  .7272,  -.4429,   .5245 },
  {  .6819,  -.4445,   .5809 },
  {  .2880,  -.4543,   .8430 },
  {  .2363,  -.5112,   .8263 },
  {  .3872,  -.3323,   .8600 },
  {  .3387,  -.3941,   .8543 },
  {  .3604,  -.4609,   .8110 },
  {  .4300,  -.4642,   .7743 },
  {  .3791,  -.5241,   .7626 },
  {  .4795,  -.2004,   .8543 },
  {  .4357,  -.2655,   .8600 },
  {  .5210,  -.1335,   .8430 },
  {  .5497,  -.2004,   .8110 },
  {  .5744,  -.2655,   .7743 },
  {  .4816,  -.4010,   .7793 },
  {  .5003,  -.4642,   .7309 },
  {  .4599,  -.3341,   .8227 },
  {  .5302,  -.3341,   .7793 },
  {  .5961,  -.3323,   .7309 },
  {  .5642,  -.4609,   .6851 },
  {  .5126,  -.5241,   .6801 },
  {  .6127,  -.3941,   .6851 },
  {  .6252,  -.4543,   .6346 },
  {  .6334,  -.5112,   .5809 },
  {  .2536,  -.5645,   .7855 },
  {  .3952,  -.5745,   .7168 },
  {  .3257,  -.5711,   .7535 },
  {  .2700,  -.6155,   .7405 },
  {  .2854,  -.6628,   .6922 },
  {  .5283,  -.5711,   .6283 },
  {  .4644,  -.5745,   .6740 },
  {  .6459,  -.5547,   .5245 },
  {  .5892,  -.5645,   .5781 },
  {  .5416,  -.6155,   .5727 },
  {  .5523,  -.6560,   .5144 },
  {  .4915,  -.6628,   .5648 },
  {  .3563,  -.6675,   .6538 },
  {  .2998,  -.7090,   .6383 },
  {  .4115,  -.6223,   .6659 },
  {  .4255,  -.6675,   .6111 },
  {  .4369,  -.7090,   .5536 },
  {  .3129,  -.7494,   .5835 },
  {  .4481,  -.7447,   .4946 },
  {  .3820,  -.7494,   .5408 },
  {  .3247,  -.7864,   .5254 },
  {  .3353,  -.8192,   .4653 },
  {  .0214,   .0658,   .9976 },
  {  .0906,   .0658,   .9937 },
  {  .0637,   .1960,   .9785 },
  {  .1333,   .1964,   .9714 },
  {  .1122,   .1317,   .9849 },
  {  .1600,   .0660,   .9849 },
  {  .2280,   .0660,   .9714 },
  {  .1048,   .3225,   .9407 },
  {  .1756,   .3249,   .9293 },
  {  .1553,   .2627,   .9523 },
  {  .1439,   .4429,   .8850 },
  {  .2146,   .4445,   .8697 },
  {  .1955,   .3859,   .9016 },
  {  .2460,   .3259,   .9128 },
  {  .3146,   .3255,   .8917 },
  {  .2944,   .2642,   .9185 },
  {  .2517,   .1326,   .9587 },
  {  .2978,   .0666,   .9523 },
  {  .2039,   .1984,   .9587 },
  {  .2740,   .1991,   .9409 },
  {  .3422,   .1983,   .9185 },
  {  .3632,   .0666,   .9293 },
  {  .4068,   .1986,   .8917 },
  {  .3860,   .1333,   .9128 },
  {  .4274,   .0667,   .9016 },
  {  .4890,   .0667,   .8697 },
  {  .1802,   .5547,   .8123 },
  {  .2536,   .5645,   .7855 },
  {  .2363,   .5112,   .8263 },
  {  .2131,   .6560,   .7240 },
  {  .2854,   .6628,   .6922 },
  {  .2700,   .6155,   .7405 },
  {  .3257,   .5711,   .7535 },
  {  .3952,   .5745,   .7168 },
  {  .3791,   .5241,   .7626 },
  {  .2420,   .7447,   .6220 },
  {  .3129,   .7494,   .5835 },
  {  .2998,   .7090,   .6383 },
  {  .2662,   .8192,   .5080 },
  {  .3353,   .8192,   .4653 },
  {  .3247,   .7864,   .5254 },
  {  .3820,   .7494,   .5408 },
  {  .4481,   .7447,   .4946 },
  {  .4369,   .7090,   .5536 },
  {  .4115,   .6223,   .6659 },
  {  .4644,   .5745,   .6740 },
  {  .3563,   .6675,   .6538 },
  {  .4255,   .6675,   .6111 },
  {  .4915,   .6628,   .5648 },
  {  .5283,   .5711,   .6283 },
  {  .5126,   .5241,   .6801 },
  {  .5523,   .6560,   .5144 },
  {  .5416,   .6155,   .5727 },
  {  .5892,   .5645,   .5781 },
  {  .6459,   .5547,   .5245 },
  {  .6334,   .5112,   .5809 },
  {  .5210,   .1335,   .8430 },
  {  .5592,   .0667,   .8263 },
  {  .4357,   .2655,   .8600 },
  {  .4795,   .2004,   .8543 },
  {  .5497,   .2004,   .8110 },
  {  .5744,   .2655,   .7743 },
  {  .6156,   .1986,   .7626 },
  {  .3387,   .3941,   .8543 },
  {  .3872,   .3323,   .8600 },
  {  .2880,   .4543,   .8430 },
  {  .3604,   .4609,   .8110 },
  {  .4300,   .4642,   .7743 },
  {  .5302,   .3341,   .7793 },
  {  .5961,   .3323,   .7309 },
  {  .4599,   .3341,   .8227 },
  {  .4816,   .4010,   .7793 },
  {  .5003,   .4642,   .7309 },
  {  .6127,   .3941,   .6851 },
  {  .6569,   .3255,   .6801 },
  {  .5642,   .4609,   .6851 },
  {  .6252,   .4543,   .6346 },
  {  .6819,   .4445,   .5809 },
  {  .6152,   .0667,   .7855 },
  {  .6685,   .1983,   .7168 },
  {  .6438,   .1333,   .7535 },
  {  .6688,   .0666,   .7405 },
  {  .7186,   .0666,   .6922 },
  {  .7064,   .3259,   .6283 },
  {  .6898,   .2642,   .6740 },
  {  .7272,   .4429,   .5245 },
  {  .7190,   .3859,   .5781 },
  {  .7527,   .3249,   .5727 },
  {  .7946,   .3225,   .5144 },
  {  .7823,   .2627,   .5648 },
  {  .7449,   .1326,   .6538 },
  {  .7669,   .0660,   .6383 },
  {  .7190,   .1991,   .6659 },
  {  .7663,   .1984,   .6111 },
  {  .8093,   .1964,   .5536 },
  {  .8094,   .0660,   .5835 },
  {  .8467,   .1960,   .4946 },
  {  .8308,   .1317,   .5408 },
  {  .8483,   .0658,   .5254 },
  {  .8827,   .0658,   .4653 },
  { -.0346,   .1065,   .9937 },
  { -.1456,   .1874,   .9714 },
  { -.0906,   .1474,   .9849 },
  { -.0133,   .1725,   .9849 },
  {  .0076,   .2372,   .9714 },
  { -.2547,   .2674,   .9293 },
  { -.2018,   .2289,   .9523 },
  { -.3564,   .3414,   .8697 },
  { -.3066,   .3052,   .9016 },
  { -.2340,   .3347,   .9128 },
  { -.2124,   .3998,   .8917 },
  { -.1603,   .3616,   .9185 },
  { -.0483,   .2803,   .9587 },
  {  .0287,   .3038,   .9523 },
  { -.1257,   .2552,   .9587 },
  { -.1047,   .3221,   .9409 },
  { -.0829,   .3867,   .9185 },
  {  .0489,   .3660,   .9293 },
  { -.0632,   .4483,   .8917 },
  { -.0075,   .4083,   .9128 },
  {  .0686,   .4271,   .9016 },
  {  .0876,   .4857,   .8697 },
  { -.4585,   .4156,   .7855 },
  { -.4132,   .3827,   .8263 },
  { -.5422,   .4762,   .6922 },
  { -.5019,   .4470,   .7405 },
  { -.4425,   .4862,   .7535 },
  { -.4242,   .5534,   .7168 },
  { -.3813,   .5225,   .7626 },
  { -.6160,   .5291,   .5835 },
  { -.5816,   .5042,   .6383 },
  { -.6755,   .5721,   .4653 },
  { -.6476,   .5519,   .5254 },
  { -.5947,   .5949,   .5408 },
  { -.5698,   .6563,   .4946 },
  { -.5393,   .6346,   .5536 },
  { -.4647,   .5837,   .6659 },
  { -.4028,   .6192,   .6740 },
  { -.5247,   .5451,   .6538 },
  { -.5034,   .6109,   .6111 },
  { -.4785,   .6723,   .5648 },
  { -.3799,   .6789,   .6283 },
  { -.3401,   .6494,   .6801 },
  { -.4532,   .7280,   .5144 },
  { -.4180,   .7052,   .5727 },
  { -.3548,   .7348,   .5781 },
  { -.3280,   .7857,   .5245 },
  { -.2904,   .7604,   .5809 },
  {  .0340,   .5368,   .8430 },
  {  .1093,   .5525,   .8263 },
  { -.1179,   .4964,   .8600 },
  { -.0424,   .5180,   .8543 },
  { -.0207,   .5847,   .8110 },
  { -.0750,   .6283,   .7743 },
  {  .0014,   .6468,   .7626 },
  { -.2702,   .4440,   .8543 },
  { -.1964,   .4709,   .8600 },
  { -.3430,   .4143,   .8430 },
  { -.3269,   .4852,   .8110 },
  { -.3086,   .5524,   .7743 },
  { -.1540,   .6075,   .7793 },
  { -.1319,   .6696,   .7309 },
  { -.1757,   .5406,   .8227 },
  { -.2325,   .5820,   .7793 },
  { -.2869,   .6192,   .7309 },
  { -.1855,   .7045,   .6851 },
  { -.1066,   .7253,   .6801 },
  { -.2640,   .6790,   .6851 },
  { -.2388,   .7350,   .6346 },
  { -.2120,   .7859,   .5809 },
  {  .1267,   .6057,   .7855 },
  {  .0179,   .6970,   .7168 },
  {  .0722,   .6535,   .7535 },
  {  .1433,   .6566,   .7405 },
  {  .1587,   .7040,   .6922 },
  { -.0917,   .7726,   .6283 },
  { -.0380,   .7377,   .6740 },
  { -.1965,   .8284,   .5245 },
  { -.1449,   .8030,   .5781 },
  { -.0764,   .8162,   .5727 },
  { -.0612,   .8553,   .5144 },
  { -.0081,   .8252,   .5648 },
  {  .1041,   .7494,   .6538 },
  {  .1742,   .7498,   .6383 },
  {  .0328,   .7453,   .6659 },
  {  .0481,   .7901,   .6111 },
  {  .0633,   .8304,   .5536 },
  {  .1874,   .7902,   .5835 },
  {  .0752,   .8659,   .4946 },
  {  .1314,   .8308,   .5408 },
  {  .1995,   .8271,   .5254 },
  {  .2102,   .8598,   .4653 },
  { -.7182,   .5721,   .3961 },
  { -.6690,   .6177,   .4135 },
  { -.6972,   .6563,   .2885 },
  { -.6454,   .7006,   .3042 },
  { -.6588,   .6609,   .3594 },
  { -.6161,   .6609,   .4285 },
  { -.5607,   .7006,   .4413 },
  { -.6628,   .7280,   .1753 },
  { -.6075,   .7718,   .1878 },
  { -.6276,   .7389,   .2454 },
  { -.6158,   .7857,   .0588 },
  { -.5576,   .8271,   .0702 },
  { -.5839,   .8015,   .1290 },
  { -.5484,   .8122,   .1990 },
  { -.4871,   .8480,   .2088 },
  { -.5100,   .8175,   .2674 },
  { -.5464,   .7435,   .3855 },
  { -.5001,   .7389,   .4516 },
  { -.5892,   .7435,   .3163 },
  { -.5294,   .7828,   .3272 },
  { -.4673,   .8175,   .3366 },
  { -.4396,   .7718,   .4593 },
  { -.4046,   .8480,   .3423 },
  { -.4232,   .8122,   .4016 },
  { -.3765,   .8015,   .4646 },
  { -.3121,   .8271,   .4673 },
  { -.5570,   .8284,  -.0588 },
  { -.4900,   .8697,  -.0588 },
  { -.5225,   .8526,   .0000 },
  { -.4875,   .8553,  -.1753 },
  { -.4172,   .8917,  -.1754 },
  { -.4547,   .8828,  -.1178 },
  { -.4195,   .9058,  -.0589 },
  { -.3469,   .9360,  -.0589 },
  { -.3827,   .9239,   .0000 },
  { -.4087,   .8659,  -.2885 },
  { -.3341,   .8968,  -.2900 },
  { -.3760,   .8964,  -.2347 },
  { -.3221,   .8598,  -.3961 },
  { -.2448,   .8850,  -.3961 },
  { -.2901,   .8929,  -.3442 },
  { -.2568,   .9219,  -.2900 },
  { -.1783,   .9407,  -.2885 },
  { -.2227,   .9462,  -.2347 },
  { -.3069,   .9444,  -.1178 },
  { -.2695,   .9612,  -.0589 },
  { -.3426,   .9227,  -.1769 },
  { -.2652,   .9478,  -.1769 },
  { -.1867,   .9666,  -.1754 },
  { -.1930,   .9794,  -.0589 },
  { -.2334,   .9724,   .0000 },
  { -.1084,   .9785,  -.1753 },
  { -.1511,   .9815,  -.1178 },
  { -.1148,   .9916,  -.0588 },
  { -.0364,   .9976,  -.0588 },
  { -.0785,   .9969,   .0000 },
  { -.2822,   .8685,   .4074 },
  { -.2337,   .8526,   .4673 },
  { -.3732,   .8848,   .2792 },
  { -.3291,   .8793,   .3442 },
  { -.2506,   .9048,   .3442 },
  { -.2181,   .9351,   .2792 },
  { -.1711,   .9239,   .3423 },
  { -.4550,   .8793,   .1404 },
  { -.4166,   .8848,   .2089 },
  { -.4906,   .8685,   .0702 },
  { -.4199,   .9048,   .0702 },
  { -.3473,   .9351,   .0703 },
  { -.2625,   .9416,   .2108 },
  { -.1830,   .9607,   .2089 },
  { -.3411,   .9161,   .2108 },
  { -.3060,   .9416,   .1405 },
  { -.2687,   .9607,   .0703 },
  { -.1488,   .9789,   .1404 },
  { -.1044,   .9724,   .2088 },
  { -.1921,   .9789,   .0702 },
  { -.1136,   .9910,   .0702 },
  { -.0351,   .9969,   .0702 },
  { -.1665,   .8697,   .4646 },
  { -.1025,   .9360,   .3366 },
  { -.1350,   .9058,   .4016 },
  { -.0980,   .8828,   .4593 },
  { -.0297,   .8917,   .4516 },
  { -.0337,   .9794,   .1990 },
  { -.0679,   .9612,   .2674 },
  {  .0363,   .9976,   .0588 },
  {  .0013,   .9916,   .1290 },
  {  .0378,   .9815,   .1878 },
  {  .1084,   .9785,   .1753 },
  {  .0734,   .9666,   .2454 },
  {  .0051,   .9227,   .3855 },
  {  .0418,   .8964,   .4413 },
  { -.0319,   .9444,   .3272 },
  {  .0397,   .9478,   .3163 },
  {  .1104,   .9462,   .3042 },
  {  .1100,   .8968,   .4285 },
  {  .1783,   .9407,   .2885 },
  {  .1445,   .9219,   .3594 },
  {  .1782,   .8929,   .4135 },
  {  .2448,   .8850,   .3961 },
  { -.7660,   .5063,   .3961 },
  { -.7596,   .5519,   .3442 },
  { -.8396,   .4603,   .2885 },
  { -.8311,   .5042,   .2347 },
  { -.7974,   .5291,   .2900 },
  { -.7497,   .5949,   .2900 },
  { -.7363,   .6346,   .2347 },
  { -.8972,   .4054,   .1753 },
  { -.8868,   .4470,   .1178 },
  { -.8617,   .4762,   .1754 },
  { -.9375,   .3428,   .0588 },
  { -.9239,   .3827,   .0000 },
  { -.9076,   .4156,   .0588 },
  { -.8718,   .4862,   .0589 },
  { -.8526,   .5225,   .0000 },
  { -.8309,   .5534,   .0589 },
  { -.7717,   .6109,   .1769 },
  { -.7192,   .6723,   .1754 },
  { -.8195,   .5451,   .1769 },
  { -.8034,   .5837,   .1178 },
  { -.7830,   .6192,   .0589 },
  { -.6991,   .7052,   .1178 },
  { -.7604,   .6494,   .0000 },
  { -.7319,   .6789,   .0589 },
  { -.6757,   .7348,   .0588 },
  { -.6494,   .7604,   .0000 },
  { -.9600,   .2737,  -.0588 },
  { -.9435,   .3052,  -.1290 },
  { -.9373,   .3414,  -.0702 },
  { -.9641,   .1993,  -.1753 },
  { -.9420,   .2289,  -.2454 },
  { -.9451,   .2674,  -.1878 },
  { -.9211,   .3347,  -.1990 },
  { -.8932,   .3616,  -.2674 },
  { -.8925,   .3998,  -.2088 },
  { -.9498,   .1211,  -.2885 },
  { -.9215,   .1474,  -.3594 },
  { -.9340,   .1874,  -.3042 },
  { -.9173,   .0406,  -.3961 },
  { -.8827,   .0658,  -.4653 },
  { -.9043,   .1065,  -.4135 },
  { -.8869,   .1725,  -.4285 },
  { -.8467,   .1960,  -.4946 },
  { -.8655,   .2372,  -.4413 },
  { -.8884,   .3221,  -.3272 },
  { -.8586,   .3867,  -.3366 },
  { -.9137,   .2552,  -.3163 },
  { -.8791,   .2803,  -.3855 },
  { -.8389,   .3038,  -.4516 },
  { -.8198,   .4083,  -.4016 },
  { -.8258,   .4483,  -.3423 },
  { -.7946,   .3225,  -.5144 },
  { -.8093,   .3660,  -.4593 },
  { -.7757,   .4271,  -.4646 },
  { -.7272,   .4429,  -.5245 },
  { -.7387,   .4857,  -.4673 },
  { -.6744,   .7350,  -.0702 },
  { -.6144,   .7859,  -.0702 },
  { -.7821,   .6192,  -.0703 },
  { -.7308,   .6790,  -.0702 },
  { -.6957,   .7045,  -.1404 },
  { -.7127,   .6696,  -.2089 },
  { -.6560,   .7253,  -.2088 },
  { -.8716,   .4852,  -.0702 },
  { -.8306,   .5524,  -.0703 },
  { -.9074,   .4143,  -.0702 },
  { -.8850,   .4440,  -.1404 },
  { -.8571,   .4709,  -.2089 },
  { -.7659,   .6075,  -.2108 },
  { -.7261,   .6283,  -.2792 },
  { -.8010,   .5820,  -.1405 },
  { -.8144,   .5406,  -.2108 },
  { -.8220,   .4964,  -.2792 },
  { -.7346,   .5847,  -.3442 },
  { -.6815,   .6468,  -.3423 },
  { -.7831,   .5180,  -.3442 },
  { -.7388,   .5368,  -.4074 },
  { -.6902,   .5525,  -.4673 },
  { -.5818,   .8030,  -.1290 },
  { -.6199,   .7377,  -.2674 },
  { -.6029,   .7726,  -.1990 },
  { -.5464,   .8162,  -.1878 },
  { -.5088,   .8252,  -.2454 },
  { -.6416,   .6535,  -.4016 },
  { -.6331,   .6970,  -.3366 },
  { -.6459,   .5547,  -.5245 },
  { -.6459,   .6057,  -.4646 },
  { -.5982,   .6566,  -.4593 },
  { -.5523,   .6560,  -.5144 },
  { -.5482,   .7040,  -.4516 },
  { -.5250,   .7901,  -.3163 },
  { -.4669,   .8304,  -.3042 },
  { -.5809,   .7453,  -.3272 },
  { -.5383,   .7494,  -.3855 },
  { -.4931,   .7498,  -.4413 },
  { -.4249,   .8308,  -.3594 },
  { -.4481,   .7447,  -.4946 },
  { -.4382,   .7902,  -.4285 },
  { -.3807,   .8271,  -.4135 },
  { -.3353,   .8192,  -.4653 },
  { -.7942,   .4454,   .4135 },
  { -.8396,   .3168,   .4413 },
  { -.8189,   .3817,   .4285 },
  { -.8321,   .4224,   .3594 },
  { -.8658,   .3974,   .3042 },
  { -.8699,   .1796,   .4593 },
  { -.8573,   .2473,   .4516 },
  { -.8831,   .0412,   .4673 },
  { -.8786,   .1104,   .4646 },
  { -.9032,   .1515,   .4016 },
  { -.9315,   .1227,   .3423 },
  { -.9219,   .1918,   .3366 },
  { -.8892,   .3306,   .3163 },
  { -.8966,   .3685,   .2454 },
  { -.8760,   .2899,   .3855 },
  { -.9081,   .2616,   .3272 },
  { -.9351,   .2325,   .2674 },
  { -.9218,   .3392,   .1878 },
  { -.9570,   .2012,   .2088 },
  { -.9419,   .2706,   .1990 },
  { -.9427,   .3077,   .1290 },
  { -.9590,   .2747,   .0702 },
  { -.8786,  -.1104,   .4646 },
  { -.8831,  -.0412,   .4673 },
  { -.8573,  -.2473,   .4516 },
  { -.8699,  -.1796,   .4593 },
  { -.9032,  -.1515,   .4016 },
  { -.9219,  -.1918,   .3366 },
  { -.9315,  -.1227,   .3423 },
  { -.8189,  -.3817,   .4285 },
  { -.8396,  -.3168,   .4413 },
  { -.7660,  -.5063,   .3961 },
  { -.7942,  -.4454,   .4135 },
  { -.8321,  -.4224,   .3594 },
  { -.8396,  -.4603,   .2885 },
  { -.8658,  -.3974,   .3042 },
  { -.9081,  -.2616,   .3272 },
  { -.9351,  -.2325,   .2674 },
  { -.8760,  -.2899,   .3855 },
  { -.8892,  -.3306,   .3163 },
  { -.8966,  -.3685,   .2454 },
  { -.9419,  -.2706,   .1990 },
  { -.9570,  -.2012,   .2088 },
  { -.8972,  -.4054,   .1753 },
  { -.9218,  -.3392,   .1878 },
  { -.9427,  -.3077,   .1290 },
  { -.9375,  -.3428,   .0588 },
  { -.9590,  -.2747,   .0702 },
  { -.9776,   .1982,   .0702 },
  { -.9724,   .2334,   .0000 },
  { -.9702,   .1228,   .2089 },
  { -.9769,   .1610,   .1404 },
  { -.9903,   .1198,   .0702 },
  { -.9967,   .0413,   .0703 },
  { -.9969,   .0785,   .0000 },
  { -.9380,   .0413,   .3442 },
  { -.9568,   .0815,   .2792 },
  { -.9132,   .0000,   .4074 },
  { -.9380,  -.0413,   .3442 },
  { -.9568,  -.0815,   .2792 },
  { -.9901,   .0000,   .1405 },
  { -.9967,  -.0413,   .0703 },
  { -.9767,   .0413,   .2108 },
  { -.9767,  -.0413,   .2108 },
  { -.9702,  -.1228,   .2089 },
  { -.9903,  -.1198,   .0702 },
  { -.9969,  -.0785,   .0000 },
  { -.9769,  -.1610,   .1404 },
  { -.9776,  -.1982,   .0702 },
  { -.9724,  -.2334,   .0000 },
  { -.9786,   .1973,  -.0588 },
  { -.9974,   .0407,  -.0589 },
  { -.9911,   .1191,  -.0589 },
  { -.9801,   .1596,  -.1178 },
  { -.9770,   .1212,  -.1754 },
  { -.9911,  -.1191,  -.0589 },
  { -.9974,  -.0407,  -.0589 },
  { -.9600,  -.2737,  -.0588 },
  { -.9786,  -.1973,  -.0588 },
  { -.9801,  -.1596,  -.1178 },
  { -.9641,  -.1993,  -.1753 },
  { -.9770,  -.1212,  -.1754 },
  { -.9834,   .0407,  -.1769 },
  { -.9687,   .0806,  -.2347 },
  { -.9930,   .0000,  -.1178 },
  { -.9834,  -.0407,  -.1769 },
  { -.9687,  -.0806,  -.2347 },
  { -.9562,   .0406,  -.2900 },
  { -.9498,  -.1211,  -.2885 },
  { -.9562,  -.0406,  -.2900 },
  { -.9389,   .0000,  -.3442 },
  { -.9173,  -.0406,  -.3961 },
  { -.7182,  -.5721,   .3961 },
  { -.7596,  -.5519,   .3442 },
  { -.6972,  -.6563,   .2885 },
  { -.7363,  -.6346,   .2347 },
  { -.7497,  -.5949,   .2900 },
  { -.7974,  -.5291,   .2900 },
  { -.8311,  -.5042,   .2347 },
  { -.6628,  -.7280,   .1753 },
  { -.6991,  -.7052,   .1178 },
  { -.7192,  -.6723,   .1754 },
  { -.6158,  -.7857,   .0588 },
  { -.6494,  -.7604,   .0000 },
  { -.6757,  -.7348,   .0588 },
  { -.7319,  -.6789,   .0589 },
  { -.7604,  -.6494,   .0000 },
  { -.7830,  -.6192,   .0589 },
  { -.8195,  -.5451,   .1769 },
  { -.8617,  -.4762,   .1754 },
  { -.7717,  -.6109,   .1769 },
  { -.8034,  -.5837,   .1178 },
  { -.8309,  -.5534,   .0589 },
  { -.8868,  -.4470,   .1178 },
  { -.8526,  -.5225,   .0000 },
  { -.8718,  -.4862,   .0589 },
  { -.9076,  -.4156,   .0588 },
  { -.9239,  -.3827,   .0000 },
  { -.5570,  -.8284,  -.0588 },
  { -.5818,  -.8030,  -.1290 },
  { -.6144,  -.7859,  -.0702 },
  { -.4875,  -.8553,  -.1753 },
  { -.5088,  -.8252,  -.2454 },
  { -.5464,  -.8162,  -.1878 },
  { -.6029,  -.7726,  -.1990 },
  { -.6199,  -.7377,  -.2674 },
  { -.6560,  -.7253,  -.2088 },
  { -.4087,  -.8659,  -.2885 },
  { -.4249,  -.8308,  -.3594 },
  { -.4669,  -.8304,  -.3042 },
  { -.3221,  -.8598,  -.3961 },
  { -.3353,  -.8192,  -.4653 },
  { -.3807,  -.8271,  -.4135 },
  { -.4382,  -.7902,  -.4285 },
  { -.4481,  -.7447,  -.4946 },
  { -.4931,  -.7498,  -.4413 },
  { -.5809,  -.7453,  -.3272 },
  { -.6331,  -.6970,  -.3366 },
  { -.5250,  -.7901,  -.3163 },
  { -.5383,  -.7494,  -.3855 },
  { -.5482,  -.7040,  -.4516 },
  { -.6416,  -.6535,  -.4016 },
  { -.6815,  -.6468,  -.3423 },
  { -.5523,  -.6560,  -.5144 },
  { -.5982,  -.6566,  -.4593 },
  { -.6459,  -.6057,  -.4646 },
  { -.6459,  -.5547,  -.5245 },
  { -.6902,  -.5525,  -.4673 },
  { -.9074,  -.4143,  -.0702 },
  { -.9373,  -.3414,  -.0702 },
  { -.8306,  -.5524,  -.0703 },
  { -.8716,  -.4852,  -.0702 },
  { -.8850,  -.4440,  -.1404 },
  { -.8571,  -.4709,  -.2089 },
  { -.8925,  -.3998,  -.2088 },
  { -.7308,  -.6790,  -.0702 },
  { -.7821,  -.6192,  -.0703 },
  { -.6744,  -.7350,  -.0702 },
  { -.6957,  -.7045,  -.1404 },
  { -.7127,  -.6696,  -.2089 },
  { -.8144,  -.5407,  -.2108 },
  { -.8220,  -.4964,  -.2792 },
  { -.8010,  -.5820,  -.1405 },
  { -.7659,  -.6075,  -.2108 },
  { -.7261,  -.6283,  -.2792 },
  { -.7831,  -.5180,  -.3442 },
  { -.8258,  -.4483,  -.3423 },
  { -.7346,  -.5847,  -.3442 },
  { -.7388,  -.5368,  -.4074 },
  { -.7387,  -.4857,  -.4673 },
  { -.9435,  -.3052,  -.1290 },
  { -.8932,  -.3616,  -.2674 },
  { -.9211,  -.3347,  -.1990 },
  { -.9451,  -.2674,  -.1878 },
  { -.9420,  -.2289,  -.2454 },
  { -.8198,  -.4083,  -.4016 },
  { -.8586,  -.3867,  -.3366 },
  { -.7272,  -.4429,  -.5245 },
  { -.7757,  -.4271,  -.4646 },
  { -.8093,  -.3660,  -.4593 },
  { -.7946,  -.3225,  -.5144 },
  { -.8389,  -.3038,  -.4516 },
  { -.9137,  -.2552,  -.3163 },
  { -.9340,  -.1874,  -.3042 },
  { -.8884,  -.3221,  -.3272 },
  { -.8791,  -.2803,  -.3855 },
  { -.8655,  -.2372,  -.4413 },
  { -.9215,  -.1474,  -.3594 },
  { -.8467,  -.1960,  -.4946 },
  { -.8869,  -.1725,  -.4285 },
  { -.9043,  -.1065,  -.4135 },
  { -.8827,  -.0658,  -.4653 },
  { -.6690,  -.6177,   .4135 },
  { -.5607,  -.7006,   .4413 },
  { -.6161,  -.6609,   .4285 },
  { -.6588,  -.6609,   .3594 },
  { -.6454,  -.7006,   .3042 },
  { -.4396,  -.7718,   .4593 },
  { -.5001,  -.7389,   .4516 },
  { -.3121,  -.8271,   .4673 },
  { -.3765,  -.8015,   .4646 },
  { -.4232,  -.8122,   .4016 },
  { -.4046,  -.8480,   .3423 },
  { -.4673,  -.8175,   .3366 },
  { -.5892,  -.7435,   .3163 },
  { -.6276,  -.7389,   .2454 },
  { -.5464,  -.7435,   .3855 },
  { -.5294,  -.7828,   .3272 },
  { -.5100,  -.8175,   .2674 },
  { -.6075,  -.7718,   .1878 },
  { -.4871,  -.8480,   .2088 },
  { -.5484,  -.8122,   .1990 },
  { -.5839,  -.8015,   .1290 },
  { -.5576,  -.8271,   .0702 },
  { -.1665,  -.8697,   .4646 },
  { -.2337,  -.8526,   .4673 },
  { -.0297,  -.8917,   .4516 },
  { -.0980,  -.8828,   .4593 },
  { -.1350,  -.9058,   .4016 },
  { -.1025,  -.9360,   .3366 },
  { -.1711,  -.9239,   .3423 },
  {  .1100,  -.8968,   .4285 },
  {  .0418,  -.8964,   .4413 },
  {  .2448,  -.8850,   .3961 },
  {  .1782,  -.8929,   .4135 },
  {  .1445,  -.9219,   .3594 },
  {  .1783,  -.9407,   .2885 },
  {  .1104,  -.9462,   .3042 },
  { -.0319,  -.9444,   .3272 },
  { -.0679,  -.9612,   .2674 },
  {  .0051,  -.9227,   .3855 },
  {  .0397,  -.9478,   .3163 },
  {  .0734,  -.9666,   .2454 },
  { -.0337,  -.9794,   .1990 },
  { -.1044,  -.9724,   .2088 },
  {  .1084,  -.9785,   .1753 },
  {  .0378,  -.9815,   .1878 },
  {  .0013,  -.9916,   .1290 },
  {  .0364,  -.9976,   .0588 },
  { -.0351,  -.9969,   .0702 },
  { -.4906,  -.8685,   .0702 },
  { -.5225,  -.8526,   .0000 },
  { -.4166,  -.8848,   .2089 },
  { -.4550,  -.8793,   .1404 },
  { -.4199,  -.9048,   .0702 },
  { -.3473,  -.9351,   .0703 },
  { -.3827,  -.9239,   .0000 },
  { -.3291,  -.8793,   .3442 },
  { -.3732,  -.8848,   .2792 },
  { -.2822,  -.8685,   .4074 },
  { -.2506,  -.9048,   .3442 },
  { -.2181,  -.9351,   .2792 },
  { -.3059,  -.9416,   .1405 },
  { -.2687,  -.9607,   .0703 },
  { -.3411,  -.9161,   .2108 },
  { -.2625,  -.9416,   .2108 },
  { -.1830,  -.9607,   .2089 },
  { -.1921,  -.9789,   .0702 },
  { -.2334,  -.9724,   .0000 },
  { -.1488,  -.9789,   .1404 },
  { -.1136,  -.9910,   .0702 },
  { -.0785,  -.9969,   .0000 },
  { -.4900,  -.8697,  -.0588 },
  { -.3469,  -.9360,  -.0589 },
  { -.4195,  -.9058,  -.0589 },
  { -.4547,  -.8828,  -.1178 },
  { -.4172,  -.8917,  -.1754 },
  { -.1930,  -.9794,  -.0589 },
  { -.2695,  -.9612,  -.0589 },
  { -.0363,  -.9976,  -.0588 },
  { -.1148,  -.9916,  -.0588 },
  { -.1511,  -.9815,  -.1178 },
  { -.1084,  -.9785,  -.1753 },
  { -.1867,  -.9666,  -.1754 },
  { -.3426,  -.9227,  -.1769 },
  { -.3760,  -.8964,  -.2347 },
  { -.3069,  -.9444,  -.1178 },
  { -.2652,  -.9478,  -.1769 },
  { -.2227,  -.9462,  -.2347 },
  { -.3341,  -.8968,  -.2900 },
  { -.1783,  -.9407,  -.2885 },
  { -.2568,  -.9219,  -.2900 },
  { -.2901,  -.8929,  -.3442 },
  { -.2448,  -.8850,  -.3961 },
  {  .3221,  -.8598,   .3961 },
  {  .2901,  -.8929,   .3442 },
  {  .4087,  -.8659,   .2885 },
  {  .3760,  -.8964,   .2347 },
  {  .3341,  -.8968,   .2900 },
  {  .2568,  -.9219,   .2900 },
  {  .2227,  -.9462,   .2347 },
  {  .4875,  -.8553,   .1753 },
  {  .4547,  -.8828,   .1178 },
  {  .4172,  -.8917,   .1754 },
  {  .5570,  -.8284,   .0588 },
  {  .5225,  -.8526,   .0000 },
  {  .4900,  -.8697,   .0588 },
  {  .4195,  -.9058,   .0589 },
  {  .3827,  -.9239,   .0000 },
  {  .3469,  -.9360,   .0589 },
  {  .2652,  -.9478,   .1769 },
  {  .1867,  -.9666,   .1754 },
  {  .3426,  -.9227,   .1769 },
  {  .3069,  -.9444,   .1178 },
  {  .2695,  -.9612,   .0589 },
  {  .1511,  -.9815,   .1178 },
  {  .2334,  -.9724,   .0000 },
  {  .1930,  -.9794,   .0589 },
  {  .1148,  -.9916,   .0588 },
  {  .0785,  -.9969,   .0000 },
  {  .6158,  -.7857,  -.0588 },
  {  .5839,  -.8015,  -.1290 },
  {  .5576,  -.8271,  -.0702 },
  {  .6628,  -.7280,  -.1753 },
  {  .6276,  -.7389,  -.2454 },
  {  .6075,  -.7718,  -.1878 },
  {  .5484,  -.8122,  -.1990 },
  {  .5100,  -.8175,  -.2674 },
  {  .4871,  -.8480,  -.2088 },
  {  .6972,  -.6563,  -.2885 },
  {  .6588,  -.6609,  -.3594 },
  {  .6454,  -.7006,  -.3042 },
  {  .7182,  -.5721,  -.3961 },
  {  .6755,  -.5721,  -.4653 },
  {  .6690,  -.6177,  -.4135 },
  {  .6161,  -.6609,  -.4285 },
  {  .5698,  -.6563,  -.4946 },
  {  .5607,  -.7006,  -.4413 },
  {  .5294,  -.7828,  -.3272 },
  {  .4673,  -.8175,  -.3366 },
  {  .5892,  -.7435,  -.3163 },
  {  .5464,  -.7435,  -.3855 },
  {  .5001,  -.7389,  -.4516 },
  {  .4232,  -.8122,  -.4016 },
  {  .4046,  -.8480,  -.3423 },
  {  .4532,  -.7280,  -.5144 },
  {  .4396,  -.7718,  -.4593 },
  {  .3765,  -.8015,  -.4646 },
  {  .3280,  -.7857,  -.5245 },
  {  .3121,  -.8271,  -.4673 },
  {  .1136,  -.9910,  -.0702 },
  {  .0351,  -.9969,  -.0702 },
  {  .2687,  -.9607,  -.0703 },
  {  .1921,  -.9789,  -.0702 },
  {  .1488,  -.9789,  -.1404 },
  {  .1830,  -.9607,  -.2089 },
  {  .1044,  -.9724,  -.2088 },
  {  .4199,  -.9048,  -.0702 },
  {  .3473,  -.9351,  -.0703 },
  {  .4906,  -.8685,  -.0702 },
  {  .4550,  -.8793,  -.1404 },
  {  .4166,  -.8848,  -.2089 },
  {  .2625,  -.9416,  -.2108 },
  {  .2181,  -.9351,  -.2792 },
  {  .3060,  -.9416,  -.1405 },
  {  .3411,  -.9161,  -.2108 },
  {  .3732,  -.8848,  -.2792 },
  {  .2506,  -.9048,  -.3442 },
  {  .1711,  -.9239,  -.3423 },
  {  .3291,  -.8793,  -.3442 },
  {  .2822,  -.8685,  -.4074 },
  {  .2337,  -.8526,  -.4673 },
  { -.0013,  -.9916,  -.1290 },
  {  .0679,  -.9612,  -.2674 },
  {  .0337,  -.9794,  -.1990 },
  { -.0378,  -.9815,  -.1878 },
  { -.0734,  -.9666,  -.2454 },
  {  .1350,  -.9058,  -.4016 },
  {  .1025,  -.9360,  -.3366 },
  {  .1965,  -.8284,  -.5245 },
  {  .1665,  -.8697,  -.4646 },
  {  .0980,  -.8828,  -.4593 },
  {  .0612,  -.8553,  -.5144 },
  {  .0297,  -.8917,  -.4516 },
  { -.0397,  -.9478,  -.3163 },
  { -.1104,  -.9462,  -.3042 },
  {  .0319,  -.9444,  -.3272 },
  { -.0051,  -.9227,  -.3855 },
  { -.0418,  -.8964,  -.4413 },
  { -.1445,  -.9219,  -.3594 },
  { -.0752,  -.8659,  -.4946 },
  { -.1100,  -.8968,  -.4285 },
  { -.1782,  -.8929,  -.4135 },
  { -.2102,  -.8598,  -.4653 },
  {  .3807,  -.8271,   .4135 },
  {  .4931,  -.7498,   .4413 },
  {  .4382,  -.7902,   .4285 },
  {  .4249,  -.8308,   .3594 },
  {  .4669,  -.8304,   .3042 },
  {  .5982,  -.6566,   .4593 },
  {  .5482,  -.7040,   .4516 },
  {  .6902,  -.5525,   .4673 },
  {  .6459,  -.6057,   .4646 },
  {  .6416,  -.6535,   .4016 },
  {  .6815,  -.6468,   .3423 },
  {  .6331,  -.6970,   .3366 },
  {  .5250,  -.7901,   .3163 },
  {  .5088,  -.8252,   .2454 },
  {  .5383,  -.7494,   .3855 },
  {  .5809,  -.7453,   .3272 },
  {  .6199,  -.7377,   .2674 },
  {  .5464,  -.8162,   .1878 },
  {  .6560,  -.7253,   .2088 },
  {  .6029,  -.7726,   .1990 },
  {  .5818,  -.8030,   .1290 },
  {  .6144,  -.7859,   .0702 },
  {  .7757,  -.4271,   .4646 },
  {  .7387,  -.4857,   .4673 },
  {  .8389,  -.3038,   .4516 },
  {  .8093,  -.3660,   .4593 },
  {  .8198,  -.4083,   .4016 },
  {  .8586,  -.3867,   .3366 },
  {  .8258,  -.4483,   .3423 },
  {  .8869,  -.1725,   .4285 },
  {  .8655,  -.2372,   .4413 },
  {  .9173,  -.0406,   .3961 },
  {  .9043,  -.1065,   .4135 },
  {  .9215,  -.1474,   .3594 },
  {  .9498,  -.1211,   .2885 },
  {  .9340,  -.1874,   .3042 },
  {  .8884,  -.3221,   .3272 },
  {  .8932,  -.3616,   .2674 },
  {  .8791,  -.2803,   .3855 },
  {  .9137,  -.2552,   .3163 },
  {  .9420,  -.2289,   .2454 },
  {  .9211,  -.3347,   .1990 },
  {  .8925,  -.3998,   .2088 },
  {  .9641,  -.1993,   .1753 },
  {  .9451,  -.2674,   .1878 },
  {  .9435,  -.3052,   .1290 },
  {  .9600,  -.2737,   .0588 },
  {  .9373,  -.3414,   .0702 },
  {  .6744,  -.7350,   .0702 },
  {  .6494,  -.7604,   .0000 },
  {  .7127,  -.6696,   .2089 },
  {  .6957,  -.7045,   .1404 },
  {  .7308,  -.6790,   .0702 },
  {  .7821,  -.6192,   .0703 },
  {  .7604,  -.6494,   .0000 },
  {  .7346,  -.5847,   .3442 },
  {  .7261,  -.6283,   .2792 },
  {  .7388,  -.5368,   .4074 },
  {  .7831,  -.5180,   .3442 },
  {  .8220,  -.4964,   .2792 },
  {  .8010,  -.5820,   .1405 },
  {  .8306,  -.5524,   .0703 },
  {  .7659,  -.6075,   .2108 },
  {  .8144,  -.5406,   .2108 },
  {  .8571,  -.4709,   .2089 },
  {  .8716,  -.4852,   .0702 },
  {  .8526,  -.5225,   .0000 },
  {  .8850,  -.4440,   .1404 },
  {  .9074,  -.4143,   .0702 },
  {  .9239,  -.3827,   .0000 },
  {  .6757,  -.7348,  -.0588 },
  {  .7830,  -.6192,  -.0589 },
  {  .7319,  -.6789,  -.0589 },
  {  .6991,  -.7052,  -.1178 },
  {  .7192,  -.6723,  -.1754 },
  {  .8718,  -.4862,  -.0589 },
  {  .8309,  -.5534,  -.0589 },
  {  .9375,  -.3428,  -.0588 },
  {  .9076,  -.4156,  -.0588 },
  {  .8868,  -.4470,  -.1178 },
  {  .8972,  -.4054,  -.1753 },
  {  .8617,  -.4762,  -.1754 },
  {  .7717,  -.6109,  -.1769 },
  {  .7363,  -.6346,  -.2347 },
  {  .8034,  -.5837,  -.1178 },
  {  .8195,  -.5451,  -.1769 },
  {  .8311,  -.5042,  -.2347 },
  {  .7497,  -.5949,  -.2900 },
  {  .8396,  -.4603,  -.2885 },
  {  .7974,  -.5291,  -.2900 },
  {  .7596,  -.5519,  -.3442 },
  {  .7660,  -.5063,  -.3961 },
  {  .9173,   .0406,   .3961 },
  {  .9389,   .0000,   .3442 },
  {  .9498,   .1211,   .2885 },
  {  .9687,   .0806,   .2347 },
  {  .9562,   .0406,   .2900 },
  {  .9562,  -.0406,   .2900 },
  {  .9687,  -.0806,   .2347 },
  {  .9641,   .1993,   .1753 },
  {  .9801,   .1596,   .1178 },
  {  .9770,   .1212,   .1754 },
  {  .9600,   .2737,   .0588 },
  {  .9724,   .2334,   .0000 },
  {  .9786,   .1973,   .0588 },
  {  .9911,   .1191,   .0589 },
  {  .9969,   .0785,   .0000 },
  {  .9974,   .0407,   .0589 },
  {  .9834,  -.0407,   .1769 },
  {  .9770,  -.1212,   .1754 },
  {  .9834,   .0407,   .1769 },
  {  .9930,   .0000,   .1178 },
  {  .9974,  -.0407,   .0589 },
  {  .9801,  -.1596,   .1178 },
  {  .9969,  -.0785,   .0000 },
  {  .9911,  -.1191,   .0589 },
  {  .9786,  -.1973,   .0588 },
  {  .9724,  -.2334,   .0000 },
  {  .9375,   .3428,  -.0588 },
  {  .9427,   .3077,  -.1290 },
  {  .9590,   .2747,  -.0702 },
  {  .8972,   .4054,  -.1753 },
  {  .8966,   .3685,  -.2454 },
  {  .9218,   .3392,  -.1878 },
  {  .9419,   .2706,  -.1990 },
  {  .9351,   .2325,  -.2674 },
  {  .9570,   .2012,  -.2088 },
  {  .8396,   .4603,  -.2885 },
  {  .8321,   .4224,  -.3594 },
  {  .8658,   .3974,  -.3042 },
  {  .7660,   .5063,  -.3961 },
  {  .7528,   .4657,  -.4653 },
  {  .7942,   .4454,  -.4135 },
  {  .8189,   .3817,  -.4285 },
  {  .8002,   .3391,  -.4946 },
  {  .8396,   .3168,  -.4413 },
  {  .9081,   .2616,  -.3272 },
  {  .9219,   .1918,  -.3366 },
  {  .8892,   .3306,  -.3163 },
  {  .8760,   .2899,  -.3855 },
  {  .8573,   .2473,  -.4516 },
  {  .9032,   .1515,  -.4016 },
  {  .9315,   .1227,  -.3423 },
  {  .8324,   .2061,  -.5144 },
  {  .8699,   .1796,  -.4593 },
  {  .8786,   .1104,  -.4646 },
  {  .8486,   .0691,  -.5245 },
  {  .8831,   .0412,  -.4673 },
  {  .9776,  -.1982,  -.0702 },
  {  .9590,  -.2747,  -.0702 },
  {  .9967,  -.0413,  -.0703 },
  {  .9903,  -.1198,  -.0702 },
  {  .9769,  -.1610,  -.1404 },
  {  .9702,  -.1228,  -.2089 },
  {  .9570,  -.2012,  -.2088 },
  {  .9903,   .1198,  -.0702 },
  {  .9967,   .0413,  -.0703 },
  {  .9776,   .1982,  -.0702 },
  {  .9769,   .1610,  -.1404 },
  {  .9702,   .1228,  -.2089 },
  {  .9767,  -.0413,  -.2108 },
  {  .9568,  -.0815,  -.2792 },
  {  .9901,   .0000,  -.1405 },
  {  .9767,   .0413,  -.2108 },
  {  .9568,   .0815,  -.2792 },
  {  .9380,  -.0413,  -.3442 },
  {  .9315,  -.1227,  -.3423 },
  {  .9380,   .0413,  -.3442 },
  {  .9132,   .0000,  -.4074 },
  {  .8831,  -.0412,  -.4673 },
  {  .9427,  -.3077,  -.1290 },
  {  .9351,  -.2325,  -.2674 },
  {  .9419,  -.2706,  -.1990 },
  {  .9218,  -.3392,  -.1878 },
  {  .8966,  -.3685,  -.2454 },
  {  .9032,  -.1515,  -.4016 },
  {  .9219,  -.1918,  -.3366 },
  {  .8486,  -.0691,  -.5245 },
  {  .8786,  -.1104,  -.4646 },
  {  .8699,  -.1796,  -.4593 },
  {  .8324,  -.2061,  -.5144 },
  {  .8573,  -.2473,  -.4516 },
  {  .8892,  -.3306,  -.3163 },
  {  .8658,  -.3974,  -.3042 },
  {  .9081,  -.2616,  -.3272 },
  {  .8760,  -.2899,  -.3855 },
  {  .8396,  -.3168,  -.4413 },
  {  .8321,  -.4224,  -.3594 },
  {  .8002,  -.3391,  -.4946 },
  {  .8189,  -.3817,  -.4285 },
  {  .7942,  -.4454,  -.4135 },
  {  .7528,  -.4657,  -.4653 },
  {  .9043,   .1065,   .4135 },
  {  .8655,   .2372,   .4413 },
  {  .8869,   .1725,   .4285 },
  {  .9215,   .1474,   .3594 },
  {  .9340,   .1874,   .3042 },
  {  .8093,   .3660,   .4593 },
  {  .8389,   .3038,   .4516 },
  {  .7387,   .4857,   .4673 },
  {  .7757,   .4271,   .4646 },
  {  .8198,   .4083,   .4016 },
  {  .8258,   .4483,   .3423 },
  {  .8586,   .3867,   .3366 },
  {  .9137,   .2552,   .3163 },
  {  .9420,   .2289,   .2454 },
  {  .8791,   .2803,   .3855 },
  {  .8884,   .3221,   .3272 },
  {  .8932,   .3616,   .2674 },
  {  .9451,   .2674,   .1878 },
  {  .8925,   .3998,   .2088 },
  {  .9211,   .3347,   .1990 },
  {  .9435,   .3052,   .1290 },
  {  .9373,   .3414,   .0702 },
  {  .6459,   .6057,   .4646 },
  {  .6902,   .5525,   .4673 },
  {  .5482,   .7040,   .4516 },
  {  .5982,   .6566,   .4593 },
  {  .6416,   .6535,   .4016 },
  {  .6331,   .6970,   .3366 },
  {  .6815,   .6468,   .3423 },
  {  .4382,   .7902,   .4285 },
  {  .4931,   .7498,   .4413 },
  {  .3221,   .8598,   .3961 },
  {  .3807,   .8271,   .4135 },
  {  .4249,   .8308,   .3594 },
  {  .4087,   .8659,   .2885 },
  {  .4669,   .8304,   .3042 },
  {  .5809,   .7453,   .3272 },
  {  .6199,   .7377,   .2674 },
  {  .5383,   .7494,   .3855 },
  {  .5250,   .7901,   .3163 },
  {  .5088,   .8252,   .2454 },
  {  .6029,   .7726,   .1990 },
  {  .6560,   .7253,   .2088 },
  {  .4875,   .8553,   .1753 },
  {  .5464,   .8162,   .1878 },
  {  .5818,   .8030,   .1290 },
  {  .5570,   .8284,   .0588 },
  {  .6144,   .7859,   .0702 },
  {  .9074,   .4143,   .0702 },
  {  .9239,   .3827,   .0000 },
  {  .8571,   .4709,   .2089 },
  {  .8850,   .4440,   .1404 },
  {  .8716,   .4852,   .0702 },
  {  .8306,   .5524,   .0703 },
  {  .8526,   .5225,   .0000 },
  {  .7831,   .5180,   .3442 },
  {  .8220,   .4964,   .2792 },
  {  .7388,   .5368,   .4074 },
  {  .7346,   .5847,   .3442 },
  {  .7261,   .6283,   .2792 },
  {  .8010,   .5820,   .1405 },
  {  .7821,   .6192,   .0703 },
  {  .8144,   .5406,   .2108 },
  {  .7659,   .6075,   .2108 },
  {  .7127,   .6696,   .2089 },
  {  .7308,   .6790,   .0702 },
  {  .7604,   .6494,   .0000 },
  {  .6957,   .7045,   .1404 },
  {  .6744,   .7350,   .0702 },
  {  .6494,   .7604,   .0000 },
  {  .9076,   .4156,  -.0588 },
  {  .8309,   .5534,  -.0589 },
  {  .8718,   .4862,  -.0589 },
  {  .8868,   .4470,  -.1178 },
  {  .8617,   .4762,  -.1754 },
  {  .7319,   .6789,  -.0589 },
  {  .7830,   .6192,  -.0589 },
  {  .6158,   .7857,  -.0588 },
  {  .6757,   .7348,  -.0588 },
  {  .6991,   .7052,  -.1178 },
  {  .6628,   .7280,  -.1753 },
  {  .7192,   .6723,  -.1754 },
  {  .8195,   .5451,  -.1769 },
  {  .8311,   .5042,  -.2347 },
  {  .8034,   .5837,  -.1178 },
  {  .7717,   .6109,  -.1769 },
  {  .7363,   .6346,  -.2347 },
  {  .7974,   .5291,  -.2900 },
  {  .6972,   .6563,  -.2885 },
  {  .7497,   .5949,  -.2900 },
  {  .7596,   .5519,  -.3442 },
  {  .7182,   .5721,  -.3961 },
  {  .2901,   .8929,   .3442 },
  {  .2227,   .9462,   .2347 },
  {  .2568,   .9219,   .2900 },
  {  .3341,   .8968,   .2900 },
  {  .3760,   .8964,   .2347 },
  {  .1511,   .9815,   .1178 },
  {  .1867,   .9666,   .1754 },
  {  .0785,   .9969,   .0000 },
  {  .1148,   .9916,   .0588 },
  {  .1930,   .9794,   .0589 },
  {  .2334,   .9724,   .0000 },
  {  .2695,   .9612,   .0589 },
  {  .3426,   .9227,   .1769 },
  {  .4172,   .8917,   .1754 },
  {  .2652,   .9478,   .1769 },
  {  .3069,   .9444,   .1178 },
  {  .3469,   .9360,   .0589 },
  {  .4547,   .8828,   .1178 },
  {  .3827,   .9239,   .0000 },
  {  .4195,   .9058,   .0589 },
  {  .4900,   .8697,   .0588 },
  {  .5225,   .8526,   .0000 },
  { -.0013,   .9916,  -.1290 },
  {  .0351,   .9969,  -.0702 },
  { -.0734,   .9666,  -.2454 },
  { -.0378,   .9815,  -.1878 },
  {  .0337,   .9794,  -.1990 },
  {  .0679,   .9612,  -.2674 },
  {  .1044,   .9724,  -.2088 },
  { -.1446,   .9219,  -.3594 },
  { -.1104,   .9462,  -.3042 },
  { -.2102,   .8598,  -.4653 },
  { -.1782,   .8929,  -.4135 },
  { -.1100,   .8968,  -.4285 },
  { -.0752,   .8659,  -.4946 },
  { -.0418,   .8964,  -.4413 },
  {  .0319,   .9444,  -.3272 },
  {  .1025,   .9360,  -.3366 },
  { -.0397,   .9478,  -.3163 },
  { -.0051,   .9227,  -.3855 },
  {  .0297,   .8917,  -.4516 },
  {  .1350,   .9058,  -.4016 },
  {  .1711,   .9239,  -.3423 },
  {  .0612,   .8553,  -.5144 },
  {  .0980,   .8828,  -.4593 },
  {  .1665,   .8697,  -.4646 },
  {  .1965,   .8284,  -.5245 },
  {  .2337,   .8526,  -.4673 },
  {  .4906,   .8685,  -.0702 },
  {  .5576,   .8271,  -.0702 },
  {  .3473,   .9351,  -.0703 },
  {  .4199,   .9048,  -.0702 },
  {  .4550,   .8793,  -.1404 },
  {  .4166,   .8848,  -.2089 },
  {  .4871,   .8480,  -.2088 },
  {  .1921,   .9789,  -.0702 },
  {  .2687,   .9607,  -.0703 },
  {  .1136,   .9910,  -.0702 },
  {  .1488,   .9789,  -.1404 },
  {  .1830,   .9607,  -.2089 },
  {  .3411,   .9161,  -.2108 },
  {  .3732,   .8848,  -.2792 },
  {  .3060,   .9416,  -.1405 },
  {  .2625,   .9416,  -.2108 },
  {  .2181,   .9351,  -.2792 },
  {  .3291,   .8793,  -.3442 },
  {  .4046,   .8480,  -.3423 },
  {  .2506,   .9048,  -.3442 },
  {  .2822,   .8685,  -.4074 },
  {  .3121,   .8271,  -.4673 },
  {  .5839,   .8015,  -.1290 },
  {  .5100,   .8175,  -.2674 },
  {  .5484,   .8122,  -.1990 },
  {  .6075,   .7718,  -.1878 },
  {  .6276,   .7389,  -.2454 },
  {  .4232,   .8122,  -.4016 },
  {  .4673,   .8175,  -.3366 },
  {  .3280,   .7857,  -.5245 },
  {  .3765,   .8015,  -.4646 },
  {  .4396,   .7718,  -.4593 },
  {  .4532,   .7280,  -.5144 },
  {  .5001,   .7389,  -.4516 },
  {  .5892,   .7435,  -.3163 },
  {  .6454,   .7006,  -.3042 },
  {  .5294,   .7828,  -.3272 },
  {  .5464,   .7435,  -.3855 },
  {  .5607,   .7006,  -.4413 },
  {  .6588,   .6609,  -.3594 },
  {  .5698,   .6563,  -.4946 },
  {  .6161,   .6609,  -.4285 },
  {  .6690,   .6177,  -.4135 },
  {  .6755,   .5721,  -.4653 },
  { -.3247,   .7864,  -.5254 },
  { -.2662,   .8192,  -.5080 },
  { -.4369,   .7090,  -.5536 },
  { -.3820,   .7494,  -.5408 },
  { -.3129,   .7494,  -.5835 },
  { -.2998,   .7090,  -.6383 },
  { -.2420,   .7447,  -.6220 },
  { -.5416,   .6155,  -.5727 },
  { -.4915,   .6628,  -.5648 },
  { -.6334,   .5112,  -.5809 },
  { -.5892,   .5645,  -.5781 },
  { -.5283,   .5711,  -.6283 },
  { -.5126,   .5241,  -.6801 },
  { -.4644,   .5745,  -.6740 },
  { -.3563,   .6675,  -.6538 },
  { -.2854,   .6628,  -.6922 },
  { -.4255,   .6675,  -.6111 },
  { -.4115,   .6223,  -.6659 },
  { -.3952,   .5745,  -.7168 },
  { -.2700,   .6155,  -.7405 },
  { -.2131,   .6560,  -.7240 },
  { -.3791,   .5241,  -.7626 },
  { -.3257,   .5711,  -.7535 },
  { -.2536,   .5645,  -.7855 },
  { -.2363,   .5112,  -.8263 },
  { -.1802,   .5547,  -.8123 },
  { -.7190,   .3859,  -.5781 },
  { -.6819,   .4445,  -.5809 },
  { -.7823,   .2627,  -.5648 },
  { -.7527,   .3249,  -.5727 },
  { -.7064,   .3259,  -.6283 },
  { -.6898,   .2642,  -.6740 },
  { -.6569,   .3255,  -.6801 },
  { -.8308,   .1317,  -.5408 },
  { -.8093,   .1964,  -.5536 },
  { -.8614,   .0000,  -.5080 },
  { -.8483,   .0658,  -.5254 },
  { -.8094,   .0660,  -.5835 },
  { -.7831,   .0000,  -.6220 },
  { -.7669,   .0660,  -.6383 },
  { -.7190,   .1991,  -.6659 },
  { -.6685,   .1983,  -.7168 },
  { -.7663,   .1984,  -.6111 },
  { -.7449,   .1326,  -.6538 },
  { -.7186,   .0666,  -.6922 },
  { -.6438,   .1333,  -.7535 },
  { -.6156,   .1986,  -.7626 },
  { -.6898,   .0000,  -.7240 },
  { -.6688,   .0666,  -.7405 },
  { -.6152,   .0667,  -.7855 },
  { -.5833,   .0000,  -.8123 },
  { -.5592,   .0667,  -.8263 },
  { -.2880,   .4543,  -.8430 },
  { -.2146,   .4445,  -.8697 },
  { -.4300,   .4642,  -.7743 },
  { -.3604,   .4609,  -.8110 },
  { -.3387,   .3941,  -.8543 },
  { -.3872,   .3323,  -.8600 },
  { -.3146,   .3255,  -.8917 },
  { -.5642,   .4609,  -.6851 },
  { -.5003,   .4642,  -.7309 },
  { -.6252,   .4543,  -.6346 },
  { -.6127,   .3941,  -.6851 },
  { -.5961,   .3323,  -.7309 },
  { -.4599,   .3341,  -.8227 },
  { -.4357,   .2655,  -.8600 },
  { -.4816,   .4010,  -.7793 },
  { -.5302,   .3341,  -.7793 },
  { -.5744,   .2655,  -.7743 },
  { -.4795,   .2004,  -.8543 },
  { -.4068,   .1986,  -.8917 },
  { -.5497,   .2004,  -.8110 },
  { -.5210,   .1335,  -.8430 },
  { -.4890,   .0667,  -.8697 },
  { -.1955,   .3859,  -.9016 },
  { -.1439,   .4429,  -.8850 },
  { -.2944,   .2642,  -.9185 },
  { -.2460,   .3259,  -.9128 },
  { -.1756,   .3249,  -.9293 },
  { -.1553,   .2627,  -.9523 },
  { -.1048,   .3225,  -.9407 },
  { -.3860,   .1333,  -.9128 },
  { -.3422,   .1983,  -.9185 },
  { -.4657,   .0000,  -.8850 },
  { -.4274,   .0667,  -.9016 },
  { -.3632,   .0666,  -.9293 },
  { -.3391,   .0000,  -.9407 },
  { -.2978,   .0666,  -.9523 },
  { -.2039,   .1984,  -.9587 },
  { -.1333,   .1964,  -.9714 },
  { -.2740,   .1991,  -.9409 },
  { -.2517,   .1326,  -.9587 },
  { -.2280,   .0660,  -.9714 },
  { -.1122,   .1317,  -.9849 },
  { -.0637,   .1960,  -.9785 },
  { -.2061,   .0000,  -.9785 },
  { -.1600,   .0660,  -.9849 },
  { -.0906,   .0658,  -.9937 },
  { -.0691,   .0000,  -.9976 },
  { -.0214,   .0658,  -.9976 },
  { -.8483,  -.0658,  -.5254 },
  { -.8093,  -.1964,  -.5536 },
  { -.8308,  -.1317,  -.5408 },
  { -.8094,  -.0660,  -.5835 },
  { -.7669,  -.0660,  -.6383 },
  { -.7527,  -.3249,  -.5727 },
  { -.7823,  -.2627,  -.5648 },
  { -.6819,  -.4445,  -.5809 },
  { -.7190,  -.3859,  -.5781 },
  { -.7064,  -.3259,  -.6283 },
  { -.6569,  -.3255,  -.6801 },
  { -.6898,  -.2642,  -.6740 },
  { -.7449,  -.1326,  -.6538 },
  { -.7186,  -.0666,  -.6922 },
  { -.7663,  -.1984,  -.6111 },
  { -.7190,  -.1991,  -.6659 },
  { -.6685,  -.1983,  -.7168 },
  { -.6688,  -.0666,  -.7405 },
  { -.6156,  -.1986,  -.7626 },
  { -.6438,  -.1333,  -.7535 },
  { -.6152,  -.0667,  -.7855 },
  { -.5592,  -.0667,  -.8263 },
  { -.5892,  -.5645,  -.5781 },
  { -.6334,  -.5112,  -.5809 },
  { -.4915,  -.6628,  -.5648 },
  { -.5416,  -.6155,  -.5727 },
  { -.5283,  -.5711,  -.6283 },
  { -.4644,  -.5745,  -.6740 },
  { -.5126,  -.5241,  -.6801 },
  { -.3820,  -.7494,  -.5408 },
  { -.4369,  -.7090,  -.5536 },
  { -.2662,  -.8192,  -.5080 },
  { -.3247,  -.7864,  -.5254 },
  { -.3129,  -.7494,  -.5835 },
  { -.2420,  -.7447,  -.6220 },
  { -.2998,  -.7090,  -.6383 },
  { -.4115,  -.6223,  -.6659 },
  { -.3952,  -.5745,  -.7168 },
  { -.4255,  -.6675,  -.6111 },
  { -.3563,  -.6675,  -.6538 },
  { -.2854,  -.6628,  -.6922 },
  { -.3257,  -.5711,  -.7535 },
  { -.3791,  -.5241,  -.7626 },
  { -.2131,  -.6560,  -.7240 },
  { -.2700,  -.6155,  -.7405 },
  { -.2536,  -.5645,  -.7855 },
  { -.1802,  -.5547,  -.8123 },
  { -.2363,  -.5112,  -.8263 },
  { -.5210,  -.1335,  -.8430 },
  { -.4890,  -.0667,  -.8697 },
  { -.5744,  -.2655,  -.7743 },
  { -.5497,  -.2004,  -.8110 },
  { -.4795,  -.2004,  -.8543 },
  { -.4357,  -.2655,  -.8600 },
  { -.4068,  -.1986,  -.8917 },
  { -.6127,  -.3941,  -.6851 },
  { -.5961,  -.3323,  -.7309 },
  { -.6252,  -.4543,  -.6346 },
  { -.5642,  -.4609,  -.6851 },
  { -.5003,  -.4642,  -.7309 },
  { -.4599,  -.3341,  -.8227 },
  { -.3872,  -.3323,  -.8600 },
  { -.5302,  -.3341,  -.7793 },
  { -.4816,  -.4010,  -.7793 },
  { -.4300,  -.4642,  -.7743 },
  { -.3387,  -.3941,  -.8543 },
  { -.3146,  -.3255,  -.8917 },
  { -.3604,  -.4609,  -.8110 },
  { -.2880,  -.4543,  -.8430 },
  { -.2146,  -.4445,  -.8697 },
  { -.4274,  -.0667,  -.9016 },
  { -.3422,  -.1983,  -.9185 },
  { -.3860,  -.1333,  -.9128 },
  { -.3632,  -.0666,  -.9293 },
  { -.2978,  -.0666,  -.9523 },
  { -.2460,  -.3259,  -.9128 },
  { -.2944,  -.2642,  -.9185 },
  { -.1439,  -.4429,  -.8850 },
  { -.1955,  -.3859,  -.9016 },
  { -.1756,  -.3249,  -.9293 },
  { -.1048,  -.3225,  -.9407 },
  { -.1553,  -.2627,  -.9523 },
  { -.2517,  -.1326,  -.9587 },
  { -.2280,  -.0660,  -.9714 },
  { -.2740,  -.1991,  -.9409 },
  { -.2039,  -.1984,  -.9587 },
  { -.1333,  -.1964,  -.9714 },
  { -.1600,  -.0660,  -.9849 },
  { -.0637,  -.1960,  -.9785 },
  { -.1122,  -.1317,  -.9849 },
  { -.0906,  -.0658,  -.9937 },
  { -.0214,  -.0658,  -.9976 },
  { -.1995,  -.8271,  -.5254 },
  { -.0633,  -.8304,  -.5536 },
  { -.1314,  -.8308,  -.5408 },
  { -.1874,  -.7902,  -.5835 },
  { -.1742,  -.7498,  -.6383 },
  {  .0764,  -.8162,  -.5727 },
  {  .0081,  -.8252,  -.5648 },
  {  .2120,  -.7859,  -.5809 },
  {  .1449,  -.8030,  -.5781 },
  {  .0917,  -.7726,  -.6283 },
  {  .1066,  -.7253,  -.6801 },
  {  .0380,  -.7377,  -.6740 },
  { -.1041,  -.7494,  -.6538 },
  { -.1587,  -.7040,  -.6922 },
  { -.0481,  -.7901,  -.6111 },
  { -.0328,  -.7453,  -.6659 },
  { -.0179,  -.6970,  -.7168 },
  { -.1433,  -.6566,  -.7405 },
  { -.0014,  -.6468,  -.7626 },
  { -.0722,  -.6535,  -.7535 },
  { -.1267,  -.6057,  -.7855 },
  { -.1093,  -.5525,  -.8263 },
  {  .3548,  -.7348,  -.5781 },
  {  .2904,  -.7604,  -.5809 },
  {  .4785,  -.6723,  -.5648 },
  {  .4180,  -.7052,  -.5727 },
  {  .3799,  -.6789,  -.6283 },
  {  .4028,  -.6192,  -.6740 },
  {  .3401,  -.6494,  -.6801 },
  {  .5947,  -.5949,  -.5408 },
  {  .5393,  -.6346,  -.5536 },
  {  .6969,  -.5063,  -.5080 },
  {  .6476,  -.5519,  -.5254 },
  {  .6160,  -.5291,  -.5835 },
  {  .6335,  -.4603,  -.6220 },
  {  .5816,  -.5042,  -.6383 },
  {  .4647,  -.5837,  -.6659 },
  {  .4242,  -.5534,  -.7168 },
  {  .5034,  -.6109,  -.6111 },
  {  .5247,  -.5451,  -.6538 },
  {  .5422,  -.4762,  -.6922 },
  {  .4425,  -.4862,  -.7535 },
  {  .3813,  -.5225,  -.7626 },
  {  .5580,  -.4054,  -.7240 },
  {  .5019,  -.4470,  -.7405 },
  {  .4585,  -.4156,  -.7855 },
  {  .4719,  -.3428,  -.8123 },
  {  .4132,  -.3827,  -.8263 },
  { -.0340,  -.5368,  -.8430 },
  { -.0876,  -.4857,  -.8697 },
  {  .0750,  -.6283,  -.7743 },
  {  .0207,  -.5847,  -.8110 },
  {  .0424,  -.5180,  -.8543 },
  {  .1179,  -.4964,  -.8600 },
  {  .0632,  -.4483,  -.8917 },
  {  .1855,  -.7045,  -.6851 },
  {  .1319,  -.6696,  -.7309 },
  {  .2388,  -.7350,  -.6346 },
  {  .2640,  -.6790,  -.6851 },
  {  .2869,  -.6192,  -.7309 },
  {  .1757,  -.5406,  -.8227 },
  {  .1964,  -.4709,  -.8600 },
  {  .1540,  -.6075,  -.7793 },
  {  .2325,  -.5820,  -.7793 },
  {  .3086,  -.5524,  -.7743 },
  {  .2702,  -.4440,  -.8543 },
  {  .2124,  -.3998,  -.8917 },
  {  .3269,  -.4852,  -.8110 },
  {  .3430,  -.4143,  -.8430 },
  {  .3564,  -.3414,  -.8697 },
  { -.0686,  -.4271,  -.9016 },
  {  .0829,  -.3867,  -.9185 },
  {  .0075,  -.4083,  -.9128 },
  { -.0489,  -.3660,  -.9293 },
  { -.0287,  -.3038,  -.9523 },
  {  .2340,  -.3347,  -.9128 },
  {  .1603,  -.3616,  -.9185 },
  {  .3767,  -.2737,  -.8850 },
  {  .3066,  -.3052,  -.9016 },
  {  .2547,  -.2674,  -.9293 },
  {  .2744,  -.1993,  -.9407 },
  {  .2018,  -.2289,  -.9523 },
  {  .0483,  -.2803,  -.9587 },
  { -.0076,  -.2372,  -.9714 },
  {  .1047,  -.3221,  -.9409 },
  {  .1257,  -.2552,  -.9587 },
  {  .1456,  -.1874,  -.9714 },
  {  .0133,  -.1725,  -.9849 },
  {  .1667,  -.1211,  -.9785 },
  {  .0906,  -.1474,  -.9849 },
  {  .0346,  -.1065,  -.9937 },
  {  .0559,  -.0406,  -.9976 },
  {  .7250,  -.4454,  -.5254 },
  {  .7702,  -.3168,  -.5536 },
  {  .7495,  -.3817,  -.5408 },
  {  .6936,  -.4224,  -.5835 },
  {  .6593,  -.3974,  -.6383 },
  {  .7999,  -.1796,  -.5727 },
  {  .7873,  -.2473,  -.5648 },
  {  .8129,  -.0412,  -.5809 },
  {  .8085,  -.1104,  -.5781 },
  {  .7631,  -.1515,  -.6283 },
  {  .7227,  -.1227,  -.6801 },
  {  .7134,  -.1918,  -.6740 },
  {  .6806,  -.3306,  -.6538 },
  {  .6205,  -.3685,  -.6922 },
  {  .7366,  -.2899,  -.6111 },
  {  .6987,  -.2616,  -.6659 },
  {  .6574,  -.2325,  -.7168 },
  {  .5802,  -.3392,  -.7405 },
  {  .6148,  -.2012,  -.7626 },
  {  .5992,  -.2706,  -.7535 },
  {  .5369,  -.3077,  -.7855 },
  {  .4916,  -.2747,  -.8263 },
  {  .8085,   .1104,  -.5781 },
  {  .8129,   .0412,  -.5809 },
  {  .7873,   .2473,  -.5648 },
  {  .7999,   .1796,  -.5727 },
  {  .7631,   .1515,  -.6283 },
  {  .7134,   .1918,  -.6740 },
  {  .7227,   .1227,  -.6801 },
  {  .7495,   .3817,  -.5408 },
  {  .7702,   .3168,  -.5536 },
  {  .6969,   .5063,  -.5080 },
  {  .7250,   .4454,  -.5254 },
  {  .6936,   .4224,  -.5835 },
  {  .6335,   .4603,  -.6220 },
  {  .6593,   .3974,  -.6383 },
  {  .6987,   .2616,  -.6659 },
  {  .6574,   .2325,  -.7168 },
  {  .7366,   .2899,  -.6111 },
  {  .6806,   .3306,  -.6538 },
  {  .6205,   .3685,  -.6922 },
  {  .5992,   .2706,  -.7535 },
  {  .6148,   .2012,  -.7626 },
  {  .5580,   .4054,  -.7240 },
  {  .5802,   .3392,  -.7405 },
  {  .5369,   .3077,  -.7855 },
  {  .4719,   .3428,  -.8123 },
  {  .4916,   .2747,  -.8263 },
  {  .5000,  -.1982,  -.8430 },
  {  .4349,  -.2334,  -.8697 },
  {  .6207,  -.1228,  -.7743 },
  {  .5625,  -.1610,  -.8110 },
  {  .5057,  -.1198,  -.8543 },
  {  .5086,  -.0413,  -.8600 },
  {  .4458,  -.0785,  -.8917 },
  {  .7273,  -.0413,  -.6851 },
  {  .6776,  -.0815,  -.7309 },
  {  .7728,   .0000,  -.6346 },
  {  .7273,   .0413,  -.6851 },
  {  .6776,   .0815,  -.7309 },
  {  .5685,   .0000,  -.8227 },
  {  .5086,   .0413,  -.8600 },
  {  .6253,  -.0413,  -.7793 },
  {  .6253,   .0413,  -.7793 },
  {  .6207,   .1228,  -.7743 },
  {  .5057,   .1198,  -.8543 },
  {  .4458,   .0785,  -.8917 },
  {  .5625,   .1610,  -.8110 },
  {  .5000,   .1982,  -.8430 },
  {  .4349,   .2334,  -.8697 },
  {  .3850,  -.1973,  -.9016 },
  {  .3934,  -.0407,  -.9185 },
  {  .3906,  -.1191,  -.9128 },
  {  .3330,  -.1596,  -.9293 },
  {  .2801,  -.1212,  -.9523 },
  {  .3906,   .1191,  -.9128 },
  {  .3934,   .0407,  -.9185 },
  {  .3767,   .2737,  -.8850 },
  {  .3850,   .1973,  -.9016 },
  {  .3330,   .1596,  -.9293 },
  {  .2744,   .1993,  -.9407 },
  {  .2801,   .1212,  -.9523 },
  {  .2815,  -.0407,  -.9587 },
  {  .2233,  -.0806,  -.9714 },
  {  .3387,   .0000,  -.9409 },
  {  .2815,   .0407,  -.9587 },
  {  .2233,   .0806,  -.9714 },
  {  .1682,  -.0406,  -.9849 },
  {  .1667,   .1211,  -.9785 },
  {  .1682,   .0406,  -.9849 },
  {  .1120,   .0000,  -.9937 },
  {  .0559,   .0406,  -.9976 },
  {  .6476,   .5519,  -.5254 },
  {  .5393,   .6346,  -.5536 },
  {  .5947,   .5949,  -.5408 },
  {  .6160,   .5291,  -.5835 },
  {  .5816,   .5042,  -.6383 },
  {  .4180,   .7052,  -.5727 },
  {  .4785,   .6723,  -.5648 },
  {  .2904,   .7604,  -.5809 },
  {  .3548,   .7348,  -.5781 },
  {  .3799,   .6789,  -.6283 },
  {  .3401,   .6494,  -.6801 },
  {  .4028,   .6192,  -.6740 },
  {  .5247,   .5451,  -.6538 },
  {  .5422,   .4762,  -.6922 },
  {  .5034,   .6109,  -.6111 },
  {  .4647,   .5837,  -.6659 },
  {  .4242,   .5534,  -.7168 },
  {  .5019,   .4470,  -.7405 },
  {  .3813,   .5225,  -.7626 },
  {  .4425,   .4862,  -.7535 },
  {  .4585,   .4156,  -.7855 },
  {  .4132,   .3827,  -.8263 },
  {  .1448,   .8030,  -.5781 },
  {  .2120,   .7859,  -.5809 },
  {  .0081,   .8252,  -.5648 },
  {  .0764,   .8162,  -.5727 },
  {  .0917,   .7726,  -.6283 },
  {  .0380,   .7377,  -.6740 },
  {  .1066,   .7253,  -.6801 },
  { -.1314,   .8308,  -.5408 },
  { -.0633,   .8304,  -.5536 },
  { -.1995,   .8271,  -.5254 },
  { -.1874,   .7902,  -.5835 },
  { -.1742,   .7498,  -.6383 },
  { -.0328,   .7453,  -.6659 },
  { -.0179,   .6970,  -.7168 },
  { -.0481,   .7901,  -.6111 },
  { -.1041,   .7494,  -.6538 },
  { -.1587,   .7040,  -.6922 },
  { -.0722,   .6535,  -.7535 },
  { -.0014,   .6468,  -.7626 },
  { -.1433,   .6566,  -.7405 },
  { -.1267,   .6057,  -.7855 },
  { -.1093,   .5525,  -.8263 },
  {  .3430,   .4143,  -.8430 },
  {  .3564,   .3414,  -.8697 },
  {  .3086,   .5524,  -.7743 },
  {  .3269,   .4852,  -.8110 },
  {  .2702,   .4440,  -.8543 },
  {  .1964,   .4709,  -.8600 },
  {  .2124,   .3998,  -.8917 },
  {  .2640,   .6790,  -.6851 },
  {  .2869,   .6192,  -.7309 },
  {  .2388,   .7350,  -.6346 },
  {  .1855,   .7045,  -.6851 },
  {  .1319,   .6696,  -.7309 },
  {  .1757,   .5406,  -.8227 },
  {  .1179,   .4964,  -.8600 },
  {  .2325,   .5820,  -.7793 },
  {  .1540,   .6075,  -.7793 },
  {  .0750,   .6283,  -.7743 },
  {  .0424,   .5180,  -.8543 },
  {  .0632,   .4483,  -.8917 },
  {  .0207,   .5847,  -.8110 },
  { -.0340,   .5368,  -.8430 },
  { -.0876,   .4857,  -.8697 },
  {  .3066,   .3052,  -.9016 },
  {  .1603,   .3616,  -.9185 },
  {  .2340,   .3347,  -.9128 },
  {  .2547,   .2674,  -.9293 },
  {  .2018,   .2289,  -.9523 },
  {  .0075,   .4083,  -.9128 },
  {  .0829,   .3867,  -.9185 },
  { -.0686,   .4271,  -.9016 },
  { -.0489,   .3660,  -.9293 },
  { -.0287,   .3038,  -.9523 },
  {  .1257,   .2552,  -.9587 },
  {  .1456,   .1874,  -.9714 },
  {  .1047,   .3221,  -.9409 },
  {  .0483,   .2803,  -.9587 },
  { -.0076,   .2372,  -.9714 },
  {  .0906,   .1474,  -.9849 },
  {  .0133,   .1725,  -.9849 },
  {  .0346,   .1065,  -.9937 }
} ;


/*---------------------------------------------------------------------------*/
IC_FACE ic4_faces[5120] =
{
  {{    1,   643,   645 }},
  {{  643,   163,   644 }},
  {{  645,   643,   644 }},
  {{  645,   644,   165 }},
  {{  163,   646,   648 }},
  {{  646,    43,   647 }},
  {{  648,   646,   647 }},
  {{  648,   647,   164 }},
  {{  165,   644,   649 }},
  {{  644,   163,   648 }},
  {{  649,   644,   648 }},
  {{  649,   648,   164 }},
  {{  165,   649,   651 }},
  {{  649,   164,   650 }},
  {{  651,   649,   650 }},
  {{  651,   650,    45 }},
  {{   43,   652,   654 }},
  {{  652,   166,   653 }},
  {{  654,   652,   653 }},
  {{  654,   653,   168 }},
  {{  166,   655,   657 }},
  {{  655,    13,   656 }},
  {{  657,   655,   656 }},
  {{  657,   656,   167 }},
  {{  168,   653,   658 }},
  {{  653,   166,   657 }},
  {{  658,   653,   657 }},
  {{  658,   657,   167 }},
  {{  168,   658,   660 }},
  {{  658,   167,   659 }},
  {{  660,   658,   659 }},
  {{  660,   659,    44 }},
  {{   45,   650,   662 }},
  {{  650,   164,   661 }},
  {{  662,   650,   661 }},
  {{  662,   661,   169 }},
  {{  164,   647,   663 }},
  {{  647,    43,   654 }},
  {{  663,   647,   654 }},
  {{  663,   654,   168 }},
  {{  169,   661,   664 }},
  {{  661,   164,   663 }},
  {{  664,   661,   663 }},
  {{  664,   663,   168 }},
  {{  169,   664,   665 }},
  {{  664,   168,   660 }},
  {{  665,   664,   660 }},
  {{  665,   660,    44 }},
  {{   45,   662,   667 }},
  {{  662,   169,   666 }},
  {{  667,   662,   666 }},
  {{  667,   666,   171 }},
  {{  169,   665,   669 }},
  {{  665,    44,   668 }},
  {{  669,   665,   668 }},
  {{  669,   668,   170 }},
  {{  171,   666,   670 }},
  {{  666,   169,   669 }},
  {{  670,   666,   669 }},
  {{  670,   669,   170 }},
  {{  171,   670,   672 }},
  {{  670,   170,   671 }},
  {{  672,   670,   671 }},
  {{  672,   671,    15 }},
  {{   13,   673,   675 }},
  {{  673,   172,   674 }},
  {{  675,   673,   674 }},
  {{  675,   674,   174 }},
  {{  172,   676,   678 }},
  {{  676,    46,   677 }},
  {{  678,   676,   677 }},
  {{  678,   677,   173 }},
  {{  174,   674,   679 }},
  {{  674,   172,   678 }},
  {{  679,   674,   678 }},
  {{  679,   678,   173 }},
  {{  174,   679,   681 }},
  {{  679,   173,   680 }},
  {{  681,   679,   680 }},
  {{  681,   680,    48 }},
  {{   46,   682,   684 }},
  {{  682,   175,   683 }},
  {{  684,   682,   683 }},
  {{  684,   683,   177 }},
  {{  175,   685,   687 }},
  {{  685,     5,   686 }},
  {{  687,   685,   686 }},
  {{  687,   686,   176 }},
  {{  177,   683,   688 }},
  {{  683,   175,   687 }},
  {{  688,   683,   687 }},
  {{  688,   687,   176 }},
  {{  177,   688,   690 }},
  {{  688,   176,   689 }},
  {{  690,   688,   689 }},
  {{  690,   689,    47 }},
  {{   48,   680,   692 }},
  {{  680,   173,   691 }},
  {{  692,   680,   691 }},
  {{  692,   691,   178 }},
  {{  173,   677,   693 }},
  {{  677,    46,   684 }},
  {{  693,   677,   684 }},
  {{  693,   684,   177 }},
  {{  178,   691,   694 }},
  {{  691,   173,   693 }},
  {{  694,   691,   693 }},
  {{  694,   693,   177 }},
  {{  178,   694,   695 }},
  {{  694,   177,   690 }},
  {{  695,   694,   690 }},
  {{  695,   690,    47 }},
  {{   48,   692,   697 }},
  {{  692,   178,   696 }},
  {{  697,   692,   696 }},
  {{  697,   696,   180 }},
  {{  178,   695,   699 }},
  {{  695,    47,   698 }},
  {{  699,   695,   698 }},
  {{  699,   698,   179 }},
  {{  180,   696,   700 }},
  {{  696,   178,   699 }},
  {{  700,   696,   699 }},
  {{  700,   699,   179 }},
  {{  180,   700,   702 }},
  {{  700,   179,   701 }},
  {{  702,   700,   701 }},
  {{  702,   701,    14 }},
  {{   15,   671,   704 }},
  {{  671,   170,   703 }},
  {{  704,   671,   703 }},
  {{  704,   703,   182 }},
  {{  170,   668,   706 }},
  {{  668,    44,   705 }},
  {{  706,   668,   705 }},
  {{  706,   705,   181 }},
  {{  182,   703,   707 }},
  {{  703,   170,   706 }},
  {{  707,   703,   706 }},
  {{  707,   706,   181 }},
  {{  182,   707,   709 }},
  {{  707,   181,   708 }},
  {{  709,   707,   708 }},
  {{  709,   708,    49 }},
  {{   44,   659,   711 }},
  {{  659,   167,   710 }},
  {{  711,   659,   710 }},
  {{  711,   710,   183 }},
  {{  167,   656,   712 }},
  {{  656,    13,   675 }},
  {{  712,   656,   675 }},
  {{  712,   675,   174 }},
  {{  183,   710,   713 }},
  {{  710,   167,   712 }},
  {{  713,   710,   712 }},
  {{  713,   712,   174 }},
  {{  183,   713,   714 }},
  {{  713,   174,   681 }},
  {{  714,   713,   681 }},
  {{  714,   681,    48 }},
  {{   49,   708,   716 }},
  {{  708,   181,   715 }},
  {{  716,   708,   715 }},
  {{  716,   715,   184 }},
  {{  181,   705,   717 }},
  {{  705,    44,   711 }},
  {{  717,   705,   711 }},
  {{  717,   711,   183 }},
  {{  184,   715,   718 }},
  {{  715,   181,   717 }},
  {{  718,   715,   717 }},
  {{  718,   717,   183 }},
  {{  184,   718,   719 }},
  {{  718,   183,   714 }},
  {{  719,   718,   714 }},
  {{  719,   714,    48 }},
  {{   49,   716,   721 }},
  {{  716,   184,   720 }},
  {{  721,   716,   720 }},
  {{  721,   720,   185 }},
  {{  184,   719,   722 }},
  {{  719,    48,   697 }},
  {{  722,   719,   697 }},
  {{  722,   697,   180 }},
  {{  185,   720,   723 }},
  {{  720,   184,   722 }},
  {{  723,   720,   722 }},
  {{  723,   722,   180 }},
  {{  185,   723,   724 }},
  {{  723,   180,   702 }},
  {{  724,   723,   702 }},
  {{  724,   702,    14 }},
  {{   15,   704,   726 }},
  {{  704,   182,   725 }},
  {{  726,   704,   725 }},
  {{  726,   725,   187 }},
  {{  182,   709,   728 }},
  {{  709,    49,   727 }},
  {{  728,   709,   727 }},
  {{  728,   727,   186 }},
  {{  187,   725,   729 }},
  {{  725,   182,   728 }},
  {{  729,   725,   728 }},
  {{  729,   728,   186 }},
  {{  187,   729,   731 }},
  {{  729,   186,   730 }},
  {{  731,   729,   730 }},
  {{  731,   730,    51 }},
  {{   49,   721,   733 }},
  {{  721,   185,   732 }},
  {{  733,   721,   732 }},
  {{  733,   732,   189 }},
  {{  185,   724,   735 }},
  {{  724,    14,   734 }},
  {{  735,   724,   734 }},
  {{  735,   734,   188 }},
  {{  189,   732,   736 }},
  {{  732,   185,   735 }},
  {{  736,   732,   735 }},
  {{  736,   735,   188 }},
  {{  189,   736,   738 }},
  {{  736,   188,   737 }},
  {{  738,   736,   737 }},
  {{  738,   737,    50 }},
  {{   51,   730,   740 }},
  {{  730,   186,   739 }},
  {{  740,   730,   739 }},
  {{  740,   739,   190 }},
  {{  186,   727,   741 }},
  {{  727,    49,   733 }},
  {{  741,   727,   733 }},
  {{  741,   733,   189 }},
  {{  190,   739,   742 }},
  {{  739,   186,   741 }},
  {{  742,   739,   741 }},
  {{  742,   741,   189 }},
  {{  190,   742,   743 }},
  {{  742,   189,   738 }},
  {{  743,   742,   738 }},
  {{  743,   738,    50 }},
  {{   51,   740,   745 }},
  {{  740,   190,   744 }},
  {{  745,   740,   744 }},
  {{  745,   744,   192 }},
  {{  190,   743,   747 }},
  {{  743,    50,   746 }},
  {{  747,   743,   746 }},
  {{  747,   746,   191 }},
  {{  192,   744,   748 }},
  {{  744,   190,   747 }},
  {{  748,   744,   747 }},
  {{  748,   747,   191 }},
  {{  192,   748,   750 }},
  {{  748,   191,   749 }},
  {{  750,   748,   749 }},
  {{  750,   749,     4 }},
  {{    1,   751,   643 }},
  {{  751,   193,   752 }},
  {{  643,   751,   752 }},
  {{  643,   752,   163 }},
  {{  193,   753,   755 }},
  {{  753,    52,   754 }},
  {{  755,   753,   754 }},
  {{  755,   754,   194 }},
  {{  163,   752,   756 }},
  {{  752,   193,   755 }},
  {{  756,   752,   755 }},
  {{  756,   755,   194 }},
  {{  163,   756,   646 }},
  {{  756,   194,   757 }},
  {{  646,   756,   757 }},
  {{  646,   757,    43 }},
  {{   52,   758,   760 }},
  {{  758,   195,   759 }},
  {{  760,   758,   759 }},
  {{  760,   759,   197 }},
  {{  195,   761,   763 }},
  {{  761,    16,   762 }},
  {{  763,   761,   762 }},
  {{  763,   762,   196 }},
  {{  197,   759,   764 }},
  {{  759,   195,   763 }},
  {{  764,   759,   763 }},
  {{  764,   763,   196 }},
  {{  197,   764,   766 }},
  {{  764,   196,   765 }},
  {{  766,   764,   765 }},
  {{  766,   765,    53 }},
  {{   43,   757,   768 }},
  {{  757,   194,   767 }},
  {{  768,   757,   767 }},
  {{  768,   767,   198 }},
  {{  194,   754,   769 }},
  {{  754,    52,   760 }},
  {{  769,   754,   760 }},
  {{  769,   760,   197 }},
  {{  198,   767,   770 }},
  {{  767,   194,   769 }},
  {{  770,   767,   769 }},
  {{  770,   769,   197 }},
  {{  198,   770,   771 }},
  {{  770,   197,   766 }},
  {{  771,   770,   766 }},
  {{  771,   766,    53 }},
  {{   43,   768,   652 }},
  {{  768,   198,   772 }},
  {{  652,   768,   772 }},
  {{  652,   772,   166 }},
  {{  198,   771,   774 }},
  {{  771,    53,   773 }},
  {{  774,   771,   773 }},
  {{  774,   773,   199 }},
  {{  166,   772,   775 }},
  {{  772,   198,   774 }},
  {{  775,   772,   774 }},
  {{  775,   774,   199 }},
  {{  166,   775,   655 }},
  {{  775,   199,   776 }},
  {{  655,   775,   776 }},
  {{  655,   776,    13 }},
  {{   16,   777,   779 }},
  {{  777,   200,   778 }},
  {{  779,   777,   778 }},
  {{  779,   778,   202 }},
  {{  200,   780,   782 }},
  {{  780,    54,   781 }},
  {{  782,   780,   781 }},
  {{  782,   781,   201 }},
  {{  202,   778,   783 }},
  {{  778,   200,   782 }},
  {{  783,   778,   782 }},
  {{  783,   782,   201 }},
  {{  202,   783,   785 }},
  {{  783,   201,   784 }},
  {{  785,   783,   784 }},
  {{  785,   784,    56 }},
  {{   54,   786,   788 }},
  {{  786,   203,   787 }},
  {{  788,   786,   787 }},
  {{  788,   787,   205 }},
  {{  203,   789,   791 }},
  {{  789,     6,   790 }},
  {{  791,   789,   790 }},
  {{  791,   790,   204 }},
  {{  205,   787,   792 }},
  {{  787,   203,   791 }},
  {{  792,   787,   791 }},
  {{  792,   791,   204 }},
  {{  205,   792,   794 }},
  {{  792,   204,   793 }},
  {{  794,   792,   793 }},
  {{  794,   793,    55 }},
  {{   56,   784,   796 }},
  {{  784,   201,   795 }},
  {{  796,   784,   795 }},
  {{  796,   795,   206 }},
  {{  201,   781,   797 }},
  {{  781,    54,   788 }},
  {{  797,   781,   788 }},
  {{  797,   788,   205 }},
  {{  206,   795,   798 }},
  {{  795,   201,   797 }},
  {{  798,   795,   797 }},
  {{  798,   797,   205 }},
  {{  206,   798,   799 }},
  {{  798,   205,   794 }},
  {{  799,   798,   794 }},
  {{  799,   794,    55 }},
  {{   56,   796,   801 }},
  {{  796,   206,   800 }},
  {{  801,   796,   800 }},
  {{  801,   800,   208 }},
  {{  206,   799,   803 }},
  {{  799,    55,   802 }},
  {{  803,   799,   802 }},
  {{  803,   802,   207 }},
  {{  208,   800,   804 }},
  {{  800,   206,   803 }},
  {{  804,   800,   803 }},
  {{  804,   803,   207 }},
  {{  208,   804,   806 }},
  {{  804,   207,   805 }},
  {{  806,   804,   805 }},
  {{  806,   805,    17 }},
  {{   13,   776,   808 }},
  {{  776,   199,   807 }},
  {{  808,   776,   807 }},
  {{  808,   807,   210 }},
  {{  199,   773,   810 }},
  {{  773,    53,   809 }},
  {{  810,   773,   809 }},
  {{  810,   809,   209 }},
  {{  210,   807,   811 }},
  {{  807,   199,   810 }},
  {{  811,   807,   810 }},
  {{  811,   810,   209 }},
  {{  210,   811,   813 }},
  {{  811,   209,   812 }},
  {{  813,   811,   812 }},
  {{  813,   812,    57 }},
  {{   53,   765,   815 }},
  {{  765,   196,   814 }},
  {{  815,   765,   814 }},
  {{  815,   814,   211 }},
  {{  196,   762,   816 }},
  {{  762,    16,   779 }},
  {{  816,   762,   779 }},
  {{  816,   779,   202 }},
  {{  211,   814,   817 }},
  {{  814,   196,   816 }},
  {{  817,   814,   816 }},
  {{  817,   816,   202 }},
  {{  211,   817,   818 }},
  {{  817,   202,   785 }},
  {{  818,   817,   785 }},
  {{  818,   785,    56 }},
  {{   57,   812,   820 }},
  {{  812,   209,   819 }},
  {{  820,   812,   819 }},
  {{  820,   819,   212 }},
  {{  209,   809,   821 }},
  {{  809,    53,   815 }},
  {{  821,   809,   815 }},
  {{  821,   815,   211 }},
  {{  212,   819,   822 }},
  {{  819,   209,   821 }},
  {{  822,   819,   821 }},
  {{  822,   821,   211 }},
  {{  212,   822,   823 }},
  {{  822,   211,   818 }},
  {{  823,   822,   818 }},
  {{  823,   818,    56 }},
  {{   57,   820,   825 }},
  {{  820,   212,   824 }},
  {{  825,   820,   824 }},
  {{  825,   824,   213 }},
  {{  212,   823,   826 }},
  {{  823,    56,   801 }},
  {{  826,   823,   801 }},
  {{  826,   801,   208 }},
  {{  213,   824,   827 }},
  {{  824,   212,   826 }},
  {{  827,   824,   826 }},
  {{  827,   826,   208 }},
  {{  213,   827,   828 }},
  {{  827,   208,   806 }},
  {{  828,   827,   806 }},
  {{  828,   806,    17 }},
  {{   13,   808,   673 }},
  {{  808,   210,   829 }},
  {{  673,   808,   829 }},
  {{  673,   829,   172 }},
  {{  210,   813,   831 }},
  {{  813,    57,   830 }},
  {{  831,   813,   830 }},
  {{  831,   830,   214 }},
  {{  172,   829,   832 }},
  {{  829,   210,   831 }},
  {{  832,   829,   831 }},
  {{  832,   831,   214 }},
  {{  172,   832,   676 }},
  {{  832,   214,   833 }},
  {{  676,   832,   833 }},
  {{  676,   833,    46 }},
  {{   57,   825,   835 }},
  {{  825,   213,   834 }},
  {{  835,   825,   834 }},
  {{  835,   834,   216 }},
  {{  213,   828,   837 }},
  {{  828,    17,   836 }},
  {{  837,   828,   836 }},
  {{  837,   836,   215 }},
  {{  216,   834,   838 }},
  {{  834,   213,   837 }},
  {{  838,   834,   837 }},
  {{  838,   837,   215 }},
  {{  216,   838,   840 }},
  {{  838,   215,   839 }},
  {{  840,   838,   839 }},
  {{  840,   839,    58 }},
  {{   46,   833,   842 }},
  {{  833,   214,   841 }},
  {{  842,   833,   841 }},
  {{  842,   841,   217 }},
  {{  214,   830,   843 }},
  {{  830,    57,   835 }},
  {{  843,   830,   835 }},
  {{  843,   835,   216 }},
  {{  217,   841,   844 }},
  {{  841,   214,   843 }},
  {{  844,   841,   843 }},
  {{  844,   843,   216 }},
  {{  217,   844,   845 }},
  {{  844,   216,   840 }},
  {{  845,   844,   840 }},
  {{  845,   840,    58 }},
  {{   46,   842,   682 }},
  {{  842,   217,   846 }},
  {{  682,   842,   846 }},
  {{  682,   846,   175 }},
  {{  217,   845,   848 }},
  {{  845,    58,   847 }},
  {{  848,   845,   847 }},
  {{  848,   847,   218 }},
  {{  175,   846,   849 }},
  {{  846,   217,   848 }},
  {{  849,   846,   848 }},
  {{  849,   848,   218 }},
  {{  175,   849,   685 }},
  {{  849,   218,   850 }},
  {{  685,   849,   850 }},
  {{  685,   850,     5 }},
  {{    1,   851,   751 }},
  {{  851,   219,   852 }},
  {{  751,   851,   852 }},
  {{  751,   852,   193 }},
  {{  219,   853,   855 }},
  {{  853,    59,   854 }},
  {{  855,   853,   854 }},
  {{  855,   854,   220 }},
  {{  193,   852,   856 }},
  {{  852,   219,   855 }},
  {{  856,   852,   855 }},
  {{  856,   855,   220 }},
  {{  193,   856,   753 }},
  {{  856,   220,   857 }},
  {{  753,   856,   857 }},
  {{  753,   857,    52 }},
  {{   59,   858,   860 }},
  {{  858,   221,   859 }},
  {{  860,   858,   859 }},
  {{  860,   859,   223 }},
  {{  221,   861,   863 }},
  {{  861,    18,   862 }},
  {{  863,   861,   862 }},
  {{  863,   862,   222 }},
  {{  223,   859,   864 }},
  {{  859,   221,   863 }},
  {{  864,   859,   863 }},
  {{  864,   863,   222 }},
  {{  223,   864,   866 }},
  {{  864,   222,   865 }},
  {{  866,   864,   865 }},
  {{  866,   865,    60 }},
  {{   52,   857,   868 }},
  {{  857,   220,   867 }},
  {{  868,   857,   867 }},
  {{  868,   867,   224 }},
  {{  220,   854,   869 }},
  {{  854,    59,   860 }},
  {{  869,   854,   860 }},
  {{  869,   860,   223 }},
  {{  224,   867,   870 }},
  {{  867,   220,   869 }},
  {{  870,   867,   869 }},
  {{  870,   869,   223 }},
  {{  224,   870,   871 }},
  {{  870,   223,   866 }},
  {{  871,   870,   866 }},
  {{  871,   866,    60 }},
  {{   52,   868,   758 }},
  {{  868,   224,   872 }},
  {{  758,   868,   872 }},
  {{  758,   872,   195 }},
  {{  224,   871,   874 }},
  {{  871,    60,   873 }},
  {{  874,   871,   873 }},
  {{  874,   873,   225 }},
  {{  195,   872,   875 }},
  {{  872,   224,   874 }},
  {{  875,   872,   874 }},
  {{  875,   874,   225 }},
  {{  195,   875,   761 }},
  {{  875,   225,   876 }},
  {{  761,   875,   876 }},
  {{  761,   876,    16 }},
  {{   18,   877,   879 }},
  {{  877,   226,   878 }},
  {{  879,   877,   878 }},
  {{  879,   878,   228 }},
  {{  226,   880,   882 }},
  {{  880,    61,   881 }},
  {{  882,   880,   881 }},
  {{  882,   881,   227 }},
  {{  228,   878,   883 }},
  {{  878,   226,   882 }},
  {{  883,   878,   882 }},
  {{  883,   882,   227 }},
  {{  228,   883,   885 }},
  {{  883,   227,   884 }},
  {{  885,   883,   884 }},
  {{  885,   884,    63 }},
  {{   61,   886,   888 }},
  {{  886,   229,   887 }},
  {{  888,   886,   887 }},
  {{  888,   887,   231 }},
  {{  229,   889,   891 }},
  {{  889,     2,   890 }},
  {{  891,   889,   890 }},
  {{  891,   890,   230 }},
  {{  231,   887,   892 }},
  {{  887,   229,   891 }},
  {{  892,   887,   891 }},
  {{  892,   891,   230 }},
  {{  231,   892,   894 }},
  {{  892,   230,   893 }},
  {{  894,   892,   893 }},
  {{  894,   893,    62 }},
  {{   63,   884,   896 }},
  {{  884,   227,   895 }},
  {{  896,   884,   895 }},
  {{  896,   895,   232 }},
  {{  227,   881,   897 }},
  {{  881,    61,   888 }},
  {{  897,   881,   888 }},
  {{  897,   888,   231 }},
  {{  232,   895,   898 }},
  {{  895,   227,   897 }},
  {{  898,   895,   897 }},
  {{  898,   897,   231 }},
  {{  232,   898,   899 }},
  {{  898,   231,   894 }},
  {{  899,   898,   894 }},
  {{  899,   894,    62 }},
  {{   63,   896,   901 }},
  {{  896,   232,   900 }},
  {{  901,   896,   900 }},
  {{  901,   900,   234 }},
  {{  232,   899,   903 }},
  {{  899,    62,   902 }},
  {{  903,   899,   902 }},
  {{  903,   902,   233 }},
  {{  234,   900,   904 }},
  {{  900,   232,   903 }},
  {{  904,   900,   903 }},
  {{  904,   903,   233 }},
  {{  234,   904,   906 }},
  {{  904,   233,   905 }},
  {{  906,   904,   905 }},
  {{  906,   905,    19 }},
  {{   16,   876,   908 }},
  {{  876,   225,   907 }},
  {{  908,   876,   907 }},
  {{  908,   907,   236 }},
  {{  225,   873,   910 }},
  {{  873,    60,   909 }},
  {{  910,   873,   909 }},
  {{  910,   909,   235 }},
  {{  236,   907,   911 }},
  {{  907,   225,   910 }},
  {{  911,   907,   910 }},
  {{  911,   910,   235 }},
  {{  236,   911,   913 }},
  {{  911,   235,   912 }},
  {{  913,   911,   912 }},
  {{  913,   912,    64 }},
  {{   60,   865,   915 }},
  {{  865,   222,   914 }},
  {{  915,   865,   914 }},
  {{  915,   914,   237 }},
  {{  222,   862,   916 }},
  {{  862,    18,   879 }},
  {{  916,   862,   879 }},
  {{  916,   879,   228 }},
  {{  237,   914,   917 }},
  {{  914,   222,   916 }},
  {{  917,   914,   916 }},
  {{  917,   916,   228 }},
  {{  237,   917,   918 }},
  {{  917,   228,   885 }},
  {{  918,   917,   885 }},
  {{  918,   885,    63 }},
  {{   64,   912,   920 }},
  {{  912,   235,   919 }},
  {{  920,   912,   919 }},
  {{  920,   919,   238 }},
  {{  235,   909,   921 }},
  {{  909,    60,   915 }},
  {{  921,   909,   915 }},
  {{  921,   915,   237 }},
  {{  238,   919,   922 }},
  {{  919,   235,   921 }},
  {{  922,   919,   921 }},
  {{  922,   921,   237 }},
  {{  238,   922,   923 }},
  {{  922,   237,   918 }},
  {{  923,   922,   918 }},
  {{  923,   918,    63 }},
  {{   64,   920,   925 }},
  {{  920,   238,   924 }},
  {{  925,   920,   924 }},
  {{  925,   924,   239 }},
  {{  238,   923,   926 }},
  {{  923,    63,   901 }},
  {{  926,   923,   901 }},
  {{  926,   901,   234 }},
  {{  239,   924,   927 }},
  {{  924,   238,   926 }},
  {{  927,   924,   926 }},
  {{  927,   926,   234 }},
  {{  239,   927,   928 }},
  {{  927,   234,   906 }},
  {{  928,   927,   906 }},
  {{  928,   906,    19 }},
  {{   16,   908,   777 }},
  {{  908,   236,   929 }},
  {{  777,   908,   929 }},
  {{  777,   929,   200 }},
  {{  236,   913,   931 }},
  {{  913,    64,   930 }},
  {{  931,   913,   930 }},
  {{  931,   930,   240 }},
  {{  200,   929,   932 }},
  {{  929,   236,   931 }},
  {{  932,   929,   931 }},
  {{  932,   931,   240 }},
  {{  200,   932,   780 }},
  {{  932,   240,   933 }},
  {{  780,   932,   933 }},
  {{  780,   933,    54 }},
  {{   64,   925,   935 }},
  {{  925,   239,   934 }},
  {{  935,   925,   934 }},
  {{  935,   934,   242 }},
  {{  239,   928,   937 }},
  {{  928,    19,   936 }},
  {{  937,   928,   936 }},
  {{  937,   936,   241 }},
  {{  242,   934,   938 }},
  {{  934,   239,   937 }},
  {{  938,   934,   937 }},
  {{  938,   937,   241 }},
  {{  242,   938,   940 }},
  {{  938,   241,   939 }},
  {{  940,   938,   939 }},
  {{  940,   939,    65 }},
  {{   54,   933,   942 }},
  {{  933,   240,   941 }},
  {{  942,   933,   941 }},
  {{  942,   941,   243 }},
  {{  240,   930,   943 }},
  {{  930,    64,   935 }},
  {{  943,   930,   935 }},
  {{  943,   935,   242 }},
  {{  243,   941,   944 }},
  {{  941,   240,   943 }},
  {{  944,   941,   943 }},
  {{  944,   943,   242 }},
  {{  243,   944,   945 }},
  {{  944,   242,   940 }},
  {{  945,   944,   940 }},
  {{  945,   940,    65 }},
  {{   54,   942,   786 }},
  {{  942,   243,   946 }},
  {{  786,   942,   946 }},
  {{  786,   946,   203 }},
  {{  243,   945,   948 }},
  {{  945,    65,   947 }},
  {{  948,   945,   947 }},
  {{  948,   947,   244 }},
  {{  203,   946,   949 }},
  {{  946,   243,   948 }},
  {{  949,   946,   948 }},
  {{  949,   948,   244 }},
  {{  203,   949,   789 }},
  {{  949,   244,   950 }},
  {{  789,   949,   950 }},
  {{  789,   950,     6 }},
  {{    1,   951,   851 }},
  {{  951,   245,   952 }},
  {{  851,   951,   952 }},
  {{  851,   952,   219 }},
  {{  245,   953,   955 }},
  {{  953,    66,   954 }},
  {{  955,   953,   954 }},
  {{  955,   954,   246 }},
  {{  219,   952,   956 }},
  {{  952,   245,   955 }},
  {{  956,   952,   955 }},
  {{  956,   955,   246 }},
  {{  219,   956,   853 }},
  {{  956,   246,   957 }},
  {{  853,   956,   957 }},
  {{  853,   957,    59 }},
  {{   66,   958,   960 }},
  {{  958,   247,   959 }},
  {{  960,   958,   959 }},
  {{  960,   959,   249 }},
  {{  247,   961,   963 }},
  {{  961,    20,   962 }},
  {{  963,   961,   962 }},
  {{  963,   962,   248 }},
  {{  249,   959,   964 }},
  {{  959,   247,   963 }},
  {{  964,   959,   963 }},
  {{  964,   963,   248 }},
  {{  249,   964,   966 }},
  {{  964,   248,   965 }},
  {{  966,   964,   965 }},
  {{  966,   965,    67 }},
  {{   59,   957,   968 }},
  {{  957,   246,   967 }},
  {{  968,   957,   967 }},
  {{  968,   967,   250 }},
  {{  246,   954,   969 }},
  {{  954,    66,   960 }},
  {{  969,   954,   960 }},
  {{  969,   960,   249 }},
  {{  250,   967,   970 }},
  {{  967,   246,   969 }},
  {{  970,   967,   969 }},
  {{  970,   969,   249 }},
  {{  250,   970,   971 }},
  {{  970,   249,   966 }},
  {{  971,   970,   966 }},
  {{  971,   966,    67 }},
  {{   59,   968,   858 }},
  {{  968,   250,   972 }},
  {{  858,   968,   972 }},
  {{  858,   972,   221 }},
  {{  250,   971,   974 }},
  {{  971,    67,   973 }},
  {{  974,   971,   973 }},
  {{  974,   973,   251 }},
  {{  221,   972,   975 }},
  {{  972,   250,   974 }},
  {{  975,   972,   974 }},
  {{  975,   974,   251 }},
  {{  221,   975,   861 }},
  {{  975,   251,   976 }},
  {{  861,   975,   976 }},
  {{  861,   976,    18 }},
  {{   20,   977,   979 }},
  {{  977,   252,   978 }},
  {{  979,   977,   978 }},
  {{  979,   978,   254 }},
  {{  252,   980,   982 }},
  {{  980,    68,   981 }},
  {{  982,   980,   981 }},
  {{  982,   981,   253 }},
  {{  254,   978,   983 }},
  {{  978,   252,   982 }},
  {{  983,   978,   982 }},
  {{  983,   982,   253 }},
  {{  254,   983,   985 }},
  {{  983,   253,   984 }},
  {{  985,   983,   984 }},
  {{  985,   984,    70 }},
  {{   68,   986,   988 }},
  {{  986,   255,   987 }},
  {{  988,   986,   987 }},
  {{  988,   987,   257 }},
  {{  255,   989,   991 }},
  {{  989,     3,   990 }},
  {{  991,   989,   990 }},
  {{  991,   990,   256 }},
  {{  257,   987,   992 }},
  {{  987,   255,   991 }},
  {{  992,   987,   991 }},
  {{  992,   991,   256 }},
  {{  257,   992,   994 }},
  {{  992,   256,   993 }},
  {{  994,   992,   993 }},
  {{  994,   993,    69 }},
  {{   70,   984,   996 }},
  {{  984,   253,   995 }},
  {{  996,   984,   995 }},
  {{  996,   995,   258 }},
  {{  253,   981,   997 }},
  {{  981,    68,   988 }},
  {{  997,   981,   988 }},
  {{  997,   988,   257 }},
  {{  258,   995,   998 }},
  {{  995,   253,   997 }},
  {{  998,   995,   997 }},
  {{  998,   997,   257 }},
  {{  258,   998,   999 }},
  {{  998,   257,   994 }},
  {{  999,   998,   994 }},
  {{  999,   994,    69 }},
  {{   70,   996,  1001 }},
  {{  996,   258,  1000 }},
  {{ 1001,   996,  1000 }},
  {{ 1001,  1000,   260 }},
  {{  258,   999,  1003 }},
  {{  999,    69,  1002 }},
  {{ 1003,   999,  1002 }},
  {{ 1003,  1002,   259 }},
  {{  260,  1000,  1004 }},
  {{ 1000,   258,  1003 }},
  {{ 1004,  1000,  1003 }},
  {{ 1004,  1003,   259 }},
  {{  260,  1004,  1006 }},
  {{ 1004,   259,  1005 }},
  {{ 1006,  1004,  1005 }},
  {{ 1006,  1005,    21 }},
  {{   18,   976,  1008 }},
  {{  976,   251,  1007 }},
  {{ 1008,   976,  1007 }},
  {{ 1008,  1007,   262 }},
  {{  251,   973,  1010 }},
  {{  973,    67,  1009 }},
  {{ 1010,   973,  1009 }},
  {{ 1010,  1009,   261 }},
  {{  262,  1007,  1011 }},
  {{ 1007,   251,  1010 }},
  {{ 1011,  1007,  1010 }},
  {{ 1011,  1010,   261 }},
  {{  262,  1011,  1013 }},
  {{ 1011,   261,  1012 }},
  {{ 1013,  1011,  1012 }},
  {{ 1013,  1012,    71 }},
  {{   67,   965,  1015 }},
  {{  965,   248,  1014 }},
  {{ 1015,   965,  1014 }},
  {{ 1015,  1014,   263 }},
  {{  248,   962,  1016 }},
  {{  962,    20,   979 }},
  {{ 1016,   962,   979 }},
  {{ 1016,   979,   254 }},
  {{  263,  1014,  1017 }},
  {{ 1014,   248,  1016 }},
  {{ 1017,  1014,  1016 }},
  {{ 1017,  1016,   254 }},
  {{  263,  1017,  1018 }},
  {{ 1017,   254,   985 }},
  {{ 1018,  1017,   985 }},
  {{ 1018,   985,    70 }},
  {{   71,  1012,  1020 }},
  {{ 1012,   261,  1019 }},
  {{ 1020,  1012,  1019 }},
  {{ 1020,  1019,   264 }},
  {{  261,  1009,  1021 }},
  {{ 1009,    67,  1015 }},
  {{ 1021,  1009,  1015 }},
  {{ 1021,  1015,   263 }},
  {{  264,  1019,  1022 }},
  {{ 1019,   261,  1021 }},
  {{ 1022,  1019,  1021 }},
  {{ 1022,  1021,   263 }},
  {{  264,  1022,  1023 }},
  {{ 1022,   263,  1018 }},
  {{ 1023,  1022,  1018 }},
  {{ 1023,  1018,    70 }},
  {{   71,  1020,  1025 }},
  {{ 1020,   264,  1024 }},
  {{ 1025,  1020,  1024 }},
  {{ 1025,  1024,   265 }},
  {{  264,  1023,  1026 }},
  {{ 1023,    70,  1001 }},
  {{ 1026,  1023,  1001 }},
  {{ 1026,  1001,   260 }},
  {{  265,  1024,  1027 }},
  {{ 1024,   264,  1026 }},
  {{ 1027,  1024,  1026 }},
  {{ 1027,  1026,   260 }},
  {{  265,  1027,  1028 }},
  {{ 1027,   260,  1006 }},
  {{ 1028,  1027,  1006 }},
  {{ 1028,  1006,    21 }},
  {{   18,  1008,   877 }},
  {{ 1008,   262,  1029 }},
  {{  877,  1008,  1029 }},
  {{  877,  1029,   226 }},
  {{  262,  1013,  1031 }},
  {{ 1013,    71,  1030 }},
  {{ 1031,  1013,  1030 }},
  {{ 1031,  1030,   266 }},
  {{  226,  1029,  1032 }},
  {{ 1029,   262,  1031 }},
  {{ 1032,  1029,  1031 }},
  {{ 1032,  1031,   266 }},
  {{  226,  1032,   880 }},
  {{ 1032,   266,  1033 }},
  {{  880,  1032,  1033 }},
  {{  880,  1033,    61 }},
  {{   71,  1025,  1035 }},
  {{ 1025,   265,  1034 }},
  {{ 1035,  1025,  1034 }},
  {{ 1035,  1034,   268 }},
  {{  265,  1028,  1037 }},
  {{ 1028,    21,  1036 }},
  {{ 1037,  1028,  1036 }},
  {{ 1037,  1036,   267 }},
  {{  268,  1034,  1038 }},
  {{ 1034,   265,  1037 }},
  {{ 1038,  1034,  1037 }},
  {{ 1038,  1037,   267 }},
  {{  268,  1038,  1040 }},
  {{ 1038,   267,  1039 }},
  {{ 1040,  1038,  1039 }},
  {{ 1040,  1039,    72 }},
  {{   61,  1033,  1042 }},
  {{ 1033,   266,  1041 }},
  {{ 1042,  1033,  1041 }},
  {{ 1042,  1041,   269 }},
  {{  266,  1030,  1043 }},
  {{ 1030,    71,  1035 }},
  {{ 1043,  1030,  1035 }},
  {{ 1043,  1035,   268 }},
  {{  269,  1041,  1044 }},
  {{ 1041,   266,  1043 }},
  {{ 1044,  1041,  1043 }},
  {{ 1044,  1043,   268 }},
  {{  269,  1044,  1045 }},
  {{ 1044,   268,  1040 }},
  {{ 1045,  1044,  1040 }},
  {{ 1045,  1040,    72 }},
  {{   61,  1042,   886 }},
  {{ 1042,   269,  1046 }},
  {{  886,  1042,  1046 }},
  {{  886,  1046,   229 }},
  {{  269,  1045,  1048 }},
  {{ 1045,    72,  1047 }},
  {{ 1048,  1045,  1047 }},
  {{ 1048,  1047,   270 }},
  {{  229,  1046,  1049 }},
  {{ 1046,   269,  1048 }},
  {{ 1049,  1046,  1048 }},
  {{ 1049,  1048,   270 }},
  {{  229,  1049,   889 }},
  {{ 1049,   270,  1050 }},
  {{  889,  1049,  1050 }},
  {{  889,  1050,     2 }},
  {{    1,   645,   951 }},
  {{  645,   165,  1051 }},
  {{  951,   645,  1051 }},
  {{  951,  1051,   245 }},
  {{  165,   651,  1053 }},
  {{  651,    45,  1052 }},
  {{ 1053,   651,  1052 }},
  {{ 1053,  1052,   271 }},
  {{  245,  1051,  1054 }},
  {{ 1051,   165,  1053 }},
  {{ 1054,  1051,  1053 }},
  {{ 1054,  1053,   271 }},
  {{  245,  1054,   953 }},
  {{ 1054,   271,  1055 }},
  {{  953,  1054,  1055 }},
  {{  953,  1055,    66 }},
  {{   45,   667,  1057 }},
  {{  667,   171,  1056 }},
  {{ 1057,   667,  1056 }},
  {{ 1057,  1056,   273 }},
  {{  171,   672,  1059 }},
  {{  672,    15,  1058 }},
  {{ 1059,   672,  1058 }},
  {{ 1059,  1058,   272 }},
  {{  273,  1056,  1060 }},
  {{ 1056,   171,  1059 }},
  {{ 1060,  1056,  1059 }},
  {{ 1060,  1059,   272 }},
  {{  273,  1060,  1062 }},
  {{ 1060,   272,  1061 }},
  {{ 1062,  1060,  1061 }},
  {{ 1062,  1061,    73 }},
  {{   66,  1055,  1064 }},
  {{ 1055,   271,  1063 }},
  {{ 1064,  1055,  1063 }},
  {{ 1064,  1063,   274 }},
  {{  271,  1052,  1065 }},
  {{ 1052,    45,  1057 }},
  {{ 1065,  1052,  1057 }},
  {{ 1065,  1057,   273 }},
  {{  274,  1063,  1066 }},
  {{ 1063,   271,  1065 }},
  {{ 1066,  1063,  1065 }},
  {{ 1066,  1065,   273 }},
  {{  274,  1066,  1067 }},
  {{ 1066,   273,  1062 }},
  {{ 1067,  1066,  1062 }},
  {{ 1067,  1062,    73 }},
  {{   66,  1064,   958 }},
  {{ 1064,   274,  1068 }},
  {{  958,  1064,  1068 }},
  {{  958,  1068,   247 }},
  {{  274,  1067,  1070 }},
  {{ 1067,    73,  1069 }},
  {{ 1070,  1067,  1069 }},
  {{ 1070,  1069,   275 }},
  {{  247,  1068,  1071 }},
  {{ 1068,   274,  1070 }},
  {{ 1071,  1068,  1070 }},
  {{ 1071,  1070,   275 }},
  {{  247,  1071,   961 }},
  {{ 1071,   275,  1072 }},
  {{  961,  1071,  1072 }},
  {{  961,  1072,    20 }},
  {{   15,   726,  1074 }},
  {{  726,   187,  1073 }},
  {{ 1074,   726,  1073 }},
  {{ 1074,  1073,   277 }},
  {{  187,   731,  1076 }},
  {{  731,    51,  1075 }},
  {{ 1076,   731,  1075 }},
  {{ 1076,  1075,   276 }},
  {{  277,  1073,  1077 }},
  {{ 1073,   187,  1076 }},
  {{ 1077,  1073,  1076 }},
  {{ 1077,  1076,   276 }},
  {{  277,  1077,  1079 }},
  {{ 1077,   276,  1078 }},
  {{ 1079,  1077,  1078 }},
  {{ 1079,  1078,    75 }},
  {{   51,   745,  1081 }},
  {{  745,   192,  1080 }},
  {{ 1081,   745,  1080 }},
  {{ 1081,  1080,   279 }},
  {{  192,   750,  1083 }},
  {{  750,     4,  1082 }},
  {{ 1083,   750,  1082 }},
  {{ 1083,  1082,   278 }},
  {{  279,  1080,  1084 }},
  {{ 1080,   192,  1083 }},
  {{ 1084,  1080,  1083 }},
  {{ 1084,  1083,   278 }},
  {{  279,  1084,  1086 }},
  {{ 1084,   278,  1085 }},
  {{ 1086,  1084,  1085 }},
  {{ 1086,  1085,    74 }},
  {{   75,  1078,  1088 }},
  {{ 1078,   276,  1087 }},
  {{ 1088,  1078,  1087 }},
  {{ 1088,  1087,   280 }},
  {{  276,  1075,  1089 }},
  {{ 1075,    51,  1081 }},
  {{ 1089,  1075,  1081 }},
  {{ 1089,  1081,   279 }},
  {{  280,  1087,  1090 }},
  {{ 1087,   276,  1089 }},
  {{ 1090,  1087,  1089 }},
  {{ 1090,  1089,   279 }},
  {{  280,  1090,  1091 }},
  {{ 1090,   279,  1086 }},
  {{ 1091,  1090,  1086 }},
  {{ 1091,  1086,    74 }},
  {{   75,  1088,  1093 }},
  {{ 1088,   280,  1092 }},
  {{ 1093,  1088,  1092 }},
  {{ 1093,  1092,   282 }},
  {{  280,  1091,  1095 }},
  {{ 1091,    74,  1094 }},
  {{ 1095,  1091,  1094 }},
  {{ 1095,  1094,   281 }},
  {{  282,  1092,  1096 }},
  {{ 1092,   280,  1095 }},
  {{ 1096,  1092,  1095 }},
  {{ 1096,  1095,   281 }},
  {{  282,  1096,  1098 }},
  {{ 1096,   281,  1097 }},
  {{ 1098,  1096,  1097 }},
  {{ 1098,  1097,    22 }},
  {{   20,  1072,  1100 }},
  {{ 1072,   275,  1099 }},
  {{ 1100,  1072,  1099 }},
  {{ 1100,  1099,   284 }},
  {{  275,  1069,  1102 }},
  {{ 1069,    73,  1101 }},
  {{ 1102,  1069,  1101 }},
  {{ 1102,  1101,   283 }},
  {{  284,  1099,  1103 }},
  {{ 1099,   275,  1102 }},
  {{ 1103,  1099,  1102 }},
  {{ 1103,  1102,   283 }},
  {{  284,  1103,  1105 }},
  {{ 1103,   283,  1104 }},
  {{ 1105,  1103,  1104 }},
  {{ 1105,  1104,    76 }},
  {{   73,  1061,  1107 }},
  {{ 1061,   272,  1106 }},
  {{ 1107,  1061,  1106 }},
  {{ 1107,  1106,   285 }},
  {{  272,  1058,  1108 }},
  {{ 1058,    15,  1074 }},
  {{ 1108,  1058,  1074 }},
  {{ 1108,  1074,   277 }},
  {{  285,  1106,  1109 }},
  {{ 1106,   272,  1108 }},
  {{ 1109,  1106,  1108 }},
  {{ 1109,  1108,   277 }},
  {{  285,  1109,  1110 }},
  {{ 1109,   277,  1079 }},
  {{ 1110,  1109,  1079 }},
  {{ 1110,  1079,    75 }},
  {{   76,  1104,  1112 }},
  {{ 1104,   283,  1111 }},
  {{ 1112,  1104,  1111 }},
  {{ 1112,  1111,   286 }},
  {{  283,  1101,  1113 }},
  {{ 1101,    73,  1107 }},
  {{ 1113,  1101,  1107 }},
  {{ 1113,  1107,   285 }},
  {{  286,  1111,  1114 }},
  {{ 1111,   283,  1113 }},
  {{ 1114,  1111,  1113 }},
  {{ 1114,  1113,   285 }},
  {{  286,  1114,  1115 }},
  {{ 1114,   285,  1110 }},
  {{ 1115,  1114,  1110 }},
  {{ 1115,  1110,    75 }},
  {{   76,  1112,  1117 }},
  {{ 1112,   286,  1116 }},
  {{ 1117,  1112,  1116 }},
  {{ 1117,  1116,   287 }},
  {{  286,  1115,  1118 }},
  {{ 1115,    75,  1093 }},
  {{ 1118,  1115,  1093 }},
  {{ 1118,  1093,   282 }},
  {{  287,  1116,  1119 }},
  {{ 1116,   286,  1118 }},
  {{ 1119,  1116,  1118 }},
  {{ 1119,  1118,   282 }},
  {{  287,  1119,  1120 }},
  {{ 1119,   282,  1098 }},
  {{ 1120,  1119,  1098 }},
  {{ 1120,  1098,    22 }},
  {{   20,  1100,   977 }},
  {{ 1100,   284,  1121 }},
  {{  977,  1100,  1121 }},
  {{  977,  1121,   252 }},
  {{  284,  1105,  1123 }},
  {{ 1105,    76,  1122 }},
  {{ 1123,  1105,  1122 }},
  {{ 1123,  1122,   288 }},
  {{  252,  1121,  1124 }},
  {{ 1121,   284,  1123 }},
  {{ 1124,  1121,  1123 }},
  {{ 1124,  1123,   288 }},
  {{  252,  1124,   980 }},
  {{ 1124,   288,  1125 }},
  {{  980,  1124,  1125 }},
  {{  980,  1125,    68 }},
  {{   76,  1117,  1127 }},
  {{ 1117,   287,  1126 }},
  {{ 1127,  1117,  1126 }},
  {{ 1127,  1126,   290 }},
  {{  287,  1120,  1129 }},
  {{ 1120,    22,  1128 }},
  {{ 1129,  1120,  1128 }},
  {{ 1129,  1128,   289 }},
  {{  290,  1126,  1130 }},
  {{ 1126,   287,  1129 }},
  {{ 1130,  1126,  1129 }},
  {{ 1130,  1129,   289 }},
  {{  290,  1130,  1132 }},
  {{ 1130,   289,  1131 }},
  {{ 1132,  1130,  1131 }},
  {{ 1132,  1131,    77 }},
  {{   68,  1125,  1134 }},
  {{ 1125,   288,  1133 }},
  {{ 1134,  1125,  1133 }},
  {{ 1134,  1133,   291 }},
  {{  288,  1122,  1135 }},
  {{ 1122,    76,  1127 }},
  {{ 1135,  1122,  1127 }},
  {{ 1135,  1127,   290 }},
  {{  291,  1133,  1136 }},
  {{ 1133,   288,  1135 }},
  {{ 1136,  1133,  1135 }},
  {{ 1136,  1135,   290 }},
  {{  291,  1136,  1137 }},
  {{ 1136,   290,  1132 }},
  {{ 1137,  1136,  1132 }},
  {{ 1137,  1132,    77 }},
  {{   68,  1134,   986 }},
  {{ 1134,   291,  1138 }},
  {{  986,  1134,  1138 }},
  {{  986,  1138,   255 }},
  {{  291,  1137,  1140 }},
  {{ 1137,    77,  1139 }},
  {{ 1140,  1137,  1139 }},
  {{ 1140,  1139,   292 }},
  {{  255,  1138,  1141 }},
  {{ 1138,   291,  1140 }},
  {{ 1141,  1138,  1140 }},
  {{ 1141,  1140,   292 }},
  {{  255,  1141,   989 }},
  {{ 1141,   292,  1142 }},
  {{  989,  1141,  1142 }},
  {{  989,  1142,     3 }},
  {{    4,  1143,  1082 }},
  {{ 1143,   293,  1144 }},
  {{ 1082,  1143,  1144 }},
  {{ 1082,  1144,   278 }},
  {{  293,  1145,  1147 }},
  {{ 1145,    78,  1146 }},
  {{ 1147,  1145,  1146 }},
  {{ 1147,  1146,   294 }},
  {{  278,  1144,  1148 }},
  {{ 1144,   293,  1147 }},
  {{ 1148,  1144,  1147 }},
  {{ 1148,  1147,   294 }},
  {{  278,  1148,  1085 }},
  {{ 1148,   294,  1149 }},
  {{ 1085,  1148,  1149 }},
  {{ 1085,  1149,    74 }},
  {{   78,  1150,  1152 }},
  {{ 1150,   295,  1151 }},
  {{ 1152,  1150,  1151 }},
  {{ 1152,  1151,   297 }},
  {{  295,  1153,  1155 }},
  {{ 1153,    23,  1154 }},
  {{ 1155,  1153,  1154 }},
  {{ 1155,  1154,   296 }},
  {{  297,  1151,  1156 }},
  {{ 1151,   295,  1155 }},
  {{ 1156,  1151,  1155 }},
  {{ 1156,  1155,   296 }},
  {{  297,  1156,  1158 }},
  {{ 1156,   296,  1157 }},
  {{ 1158,  1156,  1157 }},
  {{ 1158,  1157,    79 }},
  {{   74,  1149,  1160 }},
  {{ 1149,   294,  1159 }},
  {{ 1160,  1149,  1159 }},
  {{ 1160,  1159,   298 }},
  {{  294,  1146,  1161 }},
  {{ 1146,    78,  1152 }},
  {{ 1161,  1146,  1152 }},
  {{ 1161,  1152,   297 }},
  {{  298,  1159,  1162 }},
  {{ 1159,   294,  1161 }},
  {{ 1162,  1159,  1161 }},
  {{ 1162,  1161,   297 }},
  {{  298,  1162,  1163 }},
  {{ 1162,   297,  1158 }},
  {{ 1163,  1162,  1158 }},
  {{ 1163,  1158,    79 }},
  {{   74,  1160,  1094 }},
  {{ 1160,   298,  1164 }},
  {{ 1094,  1160,  1164 }},
  {{ 1094,  1164,   281 }},
  {{  298,  1163,  1166 }},
  {{ 1163,    79,  1165 }},
  {{ 1166,  1163,  1165 }},
  {{ 1166,  1165,   299 }},
  {{  281,  1164,  1167 }},
  {{ 1164,   298,  1166 }},
  {{ 1167,  1164,  1166 }},
  {{ 1167,  1166,   299 }},
  {{  281,  1167,  1097 }},
  {{ 1167,   299,  1168 }},
  {{ 1097,  1167,  1168 }},
  {{ 1097,  1168,    22 }},
  {{   23,  1169,  1171 }},
  {{ 1169,   300,  1170 }},
  {{ 1171,  1169,  1170 }},
  {{ 1171,  1170,   302 }},
  {{  300,  1172,  1174 }},
  {{ 1172,    80,  1173 }},
  {{ 1174,  1172,  1173 }},
  {{ 1174,  1173,   301 }},
  {{  302,  1170,  1175 }},
  {{ 1170,   300,  1174 }},
  {{ 1175,  1170,  1174 }},
  {{ 1175,  1174,   301 }},
  {{  302,  1175,  1177 }},
  {{ 1175,   301,  1176 }},
  {{ 1177,  1175,  1176 }},
  {{ 1177,  1176,    82 }},
  {{   80,  1178,  1180 }},
  {{ 1178,   303,  1179 }},
  {{ 1180,  1178,  1179 }},
  {{ 1180,  1179,   305 }},
  {{  303,  1181,  1183 }},
  {{ 1181,     9,  1182 }},
  {{ 1183,  1181,  1182 }},
  {{ 1183,  1182,   304 }},
  {{  305,  1179,  1184 }},
  {{ 1179,   303,  1183 }},
  {{ 1184,  1179,  1183 }},
  {{ 1184,  1183,   304 }},
  {{  305,  1184,  1186 }},
  {{ 1184,   304,  1185 }},
  {{ 1186,  1184,  1185 }},
  {{ 1186,  1185,    81 }},
  {{   82,  1176,  1188 }},
  {{ 1176,   301,  1187 }},
  {{ 1188,  1176,  1187 }},
  {{ 1188,  1187,   306 }},
  {{  301,  1173,  1189 }},
  {{ 1173,    80,  1180 }},
  {{ 1189,  1173,  1180 }},
  {{ 1189,  1180,   305 }},
  {{  306,  1187,  1190 }},
  {{ 1187,   301,  1189 }},
  {{ 1190,  1187,  1189 }},
  {{ 1190,  1189,   305 }},
  {{  306,  1190,  1191 }},
  {{ 1190,   305,  1186 }},
  {{ 1191,  1190,  1186 }},
  {{ 1191,  1186,    81 }},
  {{   82,  1188,  1193 }},
  {{ 1188,   306,  1192 }},
  {{ 1193,  1188,  1192 }},
  {{ 1193,  1192,   308 }},
  {{  306,  1191,  1195 }},
  {{ 1191,    81,  1194 }},
  {{ 1195,  1191,  1194 }},
  {{ 1195,  1194,   307 }},
  {{  308,  1192,  1196 }},
  {{ 1192,   306,  1195 }},
  {{ 1196,  1192,  1195 }},
  {{ 1196,  1195,   307 }},
  {{  308,  1196,  1198 }},
  {{ 1196,   307,  1197 }},
  {{ 1198,  1196,  1197 }},
  {{ 1198,  1197,    24 }},
  {{   22,  1168,  1200 }},
  {{ 1168,   299,  1199 }},
  {{ 1200,  1168,  1199 }},
  {{ 1200,  1199,   310 }},
  {{  299,  1165,  1202 }},
  {{ 1165,    79,  1201 }},
  {{ 1202,  1165,  1201 }},
  {{ 1202,  1201,   309 }},
  {{  310,  1199,  1203 }},
  {{ 1199,   299,  1202 }},
  {{ 1203,  1199,  1202 }},
  {{ 1203,  1202,   309 }},
  {{  310,  1203,  1205 }},
  {{ 1203,   309,  1204 }},
  {{ 1205,  1203,  1204 }},
  {{ 1205,  1204,    83 }},
  {{   79,  1157,  1207 }},
  {{ 1157,   296,  1206 }},
  {{ 1207,  1157,  1206 }},
  {{ 1207,  1206,   311 }},
  {{  296,  1154,  1208 }},
  {{ 1154,    23,  1171 }},
  {{ 1208,  1154,  1171 }},
  {{ 1208,  1171,   302 }},
  {{  311,  1206,  1209 }},
  {{ 1206,   296,  1208 }},
  {{ 1209,  1206,  1208 }},
  {{ 1209,  1208,   302 }},
  {{  311,  1209,  1210 }},
  {{ 1209,   302,  1177 }},
  {{ 1210,  1209,  1177 }},
  {{ 1210,  1177,    82 }},
  {{   83,  1204,  1212 }},
  {{ 1204,   309,  1211 }},
  {{ 1212,  1204,  1211 }},
  {{ 1212,  1211,   312 }},
  {{  309,  1201,  1213 }},
  {{ 1201,    79,  1207 }},
  {{ 1213,  1201,  1207 }},
  {{ 1213,  1207,   311 }},
  {{  312,  1211,  1214 }},
  {{ 1211,   309,  1213 }},
  {{ 1214,  1211,  1213 }},
  {{ 1214,  1213,   311 }},
  {{  312,  1214,  1215 }},
  {{ 1214,   311,  1210 }},
  {{ 1215,  1214,  1210 }},
  {{ 1215,  1210,    82 }},
  {{   83,  1212,  1217 }},
  {{ 1212,   312,  1216 }},
  {{ 1217,  1212,  1216 }},
  {{ 1217,  1216,   313 }},
  {{  312,  1215,  1218 }},
  {{ 1215,    82,  1193 }},
  {{ 1218,  1215,  1193 }},
  {{ 1218,  1193,   308 }},
  {{  313,  1216,  1219 }},
  {{ 1216,   312,  1218 }},
  {{ 1219,  1216,  1218 }},
  {{ 1219,  1218,   308 }},
  {{  313,  1219,  1220 }},
  {{ 1219,   308,  1198 }},
  {{ 1220,  1219,  1198 }},
  {{ 1220,  1198,    24 }},
  {{   22,  1200,  1128 }},
  {{ 1200,   310,  1221 }},
  {{ 1128,  1200,  1221 }},
  {{ 1128,  1221,   289 }},
  {{  310,  1205,  1223 }},
  {{ 1205,    83,  1222 }},
  {{ 1223,  1205,  1222 }},
  {{ 1223,  1222,   314 }},
  {{  289,  1221,  1224 }},
  {{ 1221,   310,  1223 }},
  {{ 1224,  1221,  1223 }},
  {{ 1224,  1223,   314 }},
  {{  289,  1224,  1131 }},
  {{ 1224,   314,  1225 }},
  {{ 1131,  1224,  1225 }},
  {{ 1131,  1225,    77 }},
  {{   83,  1217,  1227 }},
  {{ 1217,   313,  1226 }},
  {{ 1227,  1217,  1226 }},
  {{ 1227,  1226,   316 }},
  {{  313,  1220,  1229 }},
  {{ 1220,    24,  1228 }},
  {{ 1229,  1220,  1228 }},
  {{ 1229,  1228,   315 }},
  {{  316,  1226,  1230 }},
  {{ 1226,   313,  1229 }},
  {{ 1230,  1226,  1229 }},
  {{ 1230,  1229,   315 }},
  {{  316,  1230,  1232 }},
  {{ 1230,   315,  1231 }},
  {{ 1232,  1230,  1231 }},
  {{ 1232,  1231,    84 }},
  {{   77,  1225,  1234 }},
  {{ 1225,   314,  1233 }},
  {{ 1234,  1225,  1233 }},
  {{ 1234,  1233,   317 }},
  {{  314,  1222,  1235 }},
  {{ 1222,    83,  1227 }},
  {{ 1235,  1222,  1227 }},
  {{ 1235,  1227,   316 }},
  {{  317,  1233,  1236 }},
  {{ 1233,   314,  1235 }},
  {{ 1236,  1233,  1235 }},
  {{ 1236,  1235,   316 }},
  {{  317,  1236,  1237 }},
  {{ 1236,   316,  1232 }},
  {{ 1237,  1236,  1232 }},
  {{ 1237,  1232,    84 }},
  {{   77,  1234,  1139 }},
  {{ 1234,   317,  1238 }},
  {{ 1139,  1234,  1238 }},
  {{ 1139,  1238,   292 }},
  {{  317,  1237,  1240 }},
  {{ 1237,    84,  1239 }},
  {{ 1240,  1237,  1239 }},
  {{ 1240,  1239,   318 }},
  {{  292,  1238,  1241 }},
  {{ 1238,   317,  1240 }},
  {{ 1241,  1238,  1240 }},
  {{ 1241,  1240,   318 }},
  {{  292,  1241,  1142 }},
  {{ 1241,   318,  1242 }},
  {{ 1142,  1241,  1242 }},
  {{ 1142,  1242,     3 }},
  {{    4,  1243,  1143 }},
  {{ 1243,   319,  1244 }},
  {{ 1143,  1243,  1244 }},
  {{ 1143,  1244,   293 }},
  {{  319,  1245,  1247 }},
  {{ 1245,    85,  1246 }},
  {{ 1247,  1245,  1246 }},
  {{ 1247,  1246,   320 }},
  {{  293,  1244,  1248 }},
  {{ 1244,   319,  1247 }},
  {{ 1248,  1244,  1247 }},
  {{ 1248,  1247,   320 }},
  {{  293,  1248,  1145 }},
  {{ 1248,   320,  1249 }},
  {{ 1145,  1248,  1249 }},
  {{ 1145,  1249,    78 }},
  {{   85,  1250,  1252 }},
  {{ 1250,   321,  1251 }},
  {{ 1252,  1250,  1251 }},
  {{ 1252,  1251,   323 }},
  {{  321,  1253,  1255 }},
  {{ 1253,    25,  1254 }},
  {{ 1255,  1253,  1254 }},
  {{ 1255,  1254,   322 }},
  {{  323,  1251,  1256 }},
  {{ 1251,   321,  1255 }},
  {{ 1256,  1251,  1255 }},
  {{ 1256,  1255,   322 }},
  {{  323,  1256,  1258 }},
  {{ 1256,   322,  1257 }},
  {{ 1258,  1256,  1257 }},
  {{ 1258,  1257,    86 }},
  {{   78,  1249,  1260 }},
  {{ 1249,   320,  1259 }},
  {{ 1260,  1249,  1259 }},
  {{ 1260,  1259,   324 }},
  {{  320,  1246,  1261 }},
  {{ 1246,    85,  1252 }},
  {{ 1261,  1246,  1252 }},
  {{ 1261,  1252,   323 }},
  {{  324,  1259,  1262 }},
  {{ 1259,   320,  1261 }},
  {{ 1262,  1259,  1261 }},
  {{ 1262,  1261,   323 }},
  {{  324,  1262,  1263 }},
  {{ 1262,   323,  1258 }},
  {{ 1263,  1262,  1258 }},
  {{ 1263,  1258,    86 }},
  {{   78,  1260,  1150 }},
  {{ 1260,   324,  1264 }},
  {{ 1150,  1260,  1264 }},
  {{ 1150,  1264,   295 }},
  {{  324,  1263,  1266 }},
  {{ 1263,    86,  1265 }},
  {{ 1266,  1263,  1265 }},
  {{ 1266,  1265,   325 }},
  {{  295,  1264,  1267 }},
  {{ 1264,   324,  1266 }},
  {{ 1267,  1264,  1266 }},
  {{ 1267,  1266,   325 }},
  {{  295,  1267,  1153 }},
  {{ 1267,   325,  1268 }},
  {{ 1153,  1267,  1268 }},
  {{ 1153,  1268,    23 }},
  {{   25,  1269,  1271 }},
  {{ 1269,   326,  1270 }},
  {{ 1271,  1269,  1270 }},
  {{ 1271,  1270,   328 }},
  {{  326,  1272,  1274 }},
  {{ 1272,    87,  1273 }},
  {{ 1274,  1272,  1273 }},
  {{ 1274,  1273,   327 }},
  {{  328,  1270,  1275 }},
  {{ 1270,   326,  1274 }},
  {{ 1275,  1270,  1274 }},
  {{ 1275,  1274,   327 }},
  {{  328,  1275,  1277 }},
  {{ 1275,   327,  1276 }},
  {{ 1277,  1275,  1276 }},
  {{ 1277,  1276,    89 }},
  {{   87,  1278,  1280 }},
  {{ 1278,   329,  1279 }},
  {{ 1280,  1278,  1279 }},
  {{ 1280,  1279,   331 }},
  {{  329,  1281,  1283 }},
  {{ 1281,    10,  1282 }},
  {{ 1283,  1281,  1282 }},
  {{ 1283,  1282,   330 }},
  {{  331,  1279,  1284 }},
  {{ 1279,   329,  1283 }},
  {{ 1284,  1279,  1283 }},
  {{ 1284,  1283,   330 }},
  {{  331,  1284,  1286 }},
  {{ 1284,   330,  1285 }},
  {{ 1286,  1284,  1285 }},
  {{ 1286,  1285,    88 }},
  {{   89,  1276,  1288 }},
  {{ 1276,   327,  1287 }},
  {{ 1288,  1276,  1287 }},
  {{ 1288,  1287,   332 }},
  {{  327,  1273,  1289 }},
  {{ 1273,    87,  1280 }},
  {{ 1289,  1273,  1280 }},
  {{ 1289,  1280,   331 }},
  {{  332,  1287,  1290 }},
  {{ 1287,   327,  1289 }},
  {{ 1290,  1287,  1289 }},
  {{ 1290,  1289,   331 }},
  {{  332,  1290,  1291 }},
  {{ 1290,   331,  1286 }},
  {{ 1291,  1290,  1286 }},
  {{ 1291,  1286,    88 }},
  {{   89,  1288,  1293 }},
  {{ 1288,   332,  1292 }},
  {{ 1293,  1288,  1292 }},
  {{ 1293,  1292,   334 }},
  {{  332,  1291,  1295 }},
  {{ 1291,    88,  1294 }},
  {{ 1295,  1291,  1294 }},
  {{ 1295,  1294,   333 }},
  {{  334,  1292,  1296 }},
  {{ 1292,   332,  1295 }},
  {{ 1296,  1292,  1295 }},
  {{ 1296,  1295,   333 }},
  {{  334,  1296,  1298 }},
  {{ 1296,   333,  1297 }},
  {{ 1298,  1296,  1297 }},
  {{ 1298,  1297,    26 }},
  {{   23,  1268,  1300 }},
  {{ 1268,   325,  1299 }},
  {{ 1300,  1268,  1299 }},
  {{ 1300,  1299,   336 }},
  {{  325,  1265,  1302 }},
  {{ 1265,    86,  1301 }},
  {{ 1302,  1265,  1301 }},
  {{ 1302,  1301,   335 }},
  {{  336,  1299,  1303 }},
  {{ 1299,   325,  1302 }},
  {{ 1303,  1299,  1302 }},
  {{ 1303,  1302,   335 }},
  {{  336,  1303,  1305 }},
  {{ 1303,   335,  1304 }},
  {{ 1305,  1303,  1304 }},
  {{ 1305,  1304,    90 }},
  {{   86,  1257,  1307 }},
  {{ 1257,   322,  1306 }},
  {{ 1307,  1257,  1306 }},
  {{ 1307,  1306,   337 }},
  {{  322,  1254,  1308 }},
  {{ 1254,    25,  1271 }},
  {{ 1308,  1254,  1271 }},
  {{ 1308,  1271,   328 }},
  {{  337,  1306,  1309 }},
  {{ 1306,   322,  1308 }},
  {{ 1309,  1306,  1308 }},
  {{ 1309,  1308,   328 }},
  {{  337,  1309,  1310 }},
  {{ 1309,   328,  1277 }},
  {{ 1310,  1309,  1277 }},
  {{ 1310,  1277,    89 }},
  {{   90,  1304,  1312 }},
  {{ 1304,   335,  1311 }},
  {{ 1312,  1304,  1311 }},
  {{ 1312,  1311,   338 }},
  {{  335,  1301,  1313 }},
  {{ 1301,    86,  1307 }},
  {{ 1313,  1301,  1307 }},
  {{ 1313,  1307,   337 }},
  {{  338,  1311,  1314 }},
  {{ 1311,   335,  1313 }},
  {{ 1314,  1311,  1313 }},
  {{ 1314,  1313,   337 }},
  {{  338,  1314,  1315 }},
  {{ 1314,   337,  1310 }},
  {{ 1315,  1314,  1310 }},
  {{ 1315,  1310,    89 }},
  {{   90,  1312,  1317 }},
  {{ 1312,   338,  1316 }},
  {{ 1317,  1312,  1316 }},
  {{ 1317,  1316,   339 }},
  {{  338,  1315,  1318 }},
  {{ 1315,    89,  1293 }},
  {{ 1318,  1315,  1293 }},
  {{ 1318,  1293,   334 }},
  {{  339,  1316,  1319 }},
  {{ 1316,   338,  1318 }},
  {{ 1319,  1316,  1318 }},
  {{ 1319,  1318,   334 }},
  {{  339,  1319,  1320 }},
  {{ 1319,   334,  1298 }},
  {{ 1320,  1319,  1298 }},
  {{ 1320,  1298,    26 }},
  {{   23,  1300,  1169 }},
  {{ 1300,   336,  1321 }},
  {{ 1169,  1300,  1321 }},
  {{ 1169,  1321,   300 }},
  {{  336,  1305,  1323 }},
  {{ 1305,    90,  1322 }},
  {{ 1323,  1305,  1322 }},
  {{ 1323,  1322,   340 }},
  {{  300,  1321,  1324 }},
  {{ 1321,   336,  1323 }},
  {{ 1324,  1321,  1323 }},
  {{ 1324,  1323,   340 }},
  {{  300,  1324,  1172 }},
  {{ 1324,   340,  1325 }},
  {{ 1172,  1324,  1325 }},
  {{ 1172,  1325,    80 }},
  {{   90,  1317,  1327 }},
  {{ 1317,   339,  1326 }},
  {{ 1327,  1317,  1326 }},
  {{ 1327,  1326,   342 }},
  {{  339,  1320,  1329 }},
  {{ 1320,    26,  1328 }},
  {{ 1329,  1320,  1328 }},
  {{ 1329,  1328,   341 }},
  {{  342,  1326,  1330 }},
  {{ 1326,   339,  1329 }},
  {{ 1330,  1326,  1329 }},
  {{ 1330,  1329,   341 }},
  {{  342,  1330,  1332 }},
  {{ 1330,   341,  1331 }},
  {{ 1332,  1330,  1331 }},
  {{ 1332,  1331,    91 }},
  {{   80,  1325,  1334 }},
  {{ 1325,   340,  1333 }},
  {{ 1334,  1325,  1333 }},
  {{ 1334,  1333,   343 }},
  {{  340,  1322,  1335 }},
  {{ 1322,    90,  1327 }},
  {{ 1335,  1322,  1327 }},
  {{ 1335,  1327,   342 }},
  {{  343,  1333,  1336 }},
  {{ 1333,   340,  1335 }},
  {{ 1336,  1333,  1335 }},
  {{ 1336,  1335,   342 }},
  {{  343,  1336,  1337 }},
  {{ 1336,   342,  1332 }},
  {{ 1337,  1336,  1332 }},
  {{ 1337,  1332,    91 }},
  {{   80,  1334,  1178 }},
  {{ 1334,   343,  1338 }},
  {{ 1178,  1334,  1338 }},
  {{ 1178,  1338,   303 }},
  {{  343,  1337,  1340 }},
  {{ 1337,    91,  1339 }},
  {{ 1340,  1337,  1339 }},
  {{ 1340,  1339,   344 }},
  {{  303,  1338,  1341 }},
  {{ 1338,   343,  1340 }},
  {{ 1341,  1338,  1340 }},
  {{ 1341,  1340,   344 }},
  {{  303,  1341,  1181 }},
  {{ 1341,   344,  1342 }},
  {{ 1181,  1341,  1342 }},
  {{ 1181,  1342,     9 }},
  {{    4,   749,  1243 }},
  {{  749,   191,  1343 }},
  {{ 1243,   749,  1343 }},
  {{ 1243,  1343,   319 }},
  {{  191,   746,  1345 }},
  {{  746,    50,  1344 }},
  {{ 1345,   746,  1344 }},
  {{ 1345,  1344,   345 }},
  {{  319,  1343,  1346 }},
  {{ 1343,   191,  1345 }},
  {{ 1346,  1343,  1345 }},
  {{ 1346,  1345,   345 }},
  {{  319,  1346,  1245 }},
  {{ 1346,   345,  1347 }},
  {{ 1245,  1346,  1347 }},
  {{ 1245,  1347,    85 }},
  {{   50,   737,  1349 }},
  {{  737,   188,  1348 }},
  {{ 1349,   737,  1348 }},
  {{ 1349,  1348,   347 }},
  {{  188,   734,  1351 }},
  {{  734,    14,  1350 }},
  {{ 1351,   734,  1350 }},
  {{ 1351,  1350,   346 }},
  {{  347,  1348,  1352 }},
  {{ 1348,   188,  1351 }},
  {{ 1352,  1348,  1351 }},
  {{ 1352,  1351,   346 }},
  {{  347,  1352,  1354 }},
  {{ 1352,   346,  1353 }},
  {{ 1354,  1352,  1353 }},
  {{ 1354,  1353,    92 }},
  {{   85,  1347,  1356 }},
  {{ 1347,   345,  1355 }},
  {{ 1356,  1347,  1355 }},
  {{ 1356,  1355,   348 }},
  {{  345,  1344,  1357 }},
  {{ 1344,    50,  1349 }},
  {{ 1357,  1344,  1349 }},
  {{ 1357,  1349,   347 }},
  {{  348,  1355,  1358 }},
  {{ 1355,   345,  1357 }},
  {{ 1358,  1355,  1357 }},
  {{ 1358,  1357,   347 }},
  {{  348,  1358,  1359 }},
  {{ 1358,   347,  1354 }},
  {{ 1359,  1358,  1354 }},
  {{ 1359,  1354,    92 }},
  {{   85,  1356,  1250 }},
  {{ 1356,   348,  1360 }},
  {{ 1250,  1356,  1360 }},
  {{ 1250,  1360,   321 }},
  {{  348,  1359,  1362 }},
  {{ 1359,    92,  1361 }},
  {{ 1362,  1359,  1361 }},
  {{ 1362,  1361,   349 }},
  {{  321,  1360,  1363 }},
  {{ 1360,   348,  1362 }},
  {{ 1363,  1360,  1362 }},
  {{ 1363,  1362,   349 }},
  {{  321,  1363,  1253 }},
  {{ 1363,   349,  1364 }},
  {{ 1253,  1363,  1364 }},
  {{ 1253,  1364,    25 }},
  {{   14,   701,  1366 }},
  {{  701,   179,  1365 }},
  {{ 1366,   701,  1365 }},
  {{ 1366,  1365,   351 }},
  {{  179,   698,  1368 }},
  {{  698,    47,  1367 }},
  {{ 1368,   698,  1367 }},
  {{ 1368,  1367,   350 }},
  {{  351,  1365,  1369 }},
  {{ 1365,   179,  1368 }},
  {{ 1369,  1365,  1368 }},
  {{ 1369,  1368,   350 }},
  {{  351,  1369,  1371 }},
  {{ 1369,   350,  1370 }},
  {{ 1371,  1369,  1370 }},
  {{ 1371,  1370,    94 }},
  {{   47,   689,  1373 }},
  {{  689,   176,  1372 }},
  {{ 1373,   689,  1372 }},
  {{ 1373,  1372,   353 }},
  {{  176,   686,  1375 }},
  {{  686,     5,  1374 }},
  {{ 1375,   686,  1374 }},
  {{ 1375,  1374,   352 }},
  {{  353,  1372,  1376 }},
  {{ 1372,   176,  1375 }},
  {{ 1376,  1372,  1375 }},
  {{ 1376,  1375,   352 }},
  {{  353,  1376,  1378 }},
  {{ 1376,   352,  1377 }},
  {{ 1378,  1376,  1377 }},
  {{ 1378,  1377,    93 }},
  {{   94,  1370,  1380 }},
  {{ 1370,   350,  1379 }},
  {{ 1380,  1370,  1379 }},
  {{ 1380,  1379,   354 }},
  {{  350,  1367,  1381 }},
  {{ 1367,    47,  1373 }},
  {{ 1381,  1367,  1373 }},
  {{ 1381,  1373,   353 }},
  {{  354,  1379,  1382 }},
  {{ 1379,   350,  1381 }},
  {{ 1382,  1379,  1381 }},
  {{ 1382,  1381,   353 }},
  {{  354,  1382,  1383 }},
  {{ 1382,   353,  1378 }},
  {{ 1383,  1382,  1378 }},
  {{ 1383,  1378,    93 }},
  {{   94,  1380,  1385 }},
  {{ 1380,   354,  1384 }},
  {{ 1385,  1380,  1384 }},
  {{ 1385,  1384,   356 }},
  {{  354,  1383,  1387 }},
  {{ 1383,    93,  1386 }},
  {{ 1387,  1383,  1386 }},
  {{ 1387,  1386,   355 }},
  {{  356,  1384,  1388 }},
  {{ 1384,   354,  1387 }},
  {{ 1388,  1384,  1387 }},
  {{ 1388,  1387,   355 }},
  {{  356,  1388,  1390 }},
  {{ 1388,   355,  1389 }},
  {{ 1390,  1388,  1389 }},
  {{ 1390,  1389,    27 }},
  {{   25,  1364,  1392 }},
  {{ 1364,   349,  1391 }},
  {{ 1392,  1364,  1391 }},
  {{ 1392,  1391,   358 }},
  {{  349,  1361,  1394 }},
  {{ 1361,    92,  1393 }},
  {{ 1394,  1361,  1393 }},
  {{ 1394,  1393,   357 }},
  {{  358,  1391,  1395 }},
  {{ 1391,   349,  1394 }},
  {{ 1395,  1391,  1394 }},
  {{ 1395,  1394,   357 }},
  {{  358,  1395,  1397 }},
  {{ 1395,   357,  1396 }},
  {{ 1397,  1395,  1396 }},
  {{ 1397,  1396,    95 }},
  {{   92,  1353,  1399 }},
  {{ 1353,   346,  1398 }},
  {{ 1399,  1353,  1398 }},
  {{ 1399,  1398,   359 }},
  {{  346,  1350,  1400 }},
  {{ 1350,    14,  1366 }},
  {{ 1400,  1350,  1366 }},
  {{ 1400,  1366,   351 }},
  {{  359,  1398,  1401 }},
  {{ 1398,   346,  1400 }},
  {{ 1401,  1398,  1400 }},
  {{ 1401,  1400,   351 }},
  {{  359,  1401,  1402 }},
  {{ 1401,   351,  1371 }},
  {{ 1402,  1401,  1371 }},
  {{ 1402,  1371,    94 }},
  {{   95,  1396,  1404 }},
  {{ 1396,   357,  1403 }},
  {{ 1404,  1396,  1403 }},
  {{ 1404,  1403,   360 }},
  {{  357,  1393,  1405 }},
  {{ 1393,    92,  1399 }},
  {{ 1405,  1393,  1399 }},
  {{ 1405,  1399,   359 }},
  {{  360,  1403,  1406 }},
  {{ 1403,   357,  1405 }},
  {{ 1406,  1403,  1405 }},
  {{ 1406,  1405,   359 }},
  {{  360,  1406,  1407 }},
  {{ 1406,   359,  1402 }},
  {{ 1407,  1406,  1402 }},
  {{ 1407,  1402,    94 }},
  {{   95,  1404,  1409 }},
  {{ 1404,   360,  1408 }},
  {{ 1409,  1404,  1408 }},
  {{ 1409,  1408,   361 }},
  {{  360,  1407,  1410 }},
  {{ 1407,    94,  1385 }},
  {{ 1410,  1407,  1385 }},
  {{ 1410,  1385,   356 }},
  {{  361,  1408,  1411 }},
  {{ 1408,   360,  1410 }},
  {{ 1411,  1408,  1410 }},
  {{ 1411,  1410,   356 }},
  {{  361,  1411,  1412 }},
  {{ 1411,   356,  1390 }},
  {{ 1412,  1411,  1390 }},
  {{ 1412,  1390,    27 }},
  {{   25,  1392,  1269 }},
  {{ 1392,   358,  1413 }},
  {{ 1269,  1392,  1413 }},
  {{ 1269,  1413,   326 }},
  {{  358,  1397,  1415 }},
  {{ 1397,    95,  1414 }},
  {{ 1415,  1397,  1414 }},
  {{ 1415,  1414,   362 }},
  {{  326,  1413,  1416 }},
  {{ 1413,   358,  1415 }},
  {{ 1416,  1413,  1415 }},
  {{ 1416,  1415,   362 }},
  {{  326,  1416,  1272 }},
  {{ 1416,   362,  1417 }},
  {{ 1272,  1416,  1417 }},
  {{ 1272,  1417,    87 }},
  {{   95,  1409,  1419 }},
  {{ 1409,   361,  1418 }},
  {{ 1419,  1409,  1418 }},
  {{ 1419,  1418,   364 }},
  {{  361,  1412,  1421 }},
  {{ 1412,    27,  1420 }},
  {{ 1421,  1412,  1420 }},
  {{ 1421,  1420,   363 }},
  {{  364,  1418,  1422 }},
  {{ 1418,   361,  1421 }},
  {{ 1422,  1418,  1421 }},
  {{ 1422,  1421,   363 }},
  {{  364,  1422,  1424 }},
  {{ 1422,   363,  1423 }},
  {{ 1424,  1422,  1423 }},
  {{ 1424,  1423,    96 }},
  {{   87,  1417,  1426 }},
  {{ 1417,   362,  1425 }},
  {{ 1426,  1417,  1425 }},
  {{ 1426,  1425,   365 }},
  {{  362,  1414,  1427 }},
  {{ 1414,    95,  1419 }},
  {{ 1427,  1414,  1419 }},
  {{ 1427,  1419,   364 }},
  {{  365,  1425,  1428 }},
  {{ 1425,   362,  1427 }},
  {{ 1428,  1425,  1427 }},
  {{ 1428,  1427,   364 }},
  {{  365,  1428,  1429 }},
  {{ 1428,   364,  1424 }},
  {{ 1429,  1428,  1424 }},
  {{ 1429,  1424,    96 }},
  {{   87,  1426,  1278 }},
  {{ 1426,   365,  1430 }},
  {{ 1278,  1426,  1430 }},
  {{ 1278,  1430,   329 }},
  {{  365,  1429,  1432 }},
  {{ 1429,    96,  1431 }},
  {{ 1432,  1429,  1431 }},
  {{ 1432,  1431,   366 }},
  {{  329,  1430,  1433 }},
  {{ 1430,   365,  1432 }},
  {{ 1433,  1430,  1432 }},
  {{ 1433,  1432,   366 }},
  {{  329,  1433,  1281 }},
  {{ 1433,   366,  1434 }},
  {{ 1281,  1433,  1434 }},
  {{ 1281,  1434,    10 }},
  {{    5,  1435,  1374 }},
  {{ 1435,   367,  1436 }},
  {{ 1374,  1435,  1436 }},
  {{ 1374,  1436,   352 }},
  {{  367,  1437,  1439 }},
  {{ 1437,    97,  1438 }},
  {{ 1439,  1437,  1438 }},
  {{ 1439,  1438,   368 }},
  {{  352,  1436,  1440 }},
  {{ 1436,   367,  1439 }},
  {{ 1440,  1436,  1439 }},
  {{ 1440,  1439,   368 }},
  {{  352,  1440,  1377 }},
  {{ 1440,   368,  1441 }},
  {{ 1377,  1440,  1441 }},
  {{ 1377,  1441,    93 }},
  {{   97,  1442,  1444 }},
  {{ 1442,   369,  1443 }},
  {{ 1444,  1442,  1443 }},
  {{ 1444,  1443,   371 }},
  {{  369,  1445,  1447 }},
  {{ 1445,    28,  1446 }},
  {{ 1447,  1445,  1446 }},
  {{ 1447,  1446,   370 }},
  {{  371,  1443,  1448 }},
  {{ 1443,   369,  1447 }},
  {{ 1448,  1443,  1447 }},
  {{ 1448,  1447,   370 }},
  {{  371,  1448,  1450 }},
  {{ 1448,   370,  1449 }},
  {{ 1450,  1448,  1449 }},
  {{ 1450,  1449,    98 }},
  {{   93,  1441,  1452 }},
  {{ 1441,   368,  1451 }},
  {{ 1452,  1441,  1451 }},
  {{ 1452,  1451,   372 }},
  {{  368,  1438,  1453 }},
  {{ 1438,    97,  1444 }},
  {{ 1453,  1438,  1444 }},
  {{ 1453,  1444,   371 }},
  {{  372,  1451,  1454 }},
  {{ 1451,   368,  1453 }},
  {{ 1454,  1451,  1453 }},
  {{ 1454,  1453,   371 }},
  {{  372,  1454,  1455 }},
  {{ 1454,   371,  1450 }},
  {{ 1455,  1454,  1450 }},
  {{ 1455,  1450,    98 }},
  {{   93,  1452,  1386 }},
  {{ 1452,   372,  1456 }},
  {{ 1386,  1452,  1456 }},
  {{ 1386,  1456,   355 }},
  {{  372,  1455,  1458 }},
  {{ 1455,    98,  1457 }},
  {{ 1458,  1455,  1457 }},
  {{ 1458,  1457,   373 }},
  {{  355,  1456,  1459 }},
  {{ 1456,   372,  1458 }},
  {{ 1459,  1456,  1458 }},
  {{ 1459,  1458,   373 }},
  {{  355,  1459,  1389 }},
  {{ 1459,   373,  1460 }},
  {{ 1389,  1459,  1460 }},
  {{ 1389,  1460,    27 }},
  {{   28,  1461,  1463 }},
  {{ 1461,   374,  1462 }},
  {{ 1463,  1461,  1462 }},
  {{ 1463,  1462,   376 }},
  {{  374,  1464,  1466 }},
  {{ 1464,    99,  1465 }},
  {{ 1466,  1464,  1465 }},
  {{ 1466,  1465,   375 }},
  {{  376,  1462,  1467 }},
  {{ 1462,   374,  1466 }},
  {{ 1467,  1462,  1466 }},
  {{ 1467,  1466,   375 }},
  {{  376,  1467,  1469 }},
  {{ 1467,   375,  1468 }},
  {{ 1469,  1467,  1468 }},
  {{ 1469,  1468,   101 }},
  {{   99,  1470,  1472 }},
  {{ 1470,   377,  1471 }},
  {{ 1472,  1470,  1471 }},
  {{ 1472,  1471,   379 }},
  {{  377,  1473,  1475 }},
  {{ 1473,    11,  1474 }},
  {{ 1475,  1473,  1474 }},
  {{ 1475,  1474,   378 }},
  {{  379,  1471,  1476 }},
  {{ 1471,   377,  1475 }},
  {{ 1476,  1471,  1475 }},
  {{ 1476,  1475,   378 }},
  {{  379,  1476,  1478 }},
  {{ 1476,   378,  1477 }},
  {{ 1478,  1476,  1477 }},
  {{ 1478,  1477,   100 }},
  {{  101,  1468,  1480 }},
  {{ 1468,   375,  1479 }},
  {{ 1480,  1468,  1479 }},
  {{ 1480,  1479,   380 }},
  {{  375,  1465,  1481 }},
  {{ 1465,    99,  1472 }},
  {{ 1481,  1465,  1472 }},
  {{ 1481,  1472,   379 }},
  {{  380,  1479,  1482 }},
  {{ 1479,   375,  1481 }},
  {{ 1482,  1479,  1481 }},
  {{ 1482,  1481,   379 }},
  {{  380,  1482,  1483 }},
  {{ 1482,   379,  1478 }},
  {{ 1483,  1482,  1478 }},
  {{ 1483,  1478,   100 }},
  {{  101,  1480,  1485 }},
  {{ 1480,   380,  1484 }},
  {{ 1485,  1480,  1484 }},
  {{ 1485,  1484,   382 }},
  {{  380,  1483,  1487 }},
  {{ 1483,   100,  1486 }},
  {{ 1487,  1483,  1486 }},
  {{ 1487,  1486,   381 }},
  {{  382,  1484,  1488 }},
  {{ 1484,   380,  1487 }},
  {{ 1488,  1484,  1487 }},
  {{ 1488,  1487,   381 }},
  {{  382,  1488,  1490 }},
  {{ 1488,   381,  1489 }},
  {{ 1490,  1488,  1489 }},
  {{ 1490,  1489,    29 }},
  {{   27,  1460,  1492 }},
  {{ 1460,   373,  1491 }},
  {{ 1492,  1460,  1491 }},
  {{ 1492,  1491,   384 }},
  {{  373,  1457,  1494 }},
  {{ 1457,    98,  1493 }},
  {{ 1494,  1457,  1493 }},
  {{ 1494,  1493,   383 }},
  {{  384,  1491,  1495 }},
  {{ 1491,   373,  1494 }},
  {{ 1495,  1491,  1494 }},
  {{ 1495,  1494,   383 }},
  {{  384,  1495,  1497 }},
  {{ 1495,   383,  1496 }},
  {{ 1497,  1495,  1496 }},
  {{ 1497,  1496,   102 }},
  {{   98,  1449,  1499 }},
  {{ 1449,   370,  1498 }},
  {{ 1499,  1449,  1498 }},
  {{ 1499,  1498,   385 }},
  {{  370,  1446,  1500 }},
  {{ 1446,    28,  1463 }},
  {{ 1500,  1446,  1463 }},
  {{ 1500,  1463,   376 }},
  {{  385,  1498,  1501 }},
  {{ 1498,   370,  1500 }},
  {{ 1501,  1498,  1500 }},
  {{ 1501,  1500,   376 }},
  {{  385,  1501,  1502 }},
  {{ 1501,   376,  1469 }},
  {{ 1502,  1501,  1469 }},
  {{ 1502,  1469,   101 }},
  {{  102,  1496,  1504 }},
  {{ 1496,   383,  1503 }},
  {{ 1504,  1496,  1503 }},
  {{ 1504,  1503,   386 }},
  {{  383,  1493,  1505 }},
  {{ 1493,    98,  1499 }},
  {{ 1505,  1493,  1499 }},
  {{ 1505,  1499,   385 }},
  {{  386,  1503,  1506 }},
  {{ 1503,   383,  1505 }},
  {{ 1506,  1503,  1505 }},
  {{ 1506,  1505,   385 }},
  {{  386,  1506,  1507 }},
  {{ 1506,   385,  1502 }},
  {{ 1507,  1506,  1502 }},
  {{ 1507,  1502,   101 }},
  {{  102,  1504,  1509 }},
  {{ 1504,   386,  1508 }},
  {{ 1509,  1504,  1508 }},
  {{ 1509,  1508,   387 }},
  {{  386,  1507,  1510 }},
  {{ 1507,   101,  1485 }},
  {{ 1510,  1507,  1485 }},
  {{ 1510,  1485,   382 }},
  {{  387,  1508,  1511 }},
  {{ 1508,   386,  1510 }},
  {{ 1511,  1508,  1510 }},
  {{ 1511,  1510,   382 }},
  {{  387,  1511,  1512 }},
  {{ 1511,   382,  1490 }},
  {{ 1512,  1511,  1490 }},
  {{ 1512,  1490,    29 }},
  {{   27,  1492,  1420 }},
  {{ 1492,   384,  1513 }},
  {{ 1420,  1492,  1513 }},
  {{ 1420,  1513,   363 }},
  {{  384,  1497,  1515 }},
  {{ 1497,   102,  1514 }},
  {{ 1515,  1497,  1514 }},
  {{ 1515,  1514,   388 }},
  {{  363,  1513,  1516 }},
  {{ 1513,   384,  1515 }},
  {{ 1516,  1513,  1515 }},
  {{ 1516,  1515,   388 }},
  {{  363,  1516,  1423 }},
  {{ 1516,   388,  1517 }},
  {{ 1423,  1516,  1517 }},
  {{ 1423,  1517,    96 }},
  {{  102,  1509,  1519 }},
  {{ 1509,   387,  1518 }},
  {{ 1519,  1509,  1518 }},
  {{ 1519,  1518,   390 }},
  {{  387,  1512,  1521 }},
  {{ 1512,    29,  1520 }},
  {{ 1521,  1512,  1520 }},
  {{ 1521,  1520,   389 }},
  {{  390,  1518,  1522 }},
  {{ 1518,   387,  1521 }},
  {{ 1522,  1518,  1521 }},
  {{ 1522,  1521,   389 }},
  {{  390,  1522,  1524 }},
  {{ 1522,   389,  1523 }},
  {{ 1524,  1522,  1523 }},
  {{ 1524,  1523,   103 }},
  {{   96,  1517,  1526 }},
  {{ 1517,   388,  1525 }},
  {{ 1526,  1517,  1525 }},
  {{ 1526,  1525,   391 }},
  {{  388,  1514,  1527 }},
  {{ 1514,   102,  1519 }},
  {{ 1527,  1514,  1519 }},
  {{ 1527,  1519,   390 }},
  {{  391,  1525,  1528 }},
  {{ 1525,   388,  1527 }},
  {{ 1528,  1525,  1527 }},
  {{ 1528,  1527,   390 }},
  {{  391,  1528,  1529 }},
  {{ 1528,   390,  1524 }},
  {{ 1529,  1528,  1524 }},
  {{ 1529,  1524,   103 }},
  {{   96,  1526,  1431 }},
  {{ 1526,   391,  1530 }},
  {{ 1431,  1526,  1530 }},
  {{ 1431,  1530,   366 }},
  {{  391,  1529,  1532 }},
  {{ 1529,   103,  1531 }},
  {{ 1532,  1529,  1531 }},
  {{ 1532,  1531,   392 }},
  {{  366,  1530,  1533 }},
  {{ 1530,   391,  1532 }},
  {{ 1533,  1530,  1532 }},
  {{ 1533,  1532,   392 }},
  {{  366,  1533,  1434 }},
  {{ 1533,   392,  1534 }},
  {{ 1434,  1533,  1534 }},
  {{ 1434,  1534,    10 }},
  {{    5,   850,  1435 }},
  {{  850,   218,  1535 }},
  {{ 1435,   850,  1535 }},
  {{ 1435,  1535,   367 }},
  {{  218,   847,  1537 }},
  {{  847,    58,  1536 }},
  {{ 1537,   847,  1536 }},
  {{ 1537,  1536,   393 }},
  {{  367,  1535,  1538 }},
  {{ 1535,   218,  1537 }},
  {{ 1538,  1535,  1537 }},
  {{ 1538,  1537,   393 }},
  {{  367,  1538,  1437 }},
  {{ 1538,   393,  1539 }},
  {{ 1437,  1538,  1539 }},
  {{ 1437,  1539,    97 }},
  {{   58,   839,  1541 }},
  {{  839,   215,  1540 }},
  {{ 1541,   839,  1540 }},
  {{ 1541,  1540,   395 }},
  {{  215,   836,  1543 }},
  {{  836,    17,  1542 }},
  {{ 1543,   836,  1542 }},
  {{ 1543,  1542,   394 }},
  {{  395,  1540,  1544 }},
  {{ 1540,   215,  1543 }},
  {{ 1544,  1540,  1543 }},
  {{ 1544,  1543,   394 }},
  {{  395,  1544,  1546 }},
  {{ 1544,   394,  1545 }},
  {{ 1546,  1544,  1545 }},
  {{ 1546,  1545,   104 }},
  {{   97,  1539,  1548 }},
  {{ 1539,   393,  1547 }},
  {{ 1548,  1539,  1547 }},
  {{ 1548,  1547,   396 }},
  {{  393,  1536,  1549 }},
  {{ 1536,    58,  1541 }},
  {{ 1549,  1536,  1541 }},
  {{ 1549,  1541,   395 }},
  {{  396,  1547,  1550 }},
  {{ 1547,   393,  1549 }},
  {{ 1550,  1547,  1549 }},
  {{ 1550,  1549,   395 }},
  {{  396,  1550,  1551 }},
  {{ 1550,   395,  1546 }},
  {{ 1551,  1550,  1546 }},
  {{ 1551,  1546,   104 }},
  {{   97,  1548,  1442 }},
  {{ 1548,   396,  1552 }},
  {{ 1442,  1548,  1552 }},
  {{ 1442,  1552,   369 }},
  {{  396,  1551,  1554 }},
  {{ 1551,   104,  1553 }},
  {{ 1554,  1551,  1553 }},
  {{ 1554,  1553,   397 }},
  {{  369,  1552,  1555 }},
  {{ 1552,   396,  1554 }},
  {{ 1555,  1552,  1554 }},
  {{ 1555,  1554,   397 }},
  {{  369,  1555,  1445 }},
  {{ 1555,   397,  1556 }},
  {{ 1445,  1555,  1556 }},
  {{ 1445,  1556,    28 }},
  {{   17,   805,  1558 }},
  {{  805,   207,  1557 }},
  {{ 1558,   805,  1557 }},
  {{ 1558,  1557,   399 }},
  {{  207,   802,  1560 }},
  {{  802,    55,  1559 }},
  {{ 1560,   802,  1559 }},
  {{ 1560,  1559,   398 }},
  {{  399,  1557,  1561 }},
  {{ 1557,   207,  1560 }},
  {{ 1561,  1557,  1560 }},
  {{ 1561,  1560,   398 }},
  {{  399,  1561,  1563 }},
  {{ 1561,   398,  1562 }},
  {{ 1563,  1561,  1562 }},
  {{ 1563,  1562,   106 }},
  {{   55,   793,  1565 }},
  {{  793,   204,  1564 }},
  {{ 1565,   793,  1564 }},
  {{ 1565,  1564,   401 }},
  {{  204,   790,  1567 }},
  {{  790,     6,  1566 }},
  {{ 1567,   790,  1566 }},
  {{ 1567,  1566,   400 }},
  {{  401,  1564,  1568 }},
  {{ 1564,   204,  1567 }},
  {{ 1568,  1564,  1567 }},
  {{ 1568,  1567,   400 }},
  {{  401,  1568,  1570 }},
  {{ 1568,   400,  1569 }},
  {{ 1570,  1568,  1569 }},
  {{ 1570,  1569,   105 }},
  {{  106,  1562,  1572 }},
  {{ 1562,   398,  1571 }},
  {{ 1572,  1562,  1571 }},
  {{ 1572,  1571,   402 }},
  {{  398,  1559,  1573 }},
  {{ 1559,    55,  1565 }},
  {{ 1573,  1559,  1565 }},
  {{ 1573,  1565,   401 }},
  {{  402,  1571,  1574 }},
  {{ 1571,   398,  1573 }},
  {{ 1574,  1571,  1573 }},
  {{ 1574,  1573,   401 }},
  {{  402,  1574,  1575 }},
  {{ 1574,   401,  1570 }},
  {{ 1575,  1574,  1570 }},
  {{ 1575,  1570,   105 }},
  {{  106,  1572,  1577 }},
  {{ 1572,   402,  1576 }},
  {{ 1577,  1572,  1576 }},
  {{ 1577,  1576,   404 }},
  {{  402,  1575,  1579 }},
  {{ 1575,   105,  1578 }},
  {{ 1579,  1575,  1578 }},
  {{ 1579,  1578,   403 }},
  {{  404,  1576,  1580 }},
  {{ 1576,   402,  1579 }},
  {{ 1580,  1576,  1579 }},
  {{ 1580,  1579,   403 }},
  {{  404,  1580,  1582 }},
  {{ 1580,   403,  1581 }},
  {{ 1582,  1580,  1581 }},
  {{ 1582,  1581,    30 }},
  {{   28,  1556,  1584 }},
  {{ 1556,   397,  1583 }},
  {{ 1584,  1556,  1583 }},
  {{ 1584,  1583,   406 }},
  {{  397,  1553,  1586 }},
  {{ 1553,   104,  1585 }},
  {{ 1586,  1553,  1585 }},
  {{ 1586,  1585,   405 }},
  {{  406,  1583,  1587 }},
  {{ 1583,   397,  1586 }},
  {{ 1587,  1583,  1586 }},
  {{ 1587,  1586,   405 }},
  {{  406,  1587,  1589 }},
  {{ 1587,   405,  1588 }},
  {{ 1589,  1587,  1588 }},
  {{ 1589,  1588,   107 }},
  {{  104,  1545,  1591 }},
  {{ 1545,   394,  1590 }},
  {{ 1591,  1545,  1590 }},
  {{ 1591,  1590,   407 }},
  {{  394,  1542,  1592 }},
  {{ 1542,    17,  1558 }},
  {{ 1592,  1542,  1558 }},
  {{ 1592,  1558,   399 }},
  {{  407,  1590,  1593 }},
  {{ 1590,   394,  1592 }},
  {{ 1593,  1590,  1592 }},
  {{ 1593,  1592,   399 }},
  {{  407,  1593,  1594 }},
  {{ 1593,   399,  1563 }},
  {{ 1594,  1593,  1563 }},
  {{ 1594,  1563,   106 }},
  {{  107,  1588,  1596 }},
  {{ 1588,   405,  1595 }},
  {{ 1596,  1588,  1595 }},
  {{ 1596,  1595,   408 }},
  {{  405,  1585,  1597 }},
  {{ 1585,   104,  1591 }},
  {{ 1597,  1585,  1591 }},
  {{ 1597,  1591,   407 }},
  {{  408,  1595,  1598 }},
  {{ 1595,   405,  1597 }},
  {{ 1598,  1595,  1597 }},
  {{ 1598,  1597,   407 }},
  {{  408,  1598,  1599 }},
  {{ 1598,   407,  1594 }},
  {{ 1599,  1598,  1594 }},
  {{ 1599,  1594,   106 }},
  {{  107,  1596,  1601 }},
  {{ 1596,   408,  1600 }},
  {{ 1601,  1596,  1600 }},
  {{ 1601,  1600,   409 }},
  {{  408,  1599,  1602 }},
  {{ 1599,   106,  1577 }},
  {{ 1602,  1599,  1577 }},
  {{ 1602,  1577,   404 }},
  {{  409,  1600,  1603 }},
  {{ 1600,   408,  1602 }},
  {{ 1603,  1600,  1602 }},
  {{ 1603,  1602,   404 }},
  {{  409,  1603,  1604 }},
  {{ 1603,   404,  1582 }},
  {{ 1604,  1603,  1582 }},
  {{ 1604,  1582,    30 }},
  {{   28,  1584,  1461 }},
  {{ 1584,   406,  1605 }},
  {{ 1461,  1584,  1605 }},
  {{ 1461,  1605,   374 }},
  {{  406,  1589,  1607 }},
  {{ 1589,   107,  1606 }},
  {{ 1607,  1589,  1606 }},
  {{ 1607,  1606,   410 }},
  {{  374,  1605,  1608 }},
  {{ 1605,   406,  1607 }},
  {{ 1608,  1605,  1607 }},
  {{ 1608,  1607,   410 }},
  {{  374,  1608,  1464 }},
  {{ 1608,   410,  1609 }},
  {{ 1464,  1608,  1609 }},
  {{ 1464,  1609,    99 }},
  {{  107,  1601,  1611 }},
  {{ 1601,   409,  1610 }},
  {{ 1611,  1601,  1610 }},
  {{ 1611,  1610,   412 }},
  {{  409,  1604,  1613 }},
  {{ 1604,    30,  1612 }},
  {{ 1613,  1604,  1612 }},
  {{ 1613,  1612,   411 }},
  {{  412,  1610,  1614 }},
  {{ 1610,   409,  1613 }},
  {{ 1614,  1610,  1613 }},
  {{ 1614,  1613,   411 }},
  {{  412,  1614,  1616 }},
  {{ 1614,   411,  1615 }},
  {{ 1616,  1614,  1615 }},
  {{ 1616,  1615,   108 }},
  {{   99,  1609,  1618 }},
  {{ 1609,   410,  1617 }},
  {{ 1618,  1609,  1617 }},
  {{ 1618,  1617,   413 }},
  {{  410,  1606,  1619 }},
  {{ 1606,   107,  1611 }},
  {{ 1619,  1606,  1611 }},
  {{ 1619,  1611,   412 }},
  {{  413,  1617,  1620 }},
  {{ 1617,   410,  1619 }},
  {{ 1620,  1617,  1619 }},
  {{ 1620,  1619,   412 }},
  {{  413,  1620,  1621 }},
  {{ 1620,   412,  1616 }},
  {{ 1621,  1620,  1616 }},
  {{ 1621,  1616,   108 }},
  {{   99,  1618,  1470 }},
  {{ 1618,   413,  1622 }},
  {{ 1470,  1618,  1622 }},
  {{ 1470,  1622,   377 }},
  {{  413,  1621,  1624 }},
  {{ 1621,   108,  1623 }},
  {{ 1624,  1621,  1623 }},
  {{ 1624,  1623,   414 }},
  {{  377,  1622,  1625 }},
  {{ 1622,   413,  1624 }},
  {{ 1625,  1622,  1624 }},
  {{ 1625,  1624,   414 }},
  {{  377,  1625,  1473 }},
  {{ 1625,   414,  1626 }},
  {{ 1473,  1625,  1626 }},
  {{ 1473,  1626,    11 }},
  {{    6,  1627,  1566 }},
  {{ 1627,   415,  1628 }},
  {{ 1566,  1627,  1628 }},
  {{ 1566,  1628,   400 }},
  {{  415,  1629,  1631 }},
  {{ 1629,   109,  1630 }},
  {{ 1631,  1629,  1630 }},
  {{ 1631,  1630,   416 }},
  {{  400,  1628,  1632 }},
  {{ 1628,   415,  1631 }},
  {{ 1632,  1628,  1631 }},
  {{ 1632,  1631,   416 }},
  {{  400,  1632,  1569 }},
  {{ 1632,   416,  1633 }},
  {{ 1569,  1632,  1633 }},
  {{ 1569,  1633,   105 }},
  {{  109,  1634,  1636 }},
  {{ 1634,   417,  1635 }},
  {{ 1636,  1634,  1635 }},
  {{ 1636,  1635,   419 }},
  {{  417,  1637,  1639 }},
  {{ 1637,    31,  1638 }},
  {{ 1639,  1637,  1638 }},
  {{ 1639,  1638,   418 }},
  {{  419,  1635,  1640 }},
  {{ 1635,   417,  1639 }},
  {{ 1640,  1635,  1639 }},
  {{ 1640,  1639,   418 }},
  {{  419,  1640,  1642 }},
  {{ 1640,   418,  1641 }},
  {{ 1642,  1640,  1641 }},
  {{ 1642,  1641,   110 }},
  {{  105,  1633,  1644 }},
  {{ 1633,   416,  1643 }},
  {{ 1644,  1633,  1643 }},
  {{ 1644,  1643,   420 }},
  {{  416,  1630,  1645 }},
  {{ 1630,   109,  1636 }},
  {{ 1645,  1630,  1636 }},
  {{ 1645,  1636,   419 }},
  {{  420,  1643,  1646 }},
  {{ 1643,   416,  1645 }},
  {{ 1646,  1643,  1645 }},
  {{ 1646,  1645,   419 }},
  {{  420,  1646,  1647 }},
  {{ 1646,   419,  1642 }},
  {{ 1647,  1646,  1642 }},
  {{ 1647,  1642,   110 }},
  {{  105,  1644,  1578 }},
  {{ 1644,   420,  1648 }},
  {{ 1578,  1644,  1648 }},
  {{ 1578,  1648,   403 }},
  {{  420,  1647,  1650 }},
  {{ 1647,   110,  1649 }},
  {{ 1650,  1647,  1649 }},
  {{ 1650,  1649,   421 }},
  {{  403,  1648,  1651 }},
  {{ 1648,   420,  1650 }},
  {{ 1651,  1648,  1650 }},
  {{ 1651,  1650,   421 }},
  {{  403,  1651,  1581 }},
  {{ 1651,   421,  1652 }},
  {{ 1581,  1651,  1652 }},
  {{ 1581,  1652,    30 }},
  {{   31,  1653,  1655 }},
  {{ 1653,   422,  1654 }},
  {{ 1655,  1653,  1654 }},
  {{ 1655,  1654,   424 }},
  {{  422,  1656,  1658 }},
  {{ 1656,   111,  1657 }},
  {{ 1658,  1656,  1657 }},
  {{ 1658,  1657,   423 }},
  {{  424,  1654,  1659 }},
  {{ 1654,   422,  1658 }},
  {{ 1659,  1654,  1658 }},
  {{ 1659,  1658,   423 }},
  {{  424,  1659,  1661 }},
  {{ 1659,   423,  1660 }},
  {{ 1661,  1659,  1660 }},
  {{ 1661,  1660,   113 }},
  {{  111,  1662,  1664 }},
  {{ 1662,   425,  1663 }},
  {{ 1664,  1662,  1663 }},
  {{ 1664,  1663,   427 }},
  {{  425,  1665,  1667 }},
  {{ 1665,     7,  1666 }},
  {{ 1667,  1665,  1666 }},
  {{ 1667,  1666,   426 }},
  {{  427,  1663,  1668 }},
  {{ 1663,   425,  1667 }},
  {{ 1668,  1663,  1667 }},
  {{ 1668,  1667,   426 }},
  {{  427,  1668,  1670 }},
  {{ 1668,   426,  1669 }},
  {{ 1670,  1668,  1669 }},
  {{ 1670,  1669,   112 }},
  {{  113,  1660,  1672 }},
  {{ 1660,   423,  1671 }},
  {{ 1672,  1660,  1671 }},
  {{ 1672,  1671,   428 }},
  {{  423,  1657,  1673 }},
  {{ 1657,   111,  1664 }},
  {{ 1673,  1657,  1664 }},
  {{ 1673,  1664,   427 }},
  {{  428,  1671,  1674 }},
  {{ 1671,   423,  1673 }},
  {{ 1674,  1671,  1673 }},
  {{ 1674,  1673,   427 }},
  {{  428,  1674,  1675 }},
  {{ 1674,   427,  1670 }},
  {{ 1675,  1674,  1670 }},
  {{ 1675,  1670,   112 }},
  {{  113,  1672,  1677 }},
  {{ 1672,   428,  1676 }},
  {{ 1677,  1672,  1676 }},
  {{ 1677,  1676,   430 }},
  {{  428,  1675,  1679 }},
  {{ 1675,   112,  1678 }},
  {{ 1679,  1675,  1678 }},
  {{ 1679,  1678,   429 }},
  {{  430,  1676,  1680 }},
  {{ 1676,   428,  1679 }},
  {{ 1680,  1676,  1679 }},
  {{ 1680,  1679,   429 }},
  {{  430,  1680,  1682 }},
  {{ 1680,   429,  1681 }},
  {{ 1682,  1680,  1681 }},
  {{ 1682,  1681,    32 }},
  {{   30,  1652,  1684 }},
  {{ 1652,   421,  1683 }},
  {{ 1684,  1652,  1683 }},
  {{ 1684,  1683,   432 }},
  {{  421,  1649,  1686 }},
  {{ 1649,   110,  1685 }},
  {{ 1686,  1649,  1685 }},
  {{ 1686,  1685,   431 }},
  {{  432,  1683,  1687 }},
  {{ 1683,   421,  1686 }},
  {{ 1687,  1683,  1686 }},
  {{ 1687,  1686,   431 }},
  {{  432,  1687,  1689 }},
  {{ 1687,   431,  1688 }},
  {{ 1689,  1687,  1688 }},
  {{ 1689,  1688,   114 }},
  {{  110,  1641,  1691 }},
  {{ 1641,   418,  1690 }},
  {{ 1691,  1641,  1690 }},
  {{ 1691,  1690,   433 }},
  {{  418,  1638,  1692 }},
  {{ 1638,    31,  1655 }},
  {{ 1692,  1638,  1655 }},
  {{ 1692,  1655,   424 }},
  {{  433,  1690,  1693 }},
  {{ 1690,   418,  1692 }},
  {{ 1693,  1690,  1692 }},
  {{ 1693,  1692,   424 }},
  {{  433,  1693,  1694 }},
  {{ 1693,   424,  1661 }},
  {{ 1694,  1693,  1661 }},
  {{ 1694,  1661,   113 }},
  {{  114,  1688,  1696 }},
  {{ 1688,   431,  1695 }},
  {{ 1696,  1688,  1695 }},
  {{ 1696,  1695,   434 }},
  {{  431,  1685,  1697 }},
  {{ 1685,   110,  1691 }},
  {{ 1697,  1685,  1691 }},
  {{ 1697,  1691,   433 }},
  {{  434,  1695,  1698 }},
  {{ 1695,   431,  1697 }},
  {{ 1698,  1695,  1697 }},
  {{ 1698,  1697,   433 }},
  {{  434,  1698,  1699 }},
  {{ 1698,   433,  1694 }},
  {{ 1699,  1698,  1694 }},
  {{ 1699,  1694,   113 }},
  {{  114,  1696,  1701 }},
  {{ 1696,   434,  1700 }},
  {{ 1701,  1696,  1700 }},
  {{ 1701,  1700,   435 }},
  {{  434,  1699,  1702 }},
  {{ 1699,   113,  1677 }},
  {{ 1702,  1699,  1677 }},
  {{ 1702,  1677,   430 }},
  {{  435,  1700,  1703 }},
  {{ 1700,   434,  1702 }},
  {{ 1703,  1700,  1702 }},
  {{ 1703,  1702,   430 }},
  {{  435,  1703,  1704 }},
  {{ 1703,   430,  1682 }},
  {{ 1704,  1703,  1682 }},
  {{ 1704,  1682,    32 }},
  {{   30,  1684,  1612 }},
  {{ 1684,   432,  1705 }},
  {{ 1612,  1684,  1705 }},
  {{ 1612,  1705,   411 }},
  {{  432,  1689,  1707 }},
  {{ 1689,   114,  1706 }},
  {{ 1707,  1689,  1706 }},
  {{ 1707,  1706,   436 }},
  {{  411,  1705,  1708 }},
  {{ 1705,   432,  1707 }},
  {{ 1708,  1705,  1707 }},
  {{ 1708,  1707,   436 }},
  {{  411,  1708,  1615 }},
  {{ 1708,   436,  1709 }},
  {{ 1615,  1708,  1709 }},
  {{ 1615,  1709,   108 }},
  {{  114,  1701,  1711 }},
  {{ 1701,   435,  1710 }},
  {{ 1711,  1701,  1710 }},
  {{ 1711,  1710,   438 }},
  {{  435,  1704,  1713 }},
  {{ 1704,    32,  1712 }},
  {{ 1713,  1704,  1712 }},
  {{ 1713,  1712,   437 }},
  {{  438,  1710,  1714 }},
  {{ 1710,   435,  1713 }},
  {{ 1714,  1710,  1713 }},
  {{ 1714,  1713,   437 }},
  {{  438,  1714,  1716 }},
  {{ 1714,   437,  1715 }},
  {{ 1716,  1714,  1715 }},
  {{ 1716,  1715,   115 }},
  {{  108,  1709,  1718 }},
  {{ 1709,   436,  1717 }},
  {{ 1718,  1709,  1717 }},
  {{ 1718,  1717,   439 }},
  {{  436,  1706,  1719 }},
  {{ 1706,   114,  1711 }},
  {{ 1719,  1706,  1711 }},
  {{ 1719,  1711,   438 }},
  {{  439,  1717,  1720 }},
  {{ 1717,   436,  1719 }},
  {{ 1720,  1717,  1719 }},
  {{ 1720,  1719,   438 }},
  {{  439,  1720,  1721 }},
  {{ 1720,   438,  1716 }},
  {{ 1721,  1720,  1716 }},
  {{ 1721,  1716,   115 }},
  {{  108,  1718,  1623 }},
  {{ 1718,   439,  1722 }},
  {{ 1623,  1718,  1722 }},
  {{ 1623,  1722,   414 }},
  {{  439,  1721,  1724 }},
  {{ 1721,   115,  1723 }},
  {{ 1724,  1721,  1723 }},
  {{ 1724,  1723,   440 }},
  {{  414,  1722,  1725 }},
  {{ 1722,   439,  1724 }},
  {{ 1725,  1722,  1724 }},
  {{ 1725,  1724,   440 }},
  {{  414,  1725,  1626 }},
  {{ 1725,   440,  1726 }},
  {{ 1626,  1725,  1726 }},
  {{ 1626,  1726,    11 }},
  {{    6,   950,  1627 }},
  {{  950,   244,  1727 }},
  {{ 1627,   950,  1727 }},
  {{ 1627,  1727,   415 }},
  {{  244,   947,  1729 }},
  {{  947,    65,  1728 }},
  {{ 1729,   947,  1728 }},
  {{ 1729,  1728,   441 }},
  {{  415,  1727,  1730 }},
  {{ 1727,   244,  1729 }},
  {{ 1730,  1727,  1729 }},
  {{ 1730,  1729,   441 }},
  {{  415,  1730,  1629 }},
  {{ 1730,   441,  1731 }},
  {{ 1629,  1730,  1731 }},
  {{ 1629,  1731,   109 }},
  {{   65,   939,  1733 }},
  {{  939,   241,  1732 }},
  {{ 1733,   939,  1732 }},
  {{ 1733,  1732,   443 }},
  {{  241,   936,  1735 }},
  {{  936,    19,  1734 }},
  {{ 1735,   936,  1734 }},
  {{ 1735,  1734,   442 }},
  {{  443,  1732,  1736 }},
  {{ 1732,   241,  1735 }},
  {{ 1736,  1732,  1735 }},
  {{ 1736,  1735,   442 }},
  {{  443,  1736,  1738 }},
  {{ 1736,   442,  1737 }},
  {{ 1738,  1736,  1737 }},
  {{ 1738,  1737,   116 }},
  {{  109,  1731,  1740 }},
  {{ 1731,   441,  1739 }},
  {{ 1740,  1731,  1739 }},
  {{ 1740,  1739,   444 }},
  {{  441,  1728,  1741 }},
  {{ 1728,    65,  1733 }},
  {{ 1741,  1728,  1733 }},
  {{ 1741,  1733,   443 }},
  {{  444,  1739,  1742 }},
  {{ 1739,   441,  1741 }},
  {{ 1742,  1739,  1741 }},
  {{ 1742,  1741,   443 }},
  {{  444,  1742,  1743 }},
  {{ 1742,   443,  1738 }},
  {{ 1743,  1742,  1738 }},
  {{ 1743,  1738,   116 }},
  {{  109,  1740,  1634 }},
  {{ 1740,   444,  1744 }},
  {{ 1634,  1740,  1744 }},
  {{ 1634,  1744,   417 }},
  {{  444,  1743,  1746 }},
  {{ 1743,   116,  1745 }},
  {{ 1746,  1743,  1745 }},
  {{ 1746,  1745,   445 }},
  {{  417,  1744,  1747 }},
  {{ 1744,   444,  1746 }},
  {{ 1747,  1744,  1746 }},
  {{ 1747,  1746,   445 }},
  {{  417,  1747,  1637 }},
  {{ 1747,   445,  1748 }},
  {{ 1637,  1747,  1748 }},
  {{ 1637,  1748,    31 }},
  {{   19,   905,  1750 }},
  {{  905,   233,  1749 }},
  {{ 1750,   905,  1749 }},
  {{ 1750,  1749,   447 }},
  {{  233,   902,  1752 }},
  {{  902,    62,  1751 }},
  {{ 1752,   902,  1751 }},
  {{ 1752,  1751,   446 }},
  {{  447,  1749,  1753 }},
  {{ 1749,   233,  1752 }},
  {{ 1753,  1749,  1752 }},
  {{ 1753,  1752,   446 }},
  {{  447,  1753,  1755 }},
  {{ 1753,   446,  1754 }},
  {{ 1755,  1753,  1754 }},
  {{ 1755,  1754,   118 }},
  {{   62,   893,  1757 }},
  {{  893,   230,  1756 }},
  {{ 1757,   893,  1756 }},
  {{ 1757,  1756,   449 }},
  {{  230,   890,  1759 }},
  {{  890,     2,  1758 }},
  {{ 1759,   890,  1758 }},
  {{ 1759,  1758,   448 }},
  {{  449,  1756,  1760 }},
  {{ 1756,   230,  1759 }},
  {{ 1760,  1756,  1759 }},
  {{ 1760,  1759,   448 }},
  {{  449,  1760,  1762 }},
  {{ 1760,   448,  1761 }},
  {{ 1762,  1760,  1761 }},
  {{ 1762,  1761,   117 }},
  {{  118,  1754,  1764 }},
  {{ 1754,   446,  1763 }},
  {{ 1764,  1754,  1763 }},
  {{ 1764,  1763,   450 }},
  {{  446,  1751,  1765 }},
  {{ 1751,    62,  1757 }},
  {{ 1765,  1751,  1757 }},
  {{ 1765,  1757,   449 }},
  {{  450,  1763,  1766 }},
  {{ 1763,   446,  1765 }},
  {{ 1766,  1763,  1765 }},
  {{ 1766,  1765,   449 }},
  {{  450,  1766,  1767 }},
  {{ 1766,   449,  1762 }},
  {{ 1767,  1766,  1762 }},
  {{ 1767,  1762,   117 }},
  {{  118,  1764,  1769 }},
  {{ 1764,   450,  1768 }},
  {{ 1769,  1764,  1768 }},
  {{ 1769,  1768,   452 }},
  {{  450,  1767,  1771 }},
  {{ 1767,   117,  1770 }},
  {{ 1771,  1767,  1770 }},
  {{ 1771,  1770,   451 }},
  {{  452,  1768,  1772 }},
  {{ 1768,   450,  1771 }},
  {{ 1772,  1768,  1771 }},
  {{ 1772,  1771,   451 }},
  {{  452,  1772,  1774 }},
  {{ 1772,   451,  1773 }},
  {{ 1774,  1772,  1773 }},
  {{ 1774,  1773,    33 }},
  {{   31,  1748,  1776 }},
  {{ 1748,   445,  1775 }},
  {{ 1776,  1748,  1775 }},
  {{ 1776,  1775,   454 }},
  {{  445,  1745,  1778 }},
  {{ 1745,   116,  1777 }},
  {{ 1778,  1745,  1777 }},
  {{ 1778,  1777,   453 }},
  {{  454,  1775,  1779 }},
  {{ 1775,   445,  1778 }},
  {{ 1779,  1775,  1778 }},
  {{ 1779,  1778,   453 }},
  {{  454,  1779,  1781 }},
  {{ 1779,   453,  1780 }},
  {{ 1781,  1779,  1780 }},
  {{ 1781,  1780,   119 }},
  {{  116,  1737,  1783 }},
  {{ 1737,   442,  1782 }},
  {{ 1783,  1737,  1782 }},
  {{ 1783,  1782,   455 }},
  {{  442,  1734,  1784 }},
  {{ 1734,    19,  1750 }},
  {{ 1784,  1734,  1750 }},
  {{ 1784,  1750,   447 }},
  {{  455,  1782,  1785 }},
  {{ 1782,   442,  1784 }},
  {{ 1785,  1782,  1784 }},
  {{ 1785,  1784,   447 }},
  {{  455,  1785,  1786 }},
  {{ 1785,   447,  1755 }},
  {{ 1786,  1785,  1755 }},
  {{ 1786,  1755,   118 }},
  {{  119,  1780,  1788 }},
  {{ 1780,   453,  1787 }},
  {{ 1788,  1780,  1787 }},
  {{ 1788,  1787,   456 }},
  {{  453,  1777,  1789 }},
  {{ 1777,   116,  1783 }},
  {{ 1789,  1777,  1783 }},
  {{ 1789,  1783,   455 }},
  {{  456,  1787,  1790 }},
  {{ 1787,   453,  1789 }},
  {{ 1790,  1787,  1789 }},
  {{ 1790,  1789,   455 }},
  {{  456,  1790,  1791 }},
  {{ 1790,   455,  1786 }},
  {{ 1791,  1790,  1786 }},
  {{ 1791,  1786,   118 }},
  {{  119,  1788,  1793 }},
  {{ 1788,   456,  1792 }},
  {{ 1793,  1788,  1792 }},
  {{ 1793,  1792,   457 }},
  {{  456,  1791,  1794 }},
  {{ 1791,   118,  1769 }},
  {{ 1794,  1791,  1769 }},
  {{ 1794,  1769,   452 }},
  {{  457,  1792,  1795 }},
  {{ 1792,   456,  1794 }},
  {{ 1795,  1792,  1794 }},
  {{ 1795,  1794,   452 }},
  {{  457,  1795,  1796 }},
  {{ 1795,   452,  1774 }},
  {{ 1796,  1795,  1774 }},
  {{ 1796,  1774,    33 }},
  {{   31,  1776,  1653 }},
  {{ 1776,   454,  1797 }},
  {{ 1653,  1776,  1797 }},
  {{ 1653,  1797,   422 }},
  {{  454,  1781,  1799 }},
  {{ 1781,   119,  1798 }},
  {{ 1799,  1781,  1798 }},
  {{ 1799,  1798,   458 }},
  {{  422,  1797,  1800 }},
  {{ 1797,   454,  1799 }},
  {{ 1800,  1797,  1799 }},
  {{ 1800,  1799,   458 }},
  {{  422,  1800,  1656 }},
  {{ 1800,   458,  1801 }},
  {{ 1656,  1800,  1801 }},
  {{ 1656,  1801,   111 }},
  {{  119,  1793,  1803 }},
  {{ 1793,   457,  1802 }},
  {{ 1803,  1793,  1802 }},
  {{ 1803,  1802,   460 }},
  {{  457,  1796,  1805 }},
  {{ 1796,    33,  1804 }},
  {{ 1805,  1796,  1804 }},
  {{ 1805,  1804,   459 }},
  {{  460,  1802,  1806 }},
  {{ 1802,   457,  1805 }},
  {{ 1806,  1802,  1805 }},
  {{ 1806,  1805,   459 }},
  {{  460,  1806,  1808 }},
  {{ 1806,   459,  1807 }},
  {{ 1808,  1806,  1807 }},
  {{ 1808,  1807,   120 }},
  {{  111,  1801,  1810 }},
  {{ 1801,   458,  1809 }},
  {{ 1810,  1801,  1809 }},
  {{ 1810,  1809,   461 }},
  {{  458,  1798,  1811 }},
  {{ 1798,   119,  1803 }},
  {{ 1811,  1798,  1803 }},
  {{ 1811,  1803,   460 }},
  {{  461,  1809,  1812 }},
  {{ 1809,   458,  1811 }},
  {{ 1812,  1809,  1811 }},
  {{ 1812,  1811,   460 }},
  {{  461,  1812,  1813 }},
  {{ 1812,   460,  1808 }},
  {{ 1813,  1812,  1808 }},
  {{ 1813,  1808,   120 }},
  {{  111,  1810,  1662 }},
  {{ 1810,   461,  1814 }},
  {{ 1662,  1810,  1814 }},
  {{ 1662,  1814,   425 }},
  {{  461,  1813,  1816 }},
  {{ 1813,   120,  1815 }},
  {{ 1816,  1813,  1815 }},
  {{ 1816,  1815,   462 }},
  {{  425,  1814,  1817 }},
  {{ 1814,   461,  1816 }},
  {{ 1817,  1814,  1816 }},
  {{ 1817,  1816,   462 }},
  {{  425,  1817,  1665 }},
  {{ 1817,   462,  1818 }},
  {{ 1665,  1817,  1818 }},
  {{ 1665,  1818,     7 }},
  {{    2,  1819,  1758 }},
  {{ 1819,   463,  1820 }},
  {{ 1758,  1819,  1820 }},
  {{ 1758,  1820,   448 }},
  {{  463,  1821,  1823 }},
  {{ 1821,   121,  1822 }},
  {{ 1823,  1821,  1822 }},
  {{ 1823,  1822,   464 }},
  {{  448,  1820,  1824 }},
  {{ 1820,   463,  1823 }},
  {{ 1824,  1820,  1823 }},
  {{ 1824,  1823,   464 }},
  {{  448,  1824,  1761 }},
  {{ 1824,   464,  1825 }},
  {{ 1761,  1824,  1825 }},
  {{ 1761,  1825,   117 }},
  {{  121,  1826,  1828 }},
  {{ 1826,   465,  1827 }},
  {{ 1828,  1826,  1827 }},
  {{ 1828,  1827,   467 }},
  {{  465,  1829,  1831 }},
  {{ 1829,    34,  1830 }},
  {{ 1831,  1829,  1830 }},
  {{ 1831,  1830,   466 }},
  {{  467,  1827,  1832 }},
  {{ 1827,   465,  1831 }},
  {{ 1832,  1827,  1831 }},
  {{ 1832,  1831,   466 }},
  {{  467,  1832,  1834 }},
  {{ 1832,   466,  1833 }},
  {{ 1834,  1832,  1833 }},
  {{ 1834,  1833,   122 }},
  {{  117,  1825,  1836 }},
  {{ 1825,   464,  1835 }},
  {{ 1836,  1825,  1835 }},
  {{ 1836,  1835,   468 }},
  {{  464,  1822,  1837 }},
  {{ 1822,   121,  1828 }},
  {{ 1837,  1822,  1828 }},
  {{ 1837,  1828,   467 }},
  {{  468,  1835,  1838 }},
  {{ 1835,   464,  1837 }},
  {{ 1838,  1835,  1837 }},
  {{ 1838,  1837,   467 }},
  {{  468,  1838,  1839 }},
  {{ 1838,   467,  1834 }},
  {{ 1839,  1838,  1834 }},
  {{ 1839,  1834,   122 }},
  {{  117,  1836,  1770 }},
  {{ 1836,   468,  1840 }},
  {{ 1770,  1836,  1840 }},
  {{ 1770,  1840,   451 }},
  {{  468,  1839,  1842 }},
  {{ 1839,   122,  1841 }},
  {{ 1842,  1839,  1841 }},
  {{ 1842,  1841,   469 }},
  {{  451,  1840,  1843 }},
  {{ 1840,   468,  1842 }},
  {{ 1843,  1840,  1842 }},
  {{ 1843,  1842,   469 }},
  {{  451,  1843,  1773 }},
  {{ 1843,   469,  1844 }},
  {{ 1773,  1843,  1844 }},
  {{ 1773,  1844,    33 }},
  {{   34,  1845,  1847 }},
  {{ 1845,   470,  1846 }},
  {{ 1847,  1845,  1846 }},
  {{ 1847,  1846,   472 }},
  {{  470,  1848,  1850 }},
  {{ 1848,   123,  1849 }},
  {{ 1850,  1848,  1849 }},
  {{ 1850,  1849,   471 }},
  {{  472,  1846,  1851 }},
  {{ 1846,   470,  1850 }},
  {{ 1851,  1846,  1850 }},
  {{ 1851,  1850,   471 }},
  {{  472,  1851,  1853 }},
  {{ 1851,   471,  1852 }},
  {{ 1853,  1851,  1852 }},
  {{ 1853,  1852,   125 }},
  {{  123,  1854,  1856 }},
  {{ 1854,   473,  1855 }},
  {{ 1856,  1854,  1855 }},
  {{ 1856,  1855,   475 }},
  {{  473,  1857,  1859 }},
  {{ 1857,     8,  1858 }},
  {{ 1859,  1857,  1858 }},
  {{ 1859,  1858,   474 }},
  {{  475,  1855,  1860 }},
  {{ 1855,   473,  1859 }},
  {{ 1860,  1855,  1859 }},
  {{ 1860,  1859,   474 }},
  {{  475,  1860,  1862 }},
  {{ 1860,   474,  1861 }},
  {{ 1862,  1860,  1861 }},
  {{ 1862,  1861,   124 }},
  {{  125,  1852,  1864 }},
  {{ 1852,   471,  1863 }},
  {{ 1864,  1852,  1863 }},
  {{ 1864,  1863,   476 }},
  {{  471,  1849,  1865 }},
  {{ 1849,   123,  1856 }},
  {{ 1865,  1849,  1856 }},
  {{ 1865,  1856,   475 }},
  {{  476,  1863,  1866 }},
  {{ 1863,   471,  1865 }},
  {{ 1866,  1863,  1865 }},
  {{ 1866,  1865,   475 }},
  {{  476,  1866,  1867 }},
  {{ 1866,   475,  1862 }},
  {{ 1867,  1866,  1862 }},
  {{ 1867,  1862,   124 }},
  {{  125,  1864,  1869 }},
  {{ 1864,   476,  1868 }},
  {{ 1869,  1864,  1868 }},
  {{ 1869,  1868,   478 }},
  {{  476,  1867,  1871 }},
  {{ 1867,   124,  1870 }},
  {{ 1871,  1867,  1870 }},
  {{ 1871,  1870,   477 }},
  {{  478,  1868,  1872 }},
  {{ 1868,   476,  1871 }},
  {{ 1872,  1868,  1871 }},
  {{ 1872,  1871,   477 }},
  {{  478,  1872,  1874 }},
  {{ 1872,   477,  1873 }},
  {{ 1874,  1872,  1873 }},
  {{ 1874,  1873,    35 }},
  {{   33,  1844,  1876 }},
  {{ 1844,   469,  1875 }},
  {{ 1876,  1844,  1875 }},
  {{ 1876,  1875,   480 }},
  {{  469,  1841,  1878 }},
  {{ 1841,   122,  1877 }},
  {{ 1878,  1841,  1877 }},
  {{ 1878,  1877,   479 }},
  {{  480,  1875,  1879 }},
  {{ 1875,   469,  1878 }},
  {{ 1879,  1875,  1878 }},
  {{ 1879,  1878,   479 }},
  {{  480,  1879,  1881 }},
  {{ 1879,   479,  1880 }},
  {{ 1881,  1879,  1880 }},
  {{ 1881,  1880,   126 }},
  {{  122,  1833,  1883 }},
  {{ 1833,   466,  1882 }},
  {{ 1883,  1833,  1882 }},
  {{ 1883,  1882,   481 }},
  {{  466,  1830,  1884 }},
  {{ 1830,    34,  1847 }},
  {{ 1884,  1830,  1847 }},
  {{ 1884,  1847,   472 }},
  {{  481,  1882,  1885 }},
  {{ 1882,   466,  1884 }},
  {{ 1885,  1882,  1884 }},
  {{ 1885,  1884,   472 }},
  {{  481,  1885,  1886 }},
  {{ 1885,   472,  1853 }},
  {{ 1886,  1885,  1853 }},
  {{ 1886,  1853,   125 }},
  {{  126,  1880,  1888 }},
  {{ 1880,   479,  1887 }},
  {{ 1888,  1880,  1887 }},
  {{ 1888,  1887,   482 }},
  {{  479,  1877,  1889 }},
  {{ 1877,   122,  1883 }},
  {{ 1889,  1877,  1883 }},
  {{ 1889,  1883,   481 }},
  {{  482,  1887,  1890 }},
  {{ 1887,   479,  1889 }},
  {{ 1890,  1887,  1889 }},
  {{ 1890,  1889,   481 }},
  {{  482,  1890,  1891 }},
  {{ 1890,   481,  1886 }},
  {{ 1891,  1890,  1886 }},
  {{ 1891,  1886,   125 }},
  {{  126,  1888,  1893 }},
  {{ 1888,   482,  1892 }},
  {{ 1893,  1888,  1892 }},
  {{ 1893,  1892,   483 }},
  {{  482,  1891,  1894 }},
  {{ 1891,   125,  1869 }},
  {{ 1894,  1891,  1869 }},
  {{ 1894,  1869,   478 }},
  {{  483,  1892,  1895 }},
  {{ 1892,   482,  1894 }},
  {{ 1895,  1892,  1894 }},
  {{ 1895,  1894,   478 }},
  {{  483,  1895,  1896 }},
  {{ 1895,   478,  1874 }},
  {{ 1896,  1895,  1874 }},
  {{ 1896,  1874,    35 }},
  {{   33,  1876,  1804 }},
  {{ 1876,   480,  1897 }},
  {{ 1804,  1876,  1897 }},
  {{ 1804,  1897,   459 }},
  {{  480,  1881,  1899 }},
  {{ 1881,   126,  1898 }},
  {{ 1899,  1881,  1898 }},
  {{ 1899,  1898,   484 }},
  {{  459,  1897,  1900 }},
  {{ 1897,   480,  1899 }},
  {{ 1900,  1897,  1899 }},
  {{ 1900,  1899,   484 }},
  {{  459,  1900,  1807 }},
  {{ 1900,   484,  1901 }},
  {{ 1807,  1900,  1901 }},
  {{ 1807,  1901,   120 }},
  {{  126,  1893,  1903 }},
  {{ 1893,   483,  1902 }},
  {{ 1903,  1893,  1902 }},
  {{ 1903,  1902,   486 }},
  {{  483,  1896,  1905 }},
  {{ 1896,    35,  1904 }},
  {{ 1905,  1896,  1904 }},
  {{ 1905,  1904,   485 }},
  {{  486,  1902,  1906 }},
  {{ 1902,   483,  1905 }},
  {{ 1906,  1902,  1905 }},
  {{ 1906,  1905,   485 }},
  {{  486,  1906,  1908 }},
  {{ 1906,   485,  1907 }},
  {{ 1908,  1906,  1907 }},
  {{ 1908,  1907,   127 }},
  {{  120,  1901,  1910 }},
  {{ 1901,   484,  1909 }},
  {{ 1910,  1901,  1909 }},
  {{ 1910,  1909,   487 }},
  {{  484,  1898,  1911 }},
  {{ 1898,   126,  1903 }},
  {{ 1911,  1898,  1903 }},
  {{ 1911,  1903,   486 }},
  {{  487,  1909,  1912 }},
  {{ 1909,   484,  1911 }},
  {{ 1912,  1909,  1911 }},
  {{ 1912,  1911,   486 }},
  {{  487,  1912,  1913 }},
  {{ 1912,   486,  1908 }},
  {{ 1913,  1912,  1908 }},
  {{ 1913,  1908,   127 }},
  {{  120,  1910,  1815 }},
  {{ 1910,   487,  1914 }},
  {{ 1815,  1910,  1914 }},
  {{ 1815,  1914,   462 }},
  {{  487,  1913,  1916 }},
  {{ 1913,   127,  1915 }},
  {{ 1916,  1913,  1915 }},
  {{ 1916,  1915,   488 }},
  {{  462,  1914,  1917 }},
  {{ 1914,   487,  1916 }},
  {{ 1917,  1914,  1916 }},
  {{ 1917,  1916,   488 }},
  {{  462,  1917,  1818 }},
  {{ 1917,   488,  1918 }},
  {{ 1818,  1917,  1918 }},
  {{ 1818,  1918,     7 }},
  {{    2,  1050,  1819 }},
  {{ 1050,   270,  1919 }},
  {{ 1819,  1050,  1919 }},
  {{ 1819,  1919,   463 }},
  {{  270,  1047,  1921 }},
  {{ 1047,    72,  1920 }},
  {{ 1921,  1047,  1920 }},
  {{ 1921,  1920,   489 }},
  {{  463,  1919,  1922 }},
  {{ 1919,   270,  1921 }},
  {{ 1922,  1919,  1921 }},
  {{ 1922,  1921,   489 }},
  {{  463,  1922,  1821 }},
  {{ 1922,   489,  1923 }},
  {{ 1821,  1922,  1923 }},
  {{ 1821,  1923,   121 }},
  {{   72,  1039,  1925 }},
  {{ 1039,   267,  1924 }},
  {{ 1925,  1039,  1924 }},
  {{ 1925,  1924,   491 }},
  {{  267,  1036,  1927 }},
  {{ 1036,    21,  1926 }},
  {{ 1927,  1036,  1926 }},
  {{ 1927,  1926,   490 }},
  {{  491,  1924,  1928 }},
  {{ 1924,   267,  1927 }},
  {{ 1928,  1924,  1927 }},
  {{ 1928,  1927,   490 }},
  {{  491,  1928,  1930 }},
  {{ 1928,   490,  1929 }},
  {{ 1930,  1928,  1929 }},
  {{ 1930,  1929,   128 }},
  {{  121,  1923,  1932 }},
  {{ 1923,   489,  1931 }},
  {{ 1932,  1923,  1931 }},
  {{ 1932,  1931,   492 }},
  {{  489,  1920,  1933 }},
  {{ 1920,    72,  1925 }},
  {{ 1933,  1920,  1925 }},
  {{ 1933,  1925,   491 }},
  {{  492,  1931,  1934 }},
  {{ 1931,   489,  1933 }},
  {{ 1934,  1931,  1933 }},
  {{ 1934,  1933,   491 }},
  {{  492,  1934,  1935 }},
  {{ 1934,   491,  1930 }},
  {{ 1935,  1934,  1930 }},
  {{ 1935,  1930,   128 }},
  {{  121,  1932,  1826 }},
  {{ 1932,   492,  1936 }},
  {{ 1826,  1932,  1936 }},
  {{ 1826,  1936,   465 }},
  {{  492,  1935,  1938 }},
  {{ 1935,   128,  1937 }},
  {{ 1938,  1935,  1937 }},
  {{ 1938,  1937,   493 }},
  {{  465,  1936,  1939 }},
  {{ 1936,   492,  1938 }},
  {{ 1939,  1936,  1938 }},
  {{ 1939,  1938,   493 }},
  {{  465,  1939,  1829 }},
  {{ 1939,   493,  1940 }},
  {{ 1829,  1939,  1940 }},
  {{ 1829,  1940,    34 }},
  {{   21,  1005,  1942 }},
  {{ 1005,   259,  1941 }},
  {{ 1942,  1005,  1941 }},
  {{ 1942,  1941,   495 }},
  {{  259,  1002,  1944 }},
  {{ 1002,    69,  1943 }},
  {{ 1944,  1002,  1943 }},
  {{ 1944,  1943,   494 }},
  {{  495,  1941,  1945 }},
  {{ 1941,   259,  1944 }},
  {{ 1945,  1941,  1944 }},
  {{ 1945,  1944,   494 }},
  {{  495,  1945,  1947 }},
  {{ 1945,   494,  1946 }},
  {{ 1947,  1945,  1946 }},
  {{ 1947,  1946,   130 }},
  {{   69,   993,  1949 }},
  {{  993,   256,  1948 }},
  {{ 1949,   993,  1948 }},
  {{ 1949,  1948,   497 }},
  {{  256,   990,  1951 }},
  {{  990,     3,  1950 }},
  {{ 1951,   990,  1950 }},
  {{ 1951,  1950,   496 }},
  {{  497,  1948,  1952 }},
  {{ 1948,   256,  1951 }},
  {{ 1952,  1948,  1951 }},
  {{ 1952,  1951,   496 }},
  {{  497,  1952,  1954 }},
  {{ 1952,   496,  1953 }},
  {{ 1954,  1952,  1953 }},
  {{ 1954,  1953,   129 }},
  {{  130,  1946,  1956 }},
  {{ 1946,   494,  1955 }},
  {{ 1956,  1946,  1955 }},
  {{ 1956,  1955,   498 }},
  {{  494,  1943,  1957 }},
  {{ 1943,    69,  1949 }},
  {{ 1957,  1943,  1949 }},
  {{ 1957,  1949,   497 }},
  {{  498,  1955,  1958 }},
  {{ 1955,   494,  1957 }},
  {{ 1958,  1955,  1957 }},
  {{ 1958,  1957,   497 }},
  {{  498,  1958,  1959 }},
  {{ 1958,   497,  1954 }},
  {{ 1959,  1958,  1954 }},
  {{ 1959,  1954,   129 }},
  {{  130,  1956,  1961 }},
  {{ 1956,   498,  1960 }},
  {{ 1961,  1956,  1960 }},
  {{ 1961,  1960,   500 }},
  {{  498,  1959,  1963 }},
  {{ 1959,   129,  1962 }},
  {{ 1963,  1959,  1962 }},
  {{ 1963,  1962,   499 }},
  {{  500,  1960,  1964 }},
  {{ 1960,   498,  1963 }},
  {{ 1964,  1960,  1963 }},
  {{ 1964,  1963,   499 }},
  {{  500,  1964,  1966 }},
  {{ 1964,   499,  1965 }},
  {{ 1966,  1964,  1965 }},
  {{ 1966,  1965,    36 }},
  {{   34,  1940,  1968 }},
  {{ 1940,   493,  1967 }},
  {{ 1968,  1940,  1967 }},
  {{ 1968,  1967,   502 }},
  {{  493,  1937,  1970 }},
  {{ 1937,   128,  1969 }},
  {{ 1970,  1937,  1969 }},
  {{ 1970,  1969,   501 }},
  {{  502,  1967,  1971 }},
  {{ 1967,   493,  1970 }},
  {{ 1971,  1967,  1970 }},
  {{ 1971,  1970,   501 }},
  {{  502,  1971,  1973 }},
  {{ 1971,   501,  1972 }},
  {{ 1973,  1971,  1972 }},
  {{ 1973,  1972,   131 }},
  {{  128,  1929,  1975 }},
  {{ 1929,   490,  1974 }},
  {{ 1975,  1929,  1974 }},
  {{ 1975,  1974,   503 }},
  {{  490,  1926,  1976 }},
  {{ 1926,    21,  1942 }},
  {{ 1976,  1926,  1942 }},
  {{ 1976,  1942,   495 }},
  {{  503,  1974,  1977 }},
  {{ 1974,   490,  1976 }},
  {{ 1977,  1974,  1976 }},
  {{ 1977,  1976,   495 }},
  {{  503,  1977,  1978 }},
  {{ 1977,   495,  1947 }},
  {{ 1978,  1977,  1947 }},
  {{ 1978,  1947,   130 }},
  {{  131,  1972,  1980 }},
  {{ 1972,   501,  1979 }},
  {{ 1980,  1972,  1979 }},
  {{ 1980,  1979,   504 }},
  {{  501,  1969,  1981 }},
  {{ 1969,   128,  1975 }},
  {{ 1981,  1969,  1975 }},
  {{ 1981,  1975,   503 }},
  {{  504,  1979,  1982 }},
  {{ 1979,   501,  1981 }},
  {{ 1982,  1979,  1981 }},
  {{ 1982,  1981,   503 }},
  {{  504,  1982,  1983 }},
  {{ 1982,   503,  1978 }},
  {{ 1983,  1982,  1978 }},
  {{ 1983,  1978,   130 }},
  {{  131,  1980,  1985 }},
  {{ 1980,   504,  1984 }},
  {{ 1985,  1980,  1984 }},
  {{ 1985,  1984,   505 }},
  {{  504,  1983,  1986 }},
  {{ 1983,   130,  1961 }},
  {{ 1986,  1983,  1961 }},
  {{ 1986,  1961,   500 }},
  {{  505,  1984,  1987 }},
  {{ 1984,   504,  1986 }},
  {{ 1987,  1984,  1986 }},
  {{ 1987,  1986,   500 }},
  {{  505,  1987,  1988 }},
  {{ 1987,   500,  1966 }},
  {{ 1988,  1987,  1966 }},
  {{ 1988,  1966,    36 }},
  {{   34,  1968,  1845 }},
  {{ 1968,   502,  1989 }},
  {{ 1845,  1968,  1989 }},
  {{ 1845,  1989,   470 }},
  {{  502,  1973,  1991 }},
  {{ 1973,   131,  1990 }},
  {{ 1991,  1973,  1990 }},
  {{ 1991,  1990,   506 }},
  {{  470,  1989,  1992 }},
  {{ 1989,   502,  1991 }},
  {{ 1992,  1989,  1991 }},
  {{ 1992,  1991,   506 }},
  {{  470,  1992,  1848 }},
  {{ 1992,   506,  1993 }},
  {{ 1848,  1992,  1993 }},
  {{ 1848,  1993,   123 }},
  {{  131,  1985,  1995 }},
  {{ 1985,   505,  1994 }},
  {{ 1995,  1985,  1994 }},
  {{ 1995,  1994,   508 }},
  {{  505,  1988,  1997 }},
  {{ 1988,    36,  1996 }},
  {{ 1997,  1988,  1996 }},
  {{ 1997,  1996,   507 }},
  {{  508,  1994,  1998 }},
  {{ 1994,   505,  1997 }},
  {{ 1998,  1994,  1997 }},
  {{ 1998,  1997,   507 }},
  {{  508,  1998,  2000 }},
  {{ 1998,   507,  1999 }},
  {{ 2000,  1998,  1999 }},
  {{ 2000,  1999,   132 }},
  {{  123,  1993,  2002 }},
  {{ 1993,   506,  2001 }},
  {{ 2002,  1993,  2001 }},
  {{ 2002,  2001,   509 }},
  {{  506,  1990,  2003 }},
  {{ 1990,   131,  1995 }},
  {{ 2003,  1990,  1995 }},
  {{ 2003,  1995,   508 }},
  {{  509,  2001,  2004 }},
  {{ 2001,   506,  2003 }},
  {{ 2004,  2001,  2003 }},
  {{ 2004,  2003,   508 }},
  {{  509,  2004,  2005 }},
  {{ 2004,   508,  2000 }},
  {{ 2005,  2004,  2000 }},
  {{ 2005,  2000,   132 }},
  {{  123,  2002,  1854 }},
  {{ 2002,   509,  2006 }},
  {{ 1854,  2002,  2006 }},
  {{ 1854,  2006,   473 }},
  {{  509,  2005,  2008 }},
  {{ 2005,   132,  2007 }},
  {{ 2008,  2005,  2007 }},
  {{ 2008,  2007,   510 }},
  {{  473,  2006,  2009 }},
  {{ 2006,   509,  2008 }},
  {{ 2009,  2006,  2008 }},
  {{ 2009,  2008,   510 }},
  {{  473,  2009,  1857 }},
  {{ 2009,   510,  2010 }},
  {{ 1857,  2009,  2010 }},
  {{ 1857,  2010,     8 }},
  {{    3,  1242,  1950 }},
  {{ 1242,   318,  2011 }},
  {{ 1950,  1242,  2011 }},
  {{ 1950,  2011,   496 }},
  {{  318,  1239,  2013 }},
  {{ 1239,    84,  2012 }},
  {{ 2013,  1239,  2012 }},
  {{ 2013,  2012,   511 }},
  {{  496,  2011,  2014 }},
  {{ 2011,   318,  2013 }},
  {{ 2014,  2011,  2013 }},
  {{ 2014,  2013,   511 }},
  {{  496,  2014,  1953 }},
  {{ 2014,   511,  2015 }},
  {{ 1953,  2014,  2015 }},
  {{ 1953,  2015,   129 }},
  {{   84,  1231,  2017 }},
  {{ 1231,   315,  2016 }},
  {{ 2017,  1231,  2016 }},
  {{ 2017,  2016,   513 }},
  {{  315,  1228,  2019 }},
  {{ 1228,    24,  2018 }},
  {{ 2019,  1228,  2018 }},
  {{ 2019,  2018,   512 }},
  {{  513,  2016,  2020 }},
  {{ 2016,   315,  2019 }},
  {{ 2020,  2016,  2019 }},
  {{ 2020,  2019,   512 }},
  {{  513,  2020,  2022 }},
  {{ 2020,   512,  2021 }},
  {{ 2022,  2020,  2021 }},
  {{ 2022,  2021,   133 }},
  {{  129,  2015,  2024 }},
  {{ 2015,   511,  2023 }},
  {{ 2024,  2015,  2023 }},
  {{ 2024,  2023,   514 }},
  {{  511,  2012,  2025 }},
  {{ 2012,    84,  2017 }},
  {{ 2025,  2012,  2017 }},
  {{ 2025,  2017,   513 }},
  {{  514,  2023,  2026 }},
  {{ 2023,   511,  2025 }},
  {{ 2026,  2023,  2025 }},
  {{ 2026,  2025,   513 }},
  {{  514,  2026,  2027 }},
  {{ 2026,   513,  2022 }},
  {{ 2027,  2026,  2022 }},
  {{ 2027,  2022,   133 }},
  {{  129,  2024,  1962 }},
  {{ 2024,   514,  2028 }},
  {{ 1962,  2024,  2028 }},
  {{ 1962,  2028,   499 }},
  {{  514,  2027,  2030 }},
  {{ 2027,   133,  2029 }},
  {{ 2030,  2027,  2029 }},
  {{ 2030,  2029,   515 }},
  {{  499,  2028,  2031 }},
  {{ 2028,   514,  2030 }},
  {{ 2031,  2028,  2030 }},
  {{ 2031,  2030,   515 }},
  {{  499,  2031,  1965 }},
  {{ 2031,   515,  2032 }},
  {{ 1965,  2031,  2032 }},
  {{ 1965,  2032,    36 }},
  {{   24,  1197,  2034 }},
  {{ 1197,   307,  2033 }},
  {{ 2034,  1197,  2033 }},
  {{ 2034,  2033,   517 }},
  {{  307,  1194,  2036 }},
  {{ 1194,    81,  2035 }},
  {{ 2036,  1194,  2035 }},
  {{ 2036,  2035,   516 }},
  {{  517,  2033,  2037 }},
  {{ 2033,   307,  2036 }},
  {{ 2037,  2033,  2036 }},
  {{ 2037,  2036,   516 }},
  {{  517,  2037,  2039 }},
  {{ 2037,   516,  2038 }},
  {{ 2039,  2037,  2038 }},
  {{ 2039,  2038,   135 }},
  {{   81,  1185,  2041 }},
  {{ 1185,   304,  2040 }},
  {{ 2041,  1185,  2040 }},
  {{ 2041,  2040,   519 }},
  {{  304,  1182,  2043 }},
  {{ 1182,     9,  2042 }},
  {{ 2043,  1182,  2042 }},
  {{ 2043,  2042,   518 }},
  {{  519,  2040,  2044 }},
  {{ 2040,   304,  2043 }},
  {{ 2044,  2040,  2043 }},
  {{ 2044,  2043,   518 }},
  {{  519,  2044,  2046 }},
  {{ 2044,   518,  2045 }},
  {{ 2046,  2044,  2045 }},
  {{ 2046,  2045,   134 }},
  {{  135,  2038,  2048 }},
  {{ 2038,   516,  2047 }},
  {{ 2048,  2038,  2047 }},
  {{ 2048,  2047,   520 }},
  {{  516,  2035,  2049 }},
  {{ 2035,    81,  2041 }},
  {{ 2049,  2035,  2041 }},
  {{ 2049,  2041,   519 }},
  {{  520,  2047,  2050 }},
  {{ 2047,   516,  2049 }},
  {{ 2050,  2047,  2049 }},
  {{ 2050,  2049,   519 }},
  {{  520,  2050,  2051 }},
  {{ 2050,   519,  2046 }},
  {{ 2051,  2050,  2046 }},
  {{ 2051,  2046,   134 }},
  {{  135,  2048,  2053 }},
  {{ 2048,   520,  2052 }},
  {{ 2053,  2048,  2052 }},
  {{ 2053,  2052,   522 }},
  {{  520,  2051,  2055 }},
  {{ 2051,   134,  2054 }},
  {{ 2055,  2051,  2054 }},
  {{ 2055,  2054,   521 }},
  {{  522,  2052,  2056 }},
  {{ 2052,   520,  2055 }},
  {{ 2056,  2052,  2055 }},
  {{ 2056,  2055,   521 }},
  {{  522,  2056,  2058 }},
  {{ 2056,   521,  2057 }},
  {{ 2058,  2056,  2057 }},
  {{ 2058,  2057,    37 }},
  {{   36,  2032,  2060 }},
  {{ 2032,   515,  2059 }},
  {{ 2060,  2032,  2059 }},
  {{ 2060,  2059,   524 }},
  {{  515,  2029,  2062 }},
  {{ 2029,   133,  2061 }},
  {{ 2062,  2029,  2061 }},
  {{ 2062,  2061,   523 }},
  {{  524,  2059,  2063 }},
  {{ 2059,   515,  2062 }},
  {{ 2063,  2059,  2062 }},
  {{ 2063,  2062,   523 }},
  {{  524,  2063,  2065 }},
  {{ 2063,   523,  2064 }},
  {{ 2065,  2063,  2064 }},
  {{ 2065,  2064,   136 }},
  {{  133,  2021,  2067 }},
  {{ 2021,   512,  2066 }},
  {{ 2067,  2021,  2066 }},
  {{ 2067,  2066,   525 }},
  {{  512,  2018,  2068 }},
  {{ 2018,    24,  2034 }},
  {{ 2068,  2018,  2034 }},
  {{ 2068,  2034,   517 }},
  {{  525,  2066,  2069 }},
  {{ 2066,   512,  2068 }},
  {{ 2069,  2066,  2068 }},
  {{ 2069,  2068,   517 }},
  {{  525,  2069,  2070 }},
  {{ 2069,   517,  2039 }},
  {{ 2070,  2069,  2039 }},
  {{ 2070,  2039,   135 }},
  {{  136,  2064,  2072 }},
  {{ 2064,   523,  2071 }},
  {{ 2072,  2064,  2071 }},
  {{ 2072,  2071,   526 }},
  {{  523,  2061,  2073 }},
  {{ 2061,   133,  2067 }},
  {{ 2073,  2061,  2067 }},
  {{ 2073,  2067,   525 }},
  {{  526,  2071,  2074 }},
  {{ 2071,   523,  2073 }},
  {{ 2074,  2071,  2073 }},
  {{ 2074,  2073,   525 }},
  {{  526,  2074,  2075 }},
  {{ 2074,   525,  2070 }},
  {{ 2075,  2074,  2070 }},
  {{ 2075,  2070,   135 }},
  {{  136,  2072,  2077 }},
  {{ 2072,   526,  2076 }},
  {{ 2077,  2072,  2076 }},
  {{ 2077,  2076,   527 }},
  {{  526,  2075,  2078 }},
  {{ 2075,   135,  2053 }},
  {{ 2078,  2075,  2053 }},
  {{ 2078,  2053,   522 }},
  {{  527,  2076,  2079 }},
  {{ 2076,   526,  2078 }},
  {{ 2079,  2076,  2078 }},
  {{ 2079,  2078,   522 }},
  {{  527,  2079,  2080 }},
  {{ 2079,   522,  2058 }},
  {{ 2080,  2079,  2058 }},
  {{ 2080,  2058,    37 }},
  {{   36,  2060,  1996 }},
  {{ 2060,   524,  2081 }},
  {{ 1996,  2060,  2081 }},
  {{ 1996,  2081,   507 }},
  {{  524,  2065,  2083 }},
  {{ 2065,   136,  2082 }},
  {{ 2083,  2065,  2082 }},
  {{ 2083,  2082,   528 }},
  {{  507,  2081,  2084 }},
  {{ 2081,   524,  2083 }},
  {{ 2084,  2081,  2083 }},
  {{ 2084,  2083,   528 }},
  {{  507,  2084,  1999 }},
  {{ 2084,   528,  2085 }},
  {{ 1999,  2084,  2085 }},
  {{ 1999,  2085,   132 }},
  {{  136,  2077,  2087 }},
  {{ 2077,   527,  2086 }},
  {{ 2087,  2077,  2086 }},
  {{ 2087,  2086,   530 }},
  {{  527,  2080,  2089 }},
  {{ 2080,    37,  2088 }},
  {{ 2089,  2080,  2088 }},
  {{ 2089,  2088,   529 }},
  {{  530,  2086,  2090 }},
  {{ 2086,   527,  2089 }},
  {{ 2090,  2086,  2089 }},
  {{ 2090,  2089,   529 }},
  {{  530,  2090,  2092 }},
  {{ 2090,   529,  2091 }},
  {{ 2092,  2090,  2091 }},
  {{ 2092,  2091,   137 }},
  {{  132,  2085,  2094 }},
  {{ 2085,   528,  2093 }},
  {{ 2094,  2085,  2093 }},
  {{ 2094,  2093,   531 }},
  {{  528,  2082,  2095 }},
  {{ 2082,   136,  2087 }},
  {{ 2095,  2082,  2087 }},
  {{ 2095,  2087,   530 }},
  {{  531,  2093,  2096 }},
  {{ 2093,   528,  2095 }},
  {{ 2096,  2093,  2095 }},
  {{ 2096,  2095,   530 }},
  {{  531,  2096,  2097 }},
  {{ 2096,   530,  2092 }},
  {{ 2097,  2096,  2092 }},
  {{ 2097,  2092,   137 }},
  {{  132,  2094,  2007 }},
  {{ 2094,   531,  2098 }},
  {{ 2007,  2094,  2098 }},
  {{ 2007,  2098,   510 }},
  {{  531,  2097,  2100 }},
  {{ 2097,   137,  2099 }},
  {{ 2100,  2097,  2099 }},
  {{ 2100,  2099,   532 }},
  {{  510,  2098,  2101 }},
  {{ 2098,   531,  2100 }},
  {{ 2101,  2098,  2100 }},
  {{ 2101,  2100,   532 }},
  {{  510,  2101,  2010 }},
  {{ 2101,   532,  2102 }},
  {{ 2010,  2101,  2102 }},
  {{ 2010,  2102,     8 }},
  {{    9,  1342,  2104 }},
  {{ 1342,   344,  2103 }},
  {{ 2104,  1342,  2103 }},
  {{ 2104,  2103,   534 }},
  {{  344,  1339,  2106 }},
  {{ 1339,    91,  2105 }},
  {{ 2106,  1339,  2105 }},
  {{ 2106,  2105,   533 }},
  {{  534,  2103,  2107 }},
  {{ 2103,   344,  2106 }},
  {{ 2107,  2103,  2106 }},
  {{ 2107,  2106,   533 }},
  {{  534,  2107,  2109 }},
  {{ 2107,   533,  2108 }},
  {{ 2109,  2107,  2108 }},
  {{ 2109,  2108,   139 }},
  {{   91,  1331,  2111 }},
  {{ 1331,   341,  2110 }},
  {{ 2111,  1331,  2110 }},
  {{ 2111,  2110,   536 }},
  {{  341,  1328,  2113 }},
  {{ 1328,    26,  2112 }},
  {{ 2113,  1328,  2112 }},
  {{ 2113,  2112,   535 }},
  {{  536,  2110,  2114 }},
  {{ 2110,   341,  2113 }},
  {{ 2114,  2110,  2113 }},
  {{ 2114,  2113,   535 }},
  {{  536,  2114,  2116 }},
  {{ 2114,   535,  2115 }},
  {{ 2116,  2114,  2115 }},
  {{ 2116,  2115,   138 }},
  {{  139,  2108,  2118 }},
  {{ 2108,   533,  2117 }},
  {{ 2118,  2108,  2117 }},
  {{ 2118,  2117,   537 }},
  {{  533,  2105,  2119 }},
  {{ 2105,    91,  2111 }},
  {{ 2119,  2105,  2111 }},
  {{ 2119,  2111,   536 }},
  {{  537,  2117,  2120 }},
  {{ 2117,   533,  2119 }},
  {{ 2120,  2117,  2119 }},
  {{ 2120,  2119,   536 }},
  {{  537,  2120,  2121 }},
  {{ 2120,   536,  2116 }},
  {{ 2121,  2120,  2116 }},
  {{ 2121,  2116,   138 }},
  {{  139,  2118,  2123 }},
  {{ 2118,   537,  2122 }},
  {{ 2123,  2118,  2122 }},
  {{ 2123,  2122,   539 }},
  {{  537,  2121,  2125 }},
  {{ 2121,   138,  2124 }},
  {{ 2125,  2121,  2124 }},
  {{ 2125,  2124,   538 }},
  {{  539,  2122,  2126 }},
  {{ 2122,   537,  2125 }},
  {{ 2126,  2122,  2125 }},
  {{ 2126,  2125,   538 }},
  {{  539,  2126,  2128 }},
  {{ 2126,   538,  2127 }},
  {{ 2128,  2126,  2127 }},
  {{ 2128,  2127,    39 }},
  {{   26,  1297,  2130 }},
  {{ 1297,   333,  2129 }},
  {{ 2130,  1297,  2129 }},
  {{ 2130,  2129,   541 }},
  {{  333,  1294,  2132 }},
  {{ 1294,    88,  2131 }},
  {{ 2132,  1294,  2131 }},
  {{ 2132,  2131,   540 }},
  {{  541,  2129,  2133 }},
  {{ 2129,   333,  2132 }},
  {{ 2133,  2129,  2132 }},
  {{ 2133,  2132,   540 }},
  {{  541,  2133,  2135 }},
  {{ 2133,   540,  2134 }},
  {{ 2135,  2133,  2134 }},
  {{ 2135,  2134,   141 }},
  {{   88,  1285,  2137 }},
  {{ 1285,   330,  2136 }},
  {{ 2137,  1285,  2136 }},
  {{ 2137,  2136,   543 }},
  {{  330,  1282,  2139 }},
  {{ 1282,    10,  2138 }},
  {{ 2139,  1282,  2138 }},
  {{ 2139,  2138,   542 }},
  {{  543,  2136,  2140 }},
  {{ 2136,   330,  2139 }},
  {{ 2140,  2136,  2139 }},
  {{ 2140,  2139,   542 }},
  {{  543,  2140,  2142 }},
  {{ 2140,   542,  2141 }},
  {{ 2142,  2140,  2141 }},
  {{ 2142,  2141,   140 }},
  {{  141,  2134,  2144 }},
  {{ 2134,   540,  2143 }},
  {{ 2144,  2134,  2143 }},
  {{ 2144,  2143,   544 }},
  {{  540,  2131,  2145 }},
  {{ 2131,    88,  2137 }},
  {{ 2145,  2131,  2137 }},
  {{ 2145,  2137,   543 }},
  {{  544,  2143,  2146 }},
  {{ 2143,   540,  2145 }},
  {{ 2146,  2143,  2145 }},
  {{ 2146,  2145,   543 }},
  {{  544,  2146,  2147 }},
  {{ 2146,   543,  2142 }},
  {{ 2147,  2146,  2142 }},
  {{ 2147,  2142,   140 }},
  {{  141,  2144,  2149 }},
  {{ 2144,   544,  2148 }},
  {{ 2149,  2144,  2148 }},
  {{ 2149,  2148,   546 }},
  {{  544,  2147,  2151 }},
  {{ 2147,   140,  2150 }},
  {{ 2151,  2147,  2150 }},
  {{ 2151,  2150,   545 }},
  {{  546,  2148,  2152 }},
  {{ 2148,   544,  2151 }},
  {{ 2152,  2148,  2151 }},
  {{ 2152,  2151,   545 }},
  {{  546,  2152,  2154 }},
  {{ 2152,   545,  2153 }},
  {{ 2154,  2152,  2153 }},
  {{ 2154,  2153,    38 }},
  {{   39,  2127,  2156 }},
  {{ 2127,   538,  2155 }},
  {{ 2156,  2127,  2155 }},
  {{ 2156,  2155,   548 }},
  {{  538,  2124,  2158 }},
  {{ 2124,   138,  2157 }},
  {{ 2158,  2124,  2157 }},
  {{ 2158,  2157,   547 }},
  {{  548,  2155,  2159 }},
  {{ 2155,   538,  2158 }},
  {{ 2159,  2155,  2158 }},
  {{ 2159,  2158,   547 }},
  {{  548,  2159,  2161 }},
  {{ 2159,   547,  2160 }},
  {{ 2161,  2159,  2160 }},
  {{ 2161,  2160,   142 }},
  {{  138,  2115,  2163 }},
  {{ 2115,   535,  2162 }},
  {{ 2163,  2115,  2162 }},
  {{ 2163,  2162,   549 }},
  {{  535,  2112,  2164 }},
  {{ 2112,    26,  2130 }},
  {{ 2164,  2112,  2130 }},
  {{ 2164,  2130,   541 }},
  {{  549,  2162,  2165 }},
  {{ 2162,   535,  2164 }},
  {{ 2165,  2162,  2164 }},
  {{ 2165,  2164,   541 }},
  {{  549,  2165,  2166 }},
  {{ 2165,   541,  2135 }},
  {{ 2166,  2165,  2135 }},
  {{ 2166,  2135,   141 }},
  {{  142,  2160,  2168 }},
  {{ 2160,   547,  2167 }},
  {{ 2168,  2160,  2167 }},
  {{ 2168,  2167,   550 }},
  {{  547,  2157,  2169 }},
  {{ 2157,   138,  2163 }},
  {{ 2169,  2157,  2163 }},
  {{ 2169,  2163,   549 }},
  {{  550,  2167,  2170 }},
  {{ 2167,   547,  2169 }},
  {{ 2170,  2167,  2169 }},
  {{ 2170,  2169,   549 }},
  {{  550,  2170,  2171 }},
  {{ 2170,   549,  2166 }},
  {{ 2171,  2170,  2166 }},
  {{ 2171,  2166,   141 }},
  {{  142,  2168,  2173 }},
  {{ 2168,   550,  2172 }},
  {{ 2173,  2168,  2172 }},
  {{ 2173,  2172,   551 }},
  {{  550,  2171,  2174 }},
  {{ 2171,   141,  2149 }},
  {{ 2174,  2171,  2149 }},
  {{ 2174,  2149,   546 }},
  {{  551,  2172,  2175 }},
  {{ 2172,   550,  2174 }},
  {{ 2175,  2172,  2174 }},
  {{ 2175,  2174,   546 }},
  {{  551,  2175,  2176 }},
  {{ 2175,   546,  2154 }},
  {{ 2176,  2175,  2154 }},
  {{ 2176,  2154,    38 }},
  {{   39,  2156,  2178 }},
  {{ 2156,   548,  2177 }},
  {{ 2178,  2156,  2177 }},
  {{ 2178,  2177,   553 }},
  {{  548,  2161,  2180 }},
  {{ 2161,   142,  2179 }},
  {{ 2180,  2161,  2179 }},
  {{ 2180,  2179,   552 }},
  {{  553,  2177,  2181 }},
  {{ 2177,   548,  2180 }},
  {{ 2181,  2177,  2180 }},
  {{ 2181,  2180,   552 }},
  {{  553,  2181,  2183 }},
  {{ 2181,   552,  2182 }},
  {{ 2183,  2181,  2182 }},
  {{ 2183,  2182,   144 }},
  {{  142,  2173,  2185 }},
  {{ 2173,   551,  2184 }},
  {{ 2185,  2173,  2184 }},
  {{ 2185,  2184,   555 }},
  {{  551,  2176,  2187 }},
  {{ 2176,    38,  2186 }},
  {{ 2187,  2176,  2186 }},
  {{ 2187,  2186,   554 }},
  {{  555,  2184,  2188 }},
  {{ 2184,   551,  2187 }},
  {{ 2188,  2184,  2187 }},
  {{ 2188,  2187,   554 }},
  {{  555,  2188,  2190 }},
  {{ 2188,   554,  2189 }},
  {{ 2190,  2188,  2189 }},
  {{ 2190,  2189,   143 }},
  {{  144,  2182,  2192 }},
  {{ 2182,   552,  2191 }},
  {{ 2192,  2182,  2191 }},
  {{ 2192,  2191,   556 }},
  {{  552,  2179,  2193 }},
  {{ 2179,   142,  2185 }},
  {{ 2193,  2179,  2185 }},
  {{ 2193,  2185,   555 }},
  {{  556,  2191,  2194 }},
  {{ 2191,   552,  2193 }},
  {{ 2194,  2191,  2193 }},
  {{ 2194,  2193,   555 }},
  {{  556,  2194,  2195 }},
  {{ 2194,   555,  2190 }},
  {{ 2195,  2194,  2190 }},
  {{ 2195,  2190,   143 }},
  {{  144,  2192,  2197 }},
  {{ 2192,   556,  2196 }},
  {{ 2197,  2192,  2196 }},
  {{ 2197,  2196,   558 }},
  {{  556,  2195,  2199 }},
  {{ 2195,   143,  2198 }},
  {{ 2199,  2195,  2198 }},
  {{ 2199,  2198,   557 }},
  {{  558,  2196,  2200 }},
  {{ 2196,   556,  2199 }},
  {{ 2200,  2196,  2199 }},
  {{ 2200,  2199,   557 }},
  {{  558,  2200,  2202 }},
  {{ 2200,   557,  2201 }},
  {{ 2202,  2200,  2201 }},
  {{ 2202,  2201,    12 }},
  {{   10,  1534,  2138 }},
  {{ 1534,   392,  2203 }},
  {{ 2138,  1534,  2203 }},
  {{ 2138,  2203,   542 }},
  {{  392,  1531,  2205 }},
  {{ 1531,   103,  2204 }},
  {{ 2205,  1531,  2204 }},
  {{ 2205,  2204,   559 }},
  {{  542,  2203,  2206 }},
  {{ 2203,   392,  2205 }},
  {{ 2206,  2203,  2205 }},
  {{ 2206,  2205,   559 }},
  {{  542,  2206,  2141 }},
  {{ 2206,   559,  2207 }},
  {{ 2141,  2206,  2207 }},
  {{ 2141,  2207,   140 }},
  {{  103,  1523,  2209 }},
  {{ 1523,   389,  2208 }},
  {{ 2209,  1523,  2208 }},
  {{ 2209,  2208,   561 }},
  {{  389,  1520,  2211 }},
  {{ 1520,    29,  2210 }},
  {{ 2211,  1520,  2210 }},
  {{ 2211,  2210,   560 }},
  {{  561,  2208,  2212 }},
  {{ 2208,   389,  2211 }},
  {{ 2212,  2208,  2211 }},
  {{ 2212,  2211,   560 }},
  {{  561,  2212,  2214 }},
  {{ 2212,   560,  2213 }},
  {{ 2214,  2212,  2213 }},
  {{ 2214,  2213,   145 }},
  {{  140,  2207,  2216 }},
  {{ 2207,   559,  2215 }},
  {{ 2216,  2207,  2215 }},
  {{ 2216,  2215,   562 }},
  {{  559,  2204,  2217 }},
  {{ 2204,   103,  2209 }},
  {{ 2217,  2204,  2209 }},
  {{ 2217,  2209,   561 }},
  {{  562,  2215,  2218 }},
  {{ 2215,   559,  2217 }},
  {{ 2218,  2215,  2217 }},
  {{ 2218,  2217,   561 }},
  {{  562,  2218,  2219 }},
  {{ 2218,   561,  2214 }},
  {{ 2219,  2218,  2214 }},
  {{ 2219,  2214,   145 }},
  {{  140,  2216,  2150 }},
  {{ 2216,   562,  2220 }},
  {{ 2150,  2216,  2220 }},
  {{ 2150,  2220,   545 }},
  {{  562,  2219,  2222 }},
  {{ 2219,   145,  2221 }},
  {{ 2222,  2219,  2221 }},
  {{ 2222,  2221,   563 }},
  {{  545,  2220,  2223 }},
  {{ 2220,   562,  2222 }},
  {{ 2223,  2220,  2222 }},
  {{ 2223,  2222,   563 }},
  {{  545,  2223,  2153 }},
  {{ 2223,   563,  2224 }},
  {{ 2153,  2223,  2224 }},
  {{ 2153,  2224,    38 }},
  {{   29,  1489,  2226 }},
  {{ 1489,   381,  2225 }},
  {{ 2226,  1489,  2225 }},
  {{ 2226,  2225,   565 }},
  {{  381,  1486,  2228 }},
  {{ 1486,   100,  2227 }},
  {{ 2228,  1486,  2227 }},
  {{ 2228,  2227,   564 }},
  {{  565,  2225,  2229 }},
  {{ 2225,   381,  2228 }},
  {{ 2229,  2225,  2228 }},
  {{ 2229,  2228,   564 }},
  {{  565,  2229,  2231 }},
  {{ 2229,   564,  2230 }},
  {{ 2231,  2229,  2230 }},
  {{ 2231,  2230,   147 }},
  {{  100,  1477,  2233 }},
  {{ 1477,   378,  2232 }},
  {{ 2233,  1477,  2232 }},
  {{ 2233,  2232,   567 }},
  {{  378,  1474,  2235 }},
  {{ 1474,    11,  2234 }},
  {{ 2235,  1474,  2234 }},
  {{ 2235,  2234,   566 }},
  {{  567,  2232,  2236 }},
  {{ 2232,   378,  2235 }},
  {{ 2236,  2232,  2235 }},
  {{ 2236,  2235,   566 }},
  {{  567,  2236,  2238 }},
  {{ 2236,   566,  2237 }},
  {{ 2238,  2236,  2237 }},
  {{ 2238,  2237,   146 }},
  {{  147,  2230,  2240 }},
  {{ 2230,   564,  2239 }},
  {{ 2240,  2230,  2239 }},
  {{ 2240,  2239,   568 }},
  {{  564,  2227,  2241 }},
  {{ 2227,   100,  2233 }},
  {{ 2241,  2227,  2233 }},
  {{ 2241,  2233,   567 }},
  {{  568,  2239,  2242 }},
  {{ 2239,   564,  2241 }},
  {{ 2242,  2239,  2241 }},
  {{ 2242,  2241,   567 }},
  {{  568,  2242,  2243 }},
  {{ 2242,   567,  2238 }},
  {{ 2243,  2242,  2238 }},
  {{ 2243,  2238,   146 }},
  {{  147,  2240,  2245 }},
  {{ 2240,   568,  2244 }},
  {{ 2245,  2240,  2244 }},
  {{ 2245,  2244,   570 }},
  {{  568,  2243,  2247 }},
  {{ 2243,   146,  2246 }},
  {{ 2247,  2243,  2246 }},
  {{ 2247,  2246,   569 }},
  {{  570,  2244,  2248 }},
  {{ 2244,   568,  2247 }},
  {{ 2248,  2244,  2247 }},
  {{ 2248,  2247,   569 }},
  {{  570,  2248,  2250 }},
  {{ 2248,   569,  2249 }},
  {{ 2250,  2248,  2249 }},
  {{ 2250,  2249,    40 }},
  {{   38,  2224,  2252 }},
  {{ 2224,   563,  2251 }},
  {{ 2252,  2224,  2251 }},
  {{ 2252,  2251,   572 }},
  {{  563,  2221,  2254 }},
  {{ 2221,   145,  2253 }},
  {{ 2254,  2221,  2253 }},
  {{ 2254,  2253,   571 }},
  {{  572,  2251,  2255 }},
  {{ 2251,   563,  2254 }},
  {{ 2255,  2251,  2254 }},
  {{ 2255,  2254,   571 }},
  {{  572,  2255,  2257 }},
  {{ 2255,   571,  2256 }},
  {{ 2257,  2255,  2256 }},
  {{ 2257,  2256,   148 }},
  {{  145,  2213,  2259 }},
  {{ 2213,   560,  2258 }},
  {{ 2259,  2213,  2258 }},
  {{ 2259,  2258,   573 }},
  {{  560,  2210,  2260 }},
  {{ 2210,    29,  2226 }},
  {{ 2260,  2210,  2226 }},
  {{ 2260,  2226,   565 }},
  {{  573,  2258,  2261 }},
  {{ 2258,   560,  2260 }},
  {{ 2261,  2258,  2260 }},
  {{ 2261,  2260,   565 }},
  {{  573,  2261,  2262 }},
  {{ 2261,   565,  2231 }},
  {{ 2262,  2261,  2231 }},
  {{ 2262,  2231,   147 }},
  {{  148,  2256,  2264 }},
  {{ 2256,   571,  2263 }},
  {{ 2264,  2256,  2263 }},
  {{ 2264,  2263,   574 }},
  {{  571,  2253,  2265 }},
  {{ 2253,   145,  2259 }},
  {{ 2265,  2253,  2259 }},
  {{ 2265,  2259,   573 }},
  {{  574,  2263,  2266 }},
  {{ 2263,   571,  2265 }},
  {{ 2266,  2263,  2265 }},
  {{ 2266,  2265,   573 }},
  {{  574,  2266,  2267 }},
  {{ 2266,   573,  2262 }},
  {{ 2267,  2266,  2262 }},
  {{ 2267,  2262,   147 }},
  {{  148,  2264,  2269 }},
  {{ 2264,   574,  2268 }},
  {{ 2269,  2264,  2268 }},
  {{ 2269,  2268,   575 }},
  {{  574,  2267,  2270 }},
  {{ 2267,   147,  2245 }},
  {{ 2270,  2267,  2245 }},
  {{ 2270,  2245,   570 }},
  {{  575,  2268,  2271 }},
  {{ 2268,   574,  2270 }},
  {{ 2271,  2268,  2270 }},
  {{ 2271,  2270,   570 }},
  {{  575,  2271,  2272 }},
  {{ 2271,   570,  2250 }},
  {{ 2272,  2271,  2250 }},
  {{ 2272,  2250,    40 }},
  {{   38,  2252,  2186 }},
  {{ 2252,   572,  2273 }},
  {{ 2186,  2252,  2273 }},
  {{ 2186,  2273,   554 }},
  {{  572,  2257,  2275 }},
  {{ 2257,   148,  2274 }},
  {{ 2275,  2257,  2274 }},
  {{ 2275,  2274,   576 }},
  {{  554,  2273,  2276 }},
  {{ 2273,   572,  2275 }},
  {{ 2276,  2273,  2275 }},
  {{ 2276,  2275,   576 }},
  {{  554,  2276,  2189 }},
  {{ 2276,   576,  2277 }},
  {{ 2189,  2276,  2277 }},
  {{ 2189,  2277,   143 }},
  {{  148,  2269,  2279 }},
  {{ 2269,   575,  2278 }},
  {{ 2279,  2269,  2278 }},
  {{ 2279,  2278,   578 }},
  {{  575,  2272,  2281 }},
  {{ 2272,    40,  2280 }},
  {{ 2281,  2272,  2280 }},
  {{ 2281,  2280,   577 }},
  {{  578,  2278,  2282 }},
  {{ 2278,   575,  2281 }},
  {{ 2282,  2278,  2281 }},
  {{ 2282,  2281,   577 }},
  {{  578,  2282,  2284 }},
  {{ 2282,   577,  2283 }},
  {{ 2284,  2282,  2283 }},
  {{ 2284,  2283,   149 }},
  {{  143,  2277,  2286 }},
  {{ 2277,   576,  2285 }},
  {{ 2286,  2277,  2285 }},
  {{ 2286,  2285,   579 }},
  {{  576,  2274,  2287 }},
  {{ 2274,   148,  2279 }},
  {{ 2287,  2274,  2279 }},
  {{ 2287,  2279,   578 }},
  {{  579,  2285,  2288 }},
  {{ 2285,   576,  2287 }},
  {{ 2288,  2285,  2287 }},
  {{ 2288,  2287,   578 }},
  {{  579,  2288,  2289 }},
  {{ 2288,   578,  2284 }},
  {{ 2289,  2288,  2284 }},
  {{ 2289,  2284,   149 }},
  {{  143,  2286,  2198 }},
  {{ 2286,   579,  2290 }},
  {{ 2198,  2286,  2290 }},
  {{ 2198,  2290,   557 }},
  {{  579,  2289,  2292 }},
  {{ 2289,   149,  2291 }},
  {{ 2292,  2289,  2291 }},
  {{ 2292,  2291,   580 }},
  {{  557,  2290,  2293 }},
  {{ 2290,   579,  2292 }},
  {{ 2293,  2290,  2292 }},
  {{ 2293,  2292,   580 }},
  {{  557,  2293,  2201 }},
  {{ 2293,   580,  2294 }},
  {{ 2201,  2293,  2294 }},
  {{ 2201,  2294,    12 }},
  {{   11,  1726,  2234 }},
  {{ 1726,   440,  2295 }},
  {{ 2234,  1726,  2295 }},
  {{ 2234,  2295,   566 }},
  {{  440,  1723,  2297 }},
  {{ 1723,   115,  2296 }},
  {{ 2297,  1723,  2296 }},
  {{ 2297,  2296,   581 }},
  {{  566,  2295,  2298 }},
  {{ 2295,   440,  2297 }},
  {{ 2298,  2295,  2297 }},
  {{ 2298,  2297,   581 }},
  {{  566,  2298,  2237 }},
  {{ 2298,   581,  2299 }},
  {{ 2237,  2298,  2299 }},
  {{ 2237,  2299,   146 }},
  {{  115,  1715,  2301 }},
  {{ 1715,   437,  2300 }},
  {{ 2301,  1715,  2300 }},
  {{ 2301,  2300,   583 }},
  {{  437,  1712,  2303 }},
  {{ 1712,    32,  2302 }},
  {{ 2303,  1712,  2302 }},
  {{ 2303,  2302,   582 }},
  {{  583,  2300,  2304 }},
  {{ 2300,   437,  2303 }},
  {{ 2304,  2300,  2303 }},
  {{ 2304,  2303,   582 }},
  {{  583,  2304,  2306 }},
  {{ 2304,   582,  2305 }},
  {{ 2306,  2304,  2305 }},
  {{ 2306,  2305,   150 }},
  {{  146,  2299,  2308 }},
  {{ 2299,   581,  2307 }},
  {{ 2308,  2299,  2307 }},
  {{ 2308,  2307,   584 }},
  {{  581,  2296,  2309 }},
  {{ 2296,   115,  2301 }},
  {{ 2309,  2296,  2301 }},
  {{ 2309,  2301,   583 }},
  {{  584,  2307,  2310 }},
  {{ 2307,   581,  2309 }},
  {{ 2310,  2307,  2309 }},
  {{ 2310,  2309,   583 }},
  {{  584,  2310,  2311 }},
  {{ 2310,   583,  2306 }},
  {{ 2311,  2310,  2306 }},
  {{ 2311,  2306,   150 }},
  {{  146,  2308,  2246 }},
  {{ 2308,   584,  2312 }},
  {{ 2246,  2308,  2312 }},
  {{ 2246,  2312,   569 }},
  {{  584,  2311,  2314 }},
  {{ 2311,   150,  2313 }},
  {{ 2314,  2311,  2313 }},
  {{ 2314,  2313,   585 }},
  {{  569,  2312,  2315 }},
  {{ 2312,   584,  2314 }},
  {{ 2315,  2312,  2314 }},
  {{ 2315,  2314,   585 }},
  {{  569,  2315,  2249 }},
  {{ 2315,   585,  2316 }},
  {{ 2249,  2315,  2316 }},
  {{ 2249,  2316,    40 }},
  {{   32,  1681,  2318 }},
  {{ 1681,   429,  2317 }},
  {{ 2318,  1681,  2317 }},
  {{ 2318,  2317,   587 }},
  {{  429,  1678,  2320 }},
  {{ 1678,   112,  2319 }},
  {{ 2320,  1678,  2319 }},
  {{ 2320,  2319,   586 }},
  {{  587,  2317,  2321 }},
  {{ 2317,   429,  2320 }},
  {{ 2321,  2317,  2320 }},
  {{ 2321,  2320,   586 }},
  {{  587,  2321,  2323 }},
  {{ 2321,   586,  2322 }},
  {{ 2323,  2321,  2322 }},
  {{ 2323,  2322,   152 }},
  {{  112,  1669,  2325 }},
  {{ 1669,   426,  2324 }},
  {{ 2325,  1669,  2324 }},
  {{ 2325,  2324,   589 }},
  {{  426,  1666,  2327 }},
  {{ 1666,     7,  2326 }},
  {{ 2327,  1666,  2326 }},
  {{ 2327,  2326,   588 }},
  {{  589,  2324,  2328 }},
  {{ 2324,   426,  2327 }},
  {{ 2328,  2324,  2327 }},
  {{ 2328,  2327,   588 }},
  {{  589,  2328,  2330 }},
  {{ 2328,   588,  2329 }},
  {{ 2330,  2328,  2329 }},
  {{ 2330,  2329,   151 }},
  {{  152,  2322,  2332 }},
  {{ 2322,   586,  2331 }},
  {{ 2332,  2322,  2331 }},
  {{ 2332,  2331,   590 }},
  {{  586,  2319,  2333 }},
  {{ 2319,   112,  2325 }},
  {{ 2333,  2319,  2325 }},
  {{ 2333,  2325,   589 }},
  {{  590,  2331,  2334 }},
  {{ 2331,   586,  2333 }},
  {{ 2334,  2331,  2333 }},
  {{ 2334,  2333,   589 }},
  {{  590,  2334,  2335 }},
  {{ 2334,   589,  2330 }},
  {{ 2335,  2334,  2330 }},
  {{ 2335,  2330,   151 }},
  {{  152,  2332,  2337 }},
  {{ 2332,   590,  2336 }},
  {{ 2337,  2332,  2336 }},
  {{ 2337,  2336,   592 }},
  {{  590,  2335,  2339 }},
  {{ 2335,   151,  2338 }},
  {{ 2339,  2335,  2338 }},
  {{ 2339,  2338,   591 }},
  {{  592,  2336,  2340 }},
  {{ 2336,   590,  2339 }},
  {{ 2340,  2336,  2339 }},
  {{ 2340,  2339,   591 }},
  {{  592,  2340,  2342 }},
  {{ 2340,   591,  2341 }},
  {{ 2342,  2340,  2341 }},
  {{ 2342,  2341,    41 }},
  {{   40,  2316,  2344 }},
  {{ 2316,   585,  2343 }},
  {{ 2344,  2316,  2343 }},
  {{ 2344,  2343,   594 }},
  {{  585,  2313,  2346 }},
  {{ 2313,   150,  2345 }},
  {{ 2346,  2313,  2345 }},
  {{ 2346,  2345,   593 }},
  {{  594,  2343,  2347 }},
  {{ 2343,   585,  2346 }},
  {{ 2347,  2343,  2346 }},
  {{ 2347,  2346,   593 }},
  {{  594,  2347,  2349 }},
  {{ 2347,   593,  2348 }},
  {{ 2349,  2347,  2348 }},
  {{ 2349,  2348,   153 }},
  {{  150,  2305,  2351 }},
  {{ 2305,   582,  2350 }},
  {{ 2351,  2305,  2350 }},
  {{ 2351,  2350,   595 }},
  {{  582,  2302,  2352 }},
  {{ 2302,    32,  2318 }},
  {{ 2352,  2302,  2318 }},
  {{ 2352,  2318,   587 }},
  {{  595,  2350,  2353 }},
  {{ 2350,   582,  2352 }},
  {{ 2353,  2350,  2352 }},
  {{ 2353,  2352,   587 }},
  {{  595,  2353,  2354 }},
  {{ 2353,   587,  2323 }},
  {{ 2354,  2353,  2323 }},
  {{ 2354,  2323,   152 }},
  {{  153,  2348,  2356 }},
  {{ 2348,   593,  2355 }},
  {{ 2356,  2348,  2355 }},
  {{ 2356,  2355,   596 }},
  {{  593,  2345,  2357 }},
  {{ 2345,   150,  2351 }},
  {{ 2357,  2345,  2351 }},
  {{ 2357,  2351,   595 }},
  {{  596,  2355,  2358 }},
  {{ 2355,   593,  2357 }},
  {{ 2358,  2355,  2357 }},
  {{ 2358,  2357,   595 }},
  {{  596,  2358,  2359 }},
  {{ 2358,   595,  2354 }},
  {{ 2359,  2358,  2354 }},
  {{ 2359,  2354,   152 }},
  {{  153,  2356,  2361 }},
  {{ 2356,   596,  2360 }},
  {{ 2361,  2356,  2360 }},
  {{ 2361,  2360,   597 }},
  {{  596,  2359,  2362 }},
  {{ 2359,   152,  2337 }},
  {{ 2362,  2359,  2337 }},
  {{ 2362,  2337,   592 }},
  {{  597,  2360,  2363 }},
  {{ 2360,   596,  2362 }},
  {{ 2363,  2360,  2362 }},
  {{ 2363,  2362,   592 }},
  {{  597,  2363,  2364 }},
  {{ 2363,   592,  2342 }},
  {{ 2364,  2363,  2342 }},
  {{ 2364,  2342,    41 }},
  {{   40,  2344,  2280 }},
  {{ 2344,   594,  2365 }},
  {{ 2280,  2344,  2365 }},
  {{ 2280,  2365,   577 }},
  {{  594,  2349,  2367 }},
  {{ 2349,   153,  2366 }},
  {{ 2367,  2349,  2366 }},
  {{ 2367,  2366,   598 }},
  {{  577,  2365,  2368 }},
  {{ 2365,   594,  2367 }},
  {{ 2368,  2365,  2367 }},
  {{ 2368,  2367,   598 }},
  {{  577,  2368,  2283 }},
  {{ 2368,   598,  2369 }},
  {{ 2283,  2368,  2369 }},
  {{ 2283,  2369,   149 }},
  {{  153,  2361,  2371 }},
  {{ 2361,   597,  2370 }},
  {{ 2371,  2361,  2370 }},
  {{ 2371,  2370,   600 }},
  {{  597,  2364,  2373 }},
  {{ 2364,    41,  2372 }},
  {{ 2373,  2364,  2372 }},
  {{ 2373,  2372,   599 }},
  {{  600,  2370,  2374 }},
  {{ 2370,   597,  2373 }},
  {{ 2374,  2370,  2373 }},
  {{ 2374,  2373,   599 }},
  {{  600,  2374,  2376 }},
  {{ 2374,   599,  2375 }},
  {{ 2376,  2374,  2375 }},
  {{ 2376,  2375,   154 }},
  {{  149,  2369,  2378 }},
  {{ 2369,   598,  2377 }},
  {{ 2378,  2369,  2377 }},
  {{ 2378,  2377,   601 }},
  {{  598,  2366,  2379 }},
  {{ 2366,   153,  2371 }},
  {{ 2379,  2366,  2371 }},
  {{ 2379,  2371,   600 }},
  {{  601,  2377,  2380 }},
  {{ 2377,   598,  2379 }},
  {{ 2380,  2377,  2379 }},
  {{ 2380,  2379,   600 }},
  {{  601,  2380,  2381 }},
  {{ 2380,   600,  2376 }},
  {{ 2381,  2380,  2376 }},
  {{ 2381,  2376,   154 }},
  {{  149,  2378,  2291 }},
  {{ 2378,   601,  2382 }},
  {{ 2291,  2378,  2382 }},
  {{ 2291,  2382,   580 }},
  {{  601,  2381,  2384 }},
  {{ 2381,   154,  2383 }},
  {{ 2384,  2381,  2383 }},
  {{ 2384,  2383,   602 }},
  {{  580,  2382,  2385 }},
  {{ 2382,   601,  2384 }},
  {{ 2385,  2382,  2384 }},
  {{ 2385,  2384,   602 }},
  {{  580,  2385,  2294 }},
  {{ 2385,   602,  2386 }},
  {{ 2294,  2385,  2386 }},
  {{ 2294,  2386,    12 }},
  {{    7,  1918,  2326 }},
  {{ 1918,   488,  2387 }},
  {{ 2326,  1918,  2387 }},
  {{ 2326,  2387,   588 }},
  {{  488,  1915,  2389 }},
  {{ 1915,   127,  2388 }},
  {{ 2389,  1915,  2388 }},
  {{ 2389,  2388,   603 }},
  {{  588,  2387,  2390 }},
  {{ 2387,   488,  2389 }},
  {{ 2390,  2387,  2389 }},
  {{ 2390,  2389,   603 }},
  {{  588,  2390,  2329 }},
  {{ 2390,   603,  2391 }},
  {{ 2329,  2390,  2391 }},
  {{ 2329,  2391,   151 }},
  {{  127,  1907,  2393 }},
  {{ 1907,   485,  2392 }},
  {{ 2393,  1907,  2392 }},
  {{ 2393,  2392,   605 }},
  {{  485,  1904,  2395 }},
  {{ 1904,    35,  2394 }},
  {{ 2395,  1904,  2394 }},
  {{ 2395,  2394,   604 }},
  {{  605,  2392,  2396 }},
  {{ 2392,   485,  2395 }},
  {{ 2396,  2392,  2395 }},
  {{ 2396,  2395,   604 }},
  {{  605,  2396,  2398 }},
  {{ 2396,   604,  2397 }},
  {{ 2398,  2396,  2397 }},
  {{ 2398,  2397,   155 }},
  {{  151,  2391,  2400 }},
  {{ 2391,   603,  2399 }},
  {{ 2400,  2391,  2399 }},
  {{ 2400,  2399,   606 }},
  {{  603,  2388,  2401 }},
  {{ 2388,   127,  2393 }},
  {{ 2401,  2388,  2393 }},
  {{ 2401,  2393,   605 }},
  {{  606,  2399,  2402 }},
  {{ 2399,   603,  2401 }},
  {{ 2402,  2399,  2401 }},
  {{ 2402,  2401,   605 }},
  {{  606,  2402,  2403 }},
  {{ 2402,   605,  2398 }},
  {{ 2403,  2402,  2398 }},
  {{ 2403,  2398,   155 }},
  {{  151,  2400,  2338 }},
  {{ 2400,   606,  2404 }},
  {{ 2338,  2400,  2404 }},
  {{ 2338,  2404,   591 }},
  {{  606,  2403,  2406 }},
  {{ 2403,   155,  2405 }},
  {{ 2406,  2403,  2405 }},
  {{ 2406,  2405,   607 }},
  {{  591,  2404,  2407 }},
  {{ 2404,   606,  2406 }},
  {{ 2407,  2404,  2406 }},
  {{ 2407,  2406,   607 }},
  {{  591,  2407,  2341 }},
  {{ 2407,   607,  2408 }},
  {{ 2341,  2407,  2408 }},
  {{ 2341,  2408,    41 }},
  {{   35,  1873,  2410 }},
  {{ 1873,   477,  2409 }},
  {{ 2410,  1873,  2409 }},
  {{ 2410,  2409,   609 }},
  {{  477,  1870,  2412 }},
  {{ 1870,   124,  2411 }},
  {{ 2412,  1870,  2411 }},
  {{ 2412,  2411,   608 }},
  {{  609,  2409,  2413 }},
  {{ 2409,   477,  2412 }},
  {{ 2413,  2409,  2412 }},
  {{ 2413,  2412,   608 }},
  {{  609,  2413,  2415 }},
  {{ 2413,   608,  2414 }},
  {{ 2415,  2413,  2414 }},
  {{ 2415,  2414,   157 }},
  {{  124,  1861,  2417 }},
  {{ 1861,   474,  2416 }},
  {{ 2417,  1861,  2416 }},
  {{ 2417,  2416,   611 }},
  {{  474,  1858,  2419 }},
  {{ 1858,     8,  2418 }},
  {{ 2419,  1858,  2418 }},
  {{ 2419,  2418,   610 }},
  {{  611,  2416,  2420 }},
  {{ 2416,   474,  2419 }},
  {{ 2420,  2416,  2419 }},
  {{ 2420,  2419,   610 }},
  {{  611,  2420,  2422 }},
  {{ 2420,   610,  2421 }},
  {{ 2422,  2420,  2421 }},
  {{ 2422,  2421,   156 }},
  {{  157,  2414,  2424 }},
  {{ 2414,   608,  2423 }},
  {{ 2424,  2414,  2423 }},
  {{ 2424,  2423,   612 }},
  {{  608,  2411,  2425 }},
  {{ 2411,   124,  2417 }},
  {{ 2425,  2411,  2417 }},
  {{ 2425,  2417,   611 }},
  {{  612,  2423,  2426 }},
  {{ 2423,   608,  2425 }},
  {{ 2426,  2423,  2425 }},
  {{ 2426,  2425,   611 }},
  {{  612,  2426,  2427 }},
  {{ 2426,   611,  2422 }},
  {{ 2427,  2426,  2422 }},
  {{ 2427,  2422,   156 }},
  {{  157,  2424,  2429 }},
  {{ 2424,   612,  2428 }},
  {{ 2429,  2424,  2428 }},
  {{ 2429,  2428,   614 }},
  {{  612,  2427,  2431 }},
  {{ 2427,   156,  2430 }},
  {{ 2431,  2427,  2430 }},
  {{ 2431,  2430,   613 }},
  {{  614,  2428,  2432 }},
  {{ 2428,   612,  2431 }},
  {{ 2432,  2428,  2431 }},
  {{ 2432,  2431,   613 }},
  {{  614,  2432,  2434 }},
  {{ 2432,   613,  2433 }},
  {{ 2434,  2432,  2433 }},
  {{ 2434,  2433,    42 }},
  {{   41,  2408,  2436 }},
  {{ 2408,   607,  2435 }},
  {{ 2436,  2408,  2435 }},
  {{ 2436,  2435,   616 }},
  {{  607,  2405,  2438 }},
  {{ 2405,   155,  2437 }},
  {{ 2438,  2405,  2437 }},
  {{ 2438,  2437,   615 }},
  {{  616,  2435,  2439 }},
  {{ 2435,   607,  2438 }},
  {{ 2439,  2435,  2438 }},
  {{ 2439,  2438,   615 }},
  {{  616,  2439,  2441 }},
  {{ 2439,   615,  2440 }},
  {{ 2441,  2439,  2440 }},
  {{ 2441,  2440,   158 }},
  {{  155,  2397,  2443 }},
  {{ 2397,   604,  2442 }},
  {{ 2443,  2397,  2442 }},
  {{ 2443,  2442,   617 }},
  {{  604,  2394,  2444 }},
  {{ 2394,    35,  2410 }},
  {{ 2444,  2394,  2410 }},
  {{ 2444,  2410,   609 }},
  {{  617,  2442,  2445 }},
  {{ 2442,   604,  2444 }},
  {{ 2445,  2442,  2444 }},
  {{ 2445,  2444,   609 }},
  {{  617,  2445,  2446 }},
  {{ 2445,   609,  2415 }},
  {{ 2446,  2445,  2415 }},
  {{ 2446,  2415,   157 }},
  {{  158,  2440,  2448 }},
  {{ 2440,   615,  2447 }},
  {{ 2448,  2440,  2447 }},
  {{ 2448,  2447,   618 }},
  {{  615,  2437,  2449 }},
  {{ 2437,   155,  2443 }},
  {{ 2449,  2437,  2443 }},
  {{ 2449,  2443,   617 }},
  {{  618,  2447,  2450 }},
  {{ 2447,   615,  2449 }},
  {{ 2450,  2447,  2449 }},
  {{ 2450,  2449,   617 }},
  {{  618,  2450,  2451 }},
  {{ 2450,   617,  2446 }},
  {{ 2451,  2450,  2446 }},
  {{ 2451,  2446,   157 }},
  {{  158,  2448,  2453 }},
  {{ 2448,   618,  2452 }},
  {{ 2453,  2448,  2452 }},
  {{ 2453,  2452,   619 }},
  {{  618,  2451,  2454 }},
  {{ 2451,   157,  2429 }},
  {{ 2454,  2451,  2429 }},
  {{ 2454,  2429,   614 }},
  {{  619,  2452,  2455 }},
  {{ 2452,   618,  2454 }},
  {{ 2455,  2452,  2454 }},
  {{ 2455,  2454,   614 }},
  {{  619,  2455,  2456 }},
  {{ 2455,   614,  2434 }},
  {{ 2456,  2455,  2434 }},
  {{ 2456,  2434,    42 }},
  {{   41,  2436,  2372 }},
  {{ 2436,   616,  2457 }},
  {{ 2372,  2436,  2457 }},
  {{ 2372,  2457,   599 }},
  {{  616,  2441,  2459 }},
  {{ 2441,   158,  2458 }},
  {{ 2459,  2441,  2458 }},
  {{ 2459,  2458,   620 }},
  {{  599,  2457,  2460 }},
  {{ 2457,   616,  2459 }},
  {{ 2460,  2457,  2459 }},
  {{ 2460,  2459,   620 }},
  {{  599,  2460,  2375 }},
  {{ 2460,   620,  2461 }},
  {{ 2375,  2460,  2461 }},
  {{ 2375,  2461,   154 }},
  {{  158,  2453,  2463 }},
  {{ 2453,   619,  2462 }},
  {{ 2463,  2453,  2462 }},
  {{ 2463,  2462,   622 }},
  {{  619,  2456,  2465 }},
  {{ 2456,    42,  2464 }},
  {{ 2465,  2456,  2464 }},
  {{ 2465,  2464,   621 }},
  {{  622,  2462,  2466 }},
  {{ 2462,   619,  2465 }},
  {{ 2466,  2462,  2465 }},
  {{ 2466,  2465,   621 }},
  {{  622,  2466,  2468 }},
  {{ 2466,   621,  2467 }},
  {{ 2468,  2466,  2467 }},
  {{ 2468,  2467,   159 }},
  {{  154,  2461,  2470 }},
  {{ 2461,   620,  2469 }},
  {{ 2470,  2461,  2469 }},
  {{ 2470,  2469,   623 }},
  {{  620,  2458,  2471 }},
  {{ 2458,   158,  2463 }},
  {{ 2471,  2458,  2463 }},
  {{ 2471,  2463,   622 }},
  {{  623,  2469,  2472 }},
  {{ 2469,   620,  2471 }},
  {{ 2472,  2469,  2471 }},
  {{ 2472,  2471,   622 }},
  {{  623,  2472,  2473 }},
  {{ 2472,   622,  2468 }},
  {{ 2473,  2472,  2468 }},
  {{ 2473,  2468,   159 }},
  {{  154,  2470,  2383 }},
  {{ 2470,   623,  2474 }},
  {{ 2383,  2470,  2474 }},
  {{ 2383,  2474,   602 }},
  {{  623,  2473,  2476 }},
  {{ 2473,   159,  2475 }},
  {{ 2476,  2473,  2475 }},
  {{ 2476,  2475,   624 }},
  {{  602,  2474,  2477 }},
  {{ 2474,   623,  2476 }},
  {{ 2477,  2474,  2476 }},
  {{ 2477,  2476,   624 }},
  {{  602,  2477,  2386 }},
  {{ 2477,   624,  2478 }},
  {{ 2386,  2477,  2478 }},
  {{ 2386,  2478,    12 }},
  {{    8,  2102,  2418 }},
  {{ 2102,   532,  2479 }},
  {{ 2418,  2102,  2479 }},
  {{ 2418,  2479,   610 }},
  {{  532,  2099,  2481 }},
  {{ 2099,   137,  2480 }},
  {{ 2481,  2099,  2480 }},
  {{ 2481,  2480,   625 }},
  {{  610,  2479,  2482 }},
  {{ 2479,   532,  2481 }},
  {{ 2482,  2479,  2481 }},
  {{ 2482,  2481,   625 }},
  {{  610,  2482,  2421 }},
  {{ 2482,   625,  2483 }},
  {{ 2421,  2482,  2483 }},
  {{ 2421,  2483,   156 }},
  {{  137,  2091,  2485 }},
  {{ 2091,   529,  2484 }},
  {{ 2485,  2091,  2484 }},
  {{ 2485,  2484,   627 }},
  {{  529,  2088,  2487 }},
  {{ 2088,    37,  2486 }},
  {{ 2487,  2088,  2486 }},
  {{ 2487,  2486,   626 }},
  {{  627,  2484,  2488 }},
  {{ 2484,   529,  2487 }},
  {{ 2488,  2484,  2487 }},
  {{ 2488,  2487,   626 }},
  {{  627,  2488,  2490 }},
  {{ 2488,   626,  2489 }},
  {{ 2490,  2488,  2489 }},
  {{ 2490,  2489,   160 }},
  {{  156,  2483,  2492 }},
  {{ 2483,   625,  2491 }},
  {{ 2492,  2483,  2491 }},
  {{ 2492,  2491,   628 }},
  {{  625,  2480,  2493 }},
  {{ 2480,   137,  2485 }},
  {{ 2493,  2480,  2485 }},
  {{ 2493,  2485,   627 }},
  {{  628,  2491,  2494 }},
  {{ 2491,   625,  2493 }},
  {{ 2494,  2491,  2493 }},
  {{ 2494,  2493,   627 }},
  {{  628,  2494,  2495 }},
  {{ 2494,   627,  2490 }},
  {{ 2495,  2494,  2490 }},
  {{ 2495,  2490,   160 }},
  {{  156,  2492,  2430 }},
  {{ 2492,   628,  2496 }},
  {{ 2430,  2492,  2496 }},
  {{ 2430,  2496,   613 }},
  {{  628,  2495,  2498 }},
  {{ 2495,   160,  2497 }},
  {{ 2498,  2495,  2497 }},
  {{ 2498,  2497,   629 }},
  {{  613,  2496,  2499 }},
  {{ 2496,   628,  2498 }},
  {{ 2499,  2496,  2498 }},
  {{ 2499,  2498,   629 }},
  {{  613,  2499,  2433 }},
  {{ 2499,   629,  2500 }},
  {{ 2433,  2499,  2500 }},
  {{ 2433,  2500,    42 }},
  {{   37,  2057,  2502 }},
  {{ 2057,   521,  2501 }},
  {{ 2502,  2057,  2501 }},
  {{ 2502,  2501,   631 }},
  {{  521,  2054,  2504 }},
  {{ 2054,   134,  2503 }},
  {{ 2504,  2054,  2503 }},
  {{ 2504,  2503,   630 }},
  {{  631,  2501,  2505 }},
  {{ 2501,   521,  2504 }},
  {{ 2505,  2501,  2504 }},
  {{ 2505,  2504,   630 }},
  {{  631,  2505,  2507 }},
  {{ 2505,   630,  2506 }},
  {{ 2507,  2505,  2506 }},
  {{ 2507,  2506,   161 }},
  {{  134,  2045,  2509 }},
  {{ 2045,   518,  2508 }},
  {{ 2509,  2045,  2508 }},
  {{ 2509,  2508,   632 }},
  {{  518,  2042,  2510 }},
  {{ 2042,     9,  2104 }},
  {{ 2510,  2042,  2104 }},
  {{ 2510,  2104,   534 }},
  {{  632,  2508,  2511 }},
  {{ 2508,   518,  2510 }},
  {{ 2511,  2508,  2510 }},
  {{ 2511,  2510,   534 }},
  {{  632,  2511,  2512 }},
  {{ 2511,   534,  2109 }},
  {{ 2512,  2511,  2109 }},
  {{ 2512,  2109,   139 }},
  {{  161,  2506,  2514 }},
  {{ 2506,   630,  2513 }},
  {{ 2514,  2506,  2513 }},
  {{ 2514,  2513,   633 }},
  {{  630,  2503,  2515 }},
  {{ 2503,   134,  2509 }},
  {{ 2515,  2503,  2509 }},
  {{ 2515,  2509,   632 }},
  {{  633,  2513,  2516 }},
  {{ 2513,   630,  2515 }},
  {{ 2516,  2513,  2515 }},
  {{ 2516,  2515,   632 }},
  {{  633,  2516,  2517 }},
  {{ 2516,   632,  2512 }},
  {{ 2517,  2516,  2512 }},
  {{ 2517,  2512,   139 }},
  {{  161,  2514,  2519 }},
  {{ 2514,   633,  2518 }},
  {{ 2519,  2514,  2518 }},
  {{ 2519,  2518,   634 }},
  {{  633,  2517,  2520 }},
  {{ 2517,   139,  2123 }},
  {{ 2520,  2517,  2123 }},
  {{ 2520,  2123,   539 }},
  {{  634,  2518,  2521 }},
  {{ 2518,   633,  2520 }},
  {{ 2521,  2518,  2520 }},
  {{ 2521,  2520,   539 }},
  {{  634,  2521,  2522 }},
  {{ 2521,   539,  2128 }},
  {{ 2522,  2521,  2128 }},
  {{ 2522,  2128,    39 }},
  {{   42,  2500,  2524 }},
  {{ 2500,   629,  2523 }},
  {{ 2524,  2500,  2523 }},
  {{ 2524,  2523,   636 }},
  {{  629,  2497,  2526 }},
  {{ 2497,   160,  2525 }},
  {{ 2526,  2497,  2525 }},
  {{ 2526,  2525,   635 }},
  {{  636,  2523,  2527 }},
  {{ 2523,   629,  2526 }},
  {{ 2527,  2523,  2526 }},
  {{ 2527,  2526,   635 }},
  {{  636,  2527,  2529 }},
  {{ 2527,   635,  2528 }},
  {{ 2529,  2527,  2528 }},
  {{ 2529,  2528,   162 }},
  {{  160,  2489,  2531 }},
  {{ 2489,   626,  2530 }},
  {{ 2531,  2489,  2530 }},
  {{ 2531,  2530,   637 }},
  {{  626,  2486,  2532 }},
  {{ 2486,    37,  2502 }},
  {{ 2532,  2486,  2502 }},
  {{ 2532,  2502,   631 }},
  {{  637,  2530,  2533 }},
  {{ 2530,   626,  2532 }},
  {{ 2533,  2530,  2532 }},
  {{ 2533,  2532,   631 }},
  {{  637,  2533,  2534 }},
  {{ 2533,   631,  2507 }},
  {{ 2534,  2533,  2507 }},
  {{ 2534,  2507,   161 }},
  {{  162,  2528,  2536 }},
  {{ 2528,   635,  2535 }},
  {{ 2536,  2528,  2535 }},
  {{ 2536,  2535,   638 }},
  {{  635,  2525,  2537 }},
  {{ 2525,   160,  2531 }},
  {{ 2537,  2525,  2531 }},
  {{ 2537,  2531,   637 }},
  {{  638,  2535,  2538 }},
  {{ 2535,   635,  2537 }},
  {{ 2538,  2535,  2537 }},
  {{ 2538,  2537,   637 }},
  {{  638,  2538,  2539 }},
  {{ 2538,   637,  2534 }},
  {{ 2539,  2538,  2534 }},
  {{ 2539,  2534,   161 }},
  {{  162,  2536,  2541 }},
  {{ 2536,   638,  2540 }},
  {{ 2541,  2536,  2540 }},
  {{ 2541,  2540,   639 }},
  {{  638,  2539,  2542 }},
  {{ 2539,   161,  2519 }},
  {{ 2542,  2539,  2519 }},
  {{ 2542,  2519,   634 }},
  {{  639,  2540,  2543 }},
  {{ 2540,   638,  2542 }},
  {{ 2543,  2540,  2542 }},
  {{ 2543,  2542,   634 }},
  {{  639,  2543,  2544 }},
  {{ 2543,   634,  2522 }},
  {{ 2544,  2543,  2522 }},
  {{ 2544,  2522,    39 }},
  {{   42,  2524,  2464 }},
  {{ 2524,   636,  2545 }},
  {{ 2464,  2524,  2545 }},
  {{ 2464,  2545,   621 }},
  {{  636,  2529,  2547 }},
  {{ 2529,   162,  2546 }},
  {{ 2547,  2529,  2546 }},
  {{ 2547,  2546,   640 }},
  {{  621,  2545,  2548 }},
  {{ 2545,   636,  2547 }},
  {{ 2548,  2545,  2547 }},
  {{ 2548,  2547,   640 }},
  {{  621,  2548,  2467 }},
  {{ 2548,   640,  2549 }},
  {{ 2467,  2548,  2549 }},
  {{ 2467,  2549,   159 }},
  {{  162,  2541,  2551 }},
  {{ 2541,   639,  2550 }},
  {{ 2551,  2541,  2550 }},
  {{ 2551,  2550,   641 }},
  {{  639,  2544,  2552 }},
  {{ 2544,    39,  2178 }},
  {{ 2552,  2544,  2178 }},
  {{ 2552,  2178,   553 }},
  {{  641,  2550,  2553 }},
  {{ 2550,   639,  2552 }},
  {{ 2553,  2550,  2552 }},
  {{ 2553,  2552,   553 }},
  {{  641,  2553,  2554 }},
  {{ 2553,   553,  2183 }},
  {{ 2554,  2553,  2183 }},
  {{ 2554,  2183,   144 }},
  {{  159,  2549,  2556 }},
  {{ 2549,   640,  2555 }},
  {{ 2556,  2549,  2555 }},
  {{ 2556,  2555,   642 }},
  {{  640,  2546,  2557 }},
  {{ 2546,   162,  2551 }},
  {{ 2557,  2546,  2551 }},
  {{ 2557,  2551,   641 }},
  {{  642,  2555,  2558 }},
  {{ 2555,   640,  2557 }},
  {{ 2558,  2555,  2557 }},
  {{ 2558,  2557,   641 }},
  {{  642,  2558,  2559 }},
  {{ 2558,   641,  2554 }},
  {{ 2559,  2558,  2554 }},
  {{ 2559,  2554,   144 }},
  {{  159,  2556,  2475 }},
  {{ 2556,   642,  2560 }},
  {{ 2475,  2556,  2560 }},
  {{ 2475,  2560,   624 }},
  {{  642,  2559,  2561 }},
  {{ 2559,   144,  2197 }},
  {{ 2561,  2559,  2197 }},
  {{ 2561,  2197,   558 }},
  {{  624,  2560,  2562 }},
  {{ 2560,   642,  2561 }},
  {{ 2562,  2560,  2561 }},
  {{ 2562,  2561,   558 }},
  {{  624,  2562,  2478 }},
  {{ 2562,   558,  2202 }},
  {{ 2478,  2562,  2202 }},
  {{ 2478,  2202,    12 }}
} ;


/*---------------------------------------------------------------------------*/
MRI_SURFACE *
ICOreadOverAlloc(char *fname, double pct_over)
{
  ICOSOHEDRON *ico ;
  int         fno, vno, n1, n2, n, vn ;
  MRI_SURFACE *mris ;
  VERTEX      *v ;
  FACE        *f ;

  ico = read_icosahedron(fname) ;

  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    vno = ico->faces[fno].vno[1] ;
    ico->faces[fno].vno[1] = ico->faces[fno].vno[2] ;
    ico->faces[fno].vno[2] = vno ;
  }

#if 0
  mris = MRISoverAlloc(pct_over*ico->nvertices,pct_over*ico->nfaces,
                       ico->nvertices, ico->nfaces) ;
#endif
  mris = MRISoverAlloc((int)(pct_over*ico->nvertices),(int)(pct_over*ico->nfaces),
                       ico->nvertices, ico->nfaces) ;

  /* position vertices */
  for (vno = 0 ; vno < ico->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;

#if 0
    v->x = 100.0*ico->vertices[vno].x ;
    v->y = 100.0*ico->vertices[vno].y ;
    v->z = 100.0*ico->vertices[vno].z ;
#else
    v->x = ico->vertices[vno].x ;
    v->y = ico->vertices[vno].y ;
    v->z = ico->vertices[vno].z ;
#endif
  }

  /* fill in faces, and count # of faces each vertex is part of */
  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == 15)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      f->v[n] = ico->faces[fno].vno[n]-1 ;  /* make it zero-based */
      v = &mris->vertices[f->v[n]] ;
      v->num++ ;
      v->vnum += 2 ;   /* will remove duplicates later */
    }
  }

  for (vno = 0 ; vno < ico->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->v = (int *)calloc(v->vnum/2, sizeof(int)) ;
    if (!v->v)
      ErrorExit(ERROR_NOMEMORY, "ICread: could not allocate %dth vertex list.",
                vno) ;
    v->vnum = 0 ;
  }

  /* now build list of neighbors */
  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == 3)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;

      /* now add an edge to other 2 vertices if not already in list */
      for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
      {
        if (n1 == n)   /* don't connect vertex to itself */
          continue ;
        vn = ico->faces[fno].vno[n1]-1 ;  /* make it zero-based */

        /* now check to make sure it's not a duplicate */
        for (n2 = 0 ; n2 < v->vnum ; n2++)
        {
          if (v->v[n2] == vn)
          {
            vn = -1 ; /* mark it as a duplicate */
            break ;
          }
        }
        if (vn >= 0)
          v->v[v->vnum++] = vn ;
      }
    }
  }

  /* now allocate face arrays in vertices */
  for (vno = 0 ; vno < ico->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->f = (int *)calloc(v->num, sizeof(int)) ;
    if (!v->f)
      ErrorExit(ERROR_NOMEMORY,"ICOread: could not allocate %d faces",v->num);
    v->n = (unsigned char *)calloc(v->num,sizeof(unsigned char));
    if (!v->n)
      ErrorExit(ERROR_NOMEMORY, "ICOread: could not allocate %d nbrs", v->n);
    v->num = 0 ;   /* for use as counter in next section */
    v->dist = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist)
      ErrorExit(ERROR_NOMEMORY,
                "ICOread: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
    v->dist_orig = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist_orig)
      ErrorExit(ERROR_NOMEMORY,
                "ICOread: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
    v->vtotal = v->vnum ;
  }

  /* fill in face indices in vertex structures */
  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;
      v->n[v->num] = n ; v->f[v->num++] = fno ;
    }
  }

  MRIScomputeMetricProperties(mris) ;
#if 0
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    float dot ;
    int   ano ;

    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
      
    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
       */
    v = &mris->vertices[f->v[0]] ;
    dot = v->x * f->nx + v->y * f->ny + v->z * f->nz;
    if (dot < 0.0f)   /* not in same direction, area < 0 and reverse n */
    {
      f->area *= -1.0f ;
      f->nx *= -1.0f; f->ny *= -1.0f; f->nz *= -1.0f;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        f->angle[ano] *= -1.0f ;
    }
  }
#endif
  mris->type = MRIS_ICO_SURFACE ;
  free(ico->vertices) ;
  free(ico->faces) ;
  free(ico) ;
  return(mris) ;
}

/*---------------------------------------------------------------------------*/
MRI_SURFACE *
ICOread(char *fname)
{
  ICOSOHEDRON *ico ;
  int         fno, vno, n1, n2, n, vn ;
  MRI_SURFACE *mris ;
  VERTEX      *v ;
  FACE        *f ;

  ico = read_icosahedron(fname) ;

  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    vno = ico->faces[fno].vno[1] ;
    ico->faces[fno].vno[1] = ico->faces[fno].vno[2] ;
    ico->faces[fno].vno[2] = vno ;
  }

  mris = MRISalloc(ico->nvertices, ico->nfaces) ;

  /* position vertices */
  for (vno = 0 ; vno < ico->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;

#if 0
    v->x = 100.0*ico->vertices[vno].x ;
    v->y = 100.0*ico->vertices[vno].y ;
    v->z = 100.0*ico->vertices[vno].z ;
#else
    v->x = ico->vertices[vno].x ;
    v->y = ico->vertices[vno].y ;
    v->z = ico->vertices[vno].z ;
#endif
  }

  /* fill in faces, and count # of faces each vertex is part of */
  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == 15)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      f->v[n] = ico->faces[fno].vno[n]-1 ;  /* make it zero-based */
      v = &mris->vertices[f->v[n]] ;
      v->num++ ;
      v->vnum += 2 ;   /* will remove duplicates later */
    }
  }

  for (vno = 0 ; vno < ico->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->v = (int *)calloc(v->vnum/2, sizeof(int)) ;
    if (!v->v)
      ErrorExit(ERROR_NOMEMORY, "ICread: could not allocate %dth vertex list.",
                vno) ;
    v->vnum = 0 ;
  }

  /* now build list of neighbors */
  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (fno == 3)
      DiagBreak() ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;

      /* now add an edge to other 2 vertices if not already in list */
      for (n1 = 0 ; n1 < VERTICES_PER_FACE ; n1++)
      {
        if (n1 == n)   /* don't connect vertex to itself */
          continue ;
        vn = ico->faces[fno].vno[n1]-1 ;  /* make it zero-based */

        /* now check to make sure it's not a duplicate */
        for (n2 = 0 ; n2 < v->vnum ; n2++)
        {
          if (v->v[n2] == vn)
          {
            vn = -1 ; /* mark it as a duplicate */
            break ;
          }
        }
        if (vn >= 0)
          v->v[v->vnum++] = vn ;
      }
    }
  }

  /* now allocate face arrays in vertices */
  for (vno = 0 ; vno < ico->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->f = (int *)calloc(v->num, sizeof(int)) ;
    if (!v->f)
      ErrorExit(ERROR_NOMEMORY,"ICOread: could not allocate %d faces",v->num);
    v->n = (unsigned char *)calloc(v->num,sizeof(unsigned char));
    if (!v->n)
      ErrorExit(ERROR_NOMEMORY, "ICOread: could not allocate %d nbrs", v->n);
    v->num = 0 ;   /* for use as counter in next section */
    v->dist = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist)
      ErrorExit(ERROR_NOMEMORY,
                "ICOread: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
    v->dist_orig = (float *)calloc(v->vnum, sizeof(float)) ;
    if (!v->dist_orig)
      ErrorExit(ERROR_NOMEMORY,
                "ICOread: could not allocate list of %d "
                "dists at v=%d", v->vnum, vno) ;
    v->vtotal = v->vnum ;
  }

  /* fill in face indices in vertex structures */
  for (fno = 0 ; fno < ico->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      v = &mris->vertices[f->v[n]] ;
      v->n[v->num] = n ; v->f[v->num++] = fno ;
    }
  }

  MRIScomputeMetricProperties(mris) ;
#if 0
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    float dot ;
    int   ano ;

    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
      
    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
       */
    v = &mris->vertices[f->v[0]] ;
    dot = v->x * f->nx + v->y * f->ny + v->z * f->nz;
    if (dot < 0.0f)   /* not in same direction, area < 0 and reverse n */
    {
      f->area *= -1.0f ;
      f->nx *= -1.0f; f->ny *= -1.0f; f->nz *= -1.0f;
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
        f->angle[ano] *= -1.0f ;
    }
  }
#endif
  mris->type = MRIS_ICO_SURFACE ;
  free(ico->vertices) ;
  free(ico->faces) ;
  free(ico) ;
  return(mris) ;
}

/*---------------------------------------------------------------------------*/
static ICOSOHEDRON *
read_icosahedron(char *fname)
{
  FILE        *fp ;
  char        line[200], *cp ;
  int         vno, fno, vno1, vno2, vno3, n, nvertices, nfaces ;
  float       x, y, z ;
  ICOSOHEDRON *ico ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NOFILE, "read_icosahedron: could not open %s", fname));

  ico = (ICOSOHEDRON *)calloc(1, sizeof(ICOSOHEDRON)) ;

  fgetl(line, 150, fp) ;   /* discard # of vertices */
  sscanf(line, "%d", &nvertices) ;
  ico->nvertices = nvertices ;
  ico->vertices = 
    (IC_VERTEX *)calloc(nvertices, sizeof(IC_VERTEX)) ;
  if (!ico->vertices)
    ErrorExit(ERROR_NOMEMORY, "read_ico: could not allocate vertex list") ;


  /* first read vertices */
  n = 0 ;
  while ((cp = fgetl(line, 150, fp)) != NULL)
  {
    if (sscanf(cp, "%d %f %f %f\n", &vno, &x, &y, &z) < 4)
      break ;
    ico->vertices[vno-1].x = x ;
    ico->vertices[vno-1].y = y ;
    ico->vertices[vno-1].z = z ;
    if (++n >= ico->nvertices)
      break ;
  }
  n = 0 ; fgetl(line, 150, fp) ;   /* discard # of faces */
  sscanf(line, "%d", &nfaces) ;
  ico->nfaces = nfaces ;
  ico->faces = 
    (IC_FACE *)calloc(ico->nfaces, sizeof(IC_FACE)) ;
  if (!ico->faces)
    ErrorExit(ERROR_NOMEMORY, "read_ico: could not allocate vertex list") ;
  while ((cp = fgetl(line, 150, fp)) != NULL)
  {
    if (sscanf(cp, "%d %d %d %d\n", &fno, &vno1, &vno2, &vno3) < 4)
      break ;
    ico->faces[fno-1].vno[0] = vno1 ;
    ico->faces[fno-1].vno[1] = vno2 ;
    ico->faces[fno-1].vno[2] = vno3 ;
    if (++n >= ico->nfaces)
      break ;
  }
  fclose(fp) ;
  return(ico) ;
}


/*---------------------------------------------------------------------------*/
int
ICOreadVertexPositions(MRI_SURFACE *mris, char *fname, int which)
{
  ICOSOHEDRON *ico ;
  int         vno ;
  VERTEX      *v ;

  ico = read_icosahedron(fname) ;
  if (!ico)
    return(Gerror) ;

  if (ico->nvertices != mris->nvertices || ico->nfaces != mris->nfaces)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "ICOreadVertexPositions: different size surfaces")) ;

  /* position vertices */
  for (vno = 0 ; vno < ico->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;

    switch (which)
    {
    default:
    case CURRENT_VERTICES:
      v->x = ico->vertices[vno].x ;
      v->y = ico->vertices[vno].y ;
      v->z = ico->vertices[vno].z ;
      break ;
    case CANONICAL_VERTICES:
      v->cx = ico->vertices[vno].x ;
      v->cy = ico->vertices[vno].y ;
      v->cz = ico->vertices[vno].z ;
      break ;
    case ORIGINAL_VERTICES:
      v->origx = ico->vertices[vno].x ;
      v->origy = ico->vertices[vno].y ;
      v->origz = ico->vertices[vno].z ;
      break ;
    }
  }
  return(NO_ERROR) ;
}

/*---------------------------------------------------------------------------*/
/*------------------------------------------------------------------
  MRI_SURFACE *ReadIcoByOrder() -- loads an icosaheron given it's
  "order".  Reads the MRI_DIR environment variable, then loads the
  file "MRI_DIR/lib/bem/ic%d.tri" where %d is the "order".  As read in,
  the vertices of the icosahedron are on the unit sphere. If RescaleFactor
  is greater than zero, the coordinates of each vertex are multiplied
  by RescaleFactor (ie, the new sphere has a radius of RescaleFactor).
  Note: I (doug) don't know I'm using the term "order" properly.
  I did not know what else to call this particular property.
  -------------------------------------------------------------------*/
MRI_SURFACE *ReadIcoByOrder(int IcoOrder, float RescaleFactor)
{
  /* multiple changes - dhagler 08/25/04 */
  char funcname[STRLEN]="ReadIcoByOrder";
  char *MRI_DIR, trifile[2048];
  MRI_SURFACE *surf;
  VERTEX *v;
  int vtx;

  MRI_DIR = getenv("MRI_DIR") ;
  sprintf(trifile,"%s/lib/bem/ic%d.tri",MRI_DIR,IcoOrder);
  MsgPrintf("%s: reading icosahedron %s\n", funcname,trifile) ;
  surf = ICOread(trifile);
  if(surf == NULL)
    ErrorExit(ERROR_BADFILE, "error reading icosahedron file\n",funcname);

  /* rescale xyz to something other than the unit sphere */
  if(RescaleFactor > 0 && RescaleFactor != 1){
    for (vtx = 0 ; vtx < surf->nvertices ; vtx++){
      v = &surf->vertices[vtx] ;
      v->x *= RescaleFactor; 
      v->y *= RescaleFactor; 
      v->z *= RescaleFactor;
    }
    MRIScomputeMetricProperties(surf);
    /* this was done in ICOread, but need to recalc */
    /* because coordinates are different after rescale */
  }
  return(surf);
}

/*---------------------------------------------------------------------------*/
/*-----------------------------------------------------------------
  MRI_SURFACE *ReadIcoByNVtxs -- loads an icosaheron given the number
  of vertices.  Converts the number of vertices into the "order", then
  calls ReadIcoByOrder().
  ------------------------------------------------------------------*/
MRI_SURFACE *ReadIcoByNVtxs(int nIcoVtxs, float RescaleFactor)
{
  MRI_SURFACE *surf;
  int IcoOrder;

  IcoOrder = IcoOrderFromNVtxs(nIcoVtxs);
  if(IcoOrder == -1) return(NULL);

  surf = ReadIcoByOrder(IcoOrder,RescaleFactor);
  return(surf);
}

/*---------------------------------------------------------------------------*/
/*-------------------------------------------------------------
  IcoOrderFromNVtxs() - returns the "order" of the icosahedron
  give the number of vertices in the ico. The icosahedrons are
  stored in a file whose name is based on the order.
  Note: I (doug) don't know I'm using the term "order" properly.
  I did not know what else to call this particular property.
  --------------------------------------------------------------*/
int IcoOrderFromNVtxs(int nIcoVtxs)
{
  int IcoOrder = -1;

  switch(nIcoVtxs){
  case       12: IcoOrder = 0; break;
  case       42: IcoOrder = 1; break;
  case      162: IcoOrder = 2; break;
  case      642: IcoOrder = 3; break;
  case     2562: IcoOrder = 4; break;
  case    10242: IcoOrder = 5; break;
  case    40962: IcoOrder = 6; break;
  case   163842: IcoOrder = 7; break;
  }
  if(IcoOrder == -1)
    fprintf(stderr,"ERROR: nIcoVtxs = %d does not match an IcoOrder\n",
      nIcoVtxs);
  return(IcoOrder);
}

/*---------------------------------------------------------------------------*/
int IcoNVtxsFromOrder(int IcoOrder)
{
  int nIcoVtxs = -1;

  switch(IcoOrder){
  case       0: nIcoVtxs =     12; break;
  case       1: nIcoVtxs =     42; break;
  case       2: nIcoVtxs =    162; break;
  case       3: nIcoVtxs =    642; break;
  case       4: nIcoVtxs =   2562; break;
  case       5: nIcoVtxs =  10242; break;
  case       6: nIcoVtxs =  40962; break;
  case       7: nIcoVtxs = 163842; break;
  }
  if(nIcoVtxs == -1)
    fprintf(stderr,"ERROR: IcoOrder = %d is out of range (0-7)\n",
      IcoOrder);
  return(nIcoVtxs);
}


/*---------------------------------------------------------------*/
int GetICOOrderFromValFile(char *filename)
{
  int nIcoVtxs,IcoOrder;
  char funcname[STRLEN]="GetICOOrderFromValFile";

  nIcoVtxs = GetNVtxsFromWFile(filename);

  IcoOrder = IcoOrderFromNVtxs(nIcoVtxs);
  if(IcoOrder < 0){
    ErrorExit(ERROR_BADFILE,"%s: ### number of vertices = %d, does not mach ico\n",
      funcname, nIcoVtxs);
  }
  
  return(IcoOrder);
}
/*---------------------------------------------------------------*/


/*##########################################################################*/
/* numerical recipes */

/* global data for use with a variety of macros defined in nrutil.h */
float sqrarg;
double dsqrarg;
double dmaxarg1,dmaxarg2;
double dminarg1,dminarg2;
float maxarg1,maxarg2;
float minarg1,minarg2;
long lmaxarg1,lmaxarg2;
long lminarg1,lminarg2;
int imaxarg1,imaxarg2;
int iminarg1,iminarg2;


#define NR_END 1
#define FREE_ARG char*


float *fvector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
  if (!v) 
    ErrorExit(ERROR_NOMEMORY, "could not allocate vector(%ld, %ld)",
              nl, nh);
  return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v)
    ErrorExit(ERROR_NOMEMORY, "could not allocate int vector(%ld, %ld)",
              nl, nh);
  return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;

  v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
  if (!v)
    ErrorExit(ERROR_NOMEMORY, "could not allocate char vector(%ld, %ld)",
              nl, nh);
  return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;

  v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v)
    ErrorExit(ERROR_NOMEMORY, "could not allocate long vector(%ld, %ld)",
              nl, nh);
  return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v)
    ErrorExit(ERROR_NOMEMORY, "could not allocate double vector(%ld, %ld)",
              nl, nh);
  return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
  if (!m) 
    ErrorExit(ERROR_NOMEMORY, "could not allocate matrix(%ld, %ld)",
              nrow, ncol);
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
  if (!m[nrl])
    ErrorExit(ERROR_NOMEMORY, "could not allocate matrix(%ld, %ld) array",
              nrow, ncol);
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  if (!m) 
    ErrorExit(ERROR_NOMEMORY, "could not allocate double matrix(%ld, %ld)",
              nrow, ncol);
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl])
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate double matrix(%ld, %ld) array",
              nrow, ncol);
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
  if (!m) 
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate int matrix(%ld, %ld)", nrow, ncol);
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
  if (!m[nrl]) 
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate int matrix(%ld, %ld) array", nrow, ncol);
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
  long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;

  /* allocate array of pointers to rows */
  m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m) 
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate int submatrix(%ld, %ld)", nrow, ncol);
  m += NR_END;
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;

  /* allocate pointers to rows */
  m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
  if (!m)
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate int conversion matrix(%ld, %ld)",nrow,ncol);
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;

  /* allocate pointers to pointers to rows */
  t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) 
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate tensor(%ld, %ld, %ld)",nrow,ncol,ndep);
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl])
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate tensor(%ld, %ld, %ld) array",nrow,ncol,ndep);
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
  if (!t[nrl][ncl])
    ErrorExit(ERROR_NOMEMORY, 
              "could not allocate tensor(%ld, %ld, %ld) array",nrow,ncol,ndep);
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
  long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
}
/* (C) Copr. 1986-92 Numerical Recipes Software "25B&W#(,. */

/***************************************************************************/

void tred2(float **a, int n, float d[], float e[])
{
  int l,k,j,i;
  float scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0f;
    if (l > 1) {
      for (k=1;k<=l;k++)
        scale += (float)fabs(a[i][k]);
      if (scale == 0.0f)
        e[i]=a[i][l];
      else {
        for (k=1;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0f ? -(float)sqrt(h) : (float)sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0f;
        for (j=1;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0f;
          for (k=1;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=1;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0f;
  e[1]=0.0f;
  /* Contents of this loop can be omitted if eigenvectors not
      wanted except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
        g=0.0f;
        for (k=1;k<=l;k++)
          g += a[i][k]*a[k][j];
        for (k=1;k<=l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0f;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0f;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software "25B&W#(,. */

/***************************************************************************/

static float at,bt,ct;

/*
   Given a matrix a[1..m][1..n], the routine computes its SVD,
   A = U.W.Vt. The matrix U replaces a on output. The diagonal matrix of 
   singular values W is output as a vector w[1..n]. The matrix V 
   (not the transpose) is output as v[1..n][1..n].
*/
int
svdcmp(float **a, int m, int n, float *w, float **v)
{
  int flag,i = 1,its,j,jj,k,l = 1,nm = 1;
  float c,f,h,s,x,y,z;
  float anorm=0.0,g=0.0,scale=0.0;
  float *rv1;

  if (m < n) 
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "SVDCMP: You must augment A with extra zero rows"));

  rv1=fvector(1,n);
  for (i=1;i<=n;i++) 
  {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k=i;k<=m;k++) {
          a[k][i] /= scale;
          s += a[k][i]*a[k][i];
        }
        f=a[i][i];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][i]=f-g;
        if (i != n) {
          for (j=l;j<=n;j++) {
            for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
            f=s/h;
            for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
          }
        }
        for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k=l;k<=n;k++) {
          a[i][k] /= scale;
          s += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        a[i][l]=f-g;
        for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
        if (i != m) {
          for (j=l;j<=m;j++) {
            for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
            for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
          }
        }
        for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
        for (j=l;j<=n;j++)
          v[j][i]=(a[i][j]/a[i][l])/g;
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
          for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
        }
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=n;i>=1;i--) {
    l=i+1;
    g=w[i];
    if (i < n)
      for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      if (i != n) {
        for (j=l;j<=n;j++) {
          for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
          f=(s/a[i][i])*g;
          for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
        }
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else {
      for (j=i;j<=m;j++) a[j][i]=0.0;
    }
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
        nm=l-1;
        if (fabs(rv1[l])+anorm == anorm) {
          flag=0;
          break;
        }
        if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
        c=0.0;
        s=1.0;
        for (i=l;i<=k;i++) {
          f=s*rv1[i];
          if (fabs(f)+anorm != anorm) {
            g=w[i];
            h=PYTHAG(f,g);
            w[i]=h;
            h=1.0/h;
            c=g*h;
            s=(-f*h);
            for (j=1;j<=m;j++) {
              y=a[j][nm];
              z=a[j][i];
              a[j][nm]=y*c+z*s;
              a[j][i]=z*c-y*s;
            }
          }
        }
      }
      z=w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j=1;j<=n;j++) v[j][k]=(-v[j][k]);
        }
        break;
      }
      if (its > 5*m*n)
        ErrorReturn(ERROR_BADPARM,
                    (ERROR_BADPARM, "SVDCMP: No convergence in %d iterations",
                     5*m*n));
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=PYTHAG(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
        i=j+1;
        g=rv1[i];
        y=w[i];
        h=s*g;
        g=c*g;
        z=PYTHAG(f,h);
        rv1[j]=z;
        c=f/z;
        s=h/z;
        f=x*c+g*s;
        g=g*c-x*s;
        h=y*s;
        y=y*c;
        for (jj=1;jj<=n;jj++) {
          x=v[jj][j];
          z=v[jj][i];
          v[jj][j]=x*c+z*s;
          v[jj][i]=z*c-x*s;
        }
        z=PYTHAG(f,h);
        w[j]=z;
        if (z) {
          z=1.0/z;
          c=f*z;
          s=h*z;
        }
        f=(c*g)+(s*y);
        x=(c*y)-(s*g);
        for (jj=1;jj<=m;jj++) {
          y=a[jj][j];
          z=a[jj][i];
          a[jj][j]=y*c+z*s;
          a[jj][i]=z*c-y*s;
        }
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_vector(rv1,1,n);
  return(NO_ERROR) ;
}

/***************************************************************************/

#define MAX_ITER(n)  (n > 30 ? n : 30)

float pythag(float a, float b)
{
  float absa,absb;
  absa=(float)fabs(a);
  absb=(float)fabs(b);
  if (absa > absb) return absa*(float)sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0f ? 0.0f : absb*(float)sqrt(1.0+SQR(absa/absb)));
}

int
tqli(float d[], float e[], int n, float **z)
{
  float pythag(float a, float b);
  int m,l,iter,i,k;
  float s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0f;
  for (l=1;l<=n;l++) 
  {
    iter=0;
    do 
    {
      for (m=l;m<=n-1;m++) 
      {
        dd=(float)fabs(d[m])+(float)fabs(d[m+1]);
        if ((float)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) 
      {
        if (iter++ == MAX_ITER(n)) 
#if 1
          return(Gerror = ERROR_BADPARM) ;
#else
          ErrorReturn(ERROR_BADPARM,
                      (ERROR_BADPARM, "tqli: too many iterations"));
#endif
        g=(d[l+1]-d[l])/(2.0f*e[l]);
        r=pythag(g,1.0f);
        g=d[m]-d[l]+e[l]/(g+(float)SIGN(r,g));
        s=c=1.0f;
        p=0.0f;
        for (i=m-1;i>=l;i--) 
        {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          if (r == 0.0f) 
          {
            d[i+1] -= p;
            e[m]=0.0f;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0f*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
          for (k=1;k<=n;k++) 
          {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        if (r == 0.0 && i >= l) 
          continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0f;
      }
    } while (m != l);
  }
  return(NO_ERROR) ;
}
/* (C) Copr. 1986-92 Numerical Recipes Software "25B&W#(,. */

/*##########################################################################*/
/*---------------------------------------------------------------------------*/
/* tkregister */
/*---------------------------------------------------------------------------*/

/* from MRIio.c */

char *
lcalloc(size_t nmemb,size_t size)
{
  char *p;
  
  p = (char *)calloc(nmemb,size);
  if (p==NULL)
    ErrorExit(ERROR_NOMEMORY,"lcalloc failed\n");
  return p;
}

void
file_name(char *fpref, char *fname, int num, char *form)
{
  char ext[10];
  
  sprintf(ext,form,num);
  strcpy(fname,fpref);
  strcat(fname,ext);
}

void
buffer_to_image(unsigned char *buf, unsigned char **im,int ysize,int xsize)
{
  int i,j;
  unsigned long k;
  float sum;
  
  k=0;
  sum = 0;
  for (i=0;i<ysize;i++)
    for (j=0;j<xsize;j++)
    {
      im[i][j] = buf[k++];
      sum += im[i][j];
    }
}


/*---------------------------------------------------------------------------*/


