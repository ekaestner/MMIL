#ifndef SURFLIB_H
#define SURFLIB_H

/* surflib.h */

/*##########################################################################*/
/* includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctype.h>
#include <stdarg.h>
#include <pwd.h>
#include <errno.h>
#include <sys/time.h>
#include <time.h>
#ifdef Darwin
#  include <unistd.h>
#  include <ppc/limits.h>
#else
#  include <values.h>
#endif

/*##########################################################################*/
/* misc */
/*---------------------------------------------------------------------------*/
#define MATCH(A,B)   (!strcmp(A,B))
#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#define MTEST(ptr) \
  if((ptr)==NULL) ( fprintf(stderr,"*** Memory allocation error\n"), exit(0))
#define UNDEFSTR "unknown"
#define BIGFLOAT 1e10

/*##########################################################################*/
/* error */
#define NO_ERROR              0
#define ERROR_NONE            NO_ERROR
#define ERROR_NO_FILE         -1
#define ERROR_NOFILE          ERROR_NO_FILE
#define ERROR_NO_MEMORY       -2
#define ERROR_NOMEMORY        ERROR_NO_MEMORY
#define ERROR_UNSUPPORTED     -3
#define ERROR_BADPARM         -4
#define ERROR_BAD_PARM        ERROR_BADPARM
#define ERROR_BADFILE         -5
#define ERROR_BAD_FILE        ERROR_BADFILE
#define ERROR_SIZE            -6

#define DIAG_VERBOSE 0x10000000L   /* allows 2 levels for each type */
#define DIAG_VERBOSE_ON  (Gdiag & DIAG_VERBOSE)
#define DIAG_SHOW       0x00000040L
#define ERROR_FNAME  "error.log"
#define ErrorReturn(ret, args)  { ErrorPrintf args ; return(ret) ; }

extern int Gerror ;    /* global error value */

/*##########################################################################*/
/* matrix */
/*---------------------------------------------------------------------------*/
#define MATRIX_CELT(m,r,c)      (((COMPLEX_FLOAT **)m->rptr)[r]+c)
#define MATRIX_RELT(m,r,c)      (m->rptr[r]+c)
#define MATRIX_ELT(m,r,c)       (m->type == MATRIX_REAL ? \
                                  *MATRIX_RELT(m,r,c) : \
                                  *MATRIX_CELT(m,r,c))
#define MATRIX_PTR(m,r,c)       (m->type == MATRIX_REAL ? \
                                  MATRIX_RELT(m,r,c) : \
                                  (float *)MATRIX_CELT(m,r,c))

#define MATRIX_CELT_REAL(m,r,c)  (MATRIX_CELT(m,r,c)->real)
#define MATRIX_CELT_IMAG(m,r,c)  (MATRIX_CELT(m,r,c)->imag)

#define MATRIX_REAL        1
#define MATRIX_COMPLEX     2

#define VectorAlloc(n, type)       MatrixAlloc(n, 1, type)
#define RVectorAlloc(n, type)      MatrixAlloc(1, n, type)
#define VectorFree(pm)             MatrixFree(pm)
#define VectorAdd(v1, v2, v3)      MatrixAdd(v1, v2, v3)
#define VectorSubtract(v1,v2,v3)   MatrixSubtract(v1, v2, v3)
#define VectorScalarMul(v1,val,v2) MatrixScalarMul(v1, val, v2)
#define VectorCopy(v1, v2)         MatrixCopy(v1, v2)
#define VectorClear(v)             MatrixClear(v)
#define VectorTranspose(vsrc,vdst) MatrixTranspose(vsrc, vdst)
#define VectorAsciiWriteInto       MatrixAsciiWriteInto
#define VectorAsciiReadFrom        MatrixAsciiReadFrom

#define VECTOR_ELT(v,i)            ((v)->rptr[i][1])
#define RVECTOR_ELT(v,i)            ((v)->rptr[1][i])
#define VECTOR3_LOAD(v,x,y,z)    (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z) ;
#define VECTOR_LOAD   VECTOR3_LOAD
#define V3_LOAD       VECTOR3_LOAD

/*##########################################################################*/
/* transform */
/*---------------------------------------------------------------------------*/
typedef struct
{
  double    m[4][4];
} Transform;

#define MAX_LINE_LENGTH   10000
#define TRANSFORM_STRING "Linear_Transform"
#define FABS(x) fabs((double)x)

/*##########################################################################*/
/* icosahedron */
/*---------------------------------------------------------------------------*/
#define ICO_VERTICES_PER_FACE  3
#define ICO4_NVERTICES    2562
#define ICO4_NFACES       5120

/*---------------------------------------------------------------------------*/
typedef struct
{
  float x, y, z ;
} IC_VERTEX ;

/*---------------------------------------------------------------------------*/
typedef struct
{
  int vno[ICO_VERTICES_PER_FACE] ;
} IC_FACE ;

/*---------------------------------------------------------------------------*/
typedef struct 
{
  int       nvertices ;
  int       nfaces ;
  IC_VERTEX *vertices ;
  IC_FACE   *faces ;
} ICOSOHEDRON ;

/*---------------------------------------------------------------------------*/
extern IC_VERTEX ic2562_vertices[] ;
extern IC_FACE   ic2562_faces[] ;

/*##########################################################################*/
/* mrisurf constants */
#ifndef TRUE
#  define TRUE 1
#endif
#ifndef FALSE
#  define FALSE 0
#endif

#ifndef uchar
#  define uchar unsigned char
#endif

#define LEFT_HEMISPHERE     0
#define RIGHT_HEMISPHERE    1
#define SLICE_THICKNESS           1
#define START_Y                 (-128)

#define MRIS_BINARY_QUADRANGLE_FILE    0    /* homegrown */
#define MRIS_ASCII_TRIANGLE_FILE       1    /* homegrown */
#define MRIS_GEO_TRIANGLE_FILE         2    /* movie.byu format */
#define MRIS_ICO_SURFACE               3
#define MRIS_TRIANGULAR_SURFACE        MRIS_ICO_SURFACE
#define MRIS_ICO_FILE                  4
#define MRIS_VTK_FILE                  5

#define MRIS_BINARY_FILE    0
#define MRIS_ASCII_FILE     1
#define MRIS_GEO_FILE       2    /* movie.byu format */

#define MRIS_SURFACE               0
#define MRIS_PATCH                 1
#define MRIS_CUT                   MRIS_PATCH
#define MRIS_PLANE                 2
#define MRIS_ELLIPSOID             3
#define MRIS_SPHERE                4
#define MRIS_PARAMETERIZED_SPHERE  5
#define MRIS_RIGID_BODY            6

#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER  (-3 & 0x00ffffff)
#define NEW_VERSION_MAGIC_NUMBER  16777215

#define ORIGINAL_VERTICES   0
#define GOOD_VERTICES       1
#define TMP_VERTICES        2
#define CANONICAL_VERTICES  3
#define CURRENT_VERTICES    4
#define INFLATED_VERTICES   5
#define FLATTENED_VERTICES  6
#define PIAL_VERTICES       7

#define MAX_4_NEIGHBORS     100
#define MAX_3_NEIGHBORS     70
#define MAX_2_NEIGHBORS     20
#define MAX_1_NEIGHBORS     8
#define MAX_NEIGHBORS       (400)

#define VERTICES_PER_FACE    3
#define ANGLES_PER_TRIANGLE  3

#define MATRIX_REAL      1
#define MATRIX_COMPLEX   2

#define STRLEN      1000

#define RAN   0.001   /* one thousandth of a millimeter */

/*##########################################################################*/
/* mrisurf macros */
#define IS_QUADRANGULAR(mris)  (mris->type == MRIS_BINARY_QUADRANGLE_FILE)
#define FZERO(f)     (fabs(f) < 0.0000001F)

#define VECTOR_ELT(v,i)        ((v)->rptr[i][1])
#define VECTOR3_LOAD(v,x,y,z)    (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z) ;
#define VECTOR_LOAD              VECTOR3_LOAD
#define V3_X(v)      (VECTOR_ELT(v,1))
#define V3_Y(v)      (VECTOR_ELT(v,2))
#define V3_Z(v)      (VECTOR_ELT(v,3))
#define V3_LEN(v)    (sqrt(V3_X(v)*V3_X(v)+V3_Y(v)*V3_Y(v)+V3_Z(v)*V3_Z(v)))
#define V3_DOT(va,vb) (V3_X(va)*V3_X(vb)+V3_Y(va)*V3_Y(vb)+V3_Z(va)*V3_Z(vb))
#define V3_CROSS_PRODUCT(va,vb,vc) \
                 V3_X(vc) = V3_Y(va)*V3_Z(vb)- V3_Z(va)*V3_Y(vb),  \
                 V3_Y(vc) = V3_Z(va)*V3_X(vb)- V3_X(va)*V3_Z(vb),  \
                 V3_Z(vc) = V3_X(va)*V3_Y(vb)- V3_Y(va)*V3_X(vb) ;
#define V3_SCALAR_MUL(va,s,vb)  (V3_X(vb)=V3_X(va)*s,\
                                 V3_Y(vb)=V3_Y(va)*s,\
                                 V3_Z(vb)=V3_Z(va)*s)
#define V3_NORMALIZE(va,vb)  { float len = (V3_LEN(va)) ; \
                                  if (FZERO(len)) len = 1.0f ; \
                                  else len = 1.0f / len ; \
                                  V3_SCALAR_MUL(va,len,vb) ; }
#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y,\
                                                 v1->z-v0->z)
#define EVEN(n)      ((((n) / 2) * 2) == n)
#define ODD(n)       (!EVEN(n))
#define ISEVEN       EVEN
#define ISODD        ODD

#define VectorAlloc(n, type)       MatrixAlloc(n, 1, type)
#define VectorFree(pm)             MatrixFree(pm)

#define nint(f)   ((int)(rint((double)f)))  /* nint not on kamares! */

#ifndef UINT
#define UINT         unsigned int
#endif

/*##########################################################################*/
/* mrisurf typedefs */

/*---------------------------------------------------------------------------*/
typedef union
{
  long  l ;
  float f ;
  int   i ;
  unsigned int ui ;
  char  buf[4] ;
  short s[2] ;
} SWAP_LONG ;

/*---------------------------------------------------------------------------*/
typedef union
{
  short  s ;
  char   buf[sizeof(short)] ;
} SWAP_SHORT ;

/*---------------------------------------------------------------------------*/
typedef struct
{
  float  real ;
  float  imag ;
} COMPLEX_FLOAT, *CPTR ;

/*---------------------------------------------------------------------------*/
typedef struct
{
  short   type ;
  int     rows ;
  int     cols ;
  float **rptr;    /* pointer to an array of rows */
  float *data;     /* pointer to base of data */
} MATRIX, VECTOR ;

/*---------------------------------------------------------------------------*/
typedef struct face_type_
{
  int    v[VERTICES_PER_FACE];           /* vertex numbers of this face */
  float  nx ;
  float  ny ;
  float  nz ;
  float  area ;
  float  orig_area ;
  float  angle[ANGLES_PER_TRIANGLE] ;
  float  orig_angle[ANGLES_PER_TRIANGLE]  ;
  char   ripflag;                        /* ripped face */
} face_type, FACE ;

/*---------------------------------------------------------------------------*/
typedef struct vertex_type_
{
  float x,y,z;            /* curr position */
  float nx,ny,nz;         /* curr normal */
  float dx, dy, dz ;      /* current change in position */
  float odx, ody, odz ;   /* last change of position (for momentum) */
  float tdx, tdy, tdz ;   /* temporary storage for averaging gradient */
  float  curv;            /* curr curvature */
  float  curvbak ;
  float  val;             /* scalar data value (file: rh.val, sig2-rh.w) */
  float  imag_val ;       /* imaginary part of complex data value */
  float  cx, cy, cz ;     /* coordinates in canonical coordinate system */
  float  tx, ty, tz ;     /* tmp coordinate storage */
  float  origx, origy, origz ;  /* original coordinates */
  float  pialx, pialy, pialz ;  /* pial surface coordinates */
  float  infx, infy, infz; /* inflated coordinates */
  float  fx, fy, fz ;      /* flattened coordinates */
  float e1x, e1y, e1z ;   /* 1st basis vector for the local tangent plane */
  float e2x, e2y, e2z ;   /* 2nd basis vector for the local tangent plane */
  float nc;               /* curr length normal comp */
  float val2;             /* complex comp data value (file: sig3-rh.w) */
  float valbak;           /* scalar data stack */
  float val2bak;          /* complex comp data stack */
  float stat;             /* statistic */
  int undefval;           /* [previously dist=0] */
  int old_undefval;       /* for smooth_val_sparse */
  int fixedval;           /* [previously val=0] */
  float fieldsign;        /* fieldsign--final: -1,0,1 (file: rh.fs) */
  float fsmask;           /* significance mask (file: rh.fm) */
  uchar num;               /* number neighboring faces */
  int   *f;                /* array neighboring face numbers */
  uchar *n;                /* [0-3, num long] */
  uchar vnum;              /* number neighboring vertices */
  int   *v;                /* array neighboring vertex numbers, vnum long */
  uchar  v2num ;            /* number of 2-connected neighbors */
  uchar  v3num ;            /* number of 3-connected neighbors */
  short  vtotal ;      /* total # of neighbors, will be same as one of above*/
  float d ;              /* for distance calculations */
  uchar nsize ;            /* size of neighborhood (e.g. 1, 2, 3) */
  int   annotation;     /* area label (defunct--now from label file name!) */
  char   oripflag,origripflag;  /* cuts flags */
  float theta, phi ;     /* parameterization */
  short  marked;          /* for a variety of uses */
  char   ripflag ;
  char   border;          /* flag */
  float area,origarea ;
  float K ;             /* Gaussian curvature */
  float H ;             /* mean curvature */
  float k1, k2 ;        /* the principal curvatures */
  float *dist ;         /* original distance to neighboring vertices */
  float *dist_orig ;    /* original distance to neighboring vertices */
  char   neg ;           /* 1 if the normal vector is inverted */
  float mean ;
  float mean_imag ;    /* imaginary part of complex statistic */
  float std_error ;
} vertex_type, VERTEX ;

/*---------------------------------------------------------------------------*/
typedef struct
{
  int nvertices;
  unsigned int *vertex_indices;
} STRIP;

/*---------------------------------------------------------------------------*/
typedef struct _area_label
{
  char     name[STRLEN] ;     /* name of region */
  float    cx ;               /* centroid x */
  float    cy ;               /* centroid y */
  float    cz ;               /* centroid z */
  int      label ;            /* an identifier (used as an index) */
} MRIS_AREA_LABEL ;

/*---------------------------------------------------------------------------*/
typedef struct
{
  int          nvertices ;      /* # of vertices on surface */
  int          nfaces ;         /* # of faces on surface */ 
  int          nstrips;
  VERTEX       *vertices ;
  FACE         *faces ;
  STRIP        *strips;
  float        xctr ;
  float        yctr ;
  float        zctr ;
  float        xlo ;
  float        ylo ;
  float        zlo ;
  float        xhi ;
  float        yhi ;
  float        zhi ;
  VERTEX       *v_temporal_pole ;
  VERTEX       *v_frontal_pole ;
  VERTEX       *v_occipital_pole ;
  float        max_curv ;
  float        min_curv ;
  float        total_area ;
  float        orig_area ;
  float        neg_area ;
  float        neg_orig_area ;    /* amount of original surface in folds */
  int          zeros ;
  int          hemisphere ;       /* which hemisphere */
  int          initialized ;
/*
  General_transform transform ;
*/
  Transform    *linear_transform ;
  Transform    *inverse_linear_transform ;
  int          transform_loaded ;
  int          inverse_transform_loaded ;
  int          free_transform ;
  double       radius ;           /* radius (if status==MRIS_SPHERE) */
  float        a, b, c ;          /* ellipsoid parameters */
  char         fname[STRLEN] ;    /* file it was originally loaded from */
  float        Hmin ;             /* min mean curvature */
  float        Hmax ;             /* max mean curvature */
  float        Kmin ;             /* min Gaussian curvature */
  float        Kmax ;             /* max Gaussian curvature */
  double       Ktotal ;           /* total Gaussian curvature */
  int          status ;           /* type of surface (e.g. sphere, plane) */
  int          patch ;            /* if a patch of the surface */
  int          nlabels ;
  MRIS_AREA_LABEL *labels ;       /* nlabels of these (may be null) */
  int          nsize ;            /* size of neighborhoods */
  float        avg_nbrs ;         /* mean # of vertex neighbors */
  void         *vp ;              /* for misc. use */
  float        alpha ;            /* rotation around z-axis */
  float        beta ;             /* rotation around y-axis */
  float        gamma ;            /* rotation around x-axis */
  float        da, db, dg ;       /* old deltas */
  int          type ;             /* what type of surface was this initially*/
  int          max_vertices ;     /* may be bigger than nvertices */
  int          max_faces ;        /* may be bigger than nfaces */
  char         subject_name[STRLEN] ;/* name of the subject */
  float        canon_area ;
  int          noscale ;          /* don't scale by surface area if true */
  float        *dx2 ;       /* extra set of gradient (not always alloced) */
  float        *dy2 ;
  float        *dz2 ;
} MRI_SURFACE, MRIS ;

/*##########################################################################*/
/* prototypes */

/* file io */
short swapShort(short s);
int swapInt(int i);
unsigned int swapUInt(unsigned int ui);
float swapFloat(float f);
int fread1(int *v, FILE *fp);
int fread2(int *v, FILE *fp);
int fread3(int *v, FILE *fp);
int fread4(float *v, FILE *fp);
int freadInt(FILE *fp);
unsigned int freadUInt(FILE *fp);
short freadShort(FILE *fp);
float freadFloat(FILE *fp);
int fwrite1(int v, FILE *fp);
int fwrite2(int v, FILE *fp);
int fwrite3(int v, FILE *fp);
int fwrite4(float v, FILE *fp);
int fwriteInt(int v, FILE *fp);
int fwriteFloat(float f, FILE *fp);
char *FileNamePath(char *fname, char *pathName);
int FileExists(char *fname);
int isadir(char *path);
char *fgetl(char *s, int n, FILE *fp);

/* error */
void ErrorExit(int ecode, char *fmt, ...);
int ErrorPrintf(int ecode, char *fmt, ...);
void DiagBreak(void);
int MsgPrintf(char *fmt, ...);
void setQuiet(int qval);
int getQuiet(void);

/* random */
int setRandomSeed(long seed);
double randomNumber(double low, double hi);
void randseed();

/* bfloats */
int getNumBfloats(char *indir, char *instem, char *infix);
int readFloatHeader(char *fname, int *xsize, int *ysize, int *depth);
int readFloatHeaders(int N, char **input, char *infix, int *xsize, int *ysize, int *depth);
int readFloatImage(char *fname, float **fim, int xsize, int ysize);
int readFloatImageTS(char *fname, float ***fimts, int xsize, int ysize, int tsize);
int writeFloatImage(char *fname, float **fim, int xsize, int ysize);
int writeFloatImageTS(char *fname, float ***fimts, int xsize, int ysize, int tsize);
int writeFloatHeader(char *fname, int xsize, int ysize, int t, int num);
int writeTextImage(char *fname, float **fim, int xsize, int ysize);

/* bshorts */
int getNumBshorts(char *indir, char *instem);
int readBshortHeader(char *fname, int *xsize, int *ysize, int *tsize);
int readBshortImageTS(char *fname, float ***fimts, int xsize, int ysize, int tsize);

/* matrix */
MATRIX *MatrixAlloc(int rows, int cols, int type);
int MatrixFree(MATRIX **pmat);
MATRIX  *MatrixTranspose(MATRIX *mIn, MATRIX *mOut);
MATRIX *MatrixMultiply(MATRIX *m1, MATRIX *m2, MATRIX *m3);
float MatrixConditionNumber(MATRIX *m);
MATRIX *MatrixSVDInverse(MATRIX *m, MATRIX *m_inverse);
MATRIX *MatrixCopy(MATRIX *mIn, MATRIX *mOut);
MATRIX *MatrixEigenSystem(MATRIX *m, float *evalues, MATRIX *m_evectors);
float MatrixSVDEigenValues(MATRIX *m, float *evalues);


double Vector3Angle(VECTOR *v1, VECTOR *v2);
float VectorTripleProduct(VECTOR *v1, VECTOR *v2, VECTOR *v3);

/* surface -- from mrisurf.c */
char *StrUpper(char *str);
int MRISwriteValues(MRI_SURFACE *mris, char *sname);
int MRISwriteImagValues(MRI_SURFACE *mris, char *sname);
int MRISbuildFileName(MRI_SURFACE *mris, char *sname, char *fname);
int MRIScomputeNormals(MRI_SURFACE *mris);
double MRISaverageRadius(MRI_SURFACE *mris);
int MRIScomputeTriangleProperties(MRI_SURFACE *mris);
int MRIScomputeMetricProperties(MRI_SURFACE *mris);
int MRISsaveVertexPositions(MRI_SURFACE *mris, int which);
int MRISstoreCurrentPositions(MRI_SURFACE *mris);
int MRIScomputeCanonicalCoordinates(MRI_SURFACE *mris);
int MRISupdateEllipsoidSurface(MRI_SURFACE *mris);
double MRISaverageCanonicalRadius(MRI_SURFACE *mris);
MRI_SURFACE  *MRISalloc(int nvertices, int nfaces);
MRI_SURFACE
  *MRISoverAlloc(int max_vertices, int max_faces, int nvertices, int nfaces);
MRI_SURFACE  *MRISread(char *fname);
MRI_SURFACE  *MRISreadOverAlloc(char *fname, double pct_over);
int MRISreadValues(MRI_SURFACE *mris, char *sname);
int MRISreadImagValues(MRI_SURFACE *mris, char *fname);
float *MRISreadCurvatureVector(MRI_SURFACE *mris, char *sname);
float *MRISreadNewCurvatureVector(MRI_SURFACE *mris, char *sname);
int MRISremoveTriangleLinks(MRI_SURFACE *mris);
int MRISunrip(MRI_SURFACE *mris);
int MRISfree(MRI_SURFACE **pmris);
int MRISreadPatchNoRemove(MRI_SURFACE *mris, char *pname);
int MRISreadPatch(MRI_SURFACE *mris, char *pname);
int MRISripFaces(MRI_SURFACE *mris);
int MRISupdateSurface(MRI_SURFACE *mris);
int MRISremoveRipped(MRI_SURFACE *mris);
int MRIScomputeSecondFundamentalForm(MRI_SURFACE *mris);

/* don's surface funcs */
int MRISsmoothValues(MRIS *mris, int niter);
int MRISsmoothValuesHeat(MRIS *mris, int niter, float sigma);
int MRISsmoothComplexValues(MRIS *mris, int niter);
int MRISsmoothValuesSparse(MRIS *mris, int niter);
int MRISsmoothComplexValuesSparse(MRIS *mris, int niter);
int MRISsmoothValuesROI(MRIS *mris, int niter);
int MRISsmoothComplexValuesROI(MRIS *mris, int niter);
int MRISsmoothValuesHeatROI(MRIS *mris, int niter, float sigma);
int nverticesSurf(char *subject, char *hemi);
int maxVert(char *fname);
int readSurfVals(char *fname, float *data, int nvertices);
int writeSurfVals(char *fname, float *data, int nvertices);
int writeTextVals(char *fname, float *data, int nvertices);
MRI_SURFACE *openSurface(char *subject, char *hemi, char *surface);
int GetNVtxsFromWFile(char *wfile);
int MRISclearValues(MRI_SURFACE *mris);
int MRISclearAllValues(MRI_SURFACE *mris);
int MRISclearMarks(MRI_SURFACE *mris);
float MRISgaussFWHM(MRIS *mris, int sparseflag);
float MRISgaussCxFWHM(MRIS *mris, int sparseflag);  /* for complex data */

/* transform */
int mrisReadTransform(MRIS *mris, char *mris_fname);
int readTransform(char *filename, Transform *transform);
int invertTransform(Transform *transform, Transform *inverse);
int transformPoint(Transform *transform, double x, double y, double z,
  double *x_trans, double *y_trans, double *z_trans);
int printTransform(Transform *transform, char *transname);
void multTransform(Transform *t, Transform *t1, Transform *t2);
int selectTalairachPoint(MRIS *mris, int *vindex,
  double x_tal, double y_tal, double z_tal);
int selectVertexCoord(MRIS *mris, int *vindex,
  double x, double y, double z);
int fixMNITal(double  xmni, double  ymni, double  zmni,
        double *xtal, double *ytal, double *ztal);

/* icosahedron */
MRI_SURFACE *ICOread(char *fname);
MRI_SURFACE *ICOreadOverAlloc(char *fname, double pct_over);
int ICOreadVertexPositions(MRI_SURFACE *mris, char *fname, int which);
MRI_SURFACE *ReadIcoByOrder(int IcoOrder, float RescaleFactor);
MRI_SURFACE *ReadIcoByNVtxs(int nIcoVtxs, float RescaleFactor);
int IcoOrderFromNVtxs(int nIcoVtxs);
int IcoNVtxsFromOrder(int IcoOrder);
int GetICOOrderFromValFile(char *filename);

/* for tkregister */
char *lcalloc(size_t nmemb,size_t size) ;
void buffer_to_image(unsigned char *buf,unsigned char**im,int ysize,int xsize);
void file_name(char *fpref, char *fname, int num, char *form) ;

#endif /* SURFLIB_H */

/*##########################################################################*/
/* numerical recipes */
#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

extern float sqrarg;
#ifndef SQR
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#endif

extern double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

extern double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

extern double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

extern float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

#ifndef MAX
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
  (maxarg1) : (maxarg2))
#endif

extern float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

extern long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

extern long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

extern int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

extern int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0f ? (float)fabs(a) : (float)-fabs(a))
#endif

void nrerror(char error_text[]);
float *fvector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
  long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
  long ndl, long ndh);

#endif /* _NR_UTILS_H_ */

#ifndef _NR_FUNCS_H_
#define _NR_FUNCS_H_

#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

void tred2(float **a, int n, float d[], float e[]);
int svdcmp(float **a, int m, int n, float w[], float **v);
int tqli(float d[], float e[], int n, float **z);
float pythag(float a, float b);

#endif /* _NR_FUNCS_H_ */

