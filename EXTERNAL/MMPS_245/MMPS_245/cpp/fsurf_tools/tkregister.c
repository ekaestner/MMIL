/*============================================================================
 Copyright (c) 1996 Martin Sereno and Anders Dale
=============================================================================*/
/*   $Id: tkregister.c,v 1.10 2001/05/08 01:18:41 fischl Exp $   */
/* last modified by Don Hagler 2005/03/08 */


#ifndef lint
/*static char vcid[] = "$Id: tkregister.c,v 1.10 2001/05/08 01:18:41 fischl Exp $";*/
static char vcid[] = "$Id: tkregister.c,v 1.11 2005/03/08 17:53:50 dhagler Exp $";
#endif /* lint */

#define TCL
#if 0
#include <tcl.h>
#include <tk.h>
#else
#include "tcl-7-4.h"
#include "tk-4-0.h"
#endif
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include "surflib.h"

#if 0
#define TCL8
#endif
#ifndef OPENGL
#define OPENGL 
#endif

#ifndef RGB
#define RGB
#endif

#ifndef TCL
#define TCL
#endif

#ifdef TCL
#  define PR   {if(promptflag){fputs("% ", stdout);} fflush(stdout);}
#else
#  define PR
#endif

#ifdef OPENGL
#  include "xwindow.h"
#  include <X11/keysym.h>
#ifdef RGB
#  define swapbuffers()    pseudo_swapbuffers(); tkoSwapBuffers()
#else
#  define swapbuffers()    tkoSwapBuffers()
#endif
#  define frontbuffer(X)   glDrawBuffer(GL_FRONT)
#  define backbuffer(X)    glDrawBuffer(GL_BACK)
#  define clear()          glClear(GL_COLOR_BUFFER_BIT)
#  define getorigin(X,Y)   *(X) = w.x; *(Y) = 1024 - w.y - w.h /*WINDOW_REC w;*/
#  define getsize(X,Y)     *(X) = w.w; *(Y) = w.h
#  define Colorindex       unsigned short
#  define color(X)         glIndexs(X)
#  define mapcolor(I,R,G,B)  \
            tkoSetOneColor((int)I,(float)R/255.0,(float)G/255.0,(float)B/255.0)
#ifdef RGB
#  define rectwrite(X0,Y0,X1,Y1,P) \
            glRasterPos2i(X0,Y0); \
            glDrawPixels(X1+1,Y1+1,GL_RGB,GL_UNSIGNED_BYTE,P)
#else
#  define rectwrite(X0,Y0,X1,Y1,P) \
            glRasterPos2i(X0,Y0); \
            glDrawPixels(X1+1,Y1+1,GL_COLOR_INDEX,GL_UNSIGNED_SHORT,P)
#endif
#else
#  include <gl.h>
#  include <device.h>
#endif

/* Prototypes */

int do_one_gl_event(Tcl_Interp *interp);
void open_window(char name[]);
void  resize_window_intstep();
void move_window(int x, int y);
void reload_buffers();
void blinkbuffers();
void record_swapbuffers();
void resize_buffers(int x, int y);
void read_reg(char fname[]);
void write_reg(char fname[]);
void write_default_reg(char fname[]);
void cor2anareg(char fpref[],char ana[],char reg[]);
void set_session_subject(char *str);
void make_backup(char fname[]);
void save_rgb(char fname[]);
void scrsave_to_rgb(char fname[]);
void pix_to_rgb(char *fname);
void downslice();
void upslice();
void goto_point(char dir[]);
void write_point(char dir[]);
void rotate_brain(float a, char c);
void align_points();
void translate_brain(float a,char c);
void scale_brain(float s, char c);
void mirror_brain();
void set_cursor(float xpt, float ypt, float zpt);
void set_scale();
void redraw();
void pop_gl_window();
void mri2pix(float xpt, float ypt, float zpt, int *jpt, int *ipt, int *impt);
int imval(float px,float py,float pz);
float Error(int p,float dp);
void optimize(int maxiter);
void optimize2();
void read_images(char fpref[]);
void read_second_images(char fpref[]);
void select_pixel(short sx,short sy);
void transform(float x1,float y1,float z1,float *x2,float *y2,float *z2,float M[4][4]);
void draw_image(int imc,int ic,int jc);
void blur(float factor);
void make_filenames(char *lsubjectsdir);
void read_float_images(float ***fim,char *format,int nslices,int nperslice,int xdim,int ydim,short **buf);
void usecnap(int usec);
void initcolormap();
void pseudo_swapbuffers();
void read_AFNI_HEAD(char *fname, float *ffact, int *btype, char *byteorder);

#if 0
#ifdef USE_mriTransform

#include "mriTransformRef.h"

static mriTransformRef gRASToFuncIdxTransform = NULL;
static void InitFunctionalTransform ( char* isRegisterFile );
static int  ConvertRASToFuncIdx  ( float  ifRASX, 
           float  ifRASY, 
           float  ifRASZ,
           float* ofFuncIdxX, 
           float* ofFuncIdxY, 
           float* ofFuncIdxZ );

#endif /* USE_mriTransform */
#endif

#ifndef TRUE
#  define TRUE 1
#endif
#ifndef FALSE
#  define FALSE 0
#endif

#define APD1    1   /* imtypes */
#define BSHORT  2
#define SKIP    3
#define AFNI    4

#ifndef SQR
#define SQR(x)       ((x)*(x))
#endif
#define MATCH(A,B)   (!strcmp(A,B))
#define MATCH_STR(S) (!strcmp(str,S))

#define HEADERSIZE 80
#define NUMVALS 256
#define MAXIM 256
#define MAXPTS 10000
#define MAXPARS 10
#ifdef RGB
#define MAPOFFSET 0
#else
#define MAPOFFSET 1000
#endif
#define CORONAL 0
#define HORIZONTAL 1
#define SAGITTAL 2
#define NAME_LENGTH STRLEN
#define BLINK_DELAY 30
#define BLINK_TIME 20
#define MOTIF_XFUDGE   8
#define MOTIF_YFUDGE  32
#define TMP_DIR "tmp"
#define SLICE1_POS_UNDEF   9999.0
  /* overlay modes */
#define TARGET   1
#define MOVEABLE 2

#define BRIK_BYTE    0
#define BRIK_SHORT   1
#define BRIK_LONG    2
#define BRIK_FLOAT   3
#define BRIK_DOUBLE  4
#define BRIK_COMPLEX 5

int plane = CORONAL;
int xnum=256,ynum=256;
int ptype;
float ps,st,xx0,xx1,yy0,yy1,zz0,zz1;
int zf,ozf;
float fsf;
int xdim,ydim;
unsigned long bufsize, bufsize_2;
unsigned char *buf, *buf_2;
unsigned char **im[MAXIM];
unsigned char **sim[6]; 
int changed[MAXIM];
#ifdef RGB
GLubyte *vidbuf;
GLubyte *blinkbuft;
GLubyte *blinkbufm;
#else
Colorindex *vidbuf;
#endif
unsigned char *binbuff;
int imnr0,imnr1,numimg;
int wx0=600,wy0=100;  /* (114,302) (117,90) */
int ptsflag = FALSE;
int maxflag = FALSE;
int updateflag = FALSE;
int blinkflag = FALSE;
int blinkdelay = BLINK_DELAY;
int blinktime = BLINK_TIME;
int overlay_mode = TARGET;
int visible_mode = 0;
int last_visible_mode = 0;
int visible_plane = 0;
int last_visible_plane = 0;
int editedmatrix = FALSE;
int maskflag = FALSE;
int scrsaveflag = TRUE;
int openglwindowflag = FALSE;
int promptflag = FALSE;
int followglwinflag = TRUE;
int initpositiondoneflag = FALSE;
int npts = 0;
int prad = 0;
float TM[4][4];
float tm[4][4];
double ps_2,st_2,fscale_2; /* was float */
float xx0_2,xx1_2,yy0_2,yy1_2,zz0_2,zz1_2;
int xnum_2,ynum_2,numimg_2;
int imnr0_2,imnr1_2,xdim_2=0,ydim_2=0;
float **fim_2[MAXIM];
float ptx[MAXPTS],pty[MAXPTS],ptz[MAXPTS];
float par[MAXPARS],dpar[MAXPARS];
int nslices=0,nperslice=0;

static char *progname = NULL;
double fthresh = 0.35;
double fsquash = 12.0;
double fscale = 255;
int imc=0,ic=0,jc=0;
float xc=0,yc=0,zc=0;
float xc_old=0,yc_old=0,zc_old=0;
int impt = -1,ipt = -1,jpt = -1;

char *subjectsdir;   /* SUBJECTS_DIR */
char *srname;        /* sessiondir (funct: image) */
char *psrname;       /* parent sessiondir (funct: 970703MS) */
char *pname;         /* subject */
char *regfname;      /* register.dat */
char *afname;        /* analyse.dat */
char *targpref;      /* abs single image structural stem name */
char *movformat;     /* abs single image epi structural stem name */
char *tfname;        /* (dir!) subjectsdir/name/tmp/ */
char *sgfname;       /* (dir!) set: get from cwd: $session/rgb/ */

int blinktop = 0;    /* says whats on top while blinking */
int invalid_buffers = 1;

int pswapnext = TARGET;
char colormap[512];

#ifdef TCL
  int Register(ClientData clientData,Tcl_Interp *interp, int argc, char *argv[])
#else
  main(int argc,char **argv)
#endif
{
    char *fbuf;
    char *getenv(),*lsubjectsdir,*ext;
    char targdir[NAME_LENGTH], reltargpref[NAME_LENGTH];
    char movdir[NAME_LENGTH], movpref[NAME_LENGTH], relmovformat[NAME_LENGTH];
    FILE *fp;
    /* short dev, val; */
    /* short sx,sy; */ /* Screencoord sx,sy; */
    int i; /* ,j,k,l; */
    int regdatflag, abstargflag, skipanadatflag, subj2subjflag;
    /* long lxdim,lydim; */
    float fslice1;
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

    if (argc<2) {
      printf("%s: version %s\n", progname, vcid) ;
      printf("\n");
      printf("  Usage:\n");
      printf("    %s local                      [ mov=cwd        targ=impliedCOR ]\n",argv[0]);
      printf("    %s local <abstargdir>         [ mov=cwd        targ=absCORdir  ]\n",argv[0]);
      printf("    %s <absmovdir> <abstargdir>   [ mov=abs        targ=absCORdir  ]\n",argv[0]);
      printf("    %s <subj> <movdir> <targdir>  [ mov=relCORdir  targ=relCORdir  ]\n",argv[0]);
      printf("    %s -mkdefault\n",argv[0]);
      printf("    %s -readinfo\n",argv[0]);
      printf("      [usu. run in session subdir w/ISIR, functEPI, or alignCOR scan]\n");
      printf("      [reldirs are in $subject/mri, composite path OK]\n");
      printf("      [abstargdirs allow cross-subject reg]\n");
      printf("\n");
      printf("  Args:\n");
      printf("    local       movdir is cwd\n");
      printf("    -mkdefault  make default register.dat in cwd\n");
      printf("    -readinfo   register.dat,analyse.dat from COR-.info (COR=>COR)\n");
      printf("    <subject>   in subjects database\n");
      printf("    <movdir>    dir w/images to register (im,bshort,BRIK,brik,COR)\n");
      printf("    <targdir>   dir w/target image (COR only)\n");
      printf("\n");
      printf("  Files:\n");
      printf("    target image    COR images--subj database or session align scan\n");
      printf("    mov image       EPI,alignCOR,subjCOR--described in analyse.dat\n");
      printf("    register.dat    movdir: subj,pixsize,slicethick,transform\n");
      printf("    analyse.dat     movdir: aligndir,imformat,nslices,x/ysize\n");
      printf("    edit.dat        subj/tmp/edit.dat--for SEND/GOTO\n");
      printf("\n");
      printf("  Options: (at end, = movdir analyse.dat, required for abs => abs)\n");
      printf("    -movformat <str>   (-m) mov img format (relative printf)\n");
      printf("    -nslices <n>       (-n) num slices moveable\n");
      printf("    -reps <n>          (-e) images per slice (need BRIK offsets)\n");
      printf("    -imres <x> <y>     (-i) xsize,ysize moveable images\n");
      printf("    -regdat <file>     (-r) else default=./register.dat\n");
      printf("    -fslice1 <mm>      (-1) funct/sess offset--matrix for align\n");
      printf("\n");
      printf("  Examples:\n");
      printf("    tkregister local   (struct<=local; mk analyse.dat,register.dat)\n");
      printf("    tkregister local /fspace/990402TS/mprage  (targ must be COR)\n");
      printf("    tkregister /fspace/990402TS/mprage /subjects/test/mri/T1\n");
      printf("    tkregister /fspace/990402TS/eccen  /subjects/test/mri/T1\n");
      printf("    tkregister /fspace/990402TS/eccen  /subjects/test/mri/mprage\n");
      printf("    tkregister marty orig/im1 orig/im2   (struct2 <= struct1)\n");
      printf("\n");
#ifdef OPENGL
      printf("                                         [vers: 13c07b15--OpenGL]\n");
#else
      printf("                                         [vers: 13c07b15--GL]\n");
#endif
      exit(1);
    }
    
    lsubjectsdir = getenv("SUBJECTS_DIR");
    if (lsubjectsdir==NULL) {
      printf("register: env var SUBJECTS_DIR undefined (use setenv)\n");
      exit(1);
    }
    zf = ozf = 2;
    fsf = (float)zf;
    fslice1 = SLICE1_POS_UNDEF;
    abstargflag = FALSE;
    regdatflag = FALSE;
    skipanadatflag = FALSE;
    subj2subjflag = FALSE;

    /* default targetdir, subject, session, datfiles */
    make_filenames(lsubjectsdir);
    strcpy(movdir,".");  /* . never used */
    strcpy(targdir,"T1");
    strcpy(pname,"nobody");
    getcwd(srname, NAME_LENGTH);
    strcpy(regfname,"register.dat");
    strcpy(afname,"analyse.dat");
    strcpy(reltargpref,"COR-");   /* moveable COR only */
    strcpy(sgfname,"/tmp/tkregister.rgb");

    for (i=0;i<argc;i++) {
      if (MATCH(argv[i],"-mkdefault")) {
        write_default_reg(regfname); exit(0);
      }
      else if (MATCH(argv[i],"-readinfo")) {
        cor2anareg(reltargpref,afname,regfname); exit(0);
      }
      else if ((MATCH(argv[i],"-fslice1")||MATCH(argv[i],"-1")) && i+1<argc) {
        fslice1 = atof(argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-movformat")||MATCH(argv[i],"-m")) && i+1<argc) {
        strcpy(relmovformat,argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-nslices") || MATCH(argv[i],"-n")) && i+1<argc) {
        nslices = atoi(argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-reps") || MATCH(argv[i],"-e")) && i+1<argc) {
        nperslice = atoi(argv[i+1]); i++;
      }
      else if ((MATCH(argv[i],"-imres") || MATCH(argv[i],"-i")) && i+2<argc) {
        xdim_2 = atoi(argv[i+1]); ydim_2 = atoi(argv[i+2]); i+=2;
      }
      else if ((MATCH(argv[i],"-regdat") || MATCH(argv[i],"-r")) && i+1<argc) {
        strcpy(regfname,argv[i+1]); regdatflag = TRUE; i++;
      }
      else if (argv[i][0]=='-') {
        printf("fourier: ### unrecognized option: %s\n",argv[i]); exit(1);
      }
    }
    if (!MATCH(relmovformat,"") && nslices && nperslice && xdim_2 && ydim_2) {
      skipanadatflag = TRUE;
    }

    /*** local funct or COR => implied subj COR (name from regdat) ***/
    if (argc==2) {
      if (MATCH(argv[1],"local"))
        ; /* defaults (cmdline flags can override) */
      else if (MATCH(argv[1],"1") || MATCH(argv[1],"2") || MATCH(argv[1],"3")) {
        zf=atoi(argv[1]); fsf=(float)zf; }
      else { printf("register: ### bad arg: %s\n",argv[1]); exit(1); }
      printf("register: cwd scan => subject's %s/COR's\n",targdir);
    }
    /*** local funct or COR => absolute struct COR (overrides regdat name) ***/
    else if (argc>2 && MATCH(argv[1],"local") && argv[2][0]!='-') {
      if (argv[2][0]!='/') {printf("register: ### need abs path\n");exit(1);}
      abstargflag = TRUE;
      strcpy(targdir,argv[2]);
      printf("register: absolute targ dir (overrides subj in register.dat)\n");
      printf("register: cwd scan =>\n");
      printf("register:   %s/COR's\n",targdir);
    }
    /*** absolute funct or COR => absolute struct COR (requires all opts) ***/
    else if (argc>2 && argv[1][0]=='/' && argv[2][0]=='/') {
      if (!regdatflag) {printf("register: ### must specify -regdat\n");exit(1);}
      if (!skipanadatflag) {printf("register: ### missing options\n");exit(1);}
      abstargflag = TRUE;
      strcpy(movdir,argv[1]);
      strcpy(targdir,argv[2]);
      set_session_subject(movdir);  /* name from namefile */
      sprintf(movformat,"%s/%s",movdir,relmovformat);
      printf("register: absolute mov, targ dirs (ignores subj in regdat)\n");
    }
    /*** subj COR1 => subj COR2, CORdir's relative to $subjname/mri */
    else if (argc>3 && !(MATCH(argv[1],"local")) &&
             argv[1][0]!='-' && argv[2][0]!='-' && argv[3][0]!='-') {
      subj2subjflag = TRUE;
      strcpy(pname,argv[1]);
      strcpy(movdir,argv[2]);
      strcpy(targdir,argv[3]);
      sprintf(movpref,"%s/%s/mri/%s/COR-",lsubjectsdir,pname,movdir);
      if (!regdatflag)  /* expand abbrev */
        sprintf(regfname,"%s/%s/mri/%s/register.dat",lsubjectsdir,pname,movdir);
      if (skipanadatflag)
        sprintf(
          movformat,"%s/%s/mri/%s/%s",lsubjectsdir,pname,movdir,relmovformat);
      else {  /* expand abbrevs (movformat later read from created anadat) */
        sprintf(afname,"%s/%s/mri/%s/analyse.dat",lsubjectsdir,pname,movdir);
        sprintf(movpref,"%s/%s/mri/%s/COR-",lsubjectsdir,pname,movdir);
        cor2anareg(movpref,afname,regfname); /* created ana: relmovformat! */
      }
      printf("register: %s/COR's => %s/COR's\n",movdir,targdir);
    }
    else {  /* local plus options */
      if (!(MATCH(argv[1],"local"))) {
        printf("register: ### bad arg: %s\n",argv[1]);
        exit(1);
      }
    }
    
    /* to-be-registered xform (biff subject name only if still "nobody") */
    read_reg(regfname);   /* subj, pixsize, slicethk */

    /* read target images (COR format only), expand abbrev if needed */
    if (abstargflag)
      sprintf(targpref,"%s/COR-",targdir);  /* override register.dat */
    else
      sprintf(targpref,"%s/%s/mri/%s/COR-",lsubjectsdir,pname,targdir);
    read_images(targpref);

    /* analyse.dat for moveable images: format,nslices,nper,xydim */
    if (!skipanadatflag) {
      fp = fopen(afname,"r");
      if (fp==NULL){printf("register: ### File %s not found\n",afname);exit(1);}
      fscanf(fp,"%*s");
      fscanf(fp,"%s",relmovformat);
      fscanf(fp,"%d",&nslices);    /* hires: overlaps COR-.info */
      fscanf(fp,"%d",&nperslice);
      fscanf(fp,"%d",&xdim_2);     /* hires: overlaps COR-.info */
      fscanf(fp,"%d",&ydim_2);     /* hires: overlaps COR-.info */
      fclose(fp);
      /* reexpand overwritten subj2subj movformat for test&read below */
      if (subj2subjflag)
        sprintf(
           movformat,"%s/%s/mri/%s/%s",lsubjectsdir,pname,movdir,relmovformat);
      else
        strcpy(movformat,relmovformat);  /* local */
    }

    /* read moveable images: (1) COR vs. (2) im,bshort,skip,BRIK,brik */
    fbuf = NULL;
    ext = strstr(movformat,"COR-%03d");
    if (ext!=NULL&&strlen(ext)==strlen("COR-%03d")) { /* if moveable is COR */
      strcpy(movpref,movformat);
      movpref[strlen(movformat)-4] = '\0';  /* trim %03d */
      read_second_images(movpref);
    } else {   /* read short images using full format */
      /* huhu strange typecast necessary here.... */
      read_float_images(fim_2,movformat,nslices,nperslice,xdim_2,ydim_2,(short **)&fbuf);
      printf("register: done reading to-be-registered images\n");
      imnr0_2 = 0;
      imnr1_2 = nslices-1;
      if (fslice1 == SLICE1_POS_UNDEF) {  /* default: center z BB */
        printf("register: no slice offset (default centered z bounding box)\n");
        yy0_2 = -st_2*nslices/2.0;
        yy1_2 = st_2*nslices/2.0;
      }
      else {                              /* funct/session z offset */
        yy0_2 = -fslice1;
        yy1_2 = yy0_2 + st_2*nslices;
      }
      xx0_2 = -ps_2*xdim_2/2.0;
      xx1_2 = ps_2*xdim_2/2.0;
      zz0_2 = -ps_2*ydim_2/2.0;
      zz1_2 = ps_2*ydim_2/2.0;
    }

    sprintf(tfname,"%s/%s/%s",lsubjectsdir,pname,TMP_DIR);
#ifdef RGB
    vidbuf = (GLubyte *)lcalloc(3*bufsize*SQR(zf),sizeof(GLubyte));
    blinkbuft = (GLubyte *)lcalloc(3*bufsize*SQR(zf),sizeof(GLubyte));
    blinkbufm = (GLubyte *)lcalloc(3*bufsize*SQR(zf),sizeof(GLubyte));
#else
    vidbuf = (Colorindex *)lcalloc(bufsize*SQR(zf),sizeof(Colorindex));
#endif
    binbuff = (unsigned char *)calloc(3*xdim*ydim,sizeof(char));

    open_window(pname);

    set_scale();

    imc = zf*imnr1/2;
    ic = ydim/2;
    jc = xdim/2;
 
    /* load both buffers; sets visible_mode/plane, last_visible_mode/plane */
    overlay_mode = MOVEABLE;
    redraw();
    overlay_mode = TARGET;
    redraw();
    updateflag = FALSE;

#ifndef TCL  /* tcl: omit event loop */
#  ifdef OPENGL
    printf("register: ### non-tk register event loop not written for OpenGL\n");
    exit(1);
#  else
    while(1) {
        dev = qread(&val);  /* blocks here waiting for next event */
        switch(dev)  {
    
          case ESCKEY:     /* kill on escape */
            exit(0);
            break;
          case LEFTMOUSE:
            if (val == 0)  break;
            xc_old=xc,yc_old=yc,zc_old=zc;
            sx = getvaluator(MOUSEX);
            sy = getvaluator(MOUSEY);
            select_pixel(sx,sy);
            if (visible_mode!=overlay_mode)
              overlay_mode = visible_mode;
            updateflag = TRUE;
            break;
          case MIDDLEMOUSE:
            if (val == 0)  break;
            while (getbutton(MIDDLEMOUSE))
            {
              sx = getvaluator(MOUSEX);
              sy = getvaluator(MOUSEY);
              select_pixel(sx,sy);
              for (i= -prad;i<=prad;i++)
              for (j= -prad;j<=prad;j++)
              {
                im[imc/zf][(ydim-1-ic+i)/zf][jc/zf+j] = 255;
              }
              changed[imc/zf] = TRUE;
            }
            updateflag = TRUE;
            break;
          case RIGHTMOUSE:
            if (val == 0)  break;
            while (getbutton(RIGHTMOUSE))
            {
              sx = getvaluator(MOUSEX);
              sy = getvaluator(MOUSEY);
              select_pixel(sx,sy);
              for (i= -prad;i<=prad;i++)
              for (j= -prad;j<=prad;j++)
              {
                im[imc/zf][(ydim-1-ic+i)/zf][jc/zf+j] = 0;
              }
              changed[imc/zf] = TRUE;
            }
            updateflag = TRUE;
            break;
          case REDRAW:         /* make the window re-sizable */
            reshapeviewport();
            getsize(&lxdim,&lydim);
            xdim = (int)lxdim;
            ydim = (int)lydim;
            resize_buffers(xdim,ydim);
            zf = xdim/xnum;
            fsf = (float)zf;
            imc = (zf*imc)/ozf;
            ic = (zf*ic)/ozf;
            jc = (zf*jc)/ozf;
            redraw();
            break;
          case UPARROWKEY:
            if (val == 0)  break;
      fscale_2 *= 1.5;
            printf("fscale_2 = %f\n",fscale_2);
            updateflag = TRUE;
            break;
          case DOWNARROWKEY:
            if (val == 0)  break;
      fscale_2 /= 1.5;
            printf("fscale_2 = %f\n",fscale_2);
            updateflag = TRUE;
            break;
          case LEFTARROWKEY:
            if (val == 0)  break;
            if (plane==CORONAL)
        imc = (imc<zf)?imnr1*zf-zf+imc:imc-zf;
            else
            if (plane==HORIZONTAL)
        ic = (ic<zf)?ydim-zf+ic:ic-zf;
            else
            if (plane==SAGITTAL)
        jc = (jc<zf)?xdim-zf+jc:jc-zf;
      updateflag = TRUE;
            break;
          case RIGHTARROWKEY:
            if (val == 0)  break;
            if (plane==CORONAL)
        imc = (imc>=imnr1*zf-zf)?imc+zf-imnr1*zf:imc+zf;
            else
            if (plane==HORIZONTAL)
        ic = (ic>=ydim-zf)?ic+zf-ydim:ic+zf;
            else
            if (plane==SAGITTAL)
        jc = (jc>=xdim-zf)?jc+zf-xdim:jc+zf;
      updateflag = TRUE;
            break;
          case INSERTKEY:
            if (val == 0)  break;
            scale_brain(1/1.02,'x');
            overlay_mode = MOVEABLE;
            updateflag = TRUE;
            break;
          case DELKEY:
            if (val == 0)  break;
            scale_brain(1.02,'x');
            overlay_mode = MOVEABLE;
            updateflag = TRUE;
            break;
          case HOMEKEY:
            if (val == 0)  break;
            scale_brain(1/1.02,'y');
            overlay_mode = MOVEABLE;
            updateflag = TRUE;
            break;
          case ENDKEY:
            if (val == 0)  break;
            scale_brain(1.02,'y');
            overlay_mode = MOVEABLE;
            updateflag = TRUE;
            break;
          case PAGEUPKEY:
            if (val == 0)  break;
            scale_brain(1/1.02,'z');
            overlay_mode = MOVEABLE;
            updateflag = TRUE;
            break;
          case PAGEDOWNKEY:
            if (val == 0)  break;
            scale_brain(1.02,'z');
            overlay_mode = MOVEABLE;
            updateflag = TRUE;
            break;
          case KEYBD:
            switch((char)val)
            {
              case 'A': align_points();
                        overlay_mode = MOVEABLE;
                        updateflag = TRUE;
                        break;
              case 't': case 'T':
                        write_reg(regfname);
                        break;
              case 'x': plane=SAGITTAL;updateflag = TRUE;break;
              case 'y': plane=HORIZONTAL;updateflag = TRUE;break;
              case 'z': plane=CORONAL;updateflag = TRUE;break;
              case 'b': prad = 0;break;
              case 'B': prad = 1;break;
              case 'o': optimize(10);updateflag = TRUE;break;
              case 'O': optimize2();updateflag = TRUE;break;
              case 'M': mirror_brain();updateflag = TRUE;break;
              case 'm': maxflag = !maxflag;
                        printf("maxflag=%d",maxflag);
                        updateflag = TRUE;
                        break;
              case '*': fsquash *= 1.1; break;
              case '/': fsquash /= 1.1; break;
              case '!': fsquash = 5.0; break;
              case '+': fthresh += 0.05;
                        break;
              case '-': fthresh -= 0.05;
                        break;
              case '0': record_swapbuffers();
                        break;
              case '1': overlay_mode = TARGET;
      updateflag = TRUE;
                        break;
              case '2': overlay_mode = MOVEABLE;
      updateflag = TRUE;
                        break;
              case '[': if (plane==SAGITTAL) rotate_brain(-5.0,'x'); 
      if (plane==HORIZONTAL) rotate_brain(5.0,'z'); 
      if (plane==CORONAL) rotate_brain(-5.0,'y'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case ']': if (plane==SAGITTAL) rotate_brain(5.0,'x'); 
      if (plane==HORIZONTAL) rotate_brain(-5.0,'z'); 
      if (plane==CORONAL) rotate_brain(5.0,'y'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case '{': if (plane==SAGITTAL) rotate_brain(-50.0,'x'); 
      if (plane==HORIZONTAL) rotate_brain(50.0,'z'); 
      if (plane==CORONAL) rotate_brain(-50.0,'y'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case '}': if (plane==SAGITTAL) rotate_brain(50.0,'x'); 
      if (plane==HORIZONTAL) rotate_brain(-50.0,'z'); 
      if (plane==CORONAL) rotate_brain(50.0,'y'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case 'p': if (plane==SAGITTAL) translate_brain(-0.5,'z');
      if (plane==HORIZONTAL) translate_brain(-0.5,'y'); 
      if (plane==CORONAL) translate_brain(-0.5,'z'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case '.': if (plane==SAGITTAL) translate_brain(0.5,'z'); 
      if (plane==HORIZONTAL) translate_brain(0.5,'y'); 
      if (plane==CORONAL) translate_brain(0.5,'z'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case 'l': if (plane==SAGITTAL) translate_brain(0.5,'y'); 
      if (plane==HORIZONTAL) translate_brain(-0.5,'x'); 
      if (plane==CORONAL) translate_brain(-0.5,'x'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case ';': if (plane==SAGITTAL) translate_brain(-0.5,'y'); 
      if (plane==HORIZONTAL) translate_brain(0.5,'x'); 
      if (plane==CORONAL) translate_brain(0.5,'x'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case 'P': if (plane==SAGITTAL) translate_brain(-2.5,'z');
      if (plane==HORIZONTAL) translate_brain(-2.5,'y'); 
      if (plane==CORONAL) translate_brain(-2.5,'z'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case '>': if (plane==SAGITTAL) translate_brain(2.5,'z'); 
      if (plane==HORIZONTAL) translate_brain(2.5,'y'); 
      if (plane==CORONAL) translate_brain(2.5,'z'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case 'L': if (plane==SAGITTAL) translate_brain(2.5,'y'); 
      if (plane==HORIZONTAL) translate_brain(-2.5,'x'); 
      if (plane==CORONAL) translate_brain(-2.5,'x'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case ':': if (plane==SAGITTAL) translate_brain(-2.5,'y'); 
      if (plane==HORIZONTAL) translate_brain(2.5,'x'); 
      if (plane==CORONAL) translate_brain(2.5,'x'); 
      updateflag = TRUE;
                        overlay_mode = MOVEABLE;
                        break;
              case 'r': goto_point(tfname);
                        overlay_mode = MOVEABLE;
                        updateflag = TRUE;
                        break;
            }
            break;
          }  /* dev switch */
    set_scale();
          if (!qtest() && updateflag)  /* tcl: put in do_one_gl_event */
          {
            redraw();
            updateflag = FALSE;
          }
    }  /* while events */
#  endif
#endif  /* tcl: omit event loop */
    return GL_FALSE;
}

int do_one_gl_event(Tcl_Interp *interp)   /* tcl */
{
#ifdef OPENGL     /* derived from event.c:DoNextEvent(),tkExec() */
  XEvent current, ahead;
  char buf[1000];
  char command[NAME_LENGTH];
  KeySym ks;
  static int ctrlkeypressed = FALSE;
  static int altkeypressed = FALSE;
  static int shiftkeypressed = FALSE;
  static int button1pressed = FALSE;
  static int button2pressed = FALSE;
  static int button3pressed = FALSE;
  short sx,sy; /* Screencoord sx,sy; */
  XWindowAttributes wat;
  Window junkwin;
  int rx, ry;

  blinkbuffers();

  if (updateflag) {
    set_scale();
    redraw();
    updateflag = FALSE;
  }

  if (XPending(xDisplay)) {  /* do one if queue test */
    XNextEvent(xDisplay, &current);   /* blocks here if no event */
    switch (current.type) {
      case ConfigureNotify:
        XGetWindowAttributes(xDisplay, w.wMain, &wat);
        XTranslateCoordinates(xDisplay, w.wMain, wat.root,
                              -wat.border_width, -wat.border_width,
                              &rx, &ry, &junkwin);
        w.x = rx;  /* left orig         (current.xconfigure.x = 0 relative!) */
        w.y = ry;  /* top orig: +y down (current.xconfigure.y = 0 relative!) */
        w.w = current.xconfigure.width;
        w.h = current.xconfigure.height;

        resize_window_intstep();
#if 0
/* this causes core dump for some reason
   omitting it seems to cause no problems -- dhagler */
printf("debug: Tcl_Eval(interp,\"set zf $zf\")\n");
        Tcl_Eval(interp,"set zf $zf");  /* touch for trace */ 
#endif
        if (followglwinflag && zf != 4) {
          /* below */ 
          /*sprintf(command,"wm geometry . +%d+%d",w.x,w.y+w.h+MOTIF_YFUDGE);*/
          /* above */ 
          sprintf(command,"putaboveglwin %d %d", w.x, w.y-MOTIF_YFUDGE);
          Tcl_Eval(interp,command);
          /* Tcl_Eval(interp,"raise ."); */
        }
        if (visible_mode!=overlay_mode) /* suppress mode change */
          overlay_mode = visible_mode;
        glViewport(0, 0, xdim, ydim);
        /* updateflag = TRUE; */
        break;

      case Expose:
        if (XPending(xDisplay)) {
          XPeekEvent(xDisplay, &ahead);
          if (ahead.type==Expose) break;  /* skip extras */
        }
        if (visible_mode!=overlay_mode)
          overlay_mode = visible_mode;
        updateflag = TRUE;
        break;

      case ButtonPress:
        sx = current.xbutton.x;
        sy = current.xbutton.y;
        sx += w.x;   /* convert back to screen pos (ugh) */
        sy = 1024 - w.y - sy;
        if (current.xbutton.button == 1) {  /** left **/
          select_pixel(sx,sy);
#if 0
/* this causes core dump for some reason
   omitting it seems to cause no problems -- dhagler */
printf("debug: Tcl_Eval(interp, \"unzoomcoords $plane\")\n");
          Tcl_Eval(interp,"unzoomcoords $plane"); 
#endif
          button1pressed = TRUE;
        }
        if (current.xbutton.button == 2) {  /** middle **/
          button2pressed = TRUE;
        }
        if (current.xbutton.button == 3) {  /** right **/
          button3pressed = TRUE;
        }
        if (visible_mode!=overlay_mode)
          overlay_mode = visible_mode;
        updateflag = TRUE;
        break;

      case ButtonRelease:
        if (current.xbutton.button == 1)  button1pressed = FALSE;
        if (current.xbutton.button == 2)  button2pressed = FALSE;
        if (current.xbutton.button == 3)  button3pressed = FALSE;
        break;

      case KeyPress:
        XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
        switch (ks) {

          /* numbers */
          case XK_0: record_swapbuffers(); break;
          case XK_1: overlay_mode = TARGET;   updateflag = TRUE; break;
          case XK_2: overlay_mode = MOVEABLE; updateflag = TRUE; break;

          /* upper case */
          case XK_A: if (altkeypressed)  ; break;

          /* lower case */
          case XK_x: plane=SAGITTAL;   updateflag = TRUE; break;
          case XK_y: plane=HORIZONTAL; updateflag = TRUE; break;
          case XK_z: plane=CORONAL;    updateflag = TRUE; break;

          /* others */
          case XK_Up:
            
            Tcl_Eval(interp,
                  "set fscale_2 [expr $fscale_2 * 1.5]; set updateflag TRUE");
            break;
          case XK_Down:
             Tcl_Eval(interp,
                  "set fscale_2 [expr $fscale_2 / 1.5]; set updateflag TRUE");
            break;
          case XK_Right:
            Tcl_Eval(interp,"changeslice up $plane; set updateflag TRUE");
            break;
          case XK_Left:
            Tcl_Eval(interp,"changeslice down $plane; set updateflag TRUE");
            break;

          /* modifiers */
          case XK_Shift_L:   case XK_Shift_R:   shiftkeypressed=TRUE; break;
          case XK_Control_L: case XK_Control_R: ctrlkeypressed=TRUE;  break;
          case XK_Alt_L:     case XK_Alt_R:     altkeypressed=TRUE;   break;
        }
        break;

      case KeyRelease:   /* added this mask to xwindow.c */
        XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
        switch (ks) {
          case XK_Shift_L:   case XK_Shift_R:   shiftkeypressed=FALSE; break;
          case XK_Control_L: case XK_Control_R: ctrlkeypressed=FALSE;  break;
          case XK_Alt_L:     case XK_Alt_R:     altkeypressed=FALSE;   break;
    case XK_0:
      break;
        }
        break;
    }
    return GL_FALSE;
  }

#else  /* use gl calls */
  short dev, val;
  static int ctrlkeypressed = FALSE;
  static int altkeypressed = FALSE;
  Screencoord sx,sy;
  long lxdim,lydim;
  int i,j,k;

  blinkbuffers();

  if (updateflag) {
    set_scale();
    redraw();
    updateflag = FALSE;
  }

  if (qtest()) {  /* do one event */
    dev = qread(&val);  /* would block here if no event */
    if (dev != LEFTMOUSE && dev != RIGHTMOUSE) /* hack: mouse zeros getbutton!*/
      ctrlkeypressed = getbutton(LEFTCTRLKEY) || getbutton(RIGHTCTRLKEY);
    switch(dev) {
      case REDRAW:
        if (visible_mode!=overlay_mode)  /* suppress mode change on winrepair */
          overlay_mode = visible_mode;
        updateflag = TRUE;
        reshapeviewport();
        getsize(&lxdim,&lydim);
        xdim = (int)lxdim;
        ydim = (int)lydim;
        resize_buffers(xdim,ydim);
        ozf = zf;
        zf = xdim/xnum;
        fsf = (float)zf;
        imc = (zf*imc)/ozf;
        ic = (zf*ic)/ozf;
        jc = (zf*jc)/ozf;
        Tcl_Eval(interp,"set zf $zf");  /* touch for trace */
        if (followglwinflag && zf != 4) {
          sprintf(command,"wm geometry . +%d+%d",
                      xorig, 1024 - yorig + MOTIF_YFUDGE /*+ MOTIF_XFUDGE*/);
          Tcl_Eval(interp,command);
          Tcl_Eval(interp,"raise .");
        }
        break;
      case LEFTARROWKEY:
        if (!altkeypressed) break;
        if (val == 0)  break;
        downslice();
        updateflag = TRUE;
        break;
      case RIGHTARROWKEY:
        if (!altkeypressed) break;
        if (val == 0)  break;
        upslice();
        updateflag = TRUE;
        break;
      case LEFTMOUSE:
        if (val == 0)  break;
        xc_old=xc,yc_old=yc,zc_old=zc;
        sx = getvaluator(MOUSEX);
        sy = getvaluator(MOUSEY);
        select_pixel(sx,sy);
        Tcl_Eval(interp,"unzoomcoords $plane");
        if (visible_mode!=overlay_mode) 
          overlay_mode = visible_mode;  /* suppress mode change on select */
        updateflag = TRUE;
        break;
      case MIDDLEMOUSE:
        if (val == 0)  break;
        while (getbutton(MIDDLEMOUSE))
        {
          sx = getvaluator(MOUSEX);
          sy = getvaluator(MOUSEY);
          select_pixel(sx,sy);
          for (i= -prad;i<=prad;i++)
          for (j= -prad;j<=prad;j++)
            im[imc/zf][(ydim-1-ic+i)/zf][jc/zf+j] = 255;
          changed[imc/zf] = TRUE;
        }
        updateflag = TRUE;
        break;
      case RIGHTMOUSE:
        if (val == 0)  break;
        while (getbutton(RIGHTMOUSE))
        {
          sx = getvaluator(MOUSEX);
          sy = getvaluator(MOUSEY);
          select_pixel(sx,sy);
          for (i= -prad;i<=prad;i++)
          for (j= -prad;j<=prad;j++)
            im[imc/zf][(ydim-1-ic+i)/zf][jc/zf+j] = 0;
          changed[imc/zf] = TRUE;
        }
        updateflag = TRUE;
        break;
      case LEFTALTKEY: case RIGHTALTKEY:
        if (val)  altkeypressed = TRUE;
        else      altkeypressed = FALSE;
        break;
      case KEYBD:
        /*if (altkeypressed)*/
        switch((char)val) {
          case '0': record_swapbuffers();
                    break;
          case '1': overlay_mode = TARGET;
                    updateflag = TRUE;
                    break;
          case '2': overlay_mode = MOVEABLE;
                    updateflag = TRUE;
                    break;
          case 'x': plane=SAGITTAL;
                    updateflag = TRUE; 
                    break;
          case 'y': plane=HORIZONTAL;
                    updateflag = TRUE;
                    break;
          case 'z': plane=CORONAL;
                    updateflag = TRUE;
                    break;
        }
      }
  }
#endif
  return GL_FALSE;
}

void  open_window(char *name)
{
#ifdef OPENGL
  XSizeHints hin;

  if (openglwindowflag) {
    printf("medit: ### GL window already open: can't open second\n");PR return;}

  /* TKO_DEPTH because all OpenGL 4096 visuals have depth buffer!! */

  #ifdef RGB
  tkoInitDisplayMode(TKO_DOUBLE | TKO_RGB | TKO_DEPTH);
  #else
  tkoInitDisplayMode(TKO_DOUBLE | TKO_INDEX | TKO_DEPTH);
  #endif

  if (!initpositiondoneflag)
    tkoInitPosition(MOTIF_XFUDGE+wx0,(1024-wy0-ydim)+MOTIF_XFUDGE,xdim,ydim);
  if (!tkoInitWindow(name)) {
    printf("register: ### tkoInitWindow(name) failed\n");exit(1);}
  hin.max_width = hin.max_height = 3*xnum + xnum/2;  /* maxsize */
  hin.min_aspect.x = hin.max_aspect.x = xdim;        /* keepaspect */
  hin.min_aspect.y = hin.max_aspect.y = ydim;
  hin.flags = PMaxSize|PAspect;
  XSetWMNormalHints(xDisplay, w.wMain, &hin);

  /* TODO: bitmap drawing is slower */
  /* (did: glDisable's,4bit,RasPos:+0.5,glOrtho-1,z=0,wintitle,clrdepth) */
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glOrtho(0.0, (double)(xnum-1), 0.0, (double)(ynum-1), -1.0, 1.0);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

#else  /* use gl calls */
  if (openglwindowflag) {
    printf("medit: ### GL window already open: can't open second\n");PR return;}
  prefposition(wx0, wx0+xdim-1, wy0, wy0+ydim-1);
  foreground();  /* tcl prompt */
  winopen(name);
  keepaspect(1,1);
  stepunit((long)xnum,(long)ynum);
  maxsize((long)(3*xnum),(long)(3*ynum));
  winconstraints();
  cmode();
  doublebuffer();
  gconfig();
  ortho2(0.0, (float)(xnum-1), 0.0, (float)(ynum-1));
  zbuffer(FALSE);

  qdevice(ESCKEY);
  qdevice(REDRAW);
  qdevice(LEFTMOUSE);
  qdevice(RIGHTMOUSE);
  qdevice(MIDDLEMOUSE);
  qdevice(UPARROWKEY);
  qdevice(DOWNARROWKEY);
  qdevice(LEFTARROWKEY);
  qdevice(RIGHTARROWKEY);
  qdevice(INSERTKEY);
  qdevice(DELKEY);
  qdevice(HOMEKEY);
  qdevice(ENDKEY);
  qdevice(PAGEUPKEY);
  qdevice(PAGEDOWNKEY);
  qdevice(LEFTCTRLKEY);
  qdevice(RIGHTCTRLKEY);
  qdevice(LEFTALTKEY);
  qdevice(RIGHTALTKEY);
  qdevice(KEYBD);
#endif

  openglwindowflag = TRUE;
}

void  resize_window_intstep()
{
#ifdef OPENGL
  int tzf;

  tzf = rint((float)w.w/(float)xnum);
  tzf = (tzf<1)?1:(tzf>3)?3:tzf;
  if (w.w%xnum || w.h%ynum) {
    ozf = zf;
    zf = tzf;
    fsf = (float)zf;
    xdim = zf*xnum;
    ydim = zf*ynum;
    XResizeWindow(xDisplay, w.wMain, xdim, ydim);
    if (TKO_HAS_OVERLAY(w.type))
      XResizeWindow(xDisplay, w.wOverlay, xdim, ydim);
    resize_buffers(xdim, ydim);
    w.w = xdim;
    w.h = ydim;

    imc = (zf*imc)/ozf;
    ic = (zf*ic)/ozf;
    jc = (zf*jc)/ozf;
  }
#endif
}

void move_window(int x, int y)
{
#ifdef OPENGL
  if (openglwindowflag) {
    XMoveWindow(xDisplay, w.wMain, x, y);
    w.x = x;
    w.y = y;
  }
  else if (!initpositiondoneflag) {
    tkoInitPosition(x,y,xdim,ydim);
    initpositiondoneflag = TRUE;
  }
  else ;
#endif
}

void reload_buffers()
{
  
  /* printf("reloading buffers with %d x %d\n",xdim,ydim); */
  /*if(overlay_mode==TARGET) {
    glReadPixels(0,0,xdim,ydim,GL_RGB,GL_UNSIGNED_BYTE,blinkbuft);
    overlay_mode = MOVEABLE;
    glClear(GL_COLOR_BUFFER_BIT);
    draw_image(imc,ic,jc); */
    /* overlay_mode = TARGET; */
  /*glReadPixels(0,0,xdim,ydim,GL_RGB,GL_UNSIGNED_BYTE,blinkbufm);
    blinktop = 1;
  } else {
    glReadPixels(0,0,xdim,ydim,GL_RGB,GL_UNSIGNED_BYTE,blinkbufm);
    overlay_mode = TARGET;
    glClear(GL_COLOR_BUFFER_BIT);
    draw_image(imc,ic,jc);*/
    /* overlay_mode = MOVEABLE; */
  /*glReadPixels(0,0,xdim,ydim,GL_RGB,GL_UNSIGNED_BYTE,blinkbuft);
    blinktop = 0;
  }  
  invalid_buffers = 0;
  */
}

void blinkbuffers()
{
  if (blinkflag) {
#ifdef RGB
    /* if(blinkdelay==1) { */
      /* printf("reaload for delay 1\n"); */
    /* reload_buffers(); */
      /* printf("saving frontbuffer");
   glReadPixels(0,0,512,512,GL_RGB,GL_UNSIGNED_BYTE,blinkbuf); */
    /* } */
#endif
    if (blinkdelay<0) {
      record_swapbuffers();
      /*sginap(blinktime);*/
      /*usecnap(blinktime*10000);*/
    }
    else
      blinkdelay--;  
  }
}

void record_swapbuffers()    /* called by compare button */
{
  int swaptmp;
  
#ifdef RGB
/* if(blinkflag) {
    if(invalid_buffers == 1) {
      reload_buffers();
    } else {
      if(blinktop==1) {
        rectwrite(0,0,xdim-1,ydim-1,blinkbuft);
      } else { 
        rectwrite(0,0,xdim-1,ydim-1,blinkbufm);
      }
    }
    blinktop = (blinktop == 0) ? 1 : 0;
    } */
#endif
  swapbuffers();
  swaptmp = visible_plane;
  visible_plane = last_visible_plane;
  last_visible_plane = swaptmp;
  
  swaptmp = visible_mode;
  visible_mode = last_visible_mode;
  last_visible_mode = swaptmp;
}

void resize_buffers(int x,int y)
{
  free(vidbuf);
  free(binbuff);
  free(blinkbuft);
  free(blinkbufm);
  
#ifdef RGB
  vidbuf = (GLubyte *)lcalloc(3*(size_t)(x*y),(size_t)sizeof(GLubyte));
  blinkbuft = (GLubyte *)lcalloc(3*(size_t)(x*y),(size_t)sizeof(GLubyte));
  blinkbufm = (GLubyte *)lcalloc(3*(size_t)(x*y),(size_t)sizeof(GLubyte));
#else
  vidbuf = (Colorindex *)lcalloc((size_t)(x*y),(size_t)sizeof(Colorindex));
#endif
  binbuff = (unsigned char *)lcalloc((size_t)(3*x*y),(size_t)sizeof(char));
}

void  read_reg(char *fname)
{

#ifdef USE_mriTransform

  /* reads the register file and inits our transformation */
  InitFunctionalTransform( fname );

#else

  int i,j;
  FILE *fp;
  /* char junk[NAME_LENGTH]; */

  

  fp = fopen(fname,"r");
  if (fp==NULL) {
    printf("register: ### File %s not found (use -mkdefault option)\n",fname);
    exit(1); }

  for (i=0;i<4;i++)  /* idmat if none */
  for (j=0;j<4;j++)
    tm[i][j] = (i==j);

  if (MATCH(pname,"nobody")) fscanf(fp,"%s",pname); /* use if pname not set */
  else                       fscanf(fp,"%*s");      /* else ignore */
  fscanf(fp,"%lf",&ps_2);   /* hires: overlaps COR-.info */
  fscanf(fp,"%lf",&st_2);   /* hires: overlaps COR-.info */
  fscanf(fp,"%lf",&fscale_2);
  printf("register:   reg image transform matrix: \n");
  for (i=0;i<4;i++)
  for (j=0;j<4;j++) {
    fscanf(fp,"%f",&tm[i][j]);
    printf("%f ",tm[i][j]);
    if (j==3) printf("\n");
  }
  fclose(fp);
  fflush(stdout);

#endif /* USE_mriTransform */
}

void write_reg(char *fname)
{
  int i,j;
  FILE *fp;

  make_backup(fname);

  fp = fopen(fname,"w");
  if (fp==NULL){printf("register: ### can't create file %s\n",fname);PR return;}
  fprintf(fp,"%s\n",pname);
  fprintf(fp,"%f\n",ps_2);
  fprintf(fp,"%f\n",st_2);
  fprintf(fp,"%f\n",fscale_2);
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      fprintf(fp,"%13.6e ",tm[i][j]);
    fprintf(fp,"\n");
  }
  printf("register: file %s written\n",fname);PR
  fclose(fp);
  editedmatrix = FALSE;
}

void write_default_reg(char *fname)
{
  int i,j;
  char cwd[NAME_LENGTH];
  /* FILE *fp; */

  ps_2 = 1.562500;
  st_2 = 4.000000;
  fscale_2 = 0.200000;
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = (i==j);

  set_session_subject(getcwd(cwd,NAME_LENGTH));
  write_reg(fname);
  printf("register: ### default MGH register.dat written in cwd\n");
  printf("register:   subject name     = %s\n",pname);
  printf("register:   pixel size       = %f\n",ps_2);
  printf("register:   slice thickness  = %f\n",st_2);
  printf("register:   default bright   = %f\n",fscale_2);
  printf("register:   transform matrix = identity matrix\n");
}

void cor2anareg(char *fpref,char *ana,char *reg)
{
  FILE *fpi,*fpa,*fpr;
  char fname[NAME_LENGTH], cwd[NAME_LENGTH];
  int i,j;

  sprintf(fname,"%s.info",fpref);
  fpi = fopen(fname,"r");
  if (fpi==NULL){printf("register: ### File %s not found\n",fname);exit(1);}
  fscanf(fpi,"%*s %d",&imnr0_2);
  fscanf(fpi,"%*s %d",&imnr1_2);
  fscanf(fpi,"%*s %*d");
  fscanf(fpi,"%*s %d",&xnum_2);
  fscanf(fpi,"%*s %d",&ynum_2);
  fscanf(fpi,"%*s %*f");
  fscanf(fpi,"%*s %lf",&ps_2);
  fscanf(fpi,"%*s %lf",&st_2);
  ps_2 *= 1000;
  st_2 *= 1000;
  numimg_2 = imnr1_2-imnr0_2+1;
  fclose(fpi);

  set_session_subject(getcwd(cwd,NAME_LENGTH));
  fscale_2 = 1.0;
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = (i==j);

  fpr = fopen(reg,"r");
  if (fpr==NULL) {
    write_reg(reg);
    printf("\nregister: made register.dat from local COR-.info\n");
  }
  else {
    read_reg(reg);
    printf("\nregister: register.dat already there ...not overwritten\n");
  }
  printf("register:   subject name     = %s\n",pname);
  printf("register:   pixel size       = %f\n",ps_2);
  printf("register:   slice thickness  = %f\n",st_2);
  printf("register:   default bright   = %f\n",fscale_2);

  make_backup(ana);
  fpa = fopen(ana,"w");
  if (fpa==NULL){printf("register: ### can't create file %s\n",ana);exit(1);}
  fprintf(fpa,".\n");
  fprintf(fpa,"COR-%%03d\n");
  fprintf(fpa,"%d 1\n",numimg_2);
  fprintf(fpa,"%d %d\n",xnum_2,ynum_2);
  fclose(fpa);
  printf("\nregister: made analyse.dat from local COR-.info\n");
  printf("register:   structdir         = .\n");
  printf("register:   imageformat       = COR-%%03d\n");
  printf("register:   numslices numreps = %d 1\n",numimg_2);
  printf("register:   xpixnum ypixnum   = %d %d\n\n",xnum_2,ynum_2);
}

void set_session_subject(char *str)  /* find session to find name file to find subject */
{
  int i,j,k;
  char *word, wd[NAME_LENGTH], path[30][NAME_LENGTH];
  char fname[NAME_LENGTH];
  FILE *fp;

  strcpy(wd,str);
  word = strtok(wd,"/");
  strcpy(path[0],word);
  i = 1;
  while ((word = strtok(NULL,"/")) != NULL) { strcpy(path[i++],word); }
  if (MATCH(path[i-2],"image")) {
    printf("register: in new (APD2) format functdir\n");  j=i-1;  k=i-2;
  } else if (strspn(path[i-3],"0123456789")==5){
    printf("register: in old (APD1) format functdir\n");  j=i-2;  k=i-3;
  } else {
    printf("register: in unknown functdir type\n");       j=i;    k=i;
  }
  sprintf(srname,"/%s",path[0]);
  for(i=1;i<j;i++) sprintf(srname,"%s/%s",srname,path[i]);   /* session */
  sprintf(psrname,"/%s",path[0]);
  for(i=1;i<k;i++) sprintf(psrname,"%s/%s",psrname,path[i]); /* parent */
  sprintf(fname,"%s/name",psrname);
  if (MATCH("nobody",pname)) {  /* only if pname not set, from namefile */
    fp = fopen(fname,"r");
    if (fp==NULL) {
      printf("register: File %s not found: subj set to \"nobody\"\n",fname);
      strcpy(pname,"nobody");
    } else {
      fscanf(fp,"%s",pname);
      printf("register: subject name (from namefile): %s\n",pname);
      fclose(fp);
    }
  }
}

void make_backup(char *fname)
{
  char command[2*NAME_LENGTH];
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp!=NULL) {
    sprintf(command,"cp %s %s~",fname,fname); 
    system(command);
    fclose(fp);
  }
}

void save_rgb(char *fname)
{
  if (scrsaveflag) { scrsave_to_rgb(fname); }
  else             { pix_to_rgb(fname);     }
}

void scrsave_to_rgb(char *fname)  /* about 2X faster than pix_to_rgb */
{
  char command[2*NAME_LENGTH];
  FILE *fp;
  long xorig,xsize,yorig,ysize;
  int x0,y0,x1,y1;

  getorigin(&xorig,&yorig);
  getsize(&xsize,&ysize);

  x0 = (int)xorig;  x1 = (int)(xorig+xsize-1);
  y0 = (int)yorig;  y1 = (int)(yorig+ysize-1);
  fp = fopen(fname,"w");
  if (fp==NULL){printf("register: ### can't create file %s\n",fname);PR return;}
  fclose(fp);
  sprintf(command,"scrsave %s %d %d %d %d\n",fname,x0,x1,y0,y1);
  system(command);
  printf("register: file %s written\n",fname);PR
}

void pix_to_rgb(char *fname)
{
  /*
#ifdef OPENGL
  GLint swapbytes, lsbfirst, rowlength;
  GLint skiprows, skippixels, alignment;
  GLenum format;
  IMAGE *image;
  int y,width,height,size;
  unsigned short *r,*g,*b;
  unsigned short  *red, *green, *blue;
  FILE *fp;

  width = (int)xdim;
  height = (int)ydim;
  size = width*height;

  red = (unsigned short *)calloc(size, sizeof(unsigned short));
  green = (unsigned short *)calloc(size, sizeof(unsigned short));
  blue = (unsigned short *)calloc(size, sizeof(unsigned short));

  glGetIntegerv(GL_UNPACK_SWAP_BYTES, &swapbytes);
  glGetIntegerv(GL_UNPACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength);
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows);
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  glReadPixels(0, 0, width, height, GL_RED,GL_UNSIGNED_SHORT, (GLvoid *)red);
  glReadPixels(0, 0, width, height, GL_GREEN,GL_UNSIGNED_SHORT,(GLvoid *)green);  glReadPixels(0, 0, width, height, GL_BLUE,GL_UNSIGNED_SHORT, (GLvoid *)blue);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes);
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels);
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment);

  fp = fopen(fname,"w");
  if (fp==NULL){printf("register: ### can't create file %s\n",fname);PR return;}
  fclose(fp);
  image = iopen(fname,"w",RLE(1), 3, width, height, 3);
  for(y = 0 ; y < height; y++) {
    r = red + y * width;
    g = green + y * width;
    b = blue + y * width;
    putrow(image, r, y, 0);
    putrow(image, g, y, 1);
    putrow(image, b, y, 2);
  }
  iclose(image);
  free(red); free(green); free(blue);

  printf("register: file %s written\n",fname);PR
#else
  printf("register: ### pix_to_rgb implemented only in OpenGL version\n");PR
  return;
#endif
  */
}

void downslice()
{
  if (plane==CORONAL)
    imc = (imc<zf)?imnr1*zf-zf+imc:imc-zf;
  else if (plane==HORIZONTAL)
    ic = (ic<zf)?ydim-zf+ic:ic-zf;
  else if (plane==SAGITTAL)
    jc = (jc<zf)?xdim-zf+jc:jc-zf;
}

void upslice()
{
  if (plane==CORONAL)
    imc = (imc>=imnr1*zf-zf)?imc+zf-imnr1*zf:imc+zf;
  else if (plane==HORIZONTAL)
    ic = (ic>=ydim-zf)?ic+zf-ydim:ic+zf;
  else if (plane==SAGITTAL)
    jc = (jc>=xdim-zf)?jc+zf-xdim:jc+zf;
}

void goto_point(char *dir)
{
  char fname[NAME_LENGTH];
  FILE *fp;
  float xpt,ypt,zpt;

  sprintf(fname,"%s/edit.dat",dir);
  fp=fopen(fname,"r");
  if (fp==NULL) {printf("register: ### File %s not found\n",fname);PR return;}
  fscanf(fp,"%f %f %f",&xpt,&ypt,&zpt);
  fclose(fp);
  set_cursor(xpt,ypt,zpt);
}

void write_point(char *dir)
{
  char fname[NAME_LENGTH];
  FILE *fp;
  float xpt,ypt,zpt;

  sprintf(fname,"%s/edit.dat",dir);
  fp=fopen(fname,"w");
  if (fp==NULL){printf("register: ### can't create file %s\n",fname);PR return;}
  xpt = xx1-ps*jc/fsf;
  ypt = yy0+st*imc/fsf;
  zpt = zz1-ps*(255.0-ic/fsf);
  fprintf(fp,"%f %f %f\n",xpt,ypt,zpt);
  fclose(fp);
}

void rotate_brain(float a,char c)
{
  int i,j,k;
  float m1[4][4],m2[4][4];
  float sa,ca;

  if (c=='x')
  {
    translate_brain(yc,'y');
    translate_brain(zc,'z');
  } else
  if (c=='z')
  {
    translate_brain(xc,'x');
    translate_brain(yc,'y');
  } else
  if (c=='y')
  {
    translate_brain(xc,'x');
    translate_brain(zc,'z');
  }

  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    m1[i][j] = (i==j)?1.0:0.0;
  a = a*M_PI/1800;
  sa = sin(a);
  ca = cos(a);
  if (c=='y')
  {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  } else
  if (c=='x')
  {
    m1[1][1] = m1[2][2] = ca;
    m1[1][2] = -(m1[2][1] = sa);
  } else
  if (c=='z')
  {
    m1[0][0] = m1[1][1] = ca;
    m1[0][1] = -(m1[1][0] = sa);
  } else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  {
    m2[i][j] = 0;
    for (k=0;k<4;k++)
      m2[i][j] += tm[i][k]*m1[k][j];
  }
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = m2[i][j];

  if (c=='x')
  {
    translate_brain(-yc,'y');
    translate_brain(-zc,'z');
  } else
  if (c=='z')
  {
    translate_brain(-xc,'x');
    translate_brain(-yc,'y');
  } else
  if (c=='y')
  {
    translate_brain(-xc,'x');
    translate_brain(-zc,'z');
  }
  editedmatrix = TRUE;
}

void align_points()
{
  if (plane==SAGITTAL)
  {
    translate_brain(-(yc-yc_old),'y');
    translate_brain(-(zc-zc_old),'z');
  }
  if (plane==HORIZONTAL)
  {
    translate_brain(-(xc-xc_old),'x');
    translate_brain(-(yc-yc_old),'y');
  }
  if (plane==CORONAL)
  {
    translate_brain(-(xc-xc_old),'x');
    translate_brain(-(zc-zc_old),'z');
  }
}

void translate_brain(float a, char c)
{
  int i,j,k;
  float m1[4][4],m2[4][4];

  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    m1[i][j] = (i==j)?1.0:0.0;
  if (c=='y')
    m1[1][3] = a;
  else if (c=='x')
    m1[0][3] = a;
  else if (c=='z')
    m1[2][3] = a;
  else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  {
    m2[i][j] = 0;
    for (k=0;k<4;k++)
      m2[i][j] += tm[i][k]*m1[k][j];
  }
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = m2[i][j];
  editedmatrix = TRUE;
}

void scale_brain(float s, char c)
{
  int i,j,k;
  float m1[4][4],m2[4][4];

  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    m1[i][j] = (i==j)?1.0:0.0;
  if (c=='x')
    m1[0][0] = s;
  else if (c=='y')
    m1[1][1] = s;
  else if (c=='z')
    m1[2][2] = s;
  else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  {
    m2[i][j] = 0;
    for (k=0;k<4;k++)
      m2[i][j] += tm[i][k]*m1[k][j];
  }
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = m2[i][j];
  editedmatrix = TRUE;
}

void mirror_brain()
{
  int j;

  for (j=0;j<4;j++)
    tm[0][j] = -tm[0][j];
  editedmatrix = TRUE;
}

void set_cursor(float xpt,float ypt,float zpt)
{
  double dzf;

  if (ptype==0) /* Horizontal */
  {
    jpt = (int)((xx1-xpt)*zf/ps+0.5);
    ipt = (int)((ypt-yy0)*zf/ps+0.5);
    impt = (int)((zpt-zz0)*zf/st+0.5);
  } else if (ptype==2) /* Coronal */
  {
    jpt = (int)((xx1-xpt)*zf/ps+0.5);
    ipt = (int)((255.0-(zz1-zpt)/ps)*zf+0.5);
    impt = (int)((ypt-yy0)*zf/st+0.5);
  } else if (ptype==1) /* Sagittal */
  {
    jpt = (int)((xx1-xpt)*zf/ps+0.5);
    ipt = (int)((ypt-yy0)*zf/ps+0.5);
    impt = (int)((zpt-zz0)*zf/st+0.5);
  }
  dzf = (double)zf;
  imc = zf * (int)(rint((double)impt/dzf));  /* round to slice */
  ic =  zf * (int)(rint((double)ipt/dzf));
  jc =  zf * (int)(rint((double)jpt/dzf));
  /*jc = jpt;*/
  /*ic = ipt;*/
  /*imc = impt;*/
  jpt=ipt=impt = -1;
}

void set_scale()
{
  Colorindex i;
  float f;
  short v;

  for (i=0;i<NUMVALS;i++)
  {
/*
    f = i/fscale+fthresh;
    f = ((f<0.0)?0.0:((f>1.0)?1.0:f));
    f = pow(f,fsquash);
    v = f*fscale+0.5;
*/
    f = i/fscale;
    f = 1.0/(1.0+exp(-fsquash*(f-fthresh)));
    v = f*fscale+0.5;
#ifdef RGB
    colormap[i] = (unsigned char)v;
#else
    mapcolor(i+MAPOFFSET,v,v,v);
#endif
  }
  mapcolor(NUMVALS+MAPOFFSET,v,0.0,0.0);
}

void redraw()
{
  color(0);
  clear();
  draw_image(imc,ic,jc);

  last_visible_mode = visible_mode;
  visible_mode = overlay_mode;

  last_visible_plane = visible_plane;
  visible_plane = plane;

}

void pop_gl_window()
{
#ifdef OPENGL
  XRaiseWindow(xDisplay, w.wMain);
#else
  winpop();
#endif
}

#if 0
#include "proto.h"
#endif
void mri2pix(float xpt,float ypt,float zpt,int *jpt,int *ipt,int *impt)
{
  if (ptype==0) /* Horizontal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((ypt-yy0)/ps+0.5);
    *impt = (int)((zpt-zz0)/st+0.5);
  } else if (ptype==2) /* Coronal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((255.0-(zz1-zpt)/ps)+0.5);
    *impt = (int)((ypt-yy0)/st+0.5);
  } else if (ptype==1) /* Sagittal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((ypt-yy0)/ps+0.5);
    *impt = (int)((zpt-zz0)/st+0.5);
  }
}

int imval(float px,float py,float pz)
{
  float x,y,z;
  int j,i,imn;

  x = px*TM[0][0]+py*TM[0][1]+pz*TM[0][2]+TM[0][3]+par[0];
  y = px*TM[1][0]+py*TM[1][1]+pz*TM[1][2]+TM[1][3]+par[1];
  z = px*TM[2][0]+py*TM[2][1]+pz*TM[2][2]+TM[2][3]+par[2];
  mri2pix(x,y,z,&j,&i,&imn);
  if (imn>=0&&imn<numimg&&i>=0&&i<ynum&&j>=0&&j<xnum)
    return(im[imn][ynum-1-i][j]);
  else return 0;
}

float Error(int p,float dp)
{
  int i,num;
  float mu,error,sum;
  float mu1,mu2,sum1,sum2;

  if (p>=0)
    par[p] += dp;
  mu = mu1 = mu2 = 0;
  num = 0;
  for (i=0;i<npts;i++)
  {
    mu += imval(ptx[i],pty[i],ptz[i]);
    mu1 += imval(ptx[i]*0.9,pty[i]*0.9,ptz[i]*0.9);
    mu2 += imval(ptx[i]*1.05,pty[i]*1.05,ptz[i]*1.05);
    num ++;
  }
  mu /= num;
  mu1 /= num;
  mu2 /= num;
  sum = sum1 = sum2 = 0;
  num = 0;
  for (i=0;i<npts;i++)
  {
    error = imval(ptx[i],pty[i],ptz[i])-mu;
    sum += error*error;
    error = imval(ptx[i]*0.9,pty[i]*0.9,ptz[i]*0.9)-mu1;
    sum1 += error*error;
    error = imval(ptx[i]*1.05,pty[i]*1.05,ptz[i]*1.05)-mu2;
    sum2 += error*error;
    num ++;
  }
  sum = sqrt((sum+sum2)/num);
  if (p>=0)
    par[p] -= dp;
  return sum;
}

void optimize(int maxiter)
{
  float lambda = 0.03;
  float epsilon = 0.1;
  float momentum = 0.8;
  int iter,p;
  float dE[3];
  float error;

  for (iter=0;iter<maxiter;iter++)
  {
    error = Error(-1,0);
    printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
           iter,par[0],par[1],par[2],error);
    for (p=0;p<3;p++)
    {
      dE[p] = tanh((Error(p,epsilon/2)-Error(p,-epsilon/2))/epsilon);
    }
    for (p=0;p<3;p++)
    {
      par[p] += (dpar[p] = momentum*dpar[p] - lambda*dE[p]);
    }
  }
  error = Error(-1,0);
  printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
         iter,par[0],par[1],par[2],error);
}

void optimize2()
{
  float epsilon = 0.5;
  /* int p; */
  float p0,p1,p2,p0min,p1min,p2min;
  float error,minerror;

  error = Error(-1,0);
  minerror = error; p0min = p1min = p2min = 0 ;
  printf("%5.2f %5.2f %5.2f %7.3f\n",
         par[0],par[1],par[2],error);
  for (p0 = -10;p0 <= 10;p0 += epsilon)
  for (p1 = -10;p1 <= 10;p1 += epsilon)
  for (p2 = -10;p2 <= 10;p2 += epsilon)
  {
    par[0] = p0;
    par[1] = p1;
    par[2] = p2;
    error = Error(-1,0);
    if (error<minerror)
    {
      printf("%5.2f %5.2f %5.2f %7.3f\n",
             par[0],par[1],par[2],error);
      minerror = error;
      p0min = p0;
      p1min = p1;
      p2min = p2;
    }
  }
  par[0] = p0min;
  par[1] = p1min;
  par[2] = p2min;
  error = Error(-1,0);
  printf("%5.2f %5.2f %5.2f %7.3f\n",
         par[0],par[1],par[2],error);
}

void read_images(char *fpref)
{
  int i,j,k;                   /* loop counters */
  FILE *fptr;
  char fname[NAME_LENGTH];

        sprintf(fname,"%s.info",fpref);
        fptr = fopen(fname,"r");
        if (fptr==NULL) {
          printf("register: ### File %s not found\n",fname);exit(1);}
  fscanf(fptr,"%*s %d",&imnr0);
  fscanf(fptr,"%*s %d",&imnr1);
  fscanf(fptr,"%*s %d",&ptype);
  fscanf(fptr,"%*s %d",&xnum);
  fscanf(fptr,"%*s %d",&ynum);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %f",&ps);
  fscanf(fptr,"%*s %f",&st);
  fscanf(fptr,"%*s %*f");
        fscanf(fptr,"%*s %f",&xx0); /* strtx */
        fscanf(fptr,"%*s %f",&xx1); /* endx */
        fscanf(fptr,"%*s %f",&yy0); /* strty */
        fscanf(fptr,"%*s %f",&yy1); /* endy */
        fscanf(fptr,"%*s %f",&zz0); /* strtz */
        fscanf(fptr,"%*s %f",&zz1); /* endz */
        ps *= 1000;
        st *= 1000;
        xx0 *= 1000;
        xx1 *= 1000;
        yy0 *= 1000;
        yy1 *= 1000;
        zz0 *= 1000;
        zz1 *= 1000;
  fclose(fptr);
  numimg = imnr1-imnr0+1;
  xdim=xnum*zf;
  ydim=ynum*zf;
/*
  printf("imnr0=%d,imnr1=%d,numimg=%d,xdim=%d,ydim=%d\n",
                imnr0,imnr1,numimg,xdim,ydim);
*/

/* Allocate memory */

  bufsize = ((unsigned long)xnum)*ynum;
  buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
  for (k=0;k<numimg;k++)
  {
    im[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++)
    {
      im[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
  }
  for (k=0;k<6;k++)
  {
    sim[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++)
    {
      sim[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
        }

  for (k=0;k<numimg;k++)
  {
          changed[k] = FALSE;
    file_name(fpref,fname,k+imnr0,"%03d");
          fptr = fopen(fname,"rb");
          if (fptr==NULL) {
            printf("register: ### File %s not found\n",fname);exit(1);}
          fread(buf,sizeof(char),bufsize,fptr);
          buffer_to_image(buf,im[k],xnum,ynum);
          fclose(fptr);
/*
          printf("file %s read in\n",fname);
*/
  }
        printf("register: done reading target COR images\n");

  for (k=0;k<numimg;k++)
        for (i=0;i<ynum;i++)
        for (j=0;j<xnum;j++)
        {
          if (im[k][i][j]/2>sim[3][i][j]) sim[3][i][j] = im[k][i][j]/2;
          if (im[k][i][j]/2>sim[4][k][j]) sim[4][k][j] = im[k][i][j]/2;
          if (im[k][i][j]/2>sim[5][i][k]) sim[5][i][k] = im[k][i][j]/2;
        }
        for (i=0;i<ynum;i++)
        for (j=0;j<xnum;j++)
        for (k=0;k<3;k++)
          sim[k][i][j] = sim[k+3][i][j];
       
}

void read_second_images(char *fpref)
{
  unsigned long n;
  int i,j,k,xdim_2b,ydim_2b;
  double ps_2b,st_2b;
  FILE *fptr;
  char fname[NAME_LENGTH];

  xdim_2b=xdim_2; ydim_2b=ydim_2; ps_2b=ps_2; st_2b=st_2;

  sprintf(fname,"%s.info",fpref);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {
    printf("register: ### File %s not found\n",fname);exit(1);}
  fscanf(fptr,"%*s %d",&imnr0_2);
  fscanf(fptr,"%*s %d",&imnr1_2);
  fscanf(fptr,"%*s %d",&ptype);
  fscanf(fptr,"%*s %d",&xnum_2);
  fscanf(fptr,"%*s %d",&ynum_2);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %lf",&ps_2);
  fscanf(fptr,"%*s %lf",&st_2);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %f",&xx0_2); /* strtx */
  fscanf(fptr,"%*s %f",&xx1_2); /* endx */
  fscanf(fptr,"%*s %f",&yy0_2); /* strty */
  fscanf(fptr,"%*s %f",&yy1_2); /* endy */
  fscanf(fptr,"%*s %f",&zz0_2); /* strtz */
  fscanf(fptr,"%*s %f",&zz1_2); /* endz */
  ps_2 *= 1000;
  st_2 *= 1000;
  xx0_2 *= 1000;
  xx1_2 *= 1000;
  yy0_2 *= 1000;
  yy1_2 *= 1000;
  zz0_2 *= 1000;
  zz1_2 *= 1000;
  fclose(fptr);
  numimg_2 = imnr1_2-imnr0_2+1;
  xdim_2 = xnum_2;   /* no zoom here */
  ydim_2 = ynum_2;

  if (nslices!=numimg_2) {
    printf("register ### numslices mismatch--analyse.dat ignored \n");
    printf("  nslices: analyse.dat = %d   COR-.info = %d\n",nslices,numimg_2); }
  if (fabs(ps_2b-ps_2)>0.00001 || fabs(st_2b-st_2)>0.00001) {
    printf("register ### slicethick,pixsize mismatch--register.dat ignored\n");
    printf("  thick:   register.dat = %f   COR-.info = %f\n",st_2b,st_2);
    printf("  pixsize: register.dat = %f   COR-.info = %f\n",ps_2b,ps_2); }
  if (xdim_2b!=xdim_2 || ydim_2b!=ydim_2) {
    printf("register ### xdim,ydim mismatch--analyse.dat ignored\n");
    printf("  xdim:  analyse.dat = %d   COR-.info = %d\n",xdim_2b,xdim_2);
    printf("  ydim:  analyse.dat = %d   COR-.info = %d\n",ydim_2b,ydim_2); }

  bufsize_2 = ((unsigned long)xnum_2)*ynum_2;
  buf_2 = (unsigned char *)lcalloc(bufsize_2,sizeof(char));
  for (k=0;k<numimg_2;k++) {
    fim_2[k] = (float **)lcalloc(ynum_2,sizeof(float *));
    for (i=0;i<ynum_2;i++)
      fim_2[k][i] = (float *)lcalloc(xnum_2,sizeof(float));
  }

  for (k=0;k<numimg_2;k++) {
    file_name(fpref,fname,k+imnr0_2,"%03d");
    fptr = fopen(fname,"rb");
    if (fptr==NULL) {
      printf("register: ### File %s not found\n",fname);exit(1);}
    fread(buf_2,sizeof(float),bufsize_2,fptr);
    fclose(fptr);
    n=0;
    for (i=0;i<ynum_2;i++)
    for (j=0;j<xnum_2;j++) {
      fim_2[k][i][j] = buf_2[n++];  /* 64Meg <- 16Meg! */
    }
  }
  printf("register: done reading to-be-registered COR images\n");
}

void select_pixel(short sx, short sy)
{
  long ox,oy,lx,ly;

  getorigin(&ox,&oy);
  getsize(&lx,&ly);

  /* hack: x-y of click on tk win caught by qread when clicks buffered! */
  if (sx<ox || sx>ox+lx) sx = ox+lx/2;
  if (sy<oy || sy>oy+ly) sy = oy+ly/2;
/*
  printf("sx=%d, sy=%d, ox=%d, oy=%d, lx=%d, ly=%d\n",
         sx,sy,ox,oy,lx,ly);
*/
  if (plane==CORONAL)
  {
    ic = (sy-oy);
    jc = (sx-ox);
  } else
  if (plane==HORIZONTAL)
  {
    imc = (sy-oy);
    jc = (sx-ox);
  } else
  if (plane==SAGITTAL)
  {
    imc = (sx-ox);
    ic = (sy-oy);
  }
  xc = xx1-ps*jc/fsf; 
  yc = yy0+st*imc/fsf; 
  zc = zz1-ps*(255.0-ic/fsf);
  if (imc/zf>=0 && imc/zf<imnr1)
    printf("val=%d ",im[(int)(imc/zf)][(int)((ydim-1-ic)/zf)][(int)(jc/zf)]);
  printf("imnr=%5.1f, i=%5.1f, j=%5.1f (%5.1f, %5.1f, %5.1f)\n",
          imc/fsf,ic/fsf,jc/fsf,xc,yc,zc);
  PR
}

void transform(float x1,float y1,float z1,float *x2,float *y2,float *z2,float M[4][4])
{
  *x2 = x1*M[0][0]+y1*M[0][1]+z1*M[0][2]+M[0][3];
  *y2 = x1*M[1][0]+y1*M[1][1]+z1*M[1][2]+M[1][3];
  *z2 = x1*M[2][0]+y1*M[2][1]+z1*M[2][2]+M[2][3];
}

void draw_image(int imc,int ic,int jc)
{
  int i,j,imnr,k; /* ,p; */
  float x_1,y_1,z_1,x_2,y_2,z_2,f;
  int i_2,j_2,imnr_2;
  Colorindex v;

#ifdef RGB
  char* lvidbuf;
#endif

/*
  x_1 = xx1-ps*jc/fsf;
  y_1 = yy0+st*imc/fsf;
  z_1 = zz1-ps*(ynum-1-ic/fsf);
  printf("coord_1 : (%f,%f,%f)\n",x_1,y_1,z_1);
  transform(x_1,y_1,z_1,&x_2,&y_2,&z_2,tm);
  printf("coord_2 : (%f,%f,%f)\n",x_2,y_2,z_2);
*/

  k = 0;
  if (maxflag)
  {
    for (i=0;i<ydim;i++)
    for (j=0;j<xdim;j++)
    {
      if (plane==CORONAL) {

#ifdef RGB
        vidbuf[3*k] = colormap[sim[0][255-i/zf][j/zf]];
  vidbuf[3*k+1] = colormap[sim[0][255-i/zf][j/zf]];
  vidbuf[3*k+2] = colormap[sim[0][255-i/zf][j/zf]];
#else
  vidbuf[k] = sim[0][255-i/zf][j/zf]+MAPOFFSET;
#endif
      }
      else if (plane==HORIZONTAL) {
#ifdef RGB
        vidbuf[3*k] = colormap[sim[1][i/zf][j/zf]];
  vidbuf[3*k+1] = colormap[sim[1][i/zf][j/zf]];
  vidbuf[3*k+2] = colormap[sim[1][i/zf][j/zf]];
#else
  vidbuf[k] = sim[1][i/zf][j/zf]+MAPOFFSET;
#endif  
      }
      else if (plane==SAGITTAL) {
#ifdef RGB
        vidbuf[3*k] = colormap[sim[2][255-i/zf][j/zf]];
  vidbuf[3*k+1] = colormap[sim[2][255-i/zf][j/zf]];
  vidbuf[3*k+2] = colormap[sim[2][255-i/zf][j/zf]];
#else
  vidbuf[k] = sim[2][255-i/zf][j/zf]+MAPOFFSET;
#endif
      }
      k++;
    }
  }
  else
  {
#ifdef RGB
    if(overlay_mode==TARGET) {
      lvidbuf = blinkbuft;
      pswapnext = TARGET;
    }
    else {
      lvidbuf = blinkbufm;
      pswapnext = MOVEABLE;
    }
#endif
    if (plane==CORONAL)
    {
      for (i=0;i<ydim;i++)
      for (j=0;j<xdim;j++)
      {

  /* voxel screen to RAS */
        x_1 = xx1-ps*j/fsf;
        y_1 = yy0+st*imc/fsf;
        z_1 = zz1-ps*(ynum-1-i/fsf);

#ifdef USE_mriTransform

  ConvertRASToFuncIdx( x_1, y_1, z_1, &x_1, &y_2, &z_2 );

#else

  /* ana RAS to func RAS */
        transform(x_1,y_1,z_1,&x_2,&y_2,&z_2,tm);

  /* func RAS to func index */
        imnr_2 = floor((y_2-yy0_2)/st_2);
        i_2 = ydim_2-1-(zz1_2-z_2)/ps_2;
        j_2 = (xx1_2-x_2)/ps_2;

#endif /* USE_mriTransform */

        if (overlay_mode==TARGET)
        {
          if (imc/zf>=0&&imc/zf<imnr1)
            v = im[imc/zf][(ydim-1-i)/zf][j/zf]+MAPOFFSET;
          else v=MAPOFFSET;
          if (i==ipt||j==jpt) v=255-(v-MAPOFFSET)+MAPOFFSET;
          if (maskflag)




          if (imnr_2<0||imnr_2>=imnr1_2||i_2<0||i_2>=ydim_2||j_2<0||j_2>=xdim_2)
            v=MAPOFFSET;
          if ((i==ic&&abs(j-jc)<=2)||
              (j==jc&&abs(i-ic)<=2)) 
#ifdef RGB
    {
      lvidbuf[3*k] = 255;
      lvidbuf[3*k+1] = 0;
      lvidbuf[3*k+2] = 0;
    } else {
      lvidbuf[3*k] = colormap[v];
      lvidbuf[3*k+1] = colormap[v];
      lvidbuf[3*k+2] =colormap[v];  
    }
#else
    v=NUMVALS+MAPOFFSET;
    vidbuf[k] = v; 
#endif
          k++;
        } 
        else if (overlay_mode==MOVEABLE)
        {
          if (imnr_2>=0&&imnr_2<imnr1_2&&
              i_2>=0&&i_2<ydim_2&&j_2>=0&&j_2<xdim_2)
          {
            f = fscale_2*fim_2[imnr_2][ydim_2-1-i_2][j_2];
            v = ((f<0)?0:(f>255)?255:f)+MAPOFFSET;
          }
          else v=MAPOFFSET;
          if (i==ipt||j==jpt) v=255-(v-MAPOFFSET)+MAPOFFSET;
/*
          if (i==ic&&j==jc)
             printf("imnr_2 = %d (%f,%f,%f)\n",imnr_2,x_2,y_2,z_2);
*/

          if ((i==ic&&abs(j-jc)<=2)||
              (j==jc&&abs(i-ic)<=2)) 
#ifdef RGB
    {
      lvidbuf[3*k] = 255;
      lvidbuf[3*k+1] = 0;
      lvidbuf[3*k+2] = 0;
    } else {
      lvidbuf[3*k] = colormap[v];
      lvidbuf[3*k+1] = colormap[v];
      lvidbuf[3*k+2] = colormap[v];  
    }
#else
    v=NUMVALS+MAPOFFSET;
    vidbuf[k] = v; 
#endif
    k++;
        } 
    }
    } else
    if (plane==HORIZONTAL)
    {
      for (imnr=0;imnr<ydim;imnr++)
      for (j=0;j<xdim;j++)
      {
        x_1 = xx1-ps*j/fsf;
        y_1 = yy0+st*imnr/fsf;
        z_1 = zz1-ps*(ynum-1-ic/fsf);

#ifdef USE_mriTransform

  ConvertRASToFuncIdx( x_1, y_1, z_1, &x_1, &y_2, &z_2 );

#else

        transform(x_1,y_1,z_1,&x_2,&y_2,&z_2,tm);
        imnr_2 = floor((y_2-yy0_2)/st_2);
        i_2 = ydim_2-1-(zz1_2-z_2)/ps_2;
        j_2 = (xx1_2-x_2)/ps_2;

#endif /* USE_mriTransform */

        if (overlay_mode==TARGET)
        {
          if (imnr/zf>=0&&imnr/zf<imnr1)
            v = im[imnr/zf][(ydim-1-ic)/zf][j/zf]+MAPOFFSET;
          else v=MAPOFFSET;
          if (imnr==impt||j==jpt) v=255-(v-MAPOFFSET)+MAPOFFSET;
          if (maskflag)
          if (imnr_2<0||imnr_2>=imnr1_2||i_2<0||i_2>=ydim_2||j_2<0||j_2>=xdim_2)
            v=MAPOFFSET;
          if ((imnr==imc&&abs(j-jc)<=2)||
              (j==jc&&abs(imnr-imc)<=2))
          
#ifdef RGB
      {
        lvidbuf[3*k] = 255;
        lvidbuf[3*k+1] = 0;
        lvidbuf[3*k+2] = 0;
      } else {
        lvidbuf[3*k] = colormap[v];
        lvidbuf[3*k+1] = colormap[v];
        lvidbuf[3*k+2] = colormap[v];  
      }
#else
    v=NUMVALS+MAPOFFSET;
    vidbuf[k] = v; 
#endif
    k++;
        }
        else if (overlay_mode==MOVEABLE)
        {
          if (imnr_2>=0&&imnr_2<imnr1_2&&
              i_2>=0&&i_2<ydim_2&&j_2>=0&&j_2<xdim_2)
          {
            f = fscale_2*fim_2[imnr_2][ydim_2-1-i_2][j_2];
            v = ((f<0)?0:(f>255)?255:f)+MAPOFFSET;
          }
          else v=MAPOFFSET;
          if (imnr==impt||j==jpt) v=255-(v-MAPOFFSET)+MAPOFFSET;
/*
          if (imnr==imc&&j==jc)
             printf("imnr_2 = %d (%f,%f,%f)\n",imnr_2,x_2,y_2,z_2);
*/
    
          if ((imnr==imc&&abs(j-jc)<=2)||
              (j==jc&&abs(imnr-imc)<=2)) 
#ifdef RGB
    {
      lvidbuf[3*k] = 255;
      lvidbuf[3*k+1] = 0;
      lvidbuf[3*k+2] = 0;
    } else {
      lvidbuf[3*k] = colormap[v];
      lvidbuf[3*k+1] = colormap[v];
      lvidbuf[3*k+2] = colormap[v];  
    }
#else
    v=NUMVALS+MAPOFFSET;
    vidbuf[k] = v; 
#endif
          k++;
        } 
      }
    } else
    if (plane==SAGITTAL)
    {
      for (i=0;i<ydim;i++)
      for (imnr=0;imnr<xdim;imnr++)
      {
        x_1 = xx1-ps*jc/fsf;
        y_1 = yy0+st*imnr/fsf;
        z_1 = zz1-ps*(ynum-1-i/fsf);

#ifdef USE_mriTransform

  ConvertRASToFuncIdx( x_1, y_1, z_1, &x_1, &y_2, &z_2 );

#else

        transform(x_1,y_1,z_1,&x_2,&y_2,&z_2,tm);
        imnr_2 = floor((y_2-yy0_2)/st_2);
        i_2 = ydim_2-1-(zz1_2-z_2)/ps_2;
        j_2 = (xx1_2-x_2)/ps_2;

#endif /* USE_mriTransform */

        if (overlay_mode==TARGET)
        {
          if (imnr/zf>=0&&imnr/zf<imnr1)
            v = im[imnr/zf][(ydim-1-i)/zf][jc/zf]+MAPOFFSET;
          else v=MAPOFFSET;
          if (imnr==impt||i==ipt) v=255-(v-MAPOFFSET)+MAPOFFSET;
          if (maskflag)
          if (imnr_2<0||imnr_2>=imnr1_2||i_2<0||i_2>=ydim_2||j_2<0||j_2>=xdim_2)
            v=MAPOFFSET;
          if ((imnr==imc&&abs(i-ic)<=2)||
              (i==ic&&abs(imnr-imc)<=2)) 
#ifdef RGB
    {
      lvidbuf[3*k] = 255;
      lvidbuf[3*k+1] = 0;
      lvidbuf[3*k+2] = 0;
    } else {
      lvidbuf[3*k] = colormap[v];
      lvidbuf[3*k+1] = colormap[v];
      lvidbuf[3*k+2] = colormap[v];  
    }
#else
    v=NUMVALS+MAPOFFSET;
    vidbuf[k] = v; 
#endif
          k++;
        }
        else if (overlay_mode==MOVEABLE)
        {
          if (imnr_2>=0&&imnr_2<imnr1_2&&
              i_2>=0&&i_2<ydim_2&&j_2>=0&&j_2<xdim_2)
          {
            f = fscale_2*fim_2[imnr_2][ydim_2-1-i_2][j_2];
            v = ((f<0)?0:(f>255)?255:f)+MAPOFFSET;
          }
          else v=MAPOFFSET;
          if (imnr==impt||i==ipt) v=255-(v-MAPOFFSET)+MAPOFFSET;
/*
          if (imnr==imc&&i==ic)
             printf("imnr_2 = %d (%f,%f,%f)\n",imnr_2,x_2,y_2,z_2);
*/
          if ((imnr==imc&&abs(i-ic)<=2)||
              (i==ic&&abs(imnr-imc)<=2)) 
#ifdef RGB
    {
      lvidbuf[3*k] = 255;
      lvidbuf[3*k+1] = 0;
      lvidbuf[3*k+2] = 0;
    } else {
      lvidbuf[3*k] = colormap[v];
      lvidbuf[3*k+1] = colormap[v];
      lvidbuf[3*k+2] = colormap[v];  
    }
#else
    v=NUMVALS+MAPOFFSET;
    vidbuf[k] = v; 
#endif  
          k++;
        } 
      }
    }
  }
#ifdef RGB
  
#else
  rectwrite(0,0,xdim-1,ydim-1,vidbuf);
#endif
  swapbuffers();
  invalid_buffers=1;
}

void blur(float factor)  /* test hack */
{
  FILE *fp;
  long xorig,xsize,yorig,ysize;
  int x0,y0,x1,y1;
  int i,j,k;
  /* short r; */
  char command[2*NAME_LENGTH];

  getorigin(&xorig,&yorig);
  getsize(&xsize,&ysize);
  x0 = (int)xorig;  x1 = (int)(xorig+xsize-1);
  y0 = (int)yorig;  y1 = (int)(yorig+ysize-1);

  sprintf(command,"scrsave /tmp/tmp1.rgb %d %d %d %d\n",x0,x1,y0,y1);
  system(command);
  sprintf(command,"blur /tmp/tmp1.rgb /tmp/tmp2.rgb %d\n",
                                            (int)((float)xsize/factor));
  system(command);
  sprintf(command,"tobin /tmp/tmp2.rgb /tmp/tmp2.bin\n");
  system(command);

  fp = fopen("/tmp/tmp2.bin","r");
  fread(binbuff,3,xdim*ydim,fp);
  k = 0;
  for (i=0;i<ydim;i++)
  for (j=0;j<xdim;j++) {
    /*vidbuf[k] = binbuff[(i*xdim+j)*3+0]+MAPOFFSET;*/  /* x9 */

#ifdef RGB
    vidbuf[3*k] = binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
    vidbuf[3*k+1] = binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
    vidbuf[3*k+2] = binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
#else
    vidbuf[k] = binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
#endif    
    k++;
  }
  backbuffer(FALSE);
  frontbuffer(TRUE);
  rectwrite(0,0,xdim-1,ydim-1,vidbuf);
  backbuffer(TRUE);
  frontbuffer(FALSE);
  sprintf(command,"rm -f /tmp/tmp1.rgb /tmp/tmp2.rgb /tmp/tmp2.bin\n");
  system(command);
}

void make_filenames(char *lsubjectsdir)
{
    subjectsdir = (char *)malloc(NAME_LENGTH*sizeof(char)); /* malloc for tcl */
    srname = (char *)malloc(NAME_LENGTH*sizeof(char));
    psrname = (char *)malloc(NAME_LENGTH*sizeof(char));
    pname = (char *)malloc(NAME_LENGTH*sizeof(char));
    regfname = (char *)malloc(NAME_LENGTH*sizeof(char));
    afname = (char *)malloc(NAME_LENGTH*sizeof(char));
    targpref = (char *)malloc(NAME_LENGTH*sizeof(char));
    movformat = (char *)malloc(NAME_LENGTH*sizeof(char));
    tfname = (char *)malloc(NAME_LENGTH*sizeof(char));
    sgfname = (char *)malloc(NAME_LENGTH*sizeof(char));

    strcpy(subjectsdir,lsubjectsdir);
    strcpy(movformat,"");
    /* TODO: init others here */
}

void read_float_images(float ***fim,char *format,int nslices,int nperslice,int xdim,int ydim,short **buf) /* was MRIio.c */
{
#ifdef Linux
  int scount ;
#endif
  int i,j,k,n,t,imnr,imtype,bufsize;
  char fname[NAME_LENGTH],*ext,*BRIKext;
  float f;
  short s;
  FILE *fp = NULL;
  long offset=0;
  /* next four for AFNI BRIKs only */
  char hname[STRLEN];
  float ffact;
  int btype;
  char byteorder[STRLEN];
  int npix=0;
  double sum=0,sum2=0,avg,stdev;
  float max=-BIGFLOAT,min=BIGFLOAT;


  ext = &format[strlen(format)-1];
  while (ext>format && ext[0]!='.') ext--;
  if (!strcmp(ext,".im")) {
    imtype = APD1;
    printf("register: input file type: im (x,y: <pref><z*tcnt+t>.im)\n");
  }
  else if (!strcmp(ext,".bshort")) {
    imtype = BSHORT;
    printf("register: input file type: bshort (x,y,t: <prefix><z>.bshort)\n");
  }
  else if (!strcmp(ext,".skip")) {
    imtype = SKIP;
    printf("register: input file type: skip (x,y: <prefix>.<z>.<t>.skip)\n");
  }
  else if (!strcmp(ext,".BRIK")) {
    imtype = AFNI;
    printf("register: input file type: AFNI (x,y,z,t: <prefix>.BRIK)\n");
    strcpy(hname,format);
    /* replace .BRIK ext with .HEAD */
    BRIKext = strstr(hname,".BRIK");
    sprintf(BRIKext,".HEAD");
    /* read HEAD file to get floatfac, briktype, and byteorder */
    read_AFNI_HEAD(hname,&ffact,&btype,byteorder);
  }
  else if (!strcmp(ext,".brik")) {
    imtype = AFNI;
    printf("register: input file type: AFNI (x,y,z,t: <prefix>.brik)\n");
  }
  else {
    printf("register: ### file format %s (%s) not supported\n",format,ext);
    exit(1); }

  for (k=0;k<nslices;k++) {
    fim[k] = (float **)calloc(ydim,sizeof(float *));
    for (i=0;i<ydim;i++) {
      fim[k][i] = (float *)calloc(xdim,sizeof(float));
    }
  }
  
  t = 0;
  bufsize = xdim*ydim;
  if (imtype==AFNI) { /* one file per scan (x,y,z,t)  */
    sprintf(fname,format);
    fp = fopen(fname,"r");
    if (fp==NULL)
      ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,fname);
  }
  if (imtype==AFNI && MATCH(ext,".BRIK")) {
    for (k=0;k<nslices;k++)
    for (i=0;i<ydim;i++)
    for (j=0;j<xdim;j++) {
      if (btype==BRIK_SHORT) {
        s = freadShort(fp);
        if(MATCH(byteorder,"LSB")) s = swapShort(s);
        f = fim[k][i][j] = ffact*(float)s;
      } else if (btype==BRIK_FLOAT) {
        f = freadFloat(fp);
        if(MATCH(byteorder,"LSB")) f = swapFloat(f);
        f = fim[k][i][j] = ffact*f;
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
    fclose(fp);
  } else {
    if ((*buf)==NULL)
      *buf = (short *)calloc(bufsize,sizeof(float));
    for (k=0;k<nslices;k++) {
      if (imtype==APD1) { /* x,y file skip header */
        imnr = k*nperslice + t;
        sprintf(fname,format,imnr);
        fp = fopen(fname,"r");
        if (fp==NULL){printf("register: ### File %s not found\n",fname); return;}
        fread(*buf,sizeof(char),HEADERSIZE,fp);
      }
      if (imtype==BSHORT) { /* one file per slice (x,y,t) */
        imnr = k;
        sprintf(fname,format,imnr);
        fp = fopen(fname,"r");
        if (fp==NULL){printf("register: ### File %s not found\n",fname); return;}
      }
      if (imtype==SKIP) { /* x,y file */
        sprintf(fname,format,k,t);   /* assumes zero-based */
        fp = fopen(fname,"r");
        if (fp==NULL){printf("register: ### File %s not found\n",fname); return;}
      }
      if (imtype==AFNI) { /* brick: skip t*zcnt + z */
        offset = (t*nslices + k) * bufsize * sizeof(short);
        fseek(fp,offset,SEEK_SET);
      }
      if (k==0) printf("First slice read: %s\n",fname); fflush(stdout);
      fread(*buf,sizeof(short),bufsize,fp);
  #ifdef Linux
      for(scount=0; scount<bufsize; scount++) {
        (*buf)[scount]=swapShort((*buf)[scount]);
      }
  #endif
      if (imtype==APD1 || imtype==BSHORT || imtype==SKIP) fclose(fp);
      max = -1.0e10;
      min = 1.0e10;
      sum = sum2 = 0;
      n = 0;
      for (i=0;i<ydim;i++)
      for (j=0;j<xdim;j++) {
        f = fim[k][i][j] = (*buf)[n++];
        if (f<min) min = f;
        if (f>max) max = f;
        sum += f;
        sum2 += f*f;
      }
      sum /= xdim*ydim;
      sum2 = sqrt(sum2/(xdim*ydim)-sum*sum);
      /*printf("File %s read\n",fname);*/
      /*printf("min=%f,max=%f,avg=%f,stdev=%f\n",min,max,sum,sum2);*/
    }
    if (imtype==AFNI) fclose(fp);
    printf("Last slice read:  %s\n",fname); fflush(stdout);
  }
}

/* boilerplate wrap function defines for easier viewing */
#define WBEGIN (ClientData clientData,Tcl_Interp *interp,int argc,char *argv[]){
#define ERR(N,S)  if(argc!=N){interp->result=S;return TCL_ERROR;}
#define WEND   return TCL_OK;}
#define REND  (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL

/*=======================================================================*/
/* function wrappers and errors */
int                  W_move_window  WBEGIN
  ERR(3,"Wrong # args: move_window <x> <y>")
                       move_window(atoi(argv[1]),atoi(argv[2]));  WEND

int                  W_pop_gl_window  WBEGIN 
  ERR(1,"Wrong # args: pop_gl_window")
                       pop_gl_window();  WEND

int                  W_redraw  WBEGIN 
  ERR(1,"Wrong # args: redraw")
                       redraw();  WEND

int                  W_upslice WBEGIN 
  ERR(1,"Wrong # args: upslice")
                       upslice();  WEND

int                  W_downslice WBEGIN 
  ERR(1,"Wrong # args: downslice")
                       downslice();  WEND

int                  W_rotate_brain_x  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_x <deg>")
                       rotate_brain(atof(argv[1]),'x'); WEND

int                  W_rotate_brain_y  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_y <deg>")
                       rotate_brain(atof(argv[1]),'y'); WEND

int                  W_rotate_brain_z  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_z <deg>")
                       rotate_brain(atof(argv[1]),'z'); WEND

int                  W_translate_brain_x  WBEGIN
  ERR(2,"Wrong # args: translate_brain_x <mm>")
                       translate_brain(atof(argv[1]),'x');  WEND

int                  W_translate_brain_y  WBEGIN
  ERR(2,"Wrong # args: translate_brain_y <mm>")
                       translate_brain(atof(argv[1]),'y');  WEND

int                  W_translate_brain_z  WBEGIN
  ERR(2,"Wrong # args: translate_brain_z <mm>")
                       translate_brain(atof(argv[1]),'z');  WEND

int                  W_scale_brain_x  WBEGIN
  ERR(2,"Wrong # args: scale_brain_x <mm>")
                       scale_brain(atof(argv[1]),'x');  WEND

int                  W_scale_brain_y  WBEGIN
  ERR(2,"Wrong # args: scale_brain_y <mm>")
                       scale_brain(atof(argv[1]),'y');  WEND

int                  W_scale_brain_z  WBEGIN
  ERR(2,"Wrong # args: scale_brain_z <mm>")
                       scale_brain(atof(argv[1]),'z');  WEND

int                  W_goto_point  WBEGIN
  ERR(1,"Wrong # args: goto_point")
                       goto_point(tfname);  WEND

int                  W_write_point  WBEGIN
  ERR(1,"Wrong # args: write_point")
                       write_point(tfname);  WEND

int                  W_save_rgb  WBEGIN
  ERR(1,"Wrong # args: save_rgb")
                       save_rgb(sgfname);  WEND

int                  W_record_swapbuffers  WBEGIN
  ERR(1,"Wrong # args: record_swapbuffers")
                       record_swapbuffers();  WEND

int                  W_set_scale  WBEGIN
  ERR(1,"Wrong # args: set_scale")
                       set_scale();  WEND

int                  W_blur  WBEGIN
  ERR(2,"Wrong # args: blur <factor>")
                       blur(atof(argv[1]));  WEND

int                  W_read_reg  WBEGIN
  ERR(1,"Wrong # args: read_reg")
                       read_reg(regfname);  WEND

int                  W_write_reg  WBEGIN
  ERR(1,"Wrong # args: write_reg")
                       write_reg(regfname);  WEND

int                  W_align_points  WBEGIN
  ERR(1,"Wrong # args: align_points")
                       align_points();  WEND

int                  W_mirror_brain  WBEGIN
  ERR(1,"Wrong # args: mirror_brain")
                       mirror_brain();  WEND
/*=======================================================================*/

/* for tcl/tk */
static void StdinProc _ANSI_ARGS_((ClientData clientData, int mask));
static void Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
#ifndef TCL8
static Tk_Window mainWindow;
#endif
static Tcl_Interp *interp;
static Tcl_DString command;
static int tty;

int main(argc, argv)   /* new main */
int argc;
char **argv;
{
  int code;
#ifndef TCL8
  static char *display = NULL;
#endif
  char tkregister_tcl[NAME_LENGTH];
  /* char str[NAME_LENGTH]; */
  char *envptr;
  FILE *fp ;

  initcolormap();

  /* get tkregister tcl startup script location from environment */
  envptr = getenv("MRI_DIR");
  if (envptr==NULL) {
    printf("tkregister: env var MRI_DIR undefined (use setenv)\n");
    printf("    [dir containing mri distribution]\n");
    exit(1);
  }
  sprintf(tkregister_tcl,"%s/lib/tcl/%s",envptr,"tkregister.tcl");
  if ((fp=fopen(tkregister_tcl,"r"))==NULL) {
    printf("tkregister: startup script %s not found\n",tkregister_tcl);
    exit(1);
  }
  else fclose(fp);

  /* start program, now as function; gl window not opened yet */
  printf("tkregister: starting register\n");
  Register((ClientData) NULL, interp, argc, argv);/* event loop commented out*/

  /* start tcl/tk; first make interpreter */
  interp = Tcl_CreateInterp();

  /* make main window (not displayed until event loop starts) */

#ifdef TCL8
  /* Tk_Init(interp); */
#else
  mainWindow = Tk_CreateMainWindow(interp, display, argv[0], "Tk"); 
#endif

  /*
  if (mainWindow == NULL) {
    fprintf(stderr, "%s\n", interp->result);
    exit(1); }
  */

  /* set the "tcl_interactive" variable */
  tty = isatty(0);
  Tcl_SetVar(interp, "tcl_interactive", (tty) ? "1" : "0", TCL_GLOBAL_ONLY);
  if (tty) promptflag = TRUE;

  /* read tcl/tk internal startup scripts */
  if (Tcl_Init(interp) == TCL_ERROR) {
    fprintf(stderr, "Tcl_Init failed: %s\n", interp->result); }
  if (Tk_Init(interp)== TCL_ERROR) {
    fprintf(stderr, "Tk_Init failed: %s\n", interp->result); }

  /*=======================================================================*/
  /* register wrapped surfer functions with interpreter */
  Tcl_CreateCommand(interp, "redraw",             W_redraw,             REND);
  Tcl_CreateCommand(interp, "move_window",        W_move_window,        REND);
  Tcl_CreateCommand(interp, "pop_gl_window",      W_pop_gl_window,      REND);
  Tcl_CreateCommand(interp, "upslice",            W_upslice,            REND);
  Tcl_CreateCommand(interp, "downslice",          W_downslice,          REND);
  Tcl_CreateCommand(interp, "rotate_brain_x",     W_rotate_brain_x,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_y",     W_rotate_brain_y,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_z",     W_rotate_brain_z,     REND);
  Tcl_CreateCommand(interp, "translate_brain_x",  W_translate_brain_x,  REND);
  Tcl_CreateCommand(interp, "translate_brain_y",  W_translate_brain_y,  REND);
  Tcl_CreateCommand(interp, "translate_brain_z",  W_translate_brain_z,  REND);
  Tcl_CreateCommand(interp, "scale_brain_x",      W_scale_brain_x,      REND);
  Tcl_CreateCommand(interp, "scale_brain_y",      W_scale_brain_y,      REND);
  Tcl_CreateCommand(interp, "scale_brain_z",      W_scale_brain_z,      REND);
  Tcl_CreateCommand(interp, "goto_point",         W_goto_point,         REND);
  Tcl_CreateCommand(interp, "write_point",        W_write_point,        REND);
  Tcl_CreateCommand(interp, "save_rgb",           W_save_rgb,           REND);
  Tcl_CreateCommand(interp, "record_swapbuffers", W_record_swapbuffers, REND);
  Tcl_CreateCommand(interp, "set_scale",          W_set_scale,          REND);
  Tcl_CreateCommand(interp, "blur",               W_blur,               REND);
  Tcl_CreateCommand(interp, "read_reg",           W_read_reg,           REND);
  Tcl_CreateCommand(interp, "write_reg",          W_write_reg,          REND);
  Tcl_CreateCommand(interp, "align_points",       W_align_points,       REND);
  Tcl_CreateCommand(interp, "mirror_brain",       W_mirror_brain,       REND);
  /*=======================================================================*/
  /***** link global BOOLEAN variables to tcl equivalents */
  Tcl_LinkVar(interp,"maxflag",(char *)&maxflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"updateflag",(char *)&updateflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"blinkflag",(char *)&blinkflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"editedmatrix",(char *)&editedmatrix, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"maskflag",(char *)&maskflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"promptflag",(char *)&promptflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"followglwinflag",(char *)&followglwinflag,
                                                        TCL_LINK_BOOLEAN);
  /*=======================================================================*/
  /***** link global INT variables to tcl equivalents */
  Tcl_LinkVar(interp,"zf",(char *)&zf, TCL_LINK_INT);
  Tcl_LinkVar(interp,"xdim",(char *)&xdim, TCL_LINK_INT);
  Tcl_LinkVar(interp,"ydim",(char *)&ydim, TCL_LINK_INT);
  Tcl_LinkVar(interp,"imnr1",(char *)&imnr1, TCL_LINK_INT);
  Tcl_LinkVar(interp,"plane",(char *)&plane, TCL_LINK_INT);
  Tcl_LinkVar(interp,"imc",(char *)&imc, TCL_LINK_INT);
  Tcl_LinkVar(interp,"ic",(char *)&ic, TCL_LINK_INT);
  Tcl_LinkVar(interp,"jc",(char *)&jc, TCL_LINK_INT);
  Tcl_LinkVar(interp,"prad",(char *)&prad, TCL_LINK_INT);
  Tcl_LinkVar(interp,"blinktop",(char *)&blinktop, TCL_LINK_INT);
  Tcl_LinkVar(interp,"overlay_mode",(char *)&overlay_mode, TCL_LINK_INT);
  Tcl_LinkVar(interp,"visible_mode",(char *)&visible_mode, TCL_LINK_INT);
  Tcl_LinkVar(interp,"last_visible_mode",(char *)&last_visible_mode, 
                                                                TCL_LINK_INT);
  Tcl_LinkVar(interp,"visible_plane",(char *)&visible_plane, TCL_LINK_INT);
  Tcl_LinkVar(interp,"last_visible_plane",(char *)&last_visible_plane, 
                                                                TCL_LINK_INT);
  Tcl_LinkVar(interp,"blinkdelay",(char *)&blinkdelay, TCL_LINK_INT);
  Tcl_LinkVar(interp,"blinktime",(char *)&blinktime, TCL_LINK_INT);
  /*=======================================================================*/
  /***** link global DOUBLE variables to tcl equivalents (were float) */
  Tcl_LinkVar(interp,"fsquash",(char *)&fsquash, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fthresh",(char *)&fthresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fscale_2",(char *)&fscale_2, TCL_LINK_DOUBLE);
  /*=======================================================================*/
  /***** link global malloced STRING vars */
  Tcl_LinkVar(interp,"home",        (char *)&subjectsdir,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"session",     (char *)&srname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subject",     (char *)&pname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"registerdat", (char *)&regfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"analysedat",  (char *)&afname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subjtmpdir",  (char *)&tfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"rgb",         (char *)&sgfname,TCL_LINK_STRING);
  /*=======================================================================*/

  /* run tcl/tk startup script to set vars, make interface; no display yet */
  printf("tkregister: interface: %s\n",tkregister_tcl);
  code = Tcl_EvalFile(interp,tkregister_tcl);
  if (*interp->result != 0)  printf(interp->result);

  /* always start up command line shell too */
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  if (tty)
    Prompt(interp, 0);
  fflush(stdout);
  Tcl_DStringInit(&command);
  Tcl_ResetResult(interp);

  /*Tk_MainLoop();*/  /* standard */

  /* dual event loop (interface window made now) */
#ifdef TCL8
  while(Tk_GetNumMainWindows() > 0) {
#else
  while(tk_NumMainWindows > 0) {
#endif
    while (Tk_DoOneEvent(TK_ALL_EVENTS|TK_DONT_WAIT)) {
      /* do all the tk events; non-blocking */
    }
    do_one_gl_event(interp);
    /*sginap((long)1);*/   /* block for 10 msec */
    usecnap(10000);     /* block for 10 msec */
  }

  Tcl_Eval(interp, "exit");
  exit(0);
}

void usecnap(int usec)
{
  struct timeval delay;

  delay.tv_sec = 0;
  delay.tv_usec = (long)usec;
  select(0,NULL,NULL,NULL,&delay);
}

void initcolormap()
{
  int i;
  for(i=0; i<256; i++)
    colormap[i]=i;
  for(i=256; i<512; i++)
    colormap[i]=255; 
}

void pseudo_swapbuffers()
{
  if(pswapnext==MOVEABLE) {
    rectwrite(0,0,xdim-1,ydim-1,blinkbufm);
    pswapnext=TARGET;
  } else {
    rectwrite(0,0,xdim-1,ydim-1,blinkbuft);
    pswapnext=MOVEABLE;
  }
}
/*=== from TkMain.c ===================================================*/
static void StdinProc(clientData, mask)
  ClientData clientData;
  int mask;
{
#define BUFFER_SIZE 4000
  char input[BUFFER_SIZE+1];
  static int gotPartial = 0;
  char *cmd;
  int code, count;

  count = read(fileno(stdin), input, BUFFER_SIZE);
  if (count <= 0) {
    if (!gotPartial) {
      if (tty) {Tcl_Eval(interp, "exit"); exit(1);}
      else     {Tk_DeleteFileHandler(0);}
      return;
    }
    else count = 0;
  }
  cmd = Tcl_DStringAppend(&command, input, count);
  if (count != 0) {
    if ((input[count-1] != '\n') && (input[count-1] != ';')) {
      gotPartial = 1;
      goto prompt; }
    if (!Tcl_CommandComplete(cmd)) {
      gotPartial = 1;
      goto prompt; }
  }
  gotPartial = 0;
  Tk_CreateFileHandler(0, 0, StdinProc, (ClientData) 0);
  code = Tcl_RecordAndEval(interp, cmd, TCL_EVAL_GLOBAL);
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  Tcl_DStringFree(&command);
  if (*interp->result != 0)
    if ((code != TCL_OK) || (tty))
      puts(interp->result);
  prompt:
  if (tty)  Prompt(interp, gotPartial);
  Tcl_ResetResult(interp);
}

/*=== from TkMain.c ===================================================*/
static void Prompt(interp, partial)
  Tcl_Interp *interp;
  int partial;
{
  char *promptCmd;
  int code;

  promptCmd = Tcl_GetVar(interp,
      partial ? "tcl_prompt2" : "tcl_prompt1", TCL_GLOBAL_ONLY);
  if (promptCmd == NULL) {
    defaultPrompt:
    if (!partial)
      fputs("% ", stdout);
  }
  else {
    code = Tcl_Eval(interp, promptCmd);
    if (code != TCL_OK) {
      Tcl_AddErrorInfo(interp,
                    "\n    (script that generates prompt)");
      fprintf(stderr, "%s\n", interp->result);
      goto defaultPrompt;
    }
  }
  fflush(stdout);
}
#if 0
int __new_cfgetospeed () {
  fprintf( stderr, "__new_cfgetospeed\n" );
}
int __new_cfsetospeed () {
  fprintf( stderr, "__new_cfsetospeed\n" );
}
int __new_cfsetispeed () {
  fprintf( stderr, "__new_cfsetispeed\n" );
}
int __new_tcgetattr () {
  fprintf( stderr, "__new_tcgetattr\n" );
}
int __new_tcsetattr () {
  fprintf( stderr, "__new_tcsetattr\n" );
}
#endif


#ifdef USE_mriTransform

static void InitFunctionalTransform ( char* isRegisterFile ) {

  fMRI_REG* regInfo          = NULL;
  MATRIX*   mTmp             = NULL;
  float     ps               = 0;
  float     st               = 0;
  float     slices           = 0;
  float     rows             = 0;
  float     cols             = 0;

  /* allocate the transform */
  Trns_New( &gRASToIdxTransform );

  /* read the registration info */
  regInfo = StatReadRegistration ( fname );
  if ( NULL == regInfo ) {
    ErrorExit( ERROR_NO_FILE, "Couldn't find registration %s", fname );
  }

  /* get our stats */
  ps     = regInfo->in_plane_res;
  st     = regInfo->slice_thickness;
  slices = nslices ;
  rows   = ydim_2 ;
  cols   = xdim_2 ;
  
  // create the functional index to functional ras matrix
  mTmp = MatrixAlloc( 4, 4, MATRIX_REAL );
  MatrixClear( mTmp );
  *MATRIX_RELT(mTmp,1,1) = -ps;
  *MATRIX_RELT(mTmp,2,3) = st;
  *MATRIX_RELT(mTmp,3,2) = -ps;
  *MATRIX_RELT(mTmp,1,4) = (ps*cols) / 2.0;
  *MATRIX_RELT(mTmp,2,4) = -(st*slices) / 2.0;
  *MATRIX_RELT(mTmp,3,4) = (ps*rows) / 2.0;
  *MATRIX_RELT(mTmp,4,4) = 1.0;
  Trns_CopyBtoRAS( gRASToIdxTransform, mTmp );

  /* the ras to idx transformer gets the identity matrix here */
  MatrixIdentity( 4, mTmp );
  Trns_CopyAtoRAS( gRASToIdxTransform, mTmp );

  /* a is anatomical, b is functional, mri2fmri takes us from a_ras
     to b_ras */
  Trns_CopyARAStoBRAS( gRASToIdxTransform, regInfo->mri2fmri );

  // get rid of the registration info.
  StatFreeRegistration ( &regInfo );

  MatrixFree( &mTmp );
}

xVoxel sCoord1, sCoord2;

static int  ConvertRASToFuncIdx  ( float  ifRASX, 
           float  ifRASY, 
           float  ifRASZ,
           float* ofFuncIdxX, 
           float* ofFuncIdxY, 
           float* ofFuncIdxZ ) {

  Trns_tErr eTransform = Trns_tErr_NoErr;

  /* set coord */
  xVoxl_SetFloat( &sCoord1, ifRASX, ifRASY, ifRASZ );
  
  /* do the conversion */
  eTransform = Trns_ConvertAtoB( RASToIdxTransform, 
                                 &sCoord1, &sCoord2 );
  if( Trns_tErr_NoErr != eTransform ) {
    ErrorExit( ERROR_BADPARM, "Error %d in Trns_ConvertAtoB: %s\n",
               eTransform, Trns_GetErrorString( eTransform ) );
    return NO_ERROR;
  }

  /* return the coord */
  *ofFuncIdxX = xVoxl_GetFloatX( &sCoord2 );
  *ofFuncIdxY = xVoxl_GetFloatY( &sCoord2 );
  *ofFuncIdxZ = xVoxl_GetFloatZ( &sCoord2 );

  return(NO_ERROR) ;
}

#endif /* USE_mriTransform */

void read_AFNI_HEAD(char *fname, float *ffact, int *btype, char *byteorder)
{
  FILE *fp;
  char line[STRLEN];
  char entry[STRLEN];
  char value[STRLEN];
  int tsize;

  *ffact=-1;
  *btype=-1;
  sprintf(byteorder,"UNKNOWN");

  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### file %s not found\n",progname,fname);
  while (fgets(line,99,fp) != NULL) {
    sscanf(line,"%s %*s %s",entry,value);
    if (MATCH(value,"BRICK_FLOAT_FACS")) {
      fgets(line,99,fp);
      sscanf(line,"%s %*s %s",entry,value);
      tsize = atoi(value);
      /* read first float fac */
      fscanf(fp,"%f",ffact);
      MsgPrintf("%s: first of %d float_facs in %s: %1.4f\n",
                progname,tsize,fname,*ffact);
    }
    if (MATCH(value,"BRICK_TYPES")) {
      fgets(line,99,fp);
      sscanf(line,"%s %*s %s",entry,value);
      tsize = atoi(value);
      fscanf(fp,"%d",btype);
      MsgPrintf("%s: first of %d brick_types in %s: %d\n",
                progname,tsize,fname,*btype);
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

  switch (*btype) {
    case BRIK_SHORT:
      MsgPrintf("%s: brik data type for %s: short\n",progname,fname);
      break;
    case BRIK_FLOAT:
      MsgPrintf("%s: brik data type for %s: float\n",progname,fname);
      break;
    default:
      ErrorExit(ERROR_BADFILE,"%s: ### brik type %d is currently unsupported ...quitting\n",
              progname,btype);
  }
  if(*ffact==0)
    *ffact=1;
  else if(*ffact==-1)
    ErrorExit(ERROR_BADFILE,"%s: ### could not read float factor from %s\n",progname,fname);
  if(MATCH(byteorder,"UNKNOWN"))
    ErrorExit(ERROR_BADFILE,"%s: ### could not read byte order from %s\n",progname,fname);

  fclose(fp);
}


