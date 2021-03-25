/* fieldsign.c: read pol and ecc stats files and calculate field sign for each vertex
      created: 02/28/05 DH
     last mod: 02/16/06 DH

   purpose:
     calculating field sign of surface stats

   input:
     complex stats file (w file)

   output:
     complex stats file (w file)

   acknowledgements:
     fieldsign code from tksurfer
*/

#include "surflib.h"

#define MINARGC 7
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char polstem[STRLEN]=UNDEFSTR;
char eccstem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]="fieldsign";
char poldir[STRLEN]=".";
char eccdir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char *real_infix=NULL;
char *imag_infix=NULL;
int smooth = 50; /* pre-smoothing steps */
int revfsflag = 0;
char surf[STRLEN]="smoothwm";

/* functions */

void usage()
{
  printf("\n");
  printf("Usage: %s -polstem polstem -eccstem eccstem -name subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -polstem polstem   omit extension, infixes, hemi: <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -eccstem eccstem   omit extension, infixes, hemi: <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -name    subjname          subject name\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem   [fieldsign]     output file stem\n");
  printf("    -poldir    [.]             polar angle input dir\n");
  printf("    -eccdir    [.]             eccentricity input dir\n");
  printf("    -outdir    [.]             output dir\n");
  printf("    -hemi      [rh]            hemisphere (rh or lh)\n");
  printf("    -surf      [smoothwm]      surface used for measurements\n");
  printf("    -smooth    [50]            number of pre-smoothing steps\n");
  printf("    -infixes   [_r _i]         real,imaginary infixes\n");
  printf("    -revfs                     reverse fieldsign\n");
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

  real_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(real_infix);
  imag_infix  = (char *)malloc(STRLEN*sizeof(char));      MTEST(imag_infix);
  strcpy(real_infix,REAL_INFIX);
  strcpy(imag_infix,IMAG_INFIX);

  /* parse arguments */
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-polstem") && i+1<argc) {
        strcpy(polstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccstem") && i+1<argc) {
        strcpy(eccstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-poldir") && i+1<argc) {
        strcpy(poldir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccdir") && i+1<argc) {
        strcpy(eccdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-name") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-smooth") && i+1<argc) {
        smooth = atoi(argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-infixes")) && i+2<argc){
        strcpy(real_infix,argv[i+1]); strcpy(imag_infix,argv[i+2]); i+=2;
      } else
      if (MATCH(argv[i],"-revfs")){
        revfsflag = 1;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
      {
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
      }
    } else {
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    }
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(polstem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### polstem not specified ...quitting\n",
              progname);
  }
  if (MATCH(eccstem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccstem not specified ...quitting\n",
              progname);
  }
  if (!FileExists(poldir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### poldir %s not found ...quitting\n",
              progname,poldir);
  }
  if (!FileExists(eccdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccdir %s not found ...quitting\n",
              progname,eccdir);
  }
  if(!isadir(poldir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,poldir);
  }
  if(!isadir(eccdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,eccdir);
  }
  if (!FileExists(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  }
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (smooth<0) {
    ErrorExit(ERROR_BADPARM,"%s: ### smooth steps must be >= 0 ...quitting\n",
              progname);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}


void complex2polar(MRIS *mris)
{
  VERTEX *v;
  int k;
  float ecc_phase,ecc_amp,pol_phase,pol_amp;

  for (k=0;k<mris->nvertices;k++) {
    v = &mris->vertices[k];
    pol_phase = atan2(v->imag_val,v->val);
    pol_amp = hypot(v->imag_val,v->val);
    ecc_phase = atan2(v->val2bak,v->valbak);
    ecc_amp = hypot(v->val2bak,v->valbak);
    v->val = pol_phase;
    v->imag_val = pol_amp;
    v->valbak = ecc_phase;
    v->val2bak = ecc_amp;
  }
}

void polar2complex(MRIS *mris)
{
  VERTEX *v;
  int k;
  float ecc_phase,ecc_amp,pol_phase,pol_amp;
  
  for (k=0;k<mris->nvertices;k++) {
    v = &mris->vertices[k];
    pol_phase = v->val;
    pol_amp = v->imag_val;
    ecc_phase = v->valbak;
    ecc_amp = v->val2bak;
    v->val = pol_amp*cos(pol_phase);
    v->imag_val = pol_amp*sin(pol_phase);
    v->valbak = ecc_amp*cos(ecc_phase);
    v->val2bak = ecc_amp*sin(ecc_phase);
  }
}

void swapValues(MRIS *mris)
{
  VERTEX *v;
  int k;
  float temp;

  for (k=0;k<mris->nvertices;k++) {
    v = &mris->vertices[k];

    temp = v->val;
    v->val = v->valbak;
    v->valbak = temp;

    temp = v->imag_val;
    v->imag_val = v->val2bak;
    v->val2bak = temp;      
  }
}

float circsubtract(float a,float b)
{
  float h = a-b;
  if (h<-M_PI) h = h+2*M_PI;
  else if (h>M_PI) h = h-2*M_PI;
  return h;
}

float distance(float x1, float y1, float z1, float x2, float y2, float z2)
{
  float dist;
  
  dist = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
  
  return dist;
}

void crossprod(float x1, float y1, float z1,
               float x2, float y2, float z2,
               float *xp, float *yp, float *zp)
{
  *xp = y1*z2 - z1*y2;
  *yp = z1*x2 - x1*z2;
  *zp = x1*y2 - y1*x2;
}

float dotprod(float x1, float y1, float z1,
              float x2, float y2, float z2)
{
  float dot;
  
  dot = x1*x2 + y1*y2 + z1*z2;
  
  return dot;
}

float vectangle(float x1, float y1, float z1,
                float x2, float y2, float z2)
{
  float a,b,c,theta;

  a = distance(x1,y1,z1,x2,y2,z2);
  b = distance(0,0,0,x1,y1,z1);
  c = distance(0,0,0,x2,y2,z2);
  theta = acos((b*b + c*c - a*a)/(2*b*c));

  return theta;
}

void calcFieldSign(MRIS *mris)
{
  int j,k,m,p,q;
  VERTEX *v;
  float dv1,dv2,dx,dy,dv1dx,dv1dy,dv2dx,dv2dy;
  float m11,m12,m13,m22,m23,z1,z2,z3,z1b,z2b,z3b,denom;
  int num_nbrs,newnbr,newnbr2;
  int *nbrs,*oldnbrs,*nbrnums;
  float *x,*y,*z,*theta,*fx,*fy;
  float dist,dotp,xp,yp,zp,sumtheta,th,mintheta,maxtheta;
 
  for (k=0;k<mris->nvertices;k++) {
    v = &mris->vertices[k];
    num_nbrs = v->vnum;
    if(num_nbrs <= 1) continue;
    x=      (float *)calloc(num_nbrs,sizeof(float));  MTEST(x);
    y=      (float *)calloc(num_nbrs,sizeof(float));  MTEST(y);
    z=      (float *)calloc(num_nbrs,sizeof(float));  MTEST(z);
    theta=  (float *)calloc(num_nbrs,sizeof(float));  MTEST(theta);
    fx=     (float *)calloc(num_nbrs,sizeof(float));  MTEST(fx);
    fy=     (float *)calloc(num_nbrs,sizeof(float));  MTEST(fy);
    oldnbrs=(int   *)calloc(num_nbrs,sizeof(int));    MTEST(oldnbrs);
    nbrs=   (int   *)calloc(num_nbrs,sizeof(int));    MTEST(nbrs);
    nbrnums=(int   *)calloc(num_nbrs,sizeof(int));    MTEST(nbrnums);

    /* get x,y,z coordinates for each neighbor, offset by center vertex */
    for (m=0;m<num_nbrs;m++) {
      j = oldnbrs[m] = v->v[m];
      x[m] = mris->vertices[j].x - v->x;
      y[m] = mris->vertices[j].y - v->y;
      z[m] = mris->vertices[j].z - v->z;
    }

    for (q=0;q<num_nbrs;q++) {
      for (m=0;m<num_nbrs;m++) {
        crossprod(x[q],y[q],z[q],x[m],y[m],z[m],&xp,&yp,&zp);
        dist = distance(0,0,0,xp,yp,zp);
      }
    }

    /* reorder neighbors counterclockwise -- originally random */
    nbrs[0]=oldnbrs[0];
    nbrnums[0]=0;
    oldnbrs[0]=-1;
    for (p=1;p<num_nbrs;p++) {
      q=nbrnums[p-1];
      mintheta=2*M_PI;
      maxtheta=-2*M_PI;
      newnbr=-1;
      newnbr2=-1;
      for (m=1;m<num_nbrs;m++) {
        /* skip if this neighbor has already been reassigned */
        if(oldnbrs[m]<0) continue;
      
        /* calculate cross product between two neighbor vectors */
        crossprod(x[q],y[q],z[q],x[m],y[m],z[m],&xp,&yp,&zp);
        /* calculate angle between two neighbor vectors */
        th=vectangle(x[q],y[q],z[q],x[m],y[m],z[m]);
        dist=distance(x[q],y[q],z[q],x[m],y[m],z[m]);
        dotp = dotprod(v->nx,v->ny,v->nz,xp,yp,zp);
        /* if cross prod points out of surface and this is smallest 
           angle between neighbors, pick this as next neighbor */
        if(dotp>=0 &&
          th<mintheta) {
          mintheta=th;
          newnbr=m;
        }
        /* sometimes next neighbor will be more than 180 degrees away 
           and so cross product will be negative 
           in this case, pick the neighbor farthest away
        */
        if(dotp<=0 &&
          th>maxtheta) {
          maxtheta=th;
          newnbr2=m;
        }
      }
      if(newnbr<0 && newnbr2<0) {
        MsgPrintf("%s: unable to connect all neighbors for vertex %d\n",progname,k);
        nbrs[p]=-1;
        nbrnums[p]=-1;
      } else if(newnbr<0) {
        nbrs[p]=oldnbrs[newnbr2];
        nbrnums[p]=newnbr2;
        oldnbrs[newnbr2]=-1;
      } else {
        nbrs[p]=oldnbrs[newnbr];
        nbrnums[p]=newnbr;
        oldnbrs[newnbr]=-1;
      }
    }

    /* get x,y,z coordinates for each neighbor, offset by center vertex */
    for (m=0;m<num_nbrs;m++) {
      j = nbrs[m];
      x[m] = mris->vertices[j].x - v->x;
      y[m] = mris->vertices[j].y - v->y;
      z[m] = mris->vertices[j].z - v->z;
    }

    /* calculate angles between each neighbor */
    for (m=0;m<num_nbrs-1;m++) {
      theta[m] = vectangle(x[m],y[m],z[m],x[m+1],y[m+1],z[m+1]);
    }
    /* angle between last neighbor and first neighbor */
    theta[m] = vectangle(x[m],y[m],z[m],x[0],y[0],z[0]);

    /* scale angles so total is 2*PI */
    sumtheta = 0;
    for (m=0;m<num_nbrs;m++) sumtheta+=theta[m];
    for (m=0;m<num_nbrs;m++) theta[m]*=(2.0*M_PI/sumtheta);

    /* calculate new locally flat coordinates */
    sumtheta = 0;
    for (m=0;m<num_nbrs;m++) {
      dist = distance(0,0,0,x[m],y[m],z[m]);
      fx[m] = dist*cos(sumtheta);
      fy[m] = dist*sin(sumtheta);
      sumtheta+=theta[m];     /* counterclockwise */
    }

    dv1dx = dv1dy = dv2dx = dv2dy = 0;
    m11 = m12 = m13 = m22 = m23 = z1 = z2 = z3 = z1b = z2b = z3b = 0;
    for (m=0;m<num_nbrs;m++) {
      j = nbrs[m];
      if(j<0) continue;
      dv1 = circsubtract(v->val,mris->vertices[j].val);
      dv2 = circsubtract(v->valbak,mris->vertices[j].valbak);
      dx = fx[m];
      dy = fy[m];
      m11 += dx*dx;
      m12 += dx*dy;
      m13 += dx;
      m22 += dy*dy;
      m23 += dy;
      z1 += dx*dv1;
      z2 += dy*dv1;
      z3 += dv1;
      z1b += dx*dv2;
      z2b += dy*dv2;
      z3b += dv2;
    }
    dv1dx = (m22*z1-m23*m23*z1-m12*z2+m13*m23*z2-m13*m22*z3+m12*m23*z3);
    dv2dx = (m22*z1b-m23*m23*z1b-m12*z2b+m13*m23*z2b-m13*m22*z3b+m12*m23*z3b);
    dv1dy = (-m12*z1+m13*m23*z1+m11*z2-m13*m13*z2+m12*m13*z3-m11*m23*z3);
    dv2dy = (-m12*z1b+m13*m23*z1b+m11*z2b-m13*m13*z2b+m12*m13*z3b-m11*m23*z3b);
    denom = -m12*m12+m11*m22-m13*m13*m22+2*m12*m13*m23-m11*m23*m23;
    if (denom!=0) {
      v->fieldsign = (dv1dx*dv2dy-dv2dx*dv1dy)/(denom*denom);
    } else {
      v->fieldsign = 0;
    }
    v->fieldsign =  ((v->fieldsign<0)?-1:(v->fieldsign>0)?1:0);
    if (revfsflag)
      v->fieldsign = -(v->fieldsign);
    v->fsmask = sqrt(v->imag_val*v->val2bak);  /* geom mean of r,th power */

    free(x);
    free(y);
    free(z);
    free(theta);
    free(fx);
    free(fy);
  }
}

void writeFieldSign(MRIS *mris) {
  int k,nverts;
  float f;
  FILE *fp_fieldsign, *fp_fsmask;
  char tempstr[STRLEN];

  sprintf(tempstr,"%s/%s-%s.fs",outdir,outstem,hemi);
  MsgPrintf("fieldsign file: %s\n",tempstr);
  fp_fieldsign = fopen(tempstr,"wb");
  if (fp_fieldsign==NULL)
    ErrorExit(ERROR_NOFILE,"%s: unable to create file %s\n",progname,tempstr) ;
  sprintf(tempstr,"%s/%s-%s.fm",outdir,outstem,hemi);
  MsgPrintf("fieldsign mask file: %s\n",tempstr);
  fp_fsmask = fopen(tempstr,"wb");
  if (fp_fsmask==NULL)
    ErrorExit(ERROR_NOFILE,"%s: unable to create file %s\n",progname,tempstr) ;

  nverts = mris->nvertices;
  MsgPrintf("%s: writing fieldsign for %d vertices\n",progname,nverts);
  fwriteInt(nverts,fp_fieldsign);
  fwriteInt(nverts,fp_fsmask);
  for (k=0;k<nverts;k++) {
    f = mris->vertices[k].fieldsign;
    fwriteFloat(f,fp_fieldsign);
    f = mris->vertices[k].fsmask;
    fwriteFloat(f,fp_fsmask);
  }
  fclose(fp_fieldsign);
  fclose(fp_fsmask);
}

int main(int argc, char **argv)
{
  int ecode=NO_ERROR;
  char tempstr[STRLEN];
  MRIS *mris;
  
  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,surf);
  MsgPrintf("%s: finished opening surface\n",progname);

  /* read eccentricity */
  MsgPrintf("%s: reading eccentricity files\n",progname);
  sprintf(tempstr,"%s/%s%s-%s.w",eccdir,eccstem,real_infix,hemi);
  ecode = MRISreadValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s%s-%s.w",eccdir,eccstem,imag_infix,hemi);
  ecode = MRISreadImagValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  MsgPrintf("%s: finished reading eccentricity files\n",progname);

  MRISsmoothComplexValues(mris,smooth);
  swapValues(mris);  /* push eccen values into bak vals */

  /* read polar angle */
  MsgPrintf("%s: reading polar angle files\n",progname);
  sprintf(tempstr,"%s/%s%s-%s.w",poldir,polstem,real_infix,hemi);
  ecode = MRISreadValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s%s-%s.w",poldir,polstem,imag_infix,hemi);
  ecode = MRISreadImagValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  MsgPrintf("%s: finished reading polar angle files\n",progname);

  MRISsmoothComplexValues(mris,smooth);

  MsgPrintf("%s: converting complex to polar\n",progname);
  complex2polar(mris);

  MsgPrintf("%s: computing fieldsign\n",progname);
  calcFieldSign(mris);

  MsgPrintf("%s: writing fieldsign\n",progname);
  writeFieldSign(mris);

  MsgPrintf("%s: finished\n",progname);
  exit(0);
}
