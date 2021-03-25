/* surfgrad.cpp: read surface stats file and calculate gradient for each vertex
      created: 10/28/04 DH
     last mod: 12/30/04 DH

   purpose:
     calculating gradient of surface stats as measure of "mapness"

   input:
     complex stats file (w file)

   output:
     complex stats file (w file)
*/

#include "surflib.h"
#include "clustlib.h"
using namespace std;

#define MINARGC 3
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"

#define TYPE_GRAD    1
#define TYPE_AVGDV   2
#define TYPE_WAVGDV  3
#define TYPE_HEEGER  4

#define STR_GRAD "grad"
#define STR_AVGDV "avgdv"
#define STR_WAVGDV "wavgdv"
#define STR_HEEGER "heeger"


/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char *real_infix=NULL;
char *imag_infix=NULL;
int outcorrflag = 0;
int calctype=TYPE_GRAD;
char typestr[STRLEN]=STR_GRAD;
int noisenormflag = 0;
int amplbiasflag = 0;
int nsteps = 1; // neighbor steps
int smooth = 50; /* pre-smoothing steps */

/* functions */

void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -name subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem    instem           omit extension, infixes, hemi:\n");
  printf("                                  <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -name      subjname         subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -type      [grad]          grad, avgdv, wavgdv, or heeger\n");
  printf("                                  grad:   best fit gradient\n");
  printf("                                  avgdv:  average derivative\n");
  printf("                                  wavgdv: weighted average derivative\n");
  printf("                                  heeger: heeger's method (left-right/up-down)\n");
  printf("    -outstem   [instem-grad]   output file stem\n");
  printf("    -hemi      [rh]            hemisphere (rh or lh)\n");
  printf("    -indir     [.]             input dir\n");
  printf("    -outcorr                   output the correlation coefficient\n");
  printf("    -smooth    [50]            number of pre-smoothing steps\n");
  printf("    -nsteps    [1]             number of neighbors away to include in fits\n");
  printf("    -noisenorm                 normalize gradient by average abs(dv/dxy)\n");
  printf("    -amplbias                  bias gradient amplitudes with raw amplitudes\n");
  printf("    -infixes   [_r _i]         real,imaginary infixes\n");
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
      if (MATCH(argv[i],"-type") && i+1<argc) {
        strcpy(typestr,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
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
      if (MATCH(argv[i],"-smooth") && i+1<argc) {
        smooth = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-nsteps") && i+1<argc) {
        nsteps = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outcorr")){
        outcorrflag = 1;
      } else
      if (MATCH(argv[i],"-noisenorm")){
        noisenormflag = 1;
      } else
      if (MATCH(argv[i],"-amplbias")){
        amplbiasflag = 1;
      } else
      if ((MATCH(argv[i],"-infixes")) && i+2<argc){
        strcpy(real_infix,argv[i+1]); strcpy(imag_infix,argv[i+2]); i+=2;
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
  if (MATCH(instem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### instem not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    sprintf(outstem,"%s-grad",instem);
  }
  if (!FileExists(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",
              progname,indir);
  }
  if(!isadir(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
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
  if(!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  
  if(MATCH(typestr,STR_GRAD)) calctype=TYPE_GRAD; else
  if(MATCH(typestr,STR_AVGDV)) calctype=TYPE_AVGDV; else
  if(MATCH(typestr,STR_WAVGDV)) calctype=TYPE_WAVGDV; else
  if(MATCH(typestr,STR_HEEGER)) calctype=TYPE_HEEGER; else
    ErrorExit(ERROR_BADPARM,"%s: ### bad calc type %s\n",progname,typestr);
  
  MsgPrintf("%s: finished parsing arguments\n",progname);
}


void readPatch(MRIS *mris, char *subject, char *hemi)
{
  char fname[STRLEN];
  char *SUBJECTS_DIR = NULL;

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL)
    ErrorExit(ERROR_BADPARM,"%s: ### environment variable SUBJECTS_DIR undefined (use setenv)\n",
      progname);
  sprintf(fname,"%s/%s/surf/%s.full.patch.flat",SUBJECTS_DIR,subject,hemi);
  if(!FileExists(fname))
    ErrorExit(ERROR_BADFILE,"%s: ### file %s not found... quitting\n",
      progname,fname);
  MsgPrintf("%s: reading patch %s\n",progname,fname);
  if (mris->patch) MRISunrip(mris);
  MRISreadPatchNoRemove(mris, fname) ;
}

void complex2polar(MRIS *mris)
{
  int k;
  float phase,amplitude;
  
  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].ripflag) {
      phase = atan2(mris->vertices[k].imag_val,mris->vertices[k].val);
      amplitude = hypot(mris->vertices[k].val,mris->vertices[k].imag_val);
      mris->vertices[k].val = phase;
      mris->vertices[k].imag_val = amplitude;
    }
  }
}

void polar2complex(MRIS *mris)
{
  int k;
  float phase,amplitude;
  
  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].ripflag) {
      phase = mris->vertices[k].val;
      amplitude = mris->vertices[k].imag_val;
      mris->vertices[k].val = amplitude*cos(phase);
      mris->vertices[k].imag_val = amplitude*sin(phase);
    }
  }
}

void swapValues(MRIS *mris)
{
  VERTEX *v;
  int k;
  float temp;

  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].ripflag) {
      v = &mris->vertices[k];

      temp = v->val;
      v->val = v->valbak;
      v->valbak = temp;

      temp = v->imag_val;
      v->imag_val = v->val2bak;
      v->val2bak = temp;      
    }
  }
}


float circsubtract(float a,float b)
{
  float h = a-b;
  if (h<-M_PI) h = h+2*M_PI;
  else if (h>M_PI) h = h-2*M_PI;
  return h;
}

void defineClust(MRIS *mris, int k, Cluster *clust)
{
  int m,step,i,j,n;
  VERTEX *v;
  Cluster newclust;

//  MRISclearMarks(mris);
  v = &mris->vertices[k];

  // define list of vertices to fit
  clust->clear();
  clust->addVertex(k);
  mris->vertices[k].marked=-1; // indicates that vertex has been added
  for (step=0;step<nsteps;step++) {
    newclust.clear();
    for (i=0;i<clust->size();i++) {
      j = (*clust)[i];
      if (mris->vertices[j].marked!=1) { // skip if we already did this vertex
        for (m=0;m<mris->vertices[j].vnum;m++) {
          n = mris->vertices[j].v[m];
          if (!mris->vertices[n].ripflag && !mris->vertices[n].marked) {
            newclust.addVertex(n);
            mris->vertices[n].marked=-1;
          }
        }
      }
      mris->vertices[j].marked=1; // indicates that neighbors have been added
    }
    clust->addCluster(newclust);
  }
}

void calcGrad(MRIS *mris, int k, const Cluster *clust, float *dvdx, float *dvdy)
{
  int i,j;
  VERTEX *v;
  float dv,dx,dy;
  float m11,m12,m13,m22,m23,z1,z2,z3,denom;
  
  Cluster myclust = *clust;

  v = &mris->vertices[k];
  switch(calctype) {
    case TYPE_GRAD :
      // fit points to a plane (calculate partial derivatives of phase)
      *dvdx = *dvdy = 0;
      m11 = m12 = m13 = m22 = m23 = z1 = z2 = z3 = 0;
      for(i=0;i<myclust.size();i++) {
        j = myclust[i];
        if (j==k) continue;
        dx = mris->vertices[j].x - v->x;
        dy = mris->vertices[j].y - v->y;
        dv = circsubtract(mris->vertices[j].val,v->val);
        m11 += dx*dx;
        m12 += dx*dy;
        m13 += dx;
        m22 += dy*dy;
        m23 += dy;
        z1 += dx*dv;
        z2 += dy*dv;
        z3 += dv;
      }
      *dvdx = (m22*z1-m23*m23*z1-m12*z2+m13*m23*z2-m13*m22*z3+m12*m23*z3);
      *dvdy = (-m12*z1+m13*m23*z1+m11*z2-m13*m13*z2+m12*m13*z3-m11*m23*z3);
      denom = -m12*m12+m11*m22-m13*m13*m22+2*m12*m13*m23-m11*m23*m23;
      if (denom!=0) {
        *dvdx /= denom;
        *dvdy /= denom;
      } else {
        *dvdx = 0;
        *dvdy = 0;
      }
      break;
    case TYPE_AVGDV :
      // calculate average derviative
      *dvdx = *dvdy = 0;
      for(i=0;i<myclust.size();i++) {
        j = myclust[i];
        if (j==k) continue;
        dx = mris->vertices[j].x - v->x;
        dy = mris->vertices[j].y - v->y;
        dv = circsubtract(mris->vertices[j].val,v->val);
        if (dx!=0)
          *dvdx += dv/dx;
        if (dy!=0)
          *dvdy += dv/dy;
      }
      break;
    case TYPE_WAVGDV :
      // calculate weighted average derviative
      *dvdx = *dvdy = 0;
      for(i=0;i<myclust.size();i++) {
        j = myclust[i];
        if (j==k) continue;
        dx = mris->vertices[j].x - v->x;
        dy = mris->vertices[j].y - v->y;
        dv = circsubtract(mris->vertices[j].val,v->val) *
             sqrt(v->imag_val*mris->vertices[j].imag_val/2) / M_PI;
        // dv multiplied by geom mean of significance amplitudes
        // dv multiplied by 1/sqrt(2) so that amp of resulting vector = 1
        // dv divided by PI to normalize max difference to 1
        if (dx!=0)
          *dvdx += dv/dx;
        if (dy!=0)
          *dvdy += dv/dy;
      }
      break;
    case TYPE_HEEGER :
      // calculate average derviative
      *dvdx = *dvdy = 0;
      for(i=0;i<myclust.size();i++) {
        j = myclust[i];
        if (j==k) continue;
        dx = mris->vertices[j].x - v->x;
        dy = mris->vertices[j].y - v->y;
        dv = sin(mris->vertices[j].val) - sin(v->val);
        if      (dx > 0) *dvdx += dv;
        else if (dx < 0) *dvdx -= dv;
        if      (dy > 0) *dvdy += dv;
        else if (dy < 0) *dvdy -= dv;
      }
      break;
  }
  // reset marks for cluster vertices only
  for(i=0;i<myclust.size();i++) {
    j = myclust[i];
    mris->vertices[j].marked = 0;
  }
}

void calcCorr(MRIS *mris, int k, const Cluster *clust, float dvdx, float dvdy,
              float *corrvals)
{
  int i,j;
  VERTEX *v;
  float dv,dx,dy;
  double err, sumsqerr, sumsq;

  Cluster myclust = *clust;
  v = &mris->vertices[k];

  // calculate correlation error
  sumsqerr=sumsq=0;
  for(i=0;i<myclust.size();i++) {
    j=myclust[i];
    if (j==k) continue;
    dv = circsubtract(mris->vertices[j].val,v->val);
    dx = mris->vertices[j].x - v->x;
    dy = mris->vertices[j].y - v->y;
    err = dv - (dvdx*dx + dvdy*dy);
    sumsqerr += err*err;
    sumsq += dv*dv;
  }
  if (sumsq!=0) {
    corrvals[k] = 1 - sumsqerr/sumsq;
    if (corrvals[k]<0) corrvals[k]=0;
  } else {
    corrvals[k] = 0;
  }
}

void noiseNorm(MRIS *mris, int k, const Cluster *clust, float *dvdx, float *dvdy)
{
  int i,j;
  VERTEX *v;
  float dv,dx,dy,avgdvdx,avgdvdy;

  Cluster myclust = *clust;
  v = &mris->vertices[k];

  // calculate correlation error
  for(i=0;i<myclust.size();i++) {
    j=myclust[i];
    if (j==k) continue;
    dv = circsubtract(mris->vertices[j].val,v->val);
    dx = mris->vertices[j].x - v->x;
    dy = mris->vertices[j].y - v->y;
    if (dx!=0)
      avgdvdx += fabs(dv/dx);
    if (dy!=0)
      avgdvdy += fabs(dv/dy);
  }
  if(myclust.size()>1) {
    avgdvdx /= (float)(myclust.size()-1);
    avgdvdy /= (float)(myclust.size()-1);
  }
  if (avgdvdx==0)
    *dvdx = 0;
  else
    *dvdx /= avgdvdx;
  if (avgdvdy==0)
    *dvdy = 0;
  else
    *dvdy /= avgdvdy;
}

void compute_gradient(MRIS *mris, float *corrvals)
{
  int k,nverts,avgclustsize=0;
  VERTEX *v;
  float dvdx,dvdy;
  double avggradamp=0,amp,ang;
  Cluster clust;

  MRISclearMarks(mris);
  nverts = mris->nvertices;
  for (k=0;k<nverts;k++) {
    if (!mris->vertices[k].ripflag) {
      v = &mris->vertices[k];
      defineClust(mris,k,&clust);
      avgclustsize += clust.size();
      calcGrad(mris,k,&clust,&dvdx,&dvdy);
      if (outcorrflag) calcCorr(mris,k,&clust,dvdx,dvdy,corrvals);
      if (noisenormflag) noiseNorm(mris,k,&clust,&dvdx,&dvdy);

      amp = hypot(dvdx,dvdy);
      if(amplbiasflag) {
        ang = atan2(dvdy,dvdx);
        amp = sqrt(amp*v->imag_val); // geometric mean of grad amp and raw amp
        dvdx = amp*cos(ang);
        dvdy = amp*sin(ang);
      }
      avggradamp += amp;
      v->valbak = dvdx;
      v->val2bak = dvdy;
    }
  }
  avggradamp /= nverts;
  avgclustsize /= nverts;
  MsgPrintf("%s: average cluster size = %d vertices\n",progname,avgclustsize);
  MsgPrintf("%s: average gradient amplitude = %0.4f\n",progname,avggradamp);
}

int main(int argc, char **argv)
{
  int ecode=NO_ERROR;
  char tempstr[STRLEN];
  MRIS *mris;
  float *corrvals=NULL;
  
  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,"orig");
  MsgPrintf("%s: finished opening surface\n",progname);
  readPatch(mris,subj,hemi);
  MRISremoveTriangleLinks(mris);
  MsgPrintf("%s: finished reading patch\n",progname);
   
  corrvals = (float *)calloc(mris->nvertices,sizeof(float)); MTEST(corrvals);

  /* read input files */
  MsgPrintf("%s: reading values files\n",progname);
  sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,real_infix,hemi);
  ecode = MRISreadValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,imag_infix,hemi);
  ecode = MRISreadImagValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  MsgPrintf("%s: finished reading values files\n",progname);
  MRISsmoothComplexValues(mris,smooth);

  MsgPrintf("%s: converting complex to polar\n",progname);
  complex2polar(mris);
  MsgPrintf("%s: computing gradient\n",progname);
  compute_gradient(mris,corrvals);
  swapValues(mris);
  
  /* write output to file */
  sprintf(tempstr,"%s/%s%s-%s.w",outdir,outstem,real_infix,hemi);
  ecode = MRISwriteValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output files %s\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s%s-%s.w",outdir,outstem,imag_infix,hemi);
  ecode = MRISwriteImagValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output files %s\n",
              progname,tempstr);

  if(outcorrflag) {
    /* write corrvals to file */
    sprintf(tempstr,"%s/%s-corr-%s.w",outdir,outstem,hemi);
    if(writeSurfVals(tempstr,corrvals,mris->nvertices)==-1) {
      ErrorExit(ERROR_BADFILE,"%s: ### error writing surface values...quitting\n",
           progname);
    }
  }
  

  MsgPrintf("%s: finished\n",progname);
  exit(0);
}

