/* gradpath.cpp: calculate path of phase gradient for a surface roi
      created: 12/29/04 DH
     last mod: 02/23/06 DH

   purpose:
     calculating gradient of surface stats as measure of "mapness"

   input:
     complex stats file (w file)

   output:
     complex stats file (w file)
*/

#include <algorithm>
#include "surflib.h"
#include "clustlib.h"
using namespace std;

#define MINARGC 3
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"
#define MAX_ITERS  1000


/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char maskstem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char maskdir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char *real_infix=NULL;
char *imag_infix=NULL;
int avgdvflag = 0;
int nsteps = 1; // neighbor steps
int centmassflag = 0;
int allseedsflag = 0;
int toppaths = 1;
int fitlineflag = 0;
int ampdistflag = 0;


// class definitions
class Path {
public:
  Cluster clust;
  int seed;
  int head;
  int tail;
  float angle;
  float length;

  Path() {};
  // copy constructor
  Path(const Path &path) {
    clust = path.clust;
    seed = path.seed;
    head = path.head;
    tail = path.tail;
    angle = path.angle;
    length = path.length;
  }
  ~Path() {};
  int &operator[](int i) {
    if (i<0) {
      cout << "### path vertex index = -1 (undefined) ... quitting ###\n";
      exit(1);
    }    
    if (i>=size()) {
      cout << "### path vertices index overrun ... quitting ###\n";
      exit(1);
    }    
    return clust[i];
  }
  bool operator<(const Path &path) const {
    return size() < path.size();
  }
  Path operator=(Path path) {
    clust = path.clust;
    seed = path.seed;
    head = path.head;
    tail = path.tail;
    angle = path.angle;
    length = path.length;
    return *this;
  }
  int size() const {return clust.size();}
};

bool pathLengthComp(const Path &A, const Path &B) {
  return A.length < B.length;
}

class PathList {
private:
  vector<Path> paths; // dynamic array of paths
public:
  PathList() {};
  ~PathList() {};
  Path &operator[](int i) {
    if (i<0) {
      cout << "### pathlist index = -1 (undefined) ... quitting ###\n";
      exit(1);
    }
    if (i>=size()) {
      cout << "### pathlist index overrun ... quitting ###\n";
      exit(1);
    }
    return paths[i];
  }
  void newPath(Path &path) {
    paths.push_back(path);
  }
  int size() const {return paths.size();}
  void sortPaths() {
    stable_sort(paths.begin(),paths.end());
  }
  void sortPathsByLength() {
    stable_sort(paths.begin(),paths.end(),pathLengthComp);
  }
};



/* functions */

void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -maskstem maskstem -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem   instem           omit extension, infixes, hemi: <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -maskstem maskstem         stem for w file defining ROI\n");
  printf("    -subj     subjname         subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem   [instem-grad]   output file stem\n");
  printf("    -hemi      [rh]            hemisphere (rh or lh)\n");
  printf("    -indir     [.]             input dir\n");
  printf("    -maskdir   [.]             dir containing mask file\n");
  printf("    -outdir    [.]             output dir\n");
  printf("    -nsteps    [1]             number of neighbors away to include in fits\n");
  printf("    -avgdv                     calculate average derivative (not gradient)\n");
  printf("    -centmass                  seed path at center of mass (rather than center)\n");
  printf("    -allseeds                  seed a path at every vertex in roi\n");
  printf("    -toppaths  [1]             if allseeds, write out this many of the longest paths\n");
  printf("    -fitline                   fit line to path (rather than gradient)\n");
  printf("    -ampdist                   save distance as amplitude (rather than # of vertices)\n");
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
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-maskstem") && i+1<argc){
        strcpy(maskstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-maskdir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-nsteps") && i+1<argc) {
        nsteps = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-toppaths") && i+1<argc) {
        toppaths = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-avgdv")){
        avgdvflag = 1;
      } else
      if (MATCH(argv[i],"-centmass")){
        centmassflag = 1;
      } else
      if (MATCH(argv[i],"-allseeds")){
        allseedsflag = 1;
      } else
      if (MATCH(argv[i],"-fitline")){
        fitlineflag = 1;
      } else
      if (MATCH(argv[i],"-ampdist")){
        ampdistflag = 1;
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
  if (MATCH(maskstem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### maskstem not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    sprintf(outstem,"%s-gradpath",instem);
  }
  if (!FileExists(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",
              progname,indir);
  }
  if(!isadir(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
  }
  if (!FileExists(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### maskdir %s not found ...quitting\n",
              progname,maskdir);
  }
  if(!isadir(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,maskdir);
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
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
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

void swap_values(MRIS *mris)
{
  int k;
  float val,valbak;

  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].ripflag) {
      val = mris->vertices[k].val2;
      valbak = mris->vertices[k].imag_val;
      mris->vertices[k].val = mris->vertices[k].val2;
      mris->vertices[k].imag_val = mris->vertices[k].val2bak;
      mris->vertices[k].val2 = val;
      mris->vertices[k].val2bak = valbak;      
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

void fitline(vector<float> &x, vector<float> &y, float *intercept, float *slope)
// adapted from numerical recipes fit
{
  int i,n;
  float d,avgx,sumx=0.0,sumy=0.0,ssqd=0.0;

  if(x.size()!=y.size())
    ErrorExit(ERROR_BADPARM,"%s: ### fitline x and y vector sizes do not match... quitting\n",
              progname);
  n=x.size();
  
  *slope=0.0;
  for(i=0;i<n;i++){
    sumx+=x[i];
    sumy+=y[i];
  }
  avgx=sumx/n;
  for(i=0;i<n;i++){
    d=x[i]-avgx;
    ssqd+=d*d;
    *slope+=d*y[i];
  }
  if(ssqd!=0)
    *slope/=ssqd;
  *intercept=(sumy-sumx*(*slope))/n;
}

void defineClust(MRIS *mris, int k, Cluster *clust)
{
  int m,step,i,j,n;
  VERTEX *v;
  Cluster newclust;

  MRISclearMarks(mris);
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
  if (avgdvflag) {
    // calculate average derviative (should be zero if phase is random)
    *dvdx = *dvdy = 0;
    for(i=0;i<myclust.size();i++) {
      j = myclust[i];
      if (j==k) continue;
      dv = circsubtract(v->val,mris->vertices[j].val);
      dx = v->x - mris->vertices[j].x;
      dy = v->y - mris->vertices[j].y;
      if (dx!=0)
        *dvdx += dv/dx;
      if (dy!=0)
        *dvdy += dv/dy;
    }
    if(myclust.size()>1) {
      *dvdx /= (float)(myclust.size()-1);
      *dvdy /= (float)(myclust.size()-1);
    } else {
      *dvdx = 0;
      *dvdy = 0;
    }
  } else {
    // fit points to a plane (calculate partial derivatives of phase)
    *dvdx = *dvdy = 0;
    m11 = m12 = m13 = m22 = m23 = z1 = z2 = z3 = 0;
    for(i=0;i<myclust.size();i++) {
      j = myclust[i];
      if (j==k) continue;
      dv = circsubtract(v->val,mris->vertices[j].val);
      dx = v->x - mris->vertices[j].x;
      dy = v->y - mris->vertices[j].y;
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
  }
}


int getSeedVertex(MRIS *mris, Cluster *roi)
{
  int i,k,seed=-1;
  float x,y,cx=0,cy=0,amp,ampsum=0,dx,dy,dist,mindist=BIGFLOAT;
  
  // assumes that complex values have been converted to polar
  // with val = phase, and imag_val = amplitude

  // calculate center of mass coordinates
  MsgPrintf("%s: roi has %d vertices\n",progname,roi->size());

  for(i=0;i<roi->size();i++) {
    k=(*roi)[i];
    amp=mris->vertices[k].imag_val;
    x=mris->vertices[k].x;
    y=mris->vertices[k].y;
    if(centmassflag) {
      cx += amp*x;
      cy += amp*y;
      ampsum += amp;
    } else {
      cx += x;
      cy += y;
      ampsum++;
    }
  }
  cx /= ampsum;
  cy /= ampsum;
  if(centmassflag)
    MsgPrintf("%s: roi center of mass coordinates = (%0.3f,%0.3f)\n",
              progname,cx,cy);
  else
    MsgPrintf("%s: roi center coordinates = (%0.3f,%0.3f)\n",
              progname,cx,cy);
  
  // identify vertex closest to center of mass
  for(i=0;i<roi->size();i++) {
    k=(*roi)[i];
    x=mris->vertices[k].x;
    y=mris->vertices[k].y;
    dx=x-cx;
    dy=y-cy;
    dist=hypot(dx,dy);
    if(dist<mindist) {
      mindist=dist;
      seed=k;
    }    
  }

  MsgPrintf("%s: seed vertex = %d, distance from center = %f\n",
            progname,seed,mindist);
  
  return seed;
}

int getNextVertex(MRIS *mris, int k, float dvdx, float dvdy, int dir)
{
  float dx,dy,dv_ang,ang,diff,mindiff=BIGFLOAT;
  int j,next,m;

  next=k;

  if(dir==1)
    dv_ang = atan2(dvdy,dvdx);
  else
    dv_ang = atan2(-dvdy,-dvdx);
  
  for (m=0;m<mris->vertices[k].vnum;m++) {
    j=mris->vertices[k].v[m];
    dx = mris->vertices[k].x - mris->vertices[j].x;
    dy = mris->vertices[k].y - mris->vertices[j].y;
    ang = atan2(dy,dx);
    diff = fabs(circsubtract(dv_ang,ang));
    if(diff<mindiff) {
      mindiff=diff;
      next=j;
    }
  }
  return next;
}

void findGradPath(MRIS *mris, Cluster *roi, int seed, Path *path)
{
  int i,end;
  float dvdx,dvdy;
  Cluster clust;

  path->clust.clear();
  path->clust.addVertex(seed);
  path->seed = seed;
//  MsgPrintf("%s: findGradPath seed vertex = %d",progname,seed);

  // find plus end of path
//  MsgPrintf("%s: finding plus end of path\n",progname);
  end = seed;
  for(i=0;i<MAX_ITERS;i++) {
    path->head = end;  
    defineClust(mris,end,&clust);
    calcGrad(mris,end,&clust,&dvdx,&dvdy);
    end=getNextVertex(mris,end,dvdx,dvdy,1);
    if(path->clust.isMember(end)) {
      // if end is already on path, we are looping, so stop
//      MsgPrintf("%s: plus end of path loops back onto itself after %d steps\n",
//                progname,i);
      break;
    } else if (!roi->isMember(end)) {
      // if end is not in roi, we are finished
//     MsgPrintf("%s: plus end of path exits roi after %d steps\n",
//               progname,i);
      break;
    }
    path->clust.addVertex(end);
  }
  if(i>=MAX_ITERS)
    MsgPrintf("\n%s: plus end of path not found after %d steps\n",
              progname,i);
    
  // find minus end of path
//  MsgPrintf("%s: finding minus end of path\n",progname);
  end = seed;
  for(i=0;i<MAX_ITERS;i++) {
    path->tail = end;
    defineClust(mris,end,&clust);
    calcGrad(mris,end,&clust,&dvdx,&dvdy);
    end=getNextVertex(mris,end,dvdx,dvdy,-1);
    if(path->clust.isMember(end)) {
      // if end is already on path, we are looping, so stop
//      MsgPrintf("%s: minus end of path loops back onto itself after %d steps\n",
//                progname,i);
      break;
    } else if (!roi->isMember(end)) {
      // if end is not in roi, we are finished
//      MsgPrintf("%s: minus end of path exits roi after %d steps\n",
//                progname,i);
      break;
    }
    path->clust.addVertex(end);
  }
  if(i>=MAX_ITERS)
    MsgPrintf("\n%s: minus end of path not found after %d steps\n",
              progname,i);

//  MsgPrintf("  # vertices = %d\n",path->size());
}


int main(int argc, char **argv)
{
  int ecode=NO_ERROR;
  char tempstr[STRLEN];
  MRIS *mris;
  float *vals=NULL,*imag_vals=NULL;
  int i,j,k,nverts,seed,rank,p;
  Cluster roi;
  Path path;
  PathList pathlist;
  float dx,dy,dvdx,dvdy,amp,ang,slope,intercept;
  vector<float> x, y;

  parse_args(argc,argv);

  // load surface
  mris = openSurface(subj,hemi,"orig");
  nverts = mris->nvertices;
  MsgPrintf("%s: finished opening surface\n",progname);
  vals = new float[nverts]; MTEST(vals);

  // read input files
  MsgPrintf("%s: reading input files\n",progname);
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
  MsgPrintf("%s: reading mask file\n",progname);
  sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
  ecode = readSurfVals(tempstr,vals,nverts);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading mask file %s\n",
              progname,tempstr);
  for (k=0;k<nverts;k++)
    if(vals[k] && !mris->vertices[k].ripflag)
      roi.addVertex(k);
  delete [] vals;
  MsgPrintf("%s: %d vertices in roi\n",progname,roi.size());
  MsgPrintf("%s: finished reading input files\n",progname);
  readPatch(mris,subj,hemi);
  MsgPrintf("%s: finished reading patch\n",progname);
  MRISremoveTriangleLinks(mris);

  // calculate path of gradient
  MsgPrintf("%s: computing gradient path(s)...\n",progname);
  complex2polar(mris);
  if(!allseedsflag) {
    seed = getSeedVertex(mris,&roi);
    if (!roi.isMember(seed))
      ErrorExit(ERROR_BADPARM,"%s: ### error finding seed vertex... quitting\n",
                progname);
    findGradPath(mris,&roi,seed,&path);
    pathlist.newPath(path);
    toppaths = 1;
  } else {
    vals = new float[nverts]; MTEST(vals);
    imag_vals = new float[nverts]; MTEST(imag_vals);
    for(i=0;i<roi.size();i++) {
      seed=roi[i];
      findGradPath(mris,&roi,seed,&path);

      // calculate length of path from head to tail
      dy = mris->vertices[path.head].y - mris->vertices[path.tail].y;
      dx = mris->vertices[path.head].x - mris->vertices[path.tail].x;
      
      // todo: calculate max extent instead (max distance between any two points)
/*
      printf("### path %d: seed(%0.3f,%0.3f) head(%0.3f,%0.3f) tail(%0.3f,%0.3f)\n",
             i,
             mris->vertices[path.seed].x,
             mris->vertices[path.seed].y,
             mris->vertices[path.head].x,
             mris->vertices[path.head].y,
             mris->vertices[path.tail].x,
             mris->vertices[path.tail].y);
*/      

      // calculate angle of path and save as complex
      if(fitlineflag) {
        // calculate direction by linear fit
        x.clear();
        y.clear();
        for(j=0;j<path.size();j++) {
          k=path[j];
          x.push_back(mris->vertices[k].x);
          y.push_back(mris->vertices[k].y);
        }
        fitline(x,y,&intercept,&slope);
        dvdx=1;
        dvdy=slope;

//        printf("### path %d: fitline slope = %0.3f, intercept = %0.3f\n",
//               i,slope,intercept);
        
        if((dx==0 && dy<0) ||
           (dy==0 && dx<0) ||
           (dy*dx < 0)) {
            dvdx=-dvdx;
            dvdy=-dvdy;
        }
        // todo: find a surer way to determine which way slope should go
        //       maybe compare seed to nearest neighbors in path?
      } else
        calcGrad(mris,seed,&(path.clust),&dvdx,&dvdy);
      path.length = hypot(dx,dy);
      path.angle = atan2(dvdy,dvdx);
      ang = path.angle;
      if(ampdistflag)
        amp = path.length;
      else
        amp = (float)path.size();
      vals[seed] = amp*cos(ang);
      imag_vals[seed] = amp*sin(ang);
      pathlist.newPath(path);
    }
    sprintf(tempstr,"%s/%s-pathlengths%s-%s.w",outdir,outstem,real_infix,hemi);
    MsgPrintf("%s: writing file %s\n",progname,tempstr);
    ecode = writeSurfVals(tempstr,vals,nverts);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error writing output file %s\n",
                progname,tempstr);
    sprintf(tempstr,"%s/%s-pathlengths%s-%s.w",outdir,outstem,imag_infix,hemi);
    MsgPrintf("%s: writing file %s\n",progname,tempstr);
    ecode = writeSurfVals(tempstr,imag_vals,nverts);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error writing output file %s\n",
                progname,tempstr);
    delete [] vals;
    delete [] imag_vals;

    if(ampdistflag)
      pathlist.sortPathsByLength();
    else
      pathlist.sortPaths();
  }
  polar2complex(mris);

  if(toppaths > pathlist.size()) toppaths = pathlist.size();

  for(p=pathlist.size()-1,rank=1;p>=pathlist.size()-toppaths;p--,rank++) {
    path = pathlist[p];
    // output info to screen
    MsgPrintf("%s: path %d has %d of %d vertices\n",
              progname,rank,path.size(),roi.size());
    MsgPrintf("\tvertex\tx\ty\tamp\tphase\n");
    for(i=0;i<path.size();i++) {
      k=path[i];
      MsgPrintf("\t%d\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n",
                k,
                mris->vertices[k].x,
                mris->vertices[k].y,
                mris->vertices[k].imag_val,
                mris->vertices[k].val);
    }
    MsgPrintf("%s: seed=%d, tail=%d, head=%d, angle=%0.3f, length=%0.3f\n",
              progname,path.seed,path.tail,path.head,path.angle,path.length);

    // output to files
    vals = new float[nverts]; MTEST(vals);
    imag_vals = new float[nverts]; MTEST(imag_vals);

    for(i=0;i<path.size();i++) {
      k=path[i];
//      vals[k] = mris->vertices[k].val;
//      imag_vals[k] = mris->vertices[k].imag_val;
      if(k==path.seed) {
        ang = M_PI;
        amp = 100;
      } else if (k==path.head) {
        ang = M_PI/2;
        amp = 100;
      } else if (k==path.tail) {
        ang = M_PI*3/2;
        amp = 100;
      } else {
        ang = path.angle;
        amp = path.length;
      }
      vals[k] = amp*cos(ang);
      imag_vals[k] = amp*sin(ang);
    }

    sprintf(tempstr,"%s/%s-path%03d%s-%s.w",
            outdir,outstem,rank,real_infix,hemi);
    MsgPrintf("%s: writing file %s\n",progname,tempstr);
    ecode = writeSurfVals(tempstr,vals,nverts);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error writing output file %s\n",
                progname,tempstr);
    sprintf(tempstr,"%s/%s-path%03d%s-%s.w",
            outdir,outstem,rank,imag_infix,hemi);
    MsgPrintf("%s: writing file %s\n",progname,tempstr);
    ecode = writeSurfVals(tempstr,imag_vals,nverts);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error writing output file %s\n",
                progname,tempstr);
    delete [] vals;
    delete [] imag_vals;
  }

  MsgPrintf("%s: finished\n",progname);
  exit(0);
}

