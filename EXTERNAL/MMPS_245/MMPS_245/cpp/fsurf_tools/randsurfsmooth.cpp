/* randsurfsmooth.cpp: generate value at random vertex, smooth, calculate fwhm
      created: 02/11/05 DH
     last mod: 04/16/13 DH

   purpose:
     calculating fwhm for different smoothing steps

   input:
     options only

   output:
     stdout only

   acknowledgements:

*/

#include "surflib.h"
#include "clustlib.h"
using namespace std;

#define MINARGC 3
#define INITVALUE 1000
#define RATIO(X,Y)        ((float)X/(float)Y)
#ifdef Darwin
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),INT_MAX)*((H)-(L)))
#else
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),MAXINT)*((H)-(L)))
#endif
#define RAND_INT(H)  (int)RAND_FLOAT(0,H)
#define SMALLFLOAT 1e-10

// todo: optional pdf (probability distribution function) output

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char subj[STRLEN]=UNDEFSTR;
char instem[STRLEN]=UNDEFSTR;
char maskstem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char maskdir[STRLEN]=".";
char outdir[STRLEN]=".";
char hemi[STRLEN]="rh";
char real_infix[STRLEN]="_r";
char imag_infix[STRLEN]="_i";
char surf[STRLEN]="white";
int niter = 1;
int minsmooth = 0;
int maxsmooth = 100;
int smoothinc = 10;
int growclusterflag = 0;
int heatflag = 0;
int complexflag = 0;
int dataflag = 0;
int maskflag = 0;
int outflag = 0;
float stdev = 1.0;
float heatsigma = 0;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -subj     subjname    subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -instem   instem      omit extension, infixes, hemi:\n");
  printf("                            <instem>_{r,i}-{rh,lh}.w\n");
  printf("                            If instem not specified, will generate\n");
  printf("                            Gaussian noise instead\n");
  printf("    -indir     [.]        input dir\n");
  printf("    -maskstem maskstem    stem for w file defining ROI\n");
  printf("    -maskdir   [.]        dir containing mask file\n");
  printf("    -outstem [randsmooth] text output file stem\n");
  printf("    -outdir    [.]        output dir\n");
  printf("    -growcluster          grow cluster from non-zero seed\n");
  printf("    -hemi      [rh]       hemisphere (rh or lh)\n");
  printf("    -niter     [1]       number of iterations\n");
  printf("    -minsmooth [0]        minimum number of smoothing steps\n");
  printf("    -maxsmooth [100]      maximum number of smoothing steps\n");
  printf("    -smoothinc [10]       incremental increase in smooth steps\n");
  printf("    -surf      [white] surface used for area calculations\n");
  printf("    -stdev     [1.0]      standard deviation of random values\n");
  printf("    -heatsigma [0.0]      sigma of heat keernel smoothing\n");
  printf("    -avgdist              use average distance for noise fwhm\n");
  printf("    -complex              complex (phase-encoded) data\n");
  printf("    -quiet                suppress messages\n");
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
      if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
        dataflag=1;
      } else
      if (MATCH(argv[i],"-maskstem") && i+1<argc){
        strcpy(maskstem,argv[i+1]); i++;
        maskflag=1;
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
      if (MATCH(argv[i],"-maskdir") && i+1<argc) {
        strcpy(maskdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-niter") && i+1<argc){
        niter = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-minsmooth") && i+1<argc){
        minsmooth = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-maxsmooth") && i+1<argc){
        maxsmooth = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-smoothinc") && i+1<argc){
        smoothinc = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-growcluster")){
        growclusterflag = 1;
      } else
      if (MATCH(argv[i],"-heatsigma")){
        heatsigma = atof(argv[i+1]); i++;
        heatflag = 1;
      } else
      if (MATCH(argv[i],"-complex")){
        complexflag = 1;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-stdev") && i+1<argc) {
        stdev = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    } else
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(subj,UNDEFSTR))
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
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
  if (!MATCH(outstem,UNDEFSTR))
    outflag=1;
  if (!FileExists(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  }
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh"))
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",
              progname,hemi);
  if (maxsmooth < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### %s: maximum smooth steps (%d) must be >0\n",
              progname,maxsmooth);
  if (minsmooth > maxsmooth)
    minsmooth = maxsmooth;
  if (minsmooth < 0)
    minsmooth = 0;
  if (smoothinc <= 0)
    ErrorExit(ERROR_BADPARM,"%s: ### %s: smooth increment (%d) must be >0\n",
              progname,smoothinc);
  if (niter < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### %s: number of iterations (%d) must be >0\n",
              progname,niter);
  if (heatflag && heatsigma <= 0)
    ErrorExit(ERROR_BADPARM,"%s: ### %s: heat kernel sigma (%0.3f) must be >0\n",
              progname,heatsigma);

  MsgPrintf("%s: finished parsing arguments\n",progname);
}

float uniform ()
{
  return ( (float)drand48() );
}

void normal (float * n1, float * n2)
{
  float u1, u2;
  float r;

  u1 = 0.0;
  while (u1 <= 0.0) {u1 = uniform();}
  u2 = uniform();

#if 0
  r = sqrt(-2.0*log(u1));
#else
  r = stdev * sqrt(-2.0*log(u1));
#endif
  *n1 = r * cos(2.0*M_PI*u2);
  *n2 = r * sin(2.0*M_PI*u2);
}

void genRandSurf(MRIS *mris)
{
  int k, nverts, nvertsdiv2;
  float n1, n2;
  
  nverts = mris->nvertices;
  nvertsdiv2 = nverts/2;
  
  for (k=0;k<nvertsdiv2;k++) {
    normal(&n1,&n2);
    mris->vertices[k].val = n1;
    mris->vertices[k+nvertsdiv2].val = n2;
  }
  normal(&n1,&n2);
  mris->vertices[nverts-1].val = n1;
}

void genRandSurf(MRIS *mris, int *maskverts, int nmaskverts)
{
  int j,k, nvertsdiv2;
  float n1, n2;
  
  nvertsdiv2 = nmaskverts/2;
  for (j=0;j<nvertsdiv2;j++) {
    normal(&n1,&n2);
    k=maskverts[j];
    mris->vertices[k].val = n1;
    k=maskverts[j+nvertsdiv2];
    mris->vertices[k].val = n2;
  }  
  normal(&n1,&n2);
  k=maskverts[nmaskverts-1];
  mris->vertices[k].val = n1;
}

void genCxRandSurf(MRIS *mris)
{
  int k;
  float n1, n2;
  
  for (k=0;k<mris->nvertices;k++) {
    normal(&n1,&n2);
    mris->vertices[k].val = n1;
    mris->vertices[k].imag_val = n2;
  }  
}

void genCxRandSurf(MRIS *mris, int *maskverts, int nmaskverts)
{
  int j,k;
  float n1, n2;
  
  for (j=0;j<nmaskverts;j++) {
    normal(&n1,&n2);
    k=maskverts[j];
    mris->vertices[k].val = n1;
    mris->vertices[k].imag_val = n2;
  }  
}

int main(int argc, char **argv)
{
  int i,j,k,nverts,nmaskverts=0,seedvert=0,s;
  float *avgdiam, *stdevdiam;
  MRIS *mris;
  ClusterList clustlist;
  Cluster clust;
  int threshabsflag = 0;
  float halfmax=0,diam,area;
  float *maskvals=NULL;
  int *maskverts=NULL;
  float *vals=NULL, *imag_vals=NULL;
  FILE *fp=NULL;
  char tempstr[STRLEN];
  int ecode=NO_ERROR;

  parse_args(argc,argv);

  randseed();

  if(outflag) {
    // open output file
    sprintf(tempstr,"%s/%s.txt",outdir,outstem);
    fp = fopen(tempstr,"w");
    if(fp==NULL)
      ErrorExit(ERROR_NOFILE,"%s: ### can't create file %s\n",progname,tempstr);
  }

  // load surface
  MsgPrintf("%s: loading surface files\n",progname);
  mris = openSurface(subj,hemi,surf);
  nverts = mris->nvertices;
  MsgPrintf("%s: nverts = %d\n",progname,nverts);

  // load data
  if(dataflag && !growclusterflag) {
    if(complexflag) {
      // read complex input files
      MsgPrintf("%s: reading input files\n",progname);
      vals = new float[nverts]; MTEST(vals);
      sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,real_infix,hemi);
      ecode = readSurfVals(tempstr,vals,nverts);
      if(ecode)
        ErrorExit(ecode,"%s: ### error reading value file %s\n",
                  progname,tempstr);
      imag_vals = new float[nverts]; MTEST(imag_vals);
      sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,imag_infix,hemi);
      ecode = readSurfVals(tempstr,imag_vals,nverts);
      if(ecode!=NO_ERROR)
        ErrorExit(ecode,"%s: ### error reading value file %s\n",
                  progname,tempstr);
    } else {
      // read input file
      MsgPrintf("%s: reading input file\n",progname);
      vals = new float[nverts]; MTEST(vals);
      sprintf(tempstr,"%s/%s-%s.w",indir,instem,hemi);
      ecode = readSurfVals(tempstr,vals,nverts);
      if(ecode)
        ErrorExit(ecode,"%s: ### error reading value file %s\n",
                  progname,tempstr);
    }
  }
  if(maskflag) {
    MsgPrintf("%s: reading mask file\n",progname);
    maskvals = new float[nverts]; MTEST(maskvals);
    sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
    ecode = readSurfVals(tempstr,maskvals,nverts);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error reading mask file %s\n",
                progname,tempstr);
    MsgPrintf("%s: finished reading input files\n",progname);
    for (k=0;k<nverts;k++) {
      if(maskvals[k]>0) {
        nmaskverts++;
        mris->vertices[k].undefval=0;
      } else {
        mris->vertices[k].undefval=1;
      }
    }
    MsgPrintf("%s:   nmaskverts = %d\n",progname,nmaskverts);
    maskverts = new int[nmaskverts]; MTEST(maskverts);
    for (k=0,j=0;k<nverts;k++) {
      if(maskvals[k]>0) maskverts[j++]=k;
    }
  } else {
    nmaskverts = nverts;
    maskverts = new int[nverts]; MTEST(maskverts);
    for (k=0;k<nverts;k++) {
      maskverts[k]=k;
      mris->vertices[k].undefval=0;
    }  
  }

  // initialize arrays
  avgdiam = new float[maxsmooth+1];
  stdevdiam = new float[maxsmooth+1];
  for (s=minsmooth;s<=maxsmooth;s++) {
    avgdiam[s]=0;
    stdevdiam[s]=0;
  }

  for (i=0;i<niter;i++) {
    MsgPrintf("%s: ## iteration %d ##\n",progname,i);
    MRISclearValues(mris);
    if(growclusterflag) {
      MsgPrintf("%s: generating value at random vertex\n",progname);
      j = RAND_INT(nmaskverts);
      seedvert = maskverts[j];
      mris->vertices[seedvert].val = INITVALUE;
    } else {
      if(dataflag) {
        MsgPrintf("%s: restoring original data values at each vertex\n",progname);
        for (j=0;j<nmaskverts;j++) {
          k=maskverts[j];
          mris->vertices[k].val=vals[k];
          if(complexflag)
            mris->vertices[k].imag_val=imag_vals[k];
        }      
      } else {
        MsgPrintf("%s: generating random values at each vertex\n",progname);
        if(maskflag)
          if(complexflag)
            genCxRandSurf(mris,maskverts,nmaskverts);
          else
            genRandSurf(mris,maskverts,nmaskverts);
        else
          if(complexflag)
            genCxRandSurf(mris);
          else
            genRandSurf(mris);
      }
    }

    // ROI smoothing to skip vertices outside of mask
    if(minsmooth>0) {
      if(complexflag)
        MRISsmoothComplexValuesROI(mris,minsmooth);
        // todo: heat kernel smoothing for complex values -- why bother?
      else if(heatflag)
        MRISsmoothValuesHeatROI(mris,minsmooth,heatsigma);
      else
        MRISsmoothValuesROI(mris,minsmooth);
    }
    for (s=minsmooth;s<=maxsmooth;s+=smoothinc) {
      if(s>minsmooth) {
        if(complexflag)
          MRISsmoothComplexValuesROI(mris,smoothinc);
        else if(heatflag)
          MRISsmoothValuesHeatROI(mris,smoothinc,heatsigma);
        else
          MRISsmoothValuesROI(mris,smoothinc);
      }
      if(growclusterflag) {
        initSurfClustIndices(mris);
        halfmax = mris->vertices[seedvert].val/2;
        clust.clear();
        clust.addVertex(seedvert);
        clust.growCluster(seedvert,mris,halfmax,threshabsflag);
        clust.calcStats(mris);
        area = clust.Area();
        diam = 2*sqrt(area/M_PI);
        MsgPrintf("%s: smooth=%d, clustverts=%d, area=%0.4f, diam=%0.4f\n",
                  progname,s,clust.size(),clust.Area(),diam);
      } else {
        if(complexflag)
          diam=MRISgaussCxFWHM(mris,1);
        else
          diam=MRISgaussFWHM(mris,1);
        MsgPrintf("%s: smooth=%d, fwhm=%0.4f\n",
            progname,s,diam);
      }
      avgdiam[s]+=diam;
      stdevdiam[s]+=diam*diam;
    }
  }

  MsgPrintf("%s: iterations completed = %d\n",progname,niter);
  MsgPrintf("%s: min smooth steps = %d\n",progname,minsmooth);
  MsgPrintf("%s: max smooth steps = %d\n",progname,maxsmooth);
  if(outflag)
    MsgPrintf("%s: saving results to %s/%s.txt:\n",progname,outdir,outstem);

  if(niter>1) {
    if(outflag)
      fprintf(fp,"smooth steps ; avg diam (mm) ; stdev\n");
    printf("smooth steps ; avg diam (mm) ; stdev\n");
  } else {
    if(outflag)
      fprintf(fp,"smooth steps ; avg diam (mm)\n");
    printf("smooth steps ; avg diam (mm)\n");
  }
  for (s=minsmooth;s<=maxsmooth;s+=smoothinc) {
    if(niter<=1) {
      stdevdiam[s] = 0;
      if(outflag)
        fprintf(fp,"%d  ;  %0.10f\n",s,avgdiam[s]);
      printf("%d  ;  %0.10f\n",s,avgdiam[s]);
    } else {
      stdevdiam[s] = (stdevdiam[s] - (avgdiam[s]*avgdiam[s])/niter)/(niter-1);
      if(stdevdiam[s]<0)
        stdevdiam[s]=0;
      else
        stdevdiam[s]=sqrt(stdevdiam[s]);
      avgdiam[s]/=niter;

      if(outflag)
        fprintf(fp,"%d  ;  %0.10f ; %0.10f\n",s,avgdiam[s],stdevdiam[s]);
      printf("%d  ;  %0.10f ; %0.10f\n",s,avgdiam[s],stdevdiam[s]);
    }
  }
  
  // cleanup
  if(dataflag) {
    delete [] vals;
    if(complexflag)
      delete [] imag_vals;
  }
  if(maskflag) {
    delete [] maskvals;
  }
  delete [] maskverts;
  delete [] avgdiam;
  delete [] stdevdiam;
  if(outflag) fclose(fp);
  
  MsgPrintf("\n%s: finished\n",progname);
}


