/* randsurfclust.cpp: generate random values for a surface, make histo of cluster size
      created: 12/23/04 DH
     last mod: 12/24/04 DH

   purpose:
     generating alpha's for different cluster sizes

   input:
     options only

   output:
     stdout only

   acknowledgements:
     this code is adapted from AFNI's 3dFWHM and AlphaSim
     (B.D. Ward, Medical College of Wisconsin 1997)
*/

#include "surflib.h"
#include "cdflib/cdflib.h"
#include "clustlib.h"
using namespace std;

#define MINARGC 3
#define RATIO(X,Y)        ((float)X/(float)Y)
#ifdef Darwin
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),INT_MAX)*((H)-(L)))
#else
#  define RAND_FLOAT(L,H)  ((L)+RATIO(random(),MAXINT)*((H)-(L)))
#endif
#define EPSILON 1.0e-6

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
float pval = 0.001;
int niter = 100;
int smooth = 10;
float alpha = 0.05;
int estimate_gfw = 0;
int N = 10;
int N2 = 10;
int tstatflag = 0;
int cxfstatflag = 0;
int zeromeanflag = 0; // todo: add as option

char surf[STRLEN]="smoothwm";

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -subj     subjname   subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -niter    [100]      number of iterations\n");
  printf("    -hemi     [rh]       hemisphere (rh or lh)\n");
  printf("    -tstat               generate t-stats from N independent samples\n");
  printf("              (default assumes normal distribution, implying large N)\n");
  printf("    -N        [10]       number of independent samples (e.g. subjects)\n");
  printf("              (ignored unless tstat or cxfstat switch is used)\n");
  printf("    -cxfstat             generate complex f-stats\n");
  printf("    -zeromean            assume zero mean for cxfstat\n");
  printf("    -pval     [0.001]    probability value\n");
  printf("    -smooth   [10]       number of smoothing steps\n");
  printf("    -alpha    [0.05]     corrected probability value\n");
  printf("    -egfw                estimate gaussian filter width\n");
  printf("    -surf     [smoothwm] surface used for area calculations\n");
  printf("    -quiet               suppress messages\n");
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
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-pval") && i+1<argc){
        pval = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-alpha") && i+1<argc){
        alpha = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-niter") && i+1<argc){
        niter = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-tstat")){
        tstatflag = 1;
      } else
      if (MATCH(argv[i],"-cxfstat")){
        cxfstatflag = 1;
      } else
      if (MATCH(argv[i],"-zeromean")){
        zeromeanflag = 1;
      } else
      if (MATCH(argv[i],"-N") && i+1<argc){
        N = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-smooth") && i+1<argc){
        smooth = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-egfw")){
        estimate_gfw = 1;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
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
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh"))
    ErrorExit(ERROR_BADPARM,"%s: ### hemi (%s) must be rh or lh\n",
              progname,hemi);
  if (pval < 0 || pval > 1)
    ErrorExit(ERROR_BADPARM,"%s: ### pval (%f) must be <1 and >0\n",
              progname,pval);
  if (alpha < 0 || alpha > 1)
    ErrorExit(ERROR_BADPARM,"%s: ### alpha (%f) must be <1 and >0\n",
              progname,alpha);
  if (N < 2)
    ErrorExit(ERROR_BADPARM,"%s: ### N (%d) must be >2\n",
              progname,N);
  if (smooth < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### smooth steps (%d) must be >0\n",
              progname,smooth);
  if (niter < 0)
    ErrorExit(ERROR_BADPARM,"%s: ### number of iterations (%d) must be >0\n",
              progname,niter);

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

  r   = sqrt(-2.0*log(u1));
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

void genRandSurfImag(MRIS *mris)
{
  int k, nverts, nvertsdiv2;
  float n1, n2;
  
  nverts = mris->nvertices;
  nvertsdiv2 = nverts/2;
  
  for (k=0;k<nvertsdiv2;k++) {
    normal(&n1,&n2);
    mris->vertices[k].imag_val = n1;
    mris->vertices[k+nvertsdiv2].imag_val = n2;
  }
  normal(&n1,&n2);
  mris->vertices[nverts-1].imag_val = n1;
}

float pcalc (MRIS *mris, float zthr)
{
  int k,nverts,pcount;
  float p;
  
  nverts = mris->nvertices;
  pcount = 0;
  for (k=0;k<nverts;k++)
    if (mris->vertices[k].val > zthr) pcount++;
  p = (float)pcount/nverts;
  return (p);
}

void calcThresh(MRIS *mris, float pthr, float *zthr)
{
  int k,nverts,which, status;
  float pact;
  double p, q, z, mean, sd, bound, dof, dof2;
  static long count=0;
  static double sum=0, sumsq=0;

  nverts = mris->nvertices;

  // update sums
  if (count < 1.0e+09) {
    count += nverts;
    for (k=0;k<nverts;k++) {
      sum += mris->vertices[k].val;
      sumsq += mris->vertices[k].val * mris->vertices[k].val;
    }
  }

  if(tstatflag) {
    // calculate t-threshold
    which = 2;
    dof = N-1;
    p = 1.0 - pthr;
    q = pthr;
    cdft (&which, &p, &q, &z, &dof, &status, &bound);
    *zthr = z;
  } else if(cxfstatflag) {
    // calculate F-threshold
    which = 2;
    dof = 2;
    dof2 = 2*N;
    p = 1.0 - pthr;
    q = pthr;
    cdff (&which, &p, &q, &z, &dof, &dof2, &status, &bound);
    *zthr = z;
  } else {
    // calculate z-threshold
    which = 2;
    p = 1.0 - pthr;
    q = pthr;
    mean = sum/count;
    sd = sqrt((sumsq-(sum*sum)/count)/(count-1));
    cdfnor (&which, &p, &q, &z, &mean, &sd, &status, &bound);
    *zthr = z;
  }

  if (!getQuiet()) {
    pact = pcalc(mris,*zthr);
    MsgPrintf("%s: p(nominal)=%f  thresh=%f  p(actual)=%f\n", progname,pthr,*zthr,pact);
  }
}

int main(int argc, char **argv)
{
  int i,j,k,n,nverts;
  float totalarea, clustarea, maxarea, avgarea;
  int t_area, c_area, m_area;
  int *freqs, *maxfreqs;
  float *probs, *alphas, *cumprops;
  float *sum=NULL, *sumsq=NULL;
  float *sum_i=NULL, *sumsq_i=NULL;
  float thresh;
  MRIS *mris;
  ClusterList clustlist;
  long numclusters=0;
  float totalrands;
  float gfw, avg_gfw=0;
  float single_gfw, single_avg_gfw=0;
  int threshabsflag = 0;
  float N_inv, Nm1_inv;

  parse_args(argc,argv);

  randseed();

  // load surface
  mris = openSurface(subj,hemi,surf);
  nverts = mris->nvertices;
  
  // calculate total surface area
  totalarea=0;
  for (k=0;k<nverts;k++)
    totalarea += mris->vertices[k].area;
  t_area = int(totalarea)+1;
  
  // initialize arrays, indexed by total surface area
  freqs = new int[t_area];
  maxfreqs = new int[t_area];
  probs = new float[t_area];
  alphas = new float[t_area];
  cumprops = new float[t_area];

  if(tstatflag || cxfstatflag) {
    sum = new float[nverts];
    sumsq = new float[nverts];
    N_inv = 1.0/(double)N;
    Nm1_inv = 1.0/((double)(N-1));
  }
  if(cxfstatflag) {
    sum_i = new float[nverts];
    sumsq_i = new float[nverts];
  }

  for (i=0;i<niter;i++) {
    MsgPrintf("%s: ## iteration %d ##\n",progname,i);

    // generate random values, smooth, estimate gfw, average, calc thresh
    MsgPrintf("%s: generating random values\n",progname);
    if(tstatflag) {
      for (k=0;k<nverts;k++) {
        sum[k] = 0;
        sumsq[k] = 0;
      }
      for (n=0;n<N;n++) {
        genRandSurf(mris);
        MRISsmoothValues(mris,smooth);
        if (estimate_gfw) {
          single_gfw=MRISgaussFWHM(mris,0); // not sparse
          single_avg_gfw+=single_gfw;
          MsgPrintf("%s: single subject gaussian fwhm = %0.2f mm\n",
            progname,single_gfw);
        }
        for (k=0;k<nverts;k++) {
          sum[k] += mris->vertices[k].val;
          sumsq[k] += mris->vertices[k].val * mris->vertices[k].val;
        }
      }
      // calculate tstats
      MsgPrintf("%s: calculating t-stats\n",progname);
      for (k=0;k<nverts;k++) {
        sum[k] *= N_inv;
        sumsq[k] = sumsq[k] - N*sum[k]*sum[k];
        sumsq[k] = (sumsq[k] > 0.0) ? sqrt(sumsq[k]*Nm1_inv*N_inv) : 0.0;
        mris->vertices[k].val = (sumsq[k] > 0.0) ? \
                                sum[k]/sumsq[k] : 0.0;
      }
      if (estimate_gfw) {
        gfw=MRISgaussFWHM(mris,0); // not sparse
        avg_gfw+=gfw;
        MsgPrintf("%s: t-stat gaussian fwhm = %0.2f mm\n",progname,gfw);
      }
    } else if(cxfstatflag) {
      for (k=0;k<nverts;k++) {
        sum[k] = 0;
        sumsq[k] = 0;
        sum_i[k] = 0;
        sumsq_i[k] = 0;
      }
      for (n=0;n<N;n++) {
        genRandSurf(mris);
        genRandSurfImag(mris);
        MRISsmoothComplexValues(mris,smooth);
        if (estimate_gfw) {
          single_gfw=MRISgaussCxFWHM(mris,0); // not sparse
          single_avg_gfw+=single_gfw;
          MsgPrintf("%s: single subject complex gaussian fwhm = %0.2f mm\n",
            progname,single_gfw);
        }
        for (k=0;k<nverts;k++) {
          sum[k] += mris->vertices[k].val;
          sumsq[k] += mris->vertices[k].val * mris->vertices[k].val;
          sum_i[k] += mris->vertices[k].imag_val;
          sumsq_i[k] += mris->vertices[k].imag_val * mris->vertices[k].imag_val;
        }
      }
      // calculate cxfstats
      MsgPrintf("%s: calculating complex F-stats\n",progname);
      for (k=0;k<nverts;k++) {
        // calcuate means
        sum[k] *= N_inv;
        sum_i[k] *= N_inv;

        // square means
        sum[k] *= sum[k];
        sum_i[k] *= sum_i[k];

        // calcuate variance
        if(zeromeanflag) {
          /* null hypothesis is zero mean, so could calculate 
             variance assuming mean is zero */
          sumsq[k] *= N_inv;
          sumsq_i[k] *= N_inv;
        } else {
          sumsq[k] = sumsq[k] - (double)N*sum[k];
          sumsq_i[k] = sumsq_i[k] - (double)N*sum_i[k];
          sumsq[k] *= Nm1_inv;
          sumsq_i[k] *= Nm1_inv;
        }
        sum[k] = sum[k] + sum_i[k]; // numerator
        sumsq[k] = (sumsq[k] + sumsq_i[k])*N_inv; // denomenator
        mris->vertices[k].val = (sumsq[k] > 0.0) ? \
                                sum[k]/sumsq[k] : 0.0;
      }
      if (estimate_gfw) {
        gfw=MRISgaussCxFWHM(mris,0); // not sparse
        avg_gfw+=gfw;
        MsgPrintf("%s: Complex F-stat gaussian fwhm = %0.2f mm\n",progname,gfw);
      }
    } else {
    // generate random values, smooth, estimate gfw, calc thresh
      genRandSurf(mris);
      MRISsmoothValues(mris,smooth);
      if (estimate_gfw) {
        gfw=MRISgaussFWHM(mris,0); // not sparse
        avg_gfw+=gfw;
        MsgPrintf("%s: gaussian fwhm = %0.2f mm\n",progname,gfw);
      }
    }
    calcThresh(mris,pval,&thresh);
    
    // find clusters
    MsgPrintf("%s: finding clusters\n",progname);
    clustlist.init(nverts);
    clustlist.findClusters(mris,thresh,threshabsflag);
    clustlist.calcStats(mris);
    
    // add to histo
    MsgPrintf("%s: counting clusters\n",progname);
    maxarea = avgarea = 0;
    numclusters = clustlist.size();
    for (j=0;j<numclusters;j++) {
      clustarea = clustlist[j].Area();
      avgarea += clustarea;
      if (clustarea > maxarea) maxarea = clustarea;
      c_area = (int)rint(clustarea);
      freqs[c_area]++;
    }
    m_area = (int)rint(maxarea);
    maxfreqs[m_area]++;
    if (numclusters)
      avgarea /= numclusters;
    else
      avgarea = 0;
    MsgPrintf("%s: num clusters = %d, avg area = %0.1f mm^2, max area = %0.1f mm^2\n",
      progname,numclusters,avgarea,maxarea);
  }

  // calculate probs
  for (i=0;i<t_area;i++) {
    numclusters += freqs[i];
  }
  totalrands = (float)niter*t_area;
  for (i=0;i<t_area;i++) {
    probs[i] = (float)i*freqs[i]/totalrands;
    alphas[i] = (float)maxfreqs[i]/niter;
    cumprops[i] = (float)freqs[i]/numclusters;
  }
  for (i=1;i<t_area-1;i++) {
    j = t_area - i;
    probs[j-1] += probs[j];
    alphas[j-1] += alphas[j];
    cumprops[i+1] += cumprops[i];
  }

  printf("\n%s: clust size  frequency   cum prop   p/vertex   max freq   alpha\n",
            progname);
  c_area = 0;
  for (i=1;i<t_area;i++) {
    if (alphas[i]<EPSILON)
      break;
    else {
      if (alphas[i] <= alpha && c_area == 0) c_area = i;
      printf("%s: %6d  %10d  %12.6f %11.8f %5d %13.6f\n",
                progname,i,freqs[i],cumprops[i],
                probs[i],maxfreqs[i],alphas[i]);
    }
  }

  printf("\n%s: ---- Output Summary ----\n",progname);
  printf("%s: options used:\n",progname);
  printf("%s:   nvertices = %d\n",progname,nverts);
  printf("%s:   total surface area = %0.1f mm^2\n",progname,mris->total_area);
  printf("%s:   total vertex surface area = %d mm^2\n",progname,t_area);
  printf("%s:   iterations = %d\n",progname,niter);
  printf("%s:   pval = %0.10f\n",progname,pval);
  printf("%s:   smooth steps = %d\n",progname,smooth);
  printf("%s:   alpha = %0.4f\n",progname,alpha);
  printf("%s: result:\n",progname);
  if (estimate_gfw) {
    avg_gfw /= niter;
    if(tstatflag) {
      single_avg_gfw /= (niter*N);
      printf("%s:   single subject average gaussian fwhm = %0.2f mm\n",
             progname,single_avg_gfw);
      printf("%s:   t-stat average gaussian fwhm = %0.2f mm\n",
             progname,avg_gfw);
    } else if(cxfstatflag) {
      single_avg_gfw /= (niter*N);
      printf("%s:   single subject average complex gaussian fwhm = %0.2f mm\n",
             progname,single_avg_gfw);
      printf("%s:   complex F-stat average complex gaussian fwhm = %0.2f mm\n",
             progname,avg_gfw);
    } else {
      printf("%s:   average gaussian fwhm = %0.2f mm\n",
             progname,avg_gfw);
    }
  }
  if (c_area == 0) {
    printf("%s:   couldn't find cluster size for this alpha\n",progname);
    printf("%s:    -- need more iterations\n",progname);
  } else
    printf("%s:   cluster size for this alpha = %d\n",progname,c_area);
  MsgPrintf("\n%s: finished\n",progname);
}


