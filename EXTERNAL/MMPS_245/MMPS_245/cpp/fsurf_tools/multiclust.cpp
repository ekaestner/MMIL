/* multiclust.cpp: read surface stats file and find clusters
      created: 12/22/04 DH
     last mod: 02/28/06 DH

   purpose:
     finding clusters using a sliding threshold

   input:
     stats file (w format)

   output:
     masked stats file(s) (w format)

   acknowledgements:
     Doug Greve's mri_surfcluster
*/

/* 
 how this program works:
 1.  input w file and paired lists of thresholds and cluster sizes
 2.  optionally apply a prethreshold (and cluster exclusion)
 3.  for each threshold, starting at highest threshold, find suprathreshold
     clusters and exclude clusters based on minimum size
 4.  if a new cluster overlaps with an older cluster, don't add it
     (exclude larger overlapping cluster assuming thresholds are given
     in ascending order)
     a vertex can belong to only one cluster
 5.  for each cluster, grow cluster until it touches another cluster
     or reaches edge of some threshold
 6.  output w file with all remaining clusters, with different value for
     each cluster (or different w file for each cluster)
*/

#include "surflib.h"
#include "clustlib.h"

#define MINARGC 3
#define SMALLFLOAT 1e-10
#define DEFOUTSTEM "multiclust-output"

#define MAXGROWSTEPS 2000

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char indir[STRLEN]=".";
char outdir[STRLEN]=".";
int smooth = 0;
int nthresh = 1;
int threshabsflag = 0;
float *thresh = NULL;
float *minarea = NULL;
float growthresh = SMALLFLOAT;
char surf[STRLEN]="smoothwm";
int summaryflag = 0;
int talflag = 0;
int fixmniflag = 1;
int growflag = 1;
int outputmasksflag = 0;
int prethreshflag = 0;

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem    instem     omit extension, hemi: <instem>-rh.w\n");
  printf("    -subj      subjname   subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem [multiclust] omit extension, hemi: <outstem>-rh.w\n");
  printf("    -indir      [.]       input dir\n");
  printf("    -outdir     [.]       output dir\n");
  printf("    -hemi       [rh]      hemisphere (rh or lh)\n");
  printf("    -smooth     [0]       number of pre-smoothing iterations\n");
  printf("    -nthresh    [1]       number of thresholds (precede -thresh and -minarea)\n");
  printf("    -thresh     [10^-10]  min. value for clusters  (one for each thresh)\n");
  printf("                          order list from most liberal to most convervative\n");
  printf("    -minarea    [10^-10]  min. cluster size (mm^2) (one for each thresh)\n");
  printf("    -growthresh [10^-10]  min. value for growing clusters\n");
  printf("    -nogrow               do not grow clusters after finding\n");
  printf("    -prethresh            apply first threshold and cluster exclusion\n");
  printf("                          to data before doing sliding threshold search\n");
  printf("    -threshabs            when thresholding, use absolute values\n");
  printf("    -surf      [smoothwm] surface used for area calculations\n");
  printf("    -summary              output text summary of each cluster\n");
  printf("    -talcoords            transform coords to talairach space\n");
  printf("    -nofixmni             do not fix mni coordinates (M.Brett's txfm)\n");
  printf("    -outputmasks          output mask w file for each cluster\n");
  printf("    -quiet                suppress messages\n");
  printf("\n");
}

inline int checkopt(string *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt = argv[i];
  return(0);
}

inline int checkopt(char *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  strcpy(opt,argv[i]);
  return(0);
}

inline int checkopt(int *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt=atoi(argv[i]);
  return(0);
}

inline long checkopt(long *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt=atol(argv[i]);
  return(0);
}

inline int checkopt(float *opt, int i, int argc, char **argv)
{
  if(i>=argc) return(1);
  *opt=atof(argv[i]);
  return(0);
}

inline void nooptexit(const char *optname)
{
  printf("%s: ### missing parameter value \"%s\"\n",progname, optname);
  exit(1);
}

void parse_args(int argc, char **argv)
{
  int i,j;
  char temp[STRLEN];
  char *pch;

  progname = argv[0];
  /* strip off path */
  pch = strtok(progname,"/");
  while (pch!=NULL) {
    strcpy(temp,pch);
    pch = strtok (NULL,"/");
  }
  strcpy(progname,temp);

  if (argc<MINARGC) {usage(); exit(0);}

  thresh = new float[nthresh];
  minarea = new float[nthresh];

  /* parse arguments */
  for (i=1;i<argc;i++) {
    strcpy(temp,argv[i]);
    if (MATCH(temp,"-instem")) {
      if(checkopt(instem,i+1,argc,argv)) nooptexit("instem"); i++;
    } else
    if (MATCH(temp,"-outstem")) {
      if(checkopt(outstem,i+1,argc,argv)) nooptexit("outstem"); i++;
    } else
    if (MATCH(temp,"-subj") || MATCH(argv[i],"-name")) {
      if(checkopt(subj,i+1,argc,argv)) nooptexit("subj"); i++;
    } else
    if (MATCH(temp,"-indir")) {
      if(checkopt(indir,i+1,argc,argv)) nooptexit("indir"); i++;
    } else
    if (MATCH(temp,"-outdir")) {
      if(checkopt(outdir,i+1,argc,argv)) nooptexit("outdir"); i++;
    } else
    if (MATCH(temp,"-hemi")) {
      if(checkopt(hemi,i+1,argc,argv)) nooptexit("hemi"); i++;
    } else
    if (MATCH(temp,"-smooth")) {
      if(checkopt(&smooth,i+1,argc,argv)) nooptexit("smooth"); i++;
    } else
    if (MATCH(temp,"-nthresh")) {
      if(checkopt(&nthresh,i+1,argc,argv)) nooptexit("nthresh"); i++;
      if(thresh!=NULL)
        delete [] thresh;
      if(nthresh<=0) nthresh=1;
      thresh = new float[nthresh];
      minarea = new float[nthresh];
    } else
    if (MATCH(temp,"-thresh")){
      if(i+nthresh<argc) {
        for(j=0;j<nthresh;j++)
          thresh[j]=atof(argv[++i]);
      } else
        nooptexit("thresh");
    } else
    if (MATCH(temp,"-minarea") && i+nthresh<argc){
      if(i+nthresh<argc) {
        for(j=0;j<nthresh;j++)
          minarea[j]=atof(argv[++i]);
      } else
        nooptexit("minarea");
    } else
    if (MATCH(temp,"-growthresh")) {
      if(checkopt(&growthresh,i+1,argc,argv)) nooptexit("growthresh"); i++;
    } else
    if (MATCH(temp,"-nogrow")) {
      growflag = 0;
    } else
    if (MATCH(temp,"-prethresh")) {
      prethreshflag = 1;
    } else
    if (MATCH(temp,"-threshabs")) {
      threshabsflag = 1;
    } else
    if (MATCH(temp,"-surf")) {
      if(checkopt(surf,i+1,argc,argv)) nooptexit("surf"); i++;
    } else
    if (MATCH(temp,"-summary")) {
      summaryflag = 1;
    } else
    if (MATCH(temp,"-talcoords")){
      talflag = 1;
    } else
    if (MATCH(temp,"-nofixmni")){
      fixmniflag = 0;
    } else
    if (MATCH(temp,"-fixmni")){
      fixmniflag = 1;
    } else
    if (MATCH(temp,"-outputmasks")) {
      outputmasksflag = 1;
    } else
    if (MATCH(temp,"-quiet")){
      setQuiet(1);
    } else
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(instem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### instem not specified ...quitting\n",
              progname);
  }
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    strcpy(outstem,DEFOUTSTEM);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
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
  for (j=0;j<nthresh;j++) {
    if(thresh[j]<=0) thresh[j]=SMALLFLOAT;
    if(minarea[j]<=0) minarea[j]=SMALLFLOAT;
  }
  if(growthresh<=0) growthresh=SMALLFLOAT;
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

int main(int argc, char **argv)
{
  int i,j,k,n,m,p,ci,nverts,t;
  int ecode,addclustflag;
  int *nogrowflags=NULL;
  char valfile[STRLEN];
  float threshval;
  MRIS *mris;
  ClusterList clusters,newclusters;
  Cluster tempclust;
  int growsteps = 0;

  parse_args(argc,argv);

  // load surface
  mris = openSurface(subj,hemi,surf);
  nverts = mris->nvertices;

  // load talairach transform
  if (talflag) {
    if(mris->transform_loaded) {
      if(!getQuiet()) printTransform(mris->linear_transform,"linear_transform");
    } else {
      MsgPrintf("%s: transform NOT loaded succesfully -- cannot do talairach transform\n", progname);
      talflag = 0;
    }
  }
  
  // read values
  sprintf(valfile,"%s/%s-%s.w",indir,instem,hemi);
  ecode = MRISreadValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,valfile);
  }

  // smooth values
  if(smooth>0)
    MRISsmoothValues(mris,smooth);

  // todo: sort threshold list (but how to also sort minarea list at same time?)


  // apply first threshold to data
  if(prethreshflag) {
    t = 0;
    MsgPrintf("%s: applying prethreshold with threshold %d (%0.2f)...\n",
      progname,t+1,thresh[t]);
    newclusters.findClusters(mris,thresh[t],threshabsflag);
    newclusters.calcStats(mris);
    MsgPrintf("%s: setting excluded vals to zero\n",progname);
    for (k=0;k<mris->nvertices;k++) {
      ci = newclusters.getClusterIndex(k);
      if (newclusters.getClusterIndex(k)<0)
        mris->vertices[k].val = 0;
      else if (newclusters[ci].Area() < minarea[t])
        mris->vertices[k].val = 0;
    }
  }
  
  // find clusters for each threshold
  clusters.init(nverts);
  for (t=nthresh-1;t>=0;t--) {
    MsgPrintf("%s: finding clusters for threshold %d (%0.2f)...\n",
      progname,t+1,thresh[t]);
    newclusters.findClusters(mris,thresh[t],threshabsflag);
    newclusters.calcStats(mris);
    // add to cluster list
    n=m=p=0;
    for (i=0;i<newclusters.size();i++) {
      tempclust = newclusters[i];
      if (tempclust.size())
        n++;
      if (tempclust.Area()>=minarea[t]) {
        m++;
        addclustflag=1;
        for (j=0;j<tempclust.size();j++) {
          k=tempclust[j];
          if (clusters.getClusterIndex(k)!=-1) {
            addclustflag=0;
            break; // do not add to cluster list since this cluster
                   // contains a vertex that already belongs to a cluster
          }
        }
        if(addclustflag) {
          clusters.newCluster(tempclust);
          p++;
        }
      }
    }
    MsgPrintf("%s:   new clusters = %d\n",progname,n);
    MsgPrintf("%s:   new clusters as large as %0.2f mm^2 = %d\n",progname,minarea[t],m);
    MsgPrintf("%s:   new unique clusters = %d\n",progname,p);
  }
  MsgPrintf("%s: total number of clusters = %d\n",progname,clusters.size());

  // grow clusters with lowest threshold, keeping clusters separate
  if(!clusters.size()) {
    MsgPrintf("%s: nothing to write... quitting\n",progname);
    exit(0);
  }
  nogrowflags = new int[clusters.size()];
  for (i=0;i<clusters.size();i++)
    nogrowflags[i] = 0;
  threshval = growthresh;
  if(growflag)
    MsgPrintf("%s: growing clusters with thresh %0.2f\n",progname,growthresh);
  while(growflag) {
    growflag = 0;
    for (i=0;i<clusters.size();i++) {
      if(nogrowflags[i]) continue;
      nogrowflags[i]=1; // don't grow next time unless can grow this time
      tempclust.clear();
      for(j=0;j<clusters[i].size();j++) {
        k = clusters[i][j];
        for (m=0;m<mris->vertices[k].vnum;m++) {
          n = mris->vertices[k].v[m];
          // if neighbor does not belong to a cluster and is above growthresh, add to this one
          if (clusters.getClusterIndex(n)==-1) {
            if (applyThresh(mris->vertices[n].val,growthresh,threshabsflag)==1) {
              tempclust.addVertex(n);
              nogrowflags[i]=0; // keep growing next time
              growflag=1;
            }
          }
        }
      }
      clusters.addToCluster(&(clusters[i]),&tempclust);
    }
    growsteps++;
    MsgPrintf(".");
    if(growsteps>MAXGROWSTEPS)
      break;
  }
  MsgPrintf("\n");

  // transform coordinates and sort clusters by area
  MsgPrintf("%s: calculating cluster stats...\n",progname);
  clusters.calcStats(mris);
  if(talflag) {
    clusters.transformCoords(mris);
    if(fixmniflag) clusters.fixMNICoords();
  }
  clusters.sortByArea(mris);

  // assign different value to each cluster
  MsgPrintf("%s: assigning output values to vertices...\n",progname);
  MRISclearValues(mris);
  for(i=0;i<clusters.size();i++) {
    ci=clusters[i].Index();
    for(j=0;j<clusters[i].size();j++) {
      k=clusters[i][j];
      mris->vertices[k].val=ci+1;
    }
  }

  // write to file
  sprintf(valfile,"%s/%s-%s.w",outdir,outstem,hemi);
  MsgPrintf("%s: writing file %s...\n",progname,valfile);
  ecode = MRISwriteValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error writing value file %s\n",
              progname,valfile);
  }

  if(outputmasksflag) {
    MsgPrintf("%s: creating mask files...\n",progname);
    for(i=0;i<clusters.size();i++) {
      // create mask
      MRISclearValues(mris);
      for(j=0;j<clusters[i].size();j++) {
        k=clusters[i][j];
        mris->vertices[k].val=1;
      }

      // write to file
      sprintf(valfile,"%s/%s-%d-%s.w",outdir,outstem,i,hemi);
      MsgPrintf("%s: writing file %s...\n",progname,valfile);
      ecode = MRISwriteValues(mris,valfile);
      if(ecode!=NO_ERROR) {
        ErrorExit(ecode,"%s: ### error writing value file %s\n",
                  progname,valfile);
      }
    }
  }

  // output summary
  if(summaryflag) {
    printf("\nmulticlust output summary:\n");
    printf("ClustNo  #Vertices  Area(mm^2)  MinVal   MaxVal   AvgVal   StDev   CentX   CentY   CentZ    MaxX    MaxY    MaxZ\n");
    for (i=0;i<clusters.size();i++) {
      printf(" %3d %10d %10.2f %9.2f %8.2f %9.2f %8.2f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n",
        clusters[i].Index(),clusters[i].size(),clusters[i].Area(),
        clusters[i].Minval(),clusters[i].Maxval(),
        clusters[i].Avgval(),clusters[i].Stdev(),
        clusters[i].Cent_x(),clusters[i].Cent_y(),clusters[i].Cent_z(),
        clusters[i].Max_x(),clusters[i].Max_y(),clusters[i].Max_z());
    }
  }  
}

