#include "surflib.h"
#include "clustlib.h"
using namespace std;

/* sections in this file: */
/* utils */
/* cluster */
/* clusterlist */

/* utils */
short applyThresh(float val, float thresh, int threshabsflag)
{
  if ((threshabsflag && (fabs(val) < fabs(thresh))) ||
      (!threshabsflag && (val < thresh)))
    return -1;
  else
    return 1;
}

void initSurfClustIndices(MRIS *mris) {
  int k;
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].marked = 0;
}

bool compareClustArea(const Cluster &A, const Cluster &B) {
  return A.Area() < B.Area();
}

bool compareClustSize(const Cluster &A, const Cluster &B) {
  return A.size() < B.size();
}

bool compareClustAvgVal(const Cluster &A, const Cluster &B) {
  return A.Avgval() < B.Avgval();
}

bool compareClustMaxVal(const Cluster &A, const Cluster &B) {
  return A.Maxval() < B.Maxval();
}

/*##########################################################################*/
/* cluster */

// constructors
Cluster::Cluster() {
  clustIndex = -1;
  area = 0;
  maxdist = 0;
  minval = 0;
  maxval = 0;
  avgval = 0;
  stdev = 0;
  max_x = 0;
  max_y = 0;
  max_z = 0;
  max_v = 0; // vertex number of max vertex
  cent_x = 0;
  cent_y = 0;
  cent_z = 0;
}

Cluster::Cluster(int k) {
  vertices.push_back(k);
  clustIndex = -1;
  area = 0;
  maxdist = 0;
  minval = 0;
  maxval = 0;
  avgval = 0;
  stdev = 0;
  max_x = 0;
  max_y = 0;
  max_z = 0;
  max_v = 0; // vertex number of max vertex
  cent_x = 0;
  cent_y = 0;
  cent_z = 0;
}

Cluster::Cluster(int k, int ci) {
  vertices.push_back(k);
  clustIndex = ci;
  area = 0;
  maxdist = 0;
  minval = 0;
  maxval = 0;
  avgval = 0;
  stdev = 0;
  max_x = 0;
  max_y = 0;
  max_z = 0;
  max_v = 0; // vertex number of max vertex
  cent_x = 0;
  cent_y = 0;
  cent_z = 0;
}

// copy constructor
Cluster::Cluster(const Cluster &clust) {
  vertices = clust.vertices;
  clustIndex = clust.clustIndex;
  area = clust.area;
  maxdist = clust.maxdist;
  minval = clust.minval;
  maxval = clust.maxval;
  avgval = clust.avgval;
  stdev = clust.stdev;
  max_x = clust.max_x;
  max_y = clust.max_y;
  max_z = clust.max_z;
  max_v = clust.max_v;
  cent_x = clust.cent_x;
  cent_y = clust.cent_y;
  cent_z = clust.cent_z;
}

// destructor
Cluster::~Cluster() {}

int &Cluster::operator[](int i) {
  if (i<0) {
    cout << "### cluster vertex index = -1 (undefined) ... quitting ###\n";
    exit(1);
  }    
  if (i>=size()) {
    cout << "### cluster vertices index overrun ... quitting ###\n";
    exit(1);
  }    
  return vertices[i];
}

Cluster Cluster::operator=(Cluster clust) {
  vertices = clust.vertices;
  clustIndex = clust.clustIndex;
  area = clust.area;
  maxdist = clust.maxdist;
  minval = clust.minval;
  maxval = clust.maxval;
  avgval = clust.avgval;
  stdev = clust.stdev;
  max_x = clust.max_x;
  max_y = clust.max_y;
  max_z = clust.max_z;
  max_v = clust.max_v;
  cent_x = clust.cent_x;
  cent_y = clust.cent_y;
  cent_z = clust.cent_z;
  return *this;
}

int Cluster::isMember(int k) const {
  int i,member=0;
  for (i=0;i<size();i++)
    if (k==vertices[i]) {
      member=1; // vertex is already a member of this cluster
      break;
    }
  return member;
}

void Cluster::addVertex(int k) {
  if (!isMember(k)) vertices.push_back(k);
}

void Cluster::rmVertex(int k) {
  int i;
  std::vector<int>::iterator p = vertices.begin();
  for (i=0;i<size();i++) {
    if (k==vertices[i]) {
      vertices.erase(p);
      break;
    }
    p++;
  }
}

void Cluster::addCluster(Cluster &clust) {
  int i,k;
  for (i=0;i<clust.size();i++) {
    k=clust[i];
    addVertex(k);
  }
}

void Cluster::growCluster(int k, MRIS *mris, float thresh, int threshabsflag) {
  VERTEX *v;
  int m,i;
  
  v = &mris->vertices[k];
  mris->vertices[k].marked = 1;
  for (m=0;m<v->vnum;m++) {
    i = v->v[m];
    if(mris->vertices[i].marked) continue;
    mris->vertices[i].marked = 
      applyThresh(mris->vertices[i].val,thresh,threshabsflag);
    if(mris->vertices[i].marked == 1) {
      vertices.push_back(i); // add vertex assuming it is not already a member
      growCluster(i,mris,thresh,threshabsflag);
    }
  }
}

void Cluster::calcStats(MRIS *mris) {
  VERTEX *v;
  int i,j,k,m,nverts;
  double val,sum=0,sumsq=0,abssum=0,dist,dx,dy,dz;
  avgval=stdev=area=maxdist=0;
  minval=BIGFLOAT;
  maxval=-BIGFLOAT;

  nverts=size();
  for (i=0;i<nverts;i++) {
    k=vertices[i];
    v = &mris->vertices[k];
    val = v->val;
    sum += val;
    abssum += fabs(val);
    sumsq += val*val;
    if (val>maxval) {
      maxval=val;
      max_x = v->x;
      max_y = v->y;
      max_z = v->z;
      max_v = k;
    }
    if (val<minval) minval=val;
    area += v->area;
    cent_x += fabs(val) * v->x;
    cent_y += fabs(val) * v->y;
    cent_z += fabs(val) * v->z;

    for(j=i+1;j<nverts;j++) {
      m=vertices[j];
      dx = (v->x - mris->vertices[m].x);
      dy = (v->x - mris->vertices[m].x);
      dz = (v->x - mris->vertices[m].x);
      dist = sqrt(dx*dx + dy*dy + dz*dz);
      if (dist > maxdist) maxdist = dist;
    }
  }
  avgval = sum/nverts;
  stdev = sqrt(sumsq/nverts-avgval*avgval);
  cent_x /= abssum;
  cent_y /= abssum;
  cent_z /= abssum;
}

int Cluster::transformCoords(MRIS *mris) {
  double x,y,z;
  
  if(transformPoint(mris->linear_transform,
                  cent_x, cent_y, cent_z, &x, &y, &z)==-1) {
    return(-1);
  }
  cent_x = x;
  cent_y = y;
  cent_z = z;

  if(transformPoint(mris->linear_transform,
                  max_x, max_y, max_z, &x, &y, &z)==-1) {
    return(-1);
  }
  max_x = x;
  max_y = y;
  max_z = z;

  return 0;
}

int Cluster::fixMNICoords() {
  double x,y,z;
  
  fixMNITal(cent_x, cent_y, cent_z, &x, &y, &z);
  cent_x = x;
  cent_y = y;
  cent_z = z;

  fixMNITal(max_x, max_y, max_z, &x, &y, &z);
  max_x = x;
  max_y = y;
  max_z = z;

  return 0;
}

/*##########################################################################*/
/* clusterlist */

// constructor
ClusterList::ClusterList() {
  nvertices = 0;
  clustIndexList = NULL;
}

ClusterList::ClusterList(int nverts) {
  int i;
  nvertices = nverts;
  clustIndexList = new int[nvertices];
  for(i=0;i<nvertices;i++)
    clustIndexList[i] = -1;  
}

// copy constructor
ClusterList::ClusterList(const ClusterList &clustlist) {
  int i;
  clusters = clustlist.clusters;
  nvertices = clustlist.nvertices;
  clustIndexList = new int[nvertices];
  for (i=0;i<nvertices;i++)
    clustIndexList[i] = clustlist.clustIndexList[i];
}

// destructor
ClusterList::~ClusterList() {delete [] clustIndexList;}

Cluster &ClusterList::operator[](int i) {
  if (i<0) {
    cout << "### clusterlist index = -1 (undefined) ... quitting ###\n";
    exit(1);
  }    
  if (i>=size()) {
    cout << "### clusterlist index overrun ... quitting ###\n";
    exit(1);
  }    
  return clusters[i];
}

void ClusterList::init(int nverts) {
  int i;
  
  clusters.clear();
  nvertices = nverts;
  delete [] clustIndexList;
  clustIndexList = new int[nvertices]; 
  for(i=0;i<nvertices;i++)
    clustIndexList[i] = -1;
}

int ClusterList::size() const {
  int mysize = clusters.size();
  return mysize;
}

int ClusterList::newCluster(int k) {
  int ci;
  ci = size();
  clustIndexList[k] = ci;
  clusters.push_back(Cluster(k,ci));
  return ci;
}

int ClusterList::newCluster(Cluster &clust) {
  int i,k,ci;
  ci = size();
  clusters.push_back(clust);
  for(i=0;i<clust.size();i++) {
    k=clust[i];
    clustIndexList[k] = ci;
  }
  return ci;
}

Cluster * ClusterList::getCluster(int k) {
  if(clustIndexList[k]==-1)
    return NULL;
  else
    return &clusters[clustIndexList[k]];
}

int ClusterList::getClusterIndex(int k) const {return clustIndexList[k];}

void ClusterList::addToCluster(Cluster *clust1, Cluster *clust2) {
  int ci1, ci2, i, k;
  ci1 = clust1->Index();
  ci2 = clust2->Index();
  if (ci1 != ci2) {
    clust1->addCluster(*clust2);
    for (i=0;i<clust2->size();i++) {
      k=(*clust2)[i];
      clustIndexList[k]=ci1;
    }
    clust2->clear();
  }
}

void ClusterList::addToCluster(Cluster *clust, int k) {
  int ci;
  ci = clust->Index();
  clust->addVertex(k);
  clustIndexList[k]=ci;
}

void ClusterList::findClusters(MRIS *mris, float thresh, int threshabsflag) {
  VERTEX *v;
  Cluster *clust;
  int i,j,k,ci;

  initSurfClustIndices(mris);
  init(mris->nvertices);
  for (k=0;k<mris->nvertices;k++) {
    v = &mris->vertices[k];
    // if vertex has already been checked (added to a cluster or rejected), ignore it
    if (!v->marked) {
      v->marked = applyThresh(v->val,thresh,threshabsflag);
      if (v->marked == 1) {
        ci = newCluster(k);
        clust = getCluster(k);
        clust->growCluster(k,mris,thresh,threshabsflag);
        for (j=0;j<clust->size();j++) {
          i=(*clust)[j];
          clustIndexList[i]=ci;
        }
      }
    }
  }
}

void ClusterList::calcStats(MRIS *mris) {
  int i;
  for (i=0;i<size();i++)
    clusters[i].calcStats(mris);
}

void ClusterList::transformCoords(MRIS *mris) {
  int i;
  for (i=0;i<size();i++)
    if(clusters[i].transformCoords(mris)==-1) break;
}

void ClusterList::fixMNICoords() {
  int i;
  for (i=0;i<size();i++)
    if(clusters[i].fixMNICoords()==-1) break;
}

void ClusterList::sortByArea(MRIS *mris) {
  int i,j,k;
  Cluster clust;

  // sort clusters
  std::stable_sort(clusters.begin(),clusters.end(),compareClustArea);

  // reindex clusters
  initSurfClustIndices(mris);
  for (i=0;i<size();i++) {
    if(i!=clusters[i].Index()) {
      clusters[i].setIndex(i);
      clust = clusters[i];
      for (j=0;j<clust.size();j++) {
        k=clust[j];
        clustIndexList[k]=i;
      }
    }
  }
}
