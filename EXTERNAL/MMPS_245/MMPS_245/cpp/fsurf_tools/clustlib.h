#ifndef CLUSTLIB_H
#define CLUSTLIB_H

// clustlib.h

#include "surflib.h"
#include <iostream>
#include <new>
#include <string>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <cctype>
#include <algorithm>
using namespace std;


short applyThresh(float val, float thresh, int threshabsflag); // return 1 if val is > thresh, -1 if < thresh
void initSurfClustIndices(MRIS *mris); // set undefval=-1 and marked=0 for each vertex

class Cluster {
private:
  int clustIndex; // index of cluster within a clusterlist
  std::vector<int> vertices; // dynamic array of vertex numbers
  float area;    // surface area (mm^2) of cluster
  float maxdist; // maximum distance between members of cluster (~ diameter)
  float minval;
  float maxval;
  float avgval;
  float stdev;
  float max_x;  // x coord of vertex with maxval
  float max_y;
  float max_z;
  int max_v; // vertex number of vertex with maxval
  float cent_x;  // x coord of center of mass
  float cent_y;
  float cent_z;

public:
  Cluster(); // constructor
  Cluster(int k); // constructor
  Cluster(int k, int ci); // constructor
  Cluster(const Cluster &clust); // copy constructor
  ~Cluster(); // destructor
  int &operator[](int i); // allows access to vertices vector
  Cluster operator=(Cluster clust); // copy a cluster with = operator
  int isMember(int k) const; // returns 1 if vertex k is a member of this cluster
  void addVertex(int k); // add vertex k to this cluster
  void rmVertex(int k); // remove vertex k from this cluster
  void addCluster(Cluster &clust); // add all vertices in clust to this cluster
  void growCluster(int k, MRIS *mris, float thresh, int threshabsflag);
                            // add all vertices neighboring k to cluster
  void calcStats(MRIS *mris); // calculate cluster surface area, avg val, coords, etc.
  int transformCoords(MRIS *mris);  // apply talairach transform to coords
  int fixMNICoords();  // transform MNI to talairach

  inline void clear() {vertices.clear();clustIndex=-1;}
  inline int size() const {return vertices.size();}
  inline void setIndex(int ci) {clustIndex = ci;}
  inline int Index() const {return clustIndex;}
  inline float Area() const {return area;}
  inline float MaxDist() const {return maxdist;}
  inline float Minval() const {return minval;}
  inline float Maxval() const {return maxval;}
  inline float Avgval() const {return avgval;}
  inline float Stdev() const {return stdev;}
  inline float Max_x() const {return max_x;}
  inline float Max_y() const {return max_y;}
  inline float Max_z() const {return max_z;}
  inline int   Max_v() const {return max_v;}
  inline float Cent_x() const {return cent_x;}
  inline float Cent_y() const {return cent_y;}
  inline float Cent_z() const {return cent_z;}
};

class ClusterList {
private:
  std::vector<Cluster> clusters; // dynamic array of clusters
  int *clustIndexList; // array of cluster indices, indexed by vertex number
  int nvertices; // number of vertices in current opened surface
public:
  ClusterList();
  ClusterList(int nverts);  // constructor
  ClusterList(const ClusterList &clustlist); // copy constructor
  ~ClusterList(); // destructor
  Cluster &operator[](int i); // allows access to clusters vector
  void init(int nverts); // clears list, resets clustIndexList
  int size() const; // returns size of clusters vector
  int newCluster(int k); // creates new cluster, returns clustIndex
  int newCluster(Cluster &clust); // creates new cluster from existing one, returns clustIndex
  Cluster * getCluster(int k); // returns cluster for a vertex
  int getClusterIndex(int k) const; // returns cluster index for a vertex
  void addToCluster(Cluster *clust1, Cluster *clust2); // adds clust2 to clust1
                                                       // then clears clust2
  void addToCluster(Cluster *clust, int k); // adds vertex k to clust
  void findClusters(MRIS *mris, float thresh, int threshabsflag);
  void calcStats(MRIS *mris); // calculate stats for each cluster
  void transformCoords(MRIS *mris); // apply talairach transform to coords
  void fixMNICoords();  // transform MNI to talairach
  void sortByArea(MRIS *mris);
};


#endif

