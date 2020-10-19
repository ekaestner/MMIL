/* retinopath.cpp: find middle path through an eccentricity band of a visual area
      created: 10/06/05 DH
     last mod: 12/09/05 DH

   purpose:
     choosing surface locations for MEG/EEG source modeling in visual areas

   input:
     polar angle and eccentricity retinotopy complex stats files (w file)

   output:
     stats file (w file)
     text output     
*/

#include <algorithm>
#include "surflib.h"
#include "clustlib.h"
using namespace std;

#define MINARGC 9
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"
#define MAX_ITERS  2000
#define MASK_VAL 10

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char polstem[STRLEN]=UNDEFSTR;
char eccstem[STRLEN]=UNDEFSTR;
char maskstem[STRLEN]=UNDEFSTR;
char outstem[STRLEN]=UNDEFSTR;
char poldir[STRLEN]=".";
char eccdir[STRLEN]=".";
char maskdir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char *real_infix=NULL;
char *imag_infix=NULL;
int real_eccen_flag = 0;
float r_eccen_thresh = 1.0;
int smooth = 10; // smoothing done to data before trunc phase
int postsmooth = 10; // smoothing done to center points
float polar_offset = 0; // phase corresponding to polar angle 0
float eccen_offset = 0; // phase corresponding to min eccentricity
float eccen_minvang = 0.2;// minimum eccentricity visual angle (deg)
float eccen_maxvang = 10; // maximum eccentricity visual angle (deg)
float polar_angle0 = 0;
float polar_angle1 = 90;
float eccen_degree0 = 0;
float eccen_degree1 = 1;
int rev_eccen_phase_flag = 0;
int rev_polar_phase_flag = 0;
float search_radius = 4; // used in findPath (millimeters)
int nseg_polar = 4; // number of patches into which to divide eccen band
int output_intermed_flag = 0;
float minweight = 0.5; // which vertices of smoothed cluster to output to text file
float dist_weight = 0.1;
int eccen_logtrans_flag = 1; // for phase-encoded ecc data, assume log transform
                             // relation between phase and eccentricity

/* functions */
int areNeighbors(int v1, int v2, MRIS *mris) {
  for(int m=0;m<mris->vertices[v1].vnum;m++)
    if(v2 == mris->vertices[v1].v[m]) return 1;
  return 0;
}

float vertexDist(int v1, int v2, MRIS *mris) {
  float dx,dy,dz;

  dx = mris->vertices[v1].x - mris->vertices[v2].x;
  dy = mris->vertices[v1].y - mris->vertices[v2].y;
  dz = mris->vertices[v1].z - mris->vertices[v2].z;

  return sqrt(dx*dx + dy*dy + dz*dz);  
}

float vertexDist(int k, float x, float y, float z, MRIS *mris) {
  float dx,dy,dz;

  dx = mris->vertices[k].x - x;
  dy = mris->vertices[k].y - y;
  dz = mris->vertices[k].z - z;

  return sqrt(dx*dx + dy*dy + dz*dz);  
}

int nearestVertex(int k, MRIS *mris, Cluster *roi) {
  int i,j,n;
  float dist, minDist=BIGFLOAT;

  n = -1;
  for(i=0;i<roi->size();i++) {
    j = (*roi)[i];
    if(j==k) { // nearest vertex besides itself
      continue;
    }
    dist = vertexDist(j,k,mris);
    if(dist < minDist) {
      minDist=dist;
      n = j;
    }
  }
  return n;  // returns -1 if no vertices in roi besides k
}

int nearestVertex(float x, float y, float z, MRIS *mris, Cluster *roi) {
  int i,j,n;
  float dist, minDist=BIGFLOAT;

  n = -1;
  for(i=0;i<roi->size();i++) {
    j = (*roi)[i];
    dist = vertexDist(j,x,y,z,mris);
    if(dist < minDist) {
      minDist=dist;
      n = j;
    }
  }
  return n;
}

int farthestVertex(int k, MRIS *mris, Cluster *roi) {
  int i,j,n;
  float dist, maxDist=-BIGFLOAT;

  n = -1;
  for(i=0;i<roi->size();i++) {
    j = (*roi)[i];
    if(j==k) {
      continue;
    }
    dist = vertexDist(j,k,mris);
    if(dist > maxDist) {
      maxDist=dist;
      n = j;
    }
  }
  return n;
}

// path class definition
class Path {
public:
  Cluster pathClust;
  int head;
  int tail;
  int connected;

  Path() {
    connected = 0;
  }

  // construct path by connecting between two vertices
  Path(int v1, int v2, MRIS *mris, Cluster *roi) {
    findPath(v1,v2,mris,roi);
  }

  // copy constructor
  Path(const Path &path) {
    pathClust = path.pathClust;
    head = path.head;
    tail = path.tail;
    connected = path.connected;
  }
  ~Path() {}
  int &operator[](int i) {
    if (i<0) {
      cout << "### path vertex index = -1 (undefined) ... quitting ###\n";
      exit(1);
    }    
    if (i>=size()) {
      cout << "### path vertices index overrun ... quitting ###\n";
      exit(1);
    }
    return pathClust[i];
  }
  bool operator<(const Path &path) const {
    return size() < path.size();
  }
  Path operator=(Path path) {
    pathClust = path.pathClust;
    head = path.head;
    tail = path.tail;
    connected = path.connected;
    return *this;
  }
  int size() const {return pathClust.size();}
  void clear() {
    head=tail=connected=0;
    pathClust.clear();
  }


  // find head and tail as two farthest vertices within roi
  void findHeadTailVertex(MRIS *mris, Cluster *roi) {
    int i,j,v1,v2;
    float dist,maxDist=BIGFLOAT;

    head=tail=-1;

    MsgPrintf("%s: roi has %d vertices\n",progname,roi->size());
    maxDist=-BIGFLOAT;
    for(i=0;i<roi->size();i++){
      v1 = (*roi)[i];
      for(j=i+1;j<roi->size();j++){
        v2 = (*roi)[j];
        dist = vertexDist(v1,v2,mris);
        if(dist > maxDist) {
          maxDist = dist;
          tail = v1;
          head = v2;
        }
      }   
    }
    MsgPrintf("%s: head vertex = %d, tail vertex = %d, distance = %0.2f\n",
              progname,head,tail,maxDist);
  }

  // find head and tail as two vertices with most different phase
  void findHeadTailVertex(MRIS *mris, Cluster *roi, float *polar_phase,
                          float dweight) {
    int i,j,v1,v2;
    float dist,maxDist=BIGFLOAT;
    float phaseDiff,maxPhaseDiff,phase1,phase2;
    float score,maxScore=0;

    head=tail=-1;
    maxPhaseDiff = M_PI;

    MsgPrintf("%s: roi has %d vertices\n",progname,roi->size());

    // go through once to find max distance
    maxDist=-BIGFLOAT;
    for(i=0;i<roi->size();i++){
      v1 = (*roi)[i];
      for(j=i+1;j<roi->size();j++){
        v2 = (*roi)[j];
        dist = vertexDist(v1,v2,mris);
        if(dist > maxDist) {
          maxDist = dist;
        }
      }   
    }

    // go through a second time, calculating score based on 
    // normalized distance and normalized phase difference
    for(i=0;i<roi->size();i++){
      v1 = (*roi)[i];
      phase1 = polar_phase[i];
      for(j=i+1;j<roi->size();j++){
        v2 = (*roi)[j];
        phase2 = polar_phase[j];
        phaseDiff = fabs(phase2-phase1);
        dist = vertexDist(v1,v2,mris);

        score = (phaseDiff/maxPhaseDiff) + dweight*(dist/maxDist);

        if(score > maxScore) {
          maxScore = score;
          tail = v1;
          head = v2;
        }
      }   
    }

    // make sure tail is upper field (positive phase), head is lower field (negative phase)
    if(polar_phase[tail]<polar_phase[head]) {
      v1=tail;
      tail=head;
      head=v1;
    }

    MsgPrintf("%s: tail vertex = %d, head vertex = %d, distance = %0.2f, polar phases = %0.2f and %0.2f\n",
              progname,tail,head,maxDist,polar_phase[tail],polar_phase[head]);
  }

  // (re-)initialize path by connecting v1 and v2
  void findPath(int v1, int v2, MRIS *mris, Cluster *roi) {
    int i,k,m,n,nk;
    float dist, minDist;
    Cluster tempClust, rejects;

    // find path from v1 to v2, using any vertices
    // but use cluster vertices as guideposts to direct travel
    
//    pathClust.clear();
    connected = 0;
    k = tail = v1;
    head = v2;
    pathClust.addVertex(tail);
    i = 0;
    // if unable to connect v1 and v2 within MAX_ITERS, connected = 0
    MsgPrintf("%s: finding path from tail vertex %d to head vertex %d\n",
      progname,tail,head);

    while (!connected && i<MAX_ITERS) {
      if(!areNeighbors(k,head,mris)) {
        // select roi vertices within search radius that
        // are closer to head than current vertex
        tempClust.clear();
        for(m=0;m<roi->size();m++) {
          n = (*roi)[m];
          if(n!=k && vertexDist(n,k,mris)<search_radius &&
             vertexDist(n,head,mris)<vertexDist(k,head,mris))
            tempClust.addVertex(n);
        }

        if(tempClust.size()) {
          MsgPrintf("%s: at v=%d, selected %d out of %d vertices to guide path\n",
            progname,k,tempClust.size(),roi->size());

          // find center of mass of selected vertices
          tempClust.calcStats(mris);
        } else
          MsgPrintf("%s: at v(%d), using only head vertex to guide path\n",
            progname,k,tempClust.size(),roi->size());

        dist = 0;
        minDist = BIGFLOAT;
        nk = k;
        for(m=0;m<mris->vertices[k].vnum;m++) {
          n = mris->vertices[k].v[m];
          if(rejects.isMember(n)) {
            MsgPrintf("%s: skipping reject vertex %d\n",progname,n);
            continue;
          }
          if(tempClust.size()) {
          // pick neighbor of k closest to tempClust center
            dist = vertexDist(n,
              tempClust.Cent_x(),tempClust.Cent_y(),tempClust.Cent_z(),
              mris);
          } else {
          // for some reason, no roi vertices within search radius
          // so just pick vertex that gets us closest to head
            dist = vertexDist(n,head,mris);
          }
          if(dist < minDist) {
            minDist = dist;
            nk = n; // next k
          }
        }
        if(k==nk) {
          MsgPrintf("%s: ### no more neighbors for vertex %d!\n",progname,k);
          break;
        }
        if(pathClust.isMember(nk)) {
          MsgPrintf("%s: ### path looped back from vertex %d to vertex %d\n",
            progname,k,nk);
          rejects.addVertex(nk); // don't go here again
        } else {
          pathClust.addVertex(k);
          k = nk;
        }
      } else {
        if(pathClust.isMember(k))
          MsgPrintf("%s: vertex %d is already a path member\n",
            progname,k);
        else
          MsgPrintf("%s: now adding vertex %d to path\n",
            progname,k);
        
        if(pathClust.isMember(head))
          MsgPrintf("%s: head vertex %d is already a path member\n",
            progname,head);
        else
          MsgPrintf("%s: now adding head vertex %d to path\n",
            progname,head);

        pathClust.addVertex(k);
        pathClust.addVertex(head);
        MsgPrintf("%s: connected to head %d from %d\n",
          progname,head,k);
        connected = 1;
      }
      i++;
    }
    if(i>=MAX_ITERS)
      MsgPrintf("%s: unable to connect vertex %d with vertex %d within %d iters\n",
        progname,v1,v2,MAX_ITERS);

    if(connected)
      MsgPrintf("%s: successfully filled gaps in path\n",progname);
    else
      MsgPrintf("%s: unable to fill gaps in path\n",progname);
  }

  void findPath(MRIS *mris, Cluster *roi) {
    findHeadTailVertex(mris,roi);
    findPath(tail,head,mris,roi);
  }

  void findPath(MRIS *mris, Cluster *roi, float *polar_phase, float dweight) {
    findHeadTailVertex(mris,roi,polar_phase,dweight);
    findPath(tail,head,mris,roi);
  }
}; // end path class definition

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
  void sortpaths() {
    stable_sort(paths.begin(),paths.end());
  }
}; // end pathList class definition

void usage()
{
  printf("\n");
  printf("Usage: %s -polstem polstem -eccstem eccstem \\\n",progname);
  printf("          -maskstem maskstem -subj subjname [options]\n");
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -polstem  polstem   omit extension, infixes, hemi: <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -eccstem  eccstem   omit extension, infixes, hemi: <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -maskstem maskstem         stem for w file defining visual area\n");
  printf("    -subj     subjname         subject name\n");
  printf("\n");
  printf("  Optional parameters for \"real\" eccentricity data:\n");
  printf("    -real_eccen                 use real (not complex) eccentricity data\n");
  printf("    -r_eccen_thresh [1.0]       minimum threshold value\n");
  printf("    note: \"real\" eccentricity data can be generated from\n");
  printf("      deconvolution of block-design or event-related responses\n");
  printf("      to ring(s) at fixed eccentricities\n");
  printf("\n");
  printf("  Other optional parameters:\n");
  printf("    -outstem       [retpath]   output file stem\n");
  printf("    -poldir        [.]         polar angle input dir\n");
  printf("    -eccdir        [.]         eccentricity input dir\n");
  printf("    -maskdir       [.]         dir containing mask file\n");
  printf("    -outdir        [.]         output dir\n");
  printf("    -hemi          [rh]        hemisphere (rh or lh)\n");
  printf("    -smooth        [10]        number of pre-smoothing steps\n");
  printf("    -postsmooth    [10]        number of post-smoothing steps\n");
  printf("    -polar_offset  [0.0]       polar phase for angle 0\n");
  printf("    -eccen_offset  [0.0]       eccen phase for 0 degrees eccentricity\n");
  printf("    -eccen_minvang [0.2]       min visual angle (degrees)\n");
  printf("    -eccen_maxvang [10.0]      max visual angle (degrees)\n");
  printf("    -polar_angle0  [0.0]       starting polar angle (degrees)\n");
  printf("    -polar_angle1  [90.0]      ending polar angle (degrees)\n");
  printf("    -eccen_degree0 [0.0]       starting eccentricity (degrees)\n");
  printf("    -eccen_degree1 [1.0]       ending eccentricity (degrees)\n");
  printf("    -search_radius [4.0]       search radius (mm) for path search\n");
  printf("    -nseg_polar    [4]         number of divisions across polar angle\n");
  printf("    -minweight     [0.5]       minimum for text output of smoothed cluster\n");
  printf("    -dist_weight   [0.1]       distance weight (relative to pol-phase diff)\n");
  printf("    -rev_polar_phase           reverse polar angle phase\n");
  printf("    -rev_eccen_phase           reverse eccentricity phase\n");
  printf("    -output_intermed           write intermediate w files\n");
  printf("    -infixes       [_r _i]     real,imaginary infixes\n");
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
      if (MATCH(argv[i],"-maskstem") && i+1<argc){
        strcpy(maskstem,argv[i+1]); i++;
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
      if (MATCH(argv[i],"-maskdir") && i+1<argc) {
        strcpy(maskdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-real_eccen")){
        real_eccen_flag = 1;
      } else
      if (MATCH(argv[i],"-r_eccen_thresh") && i+1<argc) {
        r_eccen_thresh = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-smooth") && i+1<argc) {
        smooth = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-postsmooth") && i+1<argc) {
        postsmooth = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-polar_offset") && i+1<argc) {
        polar_offset = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_offset") && i+1<argc) {
        eccen_offset = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_minvang") && i+1<argc) {
        eccen_minvang = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_maxvang") && i+1<argc) {
        eccen_maxvang = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_logtrans")){
        eccen_logtrans_flag = 1;
      } else
      if (MATCH(argv[i],"-polar_angle0") && i+1<argc) {
        polar_angle0 = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-polar_angle1") && i+1<argc) {
        polar_angle1 = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_degree0") && i+1<argc) {
        eccen_degree0 = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_degree1") && i+1<argc) {
        eccen_degree1 = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-search_radius") && i+1<argc) {
        search_radius = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-nseg_polar") && i+1<argc) {
        nseg_polar = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-rev_polar_phase")){
        rev_polar_phase_flag = 1;
      } else
      if (MATCH(argv[i],"-rev_eccen_phase")){
        rev_eccen_phase_flag = 1;
      } else
      if (MATCH(argv[i],"-output_intermed")){
        output_intermed_flag = 1;
      } else
      if (MATCH(argv[i],"-minweight") && i+1<argc) {
        minweight = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-dist_weight") && i+1<argc) {
        dist_weight = atof(argv[i+1]); i++;
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
  if (!FileExists(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### maskdir %s not found ...quitting\n",
              progname,maskdir);
  }
  if(!isadir(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,maskdir);
  }
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (smooth<0) {
    ErrorExit(ERROR_BADPARM,"%s: ### smooth steps must be >= 0 ...quitting\n",
              progname);
  }
  if (postsmooth<0) {
    ErrorExit(ERROR_BADPARM,"%s: ### postsmooth steps must be >= 0 ...quitting\n",
              progname);
  }
  if (eccen_minvang >= eccen_maxvang) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccen_minvang must be < eccen_maxvang ...quitting\n",
              progname);
  }
  if (eccen_degree0 > eccen_maxvang) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccen_degree0 must be <= eccen_maxvang ...quitting\n",
              progname);
  }
  if (eccen_degree1 > eccen_maxvang) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccen_degree1 must be <= eccen_maxvang ...quitting\n",
              progname);
  }
  if (eccen_degree0 < eccen_minvang) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccen_degree0 must be <= eccen_minvang ...quitting\n",
              progname);
  }
  if (eccen_degree1 < eccen_minvang) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccen_degree1 must be <= eccen_minvang ...quitting\n",
              progname);
  }
  if (search_radius<0) {
    ErrorExit(ERROR_BADPARM,"%s: ### search radius must be >= 0 ...quitting\n",
              progname);
  }
  if (nseg_polar<=0) {
    ErrorExit(ERROR_BADPARM,"%s: ### num polar segments (nseg_polar) must be > 0 ...quitting\n",
              progname);
  }
  if (minweight<0) {
    ErrorExit(ERROR_BADPARM,"%s: ### minweight must be >= 0 ...quitting\n",
              progname);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

void complex2polar(MRIS *mris)
{
  int k;
  float phase,amplitude;
  
  for (k=0;k<mris->nvertices;k++) {
    if (!mris->vertices[k].ripflag) {
      phase = atan2(mris->vertices[k].imag_val,mris->vertices[k].val);
      amplitude = hypot(mris->vertices[k].val,mris->vertices[k].imag_val);
      if(amplitude == 0) phase = 0;
      mris->vertices[k].val = amplitude;
      mris->vertices[k].imag_val = phase;
    }
  }
}

void truncphase(float *amp, float *phase, int n, float p0, float p1)
{
  int k;
  float p;

  // p0 and p1 must be between 0 and 1 (cycles)
  p0 *= 2.0*M_PI; 
  p1 *= 2.0*M_PI;
  for (k=0;k<n;k++) {
      p = phase[k];
      if (p < 0) p+= 2.0*M_PI;
      if (((p0 < p1) && (p < p0 || p > p1)) ||
          ((p0 > p1) && (p < p0 && p > p1)))
        amp[k] = phase[k] = 0;
  }
}

void readCxData(float *amplitude,float *phase,const char *dir,const char *stem,MRIS *mris) {
  int k;
  char tempstr[STRLEN];
  int ecode=NO_ERROR;

  MsgPrintf("%s: ### reading complex data\n",progname);
  sprintf(tempstr,"%s/%s%s-%s.w",dir,stem,real_infix,hemi);
  ecode = MRISreadValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s%s-%s.w",dir,stem,imag_infix,hemi);
  ecode = MRISreadImagValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  MRISsmoothComplexValues(mris,smooth);
  MsgPrintf("%s: finished reading complex data files\n",progname);

// convert complex to amplitude and phase
  MsgPrintf("%s: convert complex to amplitude and phase\n",progname);
  complex2polar(mris);

// copy values
  for(k=0;k<mris->nvertices;k++) {
    amplitude[k] = mris->vertices[k].val;
    phase[k]     = mris->vertices[k].imag_val;
  }
}

void readRealData(float *amplitude,const char *dir,const char *stem,MRIS *mris) {
  int k;
  char tempstr[STRLEN];
  int ecode=NO_ERROR;

  MsgPrintf("%s: ### reading real data\n",progname);
  sprintf(tempstr,"%s/%s-%s.w",dir,stem,hemi);
  ecode = MRISreadValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  MRISsmoothValues(mris,smooth);
  MsgPrintf("%s: finished reading real data file\n",progname);

// copy values
  for(k=0;k<mris->nvertices;k++) {
    amplitude[k] = mris->vertices[k].val;
  }
}


void writeVals(float *vals,const char *stem,MRIS *mris) {
  int ecode=NO_ERROR;
  char tempstr[STRLEN];

// output w file with values
  sprintf(tempstr,"%s/%s-%s-%s.w",outdir,outstem,stem,hemi);
  MsgPrintf("%s: writing values to file %s\n",progname,tempstr);
  ecode = writeSurfVals(tempstr,vals,mris->nvertices);
  MsgPrintf("%s: done writing values to file %s\n",progname,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output file %s\n",
              progname,tempstr);
}

void writeClust(Cluster clust,const char *stem,MRIS *mris) {
  int i,k;
  float *vals;

// output w file with values at cluster vertices
  vals = new float[mris->nvertices]; MTEST(vals);
  for(i=0;i<clust.size();i++) {
    k=clust[i];
    if(k>mris->nvertices)
      ErrorExit(NO_ERROR,"%s: ### cluster vertex index (%d) too high (max=%d)\n",
        progname,k,mris->nvertices);
    vals[k]=MASK_VAL;
  }

  writeVals(vals,stem,mris);
  delete [] vals;
}

void writeMRISVals(MRIS *mris, const char *stem) {
  char tempstr[STRLEN];
  int ecode=NO_ERROR;

  // output w file values from mris
  sprintf(tempstr,"%s/%s-%s-%s.w",outdir,outstem,stem,hemi);
  MsgPrintf("%s: writing values to file %s\n",progname,tempstr);
  ecode = MRISwriteValues(mris,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output file %s\n",
              progname,tempstr);
}

void getBorder(Cluster roi, Cluster *border, MRIS *mris) {
  int i,k,m,n,ecode=NO_ERROR;

  // find border of area roi 
  MsgPrintf("%s: finding border of roi\n",progname);
  border->clear();
  for(i=0;i<roi.size();i++) {
    k = roi[i];
    // check if neighbors are non-zero -- then must be a border vertex
    for (m=0;m<mris->vertices[k].vnum;m++) {
      n = mris->vertices[k].v[m];
      if (mris->vertices[n].val==0) {
        border->addVertex(k);
        break;
      }
    }
  }
  if(border->size()<1)
    ErrorExit(ecode,"%s: # error # no area border vertices found\n",
              progname);
  MsgPrintf("%s: roi border has %d vertices\n",progname,border->size());
}

void avgClustValue(Cluster clust,float *values, float *avg, float *stdev) {
  int i,k;
  float val;

  *avg=*stdev=0;

  for(i=0;i<clust.size();i++) {
    k = clust[i];
    val = values[k];
    *avg += val;
    *stdev += val*val;
  }
  *avg /= clust.size();
  *stdev = sqrt((*stdev)/clust.size()-(*avg)*(*avg));
}

int main(int argc, char **argv)
{
  int ecode=NO_ERROR;
  char tempstr[STRLEN];
  MRIS *mrisSph; //, *mrisSmWM;
  float *maskvals=NULL;
  float *polar_phase=NULL, *polar_amp=NULL;
  float *eccen_phase=NULL, *eccen_amp=NULL;
  float phase0, phase1;
  float ecc_phase,ecc_angle,pol_phase,pol_angle;
  int i,j,k,nverts,maxclust,maxclustsize,s;
  int cent_v;
  float cumLength, *pathLength, segLength;
  ClusterList clusters;
  Cluster area, area_border, upper_border, lower_border, tempClust;
  Path polpath, tempPath;
  PathList Segments;
  float weight, cweight;
  FILE *fp1, *fp2;
  float avg,stdev,val;

  parse_args(argc,argv);

  // load surfaces
//  mrisSmWM = openSurface(subj,hemi,"smoothwm"); // for accurate coordinates and normals
  mrisSph = openSurface(subj,hemi,"sphere"); // for convenient distance measurements
  nverts = mrisSph->nvertices;
  MsgPrintf("%s: finished opening surfaces\n",progname);

  // allocate memory
  polar_phase = new float[nverts]; MTEST(polar_phase);
  polar_amp   = new float[nverts]; MTEST(polar_amp);
  eccen_phase = new float[nverts]; MTEST(eccen_phase);
  eccen_amp   = new float[nverts]; MTEST(eccen_amp);
  maskvals    = new float[nverts]; MTEST(maskvals);

  // read polar angle
  readCxData(polar_amp,polar_phase,poldir,polstem,mrisSph);

  if(output_intermed_flag) {
    // output values for debugging
    writeVals(polar_amp,"pol-amp",mrisSph);
    writeVals(polar_phase,"pol-phase",mrisSph);
  }

  // truncate polar angle phases
  phase0 = polar_angle0/360 + polar_offset;
  phase1 = polar_angle1/360 + polar_offset;
  MsgPrintf("%s: trunc polar phase with p0=%0.2f and p1=%0.2f\n",
    progname,phase0,phase1);
  phase0 = fmod(double(phase0),1.0); // remove any cycles > 1
  phase1 = fmod(double(phase1),1.0);
  if(phase0 < 0) phase0 += 1;
  if(phase1 < 0) phase1 += 1;
  MsgPrintf("%s: actually trunc polar phase with p0=%0.2f and p1=%0.2f\n",
    progname,phase0,phase1);
  truncphase(polar_amp,polar_phase,nverts,phase0,phase1);

  if(output_intermed_flag) {
    // output values for debugging
    writeVals(polar_amp,"trunc-pol-amp",mrisSph);
    writeVals(polar_phase,"trunc-pol-phase",mrisSph);
  }

  // read eccentricity

  if(real_eccen_flag) {
    readRealData(eccen_amp,eccdir,eccstem,mrisSph);

    // threshold eccentricity data
    for(k=0;k<nverts;k++)
      if(eccen_amp[k]<r_eccen_thresh)
        eccen_amp[k] = 0;

    if(output_intermed_flag) {
      // output values for debugging
      writeVals(eccen_amp,"trunc-ecc-amp",mrisSph);
    }
  } else {
    readCxData(eccen_amp,eccen_phase,eccdir,eccstem,mrisSph);

    if(output_intermed_flag) {
      // output values for debugging
      writeVals(eccen_amp,"ecc-amp",mrisSph);
      writeVals(eccen_phase,"ecc-phase",mrisSph);
    }

    // truncate eccentricity phases
    if(eccen_logtrans_flag) {
      phase0 = eccen_offset +
        log(eccen_degree0/eccen_minvang)/log(eccen_maxvang/eccen_minvang);
      phase1 = eccen_offset +
        log(eccen_degree1/eccen_minvang)/log(eccen_maxvang/eccen_minvang);
    } else {
      phase0 = eccen_degree0/eccen_maxvang + eccen_offset;
      phase1 = eccen_degree1/eccen_maxvang + eccen_offset;
    }
    MsgPrintf("%s: trunc eccen phase with p0=%0.2f and p1=%0.2f\n",
      progname,phase0,phase1);
    phase0 -= (int)phase0; // subtract any cycles > 1
    phase1 -= (int)phase1;
    if(phase0 < 0) phase0 += 1;
    if(phase1 < 0) phase1 += 1;
    MsgPrintf("%s: actually trunc eccen phase with p0=%0.2f and p1=%0.2f\n",
      progname,phase0,phase1);
    truncphase(eccen_amp,eccen_phase,nverts,phase0,phase1);

    if(output_intermed_flag) {
      // output values for debugging
      writeVals(eccen_amp,"trunc-ecc-amp",mrisSph);
      writeVals(eccen_phase,"trunc-ecc-phase",mrisSph);
    }
  }
  
  // read mask file
  MsgPrintf("%s: ### reading mask file\n",progname);
  sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
  ecode = readSurfVals(tempstr,maskvals,nverts);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading mask file %s\n",
              progname,tempstr);

  // copy maskvals to surface
  for (k=0;k<nverts;k++) {
    if(maskvals[k])
      mrisSph->vertices[k].val=MASK_VAL;
    else
      mrisSph->vertices[k].val=0;
  }

  if(output_intermed_flag) {
    // output w file with trunc'd area mask
    writeMRISVals(mrisSph,"area-mask");
  }

  // trim away trunc'd phases
  for (k=0;k<nverts;k++) {
    if(!eccen_amp[k] || !polar_amp[k])
      mrisSph->vertices[k].val=0;
  }

  if(output_intermed_flag) {
    // output w file with trunc'd area mask
    writeMRISVals(mrisSph,"trunc-area-mask");
  }

  // find cluster defining overlap between area mask and eccen band
  MsgPrintf("%s: searching for clusters common to area and trunc'd phases\n",progname);
  clusters.init(nverts);
  clusters.findClusters(mrisSph,0.1,0);
  MsgPrintf("%s: %d cluster(s) found\n",progname,clusters.size());
  
  if(!clusters.size()) {
    MsgPrintf("%s: no overlap found... ignoring eccentricity\n",progname);
    // no overlap between area mask and eccentricity band
    // so use entire area mask (still trunc polar)

    // copy maskvals to surface (again)
    for (k=0;k<nverts;k++) {
      if(maskvals[k])
        mrisSph->vertices[k].val=MASK_VAL;
      else
        mrisSph->vertices[k].val=0;
    }

    // trim away trunc'd phases (only polar this time)
    for (k=0;k<nverts;k++) {
      if(!polar_amp[k])
        mrisSph->vertices[k].val=0;
    }

    // find cluster defining overlap between area mask and eccen band
    MsgPrintf("%s: searching again for clusters\n",progname);
    clusters.init(nverts);
    clusters.findClusters(mrisSph,0.1,0);
    MsgPrintf("%s: %d cluster(s) found\n",progname,clusters.size());

    // if still no overlap, quit
    if(!clusters.size()) {
      ErrorExit(NO_ERROR,"%s: no vertices found in area with range of defined polar angles... quitting\n",
        progname);
    }
  }
  
  // only use biggest cluster
  maxclust=maxclustsize=0;
  for(i=0;i<clusters.size();i++) {
    if(maxclustsize < clusters[i].size()) {
      maxclustsize = clusters[i].size();
      maxclust = i;
    }
    MsgPrintf("%s: cluster %d has %d vertices\n",progname,i,clusters[i].size());
  }
  MsgPrintf("%s: cluster %d is largest with %d vertices\n",
            progname,maxclust,maxclustsize);
  area = clusters[maxclust];

  // find border of area 
  MsgPrintf("%s: finding border of area\n",progname);
  getBorder(area,&area_border,mrisSph);
  // output w file with area border
  if(output_intermed_flag) {
    writeClust(area_border,"area_border",mrisSph);
  }

  // separate area_border into two groups based on polar phase
  // get average value
  avgClustValue(area_border,polar_phase,&avg,&stdev);
  MsgPrintf("%s: avg polar phase for area border = %0.3f += %0.3f\n",
    progname,avg,stdev);

  // split area_border into upper_border and lower_border
  MsgPrintf("%s: splitting area border into upper and lower\n",progname);
  for(i=0;i<area_border.size();i++) {
    k = area_border[i];
    val = polar_phase[k];
    if(fabs(avg)>M_PI/2) {
      if(val > avg)
        lower_border.addVertex(k);
      else
        upper_border.addVertex(k);
    } else {
      if(val < avg)
        lower_border.addVertex(k);
      else
        upper_border.addVertex(k);
    }
  }

  // get average value of upper border
  avgClustValue(upper_border,polar_phase,&avg,&stdev);
  MsgPrintf("%s: avg polar phase for upper_border = %0.3f += %0.3f\n",
    progname,avg,stdev);
  // output w file with upper_border
  if(output_intermed_flag) {
    writeClust(upper_border,"upper_border",mrisSph);
  }
  
  // get average value of lower border
  avgClustValue(lower_border,polar_phase,&avg,&stdev);
  MsgPrintf("%s: avg polar phase for lower_border = %0.3f += %0.3f\n",
    progname,avg,stdev);
  // output w file with lower_border
  if(output_intermed_flag) {
    writeClust(lower_border,"lower_border",mrisSph);
  }
  
  // find center vertex of cluster (force path through this point)
  area.calcStats(mrisSph);
  cent_v = nearestVertex(area.Cent_x(),area.Cent_y(),area.Cent_z(),
                         mrisSph,&area);
  MsgPrintf("%s: center vertex %d\n",progname,cent_v);

  // find center vertex of upper border = tail
  upper_border.calcStats(mrisSph);
  polpath.tail = nearestVertex(upper_border.Cent_x(),upper_border.Cent_y(),
                               upper_border.Cent_z(),mrisSph,&upper_border);
  MsgPrintf("%s: tail vertex %d\n",progname,polpath.tail);
  
  // find center vertex of lower border = head
  lower_border.calcStats(mrisSph);
  polpath.head = nearestVertex(lower_border.Cent_x(),lower_border.Cent_y(),
                               lower_border.Cent_z(),mrisSph,&lower_border);
  MsgPrintf("%s: head vertex %d\n",progname,polpath.head);

  MsgPrintf("%s: searching for path from tail to center\n",progname);
  tempPath.findPath(polpath.tail,cent_v,mrisSph,&area);
  polpath.pathClust.addCluster(tempPath.pathClust);
  if(tempPath.connected) polpath.connected = 1;
  MsgPrintf("%s: searching for path from center to head\n",progname);
  tempPath.findPath(cent_v,polpath.head,mrisSph,&area);
  polpath.pathClust.addCluster(tempPath.pathClust);
  if(tempPath.connected && polpath.connected) {
    polpath.connected = 1;
  } else {
    polpath.connected = 0;
  }
  MsgPrintf("%s: #### polpath ####\n",progname);
  for(i=0;i<polpath.size();i++) {
    MsgPrintf("%s: polpath vertex %d: %d\n",
      progname,i,polpath[i]);
  }

//  if(output_intermed_flag) {
    // output w file with full area mid-path
    writeClust(polpath.pathClust,"polpath",mrisSph);
//  }
  
  // calculate cummulative length as function of path position
  MsgPrintf("%s: calculating path length\n",progname);
  if(!polpath.size() || !polpath.connected) {
    ErrorExit(NO_ERROR,"%s: unable to continue without connected polpath... quitting\n",
      progname);
  }
  if(polpath[0]!=polpath.tail) {
    ErrorExit(NO_ERROR,"%s: tail should be first vertex of polpath... quitting\n",
      progname);
  }
  if(polpath[polpath.size()-1]!=polpath.head) {
    ErrorExit(NO_ERROR,"%s: head should be last vertex of polpath... quitting\n",
      progname);
  }
  pathLength = new float[polpath.size()];
  cumLength = 0;
  j = polpath.tail;
  MsgPrintf("%s: polpath tail = %d, head = %d\n",
    progname,polpath.tail,polpath.head);
  pathLength[0]=0;
  for(i=1;i<polpath.size();i++) {
    k = polpath[i];
//    MsgPrintf("%s: measuring distance between %d and %d\n",progname,j,k);
    if(!areNeighbors(j,k,mrisSph)) {
      ErrorExit(NO_ERROR,"%s: polpath is not really connected between %d and %d... quitting\n",
        progname,j,k);
    }
    cumLength += vertexDist(j,k,mrisSph);
    pathLength[i] = cumLength;
    j = k;  
  }
  MsgPrintf("%s: total path length = %0.2f mm with %d vertices\n",
    progname,cumLength,polpath.size());

  // divide the path into segments
  segLength = cumLength / nseg_polar;
  s = 1;
  tempPath.clear();
  for(i=0;i<polpath.size();i++) {
    k = polpath[i];

    if(pathLength[i] > segLength*s) {
//      MsgPrintf("%s: new segment\n",progname);
      s++;
      Segments.newPath(tempPath);
      tempPath.clear();
    }

//    MsgPrintf("%s: vertex %d, segment %d, pathlength=%0.2f, seglength=%0.2f\n",
//      progname,i,s,pathLength[i],segLength);

    tempPath.pathClust.addVertex(k);  
  }
//  MsgPrintf("%s: last segment\n",progname);
  Segments.newPath(tempPath);

  MsgPrintf("%s: divided path into %d segments\n",progname,Segments.size());
  
  if(output_intermed_flag) {
    // output w files with segment
    for(s=0;s<Segments.size();s++) {
      sprintf(tempstr,"segment%d",s);
      writeClust(Segments[s].pathClust,tempstr,mrisSph);
    }
  }

  // create text files for vertex numbers
  sprintf(tempstr,"%s/%s-centers-%s.txt",outdir,outstem,hemi);
  fp1 = fopen(tempstr,"w");
  if(fp1==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### can't create file %s\n",progname,tempstr);
  fprintf(fp1,"vnum\tecc-phase\tecc-vang\tpol-phase\tpol-angle\n");

  sprintf(tempstr,"%s/%s-clusters-%s.txt",outdir,outstem,hemi);
  fp2 = fopen(tempstr,"w");
  if(fp2==NULL)
    ErrorExit(ERROR_BADFILE,"%s: ### can't create file %s\n",progname,tempstr);
  fprintf(fp2,"locnum\tvnum\tweight\tecc-phase\tecc-vang\tpol-phase\tpol-angle\n");

  // find middle point for each segment
  for(s=0;s<Segments.size();s++) {
    tempPath = Segments[s];
    i = (int)rint(0.5*tempPath.size())-1;
    if(i<0) i=0;
    MsgPrintf("%s: choosing as center of segment vertex %d of %d\n",
      progname,i+1,tempPath.size());
    j = tempPath[i];
    tempClust.clear();
    tempClust.addVertex(j);
    // output w file with center
    sprintf(tempstr,"center%d",s);
    writeClust(tempClust,tempstr,mrisSph);

    // add center vertex to text file
    // convert phase to angles
    pol_phase = polar_phase[j]/(2*M_PI);
    if(pol_phase < 0) pol_phase += 1;
    pol_angle = (pol_phase - polar_offset)*360;
    if(!real_eccen_flag) {
      ecc_phase = eccen_phase[j]/(2*M_PI);
      if(ecc_phase < 0) ecc_phase += 1;
      ecc_angle = (ecc_phase - eccen_offset)*eccen_maxvang;
      // fprintf(fp,"vnum\tecc-phase\tecc-vang\tpol-phase\tpol-angle\n");
      fprintf(fp1,"%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n",
        j,ecc_phase,ecc_angle,pol_phase,pol_angle);
    } else {
      // fprintf(fp,"vnum\tpol-phase\tpol-angle\n");
      fprintf(fp1,"%d\t%0.2f\t%0.2f\n",
        j,pol_phase,pol_angle);
    }

    // copy to surface
    for(k=0;k<nverts;k++) {
      if(k==j) mrisSph->vertices[k].val = MASK_VAL*MASK_VAL;
      else mrisSph->vertices[k].val = 0;
    }
    // apply smoothing
    MRISsmoothValues(mrisSph,postsmooth);
    sprintf(tempstr,"center%d-smooth%d",s,postsmooth);
    writeMRISVals(mrisSph,tempstr);
    
    // write vertices to text file        
    cweight = mrisSph->vertices[j].val;
    for(k=0;k<nverts;k++) {
      weight = mrisSph->vertices[k].val/cweight;
      if(weight>minweight) {
        // convert phase to angles
        pol_phase = polar_phase[k]/(2*M_PI);
        if(pol_phase < 0) pol_phase += 1;
        pol_angle = (pol_phase - polar_offset)*360;
        if(!real_eccen_flag) {
          ecc_phase = eccen_phase[k]/(2*M_PI);
          if(ecc_phase < 0) ecc_phase += 1;
          ecc_angle = (ecc_phase - eccen_offset)*eccen_maxvang;
          //  fprintf(fp2,"locnum\tvnum\tweight\tecc-phase\tecc-vang\tpol-phase\tpol-angle\n");
          fprintf(fp2,"%d\t%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n",
            s,k,weight,ecc_phase,ecc_angle,pol_phase,pol_angle);
        } else {
          //  fprintf(fp2,"locnum\tvnum\tweight\tpol-phase\tpol-angle\n");
          fprintf(fp2,"%d\t%d\t%0.2f\t%0.2f\t%0.2f\n",
            s,k,weight,pol_phase,pol_angle);
        }
      }
    }
  }

  fclose(fp1);
  fclose(fp2);

  MsgPrintf("%s: finished\n",progname);
  exit(0);
}

