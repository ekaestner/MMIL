/* retinopatch.cpp: find patch defined by polar angle and eccentricity ranges in a visual area
      created: 10/06/05 DH
     last mod: 10/12/05 DH

   purpose:
     choosing surface locations for MEG/EEG source modeling in visual areas

   input:
     complex stats file (w file)

   output:
     complex stats file (w file)
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
char outstem[STRLEN]="retpatch";
char poldir[STRLEN]=".";
char eccdir[STRLEN]=".";
char maskdir[STRLEN]=".";
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char *real_infix=NULL;
char *imag_infix=NULL;
int smooth = 10;
float polar_offset = 0; // phase corresponding to polar angle 0 (in cycles)
                        //   i.e. phase is normalized so 1 => 2*M_PI
float eccen_offset = 0; // phase corresponding to eccentricity 0 degrees
float eccen_maxvang = 10;  // visual angle (degrees) corresponding to
                           //   phase = eccen_offset + 1
int rev_eccen_phase_flag = 0;
int rev_polar_phase_flag = 0;
float targ_pwidth = 1;  // polar angle width of visual target, in degrees of visual angle
float targ_ewidth = 1;  // eccentricity width of visual target, in degrees of visual angle
float targ_ecenter = 5; // center of visual target, in degrees of visual angle
float targ_offset = 0;  // polar angle offset (not visual angle) from horizontal meridian
float num_targs = 4;    // number of target locations (per 360 degrees) for which to define patches
                        //   note: only contralateral targets will get patches defined


// functions
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
  printf("  Optional parameters:\n");
  printf("    -outstem       [retpatch]  output file stem\n");
  printf("    -poldir        [.]         polar angle input dir\n");
  printf("    -eccdir        [.]         eccentricity input dir\n");
  printf("    -maskdir       [.]         dir containing mask file\n");
  printf("    -outdir        [.]         output dir\n");
  printf("    -hemi          [rh]        hemisphere (rh or lh)\n");
  printf("    -smooth        [10]        number of pre-smoothing steps\n");
  printf("    -polar_offset  [0.0]       polar phase for angle 0\n");
  printf("    -eccen_offset  [0.0]       eccen phase for 0 degrees eccentricity\n");
  printf("    -eccen_maxvang [10.0]      max eccentricity (visual angle degrees)\n");
  printf("    -targ_pwidth   [1.0]       polar angle width of visual target\n");
  printf("                                (visual angle degrees)\n");
  printf("    -targ_ewidth   [1.0]       eccentricity width of visual target\n");
  printf("                                (visual angle degrees)\n");
  printf("    -targ_ecenter  [5.0]       center of visual target\n");
  printf("                                (visual angle degrees)\n");
  printf("    -targ_offset   [0.0]       polar angle offset (degrees)\n");
  printf("    -num_targs     [4]         number of visual targets to define patches\n");
  printf("    -rev_polar_phase           reverse polar angle phase\n");
  printf("    -rev_eccen_phase           reverse eccentricity phase\n");
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
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-smooth") && i+1<argc) {
        smooth = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-polar_offset") && i+1<argc) {
        polar_offset = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_offset") && i+1<argc) {
        eccen_offset = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-eccen_maxvang") && i+1<argc) {
        eccen_maxvang = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-targ_pwidth") && i+1<argc) {
        targ_pwidth = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-targ_ewidth") && i+1<argc) {
        targ_ewidth = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-targ_ecenter") && i+1<argc) {
        targ_ecenter = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-targ_offset") && i+1<argc) {
        targ_offset = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-num_targs") && i+1<argc) {
        num_targs = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-rev_polar_phase")){
        rev_polar_phase_flag = 1;
      } else
      if (MATCH(argv[i],"-rev_eccen_phase")){
        rev_eccen_phase_flag = 1;
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
  if (targ_pwidth <= 0) {
    ErrorExit(ERROR_BADPARM,"%s: ### targ_pwidth must be > 0 ...quitting\n",
              progname);
  }
  if (targ_ewidth <= 0) {
    ErrorExit(ERROR_BADPARM,"%s: ### targ_pwidth must be > 0 ...quitting\n",
              progname);
  }
  if (targ_ecenter > eccen_maxvang) {
    ErrorExit(ERROR_BADPARM,"%s: ### eccen_ecenter must be <= eccen_maxvang ...quitting\n",
              progname);
  }
  if (num_targs <= 0) {
    ErrorExit(ERROR_BADPARM,"%s: ### num_targs must be > 0 ...quitting\n",
              progname);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

int findCenterVertex(MRIS *mris, Cluster *roi)
{
  int i,j,center=-1;
  float x,y,z,cx=0,cy=0,cz=0,dx,dy,dz,dist,minDist=BIGFLOAT;
  
  // calculate center coordinates
  MsgPrintf("%s: roi has %d vertices\n",progname,roi->size());

  for(i=0;i<roi->size();i++) {
    j=(*roi)[i];
    x=mris->vertices[j].x;
    y=mris->vertices[j].y;
    z=mris->vertices[j].z;
    cx += x;
    cy += y;
    cz += z;
  }
  cx /= roi->size();
  cy /= roi->size();
  cz /= roi->size();
  MsgPrintf("%s: roi center coordinates = (%0.3f,%0.3f,%0.3f)\n",
            progname,cx,cy,cz);
  
  // identify vertex closest to center
  for(i=0;i<roi->size();i++) {
    j=(*roi)[i];
    x=mris->vertices[j].x;
    y=mris->vertices[j].y;
    z=mris->vertices[j].z;
    dx=x-cx;
    dy=y-cy;
    dz=z-cz;
    dist=sqrt(dx*dx + dy*dy + dz*dz);
    if(dist<minDist) {
      minDist=dist;
      center=j;
    }    
  }

  x=mris->vertices[center].x;
  y=mris->vertices[center].y;
  z=mris->vertices[center].z;
  MsgPrintf("%s: center vertex (%d) coords = (%0.3f,%0.3f,%0.3f)\n",
    progname,center,x,y,z);
  MsgPrintf("%s: distance from actual center = %f\n",progname,minDist);
  
  return center;
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

  p0 *= 2.0*M_PI;
  p1 *= 2.0*M_PI;
  for (k=0;k<n;k++) {
      p = phase[k];
      if (p < 0) p+= 2.0*M_PI;
      if (((p0 < p1) && (p < p0 || p > p1)) ||
          ((p0 > p1) && (p < p0 && p > p1)))
        amp[k] = 0;
  }
}

int main(int argc, char **argv)
{
  int ecode=NO_ERROR;
  char tempstr[STRLEN];
  MRIS *mrisPol, *mrisEcc, *mrisPatch;
  float *maskvals=NULL;
  float *polar_phase=NULL, *polar_amp=NULL;
  float *eccen_phase=NULL, *eccen_amp=NULL;
  float phase0, phase1;
  int i,j,k,nverts,maxclust,maxclustsize,targ;
  ClusterList clusters;
  Cluster roi;
  float targ_incr_angle,targ_pcenter;
  float cyc_per_deg,vang_per_deg,cyc_per_vang;

  parse_args(argc,argv);

// calculate useful conversion factors
  /*
   circumference of circle (in visual angle degrees) 
     with radius targ_ecenter, C = 2.0*M_PI*targ_ecenter

   vang_per_deg (visual angle degree per polar angle degrees)
     = C/360 = targ_ecenter*2*M_PI/360

   phase is in cycles (i.e. 1 cycle = 2*PI radians)
     cyc_per_deg (cycles per degree) = 1/360

   cyc_per_vang (cycles per visual angle degree)
     = (cyc_per_deg)/(vang_per_deg)
     = (1/360)/(targ_ecenter*2*M_PI/360)
     = 360/(360*2*M_PI*targ_ecenter)
     = 1/(2*M_PI*targ_ecenter)
  */
  cyc_per_deg = 1.0/360;
  vang_per_deg = targ_ecenter*M_PI/180;
  cyc_per_vang = 1.0/(2*M_PI*targ_ecenter);

// load surfaces
  mrisPol   = openSurface(subj,hemi,"smoothwm");
  mrisEcc   = openSurface(subj,hemi,"smoothwm");
  mrisPatch = openSurface(subj,hemi,"smoothwm");
  nverts = mrisPatch->nvertices;
  MsgPrintf("%s: finished opening surfaces\n",progname);

// allocate memory
  polar_phase = new float[nverts]; MTEST(polar_phase);
  polar_amp   = new float[nverts]; MTEST(polar_amp);
  eccen_phase = new float[nverts]; MTEST(eccen_phase);
  eccen_amp   = new float[nverts]; MTEST(eccen_amp);
  maskvals    = new float[nverts]; MTEST(maskvals);

// read polar angle
  MsgPrintf("%s: reading polar angle files\n",progname);
  sprintf(tempstr,"%s/%s%s-%s.w",poldir,polstem,real_infix,hemi);
  ecode = MRISreadValues(mrisPol,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s%s-%s.w",poldir,polstem,imag_infix,hemi);
  ecode = MRISreadImagValues(mrisPol,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  MRISsmoothComplexValues(mrisPol,smooth);
  MsgPrintf("%s: finished reading polar angle files\n",progname);

// read eccentricity
  MsgPrintf("%s: reading eccentricity files\n",progname);
  sprintf(tempstr,"%s/%s%s-%s.w",eccdir,eccstem,real_infix,hemi);
  ecode = MRISreadValues(mrisEcc,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  sprintf(tempstr,"%s/%s%s-%s.w",eccdir,eccstem,imag_infix,hemi);
  ecode = MRISreadImagValues(mrisEcc,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,tempstr);
  MRISsmoothComplexValues(mrisEcc,smooth);
  MsgPrintf("%s: finished reading eccentricity files\n",progname);

// convert complex to amplitude and phase
  MsgPrintf("%s: convert complex to amplitude and phase\n",progname);
  complex2polar(mrisPol);
  complex2polar(mrisEcc);

// output w file with ecc phase
  sprintf(tempstr,"%s/%s-ecc-phase-%s.w",outdir,outstem,hemi);
  MsgPrintf("%s: writing file %s\n",progname,tempstr);
  ecode = MRISwriteImagValues(mrisEcc,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output file %s\n",
              progname,tempstr);
  
// output w file with ecc amplitude
  sprintf(tempstr,"%s/%s-ecc-amp-%s.w",outdir,outstem,hemi);
  MsgPrintf("%s: writing file %s\n",progname,tempstr);
  ecode = MRISwriteValues(mrisEcc,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output file %s\n",
              progname,tempstr);

// output w file with pol phase
  sprintf(tempstr,"%s/%s-pol-phase-%s.w",outdir,outstem,hemi);
  MsgPrintf("%s: writing file %s\n",progname,tempstr);
  ecode = MRISwriteImagValues(mrisPol,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output file %s\n",
              progname,tempstr);
  
// output w file with pol amplitude
  sprintf(tempstr,"%s/%s-pol-amp-%s.w",outdir,outstem,hemi);
  MsgPrintf("%s: writing file %s\n",progname,tempstr);
  ecode = MRISwriteValues(mrisPol,tempstr);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error writing output file %s\n",
              progname,tempstr);

// read mask file
  MsgPrintf("%s: reading mask file\n",progname);
  sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
  ecode = readSurfVals(tempstr,maskvals,nverts);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading mask file %s\n",
              progname,tempstr);

// copy and truncate eccentricity values
  // copy values
    for(k=0;k<nverts;k++) {
      eccen_amp[k]   = mrisEcc->vertices[k].val;
      eccen_phase[k] = mrisEcc->vertices[k].imag_val;
    }
  // truncate eccentricity phases
  MsgPrintf("%s: truncating eccentricity phases\n",progname);
  phase0 = (targ_ecenter-targ_ewidth/2)/eccen_maxvang + eccen_offset;
  phase1 = (targ_ecenter+targ_ewidth/2)/eccen_maxvang + eccen_offset;
  truncphase(eccen_amp,eccen_phase,nverts,phase0,phase1);

// define patches for num_targs
  targ_incr_angle = 360.0/num_targs;
  targ_pwidth *= cyc_per_vang; // convert from visual angle to cycles
  MsgPrintf("%s: targ pwidth = %0.3f cycles\n",
    progname,targ_pwidth);
  for(targ=0;targ<num_targs;targ++) {
    MsgPrintf("%s: finding patch for target %d\n",progname,targ);
    // copy values
    for(k=0;k<nverts;k++) {
      polar_amp[k]   = mrisPol->vertices[k].val;
      polar_phase[k] = mrisPol->vertices[k].imag_val;
    }
    // truncate polar angle phases
    targ_pcenter = (targ_offset+targ_incr_angle*targ)*cyc_per_deg;
    MsgPrintf("%s: targ %d pcenter = %0.3f cycles\n",
      progname,targ,targ_pcenter);
    phase0 = (targ_pcenter-targ_pwidth/2) + polar_offset;
    phase1 = (targ_pcenter+targ_pwidth/2) + polar_offset;
    MsgPrintf("%s: targ %d truncphase(%0.3f,%0.3f)\n",
      progname,targ,phase0,phase1);
    truncphase(polar_amp,polar_phase,nverts,phase0,phase1);

    // define area with intersection of polar and eccen phase ranges
    for (k=0;k<nverts;k++) {
      if(eccen_amp[k] && polar_amp[k] && maskvals[k])
        mrisPatch->vertices[k].val=MASK_VAL;
      else
        mrisPatch->vertices[k].val=0;
    }
    
    // output w file with roi
    sprintf(tempstr,"%s/%s-roi%d-%s.w",outdir,outstem,targ,hemi);
    MsgPrintf("%s: writing file %s\n",progname,tempstr);
    ecode = MRISwriteValues(mrisPatch,tempstr);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error writing output file %s\n",
                progname,tempstr);

    // search for clusters, pick largest (if more than one) 
    MsgPrintf("%s: searching for clusters\n",progname);
    clusters.init(nverts);
    clusters.findClusters(mrisPatch,0.1,0);
    if(!clusters.size()) {
      MsgPrintf("%s: no cluster(s) found for targ %d\n",progname,targ);
      continue;
    }
    MsgPrintf("%s: %d cluster(s) found for targ %d\n",
      progname,clusters.size(),targ);
    maxclustsize=-1;
    maxclust=0;
    for(i=0;i<clusters.size();i++) {
      if(maxclustsize < clusters[i].size()) {
        maxclustsize = clusters[i].size();
        maxclust = i;
      }
      MsgPrintf("%s: cluster %d has %d vertices\n",progname,i,clusters[i].size());
    }
    MsgPrintf("%s: cluster %d is largest with %d vertices\n",
              progname,maxclust,maxclustsize);
    roi = clusters[maxclust];

    // zero values outside of cluster
    for (k=0;k<nverts;k++) {
      mrisPatch->vertices[k].val=0;
    }
    for (i=0;i<roi.size();i++) {
      j=roi[i];
      mrisPatch->vertices[j].val=MASK_VAL; 
    }

    // output w file with roi
    sprintf(tempstr,"%s/%s-patch%d-%s.w",outdir,outstem,targ,hemi);
    MsgPrintf("%s: writing file %s\n",progname,tempstr);
    ecode = MRISwriteValues(mrisPatch,tempstr);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error writing output file %s\n",
                progname,tempstr);
  
    // todo: text output with vertex numbers, center coordinates, 
    //     polar angle, eccentricity, average normal?
  }

  MsgPrintf("%s: finished\n",progname);
  exit(0);
}

