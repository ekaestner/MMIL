/* surf_fwhm.c: read surface stats file and calculate 2d gaussian smoothness
      created: 08/04/04 DH
     last mod: 01/29/06 DH

   purpose:
     calculating 2d gaussian smoothness to input to alphasim
       alphasim is used to simulate clustering, to generate
       cluster-size thresholds for multiple comparison correction

   input:
     stats file (w file)

   output:
     stdout only

   acknowledgements:
     this code is adapted from AFNI's 3dFWHM
     (B.D. Ward, Medical College of Wisconsin 1997)
*/

#include "surflib.h"

#define MINARGC 3

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char indir[STRLEN]=".";
char surf[STRLEN]="smoothwm";
int complexflag = 0;
char real_infix[STRLEN]="_r";
char imag_infix[STRLEN]="_i";

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem  instem           omit extension, hemi: <instem>-rh.w\n");
  printf("    -subj    <str>            subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -indir   <str> [.]        input dir\n");
  printf("    -surf    <str> [smoothwm] surface used for area calculations\n");
  printf("    -hemi    <str> [rh]       hemisphere (rh or lh)\n");
  printf("    -                  use average distance\n");
  printf("    -complex                  complex data (real and imaginary)\n");
  printf("    -infixes  [_r _i]    r    real,imaginary infixes\n");
  printf("    -quiet                    suppress messages\n");
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
      if (MATCH(argv[i],"-instem") && i+1<argc) {
        strcpy(instem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-complex")){
        complexflag = 1;
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
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
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
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

int main(int argc, char **argv)
{
  int i,k,m, dfcount=0, ecode;
  char tempstr[STRLEN];
  double dist=0, avgneighbors=0,
        avgdist=0, stdevdist=0, stdevneighbors=0,
        mindist=BIGFLOAT, maxdist=-BIGFLOAT;
  float fwhm;
  int minneighbors=1000,maxneighbors=-1000,numneighbors=0;
  MRIS *mris;

  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,surf);

  /* read values */
  if(complexflag) {
    MsgPrintf("%s: reading input files\n",progname);
    sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,real_infix,hemi);
    ecode = MRISreadValues(mris,tempstr);
    if(ecode)
      ErrorExit(ecode,"%s: ### error reading value file %s\n",
                progname,tempstr);
    sprintf(tempstr,"%s/%s%s-%s.w",indir,instem,imag_infix,hemi);
    ecode = MRISreadImagValues(mris,tempstr);
    if(ecode!=NO_ERROR)
      ErrorExit(ecode,"%s: ### error reading value file %s\n",
                progname,tempstr);
  } else {
    // read input file
    MsgPrintf("%s: reading input file\n",progname);
    sprintf(tempstr,"%s/%s-%s.w",indir,instem,hemi);
    ecode = MRISreadValues(mris,tempstr);
    if(ecode)
      ErrorExit(ecode,"%s: ### error reading value file %s\n",
                progname,tempstr);
  }

  /* generate some info about surface */
  for (k=0;k<mris->nvertices;k++) {
    numneighbors = mris->vertices[k].vnum;
    if(minneighbors > numneighbors) minneighbors = numneighbors;
    if(maxneighbors < numneighbors) maxneighbors = numneighbors;
    avgneighbors += numneighbors;
    stdevneighbors += numneighbors*numneighbors;
    for (m=0;m<mris->vertices[k].vnum;m++) {
      i = mris->vertices[k].v[m];
      if(i<=k) continue; /* only count neighbor pair once */
      dist = mris->vertices[k].dist[m];
      if(dist<=0) continue;
      if(mindist > dist) mindist = dist;
      if(maxdist < dist) maxdist = dist;
      avgdist += dist;
      stdevdist += dist*dist;
      dfcount++;
    }
  }
  stdevdist = sqrt((stdevdist - (avgdist*avgdist)/dfcount)/(dfcount-1));
  avgdist /= dfcount;
  stdevneighbors = sqrt((stdevneighbors - 
              (avgneighbors*avgneighbors)/mris->nvertices)/(mris->nvertices-1));
  avgneighbors /= mris->nvertices;

  /* output results */
  if(!getQuiet()) {
    printf("%s: #### RESULTS #####\n",progname);
    printf("    Subject name: %s\n", subj);
    printf("    input dir: %s\n", indir);
    printf("    input file stem: %s\n", instem);
    if(complexflag)
      printf("    complex data\n");
  }

  /* calculate fwhm */
  if(complexflag)
    fwhm = MRISgaussCxFWHM(mris,0); /* not sparse */
  else  
    fwhm = MRISgaussFWHM(mris,0); /* not sparse */

  if(!getQuiet()) {
    printf("    Number of vertices = %d\n", mris->nvertices);
    printf("    Total number of neigbor relations = %d\n", dfcount);
    printf("    Average number of neighbors per vertex = %0.4f +- %0.4f\n",
                 avgneighbors,stdevneighbors);
    printf("    Minimum/Maximum number of neighbors per vertex = %d / %d\n",
                 minneighbors,maxneighbors);
    printf("    Average distance between neighboring vertices = %0.4f +- %0.4f mm\n",
                 avgdist,stdevdist);
    printf("    Minimum/Maximum distance between neighboring vertices = %0.4f / %0.4f\n",
                 mindist,maxdist);
    printf("\n");
    printf("    Gaussian filter width:\n");
    printf("           fwhm = %5.2f mm\n", fwhm);
  } else
    printf("%0.4f\n",fwhm);

  exit(0);
}

