/* surf_info.c: read surface file and print out some info
      created: 05/21/06 DH
     last mod: 05/21/06 DH

   purpose:
     getting info about a surface

   input:
     subject name

   output:
     stdout only
*/

#include "surflib.h"

#define MINARGC 3

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char surf[STRLEN]="white";

/* functions */
void usage()
{
  printf("\n");
  printf("Usage: %s -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -subj    <str>            subject name (can be ico)\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -surf    <str> [white]    surface used\n");
  printf("    -hemi    <str> [rh]       hemisphere (rh or lh)\n");
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
      if (MATCH(argv[i],"-subj") && i+1<argc) {
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc) {
        strcpy(hemi,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-surf") && i+1<argc) {
        strcpy(surf,argv[i+1]); i++;
      } else
      {
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
      }
    } else {
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    }
  }
  /* check arguments */
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
}

int main(int argc, char **argv)
{
  int i,k,m, dfcount=0;
  double dist=0, avgdist=0, avgneighbors=0,
        stdevdist=0, stdevneighbors=0,
        mindist=BIGFLOAT, maxdist=-BIGFLOAT,
        area=0,avgarea=0,stdevarea=0,
        minarea=BIGFLOAT,maxarea=-BIGFLOAT,
        total_area=0;
  int minneighbors=1000,maxneighbors=-1000,numneighbors=0;
  MRIS *mris;

  parse_args(argc,argv);

  /* load surface */
  mris = openSurface(subj,hemi,surf);

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
    area = mris->vertices[k].area;
    if(minarea > area) minarea = area;
    if(maxarea < area) maxarea = area;
    avgarea += area;
    stdevarea += area*area;
  }
  stdevdist = sqrt((stdevdist - (avgdist*avgdist)/dfcount)/(dfcount-1));
  avgdist /= dfcount;
  stdevneighbors = sqrt((stdevneighbors - 
              (avgneighbors*avgneighbors)/mris->nvertices)/(mris->nvertices-1));
  avgneighbors /= mris->nvertices;
  total_area = avgarea;
  stdevarea = sqrt((stdevarea - (avgarea*avgarea)/mris->nvertices)/(mris->nvertices-1));
  avgarea /= mris->nvertices;

  /* output results */
  printf("%s: #### Output #####\n",progname);
  printf("    subject name: %s\n", subj);
  printf("    surface: %s\n", surf);
  printf("    hemi: %s\n", surf);
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
  printf("    Average area per vertex = %0.4f +- %0.4f mm^2\n",
               avgarea,stdevarea);
  printf("    Minimum/Maximum area per vertex = %0.4f / %0.4f\n",
               minarea,maxarea);
  printf("    Total surface area = %0.4f mm^2\n",
               total_area);

  exit(0);
}

