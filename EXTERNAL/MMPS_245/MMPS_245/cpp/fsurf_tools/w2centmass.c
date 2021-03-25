/* w2centmass.c: calculate center of mass from surface-painted values
      created: 12/08/05 DH
     last mod: 12/08/05 DH

   purpose:
     finding the vertex closest to the center of mass of non-zero vertices

   input:
     w value file

   output:
     stdout only
*/

#include "surflib.h"

#define MINARGC 2
#define BIGFLOAT 1e10

/* global variables */
static char *progname = NULL;

/* parameter defaults */
char instem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char hemi[STRLEN]="rh";
char subj[STRLEN]=UNDEFSTR;
char surf[STRLEN]="smoothwm";

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
  printf("    -surf          [smoothwm] surface used for area calculations\n");
  printf("    -hemi    <str> [rh]       hemisphere (rh or lh)\n");
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
  int k,cvert=0,ecode;
  double f,sum,x,y,z,dx,dy,dz,cx,cy,cz,dist,mindist=BIGFLOAT;
  char valfile[STRLEN];
  MRIS *mris;

  parse_args(argc,argv);

  /* set value file name */
  sprintf(valfile,"%s/%s-%s.w",indir,instem,hemi);

  /* load surface */
  mris = openSurface(subj,hemi,surf);

  /* read values */
  ecode = MRISreadValues(mris,valfile);
  if(ecode!=NO_ERROR) {
    ErrorExit(ecode,"%s: ### error reading value file %s\n",
              progname,valfile);
  }

  /* calculate center of mass */
  for(k=0;k<mris->nvertices;k++) {
    f=fabs(mris->vertices[k].val);
    if(f<=0) continue;
    x=mris->vertices[k].x;
    y=mris->vertices[k].y;
    z=mris->vertices[k].z;
    cx += f*x;
    cy += f*y;
    cz += f*z;
    sum += f;
  }
  cx /= sum;
  cy /= sum;
  cz /= sum;
  MsgPrintf("%s: center of mass coordinates = (%0.3f,%0.3f,%0.3f)\n",
            progname,cx,cy,cz);
  
  // identify vertex closest to center of mass
  for(k=0;k<mris->nvertices;k++) {
    x=mris->vertices[k].x;
    y=mris->vertices[k].y;
    z=mris->vertices[k].z;
    dx=x-cx;
    dy=y-cy;
    dz=z-cz;
    dist=sqrt(dx*dx+dy*dy+dz*dz);
    if(dist<mindist) {
      mindist=dist;
      cvert=k;
    }    
  }

  MsgPrintf("%s: center vertex = %d, distance from center = %0.1f mm\n",
            progname,cvert,mindist);


  exit(0);
}

