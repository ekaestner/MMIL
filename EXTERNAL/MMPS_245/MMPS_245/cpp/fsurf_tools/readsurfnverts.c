/* readsurfnverts.c: open surface file and read number of vertices
      created: 09/15/03 DH
     last mod: 09/22/03 DH

   purpose:
     getting number of vertices for a given surface

   input:
     w value files

   output:
     stdout
*/
#include "surflib.h"

#define MINARGC 3
#define MATCH(A,B)   (!strcmp(A,B))
#define SQR(x) ((x)*(x))

char progname[STRLEN]="readsurfnverts";

/* parameter defaults */
char hemi[STRLEN]="unknown";
char subj[STRLEN]="unknown";

/* function prototypes */
void usage();
void parse_args(int argc, char **argv);

int main(int argc, char **argv)
{
  int nvertices = -1;

  if (argc<MINARGC) {usage(); exit(0);}
  parse_args(argc,argv);

  nvertices = nverticesSurf(subj,hemi);
  printf("%d\n",nvertices);
  exit(0);
}

void usage()
{
  printf("\n");
  printf("Usage: %s subjname hemi\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    subjname                   subject name\n");
  printf("    hemi                       hemisphere (rh or lh)\n");
  printf("  Optional parameters:\n");
  printf("    -debug                     print debug statements\n");
  printf("    -h                         print this usage statement\n");
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i,a=1;

  setQuiet(1);

  /* parse arguments */
  strcpy(subj,argv[a++]);  
  strcpy(hemi,argv[a++]);  
  for (i=a;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-debug")){
        setQuiet(0);
      }
      else if (MATCH(argv[i],"-h")){
        usage();
        exit(0);
      }
      else {printf("-1\n");exit(1);}
    } else {printf("-1\n");exit(1);}
  }

  /* check arguments */
  if (MATCH(subj,"unknown")) {
    if(getQuiet()) printf("-1\n");
    MsgPrintf("%s: ### subject name not specified\n",progname);
    exit(1);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    if(getQuiet()) printf("-1\n");
    MsgPrintf("%s: ### hemi must be rh or lh\n",progname);
    exit(1);
  }
}

