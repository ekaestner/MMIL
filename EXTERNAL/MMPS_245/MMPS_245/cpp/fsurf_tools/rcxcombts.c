/* rcxcombts.c: calculates average (vector sum) or t-stats, 
   from real or complex  bfloat or wt (painted surface) data sets
      created: 08/11/03 DH
     last mod: 04/05/05 DH

   purpose:
     calculating multi-scan averages for real or complex time series
     (or fourier series) data

   input:
     bfloat volumes with depth
      or
     painted time series file (wt)

   output:
     bfloat volumes with depth
      or
     painted time series file (wt)

*/

#include "surflib.h"

#define MINARGC 10
#define REAL_INFIX "_r"
#define IMAG_INFIX "_i"


/* global variables */
static char *progname = NULL;

/* parameter defaults */
char **input;
char outdir[STRLEN]=".";
char outstem[STRLEN]=UNDEFSTR;
char *real_infix=NULL;
char *imag_infix=NULL;
char hemi[STRLEN]="rh";
int *revphase;
int nslices=0;
int xnum=0,ynum=0,depth=0;
int N=0;
int Nrev=0;
int complex = 0;
float offset=0.0;
char datatype[STRLEN]="vol";
char subj[STRLEN]=UNDEFSTR;
int voldataflag = 1;

void usage() {
  printf("\n");
  printf("Usage: %s -input indir1 instem1 indir2 instem2 ... \\\n",progname);
  printf("         [-revphase 0 1 ... -options]\n");
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -input     subsequent dirs and stems used for average\n");
  printf("      indir      each contains a data set\n");
  printf("      instem     omit infixes, suffix: <instem>_{r,i}_%%03d.bfloat\n");
  printf("                 or for surf datatype: <instem>_{r,i}-{rh,lh}.wt\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -datatype       [vol]      input and output data type: vol or surf\n");
  printf("    -subj         <subjname>   required if data type = surf\n");
  printf("    -hemi           [rh]       hemisphere (rh or lh) (surf data type only)\n");
  printf("    -nslices <int>  [0]        input slice count\n");
  printf("       if omitted, uses actual number of slices\n");
  printf("    -outdir  <str>  [.]        output dir (must exist)\n");
  printf("    -outstem <str>  [avg]      output file stem\n");
  printf("    -quiet                     suppress messages\n");
  printf("\n");
  printf("  Optional parameters for complex data only:\n");
  printf("    -complex                   for complex data sets (real + imaginary)\n");
  printf("    -revphase  subsequent flags (0 or 1) indicated phase-reversed data\n");
  printf("       if omitted, all data are treated as non-phase-reversed\n");
  printf("    -offset <float> [0.0]      phase offset value\n");
  printf("       subtraction of phase offset is done prior to phase reversal\n");
  printf("       use phase offset to correct for estimated hemodynamic delay\n");
  printf("    -infixes <r> <i>           real,imaginary infixes--default: %s %s\n",
         REAL_INFIX,IMAG_INFIX);
  printf("\n");
}

void parse_args(int argc, char **argv)
{
  int i,j;
  FILE *fp;
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
  
  /* allocate memory */
  revphase = (int *)malloc(argc*sizeof(int));            MTEST(revphase);
  input = (char **)malloc(argc*sizeof(char *));          MTEST(input);
  for(i=0;i<argc;i++)
    input[i] = (char *)malloc(STRLEN*sizeof(char));      MTEST(*input);
  real_infix = (char *)malloc(STRLEN*sizeof(char));      MTEST(real_infix);
  imag_infix = (char *)malloc(STRLEN*sizeof(char));      MTEST(imag_infix);

  strcpy(real_infix,REAL_INFIX);
  strcpy(imag_infix,IMAG_INFIX);

  /* parse arguments */
  for (i=1;i<argc;i++) {
    if (argv[i][0]=='-') {
      if (MATCH(argv[i],"-datatype") && i+1<argc) {
        strcpy(datatype,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-subj") || MATCH(argv[i],"-name")) && i+1<argc){
        strcpy(subj,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-hemi") && i+1<argc){
        strcpy(hemi,argv[i+1]); i++;
      } else
      if ((MATCH(argv[i],"-nslices") || MATCH(argv[i],"-n")) && i+1<argc) {
        nslices = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outdir") && i+1<argc) {
        strcpy(outdir,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-offset") && i+1<argc){
        offset = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-input")){
        if (++i>argc) {
          ErrorExit(ERROR_BADPARM,"%s: ### no datasets following -input ...quitting\n",
                    progname);
        }
        for (j=i;j<argc;j++) {
          if (argv[j][0]=='-') break;
          strcpy(input[N],argv[j]);
          N++;
        }
        i=j-1;
      } else
      if (MATCH(argv[i],"-complex")){
        complex = 1;
      } else
      if (MATCH(argv[i],"-revphase")){
        MsgPrintf("%s: searching for revphase flags...\n",progname);
        if (++i>argc) {
          ErrorExit(ERROR_BADPARM,"%s: ### no flags following -revphase ...quitting\n",
                    progname);
        }
        for (j=i;j<argc;j++) {
          if (argv[j][0]=='-') break;
          revphase[Nrev]=atoi(argv[j]);
          Nrev++;
        }
        i=j-1;
      } else
      if ((MATCH(argv[i],"-infixes")) && i+2<argc){
        strcpy(real_infix,argv[i+1]); strcpy(imag_infix,argv[i+2]); i+=2;
      } else
      if (MATCH(argv[i],"-quiet")){
        setQuiet(1);
      } else
        ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    } else {
      ErrorExit(ERROR_BADPARM,"%s: ### parse error opt: %s\n",progname,argv[i]);
    }
  }
  MsgPrintf("%s: starting to check arguments\n",progname);

  /* check arguments */
  if (MATCH(datatype,"vol"))
    voldataflag = 1;
  else if (MATCH(datatype,"surf")) {
    voldataflag = 0;
    if (MATCH(subj,UNDEFSTR))
      ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
                progname);
  } else
    ErrorExit(ERROR_BADPARM,"%s: ### datatype %s not supported ...quitting\n",
                progname,datatype);
  if (!voldataflag && !MATCH(hemi,"rh") && !MATCH(hemi,"lh"))
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",
              progname,hemi);

  if (MATCH(outstem,UNDEFSTR)) sprintf(outstem,"avg");
  if (N==0) {
    ErrorExit(ERROR_BADPARM,"%s: ### no input datasets specified ...quitting\n",
              progname);
  } 
  if (N%2) {
    ErrorExit(ERROR_BADPARM,"%s: ### odd input dirs+formats count ...quitting\n",
              progname);
  } 
  N /= 2;
  fp = fopen(outdir,"r");
  if (fp==NULL) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  } else fclose(fp);
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  if (Nrev!=0 && Nrev!=N) {
    ErrorExit(ERROR_BADPARM,"%s: ### number of revphase flags (%d) must equal number of datasets (%d) ...quitting\n",
              progname,Nrev,N);
  }
  for (i=0;i<Nrev;i++) {
    if (revphase[i] != 0 && revphase[i] != 1) {
      ErrorExit(ERROR_BADPARM,"%s: ### revphase flags (%d) must be 0 or 1 ...quitting\n",
                progname, revphase[i]);
    }
  }
  if (!complex) {
    free(real_infix);
    real_infix = NULL;
    free(imag_infix);
    imag_infix = NULL;
  }

  /* check that input files exist */
  if (voldataflag) {
    /* just check first real file for now */
    i = getNumBfloats(input[0],input[1],real_infix);
    if(i<=0) {
      if(real_infix==NULL)
        ErrorExit(ERROR_NOFILE,"%s: ### no bfloats found in %s/ with stem %s ...quitting\n",
               progname,input[0],input[1]);
      else
        ErrorExit(ERROR_NOFILE,"%s: ### no bfloats found in %s/ with stem %s%s ...quitting\n",
               progname,input[0],input[1],real_infix);
    }
    if(nslices==0 || nslices>i) {
      nslices = i;
      if(real_infix==NULL)
        MsgPrintf("%s: %d slices found in %s/ with stem %s\n",progname,i,input[0],input[1]);
      else
        MsgPrintf("%s: %d slices found in %s/ with stem %s%s\n",progname,i,input[0],input[1],real_infix);
    }
  } else {
    /* just check first real file for now */
    if (real_infix==NULL)
      sprintf(tempstr,"%s/%s-%s.wt",input[0],input[1],hemi);
    else
      sprintf(tempstr,"%s/%s%s-%s.wt",input[0],input[1],real_infix,hemi);
    if (!FileExists(tempstr))
      ErrorExit(ERROR_BADPARM,"%s: ### file %s not found ...quitting\n",
                progname,tempstr);
  }
  MsgPrintf("%s: finished parsing arguments\n",progname);
}

void read_bfloat_ts(float ***xyt, int slicenum, char *indir, char *instem, char *infix)
{
  char fname[STRLEN];

  if (infix==NULL)
    sprintf(fname,"%s/%s_%03d.bfloat",indir,instem,slicenum);
  else
    sprintf(fname,"%s/%s%s_%03d.bfloat",indir,instem,infix,slicenum);
  if(readFloatImageTS(fname,xyt,xnum,ynum,depth)==-1)
    ErrorExit(ERROR_BADPARM,"%s: ### error reading float image...quitting\n",
              progname);
}

void write_bfloat_ts(float ***xyt, int slicenum, char *outdir, char *outstem, char *infix)
{
  char imgname[STRLEN];
  char hdrname[STRLEN];

  if (infix == NULL) {
    sprintf(imgname,"%s/%s_%03d.bfloat",outdir,outstem,slicenum);
    sprintf(hdrname,"%s/%s_%03d.hdr",outdir,outstem,slicenum);
  } else {
    sprintf(imgname,"%s/%s%s_%03d.bfloat",outdir,outstem,infix,slicenum);
    sprintf(hdrname,"%s/%s%s_%03d.hdr",outdir,outstem,infix,slicenum);
  }
  if(writeFloatImageTS(imgname,xyt,xnum,ynum,depth)==-1) {
    ErrorExit(ERROR_BADPARM,"%s: ### error writing float image...quitting\n",
              progname);
  }
  if(writeFloatHeader(hdrname,xnum,ynum,depth,0)==-1) {
    ErrorExit(ERROR_BADPARM,"%s: ### error writing float header...quitting\n",
              progname);
  }
}

int main(int argc, char **argv)
{
  int x,y,z,t,i,v;
  float ***dat_r=NULL;
  float ***dat_i=NULL;
  float ***avg_r=NULL;
  float ***avg_i=NULL;
  float real,imag;
  float ampl,phas;
  float N_inv;
  
  float **surf_avg_r=NULL;
  float **surf_avg_i=NULL;
  int numsurfvals=0,numsurfvals_tmp=0,tpoints=0,tpoints_tmp=0;
  FILE *fp_in_r=NULL,*fp_in_i=NULL;
  FILE *fp_out_r=NULL,*fp_out_i=NULL;
  int vnum, vnum_tmp, nverts;
  char tempstr[STRLEN];

  parse_args(argc,argv);

  if (voldataflag) {
    /* get dimensions of dataset */
    if(readFloatHeaders(N,input,real_infix,&xnum,&ynum,&depth)==-1)
      ErrorExit(ERROR_BADFILE, "%s: ### error reading header files...quitting\n",
             progname);
        /* this assumes that if real files are there, so are imaginary files */

    /* allocate memory for data and calculations for one 2D+time slice*/
    MsgPrintf("%s: starting to allocate memory\n",progname);
    dat_r = (float ***)calloc(depth,sizeof(float **));         MTEST(dat_r);
    avg_r = (float ***)calloc(depth,sizeof(float **));         MTEST(avg_r);
    for (t=0;t<depth;t++) {
      dat_r[t] = (float **)calloc(ynum,sizeof(float *));       MTEST(*dat_r);
      avg_r[t] = (float **)calloc(ynum,sizeof(float *));       MTEST(*avg_r);
      for (y=0;y<ynum;y++) {
        dat_r[t][y] = (float *)calloc(xnum,sizeof(float));     MTEST(**dat_r);
        avg_r[t][y] = (float *)calloc(xnum,sizeof(float));     MTEST(**avg_r);
      }
    }
    if(complex) {
      dat_i = (float ***)calloc(depth,sizeof(float **));       MTEST(dat_i);
      avg_i = (float ***)calloc(depth,sizeof(float **));       MTEST(avg_i);
      for (t=0;t<depth;t++) {
        dat_i[t] = (float **)calloc(ynum,sizeof(float *));     MTEST(*dat_i);
        avg_i[t] = (float **)calloc(ynum,sizeof(float *));     MTEST(*avg_i);
        for (y=0;y<ynum;y++) {
          dat_i[t][y] = (float *)calloc(xnum,sizeof(float));   MTEST(**dat_i);
          avg_i[t][y] = (float *)calloc(xnum,sizeof(float));   MTEST(**avg_i);
        }
      }
    }
    MsgPrintf("%s: finished allocating memory\n",progname);

    /* calculate averages, one slice at a time */
    for (z=0;z<nslices;z++) {
      MsgPrintf("%s: calculating averages for slice %d...\n",progname,z);
      for (i=0;i<N;i++) {
        /* open file(s) for one scan */
        read_bfloat_ts(dat_r,z,input[2*i],input[2*i+1],real_infix);
        if(complex)
          read_bfloat_ts(dat_i,z,input[2*i],input[2*i+1],imag_infix);
        /* calculate sums for each pixel and timepoint */
        for (t=0;t<depth;t++) {
          for (y=0;y<ynum;y++) {
            for (x=0;x<xnum;x++) {
              real = dat_r[t][y][x];
              if(complex) {
                imag = dat_i[t][y][x];
                if (offset!=0.0) {
                  ampl = hypot(real,imag);
                  phas = atan2(imag,real) - offset*2.0*M_PI;
                  real = ampl*cos(phas);
                  imag = ampl*sin(phas);
                }
                if (Nrev!=0) if (revphase[i]) imag = -imag;
                avg_i[t][y][x] += imag;
              } /* if complex */
              avg_r[t][y][x] += real;
            } /* for x */
          } /* for y */
        } /* for t */
      } /* for i */
      N_inv = 1.0/N;
      /* calculate avgs for each pixel and timepoint */
      for (t=0;t<depth;t++) {
        for (y=0;y<ynum;y++) {
          for (x=0;x<xnum;x++) {
            avg_r[t][y][x] *= N_inv;
            if(complex)
              avg_i[t][y][x] *= N_inv;
          } /* for x */
        } /* for y */
      } /* for t */
      /* write output to file */
      MsgPrintf("%s: writing output to file...\n",progname);
      write_bfloat_ts(avg_r,z,outdir,outstem,real_infix);
      if(complex)
        write_bfloat_ts(avg_i,z,outdir,outstem,imag_infix);
    } /* for z */
  } else { /* surface data */
    nverts = nverticesSurf(subj,hemi);
    /* get dimensions of datasets and make sure they match */
    for (i=0;i<N;i++) {
      if (complex)
        sprintf(tempstr,"%s/%s%s-%s.wt",input[2*i],input[2*i+1],real_infix,hemi);
      else
        sprintf(tempstr,"%s/%s-%s.wt",input[2*i],input[2*i+1],hemi);

      MsgPrintf("%s: opening file %s...\n",progname,tempstr);
      fp_in_r = fopen(tempstr,"r");
      if (fp_in_r==NULL)
        ErrorExit(ERROR_BADFILE,"%s: ### could not open file %s\n",
                  progname,tempstr);
      fread2(&tpoints_tmp,fp_in_r);
      if (i==0) {
        tpoints = tpoints_tmp;
        MsgPrintf("%s: number of surface time points: %d\n",progname,tpoints);
      } else {
        if (tpoints != tpoints_tmp)
          ErrorExit(ERROR_BADFILE,"%s: ### number of time points (%d) in %s do not match first input (%d)\n",
                    progname,tpoints_tmp,tempstr,tpoints);
      }
      fclose(fp_in_r);

      if (complex) {
        sprintf(tempstr,"%s/%s%s-%s.wt",input[2*i],input[2*i+1],imag_infix,hemi);
        MsgPrintf("%s: opening file %s...\n",progname,tempstr);
        fp_in_i = fopen(tempstr,"r");
        if (fp_in_i==NULL)
          ErrorExit(ERROR_BADFILE,"%s: ### could not open file %s\n",
                    progname,tempstr);
        fread2(&tpoints_tmp,fp_in_i);
        if (tpoints != tpoints_tmp)
          ErrorExit(ERROR_BADFILE,"%s: ### number of time points (%d) in %s do not match first input (%d)\n",
                    progname,tpoints_tmp,tempstr,tpoints);
        if (numsurfvals != numsurfvals_tmp)
          ErrorExit(ERROR_BADFILE,"%s: ### number of values (%d) in %s do not match first input (%d)\n",
                    progname,numsurfvals_tmp,tempstr,numsurfvals);
        fclose(fp_in_i);
      }
    }

    /* allocate memory for data and calculations for one surface */
    MsgPrintf("%s: starting to allocate memory\n",progname);

    surf_avg_r = (float **)calloc(tpoints,sizeof(float *));    MTEST(surf_avg_r);
    if(complex)
      surf_avg_i = (float **)calloc(tpoints,sizeof(float *));  MTEST(surf_avg_i);
    for (t=0;t<tpoints;t++) {
      surf_avg_r[t] = (float *)calloc(nverts,sizeof(float));   MTEST(*surf_avg_r);
      if(complex)
        surf_avg_i[t] = (float *)calloc(nverts,sizeof(float)); MTEST(*surf_avg_i);
    }
    MsgPrintf("%s: finished allocating memory\n",progname);
    
    /* open and initialize output file */
    if (complex)
      sprintf(tempstr,"%s/%s%s-%s.wt",outdir,outstem,real_infix,hemi);
    else
      sprintf(tempstr,"%s/%s-%s.wt",outdir,outstem,hemi);

    MsgPrintf("%s: initializing output file %s...\n", progname,tempstr);
    fp_out_r = fopen(tempstr,"w");
    if (fp_out_r==NULL)
      ErrorExit(ERROR_BADFILE,"%s: ### cannot create file %s\n",progname,tempstr);

    MsgPrintf("%s: number of surface time points: %d\n",progname,tpoints);
    fwrite2(tpoints,fp_out_r);

    if (complex) {
      sprintf(tempstr,"%s/%s%s-%s.wt",outdir,outstem,imag_infix,hemi);
      MsgPrintf("%s: initializing output file %s...\n", progname,tempstr);
      fp_out_i = fopen(tempstr,"w");
      if (fp_out_i==NULL)
        ErrorExit(ERROR_BADFILE,"%s: ### cannot create file %s\n",progname,tempstr);
      fwrite2(tpoints,fp_out_i);
    }

    MsgPrintf("%s: reading files and calculating averages...\n", progname,tempstr);
    /* calculate averages, one file at a time */
    for (i=0;i<N;i++) {
      if (complex)
        sprintf(tempstr,"%s/%s%s-%s.wt",input[2*i],input[2*i+1],real_infix,hemi);
      else
        sprintf(tempstr,"%s/%s-%s.wt",input[2*i],input[2*i+1],hemi);

      MsgPrintf("%s: opening file %s...\n",progname,tempstr);
      fp_in_r = fopen(tempstr,"r");
      if (fp_in_r==NULL)
        ErrorExit(ERROR_BADFILE,"%s: ### could not open file %s\n",
                  progname,tempstr);
      fread2(&tpoints_tmp,fp_in_r);

      if (complex) {
        sprintf(tempstr,"%s/%s%s-%s.wt",input[2*i],input[2*i+1],imag_infix,hemi);
        MsgPrintf("%s: opening file %s...\n",progname,tempstr);
        fp_in_i = fopen(tempstr,"r");
        if (fp_in_i==NULL)
          ErrorExit(ERROR_BADFILE,"%s: ### could not open file %s\n",
                    progname,tempstr);
        fread2(&tpoints_tmp,fp_in_i);
      }

      for (t=0;t<tpoints;t++) {
        fread3(&numsurfvals,fp_in_r);
        MsgPrintf("%s: number of surface values for time point %d: %d\n",
                  progname,t,numsurfvals);
        if (complex) {
          fread3(&numsurfvals_tmp,fp_in_i);
          MsgPrintf("%s: number of surface imag values for time point %d: %d\n",
                    progname,t,numsurfvals_tmp);
          if(numsurfvals_tmp!=numsurfvals) {
            ErrorExit(ERROR_BADFILE,
            "%s: ### number of surface values imaginary (%d) do not match real...quitting\n",
            progname,numsurfvals_tmp);
          }
        }
        for (v=0;v<numsurfvals;v++) {
          fread3(&vnum,fp_in_r);
          if (vnum >= nverts)
            ErrorExit(ERROR_BADFILE,
            "%s: ### vertex number (%d) > number of vertices (%d) for this subject (%s)...quitting\n",
            progname,vnum,nverts,subj);
          fread4(&real,fp_in_r);
          if (complex) {
            fread3(&vnum_tmp,fp_in_i);
            if (vnum!=vnum_tmp)
              ErrorExit(ERROR_BADFILE,
              "%s: ### vertex number for imaginary does not match real...quitting\n",
              progname);
            fread4(&imag,fp_in_i);
            if (offset!=0.0) {
              ampl = hypot(real,imag);
              phas = atan2(imag,real) - offset*2.0*M_PI;
              real = ampl*cos(phas);
              imag = ampl*sin(phas);
            }
            if (Nrev!=0) if (revphase[i]) imag = -imag;
            surf_avg_i[t][vnum]+=imag;
          }
          surf_avg_r[t][vnum]+=real;
        }
      }

      fclose(fp_in_r);
      if (complex) fclose(fp_in_i);
    } /* for i (image)*/

    /* divide sums by N to get averages and write average vals to output file */
    N_inv = 1.0/N;
    for (t=0;t<tpoints;t++) {
      for (v=0,numsurfvals=0;v<nverts;v++)
        if (surf_avg_r[t][v]!=0 || (complex && surf_avg_i[t][v]!=0))
          numsurfvals++;
      fwrite3(numsurfvals,fp_out_r);
      if (complex)
        fwrite3(numsurfvals,fp_out_i);
      for (v=0;v<nverts;v++) {
        if (surf_avg_r[t][v]!=0 || (complex && surf_avg_i[t][v]!=0)) {
          surf_avg_r[t][v] *= N_inv;
          if (complex)
            surf_avg_i[t][v] *= N_inv;
          real = surf_avg_r[t][v];
          fwrite3(v,fp_out_r);
          fwrite4(real,fp_out_r);
          if (complex) {
            imag = surf_avg_i[t][v];
            fwrite3(v,fp_out_i);
            fwrite4(imag,fp_out_i);
          }
        }
      } /* for v */
    } /* for t */
    fclose(fp_out_r);
    if (complex) fclose(fp_out_i);
  }

  MsgPrintf("%s: finished.\n",progname);

  exit(0);
}
