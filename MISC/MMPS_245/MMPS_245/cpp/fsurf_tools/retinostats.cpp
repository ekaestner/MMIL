/* retinostats.cpp: calculate statistics about retinotopy data
      created: 02/23/06 DH
     last mod: 02/23/06 DH

   purpose:
     calculating statistics about retinotopy data within a surface ROI

   input:
     complex stats file (w file)
     mask file (w file)

   output:
     complex stats file (w file)
*/

#include "surflib.h"
using namespace std;

#define MINARGC 7

// global variables
static char *progname = NULL;

// parameter defaults
char instem[STRLEN]=UNDEFSTR;
char maskstem[STRLEN]=UNDEFSTR;
char indir[STRLEN]=".";
char maskdir[STRLEN]=".";
char outstem[STRLEN]=UNDEFSTR;
char outdir[STRLEN]=".";
char subj[STRLEN]=UNDEFSTR;
char hemi[STRLEN]="rh";
char real_infix[STRLEN]="_r";
char imag_infix[STRLEN]="_i";
char surf[STRLEN]="smoothwm";
float thresh = 0;
int smooth = 0;
int nbins_amp = 10;
int nbins_phase = 10;
float binsize_amp = 0;
int norm_histo_flag = 0;
int contra_phase_flag = 0;
int fit_Gauss_phase_flag = 0;
int fit_line_phase_flag = 0;

// possible calculations
int histo_amp_flag = 0;
int histo_phase_flag = 0;
int fwhm_amp_flag = 0;
int fwhm_phase_flag = 0;
int fwhm_cx_flag = 0;

// functions
void usage()
{
  printf("\n");
  printf("Usage: %s -instem instem -maskstem maskstem -subj subjname [options]\n",progname);
  printf("\n");
  printf("  Required parameters:\n");
  printf("    -instem   instem     omit extension, infixes, hemi: <instem>_{r,i}-{rh,lh}.w\n");
  printf("    -maskstem maskstem   stem for w file defining ROI\n");
  printf("    -subj     subjname   subject name (can be ico)\n");
  printf("\n");
  printf("  Optional calculations:\n");
  printf("    -histo_amp           calculate amplitude histogram\n");
  printf("    -histo_phase         calculate phase histogram\n");
  printf("    -fwhm_amp            calculate amplitude fwhm smoothness\n");
  printf("    -fwhm_phase          calculate phase fwhm smoothness\n");
  printf("    -fwhm_cx             calculate complex fwhm smoothness\n");
  printf("\n");
  printf("  Optional parameters:\n");
  printf("    -outstem  [retstats] text output file stem\n");
  printf("    -hemi        [rh]    hemisphere (rh or lh)\n");
  printf("    -indir       [.]     input dir\n");
  printf("    -maskdir     [.]     dir containing mask file\n");
  printf("    -outdir      [.]     output dir\n");
  printf("    -surf     [smoothwm] surface used for triangle metrics\n");
  printf("    -thresh      [0]     amplitude threshold applied before smoothing\n");
  printf("    -smooth      [0]     # smoothing steps applied before calculations\n");
  printf("    -nbins_amp   [10]    # amplitude bins\n");
  printf("    -binsize_amp [auto]  amplitude bin size (otherwise set automatically)\n");
  printf("    -nbins_phase [10]    # phase bins\n");
  printf("    -norm_histo          normalize histograms by mask surface area\n");
  printf("    -contra_phase        make phase=0 the center of contralateral space\n");
  printf("    -fit_Gauss_phase     fit Gaussian curve to phase histogram\n");
  printf("    -fit_line_phase      fit flat line to phase histogram\n");
  printf("    -infixes   [_r _i]   real,imaginary infixes\n");
  printf("    -quiet               suppress messages\n");
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
      if (MATCH(argv[i],"-maskstem") && i+1<argc){
        strcpy(maskstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-outstem") && i+1<argc){
        strcpy(outstem,argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-indir") && i+1<argc) {
        strcpy(indir,argv[i+1]); i++;
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
      if (MATCH(argv[i],"-thresh") && i+1<argc) {
        thresh = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-histo_amp")){
        histo_amp_flag = 1;
      } else
      if (MATCH(argv[i],"-histo_phase")){
        histo_phase_flag = 1;
      } else
      if (MATCH(argv[i],"-nbins_amp") && i+1<argc) {
        nbins_amp = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-nbins_phase") && i+1<argc) {
        nbins_phase = atoi(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-binsize_amp") && i+1<argc) {
        binsize_amp = atof(argv[i+1]); i++;
      } else
      if (MATCH(argv[i],"-fwhm_amp")){
        fwhm_amp_flag = 1;
      } else
      if (MATCH(argv[i],"-fwhm_phase")){
        fwhm_phase_flag = 1;
      } else
      if (MATCH(argv[i],"-fwhm_cx")){
        fwhm_cx_flag = 1;
      } else
      if (MATCH(argv[i],"-norm_histo")){
        norm_histo_flag = 1;
      } else
      if (MATCH(argv[i],"-contra_phase")){
        contra_phase_flag = 1;
      } else
      if (MATCH(argv[i],"-fit_Gauss_phase")){
        fit_Gauss_phase_flag = 1;
      } else
      if (MATCH(argv[i],"-fit_line_phase")){
        fit_line_phase_flag = 1;
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
  if (MATCH(maskstem,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### maskstem not specified ...quitting\n",
              progname);
  }
  if (MATCH(outstem,UNDEFSTR)) {
    sprintf(outstem,"retstats");
  }
  if (!FileExists(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### indir %s not found ...quitting\n",
              progname,indir);
  }
  if(!isadir(indir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,indir);
  }
  if (!FileExists(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### maskdir %s not found ...quitting\n",
              progname,maskdir);
  }
  if(!isadir(maskdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,maskdir);
  }
  if (!FileExists(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### outdir %s not found ...quitting\n",
              progname,outdir);
  }
  if(!isadir(outdir)) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s exists, but not a dir ...quitting\n",
              progname,outdir);
  }
  if (MATCH(subj,UNDEFSTR)) {
    ErrorExit(ERROR_BADPARM,"%s: ### subject name not specified ...quitting\n",
              progname);
  }
  if (!MATCH(hemi,"rh") && !MATCH(hemi,"lh")) {
    ErrorExit(ERROR_BADPARM,"%s: ### %s: hemi must be rh or lh\n",progname,hemi);
  }
  if (nbins_amp<1) nbins_amp = 1;
  if (nbins_phase<1) nbins_phase = 1;
  if (binsize_amp<0) binsize_amp = 0;

  MsgPrintf("%s: finished parsing arguments\n",progname);
}

int main(int argc, char **argv)
{
  int ecode=NO_ERROR;
  char tempstr[STRLEN];
  MRIS *mris;
  float *maskvals=NULL,*histo=NULL,*fithisto=NULL;
  int *mask=NULL,*maskverts=NULL;
  int i,j,k,m,b,nverts,nmaskverts=0,nthreshverts=0,dfcount=0;
  float amp,phase,real,imag;
  double amp_avg=0,amp_var=0,amp_min=BIGFLOAT,amp_max=-BIGFLOAT;
  double phase_avg=0,phase_var=0,phase_min=BIGFLOAT,phase_max=-BIGFLOAT;
  double r_avg=0,i_avg=0,r_var=0,i_var=0,cx_var=0;
  float bin_min,bin_max,bin_size,histo_offset;
  float dist,area,totalarea;
  double df_amp,df_amp_sum=0,df_amp_sumsq=0;
  double df_phase,df_phase_sum=0,df_phase_sumsq=0;
  double df_real,df_real_sum=0,df_real_sumsq=0;
  double df_imag,df_imag_sum=0,df_imag_sumsq=0;
  float dfvar,varratio,fwhm;
  float fit_parm,fit_sum,fit_chisq;
  FILE *fp;

  parse_args(argc,argv);

  // load surface
  MsgPrintf("%s: opening %s's %s %s surface\n",progname,subj,hemi,surf);
  mris = openSurface(subj,hemi,surf);
  nverts = mris->nvertices;
  MsgPrintf("%s: finished opening surface\n",progname);

  // read complex input files
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

  // read mask file
  MsgPrintf("%s: reading mask file\n",progname);
  maskvals = new float[nverts]; MTEST(maskvals);
  mask = new int[nverts]; MTEST(mask);
  sprintf(tempstr,"%s/%s-%s.w",maskdir,maskstem,hemi);
  ecode = readSurfVals(tempstr,maskvals,nverts);
  if(ecode!=NO_ERROR)
    ErrorExit(ecode,"%s: ### error reading mask file %s\n",
              progname,tempstr);
  MsgPrintf("%s: finished reading input files\n",progname);
  for (k=0;k<nverts;k++) {
    if(maskvals[k]>0) {
      nmaskverts++;
      mask[k]=1;
      mris->vertices[k].undefval=0;
    } else {
      mask[k]=0;
      mris->vertices[k].undefval=1;
    }
  }
  maskverts = new int[nmaskverts]; MTEST(maskverts);
  for (k=0,j=0;k<nverts;k++) {
    if(mask[k]>0) maskverts[j++]=k;
  }

  // apply amplitude threshold (to vertices inside mask only)
  if(thresh>0) {
    for (j=0;j<nmaskverts;j++) {
      k=maskverts[j];
      amp = hypot(mris->vertices[k].val,mris->vertices[k].imag_val);
      if(amp<thresh) {
        // remove from mask
        maskverts[j]=-1;
        mask[k]=0;
        mris->vertices[k].undefval=1;
      }
    }  
  }

  // apply smoothing (to suprathreshold vertices inside mask only)
  MRISsmoothComplexValuesSparse(mris,smooth);

  // convert complex to polar (amplitude and phase)
  for (j=0;j<nmaskverts;j++) {
    k=maskverts[j];
    if(k<0) continue;
    mris->vertices[k].val2 =
      hypot(mris->vertices[k].val,mris->vertices[k].imag_val);
    mris->vertices[k].val2bak = 
      atan2(mris->vertices[k].imag_val,mris->vertices[k].val)/(2.0*M_PI);
    if(contra_phase_flag) {
      if(MATCH(hemi,"rh")) mris->vertices[k].val2bak+=0.5;
      if      (mris->vertices[k].val2bak<-0.5) mris->vertices[k].val2bak+=1;
      else if (mris->vertices[k].val2bak>0.5)  mris->vertices[k].val2bak-=1;
    } else
      if(mris->vertices[k].val2bak<0) mris->vertices[k].val2bak+=1;
  }

  // open output file
  sprintf(tempstr,"%s/%s.txt",outdir,outstem);
  fp = fopen(tempstr,"w");
  if(fp==NULL)
    ErrorExit(ERROR_NOFILE,"%s: ### can't create file %s\n",progname,tempstr);

  // calculate basic stats (min,max,avg,var)
  MsgPrintf("%s: calculating basic stats...\n",progname);
  for (j=0;j<nmaskverts;j++) {
    k=maskverts[j];
    if(k<0) continue;
    real=mris->vertices[k].val;
    imag=mris->vertices[k].imag_val;
    r_avg+=real;
    i_avg+=imag;
    r_var+=real*real;
    i_var+=imag*imag;

    amp=mris->vertices[k].val2;
    phase=mris->vertices[k].val2bak;
    if(amp>amp_max) amp_max=amp;
    if(amp<amp_min) amp_min=amp;
    if(phase>phase_max) phase_max=phase;
    if(phase<phase_min) phase_min=phase;
    amp_avg+=amp;
    amp_var+=amp*amp;
    phase_avg+=phase;
    phase_var+=phase*phase;
    nthreshverts++;
  }
  r_var = (r_var - (r_avg*r_avg/nthreshverts))/(nthreshverts-1);  
  i_var = (i_var - (i_avg*i_avg/nthreshverts))/(nthreshverts-1);  
  cx_var = r_var + i_var;
  r_avg/=nthreshverts;
  i_avg/=nthreshverts;
  amp_var = (amp_var - (amp_avg*amp_avg/nthreshverts))/(nthreshverts-1);  
  phase_var = (phase_var - (phase_avg*phase_avg/nthreshverts))/(nthreshverts-1);  
  amp_avg/=nthreshverts;
  phase_avg/=nthreshverts;

  // summarize results
  fprintf(fp,"basic retinotopy stats\n");
  fprintf(fp,"nverts      \t%d\n",nverts);
  fprintf(fp,"nmaskverts  \t%d\n",nmaskverts);
  fprintf(fp,"nthreshverts\t%d\n",nthreshverts);
  fprintf(fp,"r_avg       \t%-10.4f\n",r_avg);
  fprintf(fp,"r_var       \t%-10.4f\n",r_var);
  fprintf(fp,"i_avg       \t%-10.4f\n",i_avg);
  fprintf(fp,"i_var       \t%-10.4f\n",i_var);
  fprintf(fp,"cx_var      \t%-10.4f\n",cx_var);
  fprintf(fp,"amp_min     \t%-10.4f\n",amp_min);
  fprintf(fp,"amp_max     \t%-10.4f\n",amp_max);
  fprintf(fp,"amp_avg     \t%-10.4f\n",amp_avg);
  fprintf(fp,"amp_var     \t%-10.4f\n",amp_var);
  fprintf(fp,"phase_min   \t%-10.4f\n",phase_min);
  fprintf(fp,"phase_max   \t%-10.4f\n",phase_max);
  fprintf(fp,"phase_avg   \t%-10.4f\n",phase_avg);
  fprintf(fp,"phase_var   \t%-10.4f\n",phase_var);
  fprintf(fp,"\n");

  // calculate histograms
  if(histo_amp_flag) {
    MsgPrintf("%s: calculating amplitude histogram...\n",progname);
    histo = new float[nbins_amp]; MTEST(histo);
    for(b=0;b<nbins_amp;b++) histo[b]=0;
    if(binsize_amp==0) {
      histo_offset=amp_min;
      bin_min=amp_min;
      bin_max=amp_max;
      if(bin_max-bin_min > 4*sqrt(amp_var))
        bin_max = amp_avg + 2*sqrt(amp_var);     
      bin_size = (bin_max-bin_min)/nbins_amp;
    } else {
      histo_offset = 0;
      bin_size = binsize_amp;
      bin_min = 0;
      bin_max = bin_size*nbins_amp;
    }
    MsgPrintf("%s: amplitude bin size = %0.4f\n",progname,bin_size);
    totalarea=0;
    for (j=0;j<nmaskverts;j++) {
      k=maskverts[j];
      if(k<0) continue;
      amp=mris->vertices[k].val2;
      b = (int)floor((amp-histo_offset)/bin_size);
      if(b>nbins_amp-1) b=nbins_amp-1;
      area = mris->vertices[k].area;
      if(area<=0) continue;
      histo[b]+=area;
      totalarea+=area;
    }
    // output results
    if(norm_histo_flag)
      fprintf(fp,"amplitude histogram (normalized to 100 by total area)\n");
    else
      fprintf(fp,"amplitude histogram\n");
    fprintf(fp,"total surface area \t%0.4f\n",totalarea);
    if(totalarea>0) {
      fprintf(fp,"min\tmax\tarea\n");
      for (b=0;b<nbins_amp;b++) {
        if(norm_histo_flag)
          histo[b]=100.0*histo[b]/totalarea;
        bin_min = histo_offset+bin_size*b;
        if(b==nbins_amp-1)
          bin_max = amp_max;
        bin_max = histo_offset+bin_size*(b+1);

        fprintf(fp,"%-10.4f\t%-10.4f\t%-10.4f\n",bin_min,bin_max,histo[b]);
      }
    }
    fprintf(fp,"\n");
    delete [] histo;    
  }
  
  
  if(histo_phase_flag) {
    MsgPrintf("%s: calculating phase histogram...\n",progname);
    histo = new float[nbins_phase]; MTEST(histo);
    for(b=0;b<nbins_phase;b++) histo[b]=0;
    if(contra_phase_flag) {
      histo_offset=-0.5;
      bin_min=-0.5;
      bin_max=0.5;
    } else {
      histo_offset=0;
      bin_min=0;
      bin_max=1;
    }
    bin_size = (bin_max-bin_min)/nbins_phase;
    MsgPrintf("%s: phase bin size = %0.4f\n",progname,bin_size);
    totalarea=0;
    for (j=0;j<nmaskverts;j++) {
      k=maskverts[j];
      if(k<0) continue;
      phase=mris->vertices[k].val2bak;
      b = (int)floor((phase-histo_offset)/bin_size);
      if(b>nbins_phase-1) b=nbins_phase-1; // should not happen with phase
      area = mris->vertices[k].area;
      if(area<=0) continue;
      histo[b]+=area;
      totalarea+=area;
    }
    // output results
    if(norm_histo_flag)
      fprintf(fp,"phase histogram (normalized to 100 by total area)\n");
    else
      fprintf(fp,"phase histogram\n");
    fprintf(fp,"total surface area \t%0.4f\n",totalarea);
    if(totalarea>0) {
      fprintf(fp,"min\tmax\tarea\n");
      for (b=0;b<nbins_phase;b++) {
        if(norm_histo_flag)
          histo[b]=100.0*histo[b]/totalarea;
        bin_min = histo_offset+bin_size*b;
        bin_max = histo_offset+bin_size*(b+1);
        fprintf(fp,"%-10.4f\t%-10.4f\t%-10.4f\n",bin_min,bin_max,histo[b]);
      }
    }
    fprintf(fp,"\n");

    if(fit_Gauss_phase_flag) {
      MsgPrintf("%s: fitting Gaussian to phase histogram...\n",progname);
      fithisto = new float[nbins_phase]; MTEST(fithisto);
      fit_sum = 0;
      for(b=0;b<nbins_phase;b++) {
        fit_parm = 1.0/sqrt(phase_var*2.0*M_PI);
        phase=histo_offset+bin_size*(b+0.5);
        fithisto[b] =
          fit_parm*exp(-(phase-phase_avg)*(phase-phase_avg)/(2.0*phase_var));
        fit_sum += fithisto[b];
      }
      if(norm_histo_flag)
        fit_sum = fit_sum/100.0;
      else
        fit_sum = fit_sum/totalarea;

      // output results
      if(norm_histo_flag)
        fprintf(fp,"Gaussian fit phase histogram (normalized to 100)\n");
      else
        fprintf(fp,"Gaussian fit phase histogram\n");
      fprintf(fp,"min\tmax\tarea\n");
      fit_chisq=0;
      for (b=0;b<nbins_phase;b++) {
        // normalize fitted histogram
        fithisto[b] /= fit_sum;
        // calculate difference with actual histo
        fit_chisq += (fithisto[b]-histo[b])*(fithisto[b]-histo[b]);
        // output row
        bin_min = histo_offset+bin_size*b;
        bin_max = histo_offset+bin_size*(b+1);
        fprintf(fp,"%-10.4f\t%-10.4f\t%-10.4f\n",bin_min,bin_max,fithisto[b]);
      }
      fprintf(fp,"\n");
      fprintf(fp,"Gaussian fit chi-squared = %0.6f\n",fit_chisq);
      fprintf(fp,"\n");

      delete [] fithisto;
    }

    if(fit_line_phase_flag) {
      MsgPrintf("%s: fitting flat line to phase histogram...\n",progname);
      fithisto = new float[nbins_phase]; MTEST(fithisto);
      fit_sum = 0;
      for(b=0;b<nbins_phase;b++) {
        fit_parm = 1.0;
        phase=histo_offset+bin_size*(b+0.5);
        fithisto[b] = fit_parm;
        fit_sum += fithisto[b];
      }
      if(norm_histo_flag)
        fit_sum = fit_sum/100.0;
      else
        fit_sum = fit_sum/totalarea;

      // output results
      if(norm_histo_flag)
        fprintf(fp,"line fit phase histogram (normalized to 100)\n");
      else
        fprintf(fp,"line fit phase histogram\n");
      fprintf(fp,"min\tmax\tarea\n");
      fit_chisq=0;
      for (b=0;b<nbins_phase;b++) {
        // normalize fitted histogram
        fithisto[b] /= fit_sum;
        // calculate difference with actual histo
        fit_chisq += (fithisto[b]-histo[b])*(fithisto[b]-histo[b]);
        // output row
        bin_min = histo_offset+bin_size*b;
        bin_max = histo_offset+bin_size*(b+1);
        fprintf(fp,"%-10.4f\t%-10.4f\t%-10.4f\n",bin_min,bin_max,fithisto[b]);
      }
      fprintf(fp,"\n");
      fprintf(fp,"line fit chi-squared = %0.6f\n",fit_chisq);
      fprintf(fp,"\n");

      delete [] fithisto;
    }

    delete [] histo;
  }
  
  // calculate smoothness
  if(fwhm_amp_flag || fwhm_phase_flag || fwhm_cx_flag) {
    MsgPrintf("%s: calculating fwhm smoothness\n",progname);
    fprintf(fp,"smoothness\n");
    dfcount=0;
    for (j=0;j<nmaskverts;j++) {
      k=maskverts[j];
      if(k<0) continue;
      for (m=0;m<mris->vertices[k].vnum;m++) {
        i = mris->vertices[k].v[m];
        if(!mask[i]) continue; // only include suprathreshold neighbors
        if(i<=k) continue; // only count neighbor pair once
        dist = mris->vertices[k].dist[m];
        if(dist<=0) continue;
        if(fwhm_amp_flag) {
          df_amp = (mris->vertices[i].val2 - mris->vertices[k].val2)/dist;
          df_amp_sum   += df_amp;
          df_amp_sumsq += df_amp*df_amp;
        }
        if(fwhm_phase_flag) {
          df_phase = (mris->vertices[i].val2bak - mris->vertices[k].val2bak)/dist;
          df_phase_sum   += df_phase;
          df_phase_sumsq += df_phase*df_phase;
        }
        if(fwhm_cx_flag) {
          df_real = (mris->vertices[i].val - mris->vertices[k].val)/dist;
          df_real_sum   += df_real;
          df_real_sumsq += df_real*df_real;

          df_imag = (mris->vertices[i].imag_val - mris->vertices[k].imag_val)/dist;
          df_imag_sum   += df_imag;
          df_imag_sumsq += df_imag*df_imag;
        }
        dfcount++;
      }
    }
    if(fwhm_amp_flag) {
      if(dfcount<=1) {
        fwhm=0;
      } else {
        dfvar = (df_amp_sumsq - (df_amp_sum*df_amp_sum)/dfcount)/(dfcount-1);
        varratio = 1.0 - 0.5*(dfvar/amp_var);
        if ((varratio<=0) || (dfvar<=0))
          fwhm=0;
        else
          fwhm = sqrt(-2.0*log(2.0)/log(varratio));
      }
      fprintf(fp,"ampFWHM      \t%-10.4f\n",fwhm);
      MsgPrintf("%s: amplitude FWHM\n",progname);
      MsgPrintf("    ampFWHM=%0.4f\n",fwhm);
      MsgPrintf("    amp neighbor var=%0.4f\n",dfvar);
      MsgPrintf("    amp overall var=%0.4f\n",amp_var);
      MsgPrintf("    amp sum diff=%0.4lf\n",df_amp_sum);
      MsgPrintf("    amp sumsq diff=%0.4lf\n",df_amp_sumsq);
      MsgPrintf("    num neighbors=%d\n",dfcount);
    }
    if(fwhm_phase_flag) {
      if(dfcount<=1) {
        fwhm=0;
      } else {
        dfvar = (df_phase_sumsq - (df_phase_sum*df_phase_sum)/dfcount)/(dfcount-1);
        varratio = 1.0 - 0.5*(dfvar/phase_var);
        if ((varratio<=0) || (dfvar<=0))
          fwhm=0;
        else
          fwhm = sqrt(-2.0*log(2.0)/log(varratio));
      }
      fprintf(fp,"phaseFWHM      \t%-10.4f\n",fwhm);
      MsgPrintf("%s: phase FWHM\n",progname);
      MsgPrintf("    phaseFWHM=%0.4f\n",fwhm);
      MsgPrintf("    phase neighbor var=%0.4f\n",dfvar);
      MsgPrintf("    phase overall var=%0.4f\n",phase_var);
      MsgPrintf("    phase sum diff=%0.4lf\n",df_phase_sum);
      MsgPrintf("    phase sumsq diff=%0.4lf\n",df_phase_sumsq);
      MsgPrintf("    num neighbors=%d\n",dfcount);
    }
    if(fwhm_cx_flag) {
      if(dfcount<=1) {
        fwhm=0;
      } else {
        dfvar = (df_real_sumsq - (df_real_sum*df_real_sum)/dfcount)/(dfcount-1) +
                (df_imag_sumsq - (df_imag_sum*df_imag_sum)/dfcount)/(dfcount-1);
        varratio = 1.0 - 0.5*(dfvar/cx_var);
        if ((varratio<=0) || (dfvar<=0))
          fwhm=0;
        else
          fwhm = sqrt(-2.0*log(2.0)/log(varratio));
      }
      fprintf(fp,"cxFWHM      \t%-10.4f\n",fwhm);
      MsgPrintf("%s: complex FWHM\n",progname);
      MsgPrintf("    cxFWHM=%0.4f\n",fwhm);
      MsgPrintf("    cx neighbor var=%0.4f\n",dfvar);
      MsgPrintf("    cx overall var=%0.4f\n",cx_var);
      MsgPrintf("    real sum diff=%0.4lf\n",df_real_sum);
      MsgPrintf("    real sumsq diff=%0.4lf\n",df_real_sumsq);
      MsgPrintf("    imag sum diff=%0.4lf\n",df_imag_sum);
      MsgPrintf("    imag sumsq diff=%0.4lf\n",df_imag_sumsq);
      MsgPrintf("    num neighbors=%d\n",dfcount);
    }
    fprintf(fp,"\n");
  }

  fclose(fp);

  delete [] mask;
  delete [] maskverts;

  MsgPrintf("%s: finished\n",progname);
  exit(0);
}

