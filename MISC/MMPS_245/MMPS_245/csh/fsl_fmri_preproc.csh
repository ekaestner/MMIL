#!/bin/tcsh
####################################################################
# FSL processing of BOLD data
#
# usage: fsl_fmri_preproc.csh <boldvol> <maskvol> <t1vol> <outdir> <motionpar> <tr> <ndelete> <subjects_dir> <fsID> <smooth> <subjectID>
#
# Required Arguments:
#    boldvol: full path to BOLD vol (in procfolder)
#    maskvol: full path to aseg.mgz, wmparc.mgz or brainmask.mgz
#    t1vol: full path to nu.mgz
#    outdir: full path to output directory
#    motionpar: full path to file containing motion parameters
#    tr: TR in s
#    ndelete: number of vols to delete to allow for T1 eq
#    subjects_dir: full path to FS SUBJECTS_DIR
#    fsID: subject ID in SUBJECTS_DIR
#    smooth: FWHM for BOLD data
#	
# Optional Arguments:
#    subjectID: subject ID
#      {default = subj}
#
# Created:  11/17/11 by Lars T. Westlye  l.t.westlye@psykologi.uio.no
# Last Mod: 04/03/13 by Don Hagler
#
####################################################################

# parse input
if ($#argv == 0) then
  echo " "
  echo "USAGE: fsl_fmri_preproc.csh <boldvol> <maskvol> <t1vol> <outdir> <motionpar> <tr> <ndelete> <subjects_dir> <fsID> <smooth> <subjectID>"
  echo " "
  echo " Required Arguments:"
  echo "   boldvol	  full path to BOLD vol (in procfolder)"
  echo "   maskvol	  full path to aseg.mgz or wmparc.mgz"
  echo "   t1vol	  full path to nu.mgz"
  echo "   outdir	  full path to output directory"
  echo "   motionpar	  full path to file containing motion parameters"
  echo "   tr		  TR in s"
  echo "   ndelete	  number of vols to delete to allow for T1 eq"
  echo "   subjects_dir   full path to FS SUBJECTS_DIR"
  echo "   fsID 	  subject ID in SUBJECTS_DIR"
  echo "   smooth	  FWHM for BOLD data"
  echo " "
  echo " Optional Arguments:"
  echo "   subjectID      subject ID"
  echo "    {default = subj}"
  echo " "
  echo " "
  exit 1;
else
  set boldvol=$1	  
  set maskvol=$2	  
  set t1vol=$3	  
  set outdir=$4	 
  set motionpar=$5	  
  set tr=$6		  
  set ndelete=$7	 
  set subjects_dir=$8  
  set fsID=$9	  
  set smooth=$10	  
  set s=$11
endif

if($#argv < 11) then
  set s="subj"
else
  set outstem=$11
endif

if ($#argv < 10) then
  echo " "
  echo "ERROR: Did you specifiy all arguments? Please see usage."
  echo " "
endif

####################################################################

# setup FSL
#if ( `env | grep -c FSLDIR` == 0 ) then
source /usr/pubsw/bin/SetUpFSL.csh 4.1.9_RH5_64
#endif
echo "FSLDIR is ${FSLDIR}" 
setenv FSLOUTPUTTYPE NIFTI_GZ

####################################################################

# set vars you may want to change
set hipass="150"         # in seconds
set smoothdata="1"       # 0 turns off smoothing (SUSAN)
set resample2highres="0" # 1 turns on resampling of BOLD -> highres (1mm iso)
set bbrfs="1"            # 1 turns on FS BBR (FSL_BOLD -> FS space, no resampling)
set resample2fs="0"      # 1 turns on resampling of FSL_BOLD -> FS

# set other vars
set standard_head="${FSLDIR}/data/standard/MNI152_T1_2mm"
set regstandard="${standard_head}_brain"
set regstandard_dof="12"
set regstandard_nonlinear_warpres="10"
set brain_thresh="10"
set highres_head="${s}"  
set highres_file="${s}_brain"

# set fs SUBJECTS_DIR
setenv SUBJECTS_DIR $subjects_dir

# create outdir
if (! -ed ${outdir}) then
 mkdir ${outdir}
endif

####################################################################

# prepare nu and mask
cp $t1vol ${outdir}/t1.mgz
cp $maskvol ${outdir}/mask.mgz
cd $outdir
mri_binarize --i mask.mgz --min 0.2 --dilate 1 --o mask_dil_grot.mgz
mri_mask t1.mgz mask_dil_grot.mgz t1_masked_grot.mgz 
mri_convert -i t1_masked_grot.mgz -o ${s}_brain.nii.gz
mri_convert -i t1.mgz -o ${s}.nii.gz
fslreorient2std ${s} ${s}
fslreorient2std ${s}_brain ${s}_brain 
rm *mgz

# prepare fmri
mri_convert -i $boldvol -o ${outdir}/bold.nii.gz
cd ${outdir}
fslreorient2std bold bold
set nvols="`fslnvols bold`"
set suffix="`zeropad $nvols 3`"
immv bold bold_${suffix}  
set funcdata="bold_${suffix}"			

####################################################################

# start processing BOLD data
fslmaths $funcdata prefiltered_func_data -odt float
set total_volumes="$suffix"
echo "Total original volumes = $total_volumes"

# delete images
if ($ndelete > 0) then
 echo "Deleting $ndelete volume(s) - BE WARNED for future analysis!"
 set ps="Deleting $ndelete volume(s)"
 set total_volumes=`echo "$total_volumes - $ndelete " | bc -l`
 fslroi $funcdata prefiltered_func_data $ndelete $total_volumes
 set funcdata="prefiltered_func_data"
else
 echo "Using all volumes"
 set ps="Using all volumes"
endif

# choose target image and cp to example_func 
set target_vol_number="1"
fslroi $funcdata example_func $target_vol_number 1

# run BET on BOLD
echo "non-brain removal using BET"
set ps="$ps; brain extraction using BET"
fslmaths $funcdata -Tmean mean_func
bet2 mean_func mask -f 0.3 -n -m
immv mask_mask mask
fslmaths $funcdata -mas mask prefiltered_func_data_bet
set funcdata="prefiltered_func_data_bet"

# set min intensity threshold and mask
set int_2_98=`fslstats $funcdata -p 2 -p 98`
set int2=`echo "${int_2_98}" | awk '{print $1}'` 
set int98=`echo "${int_2_98}" | awk '{print $2}'`
set intensity_threshold=`echo "${int2} + ( ${brain_thresh} * ( ${int98} - ${int2} )  / 100.0 )"  | bc -l`
set funcdata_unmasked=${funcdata}
fslmaths $funcdata -thr $intensity_threshold -Tmin -bin mask -odt char
set median_intensity=`fslstats $funcdata_unmasked -k mask -p 50`
fslmaths mask -dilF mask
fslmaths $funcdata_unmasked -mas mask prefiltered_func_data_thresh
set funcdata="prefiltered_func_data_thresh"

# spatial filtering
if (${smoothdata}>0) then
 echo "spatial smoothing using a Gaussian kernel of FWHM ${smooth} mm"
 set smoothsigma=`echo "${smooth} / 2.355 " | bc -l`
 set susan_int=`echo "( ${median_intensity} - ${int2} ) * 0.75" | bc -l` 
 fslmaths $funcdata -Tmean mean_func
 susan $funcdata $susan_int $smoothsigma 3 1 1 mean_func $susan_int prefiltered_func_data_smooth
 fslmaths prefiltered_func_data_smooth -mas mask prefiltered_func_data_smooth
 set funcdata="prefiltered_func_data_smooth"
 set ps="$ps; spatial smoothing using a Gaussian kernel of FWHM ${smooth} mm"
endif

# intensity normalization
set normmean="10000"
set scaling=`echo "${normmean} / ${median_intensity}" | bc -l`
fslmaths $funcdata -mul $scaling prefiltered_func_data_intnorm
echo "grand-mean intensity normalisation by a single multiplicative factor"
set ps="$ps; grand-mean intensity normalisation of the entire 4D dataset"
set funcdata="prefiltered_func_data_intnorm"

# Temporal filtering
set lp_sigma_vol="-1"
set hp_sigma_sec=`echo "${hipass} / 2.0 " | bc -l`
set hp_sigma_vol=`echo "${hp_sigma_sec} / ${tr}" | bc -l`
echo "hp temporal filtering (Gaussian-weighted least-squares straight line fitting, with sigma=${hp_sigma_sec} s"
set ps="${ps}; highpass temporal filtering (Gaussian-weighted least-squares straight line fitting, with sigma=${hp_sigma_sec} s"
fslmaths $funcdata -bptf $hp_sigma_vol $lp_sigma_vol prefiltered_func_data_tempfilt
set funcdata="prefiltered_func_data_tempfilt"

# set threshold
fslmaths $funcdata filtered_func_data
set ps="$ps."
set absbrainthresh=`expr "${brain_thresh} * $normmean / 100.0 " | bc -l`

#set TR in header of filtered_func_data
set IMTR="`fslval filtered_func_data pixdim4`"
fslhd -x filtered_func_data | sed 's/  dt = .*/  dt = '${tr}'/g' > tmpHeader
fslcreatehd tmpHeader filtered_func_data
rm tmpHeader

####################################################################

# prepare registration
echo "various registration steps"
mkdir -p ${outdir}/reg
cd reg
imcp ../example_func example_func

# prepare files
fslmaths ../$highres_file highres
fslmaths ../$highres_head highres_head
fslmaths $regstandard standard
fslmaths $standard_head standard_head
fslmaths $regstandard -bin -dilF -dilF standard_mask -odt char

# highres -> standard 
set in="highres"
set ref="standard"
set out="${in}2${ref}"
set interp="trilinear"
flirt -ref $ref -in $in -out $out -omat ${out}.mat -interp $interp -cost corratio -dof $regstandard_dof
immv highres2standard highres2standard_linear
set conf="T1_2_MNI152_2mm"
fnirt --in=highres_head --aff=highres2standard.mat --cout=highres2standard_warp --iout=highres2standard --jout=highres2standard_jac --config=$conf --ref=standard_head --refmask=standard_mask --warpres=$regstandard_nonlinear_warpres,$regstandard_nonlinear_warpres,$regstandard_nonlinear_warpres
convert_xfm -inverse -omat ${ref}2${in}.mat ${out}.mat
slicer $out $ref -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png ${out}1.png 
slicer $ref $out -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png ${out}2.png 
pngappend ${out}1.png - ${out}2.png ${out}.png
rm -f sl*.png

# FLIRT-BBR (incl FAST segmentation)
cd $outdir
if (! -ed reg_bbr) then
  mkdir reg_bbr
endif
cp reg/highres_head.nii.gz reg_bbr/${s}.nii.gz
cp reg/highres.nii.gz reg_bbr/${s}_brain.nii.gz
echo "run epi_reg_ltw (incl FAST segm)"
epi_reg_ltw example_func \
     reg_bbr/${s} \
     reg_bbr/example_func2highres

# apply bbr xform
echo "apply bbr xform"
flirt_bbr_ltw -in example_func \
     -ref reg_bbr/${s}_brain \
     -dof 6 \
     -applyxfm -init reg_bbr/example_func2highres.mat \
     -out reg_bbr/example_func2highres \

# copy new xform to reg 
cp -f reg_bbr/example_func* reg/

# convert/inverse xforms and prepare figs 
echo "run slicer"
cd reg
convert_xfm -omat example_func2standard.mat -concat highres2standard.mat example_func2highres.mat
convert_xfm -inverse -omat highres2example_func.mat example_func2highres.mat
slicer example_func2highres highres -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2highres1.png 
slicer highres example_func2highres -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2highres2.png 
pngappend example_func2highres1.png - example_func2highres2.png example_func2highres.png
rm -f sl*.png

# applywarp and prepare figs
echo "warp and resample example_func"
cd $outdir/reg
applywarp \
  --ref=standard \
  --in=example_func \
  --out=example_func2standard \
  --warp=highres2standard_warp \
  --premat=example_func2highres.mat

convert_xfm -inverse -omat standard2example_func.mat example_func2standard.mat
slicer example_func2standard standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2standard1.png 
slicer standard example_func2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 
pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2standard2.png 
pngappend example_func2standard1.png - example_func2standard2.png example_func2standard.png
rm -f sl*.png

# prepare html document for QC
cd ${outdir}/reg
echo "<HTML><TITLE>registration_summary</TITLE><BODY BGCOLOR=#ffffff" > index.html
echo "<a href="highres2standard.png"><img src="highres2standard.png" WIDTH=750 > ${s}_highres2standard</a><br>" >> index.html
echo "<a href="example_func2highres.png"><img src="example_func2highres.png" WIDTH=750 > ${s}_example_func2highres</a><br>" >> index.html
echo "<a href="example_func2standard.png"><img src="example_func2standard.png" WIDTH=750 > ${s}_example_func2standard</a><br>" >> index.html
echo "</BODY></HTML>" >> index.html

# prepare motion parameter files
mkdir ${outdir}/motioncorr 
cd ${outdir}/motioncorr
cp ${motionpar} motionparams_all.dat
cat motionparams_all.dat | awk '{print $2 " " $3 " " $4 " " $5 " " $6 " " $7}' | tail -n ${total_volumes} > ${outdir}/motioncorr/motionparams_pruned.dat 
fsl_tsplot -i motionparams_all.dat -t '"roll" - rotation about the I-S axis' -u 1 --start=2 --finish=2 -a degrees -w 640 -h 144 -o roll.png
fsl_tsplot -i motionparams_all.dat -t '"pitch" - rotation about the R-L axis' -u 1 --start=3 --finish=3 -a degrees -w 640 -h 144 -o pitch.png
fsl_tsplot -i motionparams_all.dat -t '"yaw" - rotation about the A-P axis' -u 1 --start=4 --finish=4 -a degrees -w 640 -h 144 -o yaw.png
fsl_tsplot -i motionparams_all.dat -t 'dS - displacement in the Superior direction' -u 1 --start=5 --finish=5 -a mm -w 640 -h 144 -o dS.png
fsl_tsplot -i motionparams_all.dat -t 'dL - displacement in the Left direction' -u 1 --start=6 --finish=6 -a mm -w 640 -h 144 -o dL.png
fsl_tsplot -i motionparams_all.dat -t 'dP - displacement in the Posterior direction' -u 1 --start=7 --finish=7 -a mm -w 640 -h 144 -o dP.png
fsl_tsplot -i motionparams_all.dat -t 'rmsold - RMS diff between input and base vol' -u 1 --start=8 --finish=8 -a rms -w 640 -h 144 -o rmsold.png
fsl_tsplot -i motionparams_all.dat -t 'rmsnew - RMS diff between output and base vol' -u 1 --start=9 --finish=9 -a rms -w 640 -h 144 -o rmsnew.png

#cat motionparams_all.dat | awk '{print $8}' > rmsold.txt
#cat motionparams_all.dat | awk '{print $9}' > rmsnew.txt
#cat rmsold.txt | awk '{ sum+=$1} END {print sum}'

# calculate summary motion measures
set motion1=`cat motionparams_all.dat | awk '{sum+=$8} END { print "mean(RMSold) = ",sum/NR}'`
set motion2=`cat motionparams_all.dat | awk '{sum+=$9} END { print "mean(RMSnew) = ",sum/NR}'`
set motion3=`cat motionparams_all.dat | awk '{sum1+=$8} {sum2+=$9} END { print "mean(RMSold)-mean(RMSnew) = ",sum1/NR - sum2/NR}'`

# generate html document for QC
echo "<HTML><TITLE>motioncorr_summary</TITLE><BODY BGCOLOR=#ffffff" > index.html
echo "<a href=roll.png""><img src="roll.png" WIDTH=750 > ${s}_motioncorr_roll</a><br>" >> index.html
echo "<a href=pitch.png""><img src="pitch.png" WIDTH=750 > ${s}_motioncorr_pitch</a><br>" >> index.html
echo "<a href=yaw.png""><img src="yaw.png" WIDTH=750 > ${s}_motioncorr_yaw</a><br>" >> index.html
echo "<a href=dS.png""><img src="dS.png" WIDTH=750 > ${s}_motioncorr_dS</a><br>" >> index.html
echo "<a href=dL.png""><img src="dL.png" WIDTH=750 > ${s}_motioncorr_dL</a><br>" >> index.html
echo "<a href=dP.png""><img src="dP.png" WIDTH=750 > ${s}_motioncorr_dP</a><br>" >> index.html
echo "<a href=rmsold.png""><img src="rmsold.png" WIDTH=750 > ${s}_motioncorr_rmsold</a><br>" >> index.html
echo "<a href=rmsnew.png""><img src="rmsnew.png" WIDTH=750 > ${s}_motioncorr_rmsnew</a><br>" >> index.html
echo "</BODY></HTML>" >> index.html


####################################################################

# resample 4D
echo "warp and resample 4D"
cd $outdir
mkdir reg_standard
applywarp --ref=reg/standard --in=filtered_func_data \
  --out=reg_standard/filtered_func_data \
  --warp=reg/highres2standard_warp \
  --premat=reg/example_func2highres.mat \
  --interp=trilinear

# clean up
foreach file (wmgrad wmsmooth wmsum)
  mv ${file}.nii.gz  reg_bbr/${file}.nii.gz
end	

# resample bold -> highres 
if (${resample2highres}>0) then
 echo "resampling 4D to highres space"
 cd $outdir
 flirt -in filtered_func_data \
    -applyxfm -init reg/example_func2highres.mat \
    -ref reg/highres -o reg/filtered_func_data_highres
endif

# FS bbreg to obtain bold->fs xform (should be ~identical to FSL bbr xform) 
if (${bbrfs}>0) then
 echo "BBR FSL->FS space"
 cd $outdir
 bbregister --s ${fsID} \
    --mov ${outdir}/example_func.nii.gz \
    --reg ${outdir}/reg/example_func2freesurfer_bbr.dat \
    --init-fsl --bold --o ${outdir}/reg/example_func2freesurfer_bbr.nii.gz
endif
   
# resample volume (omit to save space?)
if (${resample2fs}>0) then
  echo "resampling BOLD -> FS space" 
  cd $outdir
  mri_vol2vol --mov filtered_func_data.nii.gz \
     --targ ${SUBJECTS_DIR}/${fsID}/mri/nu.mgz \
     --o ${outdir}/reg/filtered_func2freesurfer.nii.gz \
     --reg ${outdir}/reg/example_func2freesurfer_bbr.dat
endif

# clean up
cd $outdir
imrm bold* example_func *thresh* *smooth_usan_size*
rm reg/*1.png reg/*2.png

####################################################################

# prepare summary output 
echo "${s} done `date`"
echo " "
echo $ps
echo "Summary of processing steps:" > $outdir/summary.txt
echo $ps >> $outdir/summary.txt
echo " " >> $outdir/summary.txt
echo " " >> $outdir/summary.txt
echo "Summary motion estimates: " >> $outdir/summary.txt
echo "$motion1" >> $outdir/summary.txt
echo "$motion2" >> $outdir/summary.txt
echo "$motion3" >> $outdir/summary.txt
echo " " >> $outdir/summary.txt
echo " " >> $outdir/summary.txt
echo "check registration manually:" >> $outdir/summary.txt
echo "firefox ${outdir}/reg/index.html" >> $outdir/summary.txt
echo " "
echo "check registration manually:" 
echo "firefox ${outdir}/reg/index.html"
echo " "
