#!/bin/bash

# Usage: makempeg.sh (container name) (output type)
# 
# Valid output types are:
# "surf, aseg, dv, petreg"
#
# Note: SUBJECTS_DIR variable must already be set
#

# check for correct number of input arguments
if [ $# -ne 2 ]; then 
   echo "Usage: $0 (container name) (output type)"
   echo "Valid output types:"
   echo "surf, aseg, dv, petreg"
   exit 1
else 
   subj=$1
   output=$2
fi

mainvol="brain.mgz"
startdir=$PWD
subjdir="$SUBJECTS_DIR/$subj"

if [ -d $subjdir ]; then # confirm specified main volume and subject dir is valid
   cd $SUBJECTS_DIR
   mkdir -p "$subj/mpeg" # maybe add error checking here in case the directory is not made
   case "$output" in 
     "surf" )
     tkmedit $subj $mainvol lh.white -aux wm.mgz -aux-surface rh.white -tcl /home/mmilrec/bin/makempeg_$output.tcl
     ;;
     
     "aseg" )
     tkmedit $subj $mainvol -segmentation mri/aseg.mgz $FREESURFER_HOME/tkmeditColorsCMA -tcl /home/mmilrec/bin/makempeg_$output.tcl
     ;;
     
     "dv" )
     echo "todo"
     exit 1
     cd $subj
     # get proper overlay settings set in tcl file
     tkmedit -f mri/affinereg_to_baseline/baseline_LIA_TP2_GlobalIntensityNorm.mgz -overlay mri/nonlinreg_to_baseline/dv.mgz -tcl /home/mmilrec/bin/makempeg_$output.tcl 
     ;;
     
     "petreg" )
     echo "todo"
     exit 1
     # change so loops through all PET_reg*.txt (e.g. GMR and Patlak for some SAX_OCD)
     cd $subj
     tkmedit -f `cat PET_reg.txt | awk '{print $3}'` -aux PET_reg.mgh -overlay PET_reg.mgh -overlay-reg-identity -tcl /home/mmilrec/bin/makempeg_$output.tcl &
     ;;
   esac
else
   echo "Error: Subject directory $subjdir not found"
   exit 1
fi

mpegoutdir="$subjdir/mpeg"
cd $mpegoutdir

# find all blank image files
echo "Finding empty slices to omit from video..."
allempty=$(find . -type f -print0 | xargs -0 md5sum | sort | uniq -w32 -d --all-repeated=separate | cut -c35-)
numempty=$(md5sum * | sort | uniq -w32 -d --all-repeated=separate | wc -l)

# loop thru each empty file and delete
for (( i = 1 ; i <= $numempty ; i++ )); do
   emptyfile=$(echo $allempty | cut -d" " -f$i)
   echo "Removing blank image $emptyfile"
   rm -f $emptyfile
done

#extract unique subject id and scan date to append to mpeg file name
subjid=$(echo $subj | sed 's/FREESURFERRECON_//' | cut -d. -f1)

#run Imagemagick to convert to mpeg
echo "Converting TIFF into MPEG for horizontal slices..."
convert -quality 100 "$output"_HOR*.tif "$subjid"_"$output"_HOR.mpeg >& /dev/null
echo "Converting TIFF into MPEG for coronal slices..."
convert -quality 100 "$output"_COR*.tif "$subjid"_"$output"_COR.mpeg >& /dev/null
echo "Converting TIFF into MPEG for sagittal slices..."
convert -quality 100 "$output"_SAG*.tif "$subjid"_"$output"_SAG.mpeg >& /dev/null
echo "Conversion complete."

#rm -f *.tif

cd "$startdir"


