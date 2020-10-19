#!/bin/csh -f

setenv SUBJECTS_DIR /space/pogo1/6/dhagler/tp_roi_analysis/subjects

set srcsubj = ben2
set subjects = (clark2 nathan2 stewart2 christopher2 josh2 stephen2 tim2)
set labels = (frontal SPL)
set hemilist = (lh rh)

# reg to other subjects
foreach subject ($subjects)
  foreach label ($labels)
    foreach hemi ($hemilist)
      mri_label2label \
        --srclabel $SUBJECTS_DIR/$srcsubj/label/$hemi-$label.label \
        --trglabel $SUBJECTS_DIR/$subject/label/$hemi-$label.label \
        --regmethod surface \
        --hemi $hemi \
        --srcsubject $srcsubj \
        --trgsubject $subject
    end
  end
end  

