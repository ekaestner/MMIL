fcfg = [];

fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
fcfg.fsr_nme = 'fsaverage';

fcfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
fcfg.prc_nme = '.aparc.annot';

fcfg.inc_reg = { 'Entorhinal' };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Alena_February';
fcfg.out_nme  = 'roi_for_alena';

ejk_roi_plot(fcfg);
