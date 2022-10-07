%%
out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Manuscript/figures/Figure1';

fcfg = [];

fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
fcfg.fsr_nme = 'fsaverage';

fcfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
fcfg.prc_nme = '.aparc.annot';

fcfg.inc_reg = { 'fusiform' 'lateralorbitofrontal' };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = out_dir;
fcfg.out_nme  = 'roi_illustration';

ejk_roi_plot(fcfg);

%%
out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Manuscript/figures/Figure2';

fcfg = [];

fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
fcfg.fsr_nme = 'fsaverage';

fcfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
fcfg.prc_nme = '.aparc.annot';

fcfg.inc_reg = { 'fusiform' };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = out_dir;
fcfg.out_nme  = 'roi_fusiform_illustration';

ejk_roi_plot(fcfg);

%%
out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/slh_atl_mem/Manuscript/figures/Figure2';

fcfg = [];

fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
fcfg.fsr_nme = 'fsaverage';

fcfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
fcfg.prc_nme = '.aparc.annot';

fcfg.inc_reg = { 'lateralorbitofrontal' };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = out_dir;
fcfg.out_nme  = 'roi_lateral_orbitofrontal_illustration';

ejk_roi_plot(fcfg);
