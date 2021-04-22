
%% Tracts
planes = {'sag' 'sag' 'cor' 'cor' 'hor' 'hor'};
views = {'L' 'R' 'A' 'P' 'I' 'S'};
fibers = [123, 117, 118, 119, 120];
subjects = {'fc050c', 'DTIPROC_fc050_fmri_120312_20120313.155403_1'};
underlay_flag = [1];
image_alpha = [0.0];
for s = 1:size(subjects,1)
    for v = 1:length(views)
        ls
    dti_render_AtlasTrack(['/space/md8/10/data/epiproj/EPDTI/Containers/' subjects{s,2} '/AtlasTrack'],...
        'fiber_suffix', 'prob_countatlas_pthresh0.08_minflen12_path.grp',...
        'fname_legend',['DTI_Fiber_Legend.csv'],...
        'fname_tif',['~/tract_tiffs/' subjects{s,1} '_' views{v} '.tif'],...
        'fibers',fibers,...
        'plane', planes{v},...
        'view', views{v},...
        'underlay_flag', 1,...
        'zflip',1)
    end
end

%% ROI
fcfg = [];

fcfg.fsr_dir = '/home/ekaestne/PROJECTS/EXTERNAL/Misc'; 
fcfg.fsr_nme = 'fsaverage';

fcfg.roi_loc = '/home/ekaestner/gitrep/MMIL/EXTERNAL/freesurfer';
fcfg.prc_nme = '.aparc.annot';

fcfg.inc_reg = { 'fusiform' };

fcfg.sph = { 'lh' 'rh' };
fcfg.sph_vew = { 'lat' 'ven' 'med' };

fcfg.out_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/pst_opr/Naming/Manuscript/Figures/Figure2';
fcfg.out_nme  = 'fusiform_roi';

ejk_roi_plot(fcfg);