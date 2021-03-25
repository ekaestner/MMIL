% 
% cfg = [];
% 
% % Folders
% prj_dir = '/home/ekaestne/PROJECTS/';
% prj_nme = 'PhenotypeConnectome';
% 
% cfg.out_dir = [ prj_dir '/' 'OUTPUT' '/' prj_nme '/' 'tracts' ];
% cfg.out_nme = 'test';
% 
% cfg.trc_nme = {};

function ejk_tract_plot(cfg)

%%
avg_sbj = {'fc050c' , 'DTIPROC_fc050_fmri_120312_20120313.155403_1'};

und_flg = 1;

pln_nme = {'sag' 'sag' 'cor' 'cor' 'hor' 'hor'};
vew_nme = {'L'   'R'   'A'   'P'   'I'   'S'};

%% Find Fibers
fib_leg = mmil_readtext('/home/ekaestne/gitrep/MMIL/ANATOMICAL/PLOT/DTI_Fiber_Legend.csv');

fib_num = zeros(1,numel(cfg.trc_nme));
for iF = 1:numel(cfg.trc_nme)
    
    fib_num(iF) = fib_leg{find(strcmpi(fib_leg(:,2),cfg.trc_nme{iF})),1};
    
end

%% Make Plot
if ~isdir(cfg.out_dir); mkdir(cfg.out_dir); end

% Define Average Subject - Hard Coded
avg_sbj_loc = ['/space/md8/10/data/epiproj/EPDTI/Containers/' avg_sbj{1,2} '/AtlasTrack'];
fib_suf     = 'prob_countatlas_pthresh0.08_minflen12_path.grp';

mcd_hme     = '/home/mmilmcd/Desktop';

% Make Plot
for iV = 1:length(vew_nme)
    dti_render_AtlasTrack_ejk( avg_sbj_loc , ...
        ...
        'fiber_suffix'  , fib_suf                                                   , ...
        'fname_legend'  , ['DTI_Fiber_Legend.csv']                                  , ...
        'fname_tif'     , [ mcd_hme '/' 'avg_sbj' '_' cfg.out_nme '_' vew_nme{iV} '_' pln_nme{iV} '.tif']   , ...
        'fname_csh'     , [ cfg.out_dir '/' 'run_tracoview' '_' 'avg_sbj' '_' cfg.out_nme '_' num2str(iV) '.csh' ] , ...
        'fibers'        , fib_num                                                   , ...
        'plane'         , pln_nme{iV}                                               , ...
        'view'          , vew_nme{iV}                                               , ...
        'underlay_flag' , und_flg                                                   , ...
        'zflip'         , 1 )
end

end