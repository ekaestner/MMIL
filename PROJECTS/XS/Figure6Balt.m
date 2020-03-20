clear; clc;

clr_loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/';

ele_loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/electrode_location_files/';

%%
krs_ele_ovr = mmil_readtext('/home/ekaestne/gitrep/MMIL/PROJECTS/XS/Figure1.csv');
    krs_ele_ovr_nme = strcat( krs_ele_ovr(:,1) , '_' , krs_ele_ovr(:,2) );

krs_ele_sel = mmil_readtext('/home/ekaestne/gitrep/MMIL/PROJECTS/XS/Figure6Balt.csv');
    krs_ele_nme = strcat( krs_ele_sel(:,1) , '_' , krs_ele_sel(:,2) );
    
sbj_nme = unique( krs_ele_ovr(:,1) );

%%
tot_ele = [];
tot_ele_lhs = [];
tot_ele_rhs = [];

for iS = 1:numel(sbj_nme)
    
    try
        frs_ele_loc_lhs = mmil_readtext( [ ele_loc '/' sbj_nme{iS} '_' 'XS' '/' 'output' '/' sbj_nme{iS} '_' 'XS' '_' 'lhs_srf_fsaverage.csv' ] );
    catch
        frs_ele_loc_lhs = [];
    end
    
    try
        frs_ele_loc_rhs = mmil_readtext( [ ele_loc '/' sbj_nme{iS} '_' 'XS' '/' 'output' '/' sbj_nme{iS} '_' 'XS' '_' 'rhs_srf_fsaverage.csv' ] );
    catch
        frs_ele_loc_rhs = [];
    end
    
    % lhs
    if ~isempty(frs_ele_loc_lhs)
        tot_ele_lhs = [ tot_ele_lhs ;  [ strcat( sbj_nme{iS} , '_' , frs_ele_loc_lhs(:,1)) frs_ele_loc_lhs(:,2:end) ] ];
    end
    
    % rhs
    if ~isempty(frs_ele_loc_rhs)
        tot_ele_rhs = [ tot_ele_rhs ;  [ strcat( sbj_nme{iS} , '_' , frs_ele_loc_rhs(:,1)) frs_ele_loc_rhs(:,2:end) ] ];
    end
    
    % total
    tot_ele_hld = [ frs_ele_loc_lhs ; frs_ele_loc_rhs ];
    
    [ ~ , srt_ind ] = sort(tot_ele_hld(:,1));
    tot_ele_hld = tot_ele_hld( srt_ind , : );
    
    tot_ele = [ tot_ele ;  [ strcat( sbj_nme{iS} , '_' , tot_ele_hld(:,1)) tot_ele_hld(:,2:end) ] ];
    
end

%%
% lhs total
cell2csv( [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_lhs_ecog'] , tot_ele_lhs );

% rhs total
cell2csv( [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_rhs_ecog'] , tot_ele_rhs );

% Over total
[ ~ , ~ , krs_ind ] = intersect( tot_ele(:,1) , krs_ele_nme );

size(krs_ind);

krs_mss_nme = setxor( krs_ele_nme , krs_ele_nme(krs_ind,1) );
    size(krs_mss_nme) % Missing 137 electrodes
    krs_mss_sbj_nme = cellfun( @(x) x(1:5) , krs_mss_nme , 'uni' , 0 );
tot_mss_nme = setxor( tot_ele(:,1) , krs_ele_nme(krs_ind,1) );

krs_ovr_nme = krs_ele_nme( krs_ind , 1 );

% Select Total
tot_col = unique( krs_ele_sel(:,end) );

for iCL = 1:numel(tot_col)
    
    krs_chs_nme_hld{iCL} = strcat( krs_ele_sel( strcmpi( krs_ele_sel(:,end) , tot_col{iCL} ) , 1) , '_' , krs_ele_sel( strcmpi( krs_ele_sel(:,end) , tot_col{iCL} ) , 2) );
        krs_chs_nme{iCL} = intersect( krs_chs_nme_hld{iCL} , krs_ovr_nme );
    krs_chs_lbl{iCL} = tot_col{iCL};
    krs_chs_col{iCL} = rgb( tot_col{iCL} );
    
end

%%
% Brain
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[ ele_loc '/' 'fsaverage' '/'  'surf' '/' 'lh.pial']                       [ ele_loc '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_lhs_ecog']  [ ele_loc '/' 'total' '/'  'electrode_location' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = krs_chs_nme;
cfg.all_ele   = krs_ovr_nme;

cfg.sel_lbl = krs_chs_lbl;
cfg.col     = krs_chs_col;
cfg.nsl_col = { rgb('light grey') };

cfg.sep_str         = ',';

cfg.sve_loc   = '/home/ekaestne/gitrep/MMIL/PROJECTS/XS';
cfg.sve_pre   = 'figure6Balt';
cfg.sve_img   = 'png';

mmil_ieeg_sensor_location_plot_v4(cfg);

% Balls
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_col;
mmil_loc_dot(pcfg)



