clear; clc;

clr_loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/';

ele_loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/XS/epoch_data/clerical/electrode_location_files/';

%%
krs_ele = mmil_readtext('/home/ekaestne/gitrep/MMIL/PROJECTS/XS/Figure1.csv');
    krs_ele_nme = strcat( krs_ele(:,1) , '_' , krs_ele(:,2) );

sbj_nme = unique( krs_ele(:,1) );

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
% cell2csv( [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_lhs_ecog'] , tot_ele_lhs );
tot_ele_lhs = mmil_readtext( [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_lhs_ecog']);

% rhs total
% cell2csv( [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_rhs_ecog'] , tot_ele_rhs );

tot_ele = [ tot_ele_lhs ; tot_ele_rhs];

% total
[ ~ , ~ , krs_ind ] = intersect( tot_ele(:,1) , krs_ele_nme );

krs_mss_nme = setxor( krs_ele_nme , krs_ele_nme(krs_ind,1) );

    krs_mss_sbj_nme = cellfun( @(x) x(1:5) , krs_mss_nme , 'uni' , 0 );
tot_mss_nme = setxor( tot_ele(:,1) , krs_ele_nme(krs_ind,1) );

krs_chs_nme = {krs_ele_nme( krs_ind , 1 )};
krs_ovr_nme = krs_ele_nme( krs_ind , 1 );

%% Report
tot_exp_ele = size( krs_ele_nme , 1);
tot_prs_ele = size( krs_ind , 1);
tot_mss_ele = size( krs_mss_nme , 1);

fprintf( [ 'FIGURE1\n ' ...
           'TOTAL EXPECTED: %i\n ' ...
           'TOTAL PRESENT: %i\n ' ...
           'TOTAL MISSING: %i\n\n\n '], ...
           tot_exp_ele , tot_prs_ele , tot_mss_ele )

%%
% Brain
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[ ele_loc '/' 'fsaverage' '/'  'surf' '/' 'lh.pial']                       [ ele_loc '/' 'fsaverage' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_lhs_ecog']  [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_rhs_ecog']};

cfg.sml_vew = 0.75;

cfg.sel_ele   = krs_chs_nme;
cfg.all_ele   = krs_ovr_nme;

cfg.sel_lbl = { 'Total' };
cfg.col     = { rgb('bluish grey') };
cfg.nsl_col = { rgb('bluish grey') };

cfg.ele_rad = 1;

cfg.sep_str         = ',';

cfg.sve_loc   = '/home/ekaestne/gitrep/MMIL/PROJECTS/XS';
cfg.sve_pre   = 'figure1';
cfg.sve_img   = 'png';

mmil_ieeg_sensor_location_plot_v4(cfg);

%% With parcellations
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[ ele_loc '/' 'fsaverage' '/'  'surf' '/' 'lh.pial']                       [ ele_loc '/' 'fsaverage' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_lhs_ecog']  [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_rhs_ecog']};
cfg.aparc     = {[ ele_loc '/' 'fsaverage' '/'  'label' '/' 'lh.aparc.annot']               [ ele_loc '/' 'fsaverage' '/'  'label' '/' 'rh.aparc.annot']};

cfg.sml_vew = 0.75;

cfg.sel_ele   = krs_chs_nme;
cfg.all_ele   = krs_ovr_nme;

cfg.sel_lbl = { 'Total' };
cfg.col     = { rgb('bluish grey') };
cfg.nsl_col = { rgb('bluish grey') };

cfg.ele_rad = 1;

cfg.sep_str         = ',';

cfg.sve_loc   = '/home/ekaestne/gitrep/MMIL/PROJECTS/XS';
cfg.sve_pre   = 'figure1_withparcellation';
cfg.sve_img   = 'png';

mmil_ieeg_sensor_location_plot_v4(cfg);

%% With parcellations & White
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[ ele_loc '/' 'fsaverage' '/'  'surf' '/' 'lh.pial']                       [ ele_loc '/' 'fsaverage' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_lhs_ecog']  [ ele_loc '/' 'fsaverage' '/'  'electrode_location' '/' 'total_rhs_ecog']};
cfg.aparc     = {[ ele_loc '/' 'fsaverage' '/'  'label' '/' 'lh.aparc.annot']               [ ele_loc '/' 'fsaverage' '/'  'label' '/' 'rh.aparc.annot']};

cfg.sml_vew = 0.75;

cfg.sel_ele   = krs_chs_nme;
cfg.all_ele   = krs_ovr_nme;

cfg.sel_lbl = { 'Total' };
cfg.col     = { rgb('white') };
cfg.nsl_col = { rgb('bluish grey') };

cfg.ele_rad = 1;

cfg.sep_str         = ',';

cfg.sve_loc   = '/home/ekaestne/gitrep/MMIL/PROJECTS/XS';
cfg.sve_pre   = 'figure1_withparcellation_whitespheres';
cfg.sve_img   = 'png';

mmil_ieeg_sensor_location_plot_v4(cfg);










