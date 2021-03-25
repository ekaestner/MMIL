clear; clc;

clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL';

sbj_nme = mmil_readtext( [clr_fld '/' 'subjects']);

sbj_col = distinguishable_colors(numel(sbj_nme)+1);
    sbj_col = flipud(sbj_col([1:7 9],:));

%%
fcfg = [];
fcfg.loc = '/home/ekaestne/PROJECTS/OUTPUT/SL/';

fcfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

%%
cfg = [];
cfg.typ = 'dir';
sbj_ele = mmil_find_file(cfg,[fcfg.loc '/' 'ele_idn']);

exs_ele{1} = mmil_readtext([fcfg.loc '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total_lhs_ecog']);
exs_ele{2} = mmil_readtext([fcfg.loc '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total_rhs_ecog']);

if ~isfield(fcfg,'sbj'); fcfg.sbj = { sbj_ele(:) }; end
if ~isfield(fcfg,'nme'); fcfg.nme = { 'overall' }; end

%%
lhs_ovr_ele_loc = mmil_readtext([fcfg.loc '/' 'manuscript' '/' 'figure1' '/' 'lhs_all_ele.csv']);
rhs_ovr_ele_loc = mmil_readtext([fcfg.loc '/' 'manuscript' '/' 'figure1' '/' 'rhs_all_ele.csv']);

for iS = 1:numel(sbj_nme)
    
    lhs_sel_ele_sbj{iS} =  lhs_ovr_ele_loc(string_find(lhs_ovr_ele_loc(:,1),sbj_nme{iS}),1);
    rhs_sel_ele_sbj{iS} =  rhs_ovr_ele_loc(string_find(rhs_ovr_ele_loc(:,1),sbj_nme{iS}),1);
    
    sel_ele_sbj{iS} = [ lhs_sel_ele_sbj{iS} ; rhs_sel_ele_sbj{iS} ];
    
    eff_lbl_sbj{iS} = sbj_nme{iS};
    
    col_sbj{iS} = sbj_col(iS,:);
    
end

ovr_ele = unique(sort([ lhs_ovr_ele_loc(:,1) ; rhs_ovr_ele_loc(:,1)]));
sel_ele_sbj{1}(1) = [];
sel_ele_sbj{3}(2) = [];

%% Electrodes on Brain
cfg = [];

cfg.clr_fld = clr_fld;

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' 'lh.pial']                           ['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' 'rh.pial']};
cfg.elec_text = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};
    
cfg.sel_lbl   = eff_lbl_sbj ;
cfg.sel_ele   = sel_ele_sbj;
cfg.all_ele   = ovr_ele;
cfg.col       = col_sbj;
cfg.nsl_col = { rgb('white')};
    
cfg.sve_img   = 'jpg';
cfg.sve_loc   = [ '/home/ekaestne/PROJECTS/OUTPUT/SL/revision/'];
cfg.sve_pre   = ['overall_electrode_location_subjects'];
cfg.sep_str   = [','];

cfg.rad    = 1.25;

cfg.ele_loc = { [fcfg.loc '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total' '_' 'split' '_' 'lhs' '_' 'ecog'] [fcfg.loc '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total' '_' 'split' '_' 'rhs' '_' 'ecog'] };
cfg.inc_reg = fcfg.inc_reg;

cfg.sml_vew = 1;

mmil_ieeg_sensor_location_plot_v4(cfg);

    
    
    
    
    
    
    
    
    
    
    
    
    