clear; clc;

%% Setup
clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL';

sbj_nme = mmil_readtext( [clr_fld '/' 'subjects']);

sbj_col = distinguishable_colors(numel(sbj_nme)+1);
    sbj_col = flipud(sbj_col([1:7 9],:));

%% Specific Electrodes
loc = clr_fld;

cfg = [];
cfg.typ = 'dir';
sbj_ele = sbj_nme; %mmil_find_file(cfg,[loc '/' 'ele_idn']);

%
fcfg.clr_fld = clr_fld;
fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_anv_1500' 'pap_anv_1500' };
fcfg.eff_clm = [ 1 2 ]; % 2
fcfg.eff_col = { rgb('purple') rgb('purple') };
fcfg.eff_lbl = { 'VisTotal'    'AudTotal' }; % 'Bigram Frequency'
fcfg.nsl_col = { rgb('white')};

% Selective Electrodes
fcfg.sbj = { sbj_ele(:) };
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

% Overall Electrodes
exs_ele{1} = mmil_readtext([ clr_fld '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total_lhs_ecog']);
exs_ele{2} = mmil_readtext([ clr_fld '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total_rhs_ecog']);

fcfg.sbj{1} = intersect(fcfg.sbj{1},mmil_readtext([clr_fld '/' 'subjects']));

ovr_ele = [];
for iS = 1:numel(fcfg.sbj{1})
    
    sbj_ele_hld = mmil_readtext([ clr_fld '/' 'ele_idn' '/' fcfg.sbj{1}{iS} '/' fcfg.sbj{1}{iS} '_include.csv']);
    
    ovr_ele = [ovr_ele ; strcat(fcfg.sbj{1}{iS},'_',sbj_ele_hld(cell2mat(sbj_ele_hld(:,3))==1 & cell2mat(sbj_ele_hld(:,4))==1,2))];
    
end

ovr_ele = ovr_ele(ismember(ovr_ele(:,1),[exs_ele{1}(:,1) ; exs_ele{2}(:,1)]),:);

%% Split into subjects
for iP = 1:numel(fcfg.sel_ele)
    for iS = 1:numel(sbj_nme)
        
        sel_ele_sbj{iP}{iS} =  fcfg.sel_ele{iP}(string_find(fcfg.sel_ele{iP},sbj_nme{iS}));
        
        eff_lbl_sbj{iP}{iS} = sbj_nme{iS};
        
        col_sbj{iP}{iS} = sbj_col(iS,:);
        
    end
end

%% Plot
sve_pre = { 'figure2_middle_pic_vis' 'figure2_middle_pic_aud' };

for iP = 1:numel(fcfg.sel_ele)
    
    cfg = [];
    
    cfg.clr_fld = clr_fld;
    
    cfg.hms = {'lhs' 'rhs'};
    cfg.hem = {'lhs' 'rhs'};
    
    cfg.pial_mat  = {['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' 'lh.pial']                           ['/home/ekaestne/PROJECTS/EXTERNAL/Misc/fsaverage/surf/' '/' 'rh.pial']};
    cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};
    
    cfg.sml_vew = 1;
    
    cfg.sel_ele   = sel_ele_sbj{iP};
    cfg.all_ele   = ovr_ele;
    
    cfg.sel_lbl = eff_lbl_sbj{iP};
    cfg.col     = col_sbj{iP};
    cfg.nsl_col = fcfg.nsl_col;
    
    cfg.sep_str         = ',';
    
    cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                    'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                    'Supramarginal' ...
                    'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                    'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
    
    cfg.sve_img   = 'eps';
    cfg.sve_loc   = [fcfg.clr_fld '/' 'revision' '/' ];
    cfg.sve_pre   = sve_pre{iP};
    cfg.sep_str   = [','];
    
    cfg.sve_img   = 'png';
    mmil_ieeg_sensor_location_plot_v4(cfg);
    
end

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_col;
mmil_loc_dot(pcfg)
