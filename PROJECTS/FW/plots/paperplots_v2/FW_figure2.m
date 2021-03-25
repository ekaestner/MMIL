clear; clc;

%% Total Electrodes
% Electrodes

% Percentages
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_rsp_600';

cfg.eff_clm = 1;
cfg.top_pct = 0.50;
cfg.eff_col = { [ rgb('dark blue')-[-.10 -.20 -.20] ; rgb('blue') ; rgb('neon blue')-[0 0.25 0] ] };
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.sve_loc = 'figure2';
mmil_include_plot(cfg)

% Colorbar
pcfg = [];
pcfg.col_map = [ rgb('medium grey')-0.075 ; rgb('dark blue')-[-.10 -.20 -.20] ; rgb('blue') ; rgb('neon blue')-[0 0.25 0] ];
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/NEW_FIGS';
pcfg.sve_pre = 'total';
mmil_color_bar(pcfg)

%% Specific Electrodes
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_rsp_600' };
fcfg.eff_clm = [ 1 ]; % 2
fcfg.eff_col = { rgb('blue') };
fcfg.eff_lbl = { 'Total' }; % 'Bigram Frequency'
fcfg.nsl_col = { rgb('white')};

fcfg.ovr_typ = 'hgp';
fcfg.ovr_nme = { 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1 2 3 ]; % 1

% Selective Electrodes
fcfg.sbj = { sbj_ele(:) };
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

% Overall Electrodes
cfg = [];
cfg.typ = 'dir';
sbj_ele = mmil_find_file(cfg,['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical' '/' 'ele_idn']);

exs_ele{1} = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical' '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total_lhs_ecog']);
exs_ele{2} = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical' '/' 'electrode_location_files' '/' 'total' '/' 'output' '/' 'total_rhs_ecog']);

fcfg.sbj{1} = intersect(fcfg.sbj{1},mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical' '/' 'subjects']));

ovr_ele = [];
for iS = 1:numel(fcfg.sbj{1})
    
    sbj_ele_hld = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical' '/' 'ele_idn' '/' fcfg.sbj{1}{iS} '/' fcfg.sbj{1}{iS} '_include.csv']);
    
    ovr_ele = [ovr_ele ; strcat(fcfg.sbj{1}{iS},'_',sbj_ele_hld(cell2mat(sbj_ele_hld(:,3))==1 & cell2mat(sbj_ele_hld(:,4))==1,2))];
    
end

ovr_ele = ovr_ele(ismember(ovr_ele(:,1),[exs_ele{1}(:,1) ; exs_ele{2}(:,1)]),:);

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.sel_ele;
cfg.all_ele   = ovr_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_col = fcfg.nsl_col;

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/'];
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_col;
mmil_loc_dot(pcfg)


