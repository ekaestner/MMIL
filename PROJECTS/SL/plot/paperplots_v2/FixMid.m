%% Figure 3
fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_lng_950' 'pap_lng_950' };
fcfg.eff_clm = [ 1 2 ]; % 2
fcfg.eff_col = { rgb('red')-[0.3 0 0] rgb('bright blue')-[0 0.13 0.3] }; % rgb('blue')
fcfg.eff_lbl = { 'Text Selective' 'Voice Selective' }; % 'Bigram Frequency'
fcfg.nsl_col = { rgb('dark purple')};

fcfg.ovr_typ = 'hgp';
fcfg.ovr_nme = { 'pap_anv_1500' 'pap_anv_1500' 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1 2 3 ]; % 1

% %%%%%%%%%%%%%%%%%%%%%%%%
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

fcfg.sel_ele{3} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{2});
fcfg.sel_ele{1} = setxor(fcfg.sel_ele{1},fcfg.sel_ele{3});
fcfg.sel_ele{2} = setxor(fcfg.sel_ele{2},fcfg.sel_ele{3});

fcfg.eff_col{3} = rgb('orange'); %rgb('bright violet')+[0.1 0.03 0];

fcfg.eff_lbl{3} = 'Bi-Modal Selective';
for iEF = 1:numel(fcfg.ovr_clm)
    ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme{iEF} '_plt']);
    
    fcfg.all_ele{iEF} = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm(iEF)+2)))+1,2);
end
fcfg.all_ele = unique(cat(1,fcfg.all_ele{:}));


% %%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};


cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.sel_ele;
cfg.all_ele   = fcfg.all_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_col = fcfg.nsl_col;

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal'  };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'NEWFIGS' '/'];

cfg.sep_str   = [','];

cfg.sve_img   = 'png';

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_lat']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_lat']};
cfg.sve_pre   = ['middle_pic_figure3_lat'];
mmil_ieeg_sensor_location_plot_v4(cfg);

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_ven']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_ven']};
cfg.sve_pre   = ['middle_pic_figure3_ven'];
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Figure 3_5
fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_phn_950' 'pap_phn_950' };
fcfg.eff_clm = [ 1 2 ]; % 2 
fcfg.eff_col = { rgb('bright red') rgb('bright blue') }; 
fcfg.eff_lbl = { 'Letter-Selective' 'Phoneme-Selective' }; 
fcfg.nsl_col = { rgb('purple') };

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = { 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 4 ]; % 1

% Middle Portion
% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end
fcfg.sel_ele{2}([20 63]) = [];

% Overall Electrodes
for iEF = 1:numel(fcfg.ovr_clm)
    ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme{iEF} '_plt']);
    
    fcfg.all_ele{iEF} = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm(iEF)+2)))+1,2);
end
fcfg.all_ele = unique(cat(1,fcfg.all_ele{:}));

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.sel_ele;
cfg.all_ele   = fcfg.all_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_col = fcfg.nsl_col;

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'NEWFIGS' '/'];
cfg.sve_pre   = ['middle_pic_figure3_5'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};
cfg.sve_pre   = ['middle_pic_figure3_5'];
mmil_ieeg_sensor_location_plot_v4(cfg);

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_lat']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_lat']};
cfg.sve_pre   = ['middle_pic_figure3_5_lat'];
mmil_ieeg_sensor_location_plot_v4(cfg);

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_ven']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_ven']};
cfg.sve_pre   = ['middle_pic_figure3_5_ven'];
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Figure 4
fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_mtc_1450'       'pap_lng_950'    'pap_lng_950'      }; %{ 'pap_mtc_1450' 'pap_mtc_1450' };
fcfg.eff_clm = [ 1                    1                2                  ];             %[ 1 2 ]; % 2 
fcfg.eff_col = { rgb('bright yellow') rgb('red')       rgb('blue')        };   %{ rgb('bright yellow') rgb('orangey yellow')-[0.35 0.2 0.08] }; % rgb('blue')
fcfg.eff_lbl = { 'StimulusMisMatch'   'VisualLanguage' 'AuditoryLanguage' }; %{ 'StimulusMisMatch' 'PostStimulusMisMatch' }; % 'Bigram Frequency'
fcfg.nsl_col = {rgb('purple')};

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = { 'pap_anv_1500' 'pap_anv_1500' 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1 2 3 ]; % 1

% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

fcfg.sel_ele{2} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{2});
fcfg.sel_ele{3} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{3}); [~,~,rmv_ind] = intersect(fcfg.sel_ele{2},fcfg.sel_ele{3}); fcfg.sel_ele{3}(rmv_ind) = [];

% Overall Electrodes
for iEF = 1:numel(fcfg.ovr_clm)
    ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme{iEF} '_plt']);
    
    fcfg.all_ele{iEF} = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm(iEF)+2)))+1,2);
end
fcfg.all_ele = unique(cat(1,fcfg.all_ele{:}));

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};

cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.sel_ele;
cfg.all_ele   = fcfg.all_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_col = fcfg.nsl_col;
cfg.rad     = 3.25;

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'NEWFIGS' '/'];

cfg.sep_str   = [','];

cfg.sve_img   = 'png';

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};
cfg.sve_pre   = ['middle_pic_Figure4'];
mmil_ieeg_sensor_location_plot_v4(cfg);

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_lat']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_lat']};
cfg.sve_pre   = ['middle_pic_Figure4_lat'];
mmil_ieeg_sensor_location_plot_v4(cfg);

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_ven']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_ven']};
cfg.sve_pre   = ['middle_pic_Figure4_ven'];
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Figure 5
fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_lng_950' 'pap_con_950' };
fcfg.eff_clm = [ 1 2 ]; % 2 
fcfg.eff_col = { rgb('red') rgb('cyan') }; % rgb('blue')
fcfg.eff_lbl = { 'Text Selective' 'Noise-Vocoded Selective' }; % 'Bigram Frequency'
fcfg.nsl_col = { rgb('purple') };

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = { 'pap_anv_1500' 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1 3 ]; % 1

% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

fcfg.sel_ele{1} = unique( [ intersect(fcfg.sel_ele{1},fcfg.sel_ele{2}) ] );

% Overall Electrodes
for iEF = 1:numel(fcfg.ovr_clm)
    ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme{iEF} '_plt']);
    
    fcfg.all_ele{iEF} = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm(iEF)+2)))+1,2);
end
fcfg.all_ele = unique(cat(1,fcfg.all_ele{:}));

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};

cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.sel_ele;
cfg.all_ele   = fcfg.all_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_col = fcfg.nsl_col;

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'NEWFIGS' '/'];
cfg.sve_pre   = ['middle_pic_figure5'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_lat']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_lat']};
cfg.sve_pre   = ['middle_pic_figure5_lat'];
mmil_ieeg_sensor_location_plot_v4(cfg);

cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_ven']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_ven']};
cfg.sve_pre   = ['middle_pic_figure5_ven'];
mmil_ieeg_sensor_location_plot_v4(cfg);


