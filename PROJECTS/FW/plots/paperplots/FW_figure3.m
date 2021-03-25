function FW_figure3

fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
fcfg.tsk     = 'FW';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = 'pap_lex_600';
fcfg.eff_clm = [ 3 1 ]; % 2 
fcfg.eff_col = { rgb('yellow') rgb('orange') }; % rgb('blue')
fcfg.eff_lbl = { 'Lexical Frequency' 'Repetition' }; % 'Bigram Frequency'
fcfg.nsl_col = {rgb('blue')};

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = 'pap_rsp_600'; % 'pap_rsp_600'
fcfg.ovr_clm = 1; % 1

%%
pcfg = [];
pcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/'];
pcfg.lne_col = {rgb('red') rgb('reddish grey') } ;
pcfg.lne_lbl = {'Word'     'False-Font'        };
mmil_lne_leg(pcfg)

%% Middle Portion
% Selective Electrodes
eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme '/' 'subjects' '/' 'total' '/' fcfg.eff_nme '_plt']);

for iEF = 1:numel(fcfg.eff_clm)  
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

% Overall Electrodes
ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme '_plt']);

fcfg.all_ele = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm+2)))+1,2);

% Restrict
% tbl{1} = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme '/' 'total' '/' 'total_' fcfg.eff_nme '_lhs_table_plot']);
% tbl{1} = tbl{1}(1:end-1,[1 (fcfg.eff_clm-1)*5+2:fcfg.eff_clm*5+1]);
% rmv_ind = [];
% for iR = 2:size(tbl{1},1)
%     if tbl{1}{iR,5}<2
%         rmv_ind = [rmv_ind ; iR];
%     end
% 
% end
% tbl{1}(rmv_ind,:) = [];
% tbl{1}(:,5:end) = [];
% tbl{1} = mmil_order_table_reverse(tbl{1});
% 
% tbl{2} = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme '/' 'total' '/' 'total_' fcfg.eff_nme '_rhs_table_plot']);
% tbl{2} = tbl{2}(1:end-1,[1 (fcfg.eff_clm-1)*5+2:fcfg.eff_clm*5+1]);
% rmv_ind = [];
% for iR = 2:size(tbl{2},1)
%     if tbl{2}{iR,5}<2
%         rmv_ind = [rmv_ind ; iR];
%     end
% end
% tbl{2}(rmv_ind,:) = [];
% tbl{2}(:,5:end) = [];
% tbl{2} = mmil_order_table_reverse(tbl{2});

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

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sep_str         = ',';

cfg.sve_img   = 'eps';
cfg.sve_loc   = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/NEW_FIGS/Percentage_Electrodes'
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_ele;
mmil_loc_dot(pcfg)

%% Middle Portion ALTERNATE
% Letter
fcfg.ovr_nme = 'pap_wrd_600'; % 'pap_rsp_600'
fcfg.ovr_clm = 1; % 1
ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme '_plt']);
fcfg.all_ele = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm+2)))+1,2);

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
cfg.nsl_ele = {rgb('purple')};

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sep_str         = ',';

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'OverlapCheck'];
cfg.sve_pre   = ['ltr_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

% Word
fcfg.ovr_nme = 'pap_wrd_600'; % 'pap_rsp_600'
fcfg.ovr_clm = 3; % 1
ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme '_plt']);
fcfg.all_ele = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm+2)))+1,2);

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
cfg.nsl_ele = {rgb('red')};

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sep_str         = ',';

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'OverlapCheck'];
cfg.sve_pre   = ['wrd_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

% False-font
fcfg.ovr_nme = 'pap_con_600'; % 'pap_rsp_600'
fcfg.ovr_clm = 1; % 1
ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme '_plt']);
fcfg.all_ele = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm+2)))+1,2);

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
cfg.nsl_ele = {rgb('reddish grey')};

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sep_str         = ',';

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'OverlapCheck'];
cfg.sve_pre   = ['fls_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

%% SECOND ALTERNATE OVERLAP
fcfg.eff_nme = {'pap_wrd_600' 'pap_con_600' 'pap_lex_600'};
fcfg.eff_clm = { [ 1 3 ]      [1]           [1] };
fcfg.eff_col = { rgb('purple') rgb('bright red') rgb('reddish grey') rgb('orange') };
fcfg.eff_lbl = { 'Letters' 'Words' 'FalseFont' 'Repetition'};

fcfg.nsl_col = {rgb('maroon')};

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = 'pap_rsp_600'; % 'pap_rsp_600'
fcfg.ovr_clm = 1; % 1

% Selective Electrodes
eff_txt_one = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{1} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{1} '_plt']);
eff_txt_two = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{2} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{2} '_plt']);
eff_txt_thr = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{3} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{3} '_plt']);

fcfg.sel_ele{1} = eff_txt_one(find(cell2mat(eff_txt_one(2:end,fcfg.eff_clm{1}(1)+2)))+1,2);
fcfg.sel_ele{2} = eff_txt_one(find(cell2mat(eff_txt_one(2:end,fcfg.eff_clm{1}(2)+2)))+1,2);
fcfg.sel_ele{3} = eff_txt_two(find(cell2mat(eff_txt_two(2:end,fcfg.eff_clm{2}(1)+2)))+1,2);
fcfg.sel_ele{4} = eff_txt_thr(find(cell2mat(eff_txt_thr(2:end,fcfg.eff_clm{3}(1)+2)))+1,2);

fcfg.sel_ele{1} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{4});
fcfg.sel_ele{2} = intersect(fcfg.sel_ele{2},fcfg.sel_ele{4});
fcfg.sel_ele{3} = intersect(fcfg.sel_ele{3},fcfg.sel_ele{4});

% Overall Electrodes
ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme '_plt']);

fcfg.all_ele = intersect(ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm+2)))+1,2),unique([fcfg.sel_ele{1} ; fcfg.sel_ele{2} ; fcfg.sel_ele{3}]));

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
cfg.nsl_ele = fcfg.nsl_col;

cfg.sep_str         = ',';

cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'RepetitionOverlap'];
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.inc_prc = { [fcfg.clr_fld '/' 'electrode_location_files' '/' 'fsaverage' '/'  'label' '/' 'lh.aparc.split.annot'] ...
                [fcfg.clr_fld '/' 'electrode_location_files' '/' 'fsaverage' '/'  'label' '/' 'rh.aparc.split.annot'] };
            
mmil_ieeg_sensor_location_plot_v4(cfg);

%% MOVE TO SL
% Language
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
cfg.sve_loc = 'figure3';
cfg.tsk     = 'SL';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_lng_950';
cfg.eff_clm = [ 1 2 ];
cfg.eff_col = { {'dark red' 'red' 'neon red'} {'dark blue' 'blue' 'bright blue'} };
cfg.top_pct = 0.50;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
mmil_include_plot(cfg)

% Control
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
cfg.sve_loc = 'figure3';
cfg.tsk     = 'SL';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_con_950';
cfg.eff_clm = [ 1 2 ];
cfg.eff_col = { {'dark magenta' 'magenta' 'bright magenta'} {'dark cyan' 'cyan' 'bright cyan'} };
cfg.top_pct = 0.50;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
mmil_include_plot(cfg)

%% Example Channels
% Example Channels
sbj = { 'NY255'            'NY226'           'NY068'            'NY090'           'NY190'           'NY090'            'NY226'              'NY255'              'NY313'                'NY077'   };
chn = { 'LTO5'             'LMT3'            'G34'              'G38'             'GA40'            'G29'              'G21'                'RF6'                'RO04'                 'G05'     };
loc = { 'lft_ven'          'lft_ven'         'lft_stg'          'lft_stg'         'lft_sup'         'lft_pos'          'lft_pre'            'rgh_ifg'            'rgh_ven'              'rgh_ven' };
typ = { [2]                [2]               [2]                [2]               [2]               [2]                [2]                  [2]                  [2]                    [2] };
ylm = { [-3*10^10 9*10^10] [-7*10^9 21*10^9] [-1*10^6 3*10^6]   [-3*10^6 9*10^6]  [-7*10^9 21*10^9] [-2*10^6 6*10^6]   [-1*10^10 3*10^10]	[-1*10^10 3*10^10]   [-1*10^11 3*10^11]     [-8*10^7 24*10^7] };


tot_sbj = unique(sbj);

% Parameters
nme             = { 'repetition'                                   'Bigram'                                                  'Frequency' };

alt_eve = { 'trialinfo'                                            'wrd_big_med'                                             'wrd_frq_med' };
eve     = { [3 4 6]                                                [141 142 6]                                               [121 122 6] };
col_ord = { { rgb('bright red') rgb('orange') rgb('reddish grey')} {rgb('dark blue') rgb('bright blue') rgb('reddish grey')} {rgb('dark yellow') rgb('bright yellow')  rgb('reddish grey')} };
                     
stt_dat = { { 'vis_stm_01'            'vis_old' }                 { 'vis_stm_01' 'vis_wrd_big' }                            {'vis_stm_01' 'vis_wrd_frq'} };
stt_col = { { rgb('black')            rgb('black') }              { ft_stt_col(rgb('black')) ft_stt_col(rgb('light blue'))} {ft_stt_col(rgb('black')) ft_stt_col(rgb('light yellow'))} };
stt_cmp = { { '0%5'                   '3vs4'                    } { '0%5' '141vs142' }                                      { '0%5' '121vs122' } };

% Put together plot
for iS = 1:numel(tot_sbj)
    
    % Load
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [fcfg.dat_fld '/' tot_sbj{iS} '_FW_overall_data.mat'];
    bcc_dat  = ft_func([],cfg);
   
    bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_big_med(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==6) = bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==6);
    
    bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_big_med_rep = bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_big_med;
    bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_big_med_rep(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==4) = bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==4);
    
    bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_frq_med(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==6) = bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==6);
    
    bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_frq_med_rep = bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_frq_med;
    bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.wrd_frq_med_rep(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==4) = bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_eve.trialinfo==4);
    
    plt_loc = find(ismember(sbj,[tot_sbj{iS}]));
    
    cfg     = [];
    cfg.stt     = { 'vis_old'    'vis_wrd_frq'    }; %
    cfg.stt_msk = { 'vis_stm_01' 'vis_stm_01' }; %
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
    
    for iP = 1:numel(plt_loc)
        
        % Put together 1-effect plot
%         for iT = 1:numel(typ{plt_loc(iP)})
%             
%             cfg = [];
%             
%             cfg.type      = 'chan';
%             cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,chn{plt_loc(iP)}))};
%             cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
%             cfg.dat_loc   = 1;
%             
%             cfg.plt_dim   = [1 1];
%             
%             cfg.y_lim     = 'maxmin';
%             
%             cfg.lgd       = 0;
%             cfg.std_err   = 1;
%             
%             cfg.alt_eve = alt_eve{typ{plt_loc(iS)}(iT)};
%             cfg.eve     = eve{typ{plt_loc(iS)}(iT)};
%             cfg.lnstyle.col_ord = col_ord{typ{plt_loc(iS)}(iT)};
%             
%             cfg.stt_dat = stt_dat{typ{plt_loc(iS)}(iT)};
%             cfg.stt_col = {stt_col{typ{plt_loc(iS)}(iT)}};
%             cfg.stt_lab = 'stt_lab';
%             cfg.stt_cmp = {stt_cmp{typ{plt_loc(iS)}(iT)}};
%             
%             cfg.v_lne       = [0 0.2 0.4];
%             cfg.v_lne_wdt    = [3 1 1];
%             cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black')};
%             
%             cfg.x_lim = [-0.2 0.8];
%             
%             cfg.print      = 1;
%             cfg.nofig      = 1;
%             cfg.print_type = 'png';
%             cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' nme{typ{plt_loc(iS)}(iT)}];
%             cfg.prefix     = [tot_sbj{iS} '_' nme{typ{plt_loc(iS)}(iT)} '_' loc{plt_loc(iP)}];
%                         
%             mmil_ieeg_sensor_plot_v5(cfg)
%                         
%             cfg.print_type = 'eps';
%             mmil_ieeg_sensor_plot_v5(cfg)
%             
%         end
        
        % Dual Effect - Lexical
%         if all(ismember(typ{plt_loc(iP)},[1 3]))
            
            cfg = [];
            
            cfg.type      = 'chan';
            cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,chn{plt_loc(iP)}))};
            cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
            cfg.dat_loc   = 1;
            
            cfg.plt_dim   = [1 1];
            
            cfg.y_lim     = ylm{plt_loc(iP)};
            
            cfg.lgd       = 0;
            cfg.std_err   = 1;
            
            cfg.alt_eve = 'wrd_frq_med_rep';
            cfg.eve     = [4 6 121 122];
            cfg.lnstyle.col_ord = {rgb('orange') rgb('reddish grey') rgb('dark yellow')-[0.35 0.27 0.01] rgb('bright yellow')+[0 0 0.30] } ;
            
            cfg.stt_dat = { 'vis_stm_01'   'vis_old_msk' 'vis_wrd_frq_msk'};
            cfg.stt_col = { { rgb('black') rgb('orange') rgb('yellow')} };
            cfg.stt_lab = 'stt_lab';
            cfg.stt_cmp = { { '0%3'        '3%6'         '6%9'   } }; % { { '0%5'                   '122vs4'          '121vs122'   } }
            
            cfg.v_lne       = [0 0.2 0.4];
            cfg.v_lne_wdt   = [12 4 4];
            cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black')};
            
            cfg.x_lim = [-0.1 0.7];
            
            cfg.print      = 1;
            cfg.nofig      = 1;
            cfg.print_type = 'png';
            cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'lexical_overlap'];
            cfg.prefix     = [tot_sbj{iS} '_' 'lexical_overlap' '_' loc{plt_loc(iP)} '_'];
                        
            mmil_ieeg_sensor_plot_v5(cfg)
                        
            cfg.print_type = 'eps';
            mmil_ieeg_sensor_plot_v5(cfg)
            
%         end
        
        % Dual Effect - Bigram
%         if all(ismember(typ{plt_loc(iP)},[1 2]))
%             
%                         cfg = [];
%             
%             cfg.type      = 'chan';
%             cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,chn{plt_loc(iP)}))};
%             cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
%             cfg.dat_loc   = 1;
%             
%             cfg.plt_dim   = [1 1];
%             
%             cfg.y_lim     = 'maxmin';
%             
%             cfg.lgd       = 0;
%             cfg.std_err   = 1;
%             
%             cfg.alt_eve = 'wrd_big_med_rep';
%             cfg.eve     = [4 6 141 142];
%             cfg.lnstyle.col_ord = {rgb('orange') rgb('reddish grey') rgb('dark blue') rgb('bright blue')} ;
%             
%             cfg.stt_dat = { 'vis_stm_01'             'vis_old' 'vis_wrd_big'};
%             cfg.stt_col = { { ft_stt_col(rgb('black')) ft_stt_col(rgb('black')) ft_stt_col(rgb('black')) } };
%             cfg.stt_lab = 'stt_lab';
%             cfg.stt_cmp = { { '0%5'                   '142vs4'          '141vs142'   } };
%             
%             cfg.v_lne       = [0 0.2 0.4];
%             cfg.v_lne_wdt    = [3 1 1];
%             cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black')};
%             
%             cfg.x_lim = [-0.2 0.8];
%             
%             cfg.print      = 1;
%             cfg.nofig      = 1;
%             cfg.print_type = 'png';
%             cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'lexical_overlap'];
%             cfg.prefix     = [tot_sbj{iS} '_' 'lexical_overlap' '_'];
%                         
%             mmil_ieeg_sensor_plot_v5(cfg)
%                         
%             cfg.print_type = 'eps';
%             mmil_ieeg_sensor_plot_v5(cfg)
%             
%         end
        
    end
    
end

% Line Legend
pcfg = [];
pcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/'];
pcfg.lne_col = {rgb('orange') rgb('reddish grey') rgb('dark blue')     rgb('bright blue')    rgb('dark yellow')-[0.35 0.27 0.01] rgb('bright yellow')+[0 0 0.30] } ;
pcfg.lne_lbl = {'Repetition'  'False-Font'        'LowFrequencyBigram' 'HighFrequencyBigram' 'LowFrequencyLexical'               'HighFrequencyLexical'};
mmil_lne_leg(pcfg)

%% REPETITION
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_lex_600';
cfg.eff_clm = 1;
cfg.top_pct = 0.30;
cfg.eff_col = {{'dark orange' 'orange' 'bright orange'}};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.sve_loc = 'figure3';
mmil_include_plot(cfg)

% Colorbar
pcfg = [];
pcfg.col_map = {'dark orange' 'orange' 'bright orange'};
pcfg.out_dir = ['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure3/repetition'];
pcfg.col_bar = [0 0.30];
mmil_color_bar(pcfg)

end