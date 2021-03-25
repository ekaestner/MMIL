clear; clc;

%% Letters
% 
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_wrd_600';

cfg.eff_clm = 1;
cfg.top_pct = 0.25;
cfg.eff_col = {{'purplish grey' 'purple' 'neon purple' }};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.sve_loc = 'figure2';
mmil_include_plot(cfg)

pcfg = [];
pcfg.col_map = [ rgb('purplish grey') ; rgb('purple') ; rgb('neon purple')];
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/NEW_FIGS';
pcfg.sve_pre = 'letters';
mmil_color_bar(pcfg)

%% Words
%
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_wrd_600';

cfg.eff_clm = 3;
cfg.top_pct = 0.25;
cfg.eff_col = {{'reddish grey' 'red' 'neon red'}};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
cfg.sve_loc = 'figure2';
mmil_include_plot(cfg)

pcfg = [];
pcfg.col_map = [ rgb('reddish grey') ; rgb('red') ; rgb('neon red') ];
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/NEW_FIGS';
pcfg.sve_pre = 'words';
mmil_color_bar(pcfg)

%% Example Waveforms
% Set-up Places to look
fcfg = [];

fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
fcfg.tsk     = 'FW';

fcfg.sbj_nme = {'NY226_FW' 'NY190_FW' 'NY086_FW'};
fcfg.chn_nme = {'LMT3'     'LPT02'    'PT05'};

fcfg.alt_eve = 'trialinfo';
fcfg.eff_typ = [3 5 6];
fcfg.col_ord = {rgb('bright red') rgb('purple') rgb('reddish grey')};
fcfg.eff_nme = {'Word' 'Letters' 'False-Font'};

fcfg.stt_dat = { 'vis_stm_01'   'vis_wrd_ffn_msk_01' 'vis_wrd_msk_01'};
fcfg.stt_col = { { rgb('black') rgb('purple')        rgb('red')} };
fcfg.stt_lab = 'stt_lab';
fcfg.stt_cmp = { { '0%3'        '3%6'         '6%9'   } };

fcfg.y_lim = [-5*10^8 18*10^9];
fcfg.x_lim = [-0.1 0.750];

% Line Plots
mmil_chk_dir(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/NEW_FIGS/LINES/']);

for iT = 1:numel(fcfg.chn_nme)

    cfg = [];
    cfg.load = 'yes';
    cfg.file = [fcfg.dat_fld '/' fcfg.sbj_nme{iT} '_overall_data.mat'];
    bcc_dat  = ft_func([],cfg);   
    
    cfg = [];
    
    cfg.type      = 'chan';
    cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,fcfg.chn_nme{iT}))};
    cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
    cfg.dat_loc   = 1;
    
    cfg.plt_dim   = [1 1];
    
    
    
    cfg.lgd       = 0;
    cfg.std_err   = 1;
    
    cfg.alt_eve         = fcfg.alt_eve;
    cfg.eve             = fcfg.eff_typ;
    cfg.lnstyle.col_ord = fcfg.col_ord;
    
    cfg.stt_dat = fcfg.stt_dat;
    cfg.stt_col = fcfg.stt_col;
    cfg.stt_lab = fcfg.stt_lab;
    cfg.stt_cmp = fcfg.stt_cmp;
    
    cfg.v_lne       = [0 0.2 0.4];
    cfg.v_lne_wdt   = [3 1 1];
    cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black')};
    
    cfg.y_lim = fcfg.y_lim;
    cfg.x_lim = fcfg.x_lim;
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/NEW_FIGS/LINES/' '/' ];
    cfg.prefix     = [fcfg.sbj_nme{iT} '_' fcfg.chn_nme{iT} '_'];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
    cfg.print_type = 'eps';
    mmil_ieeg_sensor_plot_v5(cfg)
    
end


