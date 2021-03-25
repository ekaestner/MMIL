clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

%%%%%%%%%%%%%%%%%%% FOR THE PAPER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clinical Data
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.inc_clm = {'age' 'onset' 'sex' 'handedness' 'wada' };
mmil_ovr_cln(cfg);

%% Behavior
FW_performance

%% Overall Electrodes
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
mmil_ovr_ele_loc(cfg);

%% Tables
FW_table1;

%% Included Regions
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_rsp_600';
cfg.eff_clm = 1;
cfg.eff_col = {{'dark red' 'red' 'neon red'}};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Inferior Parietal' 'Superior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
mmil_include_plot(cfg)

% Colorbar
pcfg = [];
pcfg.col_map = {'dark red' 'bright red'};
pcfg.out_dir = ['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure2'];
pcfg.col_bar = [0 0.5];
mmil_color_bar(pcfg)

%% Figure2
FW_figure2;

%% Figure3
FW_figure3;

%% Figure 4
FW_figure4;

%% Figure 5
FW_figure5;

%% Figure 6

%% Figure 7
FW_figure7;

%%%%%%%%%%%%%%%%%%% Paper Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FW_paper_stats_FisherExactTest
FW_paper_stats_Electrode_v2
FW_paper_stats_Timing_v2
FW_paper_stats_Stimulation

%%%%%%%%%%%%%%%%%%% Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Effectisze
for sbj_num = 1:numel(sbj_nme_hld);
    
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    FW_effectsize_update(cfg)
    
end

%% Connectivity
FW_plv

%% Machine Learning
for iS = 18:numel(sbj_nme_hld)
    
    cfg = [];
    cfg.sbj_nme = sbj_nme_hld{iS};

    cfg.clr_fld = sbj_clr_hld{iS};
    cfg.dat_fld = sbj_dat_hld{iS};
    
    cfg.typ = 'hgp';
    
    cfg.pre_fix = {'WRD_FFN'   'WRD_CSS'   'ORT_OVR' };
    cfg.alt_eve = {'trialinfo' 'trialinfo' 'trialinfo'};
    cfg.eve     = {[3 6]       [3 5]       [3 5 6]}; 
    cfg.frq     = {[70 170]    [70 170]    [70 170]};
    
    cfg.tme = [0.100 0.500];
    cfg.wdw = .050;
    cfg.stp = .050;
    
    mmil_machine_learning2(cfg)
    
end

cfg = [];

cfg.clr_fld = sbj_clr_hld{1};
cfg.pre_fix = {'WRD_FFN'   'WRD_CSS'   'ORT_OVR' };
cfg.eve     = {[3 6]       [3 5]       [3 5 6]};

cfg.typ = 'hgp';
cfg.eff_nme = 'pap_rsp_600';
cfg.eff_clm = 1;

mmil_total_machine_learning(cfg);

%% Original
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = {'eff_600_lap_01'     ...
               'GrainSize' ...
               'NovelN4' ...
               'UnfamiliarStimuli' ...
               'ltr_ffn_600' ...
               'wrd_nwd_600' ...
               'nov_rep_600' ...
               'lng_crr_600' ...
               'big_crr_600' ...
               'frq_crr_600' ...
               'ngh_crr_600' };
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = {'eff_600_lap_01'     ...
               'GrainSize' ...
               'NovelN4' ...
               'ltr_ffn_600' ...
               'wrd_nwd_600' ...
               'nov_rep_600' ...
               'lng_crr_600' ...
               'big_crr_600' ...
               'frq_crr_600' ...
               'ngh_crr_600' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% Initial Channels
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Stats
%     FW_stats
    
    % Effects
    FW_pap_eff_hgp_act(cfg)

    FW_pap_eff_lfp_act(cfg)
    
    % Analysis
    FW_pap_ana_hgp_act(cfg)    

    FW_pap_ana_lfp_act(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_rsp_600'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_rsp_600'};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% Word/Consonants/False-Font
% for iS = 1:numel(sbj_nme_hld)
%     
%     cfg = [];
%     
%     % Setting up Variables
%     cfg.sbj_nme     = sbj_nme_hld{sbj_num};
%     cfg.sbj_dat_hld = sbj_dat_hld{sbj_num};
%     cfg.clr_fld     = sbj_clr_hld{sbj_num};
%     
%     % Statistics
%     FW_add_stt(cfg)
%     
% end

for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask Stats
    FW_msk_wrd(cfg)
    
    % Effects
    FW_pap_eff_hgp_wrd(cfg)

    FW_pap_eff_lfp_wrd(cfg)
    
    % Analysis
    FW_pap_ana_hgp_wrd(cfg)    

    FW_pap_ana_lfp_wrd(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_wrd_600'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_wrd_600'};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

% Region Timing
cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'pap_wrd_600';
cfg.inc_eve = {'parsopercularis' 'parstriangularis' 'parsorbitalis' 'rostral-middlefrontal' 'middle-middlefrontal' 'caudal-middlefrontal'};
cfg.out_nme = 'IFG_pap_wrd_600';
mmil_region_examination(cfg)

cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'pap_wrd_600';
cfg.inc_eve = {'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-MTG' 'caudal-STG' 'middle-STG' };
cfg.out_nme = 'AudSnsMot_pap_wrd_600';
mmil_region_examination(cfg)

cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'pap_wrd_600';
cfg.inc_eve = {'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' };
cfg.out_nme = 'VntRte_pap_wrd_600';
mmil_region_examination(cfg)

%% Unfamiliar/False-Font
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Analysis
    FW_pap_ana_hgp_con(cfg)    

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_con_600'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

%% Sublexical/Lexical
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask Stats
%     FW_msk_wrd(cfg)
    
    % Effects
    FW_pap_eff_hgp_lex(cfg)

    FW_pap_eff_lfp_lex(cfg)
    
    % Analysis
    FW_pap_ana_hgp_lex(cfg)    

    FW_pap_ana_lfp_lex(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_lex_600'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_lex_600'};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

% Region Timing
cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'pap_lex_600';
cfg.inc_eve = {'parsopercularis' 'parstriangularis' 'parsorbitalis' 'rostral-middlefrontal' 'middle-middlefrontal' 'caudal-middlefrontal'};
cfg.out_nme = 'IFG_pap_lex_600';
mmil_region_examination(cfg)

cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'pap_lex_600';
cfg.inc_eve = {'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-MTG' 'caudal-STG' 'middle-STG' };
cfg.out_nme = 'AudSnsMot_pap_lex_600';
mmil_region_examination(cfg)

cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'pap_lex_600';
cfg.inc_eve = {'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'parahippocampal' 'cuneus' 'lingual'};
cfg.out_nme = 'VntRte_pap_lex_600';
mmil_region_examination(cfg)

%% Orthographic
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Add orthographic stats
%     FW_add_ort_stt(cfg)
    
    % Mask Stats
%     FW_msk_wrd(cfg)
    
    % Effects
%     FW_pap_eff_hgp_ort(cfg)

%     FW_pap_eff_lfp_ort(cfg)
    
    % Analysis
    FW_pap_ana_hgp_ort(cfg)    

    FW_pap_ana_lfp_ort(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_ort_600'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_ort_600'};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% Motor Response Control
% for sbj_num = 1:numel(sbj_nme_hld);
%            
%     cfg = [];
%     
%     % Setting up Variables
%     cfg.sbj_nme     = sbj_nme_hld{sbj_num};
%     cfg.dat_fld     = sbj_dat_hld{sbj_num};
%     cfg.clr_fld     = sbj_clr_hld{sbj_num};
%     
%     % Load
%     cfg.out_pth = '/space/seh8/1/halgdev/projects/mmilanguage/FWaddon';
%     FW_load(cfg)
%         
%     % Stats
%     FW_motor_stats(cfg)
%     
%     % Plot
%     FW_motor_plot(cfg)
%     
% end

for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = '/space/seh8/1/halgdev/projects/mmilanguage/FWaddon';
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
        
    % Mask Stats
%     FW_msk_wrd(cfg)
    
    % Effects
    FW_pap_eff_hgp_rsp(cfg)

    FW_pap_eff_lfp_rsp(cfg)
    
    % Analysis
    FW_pap_ana_hgp_rsp(cfg)    

    FW_pap_ana_lfp_rsp(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_rsp_600'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_rsp_600'};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% Depth Analysis
for sbj_num = 28:numel(sbj_nme_hld);

    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW';
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Load
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data/depth';
    FW_load_depth(cfg)
    
    % Events
    cfg.sbj_dat_hld = cfg.out_pth;
    FW_events_depth(cfg)
    
    % Stats
    
    
    % Plot
    FW_plot_depth(cfg)

end

