clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

%%%%%%%%%%%%%%%%%%% FOR THE PAPER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clinical Data
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
cfg.inc_clm = {'age' 'onset' 'sex' 'handedness' 'wada' };
mmil_ovr_cln(cfg);

%% Behavior
SL_behavior

%% Overall Electrodes
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
mmil_ovr_ele_loc(cfg);

%% Tables
SL_table;

%% Figure 2
SL_figure2;

%% Figure 3
SL_figure3;

%% Figure 3
SL_figure3_5;

%% Figure 4
SL_figure4;

%% Figure 5
SL_figure5;

%% Figure 6
SL_figure6;

%% Figure 7
SL_figure7

%% Figure 8
SL_figure8_v2

%%%%%%%%%%%%%%%%%%% Paper Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SL_paper_stats_FisherExactTest_v2
SL_paper_stats_Timing
SL_plv_stt

%%%%%%%%%%%%%%%%%%% Analyses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Connectivity
SL_connectivity;

%% Machine Learning
for iS = 7:numel(sbj_nme_hld)
    
    cfg = [];
    cfg.sbj_nme = sbj_nme_hld{iS};
    
    cfg.clr_fld = sbj_clr_hld{iS};
    cfg.dat_fld = [sbj_dat_hld{iS} '/' 'epoch_data'];
    
    cfg.typ     = 'hgp';
    
    cfg.pre_fix = {'ORT_FFN'   };
    cfg.alt_eve = {'trialinfo' };
    cfg.eve     = {[1 3]       };
    cfg.frq     = {[70 170]    };
    
    cfg.tme = [0.100 0.500];
    cfg.wdw = .050;
    cfg.stp = .050;
    
    mmil_machine_learning2(cfg)
    
end

for iS = 1:numel(sbj_nme_hld)
    
    cfg = [];
    cfg.sbj_nme = sbj_nme_hld{iS};
    
    cfg.clr_fld = sbj_clr_hld{iS};
    cfg.dat_fld = [sbj_dat_hld{iS} '/' 'epoch_data'];
    
    cfg.typ     = 'hgp';
    
    cfg.pre_fix = {'AUD_NSE'   };
    cfg.alt_eve = {'trialinfo' };
    cfg.eve     = {[1 4]       };
    cfg.frq     = {[70 170]    };
    
    cfg.tme = [0.450 0.950];
    cfg.wdw = .050;
    cfg.stp = .050;
    
    mmil_machine_learning2(cfg)
    
end

for iS = 1:numel(sbj_nme_hld)
    
    cfg = [];
    cfg.sbj_nme = sbj_nme_hld{iS};
    
    cfg.clr_fld = sbj_clr_hld{iS};
    cfg.dat_fld = [sbj_dat_hld{iS} '/' 'epoch_data'];
    
    cfg.typ     = 'hgp';
    
    cfg.pre_fix = {'MTC_MSS'   };
    cfg.alt_eve = {'trialinfo' };
    cfg.eve     = {[1 2]       };
    cfg.frq     = {[70 170]    };
    
    cfg.tme = [0.450 1.050];
    cfg.wdw = .050;
    cfg.stp = .050;
    
    mmil_machine_learning2(cfg)
    
end

cfg = [];

cfg.clr_fld = sbj_clr_hld{1};
cfg.pre_fix = {'ORT_FFN'   'AUD_NSE'   'MTC_MSS' };
cfg.eve     = {[1 3]       [1 4]       [1 2]};

cfg.typ = 'hgp';
cfg.eff_nme = 'pap_rsp_950';
cfg.eff_clm = 1;

mmil_total_machine_learning(cfg);

%% Re-run originals
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(:);
cfg.ana_nme = {'Responsitivity'     ...
               'Language_Selectivity' ...
               'Control_Selectivity' ...
               'Stimulus_Match_Mismatch' ...
               'Grapheme_Phoneme' ...
               'PostStimulus_Match_Mismatch' ...
               'Word_NonWord' };
cfg.typ     = {'hgp'};
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(:);
cfg.ana_nme = {'Responsitivity'     ...
               'Language_Selectivity' ...
               'Match_Mismatch' ...
               'Grapheme_Phoneme' ...
               'Late_Match_Mismatch' ...
               'Word_NonWord' };
cfg.typ     = {'lfp' 'hgp'};
mmil_overall_analysis(cfg);

%% Initial ANOVA
for sbj_num = 1:numel(sbj_nme_hld);
    
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
%     SL_pap_eff_hgp_rsp(cfg)
    
%     SL_pap_eff_lfp_act(cfg)
    
    % Total
    SL_pap_ana_hgp_rsp(cfg)
    
%     SL_pap_ana_lfp_rsp(cfg)
    
end

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(:);
cfg.ana_nme = { ... 'pap_rsp_1500' ...
                'pap_anv_1500' };
cfg.typ     = {'hgp' 'lfp'};
mmil_overall_analysis(cfg);

%% Initial Channels
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
%     SL_pap_eff_hgp_act(cfg)
    
    % Analysis
    SL_pap_ana_hgp_act(cfg)    

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:8);
cfg.ana_nme = { 'pap_rsp_950'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:8);
cfg.ana_nme = { 'pap_rsp_950'};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% Language Selective
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask
    SL_msk_lng(cfg);
    
    % Events
    SL_lng_eve(cfg);
    
    % Effects
    SL_pap_eff_hgp_lng(cfg);

%     SL_pap_eff_lfp_lng(cfg);
    
    % Analysis
    SL_pap_ana_hgp_lng(cfg);    

%     SL_pap_ana_lfp_lng(cfg); 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:8);
cfg.ana_nme = { 'pap_lng_950'};
cfg.typ     = { 'hgp' 'lfp' };
mmil_overall_analysis(cfg);

%% Control Selective
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};

    % Analysis
    SL_pap_ana_hgp_con(cfg);    

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:8);
cfg.ana_nme = { 'pap_con_950'};
cfg.typ     = { 'hgp' 'lfp' };
mmil_overall_analysis(cfg);

%% Match/Mismatch
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};

    % Mask
    SL_msk_mtc(cfg);
    
    % Effects
    SL_pap_eff_hgp_mtc(cfg);    
    
    % Analysis
    SL_pap_ana_hgp_mtc(cfg);    

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:8);
cfg.ana_nme = { 'pap_mtc_1450'};
cfg.typ     = { 'hgp' 'lfp' };
mmil_overall_analysis(cfg);

%% Graphemes/Phonemes
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};

    % Mask
    SL_msk_phn(cfg);
    
    % Effects
    SL_pap_eff_hgp_phn(cfg);    
    
    % Analysis
    SL_pap_ana_hgp_phn(cfg);    

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:8);
cfg.ana_nme = { 'pap_wrd_950'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

%% Word/NonWord
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};

    % Mask
    SL_msk_wrd(cfg);
    
    % Effects
    SL_pap_eff_hgp_wrd(cfg);    
    
    % Analysis
    SL_pap_ana_hgp_wrd(cfg);    

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:8);
cfg.ana_nme = { 'pap_wrd_950'};
cfg.typ     = { 'hgp'};
mmil_overall_analysis(cfg);

%% Depth Analysis
for sbj_num = 2:numel(sbj_nme_hld);

    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL';
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Load
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/depth/';
    SL_load_depth(cfg)
    
    % Events
    cfg.sbj_dat_hld = cfg.out_pth;
    cfg.tsk     = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
    SL_events_depth(cfg)
    
    % Stats
    
    
    % Plot
    SL_plot_depth(cfg)

end






