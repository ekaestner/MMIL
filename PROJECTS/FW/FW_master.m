clear; clc;

prj_dat_hld = '/home/ekaestne/PROJECTS/OUTPUT/FW';% '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

cfg = [];
cfg.clr_fld = sbj_clr_hld{1};
cfg.frs_avg_loc = 1;
mmil_mri_setup(cfg);
cfg.ch2_avg_loc = 1;
mmil_mri_setup(cfg);

%% Run Subjects
for sbj_num = 12:24 %numel(sbj_nme_hld);
    
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.sbj_dat_hld = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% Localization
    mmil_mri_setup(cfg);
    
    %% load the data
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data/';
    
%     FW_load(cfg)
%     FW_fast_load(cfg)
    %% fix the events
%     FW_events(cfg)
    
end

%% Statistics
for sbj_num = 25:numel(sbj_nme_hld);
   
    cfg = [];
        
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.sbj_dat_hld = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% Statistics
    FW_stats(cfg)
    
    %% Plot
%     fw_plot2(cfg)
    
end

%% Effects
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% Mask Stats
    FW_mask(cfg)
    
    %% Effects
    FW_effects_hgp(cfg)

    FW_effects_lfp(cfg)

    %% Effect Size
    FW_effectsize(cfg)
    
    FW_effectsize_lfp(cfg)
    
    %% Analysis
    FW_analysis_hgp(cfg)    

    FW_analysis_lfp(cfg) 

end

%% Overall Analyses
% Compile Results
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = {...'rsp_600' ...
               'eff_600_lap_01'     ...
               'ltr_ffn_600' ...
               'wrd_nwd_600' ...
               'nov_rep_600' ...
               'GrainSize' ...
               'NovelN4' ...
               };...'UnfamiliarStimuli' };
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = {'lng_crr_600' ...
               'big_crr_600'     ...
               'frq_crr_600' }; ...
               ... '' };
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = {'rsp_600' ...
               'eff_600_lap_01'     ...
               'ltr_ffn_600' ...
               'wrd_nwd_600' ...
               'nov_rep_600' ...
               'GrainSize' ...
               'NovelN4' ...
               'UnfamiliarStimuli' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = { 'lng_crr_600' ...
                'big_crr_600'     ...
                'frq_crr_600' }; ...
               ... '' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

% Compile Spreadsheet
cfg = [];
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld(1:32);
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_inc = {'eff_600'     ...
               'ltr_ffn_600' ...
               'wrd_nwd_600' ...
               'nov_rep_600' };
mmil_effect_spreadsheet(cfg);
           
% Chi-square
cfg = [];
cfg.chi_sqr_stt = {'Visual Language'      'vis_stm'     '!' 'eff_600' ; ...
                   'Letter_Selective'     'vis_ltr_msk' '5' 'ltr_ffn_600' ; ...
                   'False_Font_Selective' 'vis_ltr_msk' '6' 'ltr_ffn_600' ; ...
                   'Word_Selective'       'vis_wrd_msk' '3' 'wrd_nwd_600' ; 
                   'Non-Word_Selective'   'vis_wrd_msk' '5' 'wrd_nwd_600' ; 
                   'Novel_Selective'      'vis_old_msk' '3' 'nov_rep_600' ; 
                   'Repetition_Selective' 'vis_old_msk' '4' 'nov_rep_600' };
cfg.typ     = {'hgp'};
cfg.ele     = {'ecog'};
cfg.ele_typ = {'split'};
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld([1:24]);
cfg.tme_win = [0.050];
cfg.tme     = [0.000 0.450];
cfg.hem     = {'lhs' 'rhs'};
cfg.no_sub_fld = [];
cfg.cpy_plt    = 1;
mmil_chi_sqr(cfg)





