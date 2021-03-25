clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
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
for sbj_num = 6:numel(sbj_nme_hld);
    
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.sbj_dat_hld = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% Localization
%     mmil_mri_setup(cfg);
    
    %% load the data
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';
    
    SL_initial_load2(cfg)
    
    %% Events
    cfg.out_pth = [prj_dat_hld '/' 'epoch_data'];
    cfg.tsk     = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
    sl_events(cfg)
    
end

%% Run Stats Subjects
for sbj_num = 1:numel(sbj_nme_hld);
    
    cfg = [];
    
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.out_pth     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
            
%     wat = 1;
%     while wat
%         if exist([sbj_dat_hld{sbj_num} '/' 'epoch_data' '/' cfg.sbj_nme '_overall_data.mat']); wat = 0; end
%         pause(1800)
%     end
    
    %% Behavior
    
        
    %% Initial stats
    sl_initial_stats_2(cfg)
    
    %% Initial plots
%     SL_initial_plot_2(cfg)  
    
end
  
%% Effects
cfg = [];
for sbj_num = 1:numel(sbj_nme_hld);
    
    cfg.sbj_nme = sbj_nme_hld{sbj_num};
    cfg.dat_fld = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld = sbj_clr_hld{sbj_num};
    
    cfg.ovr_wrt = 1;
    
    %% Msk stats
    SL_mask(cfg)
    
    %% Note effects
    SL_effects_hgp(cfg)
    
    SL_effects_lfp(cfg)
    
    %% Effectsize
%     SL_effectsize_hgp(cfg)
    
%     SL_effectsize_lfp(cfg)

    %% Analysis
    SL_analysis_hgp(cfg)
      
    SL_analysis_lfp(cfg)

end

%% Overall Analyses
% Compile Results
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(:);
cfg.ana_nme = {'Responsitivity'     ...
               'Selectivity' ...
               'Language_Selectivity' ...
               'Control_Selectivity' ...
               'Stimulus_Match_Mismatch' ...
               'Match_Mismatch' };
cfg.typ     = {'lfp' 'hgp'};
mmil_overall_analysis(cfg);

% Compile Spreadsheet
cfg = [];
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld(:);
cfg.typ     = {'hgp'};
cfg.ana_inc = {'Responsitivity'     ...
               'Language_Selectivity' ...
               'Control_Selectivity' ...
               'Match_Mismatch' };
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

