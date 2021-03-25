clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%
for sbj_num = [26 28 34:37 ]; % 1:16 17 18:23 24:37
    
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.sbj_dat_hld = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
%     iSASZ_word_char(cfg)
    
    %% MRI setup
    mmil_mri_setup(cfg);
    
    %% load the data
    cfg.tsk_hld = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' cfg.sbj_nme],'tsk');
    cfg.out_pth = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
%     iSASZ_initial_load3(cfg)
    
    %% Put together initial events
%     iSASZ_initial_events3(cfg)
    
end

% STATS %%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
for sbj_num = 23:numel(sbj_nme_hld);
    
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.fle_out_pth = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    cfg.ovr_wrt = 1;
    
    wat = 1;
    if exist([cfg.fle_out_pth '/' 'epoch_data' '/' cfg.sbj_nme '_overall_data.mat']); wat = 0; end
    while wat
        if exist([cfg.fle_out_pth '/' 'epoch_data' '/' cfg.sbj_nme '_overall_data.mat']); wat = 0; end
        pause(3600)
    end
    
    %% run the initial stats
    cfg.loc = 'local';
    
    iSASZ_initial_stats3(cfg)
    
    %% Run the initial plots
    iSASZ_initial_plot_v3(cfg)
        
end

% EFFECTS %%%%%%%%%%%%%%%%%%%%%%%%
for sbj_num = 1:24 %numel(sbj_nme_hld);
    
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% Mask Stats
    iSASZ_mask(cfg)
    
    %% Effects
    iSASZ_effects_hgp2(cfg)
    
    iSASZ_effects_lfp(cfg)

    %% Effect Size
%     iSASZ_effectsize_hgp(cfg)
%     
%     iSASZ_effectsize_lfp(cfg)
    
    %% Analysis
    iSASZ_analysis_bimodal_hgp(cfg)
    iSASZ_analysis_visual_hgp(cfg)
    iSASZ_analysis_auditory_hgp(cfg)

    iSASZ_analysis_bimodal_lfp(cfg)
    iSASZ_analysis_visual_lfp(cfg)
    iSASZ_analysis_auditory_lfp(cfg)

end

%% Overall Analyses
% Compile Results
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld([2:3 5:7 9:10 14:15 17:18 20:23]);
cfg.ana_nme = {'eff_800_bim' ...
               'ani_obj_800_bim' ...
               'nov_rep_800_bim' };
cfg.typ     = {'hgp'};
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld([1:7 9:12 14:15 17:18 20:23]);
cfg.ana_nme = {'eff_800_vis' ...
               'ani_obj_800_vis' ...
               'nov_rep_800_vis'};
cfg.typ     = {'hgp'};
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld([2:3 5:10 13:23]);
cfg.ana_nme = {'eff_800_aud' ...
               'ani_obj_800_aud' ...
               'nov_rep_800_aud' };
cfg.typ     = {'hgp'};
mmil_overall_analysis(cfg);

% Compile Spreadsheet
cfg = [];
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld([2:3 5:7 9:10 14:15 17:18 20:23]);
cfg.typ     = {'hgp'};
cfg.ana_inc = {'eff_800_bim' ...
               'ani_obj_800_bim' ...
               'nov_rep_800_bim' };
cfg.sht_nme = 'bimodal';
mmil_effect_spreadsheet(cfg);

cfg = [];
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld([1:7 9:12 14:15 17:18 20:23]);
cfg.typ     = {'hgp'};
cfg.ana_inc = { 'eff_800_vis' ...
                'ani_obj_800_vis' ...
                'nov_rep_800_vis' };
cfg.sht_nme = 'visual';
mmil_effect_spreadsheet(cfg);

cfg = [];
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld([2:3 5:10 13:23]);
cfg.typ     = {'hgp'};
cfg.ana_inc = { 'eff_800_aud' ...
                'ani_obj_800_aud' ...
                'nov_rep_800_aud' };
cfg.sht_nme = 'auditory';
mmil_effect_spreadsheet(cfg);

% Chi-square
cfg = [];
cfg.chi_sqr_stt = { 'Visual_Responsive' 'vis_ovr_stt' '101' 'eff_800_vis' ; ...
                    'Visual_Object'   'vis_ani_obj_stt_msk' '122' 'ani_obj_800_vis' ; ...
                    'Visual_Animal'   'vis_ani_obj_stt_msk' '121' 'ani_obj_800_vis' ; ...
                    'Visual_Novel'        'vis_new_old_stt_msk' '111' 'nov_rep_800_vis' ; ...
                    'Visual_Repetition'   'vis_new_old_stt_msk' '112' 'nov_rep_800_vis' } ;
cfg.typ     = {'hgp'};
cfg.ele     = {'ecog'};
cfg.ele_typ = {'split'};
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld([1:7 9:12 14:15 17 20:23]);
cfg.tme_win = [0.050];
cfg.tme     = [0.000 0.450];
cfg.hem     = {'lhs' 'rhs'};
cfg.no_sub_fld = [];
cfg.cpy_plt    = 1;
cfg.sht_nme = 'visual';
mmil_chi_sqr(cfg)

cfg = [];
cfg.chi_sqr_stt = { 'Auditory_Responsive' 'aud_ovr_stt' '201' 'eff_800_aud' ; ...
                    'Auditory_Object' 'aud_ani_obj_stt_msk' '222' 'ani_obj_800_aud' ; ...
                    'Auditory_Animal' 'aud_ani_obj_stt_msk' '221' 'ani_obj_800_aud' ; ...
                    'Auditory_Novel'      'vis_new_old_stt_msk' '211' 'nov_rep_800_aud' ; ...
                    'Auditory_Repetition' 'vis_new_old_stt_msk' '212' 'nov_rep_800_aud' } ;
cfg.typ     = {'hgp'};
cfg.ele     = {'ecog'};
cfg.ele_typ = {'split'};
cfg.dat_fld = sbj_dat_hld{1};
cfg.clr_fld = sbj_clr_hld{1};
cfg.inc_sbj = sbj_nme_hld([2:3 5:10 13:23]);
cfg.tme_win = [0.050];
cfg.tme     = [0.000 0.450];
cfg.hem     = {'lhs' 'rhs'};
cfg.no_sub_fld = [];
cfg.cpy_plt    = 1;
cfg.sht_nme = 'auditory';
mmil_chi_sqr(cfg)

