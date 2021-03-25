clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

%% Overall Electrodes
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
mmil_ovr_ele_loc(cfg);

%% Behavior
SL_performance

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

%% Initial Channels
for sbj_num = 5;
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
    SL_pap_eff_hgp_act(cfg)

    SL_pap_eff_lfp_act(cfg)
    
    % Analysis
    SL_pap_ana_hgp_act(cfg)    

    SL_pap_ana_lfp_act(cfg) 

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
for sbj_num = 5;
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask
    SL_msk_lng(cfg);
    
    % Effects
    SL_pap_eff_hgp_lng(cfg);

    SL_pap_eff_lfp_lng(cfg);
    
    % Analysis
    SL_pap_ana_hgp_lng(cfg);    

    SL_pap_ana_lfp_lng(cfg); 

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
