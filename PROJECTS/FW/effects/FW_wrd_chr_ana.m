clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

for sbj_num = [1:2 4:numel(sbj_nme_hld)]
    
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.sbj_dat_hld = sbj_dat_hld{sbj_num};
    cfg.dat_fld     = sbj_dat_hld{sbj_num};
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% Statistics
    FW_stats_word_char(cfg)
   
    %% Effects
    FW_effects_word_char_hgp(cfg)
    
    FW_effects_word_char_lfp(cfg)
    
    %% Effect Size
%     FW_word_char_effectsize(cfg)
    
%     FW_word_char_effectsize_lfp(cfg)
    
    %% Analysis
    FW_word_char_analysis_hgp(cfg)
    
    FW_word_char_analysis_lfp(cfg)
    
end