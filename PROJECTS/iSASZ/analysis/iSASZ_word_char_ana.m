clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

for sbj_num = [1 3:numel(sbj_nme_hld)]
    
    cfg = [];
    
    %% Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.sbj_dat_hld = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    %% Statistics
    cfg.loc = 'local';
    
    iSASZ_stats_word_char(cfg)
   
    %% Effects
    iSASZ_effects_word_char_hgp(cfg)
    
    iSASZ_effects_word_char_lfp(cfg)
    
    %% Effect Size
%     iSASZ_word_char_effectsize(cfg)
    
%     iSASZ_word_char_effectsize_lfp(cfg)
    
    %% Analysis
    iSASZ_word_char_analysis_visual_hgp(cfg)
    iSASZ_word_char_analysis_auditory_hgp(cfg)
    
    iSASZ_word_char_analysis_visual_lfp(cfg)
    iSASZ_word_char_analysis_auditory_lfp(cfg)
    
end