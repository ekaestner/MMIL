clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

%% Initial Channels
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
    iSASZ_pap_eff_hgp_act(cfg)

    iSASZ_pap_eff_lfp_act(cfg)
    
    % Analysis
    iSASZ_pap_ana_hgp_act(cfg)    

    iSASZ_pap_ana_lfp_act(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = { ''};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = { ''};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% Repetition Effects
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask
    iSASZ_msk_rep(cfg)    
    
    % Effects
    iSASZ_pap_eff_hgp_rep(cfg)

    iSASZ_pap_eff_lfp_rep(cfg)
    
    % Analysis
    iSASZ_pap_ana_hgp_rep(cfg)    

    iSASZ_pap_ana_lfp_rep(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = { ''};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld(1:32);
cfg.ana_nme = { ''};
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
