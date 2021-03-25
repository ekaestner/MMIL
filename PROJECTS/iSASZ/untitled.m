clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

sbj_tme = [];

for sbj_num = 15:23;
    
    tic
    
    %%
    pcfg = [];
    pcfg.prj_dat_hld = sbj_dat_hld{sbj_num};
    
    pcfg.eff_typ = {'hgp'};
    pcfg.eff_nme = 'Selectivity';
    pcfg.eff_rul = {'all'};
    
    pcfg.toi = [-0.200 1.200];
    pcfg.bse_lne = [-0.200 0];
    
    
    pcfg.alt_eve = 'yes';
    pcfg.eve = [111 112 211 212];
    pcfg.eve_nme = {'VisualNew' 'VisualOld' 'AuditoryNew' 'AuditoryOld'};
    
    
    pcfg.col_bar = [0 0.6];
    
    pcfg.win = 0.050;
    
    pcfg.foi = [4 8];
    
    pcfg.mve_win = 0.010;
    pcfg.mve_tme = [0.100 1.125];
    
    pcfg.sbj_num     = sbj_num;
    
    PLV_Play3(pcfg)
    
    %% All Plot
    
    
    
    %% Focused Plot
    
    
    sbj_tme(sbj_num) = toc;
    
end