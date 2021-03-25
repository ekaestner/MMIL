function iSASZ_connectivity

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

sbj_tme = [];

for sbj_num = [19 20]; %18:numel(sbj_nme_hld);
        
        tic
        
        pcfg = [];
        
        pcfg.eff_typ = {''};
        pcfg.eff_nme = '';
        pcfg.eff_rul = {'all'};
        
        pcfg.toi = [-0.100 0.600];
        pcfg.bse_lne = [-0.100 0];
        pcfg.alt_eve = [];
        pcfg.eve     = [111         112         211           212];
        pcfg.eve_nme = {'VisualNew' 'VisualOld' 'AuditoryNew' 'AuditoryOld' };
        
        pcfg.win = 0.050;
        
        pcfg.ffp = [4 24];
        pcfg.foi = [4 8];
        
        pcfg.sig_nme = {'pap_vis_act' 'pap_aud_act'};
        pcfg.sig_typ = 'hgp';
        pcfg.sig_col = {[1] [1]};
        
        pcfg.mve_win = 0.010;
        pcfg.mve_tme = [0.050 0.600];
        
        pcfg.sbj_num     = sbj_num;
        pcfg.prj_dat_hld = prj_dat_hld;
                
        PLV_Play3(pcfg)
        
        sbj_tme(sbj_num) = toc;
                
end

end