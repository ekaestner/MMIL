function iSASZ_pap_ana_hgp_sem(fcfg)

fprintf([fcfg.sbj_nme ': Starting iSASZ_pap_ana_hgp_sem on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

%% Visual
if any(eve_typ < 9)
    
    cfg = [];
    
    cfg.typ      = { 'hgp' };
    cfg.ele_typ  = { 'ecog' };
    cfg.loc_typ  = { 'split' };
    
    cfg.cmb_nme = { 'pap_vis_sem' };
    
    cfg.sig_loc = { { 'vis_new_ovr_stt' 'vis_new_ovr_stt' 'vis_ani_obj_stt_msk_ovr' 'vis_ani_obj_stt_msk_ovr' 'vis_ani_obj_stt_msk_ovr' 'vis_ani_obj_stt_msk_ovr' } };
    
    cfg.sig_chn = { [ 1 2 1 2 3 4 ] };
    
    % Table Combination
    cfg.run_tbl     = [ 1 ];
    
    cfg.tbl_cmb_nme = { { 'VisualAnimal' 'VisualObject' 'TotalVisualAnimal' 'TotalVisualObject' 'TotalVisualObjectCorrected' } };
        
    cfg.tbl_cmb_rul = { { '1 & 3' '1 & 4' '2 & 5' '2 & 6' '~ ( 1 & 2 ) & ( 3 | 4 | 5 | 6 )' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'VisualAnimal' 'VisualObject' 'TotalVisualAnimal' 'TotalVisualObject' } };
    cfg.plt_cmb_col = { { 'Green'        'Grey'         'Dark Green'        'Dark Grey' } };
    
    cfg.plt_pct_cmb_nme = { { 'VisualAnimal' 'VisualObject' 'TotalVisualAnimal' 'TotalVisualObject' 'TotalVisualObjectCorrected' } };
    cfg.plt_pct_cmb_col = { { 'Green'        'Grey'         'Dark Green'        'Dark Grey'         'black'} };
    
    cfg.plt_cmb_rul     = { { '1 & 3' '2 & 4' '1 & 5' '2 & 6' } };
    cfg.plt_prb_cmb_rul = { { '1 & 3' '2 & 4' '1 & 5' '2 & 6' '~ ( 1 & 2 ) & ( 3 | 4 | 5 | 6 )' } };
    
    cfg.plt_prb_cmb_max = { [ 1 1 2 2 3 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 1 ];
    
    cfg.fst_cmb_nme = { { 'VisualAnimal' 'VisualObject' } };
    
    cfg.fst_cmb_col = { { 'green'        'grey' } };
    
    cfg.inc_cmb_rul = { { '2 & 5'        '2 & 6' } };
    cfg.fst_cmb_rul = { { '5'            '6' } };
    
    cfg.fst_ord     = {{}};
    cfg.fst_ord_dim = {{'avg'} };
    cfg.fst_ord_nme = {{}};
    
    % Time combinations
    cfg.run_tme     = [0];
    
    % Line Plot combinations
    cfg.run_lne     = [1];
    
    cfg.plt_dim = [ 1 1 ];
    cfg.dat_loc = [ 1 ];
    
    cfg.alt_lab = { 'stt_lab' };
    
    cfg.alt_eve = { { 'vis_ani_obj' } };
    cfg.eve     = { { [121 122] } };
    cfg.lnstyle.col_ord = { { { rgb('green') rgb('reddish gray') } } };
    
    cfg.stt_dat = { { {'vis_new_ovr_stt' 'vis_ani_obj_stt_msk_ovr' } } };
    cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('green'))} } };
    cfg.stt_cmp = { { { '0%3'   '3%6'} } };
        
    cfg.v_lne       = { [0 0.2 0.4] };
    cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} };
    
    % Run
    cfg.dat_nme     = '_overall_data';
    cfg.sbj_clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    
    cfg.fle_out_pth = fcfg.dat_fld;
    
    mmil_combo_effects2(cfg);
    
    % Save for overall analysis
    for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end
    
end

%% Auditory
if any(eve_typ > 9)
       
    cfg = [];
    
    cfg.typ      = { 'hgp' };
    cfg.ele_typ  = { 'ecog' };
    cfg.loc_typ  = { 'split' };
    
    cfg.cmb_nme = { 'pap_aud_sem' };
    
    cfg.sig_loc = { { 'aud_new_ovr_stt' 'aud_new_ovr_stt' 'aud_ani_obj_stt_msk_ovr' 'aud_ani_obj_stt_msk_ovr' 'aud_ani_obj_stt_msk_ovr' 'aud_ani_obj_stt_msk_ovr' } };
    
    cfg.sig_chn = { [ 1 2 1 2 3 4 ] };
    
    % Table Combination
    cfg.run_tbl     = [ 1 ];
    
    cfg.tbl_cmb_nme = { { 'AuditoryAnimal' 'AuditoryObject' 'TotalAuditoryAnimal' 'TotalAuditoryObject' 'TotalAuditoryObjectCorrected' } };
        
    cfg.tbl_cmb_rul = { { '1 & 3' '1 & 4' '2 & 5' '2 & 6' '~ ( 1 & 2 ) & ( 3 | 4 | 5 | 6 )' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'AuditoryAnimal' 'AuditoryObject' 'TotalAuditoryAnimal' 'TotalAuditoryObject' } };
    cfg.plt_cmb_col = { { 'Green'        'Grey'         'Dark Green'        'Dark Grey' } };
    
    cfg.plt_pct_cmb_nme = { { 'AuditoryAnimal' 'AuditoryObject' 'TotalAuditoryAnimal' 'TotalAuditoryObject' 'TotalAuditoryObjectCorrected' } };
    cfg.plt_pct_cmb_col = { { 'Green'        'Grey'         'Dark Green'        'Dark Grey'         'black'} };
    
    cfg.plt_cmb_rul     = { { '1 & 3' '2 & 4' '1 & 5' '2 & 6' } };
    cfg.plt_prb_cmb_rul = { { '1 & 3' '2 & 4' '1 & 5' '2 & 6' '~ ( 1 & 2 ) & ( 3 | 4 | 5 | 6 )' } };
    
    cfg.plt_prb_cmb_max = { [ 1 1 2 2 3 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 1 ];
    
    cfg.fst_cmb_nme = { { 'AuditoryAnimal' 'AuditoryObject' } };
    
    cfg.fst_cmb_col = { { 'green'        'grey' } };
    
    cfg.inc_cmb_rul = { { '2 & 5'        '2 & 6' } };
    cfg.fst_cmb_rul = { { '5'            '6' } };
    
    cfg.fst_ord     = {{}};
    cfg.fst_ord_dim = {{'avg'} };
    cfg.fst_ord_nme = {{}};
    
    % Time combinations
    cfg.run_tme     = [0];
    
    % Line Plot combinations
    cfg.run_lne     = [1];
    
    cfg.plt_dim = [ 1 1 ];
    cfg.dat_loc = [ 1 ];
    
    cfg.alt_lab = { 'stt_lab' };
    
    cfg.alt_eve = { { 'aud_ani_obj' } };
    cfg.eve     = { { [221 222] } };
    cfg.lnstyle.col_ord = { { { rgb('bluish green') rgb('bluish grey') } } };
    
    cfg.stt_dat = { { {'aud_new_ovr_stt' 'aud_ani_obj_stt_msk_ovr' } } };
    cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('green'))} } };
    cfg.stt_cmp = { { { '0%3'   '3%6'} } };
        
    cfg.v_lne       = { [0 0.2 0.4] };
    cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} };
    
    % Run
    cfg.dat_nme     = '_overall_data';
    cfg.sbj_clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    
    cfg.fle_out_pth = fcfg.dat_fld;
    
    mmil_combo_effects2(cfg);
    
    % Save for overall analysis
    for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end
    
end

%% BiModal
if any(eve_typ < 9) && any(eve_typ > 9)
    
    cfg = [];
    
    cfg.typ      = { 'hgp' };
    cfg.ele_typ  = { 'ecog' };
    cfg.loc_typ  = { 'split' };
    
    cfg.cmb_nme = { 'pap_bim_sem' };
    
    cfg.sig_loc = { { 'vis_new_ovr_stt' 'vis_new_ovr_stt' 'vis_ani_obj_stt_msk_ovr' 'vis_ani_obj_stt_msk_ovr' 'vis_ani_obj_stt_msk_ovr' 'vis_ani_obj_stt_msk_ovr' 'aud_new_ovr_stt' 'aud_new_ovr_stt' 'aud_ani_obj_stt_msk_ovr' 'aud_ani_obj_stt_msk_ovr' 'aud_ani_obj_stt_msk_ovr' 'aud_ani_obj_stt_msk_ovr' } };
    
    cfg.sig_chn = { [ 1 2 1 2 3 4 1 2 1 2 3 4] };
    
    % Table Combination
    cfg.run_tbl     = [ 1 ];
    
    cfg.tbl_cmb_nme = { { 'BimodalAnimal' 'BimodalObject' 'TotalBimodalAnimal' 'TotalBimodalObject' 'BimodalCorrected' } };
        
    cfg.tbl_cmb_rul = { { '1 & 3 & 7 & 9' '1 & 4 & 7 & 10' '2 & 5 & 8 & 11' '2 & 6 & 8 & 12' '( 5 | 6 | 11 | 12 ) & ~ ( 2 & 8 )' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'BimodalAnimal' 'BimodalObject' 'TotalBimodalAnimal' 'TotalBimodalObject' } };
    cfg.plt_cmb_col = { { 'green'         'grey'          'dark green'          'dark grey'         } };
    
    cfg.plt_pct_cmb_nme = { { 'BimodalAnimal' 'BimodalObject' 'TotalBimodalAnimal' 'TotalBimodalObject' 'BimodalCorrected' } };
    cfg.plt_pct_cmb_col = { { 'green'         'grey'          'dark green'          'dark grey'         'black'} };
    
    cfg.plt_cmb_rul     = { { '1 & 3 & 7 & 9' '1 & 4 & 7 & 10' '2 & 5 & 8 & 11' '2 & 6 & 8 & 12' } };
    cfg.plt_prb_cmb_rul = { { '1 & 3 & 7 & 9' '1 & 4 & 7 & 10' '2 & 5 & 8 & 11' '2 & 6 & 8 & 12' '( 5 | 6 | 11 | 12 ) & ~ ( 2 & 8 )' } };
    
    cfg.plt_prb_cmb_max = { [ 1 1 2 2 3 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 0 ];
    
    
    % Time combinations
    cfg.run_tme     = [0];
    
    % Line Plot combinations
    cfg.run_lne     = [1];
    
    cfg.plt_dim = [ 2 2 ];
    cfg.dat_loc = [ 1 1 ; 0 0];
    
    cfg.alt_lab = { 'stt_lab' };
    
    cfg.alt_eve = { { 'vis_ani_obj' 'aud_ani_obj' ; '' ''} };
    cfg.eve     = { { [121 122] [221 222] ; [] [] } };
    cfg.lnstyle.col_ord = { { { rgb('green') rgb('reddish gray') } { rgb('bluish green') rgb('bluish grey') } ; {} {} } };
    
    cfg.stt_dat = { { {'vis_new_ovr_stt' 'vis_ani_obj_stt_msk_ovr' } {'aud_new_ovr_stt' 'aud_ani_obj_stt_msk_ovr'} ; {} {}} };
    cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('green'))} { ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('bluish green')) }; {} {} } };
    cfg.stt_cmp = { { { '0%3'   '6%9'} { '3%6'            '9%12'} ; {} {}} };
        
    cfg.v_lne       = { [0 0.2 0.4] };
    cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} };
    
    % Run
    cfg.dat_nme     = '_overall_data';
    cfg.sbj_clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    
    cfg.fle_out_pth = fcfg.dat_fld;
    
    mmil_combo_effects2(cfg);
    
    % Save for overall analysis
    for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end
        
end