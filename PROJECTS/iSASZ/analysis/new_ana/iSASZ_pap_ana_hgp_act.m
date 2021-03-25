function iSASZ_pap_ana_hgp_act(fcfg)

fprintf([fcfg.sbj_nme ': Starting iSASZ_pap_ana_hgp_act on %s \n'],fcfg.sbj_nme)

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
    
    cfg.cmb_nme = { 'pap_vis_act' };
    
    cfg.sig_loc = { { 'vis_new_ovr_stt' 'vis_new_ovr_stt' } };
    
    cfg.sig_chn = { [ 1 2 ] };
    
    % Table Combination
    cfg.run_tbl     = [1];
    
    cfg.tbl_cmb_nme = { { 'VisualResponsive' 'TotalVisualResponsive' } };
    
    cfg.tbl_cmb_rul = { { '1' '2' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'VisualResponsive' 'LateVisualResponsive' } };
    cfg.plt_cmb_col = { { 'reddish grey' 'dark red' } };
    
    cfg.plt_pct_cmb_nme = { { 'Visual Responsive' 'LateVisualResponsive' } };
    cfg.plt_pct_cmb_col = { { 'reddish grey' 'dark red' } };
    
    cfg.plt_cmb_rul     = { { '1' '2'} };
    cfg.plt_prb_cmb_rul = { { '1' '2'} };
    
    cfg.plt_prb_cmb_max = { [ 1 2 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 1 ];
    
    cfg.fst_cmb_nme = { { 'VisualResponsive' } };
    
    cfg.fst_cmb_col = { { 'red' } };
    
    cfg.inc_cmb_rul = { { '2' } };
    cfg.fst_cmb_rul = { { '2' } };
    
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
    
    cfg.alt_eve = { { 'vis_new_bse' } };
    cfg.eve     = { { [102] } };
    cfg.lnstyle.col_ord = { { { rgb('reddish grey') } } };
    
    cfg.stt_dat = { { { 'vis_new_ovr_stt'} } };
    cfg.stt_col = { { { ft_stt_col(rgb('red')) } } };
    cfg.stt_cmp = { { { '0%5'              } } };
        
    cfg.v_lne       = { [0 0.2 0.4] };
    cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} };
    
    cfg.dat_use = 2;
    
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
    
    cfg.cmb_nme = { 'pap_aud_act' };
    
    cfg.sig_loc = { { 'aud_new_ovr_stt' 'aud_new_ovr_stt' } };
    
    cfg.sig_chn = { [ 1 2 ] };
    
    % Table Combination
    cfg.run_tbl     = [1];
    
    cfg.tbl_cmb_nme = { { 'AuditoryResponsive' 'AuditoryVisualResponsive' } };
    
    cfg.tbl_cmb_rul = { { '1' '2' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'AuditoryResponsive' 'LateAuditoryResponsive' } };
    cfg.plt_cmb_col = { { 'bluish grey' 'dark blue'} };
    
    cfg.plt_pct_cmb_nme = { { 'AuditoryResponsive' 'LateAuditoryResponsive' } };
    cfg.plt_pct_cmb_col = { { 'bluish grey' 'dark blue' } };
    
    cfg.plt_cmb_rul     = { { '1' '2' } };
    cfg.plt_prb_cmb_rul = { { '1' '2' } };
    
    cfg.plt_prb_cmb_max = { [ 1 2 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 1 ];
    
    cfg.fst_cmb_nme = { { 'AuditoryResponsive' } };
    
    cfg.fst_cmb_col = { { 'bluish grey' } };
    
    cfg.inc_cmb_rul = { { '2' } };
    cfg.fst_cmb_rul = { { '2' } };
    
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
    
    cfg.alt_eve = { { 'aud_new_bse' } };
    cfg.eve     = { { [202] } };
    cfg.lnstyle.col_ord = { { { rgb('bluish grey') } } };
    
    cfg.stt_dat = { { {'aud_new_ovr_stt'} } };
    cfg.stt_col = { { { ft_stt_col(rgb('blue')) } } };
    cfg.stt_cmp = { { { '0%5'                } } };
    
    cfg.v_lne       = { [0 0.2 0.4] };
    cfg.v_lne_col   = { {rgb('blue') rgb('black') rgb('black')} };
    
    cfg.dat_use = 2;
    
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
    
    cfg.cmb_nme = { 'pap_bim_act' };
    
    cfg.sig_loc = { { 'vis_new_ovr_stt' 'vis_new_ovr_stt' 'aud_new_ovr_stt' 'aud_new_ovr_stt' } };
    
    cfg.sig_chn = { [ 1 2 1 2] };
    
    % Table Combination
    cfg.run_tbl     = [1];
    
    cfg.tbl_cmb_nme = { { 'BimodalResponsive' 'LateBimodalResponsive' } };
    
    cfg.tbl_cmb_rul = { { '1 & 3' '2 & 4' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'BimodalResponsive' 'LateBimodalResponsive' } };
    
    cfg.plt_cmb_col = { { 'purple grey' 'dark purple' } };
    
    cfg.plt_pct_cmb_nme = { { 'BimodalResponsive' 'LateBimodalResponsive' } };
    cfg.plt_pct_cmb_col = { { 'purple grey' 'dark purple' } };
    
    cfg.plt_cmb_rul     = { { '1 & 3' '2 & 4' } };
    cfg.plt_prb_cmb_rul = { { '1 & 3' '2 & 4' } };
    
    cfg.plt_prb_cmb_max = { [ 1 2 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 0 ];
    
    % Time combinations
    cfg.run_tme     = [0];
    
    % Line Plot combinations
    cfg.run_lne     = [1];
    
    cfg.plt_dim = [ 2 2 ];
    cfg.dat_loc = [ 1 1 ; 0 0];
    
    cfg.alt_lab = { 'stt_lab' };
    
    cfg.alt_eve = { { 'vis_new_bse' 'aud_new_bse' ; '' '' } };
    cfg.eve     = { { [102] [202] ; [] [] } };
    cfg.lnstyle.col_ord = { { { rgb('reddish grey')} {rgb('bluish grey')} ; {} {} } };
    
    cfg.stt_dat = { { { 'vis_new_ovr_stt'}                   { 'aud_new_ovr_stt' }                  ; {} {} } };
    cfg.stt_col = { { { ft_stt_col(rgb('red')) } { ft_stt_col(rgb('blue')) } ; {} {} } };
    cfg.stt_cmp = { { { '0%5' }                           { '0%5' }                         ; {} {} } };
    
    cfg.v_lne       = { [0 0.2 0.4] };
    cfg.v_lne_col   = { {rgb('blue') rgb('black') rgb('black')} };
    
    cfg.dat_use = 2;
    
    % Run
    cfg.dat_nme     = '_overall_data';
    cfg.sbj_clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    
    cfg.fle_out_pth = fcfg.dat_fld;
    
    mmil_combo_effects2(cfg);
    
    % Save for overall analysis
    for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end
        
end