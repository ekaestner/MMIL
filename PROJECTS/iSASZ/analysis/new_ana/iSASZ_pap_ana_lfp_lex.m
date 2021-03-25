function iSASZ_pap_ana_lfp_lex(fcfg)

fprintf([fcfg.sbj_nme ': Starting iSASZ_pap_ana_lfp_lex on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

%% Visual
if any(eve_typ < 9)
    
    cfg = [];
    
    cfg.typ      = { 'lfp' };
    cfg.ele_typ  = { 'ecog' };
    cfg.loc_typ  = { 'split' };
    
    cfg.cmb_nme = { 'pap_vis_lex' };
    
    cfg.sig_loc = { { 'vis_new_ovr_stt' 'vis_new_ovr_stt' 'vis_frq_stt_msk_ovr' 'vis_frq_stt_msk_ovr'} };
    
    cfg.sig_chn = { [ 1 2 1 2 ] };
    
    % Table Combination
    cfg.run_tbl     = [ 1 ];
    
    cfg.tbl_cmb_nme = { { 'VisualLexical' 'TotalVisualLexical' 'LexicalCorrected'} };
        
    cfg.tbl_cmb_rul = { { '1 & 3' '2 & 4' '~ ( 1 & 2 ) & ( 3 | 4 )' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'VisualLexical' 'TotalVisualLexical' } };
    cfg.plt_cmb_col = { { 'dark yellow' 'yellowish tan' } };
    
    cfg.plt_pct_cmb_nme = { { 'VisualLexical' 'TotalVisualLexical' 'LexicalCorrected'} };
    cfg.plt_pct_cmb_col = { { 'dark yellow' 'yellowish tan' 'black'} };
    
    cfg.plt_cmb_rul     = { { '1 & 3' '2 & 4' } };
    cfg.plt_prb_cmb_rul = { { '1 & 3' '2 & 4' '~ ( 1 & 2 ) & ( 3 | 4 )' } };
    
    cfg.plt_prb_cmb_max = { [ 1 1 2 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 1 ];
    
    cfg.fst_cmb_nme = { { 'VisualResponsive' 'VisualLexical' } };
    
    cfg.fst_cmb_col = { { 'yellowish tan' } };
    
    cfg.inc_cmb_rul = { { '2 & 4' } };
    cfg.fst_cmb_rul = { { '4' } };
    
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
    
    cfg.alt_eve = { { 'vis_frq_med' } };
    cfg.eve     = { { [151 152] } };
    cfg.lnstyle.col_ord = { { { rgb('yellowish tan')-[0.4 0.4 .2]  rgb('yellowish tan') } } };
    
    cfg.stt_dat = { { {'vis_new_ovr_stt' 'vis_frq_stt_msk_ovr' } } };
    cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('yellowish tan')) } } };
    cfg.stt_cmp = { { { '0%3'               '3%6'} } };
    
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
    
    cfg.typ      = { 'lfp' };
    cfg.ele_typ  = { 'ecog' };
    cfg.loc_typ  = { 'split' };
    
    cfg.cmb_nme = { 'pap_aud_lex' };
    
    cfg.sig_loc = { { 'aud_new_ovr_stt' 'aud_new_ovr_stt' 'aud_frq_stt_msk_ovr' 'aud_frq_stt_msk_ovr'} };
    
    cfg.sig_chn = { [ 1 2 1 2 ] };
    
    % Table Combination
    cfg.run_tbl     = [ 1 ];
    
    cfg.tbl_cmb_nme = { { 'AuditoryLexical' 'TotalAuditoryLexical' 'LexicalCorrected'} };
        
    cfg.tbl_cmb_rul = { { '1 & 3' '2 & 4' '~ ( 1 & 2 ) & ( 3 | 4 )' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'AuditoryLexical' 'TotalAuditoryLexical' } };
    cfg.plt_cmb_col = { { 'dark green' 'bluish green' } };
    
    cfg.plt_pct_cmb_nme = { { 'AuditoryLexical' 'TotalAuditoryLexical' 'LexicalCorrected'} };
    cfg.plt_pct_cmb_col = { { 'dark green' 'bluish green' 'black'} };
    
    cfg.plt_cmb_rul     = { { '1 & 3' '2 & 4' } };
    cfg.plt_prb_cmb_rul = { { '1 & 3' '2 & 4' '~ ( 1 & 2 ) & ( 3 | 4 )' } };
    
    cfg.plt_prb_cmb_max = { [ 1 2 3 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 1 ];
    
    cfg.fst_cmb_nme = { { 'AuditoryLexical' } };
    
    cfg.fst_cmb_col = { { 'bluish green' } };
    
    cfg.inc_cmb_rul = { { '2 & 4' } };
    cfg.fst_cmb_rul = { { '4' } };
    
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
    
    cfg.alt_eve = { { 'aud_frq_med' } };
    cfg.eve     = { { [251 252] } };
    cfg.lnstyle.col_ord = { { { rgb('bluish green')-[0 0.4 0.25]  rgb('bluish green') } } };
    
    cfg.stt_dat = { { {'aud_new_ovr_stt' 'aud_frq_stt_msk_ovr' } } };
    cfg.stt_col = { { { ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('bluish green')) } } };
    cfg.stt_cmp = { { { '0%3'               '3%6'} } };
        
    cfg.v_lne       = { [0 0.2 0.4] };
    cfg.v_lne_col   = { {rgb('blue') rgb('black') rgb('black')} };
    
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
    
    cfg.typ      = { 'lfp' };
    cfg.ele_typ  = { 'ecog' };
    cfg.loc_typ  = { 'split' };
    
    cfg.cmb_nme = { 'pap_bim_lex' };
    
    cfg.sig_loc = { { 'vis_new_ovr_stt' 'vis_new_ovr_stt' 'aud_new_ovr_stt' 'aud_new_ovr_stt' 'vis_frq_stt_msk_ovr' 'vis_frq_stt_msk_ovr' 'aud_frq_stt_msk_ovr' 'aud_frq_stt_msk_ovr' } };
    
    cfg.sig_chn = { [ 1 2 1 2 1 2 1 2 ] };
    
    % Table Combination
    cfg.run_tbl     = [ 1 ];
    
    cfg.tbl_cmb_nme = { { 'BimodalLexical' 'TotalBimodalLexical' 'BimodalCorrected' } };
        
    cfg.tbl_cmb_rul = { { '1 & 5 & 3 & 7' '2 & 6 & 4 & 8' '( ( 5 | 6 ) & ~ ( 1 & 2 ) ) & ( ( 7 | 8 ) & ~ ( 3 & 4 ) )' } };
    
    cfg.tbl_sum = { [] };
    
    % Plot combination
    cfg.run_plt     = [1];
    
    cfg.plt_cmb_nme = { { 'BimodalLexical' 'TotalBimodalLexical'  } };
    cfg.plt_cmb_col = { { 'yellow green' 'brownish yellow' } };
    
    cfg.plt_pct_cmb_nme = { { 'BimodalLexical' 'TotalBimodalLexical' 'BimodalCorrected' } };
    cfg.plt_pct_cmb_col = { { 'yellow green' 'brownish yellow' 'black'} };
    
    cfg.plt_cmb_rul     = { { '1 & 5 & 3 & 7' '2 & 6 & 4 & 8' } };
    cfg.plt_prb_cmb_rul = { { '1 & 5 & 3 & 7' '2 & 6 & 4 & 8' '( ( 5 | 6 ) & ~ ( 1 & 2 ) ) & ( ( 7 | 8 ) & ~ ( 3 & 4 ) )' } };
    
    cfg.plt_prb_cmb_max = { [ 1 2 3 ] };
    
    % First Pass combinations
    cfg.run_fst     = [ 0 ];
    
    
    % Time combinations
    cfg.run_tme     = [0];
    
    % Line Plot combinations
    cfg.run_lne     = [1];
    
    cfg.plt_dim = [ 2 2 ];
    cfg.dat_loc = [ 1 1 ; 0 0];
    
    cfg.alt_lab = { 'stt_lab' };
    
    cfg.alt_eve = { { 'vis_frq_med' 'aud_frq_med' ; '' ''} };
    cfg.eve     = { { [151 152] [251 252] ; [] [] } };
    cfg.lnstyle.col_ord = { { { rgb('yellowish tan')-[0.4 0.4 .2]  rgb('yellowish tan') } { rgb('bluish green')-[0 0.4 0.25]  rgb('bluish green') } ; {} {} } };
    
    cfg.stt_dat = { { {'vis_new_ovr_stt' 'vis_frq_stt_msk_ovr' } {'aud_new_ovr_stt' 'aud_frq_stt_msk_ovr'} ; {} {}} };
    cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('yellowish tan'))} { ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('bluish green')) }; {} {} } };
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