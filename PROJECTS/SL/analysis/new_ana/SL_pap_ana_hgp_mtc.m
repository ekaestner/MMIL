function SL_pap_ana_hgp_mtc(fcfg)

%%
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_mtc_1450' };

cfg.sig_loc = { { 'vis_mtc_stt_msk_anv' 'vis_mtc_stt_msk_anv' 'vis_mtc_stt_msk_anv' 'pap_ovr_stt' 'pap_anv_stt' 'vis_mtc_stt' } };

cfg.sig_chn = { [ 1 2 3 4 4 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Mismatch_sensory' 'Mismatch_post' 'Mismatch_total' 'Corrected' } };

cfg.tbl_cmb_rul = { { '1 & 4 & 5' '2 & 4 & 5' '3 & 4 & 5' '( 6 & 4 & 5 ) | ( ( 1 | 2 | 3 ) & ~ ( 4 & 5 ) )' } };

cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Mismatch_sensory' 'Mismatch_post' 'Mismatch_total' } };
cfg.plt_cmb_col = { { 'bright yellow'    'dark yellow'   'yellow' } };

cfg.plt_pct_cmb_nme = { { 'Mismatch_sensory' 'Mismatch_post' 'Mismatch_total' 'Corrected' } };
cfg.plt_pct_cmb_col = { { 'bright yellow'    'dark yellow'   'yellow'         'grey' } };

cfg.plt_cmb_rul     = { { '1 & 4 & 5' '2 & 4 & 5' '3 & 4 & 5' } };
cfg.plt_prb_cmb_rul = { { '1 & 4 & 5' '2 & 4 & 5' '3 & 4 & 5' '( 6 & 4 & 5 ) | ( ( 1 | 2 | 3 ) & ~ ( 4 & 5 ) )'  } };

cfg.plt_prb_cmb_max = { [ 1 1 1 2 ] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Mismatch_sensory' 'Mismatch_post' 'Mismatch_total' } };

cfg.fst_cmb_col = { { 'bright yellow'    'dark yellow'   'yellow' } };

cfg.inc_cmb_rul = { { '1 & 4 & 5' '2 & 4 & 5' '3 & 4 & 5' } };
cfg.fst_cmb_rul = { { '1'         '2'         '3' } };

cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };
cfg.fst_ord_nme = {{}};

% Time combinations
cfg.run_tme     = [0];

% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 1 ; 1 0];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'vis_tot_nse' 'aud_tot_nse' ; 'trialinfo' ''} };
cfg.eve     = { { [101 102] [111 112] ; [1 2] []} };
cfg.lnstyle.col_ord = { { { rgb('red') rgb('reddish grey') } {rgb('blue') rgb('bluish grey')} ; {rgb('green') rgb('yellow')} {} } };

cfg.stt_dat = { { {'vis_nse_stt'}                     {'aud_nse_stt'}                    ; {'vis_mtc_stt_msk_anv' 'vis_mtc_stt'}                          {''} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) } { ft_stt_col(rgb('bluish grey')) } ; { ft_stt_col(rgb('yellow')) ft_stt_col(rgb('bright yellow')) } {} } };
cfg.stt_cmp = { { { '0%5' }                           { '0%5' }                          ; { '0%3'                     '3%6' }                            {} } };

cfg.v_lne       = { [0 0.450 0.900] };
cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };

% Run
cfg.dat_nme     = '_overall_data';
cfg.sbj_clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;

cfg.fle_out_pth = fcfg.dat_fld;

mmil_combo_effects2(cfg);

% Save for overall analysis
for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

end