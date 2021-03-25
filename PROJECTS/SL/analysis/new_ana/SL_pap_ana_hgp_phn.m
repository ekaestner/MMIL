function SL_pap_ana_hgp_mtc(fcfg)

%%
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_phn_950' };

cfg.sig_loc = { { 'vis_con_stt_msk_anv' 'aud_con_stt_msk_anv' 'pap_ovr_stt' 'pap_anv_stt' 'vis_con_stt' 'aud_con_stt' } };

cfg.sig_chn = { [ 3 3 4 4 3 3 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'VisualLetter' 'AuditoryPhoneme' 'LetterCorrected' 'PhonemeCorrected' } };

cfg.tbl_cmb_rul = { { '1 & 3 & 4' '2 & 3 & 4' '( 1 & ~ ( 3 & 4 ) ) | ( 5 & ~ 1 )' '( 2 & ~ ( 3 & 4 ) ) | ( 6 & ~ 2 )' } };

cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'VisualLetter' 'AuditoryPhoneme' } };
cfg.plt_cmb_col = { { 'bright red'   'bright blue'     } };

cfg.plt_pct_cmb_nme = { { 'VisualLetter' 'AuditoryPhoneme' 'LetterCorrected' 'PhonemeCorrected' } };
cfg.plt_pct_cmb_col = { { 'bright red'   'bright blue'     'black' 'black'} };

cfg.plt_cmb_rul     = { { '1 & 3 & 4' '2 & 3 & 4' } };
cfg.plt_prb_cmb_rul = { { '1 & 3 & 4' '2 & 3 & 4' '( 1 & ~ ( 3 & 4 ) ) | ( 5 & ~ 1 )' '( 2 & ~ ( 3 & 4 ) ) | ( 6 & ~ 2 )'  } };

cfg.plt_prb_cmb_max = { [ 1 1 2 2] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'VisualLetter' 'AuditoryPhoneme' } };

cfg.fst_cmb_col = { { 'bright red'   'bright blue' } };

cfg.inc_cmb_rul = { { '1 & 3 & 4' '2 & 3 & 4' } };
cfg.fst_cmb_rul = { { '1'         '2'         } };

cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };
cfg.fst_ord_nme = {{}};

% Time combinations
cfg.run_tme     = [0];

% Line Plot combinations
phn_col = distinguishable_colors(24);

cfg.run_lne     = [1];

cfg.dat_use = 2;

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 1 ; 1 1];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'vis_tot_nse' 'aud_tot_nse' ; 'vis_con' 'aud_con'} };
cfg.eve     = { { [101 102]     [111 112]     ; [1:12]    [1:12]} };
cfg.lnstyle.col_ord = { { { rgb('red') rgb('reddish grey') } {rgb('blue') rgb('bluish grey')} ; mat2cell(phn_col(1:12,:),ones(1,size(phn_col(1:12,:),1)),3)' mat2cell(phn_col(13:24,:),ones(1,size(phn_col(13:24,:),1)),3)' } };

cfg.stt_dat = { { {'vis_nse_stt'}                     {'aud_nse_stt'}                    ; { 'vis_con_stt_msk_anv' 'vis_con_stt' }                  { 'aud_con_stt_msk_anv' 'aud_con_stt' } } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) } { ft_stt_col(rgb('bluish grey')) } ; { ft_stt_col(rgb('red')) ft_stt_col(rgb('bright red')) } {ft_stt_col(rgb('blue')) ft_stt_col(rgb('bright blue'))} } };
cfg.stt_cmp = { { { '0%5' }                           { '0%5' }                          ; { '0%3' '3%6' }                                          { '0%3' '3%6' } } };

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