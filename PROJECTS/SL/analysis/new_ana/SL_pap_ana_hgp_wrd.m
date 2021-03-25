function SL_pap_ana_hgp_wrd(fcfg)

%%
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_wrd_950' };

cfg.sig_loc = { { 'ovr_stt' 'ovr_stt' 'pap_anv_stt' 'pap_anv_stt' 'vis_wrd_stt_msk_anv' 'vis_wrd_stt_msk_anv' 'aud_wrd_stt_msk_anv' 'aud_wrd_stt_msk_anv' } };

cfg.sig_chn = { [ 1 2 1 2 1 2 1 2 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'VisualWord' 'VisualNonWord' 'AuditoryWord' 'AuditoryNonWord' } };

cfg.tbl_cmb_rul = { { '1 & 3 & 5'  '1 & 3 & 6'     '2 & 4 & 7'    '2 & 4 & 8' } };

cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'VisualWord' 'VisualNonWord' 'AuditoryWord' 'AuditoryNonWord' } };
cfg.plt_cmb_col = { { 'bright red' 'dark red'      'bright blue'  'dark blue' } };

cfg.plt_pct_cmb_nme = { { 'VisualWord' 'VisualNonWord' 'AuditoryWord' 'AuditoryNonWord' } };
cfg.plt_pct_cmb_col = { { 'bright red' 'dark red'      'bright blue'  'dark blue' } };

cfg.plt_cmb_rul     = { { '1 & 3 & 5'  '1 & 3 & 6'     '2 & 4 & 7'    '2 & 4 & 8' } };
cfg.plt_prb_cmb_rul = { { '1 & 3 & 5'  '1 & 3 & 6'     '2 & 4 & 7'    '2 & 4 & 8' } };

cfg.plt_prb_cmb_max = { [ 1 1 2 2 ] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'VisualWord' 'VisualNonWord' 'AuditoryWord' 'AuditoryNonWord' } };

cfg.fst_cmb_col = { { 'bright red' 'dark red'      'bright blue'  'dark blue' } };

cfg.inc_cmb_rul = { { '1 & 3 & 5'  '1 & 3 & 6'     '2 & 4 & 7'    '2 & 4 & 8' } };
cfg.fst_cmb_rul = { { '5'          '6'             '7'            '8' } };

cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };
cfg.fst_ord_nme = {{}};

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { '' } };
               
cfg.tme_cmb_col = { { '' } };
               
cfg.tme_cmb_rul = { { '' } };

cfg.tme_ord     = { { [] } };
cfg.tme_ord_nme = { { '' } };  

% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 1 ; 1 0];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'vis_wrd' 'aud_wrd' ; 'trialinfo' ''} };
cfg.eve     = { { [101 102] [201 202] ; [1 2 3 4] []} };
cfg.lnstyle.col_ord = { { { rgb('bright red') rgb('dark red') } {rgb('bright blue') rgb('dark blue')} ; {rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey')} {} } };
                     
cfg.stt_dat = { { {'ovr_stt' 'vis_anv_stt' 'vis_wrd_stt_msk_anv'}                              {'ovr_stt' 'vis_anv_stt' 'aud_wrd_stt_msk_anv'}                              ; {'vis_nse_stt'                    'aud_nse_stt'                  'vis_mtc_stt_msk_anv'} { ''} } };
cfg.stt_col = { { { ft_stt_col(rgb('black')) ft_stt_col(rgb('green')) ft_stt_col(rgb('red')) } { ft_stt_col(rgb('black')) ft_stt_col(rgb('green')) ft_stt_col(rgb('blue'))} ; { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('bluish grey')) ft_stt_col(rgb('yellow')) } {} } };
cfg.stt_cmp = { { { '0%3' '3%6' '6%9'}                                                         { '0%3' '3%6' '6%9'}                                                         ; { '0%3'                           '3%6'                          '6%9'}                            {} } };

cfg.v_lne       = { [0 0.450 0.900] };
cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
                
% Run
cfg.dat_nme     = '_overall_data';
cfg.sbj_clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;

cfg.fle_out_pth = fcfg.dat_fld;

mmil_combo_effects(cfg);

% Save for overall analysis
for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

end