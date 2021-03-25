function SL_pap_ana_hgp_lng(fcfg)

fprintf([fcfg.sbj_nme ': Starting SL_pap_ana_hgp_lng on %s \n'],fcfg.sbj_nme)

%%
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_lng_950' };

cfg.sig_loc = { { 'ovr_stt' 'ovr_stt' 'pap_anv_stt' 'pap_anv_stt' 'vis_nse_stt_msk_anv' 'aud_nse_stt_msk_anv' } };
           
cfg.sig_chn = { [ 1 2 1 2 1 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'VisualLanguage' 'AuditoryLanguage' 'LanguageNonSpecific'} };
                
cfg.tbl_cmb_rul = { { '1 & 3 & 5'          '2 & 4 & 6'            '( 5 | 6 ) & ( ~ ( 1 & 3 & 5 ) & ~ ( 2 & 4 & 6 ) )'} };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'VisualLanguage' 'AuditoryLanguage' } };            
cfg.plt_cmb_col = { { 'red'            'blue'             } };

cfg.plt_pct_cmb_nme = { { 'VisualLanguage'   'AuditoryLanguage' 'LanguageNonSpecific'} };
cfg.plt_pct_cmb_col = { { 'red'              'blue'             'black'} };
               
cfg.plt_cmb_rul     = { { '1 & 3 & 5'          '2 & 4 & 6'            } };
cfg.plt_prb_cmb_rul = { { '1 & 3 & 5'          '2 & 4 & 6'            '( 5 | 6 ) & ( ~ ( 1 & 3 & 5 ) & ~ ( 2 & 4 & 6 ) )' } };

cfg.plt_prb_cmb_max = { [ 1 2 3] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'VisualLanguage' 'AuditoryLanguage'} };
               
cfg.fst_cmb_col = { { 'red' 'blue' } };
           
cfg.inc_cmb_rul = { { '1 & 3 & 5' '2 & 4 & 6' } };
cfg.fst_cmb_rul = { { '5'         '6' } };
               
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

cfg.alt_eve = { { 'lng_tot_nse' } };
cfg.eve     = { { [601 603 604] } };
cfg.lnstyle.col_ord = { { { rgb('greenish yellow') rgb('reddish grey') rgb('bluish grey') } } };
                     
cfg.stt_dat = { { { 'ovr_stt'                'vis_anv_stt'            'vis_nse_stt_msk_anv'  'aud_nse_stt_msk_anv'} } };
cfg.stt_col = { { { ft_stt_col(rgb('black')) ft_stt_col(rgb('green')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue'))  } } };
cfg.stt_cmp = { { { '0%3'                    '3%6'                    '6%9'                 '9%12' } } };

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