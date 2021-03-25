function SL_pap_ana_lfp_act(fcfg)

fprintf([fcfg.sbj_nme ': Starting SL_pap_ana_hgp_act on %s \n'],fcfg.sbj_nme)

%%
cfg = [];

cfg.typ      = { 'lfp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_rsp_950' };

cfg.sig_loc = { { 'ovr_stt' 'ovr_stt' 'vis_anv_stt' 'vis_anv_stt' } };
           
cfg.sig_chn = { [ 1 2 1 2 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Specific'             'EventResponsive'   'SpecificNonSpecific'} };
               
cfg.tbl_cmb_rul = { { '( 1 & 3 ) | ( 2 & 4 )' '1 | 2' '~ ( 1 & 2 ) & ( 3 | 4 )'} };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Specific'             'EventResponsive' } };            
cfg.plt_cmb_col = { { 'red'                  'Black'           } };

cfg.plt_pct_cmb_nme = { { 'Specific'         'EventResponsive' 'SpecificNonSpecific'} };
cfg.plt_pct_cmb_col = { { 'red'              'Black'           'blue'} };
               
cfg.plt_cmb_rul     = { { '( 1 & 3 ) | ( 2 & 4 )' '1 | 2' } };
cfg.plt_prb_cmb_rul = { { '( 1 & 3 ) | ( 2 & 4 )' '1 | 2' '~ ( 1 & 2 ) & ( 3 | 4 )' } };

cfg.plt_prb_cmb_max = { [ 1 2 3] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'VisualLanguage' 'AuditoryLanguage'} };
               
cfg.fst_cmb_col = { { 'red' 'blue' } };
           
cfg.inc_cmb_rul = { { '1 & 3' '2 & 4'} };
cfg.fst_cmb_rul = { { '3'     '4' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Visual Language' } };
               
cfg.tme_cmb_col = { { 'maroon' } };
               
cfg.tme_cmb_rul = { { '1' } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 1 1 ];
cfg.dat_loc = [ 1 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'trialinfo' } };
cfg.eve     = { { [1 2 3 4] } };
cfg.lnstyle.col_ord = { { { rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey') } } };
                     
cfg.stt_dat = { { {'ovr_stt' 'vis_anv_stt' 'vis_anv_stt_01'} } };
cfg.stt_col = { { { ft_stt_col(rgb('grey')) ft_stt_col(rgb('green')) ft_stt_col(rgb('greenish grey')) } } };
cfg.stt_cmp = { { { '0%5' '5%10' '10%15' } } };

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