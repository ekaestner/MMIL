function FW_pap_ana_hgp_act(fcfg)

fprintf([fcfg.sbj_nme ': Starting FW_pap_ana_hgp_act on %s \n'],fcfg.sbj_nme)

%%
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_rsp_600' };

cfg.sig_loc = { { 'vis_stm_01' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Specific'             'EventResponsive'   'SpecificNonSpecific'} };
               
cfg.tbl_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 )' '( 2 | 3 | 4 | 5 )' '1 & ~ ( 2 | 3 | 4 | 5 )'} };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Specific'             'EventResponsive' } };            
cfg.plt_cmb_col = { { 'red'                  'Black'           } };

cfg.plt_pct_cmb_nme = { { 'Specific'         'EventResponsive' 'SpecificNonSpecific'} };
cfg.plt_pct_cmb_col = { { 'red'              'Black'           'blue'} };
               
cfg.plt_cmb_rul     = { { '1 & ( 2 | 3 | 4 | 5 )' '( 2 | 3 | 4 | 5 )' } };
cfg.plt_prb_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 )' '( 2 | 3 | 4 | 5 )' '1 & ~ ( 2 | 3 | 4 | 5 )'} };

cfg.plt_prb_cmb_max = { [ 1 2 3] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'VisualLanguage' } };
               
cfg.fst_cmb_col = { { 'maroon' } };
           
cfg.inc_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 )' } };
cfg.fst_cmb_rul = { { '1' } };
               
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
cfg.eve     = { { [3 4 5 6] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} } };
                     
cfg.stt_dat = { { { 'vis_stm_01'             'fw_ffn_stt'                    'fw_ltr_stt'                   'fw_wrd_stt'                  'fw_rep_stt'       } } };
cfg.stt_col = { { { ft_stt_col(rgb('black')) ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('dark purple')) ft_stt_col(rgb('bright red')) ft_stt_col(rgb('orange')) } } };

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