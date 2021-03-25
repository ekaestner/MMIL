function FW_pap_ana_hgp_lex(fcfg)

fprintf([fcfg.sbj_nme ': Starting FW_pap_ana_hgp_lex on %s \n'],fcfg.sbj_nme)

%%
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_lex_600' };

cfg.sig_loc = { { 'vis_old' 'vis_wrd_big' 'vis_wrd_frq' 'vis_stm_01' 'fw_wrd_stt' 'fw_rep_stt' 'vis_old' 'vis_wrd_big' 'vis_wrd_frq' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 1 2 2 2 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'VisualNovelPriming' 'VisualBigram'      'VisualFrq'         'NonSpecific'} };
               
cfg.tbl_cmb_rul = { { '1 & 4 & ( 5 | 6 )'  '2 & 4 & ( 5 | 6 )' '3 & 4 & ( 5 | 6 )' '( ( 1 | 2 | 3 ) & ( ~ 4 | ~ ( 5 | 6 ) ) ) | ( ( 7 | 8 | 9 ) & 4 & ( 5 | 6 ) )'} };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'VisualNovelPriming' 'VisualBigram' 'VisualFrq'} };            
cfg.plt_cmb_col = { { 'orange'             'blue'         'yellow'    } };

cfg.plt_pct_cmb_nme = { { 'VisualNovelPriming' 'VisualBigram' 'VisualFrq' 'NonSpecific'} };
cfg.plt_pct_cmb_col = { { 'orange'             'blue'         'yellow'       'black'} };
               
cfg.plt_cmb_rul     = { { '1 & 4 & ( 5 | 6 )'  '2 & 4 & ( 5 | 6 )' '3 & 4 & ( 5 | 6 )' } };
cfg.plt_prb_cmb_rul = { { '1 & 4 & ( 5 | 6 )'  '2 & 4 & ( 5 | 6 )' '3 & 4 & ( 5 | 6 )' '( ( 1 | 2 | 3 ) & ( ~ 4 | ~ ( 5 | 6 ) ) ) | ( ( 7 | 8 | 9 ) & 4 & ( 5 | 6 ) )'} };

cfg.plt_prb_cmb_max = { [ 1 2 3 4] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'VisualNovelPriming' 'VisualBigram' 'VisualFrq' } };
               
cfg.fst_cmb_col = { { 'orange'             'blue'         'yellow'    } };
           
cfg.inc_cmb_rul = { { '1 & 4 & ( 5 | 6 )'  '2 & 4 & ( 5 | 6 )' '3 & 4 & ( 5 | 6 )' } };
cfg.fst_cmb_rul = { { '1'                  '2'                 '3'} };
               
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

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 0 ; 1 1 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'trialinfo' '' ; 'wrd_big_med' 'wrd_frq_med' } };
cfg.eve     = { { [3 4 6] [] ; [141 142] [121 122] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange')   rgb('reddish grey')} {} ; {rgb('dark blue')   rgb('bright blue')} {rgb('dark yellow')   rgb('bright yellow')} } };
                     
cfg.stt_dat = { { { 'vis_stm_01'             'vis_old' }                       {} ; { 'vis_stm_01' 'vis_wrd_big' }                            {'vis_stm_01' 'vis_wrd_frq'} } };
cfg.stt_col = { { { ft_stt_col(rgb('black')) ft_stt_col(rgb('orange')) } {} ; { ft_stt_col(rgb('black')) ft_stt_col(rgb('light blue'))} {ft_stt_col(rgb('black')) ft_stt_col(rgb('light yellow'))} } };
cfg.stt_cmp = { { { '0%5'                   '3vs4'                     } {} ; { '0%5' '141vs142' }                                      { '0%5' '121vs122' } } };

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