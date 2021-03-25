function FW_analysis_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial analysis work on %s \n'],fcfg.sbj_nme)

%% Responsive Effects
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'eff_600_lap_01' };

cfg.sig_loc = { { 'vis_stm' 'vis_stm_01' 'vis_ltr_msk' 'vis_ltr_msk' 'vis_wrd_msk' 'vis_wrd_msk' 'vis_old_msk' 'vis_old_msk' } };
           
cfg.sig_chn = { [ 1 1 1 2 1 2 1 2] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Non-Specific' 'Specific .05' 'Specific .01' 'Specific' } };
               
cfg.tbl_cmb_rul = { { '( 1 | 2 ) & ~ ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & ~ 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( ~ 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Non-Specific' 'Specific .05' 'Specific .01' 'Specific' } };
               
cfg.plt_cmb_col = { { 'black' 'light yellow' 'light orange' 'bright red'} };
               
cfg.plt_cmb_rul     = { { '( 1 | 2 ) & ~ ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & ~ 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( ~ 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' } };
cfg.plt_prb_cmb_rul = { { '( 1 | 2 ) & ~ ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & ~ 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( ~ 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' } };

cfg.plt_prb_cmb_max = { [ 1 1 1 1 ] };

% First Pass combinations
cfg.run_fst     = [ 0 ];

cfg.fst_cmb_nme = { { 'Visual Language' } };
               
cfg.fst_cmb_col = { { 'maroon' } };
               
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
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark red')   rgb('reddish grey')} } };
                     
cfg.stt_dat = { { {'vis_stm'                        'vis_ltr_msk'               'vis_wrd_msk'                 'vis_old_msk'} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('dark red')) ft_stt_col(rgb('bright red')) ft_stt_col(rgb('orange')) } } };

cfg.v_lne       = { [0 0.2 0.4] };
cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} };   
                
% Run
cfg.dat_nme     = '_overall_data';
cfg.sbj_clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;

cfg.fle_out_pth = fcfg.dat_fld;

mmil_combo_effects(cfg);

% Save for overall analysis
for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

end