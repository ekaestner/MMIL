function iSASZ_analysis_visual_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial visual analysis work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

if any(eve_typ < 9)

%% Responsive Effect
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'eff_800_vis' };

cfg.sig_loc = { { 'vis_ovr_stt' } };
           
cfg.sig_chn = { [ 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Visual Responsive' } };
               
cfg.tbl_cmb_rul = { { '1' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Visual Responsive' } };
               
cfg.plt_cmb_col = { { 'red' } };
               
cfg.plt_cmb_rul     = { { '1' } };
cfg.plt_prb_cmb_rul = { { '1' } };

cfg.plt_prb_cmb_max = { [ 1 ] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Visual Responsive' } };
               
cfg.fst_cmb_col = { { 'red' } };
               
cfg.fst_cmb_rul = { { '1' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Visual Responsive' } };
               
cfg.tme_cmb_col = { { 'red' } };
               
cfg.tme_cmb_rul = { { '1' } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 1 1 ];
cfg.dat_loc = [ 1 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'vis_ovr' } };
cfg.eve     = { { [101] } };
cfg.lnstyle.col_ord = { { { rgb('reddish grey') } } };
                     
cfg.stt_dat = { { {'vis_ovr_stt'} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) } } };

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


%% Semantic Effect
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'ani_obj_800_vis' };

cfg.sig_loc = { { 'vis_ani_obj_stt_msk' 'vis_ani_obj_stt_msk' 'vis_ovr_stt' } };
           
cfg.sig_chn = { [ 1 2 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Visual Animal' 'Visual Object' } };
               
cfg.tbl_cmb_rul = { { '1 & 3' '2 & 3' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Visual Animal' 'Visual Object' } };
               
cfg.plt_cmb_col = { { 'bright red' 'rust' } };
               
cfg.plt_cmb_rul     = { { '1 & 3' '2 & 3' } };
cfg.plt_prb_cmb_rul = { { '1 & 3' '2 & 3' } };

cfg.plt_prb_cmb_max = { [ 1 1] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Visual Novel' 'Visual Repetition' } };
               
cfg.fst_cmb_col = { { 'lime' 'reddish grey' } };
               
cfg.fst_cmb_rul = { { '1' '2' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Visual Responsive' } };
               
cfg.tme_cmb_col = { { 'red' } };
               
cfg.tme_cmb_rul = { { '1'  } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 1 1 ];
cfg.dat_loc = [ 1 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'vis_ani_obj' } };
cfg.eve     = { { [121 122]  } };
cfg.lnstyle.col_ord = { { { rgb('lime') rgb('reddish grey') } } };
                     
cfg.stt_dat = { { {'vis_ovr_stt'                   'vis_new_old_stt_msk'} } };
cfg.stt_col = { { { ft_stt_col(rgb('bluish grey')) ft_stt_col(rgb('yellowish orange')) } } };

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

%% Repetition Effect
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'nov_rep_800_vis' };

cfg.sig_loc = { { 'vis_new_old_stt_msk' 'vis_new_old_stt_msk' 'vis_ovr_stt'} };
           
cfg.sig_chn = { [ 1 2 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Visual Novel' 'Visual Repetition'} };
               
cfg.tbl_cmb_rul = { { '1 & 3' '2 & 3' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Visual Novel' 'Visual Repetition' } };
               
cfg.plt_cmb_col = { { 'bright red' 'rust'  } };
               
cfg.plt_cmb_rul     = { { '1 & 3' '2 & 3' } };
cfg.plt_prb_cmb_rul = { { '1 & 3' '2 & 3' } };

cfg.plt_prb_cmb_max = { [ 1 1] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Visual Novel' 'Visual Repetition' } };
               
cfg.fst_cmb_col = { { 'bright red' 'rust' } };
               
cfg.fst_cmb_rul = { { '1' '2' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Visual Responsive' } };
               
cfg.tme_cmb_col = { { 'red' } };
               
cfg.tme_cmb_rul = { { '1' } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 1 1 ];
cfg.dat_loc = [ 1 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'vis_new_old' } };
cfg.eve     = { { [111 112] } };
cfg.lnstyle.col_ord = { { { rgb('bright red') rgb('rust') } } };
                     
cfg.stt_dat = { { { 'vis_ovr_stt' 'vis_new_old_stt_msk'} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('reddish orange')) } } };

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

end