function iSASZ_analysis_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial analysis work on %s \n'],fcfg.sbj_nme)

%% Responsive Effect
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'eff_800' };

cfg.sig_loc = { { 'vis_ovr_stt' 'aud_ovr_stt' } };
           
cfg.sig_chn = { [ 1 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Visual Responsive' 'Auditory Responsive' 'Bi-Modal Responsive'} };
               
cfg.tbl_cmb_rul = { { '1' '2' '1 & 2' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Visual Responsive' 'Auditory Responsive' 'Bi-Modal Responsive' } };
               
cfg.plt_cmb_col = { { 'red' 'blue' 'purple' } };
               
cfg.plt_cmb_rul     = { { '1 & ~ 2' '~ 1 & 2' '1 & 2' } };
cfg.plt_prb_cmb_rul = { { '1'       '2'       '1 & 2'} };

cfg.plt_prb_cmb_max = { [ 1 1 2] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Visual Responsive' 'Auditory Responsive' } };
               
cfg.fst_cmb_col = { { 'red' 'blue' } };
               
cfg.fst_cmb_rul = { { '1' '2' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Visual Responsive' 'Auditory Responsive' } };
               
cfg.tme_cmb_col = { { 'red' 'blue' } };
               
cfg.tme_cmb_rul = { { '1' '2' } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 1 ; 0 0];

cfg.alt_lab = { 'label' };

cfg.alt_eve = { { 'vis_ovr' 'aud_ovr' ; '' '' } };
cfg.eve     = { { [101]  [201]; [] []} };
cfg.lnstyle.col_ord = { { { rgb('reddish grey') } { rgb('bluish grey') } ; {} {} } };
                     
cfg.stt_dat = { { {'vis_ovr_stt'} {'aud_ovr_stt'} ; {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) } { ft_stt_col(rgb('bluish gray')) } ; {} {} } };

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

cfg.cmb_nme = { 'ani_obj_800' };

cfg.sig_loc = { { 'vis_ani_obj_stt_msk' 'vis_ani_obj_stt_msk' 'aud_ani_obj_stt_msk' 'aud_ani_obj_stt_msk' 'vis_ovr_stt' 'aud_ovr_stt'} };
           
cfg.sig_chn = { [ 1 2 1 2 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Visual Animal' 'Visual Object' 'Auditory Animal' 'Auditory Object' 'Bi-Modal Animal' 'Bi-Modal Object'} };
               
cfg.tbl_cmb_rul = { { '1 & 5' '2 & 5' '3 & 6' '4 & 6' '1 & 3 & 5 & 6' '2 & 4 & 5 & 6' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Visual Animal' 'Visual Object' 'Auditory Animal' 'Auditory Object' 'Bi-Modal Animal' 'Bi-Modal Object' } };
               
cfg.plt_cmb_col = { { 'bright red' 'rust' 'bright blue' 'yellow' 'indigo' 'black' } };
               
cfg.plt_cmb_rul     = { { '1 & ~ 3 & 5' '2 & ~ 4 & 5' '~ 1 & 3 & 6' '~ 2 & 4 & 6' '1 & 3 & 5 & 6' '2 & 4 & 5 & 6' } };
cfg.plt_prb_cmb_rul = { { '1 & 5'       '2 & 5'       '3 & 6'       '4 & 6'       '1 & 3 & 5 & 6' '2 & 4 & 5 & 6'} };

cfg.plt_prb_cmb_max = { [ 1 1 1 1 2 2] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Visual Novel' 'Visual Repetition' 'Auditory Novel' 'Auditory Repetition' } };
               
cfg.fst_cmb_col = { { 'lime' 'reddish grey' 'aqua' 'bluish grey' } };
               
cfg.fst_cmb_rul = { { '1' '2' '3' '4'} };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Visual Responsive' 'Auditory Responsive' } };
               
cfg.tme_cmb_col = { { 'red' 'blue' } };
               
cfg.tme_cmb_rul = { { '1' '2' } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 1 ; 0 0];

cfg.alt_lab = { 'label' };

cfg.alt_eve = { { 'vis_ani_obj' 'aud_ani_obj' ; '' '' } };
cfg.eve     = { { [121 122]  [221 222]; [] []} };
cfg.lnstyle.col_ord = { { { rgb('lime') rgb('reddish grey') } { rgb('aqua') rgb('bluish grey') } ; {} {} } };
                     
cfg.stt_dat = { { { 'vis_ovr_stt' 'vis_new_old_stt_msk'}                                    {'aud_ovr_stt'                       'aud_new_old_stt_msk'} ; {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('reddish orange')) } { ft_stt_col(rgb('bluish grey')) ft_stt_col(rgb('yellowish orange')) } ; {} {} } };

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

cfg.cmb_nme = { 'nov_rep_800' };

cfg.sig_loc = { { 'vis_new_old_stt_msk' 'vis_new_old_stt_msk' 'aud_new_old_stt_msk' 'aud_new_old_stt_msk' 'vis_ovr_stt' 'aud_ovr_stt' } };
           
cfg.sig_chn = { [ 1 2 1 2 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Visual Novel' 'Visual Repetition' 'Auditory Novel' 'Auditory Repetition' 'Bi-Modal Novel' 'Bi-Modal Repetition'} };
               
cfg.tbl_cmb_rul = { { '1 & 5' '2 & 5' '3 & 6' '4 & 6' '1 & 3 & 5 & 6' '2 & 4 & 5 & 6' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Visual Novel' 'Visual Repetition' 'Auditory Novel' 'Auditory Repetition' 'Bi-Modal Novel' 'Bi-Modal Repetition' } };
               
cfg.plt_cmb_col = { { 'bright red' 'rust' 'bright blue' 'yellow' 'purple' 'bright orange' } };
               
cfg.plt_cmb_rul     = { { '1 & ~ 3' '2 & ~ 4' '~ 1 & 3' '~ 2 & 4' '1 & 3' '2 & 4' } };
cfg.plt_prb_cmb_rul = { { '1'       '2'       '3'       '4'       '1 & 3' '2 & 4'} };

cfg.plt_prb_cmb_max = { [ 1 1 1 1 2 2] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Visual Novel' 'Visual Repetition' 'Auditory Novel' 'Auditory Repetition' } };
               
cfg.fst_cmb_col = { { 'bright red' 'rust' 'bright blue' 'yellow' } };
               
cfg.fst_cmb_rul = { { '1' '2' '3' '4'} };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Visual Responsive' 'Auditory Responsive' } };
               
cfg.tme_cmb_col = { { 'red' 'blue' } };
               
cfg.tme_cmb_rul = { { '1' '2' } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 1 ; 0 0];

cfg.alt_lab = { 'label' };

cfg.alt_eve = { { 'vis_new_old' 'aud_new_old' ; '' '' } };
cfg.eve     = { { [111 112]  [211 212]; [] []} };
cfg.lnstyle.col_ord = { { { rgb('bright red') rgb('rust') } { rgb('bright blue') rgb('yellow') } ; {} {} } };
                     
cfg.stt_dat = { { { 'vis_ovr_stt' 'vis_new_old_stt_msk'}                                    {'aud_ovr_stt'                       'aud_new_old_stt_msk'} ; {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('reddish orange')) } { ft_stt_col(rgb('bluish grey')) ft_stt_col(rgb('yellowish orange')) } ; {} {} } };

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