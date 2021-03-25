function FW_analysis_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial analysis work on %s \n'],fcfg.sbj_nme)

%% Responsive Effects
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'rsp_600' };

cfg.sig_loc = { { 'vis_ovr_stt' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'OverallResponsive' 'EventResponsive' 'Both' } };
               
cfg.tbl_cmb_rul = { { '1' '2 | 3 | 4 | 5' '1 & ( 2 | 3 | 4 | 5 )' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'OverallResponsive' 'EventResponsive' 'Both' } };
               
cfg.plt_cmb_col = { { 'black' 'Blue' 'red' } };
               
cfg.plt_cmb_rul     = { { '1 & ~ ( 2 | 3 | 4 | 5 )' '~ 1 & ( 2 | 3 | 4 | 5 )' '1 & ( 2 | 3 | 4 | 5 )' } };
cfg.plt_prb_cmb_rul = { { '1'                       '( 2 | 3 | 4 | 5 )'       '1 & ( 2 | 3 | 4 | 5 )' } };

cfg.plt_prb_cmb_max = { [ 1 1 2] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'OverallResponsive' 'FalseFont' 'Letters' 'Words' 'Repetition' } };
               
cfg.fst_cmb_col = { { 'black' 'reddish gray' 'dark purple' 'bright red' 'orange' } };
               
cfg.fst_cmb_rul = { { '1' '2' '3' '4' '5'} };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { '' } };
               
cfg.tme_cmb_col = { { 'black' 'reddish gray' 'dark purple' 'bright red' 'orange' } };
               
cfg.tme_cmb_rul = { { '1' '2' '3' '4' '5' } };

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
                     
cfg.stt_dat = { { { 'vis_ovr_stt'           'fw_ffn_stt'                    'fw_ltr_stt'                'fw_wrd_stt'                  'fw_rep_stt' } } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('dark purple')) ft_stt_col(rgb('bright red')) ft_stt_col(rgb('orange')) } } };

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

%% Selective Effects
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'eff_600_lap_01' };

cfg.sig_loc = { { 'vis_stm' 'vis_stm_01' 'vis_ltr_msk' 'vis_ltr_msk' 'vis_wrd_msk' 'vis_wrd_msk' 'vis_old_msk' 'vis_old_msk' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' } };
           
cfg.sig_chn = { [ 1 1 1 2 1 2 1 2 1 1 1 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Non-Specific' 'Specific .05' 'Specific .01' 'Specific' 'EventResponsive'} };
               
cfg.tbl_cmb_rul = { { '( 1 | 2 ) & ~ ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & ~ 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( ~ 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '9 | 10 | 11 | 12' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Non-Specific' 'Specific .05' 'Specific .01' 'Specific' 'EventResponsive'} };
               
cfg.plt_cmb_col = { { 'black' 'light yellow' 'light orange' 'bright red' 'Blue' } };
               
cfg.plt_cmb_rul     = { { '( 1 | 2 ) & ~ ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & ~ 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( ~ 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '9 | 10 | 11 | 12' } };
cfg.plt_prb_cmb_rul = { { '( 1 | 2 ) & ~ ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & ~ 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( ~ 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '( 1 & 2 ) & ( 3 | 4 | 5 | 6 | 7 | 8 )' '9 | 10 | 11 | 12' } };

cfg.plt_prb_cmb_max = { [ 1 1 1 1 2] };

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
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} } };
                     
cfg.stt_dat = { { { 'vis_stm'                       'vis_stm_01'                    'vis_ltr_msk'                  'vis_wrd_msk'                'vis_old_msk'} } };
cfg.stt_col = { { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('bright yellow')) ft_stt_col(rgb('dark purple')) ft_stt_col(rgb('bright red')) ft_stt_col(rgb('orange')) } } };

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

%% Letter/False Font Effects
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'ltr_ffn_600' };

cfg.sig_loc = { { 'vis_stm_01' 'vis_ltr_msk' 'vis_ltr_msk' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' } };
           
cfg.sig_chn = {[1 1 2 1 1 1 1]};

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Letter Selective' 'False Font Selective' 'Corrected'} };
               
cfg.tbl_cmb_rul = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };
               
cfg.tbl_sum = {[]};

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Letter Selective' 'False Font Selective' 'Corrected'} };
               
cfg.plt_cmb_col = { { 'dark purple'   'reddish grey' 'Black'} };
               
cfg.plt_cmb_rul     = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )' } };
cfg.plt_prb_cmb_rul = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )' } };

cfg.plt_prb_cmb_max = {[1 1 2]};

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Letter Selective' 'Voice Selective' } };
               
cfg.fst_cmb_col = { { 'dark purple' 'reddish grey' } };
               
cfg.fst_cmb_rul = { { '2' '3' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Letter Selective' 'Voice Selective' } };
               
cfg.tme_cmb_col = { { 'dark purple' 'reddish grey' } };
               
cfg.tme_cmb_rul = {{'1' '2'} };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [2 2];
cfg.dat_loc = [1 1 ; 1 1];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { {'trialinfo' 'trialinfo' ; 'trialinfo' 'trialinfo' } };
cfg.eve     = { {[3 4 5 6] [5 6]; [3 5] [3 4]} };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark purple')   rgb('reddish grey')} ; {rgb('bright red')   rgb('dark purple')} {rgb('bright red')   rgb('orange')} } };
                     
cfg.stt_dat = { {{'vis_stm'} {'vis_ltr_msk'} ; {'vis_wrd_msk'} {'vis_old_msk'} } };
cfg.stt_col = { {{ft_stt_col(rgb('reddish gray'))} {ft_stt_col(rgb('dark purple'))} ; {ft_stt_col(rgb('bright red'))} {ft_stt_col(rgb('orange'))} } };

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

%% Word/Non-Word Effects
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'wrd_nwd_600' };

cfg.sig_loc = { { 'vis_stm_01' 'vis_wrd_msk' 'vis_wrd_msk' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' } };
           
cfg.sig_chn = {[1 1 2 1 1 1 1 ]};

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Word Selective' 'Non-Word Selective' 'Corrected'} };
               
cfg.tbl_cmb_rul = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };
               
cfg.tbl_sum = {[]};

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Word Selective' 'Non-Word Selective' 'Corrected'} };
               
cfg.plt_cmb_col = { { 'bright red'   'dark purple' 'black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };
cfg.plt_prb_cmb_rul = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };

cfg.plt_prb_cmb_max = {[1 1 2]};

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Word Selective' 'Non-Word Selective' } };
               
cfg.fst_cmb_col = { { 'bright red'   'dark purple' } };
               
cfg.fst_cmb_rul = { { '2' '3' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Word Selective' 'Non-Word Selective' } };
               
cfg.tme_cmb_col = { { 'bright red'   'dark purple' } };
               
cfg.tme_cmb_rul = {{'1' '2'} };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [2 2];
cfg.dat_loc = [1 1 ; 1 1];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { {'trialinfo' 'trialinfo' ; 'trialinfo' 'trialinfo' } };
cfg.eve     = { {[3 4 5 6] [5 6]; [3 5] [3 4]} };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark purple')   rgb('reddish grey')} ; {rgb('bright red')   rgb('dark purple')} {rgb('bright red')   rgb('orange')} } };
                     
cfg.stt_dat = { {{'vis_stm'} {'vis_ltr_msk'} ; {'vis_wrd_msk'} {'vis_old_msk'} } };
cfg.stt_col = { {{ft_stt_col(rgb('reddish gray'))} {ft_stt_col(rgb('dark purple'))} ; {ft_stt_col(rgb('bright red'))} {ft_stt_col(rgb('orange'))} } };

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

%% Novel/Repetition Effects
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'nov_rep_600' };

cfg.sig_loc = { { 'vis_stm_01' 'vis_old_msk' 'vis_old_msk' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' } };
           
cfg.sig_chn = { [ 1 1 2 1 1 1 1 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Novel Selective' 'Repetition Selective' 'Corrected'} };
               
cfg.tbl_cmb_rul = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };
               
cfg.tbl_sum = {[]};

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Novel Selective' 'Repetition Selective' 'Corrected'} };
               
cfg.plt_cmb_col = { { 'bright red'   'orange' 'black'} };
               
cfg.plt_cmb_rul     = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };
cfg.plt_prb_cmb_rul = { { '1 & 2 & ( 4 | 5 | 6 | 7 )' '1 & 3 & ( 4 | 5 | 6 | 7 )' '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };

cfg.plt_prb_cmb_max = {[1 1 2]};

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Novel Selective' 'Repetition Selective' } };
               
cfg.fst_cmb_col = { { 'bright red'   'orange' } };
               
cfg.fst_cmb_rul = { { '2' '3' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Novel Selective' 'Repetition Selective' } };
               
cfg.tme_cmb_col = { { 'bright red'   'orange' } };
               
cfg.tme_cmb_rul = {{'1' '2'} };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [2 2];
cfg.dat_loc = [1 1 ; 1 1];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { {'trialinfo' 'trialinfo' ; 'trialinfo' 'trialinfo' } };
cfg.eve     = { {[3 4 5 6] [5 6]; [3 5] [3 4]} };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark purple')   rgb('reddish grey')} ; {rgb('bright red')   rgb('dark purple')} {rgb('bright red')   rgb('orange')} } };
                     
cfg.stt_dat = { {{'vis_stm'} {'vis_ltr_msk'} ; {'vis_wrd_msk'} {'vis_old_msk'} } };
cfg.stt_col = { {{ft_stt_col(rgb('reddish gray'))} {ft_stt_col(rgb('dark purple'))} ; {ft_stt_col(rgb('bright red'))} {ft_stt_col(rgb('orange'))} } };

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

%% GrainSize
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'GrainSize' };

cfg.sig_loc = { { 'vis_stm_01' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' 'vis_ltr_msk' 'vis_wrd_msk' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'LetterOnly' 'WordOnly' 'Both' } };
               
cfg.tbl_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 ) & 6' '1 & ( 2 | 3 | 4 | 5 ) & 7' '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7'} };
               
cfg.tbl_sum = {[]};

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'LetterOnly' 'WordOnly' 'Both' } };
               
cfg.plt_cmb_col = { { 'dark purple'   'bright red' 'bright purple'} };
               
cfg.plt_cmb_rul     = { { '1 & ( 2 | 3 | 4 | 5 ) & 6 & ~ 7' '1 & ( 2 | 3 | 4 | 5 ) & ~ 6 & 7' '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7'} };
cfg.plt_prb_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 ) & 6'       '1 & ( 2 | 3 | 4 | 5 ) & 7'       '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7'} };

cfg.plt_prb_cmb_max = {[1 1 2]};

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'LetterOnly' 'WordOnly' } };
               
cfg.fst_cmb_col = { { 'dark purple' 'bright red' } };
               
cfg.fst_cmb_rul = { { '6' '7' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Novel Selective' 'Repetition Selective' } };
               
cfg.tme_cmb_col = { { 'bright red'   'orange' } };
               
cfg.tme_cmb_rul = {{'1' '2'} };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [2 2];
cfg.dat_loc = [1 1 ; 1 1];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { {'trialinfo' 'trialinfo' ; 'trialinfo' 'trialinfo' } };
cfg.eve     = { {[3 4 5 6] [5 6]; [3 5] [3 4]} };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark purple')   rgb('reddish grey')} ; {rgb('bright red')   rgb('dark purple')} {rgb('bright red')   rgb('orange')} } };
                     
cfg.stt_dat = { {{'vis_stm'} {'vis_ltr_msk'} ; {'vis_wrd_msk'} {'vis_old_msk'} } };
cfg.stt_col = { {{ft_stt_col(rgb('reddish gray'))} {ft_stt_col(rgb('dark purple'))} ; {ft_stt_col(rgb('bright red'))} {ft_stt_col(rgb('orange'))} } };

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

%% NovelN4
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'NovelN4' };

cfg.sig_loc = { { 'vis_stm_01' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' 'vis_ltr_msk' 'vis_wrd_msk' 'vis_old_msk' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 1 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'LetterNovel' 'WordNovel' 'BothNovel' 'NonSpecificNovel' } };
               
cfg.tbl_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 ) & 6 & 8' '1 & ( 2 | 3 | 4 | 5 ) & 7 & 8' '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7 & 8' '1 & ( 2 | 3 | 4 | 5 ) & 8' } };
               
cfg.tbl_sum = {[]};

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'LetterNovel' 'WordNovel' 'BothNovel' 'NonSpecificNovel' } };
               
cfg.plt_cmb_col = { { 'dark purple' 'bright red' 'bright purple' 'orange'} };
               
cfg.plt_cmb_rul     = { { '1 & ( 2 | 3 | 4 | 5 ) & 6 & ~ 7 & 8' '1 & ( 2 | 3 | 4 | 5 ) & ~ 6 & 7 & 8' '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7 & 8' '1 & ( 2 | 3 | 4 | 5 ) & ~ 6 & ~ 7 & 8'} };
cfg.plt_prb_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 ) & 6 & 8'       '1 & ( 2 | 3 | 4 | 5 ) & 7 & 8'       '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7 & 8' '1 & ( 2 | 3 | 4 | 5 ) & 8'} };

cfg.plt_prb_cmb_max = {[1 1 2 3]};

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'LetterNovel' 'WordNovel'  'NovelN4' } };
               
cfg.fst_cmb_col = { { 'dark purple' 'bright red' 'orange' } };
               
cfg.fst_cmb_rul = { { '6' '7' '8' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Novel Selective' 'Repetition Selective' } };
               
cfg.tme_cmb_col = { { 'bright red'   'orange' } };
               
cfg.tme_cmb_rul = {{'1' '2'} };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [2 2];
cfg.dat_loc = [1 1 ; 1 1];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { {'trialinfo' 'trialinfo' ; 'trialinfo' 'trialinfo' } };
cfg.eve     = { {[3 4 5 6] [5 6]; [3 5] [3 4]} };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark purple')   rgb('reddish grey')} ; {rgb('bright red')   rgb('dark purple')} {rgb('bright red')   rgb('orange')} } };
                     
cfg.stt_dat = { {{'vis_stm'} {'vis_ltr_msk'} ; {'vis_wrd_msk'} {'vis_old_msk'} } };
cfg.stt_col = { {{ft_stt_col(rgb('reddish gray'))} {ft_stt_col(rgb('dark purple'))} ; {ft_stt_col(rgb('bright red'))} {ft_stt_col(rgb('orange'))} } };

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

%% UnfamiliarStimuli
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'UnfamiliarStimuli' };

cfg.sig_loc = { { 'vis_stm_01' 'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt' 'vis_ltr_msk' 'vis_wrd_msk' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 2 2 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'FalseFont' 'NonWord' 'Both' } };
               
cfg.tbl_cmb_rul = { { '1 & ( 2 | 3 | 4 | 5 ) & 6' '1 & ( 2 | 3 | 4 | 5 ) & 7' '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7'} };
               
cfg.tbl_sum = {[]};

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'FalseFont' 'NonWord' 'Both' } };
               
cfg.plt_cmb_col = { { 'reddish grey'   'dark purple' 'bright magenta'} };
               
cfg.plt_cmb_rul     = { { '1 & ( 2 | 3 | 4 | 5 ) & 6 & ~ 7' '1 & ( 2 | 3 | 4 | 5 ) & ~ 6 & 7' '1 & ( 2 | 3 | 4 | 5 ) & 6 & 7'} };
cfg.plt_prb_cmb_rul = { { '1 & 2 & ( 4 | 5 | 6 | 7 )'       '1 & 3 & ( 4 | 5 | 6 | 7 )'       '( 2 | 3 ) & ( ~ 4 & ~ 5 & ~ 6 & ~ 7 )'} };

cfg.plt_prb_cmb_max = {[1 1 2]};

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'FalseFont' 'NonWord' } };
               
cfg.fst_cmb_col = { { 'reddish grey' 'dark purple' } };
               
cfg.fst_cmb_rul = { { '6' '7' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Novel Selective' 'Repetition Selective' } };
               
cfg.tme_cmb_col = { { 'bright red'   'orange' } };
               
cfg.tme_cmb_rul = {{'1' '2'} };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [2 2];
cfg.dat_loc = [1 1 ; 1 1];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { {'trialinfo' 'trialinfo' ; 'trialinfo' 'trialinfo' } };
cfg.eve     = { {[3 4 5 6] [5 6]; [3 5] [3 4]} };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark purple')   rgb('reddish grey')} ; {rgb('bright red')   rgb('dark purple')} {rgb('bright red')   rgb('orange')} } };
                     
cfg.stt_dat = { {{'vis_stm'} {'vis_ltr_msk'} ; {'vis_wrd_msk'} {'vis_old_msk'} } };
cfg.stt_col = { {{ft_stt_col(rgb('reddish gray'))} {ft_stt_col(rgb('dark purple'))} ; {ft_stt_col(rgb('bright red'))} {ft_stt_col(rgb('orange'))} } };

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