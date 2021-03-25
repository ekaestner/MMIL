function SL_analysis_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial analysis work on %s \n'],fcfg.sbj_nme)

%% Responsive Effects
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'Responsitivity' };
% 
% cfg.sig_loc = { { 'ovr_stt' 'ovr_stt' } };
%            
% cfg.sig_chn = { [ 1 2 ] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'Early Visual' 'Early Auditory' 'Sustained Activity' } };
%                
% cfg.tbl_cmb_rul = { { '1 & ~ 2' '~ 1 & 2' '1 & 2' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [1];
% 
% cfg.plt_cmb_nme = { { 'Early Visual' 'Early Auditory' 'Sustained Activity' } };
%                
% cfg.plt_cmb_col = { { 'dark red' 'dark blue' 'purple' } };
%                
% cfg.plt_cmb_rul     = { { '1 & ~ 2' '~ 1 & 2' '1 & 2' } };
% cfg.plt_prb_cmb_rul = { { '1 & ~ 2' '~ 1 & 2' '1 & 2' } };
% 
% cfg.plt_prb_cmb_max = { [ 1 1 1 ] };
% 
% % First Pass combinations
% cfg.run_fst     = [ 0 ];
% 
% cfg.fst_cmb_nme = { { 'Early Visual' 'Early Auditory' } };
%                
% cfg.fst_cmb_col = { { 'dark red' 'dark blue' } };
%                
% cfg.fst_cmb_rul = { { '1' '2' } };
%                
% cfg.fst_ord     = {{}};
% cfg.fst_ord_dim = {{'avg'} };   
% cfg.fst_ord_nme = {{}};   
% 
% % Time combinations
% cfg.run_tme     = [0];
% 
% cfg.tme_cmb_nme = { { 'Early Visual' 'Early Auditory' } };
%                
% cfg.tme_cmb_col = { { 'dark red' 'dark blue' } };
%                
% cfg.tme_cmb_rul = { { '1' '2' } };
% 
% cfg.tme_ord     = { { 1 } };
% cfg.tme_ord_nme = { { 'vis' } };   
%     
% % Line Plot combinations
% cfg.run_lne     = [0];
% 
% cfg.plt_dim = [ 1 1 ];
% cfg.dat_loc = [ 1 ];
% 
% cfg.alt_lab = { 'stt_lab' };
% 
% cfg.alt_eve = { { 'ovr_all_eve' } };
% cfg.eve     = { { [101] } };
% cfg.lnstyle.col_ord = { { { rgb('black') } } };
%                      
% cfg.stt_dat = { { {'ovr_stt'} } };
% cfg.stt_col = { { { ft_stt_col(rgb('grey')) } } };
% 
% cfg.v_lne       = { [0 0.450 0.900] };
% cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects(cfg);
% 
% % Save for overall analysis
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% Sensitive Electrodes
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'Selectivity' };
% 
% cfg.sig_loc = { { 'vis_anv_stt' 'vis_anv_stt' 'ovr_stt' 'ovr_stt' } };
%            
% cfg.sig_chn = { [ 1 2 1 2] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'Corrected'  'Early Visual Sensitivity' 'Early Auditory Sensitivity' 'Sustained Sensitivity' 'Early Visual Activity' 'Early Auditory Activity' 'Sustained Activity' } };
%                
% cfg.tbl_cmb_rul = { { '( 1 | 2 ) & ( ~ 3 & ~ 4 )' '1 & ~ 2 & ( 3 | 4 )' '~ 1 & 2 & ( 3 | 4 )' '1 & 2  & ( 3 | 4 )' '3 & ~ 4' '~ 3 & 4' '3 & 4' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [1];
% 
% cfg.plt_cmb_nme = { { 'Corrected'  'Early Visual Sensitivity' 'Early Auditory Sensitivity' 'Sustained Sensitivity' 'Early Visual Activity' 'Early Auditory Activity' 'Sustained Activity' } };
%                
% cfg.plt_cmb_col = { { 'yellow'     'dark green'                  'bright green'                  'cyan'                     'bright red'            'bright blue'             'bright purple' } };
%                
% cfg.plt_cmb_rul     = { { '( 1 | 2 ) & ( ~ 3 & ~ 4 )' '1 & ~ 2 & ( 3 | 4 )' '~ 1 & 2 & ( 3 | 4 )' '1 & 2 & ( 3 | 4 )' '3 & ~ 4' '~ 3 & 4' '3 & 4' } };
% cfg.plt_prb_cmb_rul = { { '( 1 | 2 ) & ( ~ 3 & ~ 4 )' '1 & ~ 2 & ( 3 | 4 )' '~ 1 & 2 & ( 3 | 4 )' '1 & 2 & ( 3 | 4 )' '3 & ~ 4' '~ 3 & 4' '3 & 4' } };
% 
% cfg.plt_prb_cmb_max = { [ 1 2 2 2 3 3 3] };
% 
% % First Pass combinations
% cfg.run_fst     = [ 0 ];
% 
% cfg.fst_cmb_nme = { { 'Early Visual' 'Early Auditory' } };
%                
% cfg.fst_cmb_col = { { 'dark red' 'dark blue' } };
%                
% cfg.fst_cmb_rul = { { '1' '2' } };
%                
% cfg.fst_ord     = {{}};
% cfg.fst_ord_dim = {{'avg'} };   
% cfg.fst_ord_nme = {{}};   
% 
% % Time combinations
% cfg.run_tme     = [0];
% 
% cfg.tme_cmb_nme = { { 'Early Visual' 'Early Auditory' } };
%                
% cfg.tme_cmb_col = { { 'dark red' 'dark blue' } };
%                
% cfg.tme_cmb_rul = { { '1' '2' } };
% 
% cfg.tme_ord     = { { 1 } };
% cfg.tme_ord_nme = { { 'vis' } };   
%     
% % Line Plot combinations
% cfg.run_lne     = [1];
% 
% cfg.plt_dim = [ 1 1 ];
% cfg.dat_loc = [ 1 ];
% 
% cfg.alt_lab = { 'stt_lab' };
% 
% cfg.alt_eve = { { 'trialinfo' } };
% cfg.eve     = { { [1 2 3 4] } };
% cfg.lnstyle.col_ord = { { { rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey') } } };
%                      
% cfg.stt_dat = { { {'ovr_stt' 'vis_anv_stt' 'vis_anv_stt_01'} } };
% cfg.stt_col = { { { ft_stt_col(rgb('grey')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('green')) } } };
% 
% cfg.v_lne       = { [0 0.450 0.900] };
% cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects(cfg);
% 
% % Save for overall analysis 
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% Linguistic Effects
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'Language_Selectivity' };
% 
% cfg.sig_loc = { { 'vis_nse_stt_msk' 'aud_nse_stt_msk' 'ovr_stt' 'ovr_stt' 'vis_anv_stt' 'vis_anv_stt' } };
%            
% cfg.sig_chn = { [ 1 1 1 2 1 2] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'Early Visual Language' 'Early Auditory Language' 'Bi-modal Activity Language' 'Corrected'} };
%                
% cfg.tbl_cmb_rul = { { '1 & 3 & 5' '2 & 4 & 6' '1 & 2 & 3 & 4 & 5 & 6' '( 1 & ~ ( 1 & ~ 2 & 3 & 5 ) ) | ( 2 & ~ ( ~ 1 & 2 & 4 & 6 ) )' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [1];
% 
% cfg.plt_cmb_nme = { { 'Early Visual Language' 'Early Auditory Language' 'Bi-modal Activity Language' 'Corrected'} };
%                
% cfg.plt_cmb_col = { { 'bright red' 'bright blue' 'bright purple' 'black'} };
%                
% cfg.plt_cmb_rul     = { { '1 & ~ 2 & 3 & 5' '~ 1 & 2 & 4 & 6' '1 & 2 & 3 & 4 & 5 & 6' '( 1 & ~ ( 1 & ~ 2 & 3 & 5 ) ) | ( 2 & ~ ( ~ 1 & 2 & 4 & 6 ) )' } };
% cfg.plt_prb_cmb_rul = { { '1 & 3 & 5'       '2 & 4 & 6'       '1 & 2 & 3 & 4 & 5 & 6' '( 1 & ~ ( 1 & ~ 2 & 3 & 5 ) ) | ( 2 & ~ ( ~ 1 & 2 & 4 & 6 ) )' } };
% 
% cfg.plt_prb_cmb_max = { [ 1 1 2 2] };
% 
% % First Pass combinations
% cfg.run_fst     = [ 1 ];
% 
% cfg.fst_cmb_nme = { { 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.fst_cmb_col = { { 'bright red' 'bright blue' } };
%                
% cfg.fst_cmb_rul = { { '1' '2' } };
%                
% cfg.fst_ord     = {{}};
% cfg.fst_ord_dim = {{'avg'} };   
% cfg.fst_ord_nme = {{}};   
% 
% % Time combinations
% cfg.run_tme     = [0];
% 
% cfg.tme_cmb_nme = { { 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.tme_cmb_col = { { 'bright red' 'bright blue' } };
%                
% cfg.tme_cmb_rul = { { '1' '2' } };
% 
% cfg.tme_ord     = { { 1 } };
% cfg.tme_ord_nme = { { 'vis' } };   
%     
% % Line Plot combinations
% cfg.run_lne     = [1];
% 
% cfg.plt_dim = [ 2 2 ];
% cfg.dat_loc = [ 1 1 ; 0 0];
% 
% cfg.alt_lab = { 'stt_lab' };
% 
% cfg.alt_eve = { { 'vis_tot_nse' 'aud_tot_nse' ; '' ''} };
% cfg.eve     = { { [101 102] [111 112] ; [] []} };
% cfg.lnstyle.col_ord = { { { rgb('red') rgb('reddish grey') } {rgb('blue') rgb('bluish grey')} ; {} {} } };
%                      
% cfg.stt_dat = { { {'vis_nse_stt'} {'aud_nse_stt'} ; {''} {''} } };
% cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) } { ft_stt_col(rgb('bluish grey')) } ; {} {} } };
% 
% cfg.v_lne       = { [0 0.450 0.900] };
% cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects(cfg);
% 
% % Save for overall analysis
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% Unfamiliar Effects
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'Control_Selectivity' };
% 
% cfg.sig_loc = { { 'vis_nse_stt_msk' 'aud_nse_stt_msk' 'vis_nse_stt_msk' 'aud_nse_stt_msk' 'ovr_stt' 'ovr_stt' 'vis_anv_stt' 'vis_anv_stt' } };
%            
% cfg.sig_chn = { [ 2 2 1 1 1 2 1 2] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'Early Visual Control' 'Early Auditory Control' 'Bi-modal Activity Control' 'Early Visual Language' 'Early Auditory Language'   'Bi-modal Activity Language'     'Corrected'} };
%                
% cfg.tbl_cmb_rul = { { '1 & 5 & 7'            '2 & 6 & 8'              '1 & 2 & 5 & 6 & 7 & 8'     '( 1 | 2 ) & 3 & 5 & 7' '( 1 | 2 ) & 4 & 6 & 8' '( 1 | 2 ) & 3 & 4 & 5 & 6 & 7 & 8'  '( 1 & ~ 5 & ~ 7 ) | ( 2 & ~ 6 & ~ 8 )' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [1];
% 
% cfg.plt_cmb_nme = { { 'Early Visual Control' 'Early Auditory Control' 'Bi-modal Activity Control' 'Early Visual Language' 'Early Auditory Language' 'Bi-modal Activity Language' 'Corrected' } };
%                
% cfg.plt_cmb_col = { { 'yellow' 'cyan' 'orange' 'bright red' 'bright blue' 'bright purple' 'black' } };
%                
% cfg.plt_cmb_rul     = { { '1 & ~ 2 & 5 & 7'      '~ 1 & 2 & 6 & 8' '1 & 2 & 5 & 6 & 7 & 8'     '( 1 | 2 ) & 3 & ~ 4 & 5 & 7' '( 1 | 2 ) & ~ 3 & 4 & 6 & 8' '( 1 | 2 ) & 3 & 4 & 5 & 6 & 7 & 8' '( 1 & ~ 5 & ~ 7 ) | ( 2 & ~ 6 & ~ 8 )' } };
% cfg.plt_prb_cmb_rul = { { '1 & 5 & 7'            '2 & 6 & 8'       '1 & 2 & 5 & 6 & 7 & 8'     '( 1 | 2 ) & 3 & 5 & 7'       '( 1 | 2 ) & 4 & 6 & 8'       '( 1 | 2 ) & 3 & 4 & 5 & 6 & 7 & 8' '( 1 & ~ 5 & ~ 7 ) | ( 2 & ~ 6 & ~ 8 )' } };
% 
% cfg.plt_prb_cmb_max = { [ 1 1 1 2 2 2 3] };
% 
% % First Pass combinations
% cfg.run_fst     = [ 1 ];
% 
% cfg.fst_cmb_nme = { { 'Early Visual Control' 'Early Auditory Control' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.fst_cmb_col = { { 'yellow' 'cyan' 'bright red' 'bright blue' } };
%                
% cfg.fst_cmb_rul = { { '1' '2' '3' '4' } };
%                
% cfg.fst_ord     = {{}};
% cfg.fst_ord_dim = {{'avg'} };   
% cfg.fst_ord_nme = {{}};   
% 
% % Time combinations
% cfg.run_tme     = [0];
% 
% cfg.tme_cmb_nme = { { 'Early Visual Control' 'Early Auditory Control' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.tme_cmb_col = { { 'yellow' 'cyan' 'bright red' 'bright blue' } };
%                
% cfg.tme_cmb_rul = { { '1' '2' '3' '4' } };
% 
% cfg.tme_ord     = { { 1 } };
% cfg.tme_ord_nme = { { 'vis' } };   
%     
% % Line Plot combinations
% cfg.run_lne     = [1];
% 
% cfg.plt_dim = [ 2 2 ];
% cfg.dat_loc = [ 1 1 ; 0 0];
% 
% cfg.alt_lab = { 'stt_lab' };
% 
% cfg.alt_eve = { { 'vis_tot_nse' 'aud_tot_nse' ; '' ''} };
% cfg.eve     = { { [101 102] [111 112] ; [] []} };
% cfg.lnstyle.col_ord = { { { rgb('red') rgb('reddish grey') } {rgb('blue') rgb('bluish grey')} ; {} {} } };
%                      
% cfg.stt_dat = { { {'vis_nse_stt'} {'aud_nse_stt'} ; {''} {''} } };
% cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) } { ft_stt_col(rgb('bluish grey')) } ; {} {} } };
% 
% cfg.v_lne       = { [0 0.450 0.900] };
% cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects(cfg);
% 
% % Save for overall analysis
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% Stimulus Driven Mismatch/Match Effects  - Corrected
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'Stimulus_Match_Mismatch' };
% 
% cfg.sig_loc = { { 'vis_mtc_stt_msk' 'vis_mtc_stt_msk' 'vis_nse_stt_msk' 'aud_nse_stt_msk' 'ovr_stt' } };
%            
% cfg.sig_chn = { [ 1 3 1 1 2] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'Match' 'Mismatch'  'Early Visual Language' 'Early Auditory Language' 'Bi-modal Activity Language' 'Corrected' } };
%                
% cfg.tbl_cmb_rul = { { '1 & 5' '2 & 5' '( 1 | 2 ) & 3 & 5' '( 1 | 2 ) & 4 & 5' '( 1 | 2 ) & 3 & 4 & 5' '( 1 | 2 ) & ~ 5' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [1];
% 
% cfg.plt_cmb_nme = { { 'Match' 'Mismatch'  'Early Visual Language' 'Early Auditory Language' 'Bi-modal Activity Language' 'Corrected'} };
%                
% cfg.plt_cmb_col = { { 'green' 'yellow' 'bright red' 'bright blue' 'bright purple' 'black' } };
%                
% cfg.plt_cmb_rul     = { { '1 & 5' '2 & 5' '( 1 | 2 ) & 3 & ~ 4 & 5' '( 1 | 2 ) & ~ 3 & 4 & 5' '( 1 | 2 ) & 3 & 4 & 5' '( 1 | 2 ) & ~ 5' } };
% cfg.plt_prb_cmb_rul = { { '1 & 5' '2 & 5' '( 1 | 2 ) & 3 & 5'       '( 1 | 2 ) & 4 & 5'       '( 1 | 2 ) & 3 & 4 & 5' '( 1 | 2 ) & ~ 5'  } };
% 
% cfg.plt_prb_cmb_max = { [ 1 1 2 2 2 3 ] };
% 
% % First Pass combinations
% cfg.run_fst     = [ 1 ];
% 
% cfg.fst_cmb_nme = { { 'Match' 'Mismatch' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.fst_cmb_col = { { 'green' 'yellow' 'bright red' 'bright blue' } };
%                
% cfg.fst_cmb_rul = { { '1' '2' '3' '4' } };
%                
% cfg.fst_ord     = {{}};
% cfg.fst_ord_dim = {{'avg'} };   
% cfg.fst_ord_nme = {{}};   
% 
% % Time combinations
% cfg.run_tme     = [0];
% 
% cfg.tme_cmb_nme = { { 'Match' 'Mismatch' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.tme_cmb_col = { { 'yellow' 'cyan' 'bright red' 'bright blue' } };
%                
% cfg.tme_cmb_rul = { { '1' '2' '3' '4' } };
% 
% cfg.tme_ord     = { { 1 } };
% cfg.tme_ord_nme = { { 'vis' } };   
%     
% % Line Plot combinations
% cfg.run_lne     = [1];
% 
% cfg.plt_dim = [ 2 2 ];
% cfg.dat_loc = [ 1 1 ; 1 0];
% 
% cfg.alt_lab = { 'stt_lab' };
% 
% cfg.alt_eve = { { 'vis_tot_nse' 'aud_tot_nse' ; 'trialinfo' ''} };
% cfg.eve     = { { [101 102] [111 112] ; [1 2] []} };
% cfg.lnstyle.col_ord = { { { rgb('red') rgb('reddish grey') } {rgb('blue') rgb('bluish grey')} ; {rgb('green') rgb('yellow')} {} } };
%                      
% cfg.stt_dat = { { {'vis_nse_stt'} {'aud_nse_stt'} ; {'vis_mtc_stt_msk'} {''} } };
% cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) } { ft_stt_col(rgb('bluish grey')) } ; { ft_stt_col(rgb('yellow')) } {} } };
% 
% cfg.v_lne       = { [0 0.450 0.900] };
% cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects(cfg);
% 
% % Save for overall analysis
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% Post-Stimulus Match/Mismatch  - Corrected
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'PostStimulus_Match_Mismatch' };
% 
% cfg.sig_loc = { { 'vis_mtc_stt_msk' 'vis_mtc_stt_msk' 'vis_nse_stt_msk' 'aud_nse_stt_msk' 'ovr_stt' } };
%            
% cfg.sig_chn = { [ 2 4 1 1 2] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'Match' 'Mismatch'  'Early Visual Language' 'Early Auditory Language' 'Bi-modal Activity Language' 'Corrected' } };
%                
% cfg.tbl_cmb_rul = { { '1 & 5' '2 & 5' '( 1 | 2 ) & 3 & 5' '( 1 | 2 ) & 4 & 5' '( 1 | 2 ) & 3 & 4 & 5' '( 1 | 2 ) & ~ 5' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [1];
% 
% cfg.plt_cmb_nme = { { 'Match' 'Mismatch'  'Early Visual Language' 'Early Auditory Language' 'Bi-modal Activity Language' 'Corrected'} };
%                
% cfg.plt_cmb_col = { { 'green' 'yellow' 'bright red' 'bright blue' 'bright purple' 'black' } };
%                
% cfg.plt_cmb_rul     = { { '1 & 5' '2 & 5' '( 1 | 2 ) & 3 & ~ 4 & 5' '( 1 | 2 ) & ~ 3 & 4 & 5' '( 1 | 2 ) & 3 & 4 & 5' '( 1 | 2 ) & ~ 5' } };
% cfg.plt_prb_cmb_rul = { { '1 & 5' '2 & 5' '( 1 | 2 ) & 3 & 5'       '( 1 | 2 ) & 4 & 5'       '( 1 | 2 ) & 3 & 4 & 5' '( 1 | 2 ) & ~ 5'  } };
% 
% cfg.plt_prb_cmb_max = { [ 1 1 2 2 2 3 ] };
% 
% % First Pass combinations
% cfg.run_fst     = [ 1 ];
% 
% cfg.fst_cmb_nme = { { 'Match' 'Mismatch' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.fst_cmb_col = { { 'green' 'yellow' 'bright red' 'bright blue' } };
%                
% cfg.fst_cmb_rul = { { '1' '2' '3' '4' } };
%                
% cfg.fst_ord     = {{}};
% cfg.fst_ord_dim = {{'avg'} };   
% cfg.fst_ord_nme = {{}};   
% 
% % Time combinations
% cfg.run_tme     = [0];
% 
% cfg.tme_cmb_nme = { { 'Match' 'Mismatch' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.tme_cmb_col = { { 'yellow' 'cyan' 'bright red' 'bright blue' } };
%                
% cfg.tme_cmb_rul = { { '1' '2' '3' '4' } };
% 
% cfg.tme_ord     = { { 1 } };
% cfg.tme_ord_nme = { { 'vis' } };   
%     
% % Line Plot combinations
% cfg.run_lne     = [1];
% 
% cfg.plt_dim = [ 2 2 ];
% cfg.dat_loc = [ 1 1 ; 1 0];
% 
% cfg.alt_lab = { 'stt_lab' };
% 
% cfg.alt_eve = { { 'vis_tot_nse' 'aud_tot_nse' ; 'trialinfo' ''} };
% cfg.eve     = { { [101 102] [111 112] ; [1 2] []} };
% cfg.lnstyle.col_ord = { { { rgb('red') rgb('reddish grey') } {rgb('blue') rgb('bluish grey')} ; {rgb('green') rgb('yellow')} {} } };
%                      
% cfg.stt_dat = { { {'vis_nse_stt'} {'aud_nse_stt'} ; {'vis_mtc_stt_msk'} {''} } };
% cfg.stt_col = { { { ft_stt_col(rgb('reddish grey')) } { ft_stt_col(rgb('bluish grey')) } ; { ft_stt_col(rgb('yellow')) } {} } };
% 
% cfg.v_lne       = { [0 0.450 0.900] };
% cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects(cfg);
% 
% % Save for overall analysis
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% Word/Non-Word - corrected
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'Word_NonWord' };
% 
% cfg.sig_loc = { { 'vis_wrd_stt_msk' 'vis_wrd_stt_msk' 'aud_wrd_stt_msk' 'aud_wrd_stt_msk' 'ovr_stt' 'ovr_stt' } };
%            
% cfg.sig_chn = { [ 1 2 1 2 1 2] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'Visual Word' 'Visual NonWord'  'Auditory Word' 'Auditory NonWord' 'Corrected' } };
%                
% cfg.tbl_cmb_rul = { { '1 & 5' '2 & 5' '3 & 6' '4 & 6' '( ( 1 | 2 ) & ~ 5 ) | ( ( 3 | 4 ) & ~ 6 )' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [1];
% 
% cfg.plt_cmb_nme = { { 'Visual Word' 'Visual NonWord'  'Auditory Word' 'Auditory NonWord' 'Corrected'} };
%                
% cfg.plt_cmb_col = { { 'bright red' 'dark red' 'bright blue' 'dark blue' 'black' } };
%                
% cfg.plt_cmb_rul     = { { '1 & 5' '2 & 5' '3 & 6' '4 & 6' '( ( 1 | 2 ) & ~ 5 ) | ( ( 3 | 4 ) & ~ 6 )' } };
% cfg.plt_prb_cmb_rul = { { '1 & 5' '2 & 5' '3 & 6' '4 & 6' '( ( 1 | 2 ) & ~ 5 ) | ( ( 3 | 4 ) & ~ 6 )'  } };
% 
% cfg.plt_prb_cmb_max = { [ 1 1 1 1 2 ] };
% 
% % First Pass combinations
% cfg.run_fst     = [ 1 ];
% 
% cfg.fst_cmb_nme = { { 'Match' 'Mismatch' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.fst_cmb_col = { { 'green' 'yellow' 'bright red' 'bright blue' } };
%                
% cfg.fst_cmb_rul = { { '1' '2' '3' '4' } };
%                
% cfg.fst_ord     = {{}};
% cfg.fst_ord_dim = {{'avg'} };   
% cfg.fst_ord_nme = {{}};   
% 
% % Time combinations
% cfg.run_tme     = [0];
% 
% cfg.tme_cmb_nme = { { 'Match' 'Mismatch' 'Early Visual Language' 'Early Auditory Language' } };
%                
% cfg.tme_cmb_col = { { 'yellow' 'cyan' 'bright red' 'bright blue' } };
%                
% cfg.tme_cmb_rul = { { '1' '2' '3' '4' } };
% 
% cfg.tme_ord     = { { 1 } };
% cfg.tme_ord_nme = { { 'vis' } };   
%     
% % Line Plot combinations
% cfg.run_lne     = [1];
% 
% cfg.plt_dim = [ 2 2 ];
% cfg.dat_loc = [ 1 1 ; 1 0];
% 
% cfg.alt_lab = { 'stt_lab' };
% 
% cfg.alt_eve = { { 'vis_wrd' 'aud_wrd' ; 'trialinfo' ''} };
% cfg.eve     = { { [101 102] [201 202] ; [1 2 3 4] []} };
% cfg.lnstyle.col_ord = { { { rgb('bright red') rgb('dark red') } {rgb('bright blue') rgb('dark blue')} ; {rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey')} {} } };
%                      
% cfg.stt_dat = { { {'vis_wrd_stt'} {'aud_wrd_stt'} ; {'vis_nse_stt' 'aud_nse_stt' 'vis_mtc_stt_msk'} { ''} } };
% cfg.stt_col = { { { ft_stt_col(rgb('dark red')) } { ft_stt_col(rgb('dark blue')) } ; { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('bluish grey')) ft_stt_col(rgb('yellow')) } {} } };
% 
% cfg.v_lne       = { [0 0.450 0.900] };
% cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects(cfg);
% 
% % Save for overall analysis
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% Letter/Phoneme
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'Grapheme_Phoneme' };

cfg.sig_loc = { { 'vis_con_stt_msk' 'vis_con_stt_msk' 'aud_con_stt_msk' 'aud_con_stt_msk' 'ovr_stt' 'ovr_stt' } };
           
cfg.sig_chn = { [ 1 2 1 2 1 2 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'EarlyVisualLetter' 'LateVisualLetter' 'EarlyAuditoryPhoneme' 'Transition' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 5' '~ 1 & 2 & 5' '3 & 6' '2 & 3 & 5 & 6' '( ( 1 | 2 ) & ~ 5 ) | ( ( 3 | 4 ) & ~ 6 )' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'EarlyVisualLetter' 'LateVisualLetter' 'EarlyAuditoryPhoneme' 'Transition' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'red' 'magenta' 'blue' 'purple' 'black' } };
               
cfg.plt_cmb_rul     = { { '1 & 5' '~ 1 & 2 & 5' '3 & 6' '2 & 3 & 5 & 6' '( ( 1 | 2 ) & ~ 5 ) | ( ( 3 | 4 ) & ~ 6 )' } };
cfg.plt_prb_cmb_rul = { { '1 & 5' '~ 1 & 2 & 5' '3 & 6' '2 & 3 & 5 & 6' '( ( 1 | 2 ) & ~ 5 ) | ( ( 3 | 4 ) & ~ 6 )' } };

cfg.plt_prb_cmb_max = { [ 1 1 1 2 3 ] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'ElyVisualLetter' 'LteVisualLetter' 'ElyAuditoryPhoneme' 'LteAuditoryPhoneme' } };
               
cfg.fst_cmb_col = { { 'red' 'dark red' 'blue' 'dark blue' } };
               
cfg.fst_cmb_rul = { { '1' '2' '3' '4' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { { 'Match' 'Mismatch' 'Early Visual Language' 'Early Auditory Language' } };
               
cfg.tme_cmb_col = { { 'yellow' 'cyan' 'bright red' 'bright blue' } };
               
cfg.tme_cmb_rul = { { '1' '2' '3' '4' } };

cfg.tme_ord     = { { 1 } };
cfg.tme_ord_nme = { { 'vis' } };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 2 2 ];
cfg.dat_loc = [ 1 1 ; 1 0];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'vis_con' 'aud_con' ; 'trialinfo' ''} };
cfg.eve     = { { [1:12] [1:12] ; [1 2 3 4] []} };
eme_col = distinguishable_colors(12);
cfg.lnstyle.col_ord = { { mat2cell(eme_col,ones(1,size(eme_col,1)),3)' mat2cell(eme_col,ones(1,size(eme_col,1)),3)' ; {rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey')} {} } };
                     
cfg.stt_dat = { { {'vis_con_stt_msk'} {'aud_con_stt_msk'} ; {'vis_nse_stt' 'aud_nse_stt' 'vis_mtc_stt_msk'} { ''} } };
cfg.stt_col = { { { ft_stt_col(rgb('magenta')) } { ft_stt_col(rgb('purple')) } ; { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('bluish grey')) ft_stt_col(rgb('yellow')) } {} } };

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