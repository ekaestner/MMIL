function FW_word_char_analysis_lfp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial analysis work on %s \n'],fcfg.sbj_nme)

%% Length
cfg = [];

cfg.typ      = { 'lfp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'lng_crr_600' };

cfg.sig_loc = { { 'fw_wrd_stt' 'vis_stm_01' 'vis_wrd_lng' 'vis_wrd_big_lng' 'vis_wrd_frq_lng' 'vis_wrd_ngh_lng' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'LengthSensitive' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3 ' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'LengthSensitive' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Red' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3 ' } };
cfg.plt_prb_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3 ' } };

cfg.plt_prb_cmb_max = { [ 1 2 3] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Longer' } };
               
cfg.fst_cmb_col = { { 'Bright Red' } };
               
cfg.fst_cmb_rul = { { '3' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { {  } };
               
cfg.tme_cmb_col = { { } };
               
cfg.tme_cmb_rul = { {  } };

cfg.tme_ord     = { [] };
cfg.tme_ord_nme = { [] };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 3 3 ];
cfg.dat_loc = [ 1 1 1 ; 0 1 1 ; 0 0 0 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'trialinfo' 'wrd_lng_med' 'wrd_big_med' ; '' 'wrd_frq_med' 'wrd_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [3 4 5 6]                                                                    [111 112]                           [141 142] ; []                             [121 122]                                 [131 132] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'vis_stm_01'            'vis_ltr_msk'                   'vis_wrd_msk'                  'vis_old_msk'                 } {'vis_wrd_lng'          'vis_wrd_lng_big'       'vis_wrd_lng_frq'         'vis_wrd_lng_ngh'}         {'vis_wrd_big'           'vis_wrd_big_lng'      'vis_wrd_big_frq'       'vis_wrd_big_ngh'}         ; {} {'vis_wrd_frq'             'vis_wrd_frq_lng'      'vis_wrd_frq_big'       'vis_wrd_frq_ngh'}         {'vis_wrd_ngh'             'vis_wrd_ngh_lng'      'vis_wrd_ngh_big'       'vis_wrd_ngh_frq'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('dark purple')) ft_stt_col(rgb('bright red')) } {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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

%% Bigram
cfg = [];

cfg.typ      = { 'lfp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'big_crr_600' };

cfg.sig_loc = { { 'fw_wrd_stt' 'vis_stm_01' 'vis_wrd_big' 'vis_wrd_lng_big' 'vis_wrd_frq_big' 'vis_wrd_ngh_big' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'BigramSensitive' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3 ' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'BigramSensitive' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Blue' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3 ' } };
cfg.plt_prb_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3 ' } };

cfg.plt_prb_cmb_max = { [ 1 2 3] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Higher' } };
               
cfg.fst_cmb_col = { { 'Bright Blue' } };
               
cfg.fst_cmb_rul = { { '3' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { {  } };
               
cfg.tme_cmb_col = { { } };
               
cfg.tme_cmb_rul = { {  } };

cfg.tme_ord     = { [] };
cfg.tme_ord_nme = { [] };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 3 3 ];
cfg.dat_loc = [ 1 1 1 ; 0 1 1 ; 0 0 0 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'trialinfo' 'wrd_lng_med' 'wrd_big_med' ; '' 'wrd_frq_med' 'wrd_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [3 4 5 6]                                                                    [111 112]                           [141 142] ; []                             [121 122]                                 [131 132] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'vis_stm_01'            'vis_ltr_msk'                   'vis_wrd_msk'                  'vis_old_msk'                 } {'vis_wrd_lng'          'vis_wrd_lng_big'       'vis_wrd_lng_frq'         'vis_wrd_lng_ngh'}         {'vis_wrd_big'           'vis_wrd_big_lng'      'vis_wrd_big_frq'       'vis_wrd_big_ngh'}         ; {} {'vis_wrd_frq'             'vis_wrd_frq_lng'      'vis_wrd_frq_big'       'vis_wrd_frq_ngh'}         {'vis_wrd_ngh'             'vis_wrd_ngh_lng'      'vis_wrd_ngh_big'       'vis_wrd_ngh_frq'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('dark purple')) ft_stt_col(rgb('bright red')) } {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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

%% Frequency
cfg = [];

cfg.typ      = { 'lfp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'frq_crr_600' };

cfg.sig_loc = { { 'fw_wrd_stt' 'vis_stm_01' 'vis_wrd_frq' 'vis_wrd_lng_frq' 'vis_wrd_big_frq' 'vis_wrd_ngh_frq' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'FrequencySensitive' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'FrequencySensitive' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Yellow' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3' } };
cfg.plt_prb_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3' } };

cfg.plt_prb_cmb_max = { [ 1 2 3] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'FrequencySensitive' } };
               
cfg.fst_cmb_col = { { 'Bright Yellow' } };
               
cfg.fst_cmb_rul = { { '3' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { {  } };
               
cfg.tme_cmb_col = { { } };
               
cfg.tme_cmb_rul = { {  } };

cfg.tme_ord     = { [] };
cfg.tme_ord_nme = { [] };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 3 3 ];
cfg.dat_loc = [ 1 1 1 ; 0 1 1 ; 0 0 0 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'trialinfo' 'wrd_lng_med' 'wrd_big_med' ; '' 'wrd_frq_med' 'wrd_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [3 4 5 6]                                                                    [111 112]                           [141 142] ; []                             [121 122]                                 [131 132] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'vis_stm_01'            'vis_ltr_msk'                   'vis_wrd_msk'                  'vis_old_msk'                 } {'vis_wrd_lng'          'vis_wrd_lng_big'       'vis_wrd_lng_frq'         'vis_wrd_lng_ngh'}         {'vis_wrd_big'           'vis_wrd_big_lng'      'vis_wrd_big_frq'       'vis_wrd_big_ngh'}         ; {} {'vis_wrd_frq'             'vis_wrd_frq_lng'      'vis_wrd_frq_big'       'vis_wrd_frq_ngh'}         {'vis_wrd_ngh'             'vis_wrd_ngh_lng'      'vis_wrd_ngh_big'       'vis_wrd_ngh_frq'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('dark purple')) ft_stt_col(rgb('bright red')) } {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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

%% Neighborhood
cfg = [];

cfg.typ      = { 'lfp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'ngh_crr_600' };

cfg.sig_loc = { { 'fw_wrd_stt' 'vis_stm_01' 'vis_wrd_ngh' 'vis_wrd_lng_ngh' 'vis_wrd_big_ngh' 'vis_wrd_frq_ngh' } };
           
cfg.sig_chn = { [ 1 1 1 1 1 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'NeighborhoodSensitive' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'NeighborhoodSensitive' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Orange' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3' } };
cfg.plt_prb_cmb_rul = { { '1 & 2 & 3' '1 & 2 & 3 & ( 4 | 5 | 6 )' '~ ( 1 & 2 ) & 3' } };

cfg.plt_prb_cmb_max = { [ 1 2 3] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'NeighborhoodSensitive' } };
               
cfg.fst_cmb_col = { { 'Bright Orange' } };
               
cfg.fst_cmb_rul = { { '3' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};   

% Time combinations
cfg.run_tme     = [0];

cfg.tme_cmb_nme = { {  } };
               
cfg.tme_cmb_col = { { } };
               
cfg.tme_cmb_rul = { {  } };

cfg.tme_ord     = { [] };
cfg.tme_ord_nme = { [] };   
    
% Line Plot combinations
cfg.run_lne     = [1];

cfg.plt_dim = [ 3 3 ];
cfg.dat_loc = [ 1 1 1 ; 0 1 1 ; 0 0 0 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'trialinfo' 'wrd_lng_med' 'wrd_big_med' ; '' 'wrd_frq_med' 'wrd_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [3 4 5 6]                                                                    [111 112]                           [141 142] ; []                             [121 122]                                 [131 132] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange') rgb('dark purple')   rgb('reddish grey')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'vis_stm_01'            'vis_ltr_msk'                   'vis_wrd_msk'                  'vis_old_msk'                 } {'vis_wrd_lng'          'vis_wrd_lng_big'       'vis_wrd_lng_frq'         'vis_wrd_lng_ngh'}         {'vis_wrd_big'           'vis_wrd_big_lng'      'vis_wrd_big_frq'       'vis_wrd_big_ngh'}         ; {} {'vis_wrd_frq'             'vis_wrd_frq_lng'      'vis_wrd_frq_big'       'vis_wrd_frq_ngh'}         {'vis_wrd_ngh'             'vis_wrd_ngh_lng'      'vis_wrd_ngh_big'       'vis_wrd_ngh_frq'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('dark purple')) ft_stt_col(rgb('bright red')) } {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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