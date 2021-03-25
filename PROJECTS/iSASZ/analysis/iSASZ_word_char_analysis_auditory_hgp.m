function iSASZ_word_char_analysis_auditory_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial analysis work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

if any(eve_typ > 9)

%% Length
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'aud_lng_crr_600' };

cfg.sig_loc = { { 'aud_ovr_stt' 'aud_lng_stt' 'aud_lng_stt' 'aud_big_lng_stt' 'aud_big_lng_stt' 'aud_frq_lng_stt' 'aud_frq_lng_stt' 'aud_ngh_lng_stt' 'aud_ngh_lng_stt' } };
           
cfg.sig_chn = { [ 1 2 1 2 1 2 1 2 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Longer' 'Shorter' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Longer' 'Shorter' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Red' 'Dark Red' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
cfg.plt_prb_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };

cfg.plt_prb_cmb_max = { [ 1 1 3 4] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Longer' 'Shorter' } };
               
cfg.fst_cmb_col = { { 'Bright Red' 'Dark Red' } };
               
cfg.fst_cmb_rul = { { '2' '3' } };
               
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

cfg.alt_eve = { {           'aud_new_old'                      'aud_lng_med'                       'aud_big_med' ; ''                         'aud_frq_med'                             'aud_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [211 212]                          [231 232]                           [241 242] ; []                             [251 252]                                 [261 262] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'aud_ovr_stt'           'aud_new_old_stt_msk'}     {'aud_lng_stt'          'aud_lng_big_stt'       'aud_lng_frq_stt'         'aud_lng_ngh_stt'}         {'aud_big_stt'           'aud_big_lng_stt'      'aud_big_frq_stt'       'aud_big_ngh_stt'}         ; {} {'aud_frq_stt'             'aud_frq_lng_stt'      'aud_frq_big_stt'       'aud_frq_ngh_stt'}         {'aud_ngh_stt'             'aud_ngh_lng_stt'      'aud_ngh_big_stt'       'aud_ngh_frq_stt'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'aud_big_crr_600' };

cfg.sig_loc = { { 'aud_ovr_stt' 'aud_big_stt' 'aud_big_stt' 'aud_lng_big_stt' 'aud_lng_big_stt' 'aud_frq_big_stt' 'aud_frq_big_stt' 'aud_ngh_big_stt' 'aud_ngh_big_stt' } };
           
cfg.sig_chn = { [ 1 2 1 2 1 2 1 2 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Higher' 'Lower' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Higher' 'Lower' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Blue' 'Dark Blue' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
cfg.plt_prb_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };

cfg.plt_prb_cmb_max = { [ 1 1 3 4] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Higher' 'Lower' } };
               
cfg.fst_cmb_col = { { 'Bright Blue' 'Dark Blue' } };
               
cfg.fst_cmb_rul = { { '2' '3' } };
               
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

cfg.alt_eve = { {           'aud_new_old'                      'aud_lng_med'                       'aud_big_med' ; ''                         'aud_frq_med'                             'aud_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [211 212]                          [231 232]                           [241 242] ; []                             [251 252]                                 [261 262] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'aud_ovr_stt'           'aud_new_old_stt_msk'}     {'aud_lng_stt'          'aud_lng_big_stt'       'aud_lng_frq_stt'         'aud_lng_ngh_stt'}         {'aud_big_stt'           'aud_big_lng_stt'      'aud_big_frq_stt'       'aud_big_ngh_stt'}         ; {} {'aud_frq_stt'             'aud_frq_lng_stt'      'aud_frq_big_stt'       'aud_frq_ngh_stt'}         {'aud_ngh_stt'             'aud_ngh_lng_stt'      'aud_ngh_big_stt'       'aud_ngh_frq_stt'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'aud_frq_crr_600' };

cfg.sig_loc = { { 'aud_ovr_stt' 'aud_frq_stt' 'aud_frq_stt' 'aud_lng_frq_stt' 'aud_lng_frq_stt' 'aud_big_frq_stt' 'aud_big_frq_stt' 'aud_ngh_frq_stt' 'aud_ngh_frq_stt' } };
           
cfg.sig_chn = { [ 1 2 1 2 1 2 1 2 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Higher' 'Lower' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Higher' 'Lower' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Yellow' 'Dark Yellow' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
cfg.plt_prb_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };

cfg.plt_prb_cmb_max = { [ 1 1 3 4] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Higher' 'Lower' } };
               
cfg.fst_cmb_col = { { 'Bright Yellow' 'Dark Yellow' } };
               
cfg.fst_cmb_rul = { { '2' '3' } };
               
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

cfg.alt_eve = { {           'aud_new_old'                      'aud_lng_med'                       'aud_big_med' ; ''                         'aud_frq_med'                             'aud_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [211 212]                          [231 232]                           [241 242] ; []                             [251 252]                                 [261 262] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'aud_ovr_stt'           'aud_new_old_stt_msk'}     {'aud_lng_stt'          'aud_lng_big_stt'       'aud_lng_frq_stt'         'aud_lng_ngh_stt'}         {'aud_big_stt'           'aud_big_lng_stt'      'aud_big_frq_stt'       'aud_big_ngh_stt'}         ; {} {'aud_frq_stt'             'aud_frq_lng_stt'      'aud_frq_big_stt'       'aud_frq_ngh_stt'}         {'aud_ngh_stt'             'aud_ngh_lng_stt'      'aud_ngh_big_stt'       'aud_ngh_frq_stt'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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

%% Length
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'aud_ngh_crr_600' };

cfg.sig_loc = { { 'aud_ovr_stt' 'aud_ngh_stt' 'aud_ngh_stt' 'aud_lng_ngh_stt' 'aud_lng_ngh_stt' 'aud_big_ngh_stt' 'aud_big_ngh_stt' 'aud_frq_ngh_stt' 'aud_frq_ngh_stt'  } };
           
cfg.sig_chn = { [ 1 2 1 2 1 2 1 2 1] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'Higher' 'Lower' 'Overlap' 'Corrected' } };
               
cfg.tbl_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [1];

cfg.plt_cmb_nme = { { 'Higher' 'Lower' 'Overlap' 'Corrected' } };
               
cfg.plt_cmb_col = { { 'Bright Orange' 'Dark Orange' 'Bright Purple' 'Black' } };
               
cfg.plt_cmb_rul     = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };
cfg.plt_prb_cmb_rul = { { '1 & 2' '1 & 3' '1 & ( 2 | 3 ) & ( 4 | 5 | 6 | 7 | 8 | 9 )' '~ 1 & ( 2 | 3 )' } };

cfg.plt_prb_cmb_max = { [ 1 1 3 4] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'Higher' 'Lower' } };
               
cfg.fst_cmb_col = { { 'Bright Orange' 'Dark Orange' } };
               
cfg.fst_cmb_rul = { { '2' '3' } };
               
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

cfg.alt_eve = { {           'aud_new_old'                      'aud_lng_med'                       'aud_big_med' ; ''                         'aud_frq_med'                             'aud_ngh_med' ; '' '' '' } };
cfg.eve     = { {           [211 212]                          [231 232]                           [241 242] ; []                             [251 252]                                 [261 262] ; [] [] [] } };
cfg.lnstyle.col_ord = { { { rgb('bright red')   rgb('orange')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} ; {} {rgb('dark yellow') rgb('bright yellow')} {rgb('dark orange') rgb('bright orange')} ; {} {} {} } };
                     
cfg.stt_dat = { { { 'aud_ovr_stt'           'aud_new_old_stt_msk'}     {'aud_lng_stt'          'aud_lng_big_stt'       'aud_lng_frq_stt'         'aud_lng_ngh_stt'}         {'aud_big_stt'           'aud_big_lng_stt'      'aud_big_frq_stt'       'aud_big_ngh_stt'}         ; {} {'aud_frq_stt'             'aud_frq_lng_stt'      'aud_frq_big_stt'       'aud_frq_ngh_stt'}         {'aud_ngh_stt'             'aud_ngh_lng_stt'      'aud_ngh_big_stt'       'aud_ngh_frq_stt'} ; {} {} {} } };
cfg.stt_col = { { { ft_stt_col(rgb('gray')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} ; {} {ft_stt_col(rgb('yellow')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('orange'))} {ft_stt_col(rgb('orange')) ft_stt_col(rgb('red')) ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))} ; {} {} {} } };

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