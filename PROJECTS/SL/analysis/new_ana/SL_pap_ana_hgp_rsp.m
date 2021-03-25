function SL_pap_ana_hgp_rsp(fcfg)

fprintf([fcfg.sbj_nme ': Starting SL_pap_ana_hgp_rsp on %s \n'],fcfg.sbj_nme)

%% OVR_STT
% cfg = [];
% 
% cfg.typ      = { 'hgp' };
% cfg.ele_typ  = { 'ecog' };
% cfg.loc_typ  = { 'split' };
% 
% cfg.cmb_nme = { 'pap_rsp_1500' };
% 
% cfg.sig_loc = { { 'pap_ovr_stt' 'pap_ovr_stt' 'pap_ovr_stt' 'pap_ovr_stt' } };
%            
% cfg.sig_chn = { [ 1 2 3 4] };
% 
% % Table Combination
% cfg.run_tbl     = [1];
% 
% cfg.tbl_cmb_nme = { { 'VisualActive' 'AuditoryActive' 'ResponseActive' 'OverallActive' } };
%                
% cfg.tbl_cmb_rul = { { '1' '2' '3' '4' } };
%                
% cfg.tbl_sum = { [] };
% 
% % Plot combination
% cfg.run_plt     = [ 0 ];
% 
% % First Pass combinations
% cfg.run_fst     = [ 0 ];
% 
% % Time combinations
% cfg.run_tme     = [ 0 ];
%     
% % Line Plot combinations
% cfg.run_lne     = [ 0 ];
%                 
% % Run
% cfg.dat_nme     = '_overall_data';
% cfg.sbj_clr_fld = fcfg.clr_fld;
% cfg.sbj_nme = fcfg.sbj_nme;
% 
% cfg.fle_out_pth = fcfg.dat_fld;
% 
% mmil_combo_effects2(cfg);
% 
% % Save for overall analysis
% for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% OVR_ANV
cfg = [];

cfg.typ      = { 'hgp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_anv_1500' };

cfg.sig_loc = { { 'pap_ovr_stt' 'pap_ovr_stt' 'pap_ovr_stt' 'pap_ovr_stt' 'pap_anv_stt' 'pap_anv_stt' 'pap_anv_stt' 'pap_anv_stt' } };
           
cfg.sig_chn = { [ 1 2 3 4 1 2 3 4 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'VisualSelective' 'AuditorySelective'   'ResponseSelective' 'OverallSelective' } };
               
cfg.tbl_cmb_rul = { { '1 & 5'           '2 & 6'               '3 & 7'             '4 & 8' } };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [ 0 ];

cfg.plt_cmb_nme = { { 'VisualSelective' 'AuditorySelective'   'ResponseSelective' 'OverallSelective' } };            
cfg.plt_cmb_col = { { 'red'                  'blue'           'green' 'black' } };

cfg.plt_pct_cmb_nme = { { 'VisualSelective' 'AuditorySelective'   'ResponseSelective' 'OverallSelective' } };
cfg.plt_pct_cmb_col = { { 'red'              'blue'           'green' 'black'} };
               
cfg.plt_cmb_rul     = { { '1 & 5' '2 & 6' '3 & 7' '4 & 8' } };
cfg.plt_prb_cmb_rul = { { '1 & 5' '2 & 6' '3 & 7' '4 & 8' } };

cfg.plt_prb_cmb_max = { [ 1 1 1 1 ] };

% First Pass combinations
cfg.run_fst     = [ 1 ];

cfg.fst_cmb_nme = { { 'VisualSelective' 'AuditorySelective'   'ResponseSelective' 'OverallSelective' } };
               
cfg.fst_cmb_col = { { 'red'             'blue'                'green'             'black' } };
           
cfg.inc_cmb_rul = { { '1 & 5'           '2 & 6'               '3 & 7'             '4 & 8' } };
cfg.fst_cmb_rul = { { '5'               '6'                   '7'                 '8' } };
               
cfg.fst_ord     = {{}};
cfg.fst_ord_dim = {{'avg'} };   
cfg.fst_ord_nme = {{}};  

% Time combinations
cfg.run_tme     = [ 0 ];
    
% Line Plot combinations
cfg.run_lne     = [0];

cfg.plt_dim = [ 1 1 ];
cfg.dat_loc = [ 1 ];

cfg.alt_lab = { 'stt_lab' };

cfg.alt_eve = { { 'trialinfo' } };
cfg.eve     = { { [1 2 3 4] } };
cfg.lnstyle.col_ord = { { { rgb('green') rgb('yellow') rgb('reddish grey') rgb('bluish grey') } } };
                     
cfg.stt_dat = { { {'ovr_stt'                'vis_anv_stt' } } };
cfg.stt_col = { { { ft_stt_col(rgb('grey')) ft_stt_col(rgb('green')) } } };
cfg.stt_cmp = { { { '0%5'                   '5%10' } } };

cfg.v_lne       = { [0 0.450 0.900] };
cfg.v_lne_col   = { {rgb('red') rgb('blue') rgb('black')} };   

cfg.xlm         = [-0.2 1.5];

% Run
cfg.dat_nme     = '_overall_data';
cfg.sbj_clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;

cfg.fle_out_pth = fcfg.dat_fld;

mmil_combo_effects2(cfg);

% Save for overall analysis
for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

end