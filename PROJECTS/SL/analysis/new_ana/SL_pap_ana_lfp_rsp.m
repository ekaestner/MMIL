function SL_pap_ana_lfp_rsp(fcfg)

fprintf([fcfg.sbj_nme ': Starting SL_pap_ana_hgp_rsp on %s \n'],fcfg.sbj_nme)

%% OVR_STT
cfg = [];

cfg.typ      = { 'lfp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_rsp_1500' };

cfg.sig_loc = { { 'pap_ovr_stt' 'pap_ovr_stt' 'pap_ovr_stt' } };
           
cfg.sig_chn = { [ 1 2 3 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'VisualActive' 'AuditoryActive'   'ResponseActive'} };
               
cfg.tbl_cmb_rul = { { '1' '2' '3'} };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [ 0 ];

% First Pass combinations
cfg.run_fst     = [ 0 ];

% Time combinations
cfg.run_tme     = [ 0 ];
    
% Line Plot combinations
cfg.run_lne     = [ 0 ];
                
% Run
cfg.dat_nme     = '_overall_data';
cfg.sbj_clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;

cfg.fle_out_pth = fcfg.dat_fld;

mmil_combo_effects2(cfg);

% Save for overall analysis
for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

%% OVR_ANV
cfg = [];

cfg.typ      = { 'lfp' };
cfg.ele_typ  = { 'ecog' };
cfg.loc_typ  = { 'split' };

cfg.cmb_nme = { 'pap_anv_1500' };

cfg.sig_loc = { { 'pap_ovr_stt' 'pap_ovr_stt' 'pap_ovr_stt' 'pap_anv_stt' 'pap_anv_stt' 'pap_anv_stt' } };
           
cfg.sig_chn = { [ 1 2 3 1 2 3 ] };

% Table Combination
cfg.run_tbl     = [1];

cfg.tbl_cmb_nme = { { 'VisualSelective' 'AuditorySelective'   'ResponseSelective'} };
               
cfg.tbl_cmb_rul = { { '1 & 4' '2 & 5' '3 & 6'} };
               
cfg.tbl_sum = { [] };

% Plot combination
cfg.run_plt     = [ 0 ];

% First Pass combinations
cfg.run_fst     = [ 0 ];

% Time combinations
cfg.run_tme     = [ 0 ];
    
% Line Plot combinations
cfg.run_lne     = [ 0 ];
                
% Run
cfg.dat_nme     = '_overall_data';
cfg.sbj_clr_fld = fcfg.clr_fld;
cfg.sbj_nme = fcfg.sbj_nme;

cfg.fle_out_pth = fcfg.dat_fld;

mmil_combo_effects2(cfg);

% Save for overall analysis
for iC = 1:numel(cfg.cmb_nme); cfg.clr_fld = fcfg.clr_fld; cfg.iC = iC; mmil_ovr_ana_sve(cfg) ; end

end