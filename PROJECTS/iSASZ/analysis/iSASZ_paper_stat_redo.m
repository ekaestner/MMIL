% Activation
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_act' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_act' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_act' };
cfg.typ     = { 'lfp' }; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.typ     = { 'hgp' }; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Repetition
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_rep' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_rep' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_rep' };
cfg.typ     = { 'lfp' };  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.typ     = { 'hgp' }; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N400
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_n400' };
cfg.typ     = { 'lfp' }; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.typ     = { 'hgp' }; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lexical
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_lex' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_lex' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_lex' };
cfg.typ     = { 'lfp' }; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg.typ     = { 'hgp' }; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmil_overall_analysis(cfg); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


