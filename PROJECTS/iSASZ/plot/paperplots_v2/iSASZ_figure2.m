ejk_vis = hsv2rgb([ .3  .03 .45   ; ...
                    .05 .1  .50   ; ...
                    .04 .2  .6 ; ...
                    .03 .45 .70 ; ...
                    .01 .7  .80 ; ...
                    .97 .95 .95 ]);

ejk_aud = hsv2rgb([ .3  .03 .45   ; ...
                    .66 .1  .50   ; ...
                    .65 .2  .6 ; ...
                    .64 .45 .70 ; ...
                    .62 .7  .80 ; ...
                    .58 .95 .95 ]);
                
ejk_bim = hsv2rgb([ .3  .03 .4   ; ...
                    .81 .1  .5   ; ...
                    .80 .2  .5 ; ...
                    .79 .45 .65 ; ...
                    .77 .7  .85 ; ...
                    .73 .95 .95 ]);
                
%% Auditory Responsive HGP
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.sve_loc = 'NEW_FIGS';
cfg.tsk     = 'iSASZ';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_aud_act';
cfg.eff_clm = [1];
cfg.eff_col = { ejk_aud };
cfg.top_pct = 0.60;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_include_plot(cfg)

% % Colorbar
% pcfg = [];
% pcfg.col_map = {'dark blue' 'blue' 'bright blue'};
% pcfg.out_dir = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'manuscript' '/' 'NEW_FIGS' '/' ];
% pcfg.sve_pre = 'auditory';
% pcfg.col_bar = [0 0.60];
% mmil_color_bar(pcfg)

%% Visual Responsive HGP
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.sve_loc = 'NEW_FIGS';
cfg.tsk     = 'iSASZ';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_vis_act';
cfg.eff_clm = [1];
cfg.eff_col = { ejk_vis };
cfg.top_pct = 0.60;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_include_plot(cfg)

% % Colorbar
% pcfg = [];
% pcfg.col_map = {'dark red' 'red' 'bright red'};
% pcfg.out_dir = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'manuscript' '/' 'NEW_FIGS' '/' ];
% pcfg.sve_pre = 'visual';
% pcfg.col_bar = [0 0.60];
% mmil_color_bar(pcfg)

%% Bi-modal Responsive HGP
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.sve_loc = 'NEW_FIGS';
cfg.tsk     = 'iSASZ';

cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_bim_act';
cfg.eff_clm = [1];

% cfg.scl_typ = 'hgp';
% cfg.scl_nme = {'pap_vis_act' 'pap_aud_act'};
% cfg.scl_clm = [1             1];

cfg.eff_col = { ejk_bim };
cfg.top_pct = 0.30;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' 'Inferior Parietal' ......
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_include_plot(cfg)

% Colorbar
% pcfg = [];
% pcfg.col_map = {'dark purple' 'purple' 'bright purple'};
% pcfg.out_dir = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'manuscript' '/' 'NEW_FIGS' '/' ];
% pcfg.sve_pre = 'bimodal';
% pcfg.col_bar = [0 0.75];
% mmil_color_bar(pcfg)

%% Auditory Responsive LFP
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.sve_loc = 'NEW_FIGS';
cfg.tsk     = 'iSASZ';
cfg.eff_typ = 'lfp';
cfg.eff_nme = 'pap_aud_act';
cfg.eff_clm = [1];
cfg.eff_col = { ejk_aud };
cfg.top_pct = 1.00;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_include_plot(cfg)

%% Visual Responsive LFP
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.sve_loc = 'NEW_FIGS';
cfg.tsk     = 'iSASZ';
cfg.eff_typ = 'lfp';
cfg.eff_nme = 'pap_vis_act';
cfg.eff_clm = [1];
cfg.eff_col = { ejk_vis };
cfg.top_pct = 1.00;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_include_plot(cfg)

%% Bi-modal Responsive LFP
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.sve_loc = 'NEW_FIGS';
cfg.tsk     = 'iSASZ';

cfg.eff_typ = 'lfp';
cfg.eff_nme = 'pap_bim_act';
cfg.eff_clm = [1];

%cfg.scl_typ = 'lfp';
%cfg.scl_nme = {'pap_vis_act' 'pap_aud_act'};
%cfg.scl_clm = [1             1];

cfg.eff_col = { ejk_bim };
cfg.top_pct = 0.50;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_include_plot(cfg)
