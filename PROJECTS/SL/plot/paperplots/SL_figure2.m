function SL_figure2

%% Make Plot
% rsp
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
cfg.sve_loc = 'figure2';
cfg.tsk     = 'SL';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_rsp_1500';
cfg.eff_clm = [1 2];
cfg.eff_col = { {'dark purple' 'purple' 'neon purple'} {'dark purple' 'purple' 'neon purple'} };
cfg.top_pct = 0.50;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_include_plot(cfg)

% anv
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
cfg.sve_loc = 'figure2';
cfg.tsk     = 'SL';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_anv_1500';
cfg.eff_clm = [ 1 2 ];
cfg.eff_col = { {'dark purple' 'purple' 'neon purple'} {'dark purple' 'purple' 'neon purple'} };
cfg.top_pct = 0.50;
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };
mmil_include_plot(cfg)

% Colorbar
pcfg = [];
pcfg.col_map = {'dark purple' 'purple' 'neon purple'};
pcfg.out_dir = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical' '/' 'manuscript' '/' 'figure2' '/' ];
pcfg.col_bar = [0 0.5];
mmil_color_bar(pcfg)

end





