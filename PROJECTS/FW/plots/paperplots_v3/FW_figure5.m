clear; clc;

%% Repetition
% 
cfg = [];
cfg.chn_crr = 1;
cfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/FW/clerical';
cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = 'pap_lex_600';

cfg.eff_clm = 1;
cfg.top_pct = 0.25;
cfg.eff_col = {[ rgb('medium grey')-0.075 ; rgb('medium grey')-0.15 ; rgb('dark orange') ; rgb('bright orange')-[.15 -.075 0] ; rgb('yellowish orange')-[0 0.075 0]]};
cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' ...
                'Inferior Precentral' 'Middle Precentral' ...
                'Inferior Parietal' 'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' };
cfg.sve_loc = 'NewFigs_v3';
mmil_include_plot(cfg)

pcfg = [];
pcfg.col_map = [ rgb('medium grey')-0.075 ; rgb('medium grey')-0.15 ; rgb('dark orange') ; rgb('bright orange')-[.15 -.075 0] ; rgb('yellowish orange')-[0 0.075 0]];
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/NEW_FIGS';
pcfg.sve_pre = 'repetition';
mmil_color_bar(pcfg)