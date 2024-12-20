clear; clc;

%%
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Caudal ITG' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

mmil_ovr_ele_loc(cfg);

% Colorbar
ttt = summer(16);

pcfg = [];
pcfg.col_map = [ rgb('greenish grey') ; ttt([10 14],:)]; %[ rgb('light brown') ; [0.65 0.47 0.20]  ; [0.85 0.60 0.08] ]; %{ 'light brown' 'copper'};
pcfg.col_bar = [0 0.5];
pcfg.out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/manuscript/NEW_FIGS';
pcfg.sve_pre = 'total_electrodes';
mmil_color_bar(pcfg)