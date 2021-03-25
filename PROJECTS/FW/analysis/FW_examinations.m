%% Overall Location Plot Figure
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
mmil_ovr_ele_loc(cfg);

%% Stimulation Results
cfg = [];
cfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/FW/clerical/';
cfg.tsk = 'FW';
cfg.grp = { 'Motor' 	  'Language'};
cfg.col = { rgb('orange') rgb('blue') };
mmil_stimulation_plot(cfg);

%% IFG Grainsize
cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'GrainSize';
cfg.inc_eve = {'parsopercularis' 'parstriangularis' 'parsorbitalis'};
cfg.out_nme = 'IFG_GrainSize';
mmil_region_examination(cfg)

%% Sensory-Motor Integration Areas - GrainSize
cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'GrainSize';
cfg.inc_eve = {'inferior-postcentral' 'middle-postcentral' 'inferior-precentral' 'middle-precentral' 'supramarginal' 'inferiorparietal' 'parsopercularis'};
cfg.out_nme = 'SnsMot_GrainSize';
mmil_region_examination(cfg)

%% Ventral Route Areas - GrainSize
cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'GrainSize';
cfg.inc_eve = {'lingual' 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' 'parahippocampal' 'entorhinal' 'temporalpole' 'parstriangularis' 'parsorbitalis'};
cfg.out_nme = 'VntRte_GrainSize';
mmil_region_examination(cfg)

%% Auditory Route Areas - GrainSize
cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp' 'lfp'};
cfg.ana_nme = 'GrainSize';
cfg.inc_eve = { 'caudal-MTG' 'middle-MTG' 'rostral-MTG' 'caudal-STG' 'middle-STG' 'rostral-STG'};
cfg.out_nme = 'AudRte_GrainSize';
mmil_region_examination(cfg)

%% Sensory-Motor Integration Areas - Unfamiliar
cfg = [];
cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
cfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
cfg.typ     = {'hgp'};
cfg.ana_nme = 'UnfamiliarStimuli';
cfg.inc_eve = {'inferior-postcentral' 'middle-postcentral' 'inferior-precentral' 'middle-precentral' 'supramarginal' 'inferiorparietal' 'parsopercularis'};
cfg.out_nme = 'SnsMot_Unfamiliar';
mmil_region_examination(cfg)
