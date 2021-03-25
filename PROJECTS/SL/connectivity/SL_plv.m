clear; clc;

pcfg = [];
pcfg.prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';

pcfg.eff_typ = {'hgp'};
pcfg.eff_nme = 'Selectivity';
pcfg.eff_rul = {'all'};

pcfg.toi = [-0.2 1.2];
pcfg.bse_lne = [-0.200 0];
pcfg.eve = [1 2 3 4];
pcfg.eve_nme = {'Match' 'MisMatch' 'False-Font' 'Noise-Vocoding'};
pcfg.col_bar = [0 0.6];

pcfg.win = 0.050;

pcfg.foi = [4 8];

pcfg.mve_win = 0.010;
pcfg.mve_tme = [0.100 1.125];

pcfg.lbl     = {'Visual'      'Auditory'    ''};
pcfg.lbl_tme = {[0.000 0.450] [0.450 0.900] [0.900 1.150]};

pcfg.nde_nme = { { 'caudal-fusiform_L'     'caudal-fusiform_R'     'lateraloccipital_L'  'lateraloccipital_R' } ...
                { 'inferior-precentral_L' 'inferior-precentral_R' 'middle-precentral_L' 'middle-precentral_R' } ...
                { 'caudal-STG_L'          'caudal-STG_R'          'middle-STG_L'        'middle-STG_R' 'supramarginal_L' 'supramarginal_R'} };
pcfg.nde_col = { 'red'    ...
                'yellow' ...
                'blue'   };

pcfg.edg_nme = { [1 2]    [1 3]    [2 3] };
pcfg.edg_col = { 'orange' 'purple' 'aquamarine'};

for iPL = 1;
  pcfg.sbj_num     = iPL;
  PLV_Play3(pcfg)
end