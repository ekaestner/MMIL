%% Fusiform - SL
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
fcfg.sbj_nme = 'NY226_FW';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 1.250];
fcfg.stt_toi = [ 0.050 0.500];

fcfg.bse_lne = [-0.300 0];
fcfg.sig__ms = 100;

fcfg.eve 	 = [ 1 ];
fcfg.eve_nme = {'Match'};
fcfg.eve_col = {'red'};

fcfg.frm_loc = { 'caudal-fusiform_L' };
fcfg.frm_nme = { 'Caudal Fusiform'};
fcfg.end_loc = { {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-STG_L' 'middle-STG_L'} {'parsopercularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'Precentral'                                    'STG'                           'Pars Opercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Fusiform_Text'];
PLV_sig(fcfg)

fcfg.stt_toi = [ 0.450 0.950];
fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Fusiform_Voice'];
PLV_sig(fcfg)

%% Precentral - SL
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
fcfg.sbj_nme = 'NY226_FW';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 1.250];
fcfg.stt_toi = [ 0.050 0.500];

fcfg.bse_lne = [-0.300 0];
fcfg.sig__ms = 100;

fcfg.eve 	 = [ 1 ];
fcfg.eve_nme = {'Match'};
fcfg.eve_col = {'red'};

fcfg.frm_loc = { 'middle-precentral_L' 'inferior-precentral_L' };
fcfg.frm_nme = { 'Precentral'};
fcfg.end_loc = { {'caudal-fusiform_L'} {'caudal-STG_L' 'middle-STG_L'} {'parsopercularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'Caudal Fusiform'      'STG'                           'Pars Opercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Precentral_Text'];
PLV_sig(fcfg)

fcfg.stt_toi = [ 0.450 0.950];
fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Precentral_Voice'];
PLV_sig(fcfg)

%% STG - SL
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
fcfg.sbj_nme = 'NY226_FW';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 1.250];
fcfg.stt_toi = [ 0.050 0.500];

fcfg.bse_lne = [-0.300 0];
fcfg.sig__ms = 100;

fcfg.eve 	 = [ 1 ];
fcfg.eve_nme = {'Match'};
fcfg.eve_col = {'red'};

fcfg.frm_loc = {'caudal-STG_L' 'middle-STG_L' };
fcfg.frm_nme = { 'STG'};
fcfg.end_loc = {{'caudal-fusiform_L'} {'middle-precentral_L' 'inferior-precentral_L'}  {'parsopercularis_L'} {'supramarginal_L'}};
fcfg.end_nme = { 'Caudal Fusiform'     'Precentral'                                    'Pars Opercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'STG_Text'];

fcfg.stt_toi = [ 0.450 0.950];
fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'STG_Voice'];
PLV_sig(fcfg)

%% Pars Opercularis - SL
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
fcfg.sbj_nme = 'NY226_FW';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 1.250];
fcfg.stt_toi = [ 0.050 0.500];

fcfg.bse_lne = [-0.300 0];
fcfg.sig__ms = 100;

fcfg.eve 	 = [ 1 ];
fcfg.eve_nme = {'Match'};
fcfg.eve_col = {'red'};

fcfg.frm_loc = { 'parsopercularis_L'};
fcfg.frm_nme = { 'Pars Opercularis' };
fcfg.end_loc = {{'caudal-fusiform_L'} {'middle-precentral_L' 'inferior-precentral_L'}  {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'}};
fcfg.end_nme = { 'Caudal Fusiform'     'Precentral'                                    'STG'                            'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'ParsOpercularis_Text'];
PLV_sig(fcfg)

fcfg.stt_toi = [ 0.450 0.950];
fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'ParsOpercularis_Voice'];
PLV_sig(fcfg)

%% Supramarginal - SL
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';
fcfg.sbj_nme = 'NY226_FW';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 1.250];
fcfg.stt_toi = [ 0.050 0.500];

fcfg.bse_lne = [-0.300 0];
fcfg.sig__ms = 100;

fcfg.eve 	 = [ 1 ];
fcfg.eve_nme = {'Match'};
fcfg.eve_col = {'red'};

fcfg.frm_loc = {'caudal-STG_L' 'middle-STG_L' };
fcfg.frm_nme = { 'STG'};
fcfg.end_loc = {{'caudal-fusiform_L'} {'middle-precentral_L' 'inferior-precentral_L'}  {'parsopercularis_L'} {'supramarginal_L'}};
fcfg.end_nme = { 'Caudal Fusiform'     'Precentral'                                    'Pars Opercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'STG_Text'];

fcfg.stt_toi = [ 0.450 0.950];
fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'STG_Voice'];
PLV_sig(fcfg)
