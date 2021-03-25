%% SL Vis - Stats
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 1.200];
fcfg.stt_toi = [ 0.050 0.500];

fcfg.bse_lne = [-0.200 0];
fcfg.sig__ms = 50;

fcfg.eve 	 = [ 1              2               3 ];
fcfg.sig_nme = { 'pap_anv_1500' 'pap_anv_1500'  'pap_anv_1500'};
fcfg.sig_col = [ 1              1               1  ];
fcfg.eve_nme = { 'Match'        'Mismatch'      'FalseFont'};
fcfg.eve_col = { 'green'        'dark yellow'   'reddish grey' };

% Fusiform
fcfg.frm_loc = { 'caudal-fusiform_L' };
fcfg.frm_nme = { 'CaudalFusiform'};
fcfg.end_loc = { {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Precentral'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'CaudalFusiform'];
PLV_sig(fcfg)

% Precentral
fcfg.frm_loc = { 'middle-precentral_L' 'inferior-precentral_L' };
fcfg.frm_nme = { 'Precentral'};
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'Precentral'];
PLV_sig(fcfg)

% STG
fcfg.frm_loc = { 'caudal-STG_L' 'middle-STG_L' };
fcfg.frm_nme = { 'STG' };
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                 'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'STG'];
PLV_sig(fcfg)

% Lateral Occipital
fcfg.frm_loc = { 'lateraloccipital_L' };
fcfg.frm_nme = { 'LateralOccipital' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'LateralOccipital'];
PLV_sig(fcfg)

% MTG
fcfg.frm_loc = { 'caudal-MTG_L' 'middle-MTG_L' };
fcfg.frm_nme = { 'MTG' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'MTG'];
PLV_sig(fcfg)

% Supramarginal
fcfg.frm_loc = { 'supramarginal_L' };
fcfg.frm_nme = { 'Supramarginal' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'Supramarginal'];
PLV_sig(fcfg)

% ParsTriangularis
fcfg.frm_loc = { 'parstriangularis_L'  };
fcfg.frm_nme = { 'ParsTrinagularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'ParsTriangularis'];
PLV_sig(fcfg)

% ParsOpercularis
fcfg.frm_loc = { 'parsopercularis_L' };
fcfg.frm_nme = { 'ParsOpercularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parstriangularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsTrinagularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/'  'ParsOpercularis'];
PLV_sig(fcfg)

%% SL Vis - Aud
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 1.200];
fcfg.stt_toi = [ 0.500 0.950];

fcfg.bse_lne = [-0.200 0];
fcfg.sig__ms = 50;

fcfg.eve 	 = [ 1              2               3 ];
fcfg.sig_nme = { 'pap_anv_1500' 'pap_anv_1500'  'pap_anv_1500'};
fcfg.sig_col = [ 1              1               1  ];
fcfg.eve_nme = { 'Match'        'Mismatch'      'FalseFont'};
fcfg.eve_col = { 'green'        'dark yellow'   'reddish grey' };

% Fusiform
fcfg.frm_loc = { 'caudal-fusiform_L' };
fcfg.frm_nme = { 'CaudalFusiform'};
fcfg.end_loc = { {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Precentral'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'CaudalFusiform'];
PLV_sig(fcfg)

% Precentral
fcfg.frm_loc = { 'middle-precentral_L' 'inferior-precentral_L' };
fcfg.frm_nme = { 'Precentral'};
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'Precentral'];
PLV_sig(fcfg)

% STG
fcfg.frm_loc = { 'caudal-STG_L' 'middle-STG_L' };
fcfg.frm_nme = { 'STG' };
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                 'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'STG'];
PLV_sig(fcfg)

% Lateral Occipital
fcfg.frm_loc = { 'lateraloccipital_L' };
fcfg.frm_nme = { 'LateralOccipital' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'LateralOccipital'];
PLV_sig(fcfg)

% MTG
fcfg.frm_loc = { 'caudal-MTG_L' 'middle-MTG_L' };
fcfg.frm_nme = { 'MTG' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'MTG'];
PLV_sig(fcfg)

% Supramarginal
fcfg.frm_loc = { 'supramarginal_L' };
fcfg.frm_nme = { 'Supramarginal' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'Supramarginal'];
PLV_sig(fcfg)

% ParsTriangularis
fcfg.frm_loc = { 'parstriangularis_L'  };
fcfg.frm_nme = { 'ParsTrinagularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'ParsTriangularis'];
PLV_sig(fcfg)

% ParsOpercularis
fcfg.frm_loc = { 'parsopercularis_L' };
fcfg.frm_nme = { 'ParsOpercularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parstriangularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsTrinagularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Aud' '/'  'ParsOpercularis'];
PLV_sig(fcfg)

%% Put together Stats
clear plv_hld

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/';

nme     = { 'CaudalFusiform' 'Precentral' }; % { 'LateralOccipital' 'Fusiform' 'Supramarginal' 'MTG' 'STG' 'Precentral' 'ParsOpercularis' }; %'ParsTriangularis' }; % 
eve_nme = { 'VisualMatch' 'AuditoryMatch' };

dat_hld = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'Stats_Vis' '/' ];

for iN = 1:numel(nme)
    for iEN = 1:numel(eve_nme)
    
        sig_hld = mmil_readtext([ dat_hld '/' nme{iN} '/' eve_nme{iEN} '/' eve_nme{iEN} '_total' ]);
    
        if size(sig_hld,2)>7
            for iR = 1:size(sig_hld,1)
                
                plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){iEN} = [];
                
                ele_nme = strsplit(sig_hld{iR,end},' /');
                
                %
                sbj_loc = string_find(ele_nme,'SL');
                for iEL = 1:numel(sbj_loc)-1
                    
                    sbj_nme = strfind(ele_nme{sbj_loc(iEL)},'SL');
                    sbj_nme = ele_nme{sbj_loc(iEL)}(sbj_nme-6:sbj_nme-2);
                    
                    ele_loc{iEL} = sbj_loc(iEL)+1:sbj_loc(iEL+1)-1;
                    if ~isempty(ele_loc{iEL})
                        plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){iEN} = [plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){iEN} strcat([sbj_nme '-'],ele_nme(ele_loc{iEL}))];
                    end
                end
                
                %
                if isempty(iEL)
                    iEL = 1;
                else
                    iEL = iEL + 1;
                end
                sbj_nme = strfind(ele_nme{sbj_loc(iEL)},'SL');
                sbj_nme = ele_nme{sbj_loc(iEL)}(sbj_nme-6:sbj_nme-2);
                
                ele_loc{iEL} = sbj_loc(iEL)+1:numel(ele_nme);
                if ~isempty(ele_loc{iEL})
                    plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){iEN} = [plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){iEN} strcat([sbj_nme '-'],ele_nme(ele_loc{iEL}))];
                end
                
            end
        else
            for iR = 1:size(sig_hld,1)
                plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){iEN} = [];
            end
        end
    end  
    for iR = 1:size(sig_hld,1)
        plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){4} = sig_hld{iR,7};
    end
end

%% Plot
clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';

%
vis_cnt = [];
aud_cnt = [];

fst_nme = fieldnames(plv_hld);
for iN = 1:numel(fst_nme)
    scd_nme = fieldnames(plv_hld.(fst_nme{iN})); scd_nme(end) = [];
    for iS = 1:numel(scd_nme)
        vis_cnt = [vis_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){1}];
        aud_cnt = [aud_cnt plv_hld.(fst_nme{iN}).(scd_nme{iS}){2}];       
    end    
end

%
vis_cnt = unique(vis_cnt);

vis_cnt_nme = cellfun(@(x) x(1:5),vis_cnt,'uni',0);

vis_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),vis_cnt,'uni',0),'uni',0);
for iE = 1:numel(vis_cnt_par); 
    vis_cnt_par{iE}{1} = [vis_cnt_nme{iE} '_SL_' vis_cnt_par{iE}{1}]; vis_cnt_par{iE}{2} = [vis_cnt_nme{iE} '_SL_' vis_cnt_par{iE}{2}]; 
end

vis_cnt_ele = unique([vis_cnt_par{:}]);

vis_ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_anv_1500' '/' 'subjects' '/' 'total' '/' 'pap_anv_1500' '_plt']);
vis_ovr_txt = vis_ovr_txt(2:end,:);
vis_ovr_txt = vis_ovr_txt(logical(cell2mat(vis_ovr_txt(:,3))),2);

%
aud_cnt = unique(aud_cnt);

aud_cnt_nme = cellfun(@(x) x(1:5),aud_cnt,'uni',0);

aud_cnt_par = cellfun(@(x) strsplit(x,'--'),cellfun(@(x) x(7:end),aud_cnt,'uni',0),'uni',0);
for iE = 1:numel(aud_cnt_par); aud_cnt_par{iE}{1} = [aud_cnt_nme{iE} '_SL_' aud_cnt_par{iE}{1}]; aud_cnt_par{iE}{2} = [aud_cnt_nme{iE} '_SL_' aud_cnt_par{iE}{2}]; end

aud_cnt_ele = unique([aud_cnt_par{:}]);

aud_ovr_txt = mmil_readtext([clr_fld '/' 'sig_chn' '/' 'hgp' '/' 'ecog' '/' 'split' '/' 'pap_anv_1500' '/' 'subjects' '/' 'total' '/' 'pap_anv_1500' '_plt']);
aud_ovr_txt = aud_ovr_txt(2:end,:);
aud_ovr_txt = aud_ovr_txt(logical(cell2mat(aud_ovr_txt(:,4))),2);

%
ovr_txt = unique([vis_ovr_txt ; aud_ovr_txt]);

% PLOT
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']          [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog'] [clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = {aud_cnt_ele vis_cnt_ele};
cfg.sel_lbl = {'Visual'   'Auditory'};
cfg.col     = {rgb('blue') rgb('red')};
cfg.nsl_col = {rgb('dark purple')};

cfg.par_ele   = {aud_cnt_par vis_cnt_par};
cfg.par_ele_col = {rgb('bluish grey')+[0.1 0.1 0.2] rgb('reddish grey')+[0.15 0.075 0.075]};
cfg.par_ele_lne = [0.55                0.55];

cfg.all_ele   = ovr_txt; %ovr_txt;



cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' 'parahippocampal' 'entorhinal' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' 'bankssts' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [clr_fld '/' 'manuscript' '/' 'figure8' '/' 'Mixed'];
cfg.sve_pre   = ['middle_pic_visual'];
cfg.sep_str   = [','];

cfg.rad       = 2.35;

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);


























