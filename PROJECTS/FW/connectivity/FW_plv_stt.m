%% iSASZ - Stats
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 0.850];
fcfg.stt_toi = [ 0.050 0.650];

fcfg.bse_lne = [-0.200 0];
fcfg.sig__ms = 50;

fcfg.eve 	 = [ 1             4];
fcfg.sig_nme = { 'pap_wrd_600' 'pap_con_600'};
fcfg.sig_col = [ 3             1 ];
fcfg.eve_nme = { 'Word'        'FalseFont'};
fcfg.eve_col = { 'red'         'reddish grey' };

% Fusiform
fcfg.frm_loc = { 'caudal-fusiform_L' 'middle-fusiform_L' };
fcfg.frm_nme = { 'Fusiform'};
fcfg.end_loc = { {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Precentral'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'Fusiform'];
PLV_sig(fcfg)

% Precentral
fcfg.frm_loc = { 'middle-precentral_L' 'inferior-precentral_L' };
fcfg.frm_nme = { 'Precentral'};
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'Precentral'];
PLV_sig(fcfg)

% STG
fcfg.frm_loc = { 'caudal-STG_L' 'middle-STG_L' };
fcfg.frm_nme = { 'STG' };
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                 'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'STG'];
PLV_sig(fcfg)

% Lateral Occipital
fcfg.frm_loc = { 'lateraloccipital_L' };
fcfg.frm_nme = { 'LateralOccipital' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'LateralOccipital'];
PLV_sig(fcfg)

% MTG
fcfg.frm_loc = { 'caudal-MTG_L' 'middle-MTG_L' };
fcfg.frm_nme = { 'MTG' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'MTG'];
PLV_sig(fcfg)

% Supramarginal
fcfg.frm_loc = { 'supramarginal_L' };
fcfg.frm_nme = { 'Supramarginal' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'Supramarginal'];
PLV_sig(fcfg)

% ParsTriangularis
fcfg.frm_loc = { 'parstriangularis_L'  };
fcfg.frm_nme = { 'ParsTrinagularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'ParsTriangularis'];
PLV_sig(fcfg)

% ParsOpercularis
fcfg.frm_loc = { 'parsopercularis_L' };
fcfg.frm_nme = { 'ParsOpercularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parstriangularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsTrinagularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/'  'ParsOpercularis'];
PLV_sig(fcfg)

%% Put together Stats
clear plv_hld

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/';

sig_nme = 'act_dis_650_hlf_sig_00005';
nme     = { 'LateralOccipital' 'Fusiform' 'Supramarginal' 'MTG' 'STG' 'Precentral' 'ParsOpercularis' 'ParsTriangularis' }; % 
eve_nme = { 'Word' };

dat_hld = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure7' '/' 'Stats' '/' sig_nme];

for iN = 1:numel(nme)
    for iEN = 1:numel(eve_nme)
    
        sig_hld = mmil_readtext([ dat_hld '/' nme{iN} '/' eve_nme{iEN} '/' eve_nme{iEN} '_total' ]);
    
        if size(sig_hld,2)>7
            for iR = 1:size(sig_hld,1)
                
                plv_hld.(sig_hld{iR,1}).(sig_hld{iR,2}){iEN} = [];
                
                ele_nme = strsplit(sig_hld{iR,end},' /');
                
                %
                sbj_loc = string_find(ele_nme,'SA_SZ');
                for iEL = 1:numel(sbj_loc)-1
                    
                    sbj_nme = strfind(ele_nme{sbj_loc(iEL)},'SA_SZ');
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
                sbj_nme = strfind(ele_nme{sbj_loc(iEL)},'SA_SZ');
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


fst_nme = fieldnames(plv_hld);
for iFS = 1:numel(fst_nme)
    
    scd_nme = fieldnames(plv_hld.(fst_nme{iFS}));
    for iSC = 1:numel(scd_nme)
        plv_hld.(fst_nme{iFS}).(scd_nme{iSC}){3} = intersect(plv_hld.(fst_nme{iFS}).(scd_nme{iSC}){1},plv_hld.(fst_nme{iFS}).(scd_nme{iSC}){2});
        plv_hld.(fst_nme{iFS}).tot{iSC,1} = scd_nme{iSC};
        plv_hld.(fst_nme{iFS}).tot{iSC,2} = numel(plv_hld.(fst_nme{iFS}).(scd_nme{iSC}){1});
        plv_hld.(fst_nme{iFS}).tot{iSC,3} = numel(plv_hld.(fst_nme{iFS}).(scd_nme{iSC}){2});
        plv_hld.(fst_nme{iFS}).tot{iSC,4} = numel(plv_hld.(fst_nme{iFS}).(scd_nme{iSC}){3});
        plv_hld.(fst_nme{iFS}).tot{iSC,5} = plv_hld.(fst_nme{iFS}).(scd_nme{iSC}){4};
    end
end

plv_hld.LateralOccipital.tot
plv_hld.Fusiform.tot
plv_hld.Supramarginal.tot
plv_hld.STG.tot
plv_hld.MTG.tot
plv_hld.Precentral.tot
plv_hld.ParsOpercularis.tot
plv_hld.ParsTrinagularis.tot
