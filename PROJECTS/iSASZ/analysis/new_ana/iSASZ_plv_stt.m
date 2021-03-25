%% iSASZ - 
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 0.850];
fcfg.stt_toi = [ 0.050 0.650];

fcfg.bse_lne = [-0.300 0];
fcfg.sig__ms = 50;

fcfg.eve 	 = [ 1        3];
fcfg.eve_nme = { 'Visual' 'Auditory'};
fcfg.eve_col = { 'red'    'blue' };

% Fusiform
fcfg.frm_loc = { 'caudal-fusiform_L' 'middle-fusiform_L' };
fcfg.frm_nme = { 'Fusiform'};
fcfg.end_loc = { {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Precentral'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Fusiform'];
PLV_sig(fcfg)

% Precentral - iSASZ
fcfg.frm_loc = { 'middle-precentral_L' 'inferior-precentral_L' };
fcfg.frm_nme = { 'Precentral'};
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Precentral'];
PLV_sig(fcfg)

% STG - iSASZ
fcfg.frm_loc = { 'caudal-STG_L' 'middle-STG_L' };
fcfg.frm_nme = { 'STG' };
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                 'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'STG'];
PLV_sig(fcfg)

% Lateral Occipital - iSASZ
fcfg.frm_loc = { 'lateraloccipital_L' };
fcfg.frm_nme = { 'LateralOccipital' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'LateralOccipital'];
PLV_sig(fcfg)

% MTG - iSASZ
fcfg.frm_loc = { 'caudal-MTG_L' 'middle-MTG_L' };
fcfg.frm_nme = { 'MTG' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'MTG'];
PLV_sig(fcfg)

% Supramarginal - iSASZ
fcfg.frm_loc = { 'supramarginal_L' };
fcfg.frm_nme = { 'Supramarginal' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Supramarginal'];
PLV_sig(fcfg)

% ParsTriangularis - iSASZ
fcfg.frm_loc = { 'parstriangularis_L'  };
fcfg.frm_nme = { 'ParsTrinagularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'ParsTriangularis'];
PLV_sig(fcfg)

% ParsOpercularis - iSASZ
fcfg.frm_loc = { 'parsopercularis_L' };
fcfg.frm_nme = { 'ParsOpercularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parstriangularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsTrinagularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'ParsOpercularis'];
PLV_sig(fcfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iSASZ - Stats
fcfg = [];

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';

fcfg.foi     = [6 12];
fcfg.plt_toi = [ 0.050 0.850];
fcfg.stt_toi = [ 0.050 0.650];

fcfg.bse_lne = [-0.300 0];
fcfg.sig__ms = 50;

fcfg.eve 	 = [ 1        3];
fcfg.sig_nme = { 'pap_vis_act' 'pap_aud_act'};
fcfg.sig_col = [ 1             1 ];
fcfg.eve_nme = { 'Visual' 'Auditory'};
fcfg.eve_col = { 'red'    'blue' };

% Fusiform
fcfg.frm_loc = { 'caudal-fusiform_L' 'middle-fusiform_L' };
fcfg.frm_nme = { 'Fusiform'};
fcfg.end_loc = { {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Precentral'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Fusiform'];
PLV_sig(fcfg)

% Precentral - iSASZ
fcfg.frm_loc = { 'middle-precentral_L' 'inferior-precentral_L' };
fcfg.frm_nme = { 'Precentral'};
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'caudal-STG_L' 'middle-STG_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                    'MTG'                           'STG'                           'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Precentral'];
PLV_sig(fcfg)

% STG - iSASZ
fcfg.frm_loc = { 'caudal-STG_L' 'middle-STG_L' };
fcfg.frm_nme = { 'STG' };
fcfg.end_loc = { {'lateraloccipital_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'LateralOccipital'     'Fusiform'                                 'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'STG'];
PLV_sig(fcfg)

% Lateral Occipital - iSASZ
fcfg.frm_loc = { 'lateraloccipital_L' };
fcfg.frm_nme = { 'LateralOccipital' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'MTG'                           'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'LateralOccipital'];
PLV_sig(fcfg)

% MTG - iSASZ
fcfg.frm_loc = { 'caudal-MTG_L' 'middle-MTG_L' };
fcfg.frm_nme = { 'MTG' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'supramarginal_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'Supramarginal'     'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'MTG'];
PLV_sig(fcfg)

% Supramarginal - iSASZ
fcfg.frm_loc = { 'supramarginal_L' };
fcfg.frm_nme = { 'Supramarginal' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'parstriangularis_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'ParsTrinagularis' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'Supramarginal'];
PLV_sig(fcfg)

% ParsTriangularis - iSASZ
fcfg.frm_loc = { 'parstriangularis_L'  };
fcfg.frm_nme = { 'ParsTrinagularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parsopercularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsOpercularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'ParsTriangularis'];
PLV_sig(fcfg)

% ParsOpercularis - iSASZ
fcfg.frm_loc = { 'parsopercularis_L' };
fcfg.frm_nme = { 'ParsOpercularis' };
fcfg.end_loc = { {'caudal-STG_L' 'middle-STG_L'} {'caudal-fusiform_L' 'middle-fusiform_L'} {'lateraloccipital_L'} {'middle-precentral_L' 'inferior-precentral_L'} {'caudal-MTG_L' 'middle-MTG_L'} {'parstriangularis_L'} {'supramarginal_L'} };
fcfg.end_nme = { 'STG'                            'Fusiform'                                'LateralOccipital'     'Precentral'                                   'MTG'                            'ParsTrinagularis'    'Supramarginal' };

fcfg.out_dir = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/'  'ParsOpercularis'];
PLV_sig(fcfg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear plv_hld

fcfg.dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';

% act_550_ful_sig_00005
% act_650_hlf_sig_00005
% act_650_hlf_sig_0001
% act_dis_650_hlf_sig_00005
% act_850_ful_sig_00005
% act_850_hlf_sig_00005

% rep_650_hlf_sig_00005


sig_nme = 'act_dis_650_hlf_sig_00005';
nme     = { 'LateralOccipital' 'Fusiform' 'Supramarginal' 'MTG' 'STG' 'Precentral' 'ParsOpercularis' 'ParsTriangularis' }; % 
eve_nme = { 'Visual' 'Auditory' };

dat_hld = [fcfg.dat_hld '/' 'clerical' '/' 'manuscript' '/' 'figure8' '/' 'NetworkWide' '/' sig_nme];

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

%% PUT TOGETHER ACTUAL NUMBERS
% tot_cnt = 0;

fst_nme = fieldnames(plv_hld);

tot_cnt.visual(1,1:9) = [ cell(1,1) fst_nme'];
tot_cnt.visual(1:9,1) = [ cell(1,1) ; fst_nme];

tot_cnt.auditory(1,1:9) = [ cell(1,1) fst_nme'];
tot_cnt.auditory(1:9,1) = [ cell(1,1) ; fst_nme];

tot_cnt.bimodal(1,1:9) = [ cell(1,1) fst_nme'];
tot_cnt.bimodal(1:9,1) = [ cell(1,1) ; fst_nme];

tot_cnt.tot(1,1:9) = [ cell(1,1) fst_nme'];
tot_cnt.tot(1:9,1) = [ cell(1,1) ; fst_nme];

for iONE = 1:numel(fst_nme)
    
    scd_nme = fieldnames(plv_hld.(fst_nme{iONE}));
    for iTWO = 1:numel(scd_nme)-1
        
        iR = find(strcmpi(tot_cnt.visual(1:9,1),fst_nme{iONE}));
        iC = find(strcmpi(tot_cnt.visual(1,1:9),scd_nme{iTWO}));
        
        if iR == iC
            tot_cnt.visual{iR,iC}   = NaN;
            tot_cnt.auditory{iR,iC} = NaN;
            tot_cnt.bimodal{iR,iC}  = NaN;
            tot_cnt.tot{iR,iC}      = NaN;
        else
            tot_cnt.visual{iR,iC}   = numel(plv_hld.(fst_nme{iONE}).(scd_nme{iTWO}){1});
            tot_cnt.auditory{iR,iC} = numel(plv_hld.(fst_nme{iONE}).(scd_nme{iTWO}){2});
            tot_cnt.bimodal{iR,iC}  = numel(plv_hld.(fst_nme{iONE}).(scd_nme{iTWO}){3});
            tot_cnt.tot{iR,iC}      = plv_hld.(fst_nme{iONE}).(scd_nme{iTWO}){4};
        end
        
    end
end


vis_num = 0;
aud_num = 0;
bim_num = 0;
tot_num = 0;
for iR = 2:8
    for iC = iR+1:9
        
        vis_num = vis_num + tot_cnt.visual{iR,iC};
        aud_num = aud_num + tot_cnt.auditory{iR,iC};
        bim_num = bim_num + tot_cnt.bimodal{iR,iC};
        tot_num = tot_num + tot_cnt.tot{iR,iC};
        
    end
end

vis_num / tot_num
aud_num / tot_num
bim_num / tot_num

%% Overlap
myBinomTest(bim_num,vis_num,aud_num / tot_num(1),'lesser')
myBinomTest(bim_num,vis_num,aud_num / tot_num(1),'greater')

%% Put together numbers
% Total Number per Region - Visual
for iR = 2:9
    reg_num(iR,1) = 0;
    tot_num(iR,1) = 0;
    for iC = 2:9
        if iR ~= iC
            reg_num(iR,1) = reg_num(iR,1)+tot_cnt.visual{iR,iC};
            tot_num(iR,1) = tot_num(iR,1)+tot_cnt.tot{iR,iC};
        end
    end
end
[tot_cnt.visual(2:end,1) num2cell( reg_num(2:end) ) num2cell( reg_num(2:end) ./ tot_num(2:end) )]

% Total Number per Region - Auditory
for iR = 2:9
    reg_num(iR,1) = 0;
    tot_num(iR,1) = 0;
    for iC = 2:9
        if iR ~= iC
            reg_num(iR,1) = reg_num(iR,1)+tot_cnt.auditory{iR,iC};
            tot_num(iR,1) = tot_num(iR,1)+tot_cnt.tot{iR,iC};
        end
    end
end
[tot_cnt.auditory(2:end,1) num2cell( reg_num(2:end) ) num2cell( reg_num(2:end) ./ tot_num(2:end) )]

% Total Number per Region - Bi-Modal
for iR = 2:9
    reg_num(iR,1) = 0;
    tot_num(iR,1) = 0;
    for iC = 2:9
        if iR ~= iC
            reg_num(iR,1) = reg_num(iR,1)+tot_cnt.bimodal{iR,iC};
            tot_num(iR,1) = tot_num(iR,1)+tot_cnt.tot{iR,iC};
        end
    end
end
[tot_cnt.bimodal(2:end,1) num2cell( reg_num(2:end) ) num2cell( reg_num(2:end) ./ tot_num(2:end) )]

% Total Number per Connection - Bi-Modal
tot_cnt.pro = tot_cnt.tot;
row_cnt = 1;
for iR = 2:9
    for iC = 2:9
        if iR ~= iC
            tot_cnt.pro{iR,iC} = tot_cnt.bimodal{iR,iC} / tot_cnt.tot{iR,iC};
            hld_cnt{row_cnt,1} = [tot_cnt.pro{iR,1} '--' tot_cnt.pro{1,iC}];
            hld_cnt{row_cnt,2} = tot_cnt.bimodal{iR,iC};
            hld_cnt{row_cnt,3} = tot_cnt.bimodal{iR,iC} / tot_cnt.tot{iR,iC};
            row_cnt = row_cnt + 1;
        end
    end
end

[~,mix] = sort(cell2mat(hld_cnt(:,3)));
hld_cnt(mix,:)






















