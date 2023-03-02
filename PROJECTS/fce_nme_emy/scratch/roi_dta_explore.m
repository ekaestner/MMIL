clear; clc;

%% Constants
prj_dir = '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/Imaging/fce_nme_emy';

out_dir = [ prj_dir '/' 'Imaging/plots/' ];

dta_dir = [ prj_dir '/' 'Imaging/ROI/' ];
csv_dir = [ prj_dir '/' 'RedcapData/' ];

dta_bta_fle = 'sub_all_block_betas_eventTent_roi_subbrick_extracted_repulled_12_15_2022.csv';
dta_lat_fle = 'grand_lat.csv';
csv_fle     = 'UW_New_Total_Database.csv';

%% Load Data
fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_bta_fle ];
[ dta_bta, dta_bta_sbj, dta_bta_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ dta_dir '/' dta_lat_fle ];
[ dta_lat, dta_lat_sbj, dta_lat_col ] = ejk_dta_frm(fcfg);

fcfg = [];
fcfg.dta_loc = [ csv_dir '/' csv_fle ];
[ csv_dta, csv_dta_sbj, csv_dta_col ] = ejk_dta_frm(fcfg);

%% Split Data
clear plt_dta

% Betas
cts_int     = unique(dta_bta(:,strcmpi(dta_bta_col,'Contrast')));
cts_int_nme = { 'Nov_plus_Rep' 'Nov_minus_Rep' 'Nov_vs_Base' 'Rep_vs_Base' };

mes_int     = { 'NZMean_1' };
mes_int_nme = { 'nzmean' };

sbj_nme = unique(dta_bta_sbj);

roi_int     = unique(dta_bta(:,strcmpi(dta_bta_col,'ROI')));
roi_int_nme = cellfun(@(x) mmil_spec_char(x,{' '},{'_'}),lower(roi_int),'uni',0);

for iC = 1:numel(cts_int)
    cts_ind = strcmpi(dta_bta(:,strcmpi(dta_bta_col,'Contrast')),cts_int{iC});
    for iM = 1:numel(mes_int)
        mes_col = strcmpi(dta_bta_col,mes_int{iM});
        for iS = 1:numel(sbj_nme)
            col_ind = 2;
            plt_dta.bta_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]){iS,1} = sbj_nme{iS};
            sbj_ind = strcmpi(dta_bta_sbj,sbj_nme{iS});
            for iRO = 1:numel(roi_int)
                roi_ind = strcmpi(dta_bta(:,strcmpi(dta_bta_col,'ROI')),roi_int{iRO});
                if sum(cts_ind & sbj_ind & roi_ind)==0
                    plt_dta.bta_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]){iS,col_ind} = NaN;
                else
                    plt_dta.bta_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]){iS,col_ind} = dta_bta{cts_ind & sbj_ind & roi_ind,mes_col};
                end
                col_ind = col_ind + 1;
            end
        end
        plt_dta.bta_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]) = [ {'sbj_ind'}  roi_int_nme' ; plt_dta.bta_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]) ];
    end
end

% Laterality
cts_int     = unique(dta_lat(:,strcmpi(dta_lat_col,'Contrast')));
cts_int_nme = { 'Nov_plus_Rep' 'Nov_minus_Rep' 'Nov_vs_Base' 'Rep_vs_Base' };

mes_int     = { 'LI..overall.' };
mes_int_nme = { 'LI' };

sbj_nme = unique(dta_lat_sbj);

roi_int     = unique(dta_lat(:,strcmpi(dta_lat_col,'ROI')));
roi_int_nme = cellfun(@(x) mmil_spec_char(x,{' '},{'_'}),lower(roi_int),'uni',0);

clear lat_dat
for iC = 1:numel(cts_int)
    cts_ind = strcmpi(dta_lat(:,strcmpi(dta_lat_col,'Contrast')),cts_int{iC});
    for iM = 1:numel(mes_int)
        mes_col = strcmpi(dta_lat_col,mes_int{iM});
        for iS = 1:numel(sbj_nme)
            col_ind = 2;
            plt_dta.lat_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]){iS,1} = sbj_nme{iS};
            sbj_ind = strcmpi(dta_lat_sbj,sbj_nme{iS});
            for iRO = 1:numel(roi_int)
                roi_ind = strcmpi(dta_lat(:,strcmpi(dta_lat_col,'ROI')),roi_int{iRO});
                if sum(cts_ind & sbj_ind & roi_ind)==0
                    plt_dta.lat_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]){iS,col_ind} = NaN;
                else
                    plt_dta.lat_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]){iS,col_ind} = dta_lat{cts_ind & sbj_ind & roi_ind,mes_col};
                end
                col_ind = col_ind + 1;
            end
        end
        plt_dta.lat_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]) = [ {'sbj_ind'}  roi_int_nme' ; plt_dta.lat_dat.([cts_int_nme{iC} '_' mes_int_nme{iM}]) ];
    end
end

% CSV Data - Beta
sbj_nme = unique(dta_bta_sbj);
csv_nme = { 'sze_sde'            'VPAII'                                  'LMII' };
csv_col = { 'side_seizure_onset' 'WMS_IV_Verbal_Paired_Associates_II_raw' 'WMS_IV_Logical_Memory_II_raw' };
bta_csv_hld = cell( size(plt_dta.bta_dat.Nov_minus_Rep_nzmean,1), numel(csv_nme) );
for iS = 1:numel(sbj_nme)
    csv_sbj_ind = strcmpi(csv_dta_sbj,[ sbj_nme{iS}(5:end) '-' 'T1' ]);
    for iC = 1:numel(csv_nme)
        col_ind = strcmpi(csv_dta_col,csv_col{iC});
        if ~(sum(csv_sbj_ind)==0)
            bta_csv_hld{iS,iC} = csv_dta{csv_sbj_ind,col_ind};
            if isnumeric(bta_csv_hld{iS,iC}) && bta_csv_hld{iS,iC}==-9999
                bta_csv_hld{iS,iC} = NaN;
            end
        else
            if ischar(csv_dta{1,col_ind})
                bta_csv_hld{iS,iC} = '';
            elseif isnumeric(csv_dta{1,col_ind})
                bta_csv_hld{iS,iC} = NaN;
            end
        end
    end
end

% CSV Data - Laterality
sbj_nme = unique(dta_lat_sbj);
csv_nme = { 'sze_sde'            'VPAII'                                  'LMII' };
csv_col = { 'side_seizure_onset' 'WMS_IV_Verbal_Paired_Associates_II_raw' 'WMS_IV_Logical_Memory_II_raw' };
lat_csv_hld = cell( size(plt_dta.lat_dat.Nov_minus_Rep_LI,1), numel(csv_nme) );
for iS = 1:numel(sbj_nme)
    csv_sbj_ind = strcmpi(csv_dta_sbj,[ sbj_nme{iS} '-' 'T1' ]);
    for iC = 1:numel(csv_nme)
        col_ind = strcmpi(csv_dta_col,csv_col{iC});
        if ~(sum(csv_sbj_ind)==0)
            lat_csv_hld{iS,iC} = csv_dta{csv_sbj_ind,col_ind};
            if isnumeric(lat_csv_hld{iS,iC}) && lat_csv_hld{iS,iC}==-9999
                lat_csv_hld{iS,iC} = NaN;
            end
        else
            if ischar(csv_dta{1,col_ind})
                lat_csv_hld{iS,iC} = '';
            elseif isnumeric(csv_dta{1,col_ind})
                lat_csv_hld{iS,iC} = NaN;
            end
        end
    end
end

%% Groups
% Betas %%%%%%%%%%%%%%%
grp.bta_dat.disease.HC  = string_find(plt_dta.bta_dat.Nov_minus_Rep_nzmean(2:end,1),'CDUW');
grp.bta_dat.disease.TLE = string_find(plt_dta.bta_dat.Nov_minus_Rep_nzmean(2:end,1),'-DUW');

grp.bta_dat.lateralized.HC   = grp.bta_dat.disease.HC;
grp.bta_dat.lateralized.LTLE = find(strcmpi(bta_csv_hld(:,strcmpi(csv_nme,'sze_sde')),'left temporal'));
grp.bta_dat.lateralized.RTLE = find(strcmpi(bta_csv_hld(:,strcmpi(csv_nme,'sze_sde')),'right temporal'));

% Laterality %%%%%%%%%%%%%%%
grp.lat_dat.disease.HC  = string_find(plt_dta.lat_dat.Nov_minus_Rep_LI(2:end,1),'CDUW');
grp.lat_dat.disease.TLE = string_find(plt_dta.lat_dat.Nov_minus_Rep_LI(2:end,1),'^DUW');

grp.lat_dat.lateralized.HC   = grp.lat_dat.disease.HC;
grp.lat_dat.lateralized.LTLE = find(strcmpi(lat_csv_hld(:,strcmpi(csv_nme,'sze_sde')),'left temporal'));
grp.lat_dat.lateralized.RTLE = find(strcmpi(lat_csv_hld(:,strcmpi(csv_nme,'sze_sde')),'right temporal'));

%% Group Plots
grp_ovr = { 'lat_dat' 'bta_dat' };

grp_typ = { 'disease' 'lateralized' };

cts_int_nme = { 'Nov_plus_Rep' 'Nov_minus_Rep' 'Nov_vs_Base' 'Rep_vs_Base' };

% colors
grp_nme = { { 'HC'             'TLE' }         { 'HC'             'LTLE'            'RTLE'} };
grp_col = { { rgb('dark grey') rgb('orange') } { rgb('dark grey') rgb('light teal') rgb('light red')} };

% Plot
for iT = 1:numel(grp_ovr)
    
    fld_nme_hld = fieldnames(plt_dta.(grp_ovr{iT}));
    mes_int_nme = strrep(fld_nme_hld(string_find(fld_nme_hld,cts_int_nme{1})),[cts_int_nme{1} '_'],'');
    
    for iGT = 1:numel(grp_typ)
        
        for iM = 1:numel(mes_int_nme)
            plt_out_dir = [ out_dir '/' 'ROIs' '/' mes_int_nme{iM} '_' grp_typ{iGT}];
            ejk_chk_dir(plt_out_dir);
            
            roi_int_nme = plt_dta.(grp_ovr{iT}).(fld_nme_hld{1})(1,2:end);
            cts_int_nme = strrep(fld_nme_hld(string_find(fld_nme_hld,mes_int_nme{iM})),[ '_' mes_int_nme{iM}],'');
            
            for iRO = 1:numel(roi_int_nme)
                
                clear fcfg
                
                fcfg = [];
                
                cts_ind = 1;
                cts_xps = 1;
                for iC = 1:numel(cts_int_nme)
                    for iG = 1:numel(grp_nme{iGT})
                        fcfg.xdt{cts_ind} = cts_xps;
                        fcfg.ydt{cts_ind} = cell2mat(plt_dta.(grp_ovr{iT}).([cts_int_nme{iC} '_' mes_int_nme{iM}])( grp.(grp_ovr{iT}).(grp_typ{iGT}).(grp_nme{iGT}{iG})+1, strcmpi(plt_dta.(grp_ovr{iT}).([cts_int_nme{iC} '_' mes_int_nme{iM}])(1,:),roi_int_nme{iRO})));
                        
                        fcfg.fce_col(cts_ind)     = { grp_col{iGT}{strcmpi(grp_nme{iGT},grp_nme{iGT}{iG})} };
                        fcfg.box_plt_col(cts_ind) = { grp_col{iGT}{strcmpi(grp_nme{iGT},grp_nme{iGT}{iG})} };
                                                
                        cts_xps = cts_xps + 1;
                        cts_ind = cts_ind + 1;
                    end
                    fcfg.xlb{cts_ind-1} = cts_int_nme{iC};
                    cts_xps = cts_xps + 1;
                end
                                
                fcfg.edg_col     = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];
                fcfg.box_plt = ones(1,numel(fcfg.xdt));
                
                fcfg.xlm = [ 0.5 max(cell2mat(fcfg.xdt))+0.5 ];
                
                fcfg.ylb = mes_int_nme(iM);
                
                fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
                fcfg.aph_val = 0.45;
                
                fcfg.hln = 0;
                fcfg.hln_col = rgb('black');
                
                fcfg.ttl = roi_int_nme{iRO};
                
                fcfg.out_dir = plt_out_dir;
                fcfg.out_nme = [ roi_int_nme{iRO} '_'  mes_int_nme{iM} ];
                
                ejk_scatter(fcfg)
                
            end
        end
    end
end

%% Scatterplots
grp_ovr = { 'lat_dat' 'bta_dat' };

grp_typ = { 'disease' 'lateralized' };

cts_int_nme = { 'Nov_plus_Rep' 'Nov_minus_Rep' 'Nov_vs_Base' 'Rep_vs_Base' };

% scores
scr_nme = { 'VPAII' 'LMII' };

scr_hld.lat_dat = lat_csv_hld;
scr_hld.bta_dat = bta_csv_hld;

% colors
grp_nme = { { 'TLE' }         { 'LTLE'            'RTLE'} };
grp_col = { { rgb('orange') } { rgb('light teal') rgb('light red')} };

% Plot
for iSC = 1:numel(scr_nme)
    for iT = 1:numel(grp_ovr)
        
        fld_nme_hld = fieldnames(plt_dta.(grp_ovr{iT}));
        mes_int_nme = strrep(fld_nme_hld(string_find(fld_nme_hld,cts_int_nme{1})),[cts_int_nme{1} '_'],'');
        
        for iGT = 1:numel(grp_typ)
            
            for iM = 1:numel(mes_int_nme)
                plt_out_dir = [ out_dir '/' 'ROIs' '/' mes_int_nme{iM} '_' grp_typ{iGT} '_scatter'];
                ejk_chk_dir(plt_out_dir);
                
                roi_int_nme = plt_dta.(grp_ovr{iT}).(fld_nme_hld{1})(1,2:end);
                cts_int_nme = strrep(fld_nme_hld(string_find(fld_nme_hld,mes_int_nme{iM})),[ '_' mes_int_nme{iM}],'');
                
                for iRO = 1:numel(roi_int_nme)
                      
                    figure('Visible','off');
                    
                    for iC = 1:numel(cts_int_nme)
                        
                        sub_plt = subplot(sqrt(numel(cts_int_nme)),sqrt(numel(cts_int_nme)),iC);
                        
                        clear fcfg
                        
                        fcfg = [];
                        
                        for iG = 1:numel(grp_nme{iGT})
                            fcfg.xdt{iG} = cell2mat(scr_hld.(grp_ovr{iT})(grp.(grp_ovr{iT}).(grp_typ{iGT}).(grp_nme{iGT}{iG}),strcmpi(csv_nme,scr_nme{iSC})));
                            fcfg.ydt{iG} = cell2mat(plt_dta.(grp_ovr{iT}).([cts_int_nme{iC} '_' mes_int_nme{iM}])( grp.(grp_ovr{iT}).(grp_typ{iGT}).(grp_nme{iGT}{iG})+1, strcmpi(plt_dta.(grp_ovr{iT}).([cts_int_nme{iC} '_' mes_int_nme{iM}])(1,:),roi_int_nme{iRO})));
                            
                            fcfg.fce_col(iG)     = { grp_col{iGT}{strcmpi(grp_nme{iGT},grp_nme{iGT}{iG})} };
                            fcfg.box_plt_col(iG) = { grp_col{iGT}{strcmpi(grp_nme{iGT},grp_nme{iGT}{iG})} };                            
                        end
                        
                        fcfg.edg_col = [ repmat({[0 0 0]},1,numel(fcfg.xdt)) ];
                        
                        fcfg.trd_lne = ones(1,numel(grp_nme{iGT}));
                        
                        fcfg.xlb = scr_nme(iSC);
                        fcfg.ylb = [ roi_int_nme{iRO} '_' mes_int_nme(iM)];
                        
                        fcfg.mkr_sze = repmat(20,1,numel(fcfg.xdt));
                        fcfg.aph_val = 0.45;
                        
                        fcfg.hln = 0;
                        fcfg.hln_col = rgb('black');
                        
                        fcfg.ttl = cts_int_nme{iC};
                        
                        fcfg.sbp = sub_plt;
                        
                        ejk_scatter(fcfg)
                        
                    end
                   
                    fcfg.out_dir = plt_out_dir;
                    fcfg.out_nme = [ roi_int_nme{iRO} '_'  mes_int_nme{iM} ];
                    
                    print(gcf,[ plt_out_dir '/' scr_nme{iSC} '_' roi_int_nme{iRO} '_'  mes_int_nme{iM} '.png'],'-dpng')
                    close all
                    
                end
            end
        end
    end
end



