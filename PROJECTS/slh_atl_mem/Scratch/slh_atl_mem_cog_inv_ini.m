%% Load Data
out_dir = [ prj_dir '/' prj_nme '/' 'InitialAnalysis_v2' '/' 'Cognitive' '/' ];
    ejk_chk_dir(out_dir);

load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp.mat'])

fcfg = [];
fcfg.dta_loc = [ prj_dir '/' prj_nme '/' 'Data' '/' 'Cognitive.csv'];
fcfg.dta_col = 2;
[ cog_dta, cog_dta_sbj, cog_dta_col] = ejk_dta_frm( fcfg );

%% WMS report
out_dir_wms = [ out_dir '/' 'WMS' '/' ];
    ejk_chk_dir(out_dir_wms);

% Patients
epd_num = [ grp.site.ucsd ; grp.site.ucsf ];
    
% Table
[tbl, ~, ~, lbl] = crosstab( categorical(grp.grp_cat.site), categorical(grp.grp_cat.lm2_wms_ver) );
cell2csv([ out_dir_wms '/' 'lm2_cross_table.csv'],[ {'LM2'} lbl(1:size(tbl,2),2)'; lbl(1:size(tbl,1),1) num2cell(tbl) ])

[tbl, ~, ~, lbl] = crosstab( categorical(grp.grp_cat.site), categorical(grp.grp_cat.vp2_wms_ver) );
cell2csv([ out_dir_wms '/' 'vp2_cross_table.csv'],[ {'VP2'} lbl(1:size(tbl,2),2)'; lbl(1:size(tbl,1),1) num2cell(tbl) ])

% LM2 Plot
fcfg = [];

fcfg.xdt = { 1 2 4 5};
fcfg.ydt = { cell2mat(cog_dta(intersect(grp.lm2_wms_ver.III,epd_num),strcmpi(cog_dta_col,'lm2_pre'))) cell2mat(cog_dta(grp.lm2_wms_ver.IV,strcmpi(cog_dta_col,'lm2_pre'))) ...
             cell2mat(cog_dta(grp.lm2_wms_ver.III,strcmpi(cog_dta_col,'lm2_chg'))) cell2mat(cog_dta(grp.lm2_wms_ver.IV,strcmpi(cog_dta_col,'lm2_chg'))) };

fcfg.fce_col     = { rgb('orange') rgb('blue')  rgb('orange') rgb('blue') };
fcfg.edg_col     = { [0 0 0] [0 0 0] [0 0 0] [0 0 0] };
fcfg.box_plt_col = { rgb('orange') rgb('blue')  rgb('orange') rgb('blue') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Pre III' 'Pre IV' 'Change III' 'Change IV' };
fcfg.xlm = [ 0.5 fcfg.xdt{end}+0.5 ];
fcfg.ylb = {'LM2'};
fcfg.ylm = [ min(cat(1,fcfg.ydt{:})) max(cat(1,fcfg.ydt{:})) ];
fcfg.out_dir = out_dir_wms;
fcfg.out_nme = 'LM2.png';
ejk_scatter(fcfg);
            
% VP2 Plot
fcfg = [];

fcfg.xdt = { 1 2 3 5 6 7};
fcfg.ydt = { cell2mat(cog_dta(grp.vp2_wms_ver.II,strcmpi(cog_dta_col,'vp2_pre'))) cell2mat(cog_dta(intersect(grp.vp2_wms_ver.III,epd_num),strcmpi(cog_dta_col,'vp2_pre')))  cell2mat(cog_dta(grp.vp2_wms_ver.IV,strcmpi(cog_dta_col,'vp2_pre'))) ...
             cell2mat(cog_dta(grp.vp2_wms_ver.II,strcmpi(cog_dta_col,'vp2_chg'))) cell2mat(cog_dta(grp.vp2_wms_ver.III,strcmpi(cog_dta_col,'vp2_chg')))                     cell2mat(cog_dta(grp.vp2_wms_ver.IV,strcmpi(cog_dta_col,'vp2_chg'))) };

fcfg.fce_col     = { rgb('black') rgb('orange') rgb('blue')  rgb('black') rgb('orange') rgb('blue') };
fcfg.edg_col     = { [0 0 0] [0 0 0] [0 0 0] [0 0 0]  [0 0 0] [0 0 0] };
fcfg.box_plt_col = { rgb('black') rgb('orange') rgb('blue')  rgb('black') rgb('orange') rgb('blue') };

fcfg.box_plt = ones(1,numel(fcfg.xdt));
fcfg.xlb = { 'Pre II' 'Pre III' 'Pre IV' 'Change II' 'Change III' 'Change IV' };
fcfg.xlm = [ 0.5 fcfg.xdt{end}+0.5 ];
fcfg.ylb = {'VP2'};
fcfg.ylm = [ min(cat(1,fcfg.ydt{:})) max(cat(1,fcfg.ydt{:})) ];
fcfg.out_dir = out_dir_wms;
fcfg.out_nme = 'VP2.png';
ejk_scatter(fcfg);

%% Categorize Patients - RAW CHANGE
cog_cat_nme = { 'lm2_chg' 'vp2_chg' 'cv2_chg' };

cog_chg_cat = { 'chg_thr' 'chg_for' 'chg_fve' };
cog_chg_num = [ 3         4         5 ];

% Categorize based on change
for iCO = 1:numel(cog_cat_nme)
    for iNC = 1:numel(cog_chg_cat)
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC}) = repmat({''},size(cog_dta,1),1);
        dta_use = cog_dta(:,strcmpi(cog_dta_col,cog_cat_nme{iCO}));
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC})(cell2mat(dta_use <= -cog_chg_num(iNC))) = {'Impaired'};
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC})(cell2mat(dta_use > -cog_chg_num(iNC)) & cell2mat(dta_use < cog_chg_num(iNC))) = {'NoChange'};
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC})(cell2mat(dta_use >= cog_chg_num(iNC))) = {'Improved'};
    end
end

% Put into table
tbl_col_cnt = 1;
for iCO = 1:numel(cog_cat_nme)
    for iNC = 1:numel(cog_chg_cat)
        tbl_dta(:,tbl_col_cnt) = cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC});
        tbl_nme(1,tbl_col_cnt) = { [cog_cat_nme{iCO} '_' cog_chg_cat{iNC}] };
        tbl_dsg(tbl_col_cnt,:) = { 'count' tbl_nme{1,tbl_col_cnt} 'Impaired/NoChange/Improved' 1 };
        tbl_col_cnt = tbl_col_cnt+1;
    end
end

%% Save tables out - RAW CHANGE
ejk_chk_dir([ out_dir '/' 'tables' '/' ]);

grp_nme = fieldnames(grp); grp_nme(end) = [];

for iG = 1:numel(grp_nme)

    grp_lvl = fieldnames(grp.(grp_nme{iG}));
    
    for iR = 1:size(tbl_dsg,1)
        for iC = 1:numel(grp_lvl)
           tbl_dsg_use{iR,iC} = [ tbl_dsg{iR,1} ',' num2str(tbl_dsg{iR,4}) ',' grp_lvl{iC} ',' tbl_dsg{iR,2} ',' tbl_dsg{iR,3}]; 
        end
    end
    
    fcfg = [];
    fcfg.tbl = tbl_dsg_use;
    fcfg.dta = {[ tbl_nme ; tbl_dta ]};
    fcfg.grp = grp.(grp_nme{iG});
    tbl_out = ejk_create_table( fcfg );

    cell2csv( [ out_dir '/' 'tables' '/' grp_nme{iG} '_postop_category_rawchange.csv'], [ [ {''} grp_lvl' ] ; tbl_nme' tbl_out ] )
    
    clear tbl_dsg_use
    
end

%% Categorize Patients - RCI CHANGE
cog_cat_nme = { 'lm2_rci' 'vp2_rci' 'cv2_rci' };

cog_chg_cat = { 'rci_one' 'rci_one_half' 'rci_two' };
cog_chg_num = [ 1         1.5        2 ];

% Categorize based on change
for iCO = 1:numel(cog_cat_nme)
    for iNC = 1:numel(cog_chg_cat)
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC}) = repmat({''},size(cog_dta,1),1);
        dta_use = cog_dta(:,strcmpi(cog_dta_col,cog_cat_nme{iCO}));
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC})(cell2mat(dta_use <= -cog_chg_num(iNC))) = {'Impaired'};
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC})(cell2mat(dta_use > -cog_chg_num(iNC)) & cell2mat(dta_use < cog_chg_num(iNC))) = {'NoChange'};
        cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC})(cell2mat(dta_use >= cog_chg_num(iNC))) = {'Improved'};
    end
end

% Put into table
tbl_col_cnt = 1;
for iCO = 1:numel(cog_cat_nme)
    for iNC = 1:numel(cog_chg_cat)
        tbl_dta(:,tbl_col_cnt) = cog_cat.(cog_cat_nme{iCO}).(cog_chg_cat{iNC});
        tbl_nme(1,tbl_col_cnt) = { [cog_cat_nme{iCO} '_' cog_chg_cat{iNC}] };
        tbl_dsg(tbl_col_cnt,:) = { 'count' tbl_nme{1,tbl_col_cnt} 'Impaired/NoChange/Improved' 1 };
        tbl_col_cnt = tbl_col_cnt+1;
    end
end

%% Save tables out - RAW CHANGE
ejk_chk_dir([ out_dir '/' 'tables' '/' ]);

grp_nme = fieldnames(grp); grp_nme(end) = [];

for iG = 1:numel(grp_nme)

    grp_lvl = fieldnames(grp.(grp_nme{iG}));
    
    for iR = 1:size(tbl_dsg,1)
        for iC = 1:numel(grp_lvl)
           tbl_dsg_use{iR,iC} = [ tbl_dsg{iR,1} ',' num2str(tbl_dsg{iR,4}) ',' grp_lvl{iC} ',' tbl_dsg{iR,2} ',' tbl_dsg{iR,3}]; 
        end
    end
    
    fcfg = [];
    fcfg.tbl = tbl_dsg_use;
    fcfg.dta = {[ tbl_nme ; tbl_dta ]};
    fcfg.grp = grp.(grp_nme{iG});
    tbl_out = ejk_create_table( fcfg );

    cell2csv( [ out_dir '/' 'tables' '/' grp_nme{iG} '_postop_category_rci.csv'], [ [ {''} grp_lvl' ] ; tbl_nme' tbl_out ] )
    
    clear tbl_dsg_use
    
end

%% Cognitive group level effects: Significant declines (chg, pct, rci)
out_dir_stt = [ out_dir '/' 'stats' '/' ]; ejk_chk_dir(out_dir_stt);

grp_nme = fieldnames(grp); grp_nme(end) = [];

cog_tst_nme = { 'lm2' 'vp2' 'cv2' };
scr_typ_nme = { 'chg' 'rci' 'pct' };

scr_col_ind = 1;
for iCO = 1:numel(cog_tst_nme)
    for iST = 1:numel(scr_typ_nme)
        scr_col{scr_col_ind} = [ cog_tst_nme{iCO} '_' scr_typ_nme{iST} ];
        scr_col_ind = scr_col_ind+1;
    end
end

for iG = 1:numel(grp_nme)
    
    grp_nme_stt = fieldnames(grp.(grp_nme{iG}));
    for iGN = 1:numel(grp_nme_stt)
        
        [~, use_dta_col ] = intersect( cog_dta_col, scr_col );
        
        fcfg = [];
        fcfg.grp     = grp.(grp_nme{iG});
        fcfg.grp_inc = {grp_nme_stt(iGN)};
        fcfg.grp_nme = {grp_nme_stt(iGN)};
        fcfg.dta = cog_dta(:,use_dta_col);
        fcfg.sbj = cog_dta_sbj;
        [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );

        use_dta_col(sum(cell2mat(~isnan(grp_dta{1})))<5) = [];
        grp_dta{1}(:,sum(cell2mat(~isnan(grp_dta{1})))<5) = [];        
        
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = cog_dta_col(:,use_dta_col);
        fcfg.men     = 0;
        fcfg.out_dir = [ out_dir_stt '/' grp_nme{iG} '/' 'versus0' '/' grp_nme_stt{iGN} ];
        ejk_ttest1( fcfg );
        
    end
end

%% Cognitive group level effects: Significant differences (pre, pst, chg, pct, rci)
grp_nme = fieldnames(grp); grp_nme(end) = [];

cog_tst_nme = { 'lm2' 'vp2' 'cv2' };
scr_typ_nme = { 'pre' 'chg' 'rci' 'pct' };

scr_col_ind = 1;
for iCO = 1:numel(cog_tst_nme)
    for iST = 1:numel(scr_typ_nme)
        scr_col{scr_col_ind} = [ cog_tst_nme{iCO} '_' scr_typ_nme{iST} ];
        scr_col_ind = scr_col_ind+1;
    end
end

for iG = 1:numel(grp_nme)
    
    [~, use_dta_col ] = intersect( cog_dta_col, scr_col );
    
    % get group data
    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    fcfg = [];
    fcfg.grp     = grp.(grp_nme{iG});
    fcfg.grp_inc = {fld_nme};
    fcfg.grp_nme = {fld_nme};
    fcfg.dta = cog_dta(:,use_dta_col);
    fcfg.sbj = cog_dta_sbj;
    [ grp_dta, grp_typ, grp_sbj ] = ejk_group_create( fcfg );
    
    if numel(fld_nme)==2
        rmv_ind = [];
        for iC = 1:size(grp_dta{1},2)
            cog_tst_hld = zeros(1,numel(fld_nme));
            for iT = 1:numel(fld_nme)
                cog_tst_hld(iT) = sum(~isnan(cell2mat(grp_dta{1}(strcmpi(grp_typ{1},fld_nme{iT}),iC))));
            end
            if sum(cog_tst_hld>5)<2
                rmv_ind = [ rmv_ind iC ];
            end
        end
        use_dta_col(rmv_ind) = [];
        grp_dta{1}(:,rmv_ind) = [];
    end
    
    if numel(fld_nme)>2
        % ANOVA
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = cog_dta_col(:,use_dta_col);
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = grp_nme(iG);
        fcfg.out_dir = [ out_dir_stt '/' grp_nme{iG} '/' 'groupcomparison' '/' ];
        ejk_1way_anova( fcfg )
    elseif numel(fld_nme)==2
        % t-test
        fcfg = [];
        fcfg.sbj_nme = grp_sbj{1};
        fcfg.dta     = grp_dta{1};
        fcfg.dta_nme = cog_dta_col(:,use_dta_col);
        fcfg.grp     = grp_typ{1};
        fcfg.grp_nme = grp_nme(iG);
        fcfg.out_dir = [ out_dir_stt '/' grp_nme{iG} '/' 'groupcomparison' '/' ];
        ejk_ttest2_independent( fcfg );
        
    end
    
end

%% Cognitive table
out_dir_tbl = [ out_dir '/' 'stats' '/' ];

grp_nme = fieldnames(grp); grp_nme(end) = [];

cog_tst_nme = { 'lm2' 'vp2' 'cv2' };
scr_typ_nme = { 'pre' 'chg' 'rci' 'pct' };

scr_col_ind = 1;
for iCO = 1:numel(cog_tst_nme)
    for iST = 1:numel(scr_typ_nme)
        scr_col{scr_col_ind} = [ cog_tst_nme{iCO} '_' scr_typ_nme{iST} ];
        scr_col_ind = scr_col_ind+1;
    end
end

vrs_zro = scr_col(cellfun(@isempty,strfind(scr_col,'pre')));

clear tbl_dta
for iG = 1:numel(grp_nme)

    fld_nme = fieldnames(grp.(grp_nme{iG}));
    
    tbl_dta{1} = [ cog_dta_col ; cog_dta ];
    tbl_dta{2} = mmil_readtext([ out_dir '/' 'stats' '/' grp_nme{iG} '/' 'groupcomparison' '/' mmil_spec_char(grp_nme{iG},{'_'},{'.'}) '/' 'output_table.csv' ]);
    for iFN = 1:numel(fld_nme)
        if exist([ out_dir '/' 'stats' '/' grp_nme{iG} '/' 'versus0' '/' fld_nme{iFN} '/' 'output_table.csv' ])
            tbl_dta{iFN+2} = mmil_readtext([ out_dir '/' 'stats' '/' grp_nme{iG} '/' 'versus0' '/' fld_nme{iFN} '/' 'output_table.csv' ]);
        else
            tbl_dta{iFN+2} = [];
        end
    end
    
    % Setup table
    for iR = 1:numel(scr_col)
        for iC = 1:numel(fld_nme)
            tbl_dsg{iR,iC} = [ 'mean/std' ','  '1' ',' fld_nme{iC} ',' scr_col{iR} ];
        end
        if any(strcmpi(tbl_dta{2}(:,1),mmil_spec_char(scr_col{iR},{'_'},{'.'})))
            tbl_dsg{iR,iC+1} = [ 'copy' ','  '2' ',' scr_col{iR} ',' 'report'  ];
        else
            tbl_dsg{iR,iC+1} = 'empty';
        end
    end
    
    fcfg = [];
    fcfg.tbl = tbl_dsg;
    fcfg.dta = tbl_dta;
    fcfg.grp = grp.(grp_nme{iG});
    tbl_out = ejk_create_table( fcfg );
    
    % Note versus0 comparisons
    for iR = 1:numel(vrs_zro)
        for iC = 1:numel(fld_nme)            
            if ~isempty(tbl_dta{iC+2}) && sum(strcmpi(tbl_dta{iC+2}(:,1),mmil_spec_char(vrs_zro{iR},{'_'},{'.'})))>0
                if tbl_dta{iC+2}{ strcmpi(tbl_dta{iC+2}(:,1),mmil_spec_char(vrs_zro{iR},{'_'},{'.'})), strcmpi(tbl_dta{iC+2}(1,:),'pvalue') }<.05
                    tbl_out{strcmpi(scr_col,vrs_zro{iR}),iC} = [tbl_out{strcmpi(scr_col,vrs_zro{iR}),iC} '*'];
                end
            end              
        end
    end
    
    % Save out
    cell2csv( [ out_dir_tbl '/' grp_nme{iG} '/' grp_nme{iG} '_stat_table.csv'], [ {''} fld_nme' {'test'}; scr_col'  tbl_out ] )
    
    clear tbl_dta tbl_dsg
    
end

%% Scatter plots: Scores -BY- Group
out_dir_plt = [ out_dir '/' 'plots' '/' ];

grp_nme = fieldnames(grp); grp_nme(end) = [];

cog_tst_nme = { 'lm2' 'vp2' 'cv2' };
scr_typ_nme = { 'pre' 'pst' 'chg' 'rci' 'pct' };

dst_clr = distinguishable_colors(25);
dst_clr_cnt = 1;
for iG = 1:numel(grp_nme)
    grp_nme_plt = fieldnames(grp.(grp_nme{iG}));
    for iGP = 1:numel(grp_nme_plt)
        col_hld{iG}{iGP} = dst_clr(dst_clr_cnt,:);
        dst_clr_cnt = dst_clr_cnt+1;
    end
end    
    
for iST = 1:numel(scr_typ_nme)
    for iG = 1:numel(grp_nme)
        
        out_dir_plt_hld = [ out_dir_plt '/'  grp_nme{iG} ]; ejk_chk_dir(out_dir_plt_hld);
        plt_hld = figure('Visible','off');
                
        grp_nme_plt = fieldnames(grp.(grp_nme{iG}));
        dta_hld = [];
        for iD = 1:numel(grp_nme_plt)
            for iCT = 1:numel(cog_tst_nme)
                dta_hld = [ dta_hld ; cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' scr_typ_nme{iST}])))];
            end
        end
        dta_hld(isinf(dta_hld)) = [];
        
        for iCT = 1:numel(cog_tst_nme)
            sbp_plt = subplot(1,numel(cog_tst_nme),iCT);
            fcfg = [];
            for iD = 1:numel(grp_nme_plt)
                fcfg.xdt{iD} = iD;
                fcfg.ydt{iD} = cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' scr_typ_nme{iST}])));
            
                fcfg.fce_col{iD}     = col_hld{iG}{iD};
                fcfg.edg_col{iD}     = [0 0 0];
                fcfg.box_plt_col{iD} = col_hld{iG}{iD};
                
            end
            
            fcfg.box_plt = ones(1,numel(grp_nme_plt));
            fcfg.xlb = grp_nme_plt;
            fcfg.xlm = [ 0.5 numel(grp_nme_plt)+0.5 ];
            fcfg.ylb = {[ cog_tst_nme{iCT} '_' scr_typ_nme{iST}]};
            fcfg.ylm = [ min(dta_hld) max(dta_hld) ];
            fcfg.sbp = sbp_plt;
            ejk_scatter(fcfg);
        end
        print(gcf,[ out_dir_plt_hld '/' 'scores' '_' scr_typ_nme{iST} '.png'],'-dpng')
        close all
    end
end

%% Scatter plots: Pre-op (X) -BY- Change (Y)
for iG = 1:numel(grp_nme)
    out_dir_plt_hld = [ out_dir_plt '/'  grp_nme{iG} ]; ejk_chk_dir(out_dir_plt_hld);
    grp_nme_plt = fieldnames(grp.(grp_nme{iG}));
    for iCT = 1:numel(cog_tst_nme)
        fcfg = [];
        for iD = 1:numel(grp_nme_plt)
            fcfg.xdt{iD} = cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' 'chg'])));
            fcfg.ydt{iD} = cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' 'pre'])));
            
            fcfg.fce_col{iD}     = col_hld{iG}{iD};
            fcfg.edg_col{iD}     = [0 0 0];
            fcfg.box_plt_col{iD} = col_hld{iG}{iD};            
        end
        
        fcfg.xlb = {'Change Score'};
        fcfg.ylb = {'Pre-operative Score'};
        fcfg.ttl = cog_tst_nme{iCT};
        fcfg.out_dir = out_dir_plt_hld;
        fcfg.out_nme = [ 'scatters' '_' 'pre' '_' 'change' '_' cog_tst_nme{iCT} '.png'];
        ejk_scatter(fcfg);
        
    end
end

%% Scatter plots: Pct (X) -BY- Change (Y)
for iG = 1:numel(grp_nme)
    out_dir_plt_hld = [ out_dir_plt '/'  grp_nme{iG} ]; ejk_chk_dir(out_dir_plt_hld);
    grp_nme_plt = fieldnames(grp.(grp_nme{iG}));
    for iCT = 1:numel(cog_tst_nme)
        fcfg = [];
        for iD = 1:numel(grp_nme_plt)
            fcfg.xdt{iD} = cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' 'chg'])));
            fcfg.ydt{iD} = cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' 'pct'])));
            
            fcfg.fce_col{iD}     = col_hld{iG}{iD};
            fcfg.edg_col{iD}     = [0 0 0];
            fcfg.box_plt_col{iD} = col_hld{iG}{iD};            
        end
        
        fcfg.xlb = {'Change Score'};
        fcfg.ylb = {'Percent Change (%)'};
        fcfg.ttl = cog_tst_nme{iCT};
        fcfg.out_dir = out_dir_plt_hld;
        fcfg.out_nme = [ 'scatters' '_' 'pct' '_' 'change' '_' cog_tst_nme{iCT} '.png'];
        ejk_scatter(fcfg);
        
    end
end

%% Scatter plots: RCI (X) -BY- Change (Y)
for iG = 1:numel(grp_nme)
    out_dir_plt_hld = [ out_dir_plt '/'  grp_nme{iG} ]; ejk_chk_dir(out_dir_plt_hld);
    grp_nme_plt = fieldnames(grp.(grp_nme{iG}));
    for iCT = 1:numel(cog_tst_nme)
        fcfg = [];
        for iD = 1:numel(grp_nme_plt)
            fcfg.xdt{iD} = cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' 'chg'])));
            fcfg.ydt{iD} = cell2mat(cog_dta(grp.(grp_nme{iG}).(grp_nme_plt{iD}),strcmpi(cog_dta_col,[ cog_tst_nme{iCT} '_' 'rci'])));
            
            fcfg.fce_col{iD}     = col_hld{iG}{iD};
            fcfg.edg_col{iD}     = [0 0 0];
            fcfg.box_plt_col{iD} = col_hld{iG}{iD};            
        end
        
        fcfg.xlb = {'Change Score'};
        fcfg.ylb = {'RCI'};
        fcfg.ttl = cog_tst_nme{iCT};
        fcfg.out_dir = out_dir_plt_hld;
        fcfg.out_nme = [ 'scatters' '_' 'rci' '_' 'change' '_' cog_tst_nme{iCT} '.png'];
        ejk_scatter(fcfg);
        
    end
end

%% examine for outliers