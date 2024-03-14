tbl_out_dir = [ out_dir '/' 'Tables' '/' ]; ejk_chk_dir(tbl_out_dir);

load([ dta_dir '/' 'performance_flip_reihaneh.mat' ])

org_mdl = { 'fcn_3dm_org'                  'fcn_3dm_fll'   'fcn_3dm_fll'   'fcn_3dm_750'   'fcn_3dm_300' };

cmp_grp = { {'fcn_3dm_cbt' 'fcn_2dm_org' } {'fcn_3dm_shf'} {'fcn_2dm_fll'} {'fcn_2dm_750'} {'fcn_2dm_300'}};

%%
% Performance Table
avg_out_tbl = cell(numel(mdl_int_nme),numel(mes_int_nme));
med_out_tbl = cell(numel(mdl_int_nme),numel(mes_int_nme));
for iC = 1:numel(mes_int_nme)
    for iR = 1:numel(mdl_int_nme)
        med_val = nanmedian(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        avg_val = nanmean(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        std_val = nanstd(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        iqr_val = iqr(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        max_val = nanmax(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        min_val = nanmin(cell2mat(dta_hld.(mdl_int_nme{iR}).(mes_int_nme{iC})(:,3)));
        avg_out_tbl{iR,iC} = sprintf('%.1f (%.1f); %.1f',avg_val,std_val, max_val);
        med_out_tbl{iR,iC} = sprintf('%.1f (%.1f); %.1f',med_val,iqr_val, max_val);
    end
end
cell2csv([ tbl_out_dir '/' 'PerformanceTable_mean.csv'],   [ {''}  mes_int_nme; mdl_int_nme' avg_out_tbl ])
cell2csv([ tbl_out_dir '/' 'PerformanceTable_median.csv'], [ {''}  mes_int_nme; mdl_int_nme' med_out_tbl ])

% Diff table
avg_out_tbl = cell(numel(cmp_grp),numel(mes_int_nme));
med_out_tbl = cell(numel(cmp_grp),numel(mes_int_nme));
fdc_out_tbl = cell(numel(cmp_grp),numel(mes_int_nme));

cmp_nme     = cell(numel(cmp_grp),1);

for iC = 1:numel(mes_int_nme)
cnt = 1;
    for iD = 1:numel(org_mdl)
        for iDF = 1:numel(cmp_grp{iD})
            cmp_nme{cnt,1} = [ org_mdl{iD} ' - ' cmp_grp{iD}{iDF}];
            
            med_val = nanmedian(cell2mat(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3))) - ...
                nanmedian(cell2mat(dta_hld.(cmp_grp{iD}{iDF}).(mes_int_nme{iC})(:,3)));
            avg_val = nanmean(cell2mat(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3))) - ...
                nanmean(cell2mat(dta_hld.(cmp_grp{iD}{iDF}).(mes_int_nme{iC})(:,3)));
            std_val = nanstd(cell2mat(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3)) - ...
                cell2mat(dta_hld.(cmp_grp{iD}{iDF}).(mes_int_nme{iC})(:,3)));
            iqr_val = iqr(cell2mat(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3)) - ...
                cell2mat(dta_hld.(cmp_grp{iD}{iDF}).(mes_int_nme{iC})(:,3)));
            max_val = nanmax(cell2mat(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3)) - ...
                cell2mat(dta_hld.(cmp_grp{iD}{iDF}).(mes_int_nme{iC})(:,3)));
            min_val = nanmin(cell2mat(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3)) - ...
                cell2mat(dta_hld.(cmp_grp{iD}{iDF}).(mes_int_nme{iC})(:,3)));
            fdc_val = (sum( cell2mat(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3)) > ...
                cell2mat(dta_hld.(cmp_grp{iD}{iDF}).(mes_int_nme{iC})(:,3))) / ...
                numel(dta_hld.(org_mdl{iD}).(mes_int_nme{iC})(:,3))) * 100;
            
            avg_out_tbl{cnt,iC} = sprintf('%.1f (%.1f); %.1f',avg_val,std_val, max_val);
            med_out_tbl{cnt,iC} = sprintf('%.1f (%.1f); %.1f',med_val,iqr_val, max_val);
            fdc_out_tbl{cnt,iC} = sprintf('%.1f',fdc_val);
            
            cnt = cnt + 1;
        end
    end
end
cell2csv([ tbl_out_dir '/' 'DiffTable_mean.csv'],   [ {''}  mes_int_nme; cmp_nme avg_out_tbl ])
cell2csv([ tbl_out_dir '/' 'DiffTable_median.csv'], [ {''}  mes_int_nme; cmp_nme med_out_tbl ])
cell2csv([ tbl_out_dir '/' 'DiffTable_FDIC.csv'],   [ {''}  mes_int_nme; cmp_nme fdc_out_tbl ])

