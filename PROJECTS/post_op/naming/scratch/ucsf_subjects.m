load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

dta_dir = [ prj_dir '/' prj_nme '/' 'Data' ];
cog_dta_nme = [ dta_dir '/' 'Cognitive'                '_' 'QC' '.csv'];

cog_dta = mmil_readtext(cog_dta_nme);
cog_dta_col = ejk_fix_column_names(cog_dta(1,2:end));
cog_dta_sbj = cog_dta(2:end,1);
cog_dta     = cog_dta(2:end,2:end);

sbj_out = cell(size(cog_dta,1),5);
for iS = 1:size(cog_dta,1)
    
    
    sbj_out{iS,1} = cog_dta_sbj{iS,1};
    sbj_out{iS,2} = ~isnan(cog_dta{iS,2});
    sbj_out{iS,3} = ~isnan(cog_dta{iS,3});
    sbj_out{iS,4} = ~isnan(cog_dta{iS,5});
    sbj_out{iS,5} = ~isnan(cog_dta{iS,6});
    
end

inc_ind = find(sum(cell2mat(sbj_out(:,2:5)),2));
usf_ind = string_find(sbj_out(:,1),'ucsf');

sbj_out = sbj_out(intersect(inc_ind,usf_ind),:);

col_nme = { 'SubjID' 'PRE BNT' 'PRE ANT' 'POST BNT' 'POST ANT' };

cell2csv( '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming/Data/UCSF_subjects.csv', [ col_nme ; sbj_out ] )