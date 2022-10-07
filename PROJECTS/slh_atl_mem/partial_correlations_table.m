out_dir = [ prj_dir '/' prj_nme '/' 'FinalAnalysis' '/' 'stats' '/' 'partial_correlation__surgery_pst_cog_dti' '/' ];

load([ prj_dir '/' prj_nme '/' 'Data' '/' 'grp_img_qal.mat'])

neu_var_nme = { 'Unc' 'fusiform' 'lateralorbitofrontal' };
fld_nme = fieldnames(grp.surgery.pst_cog_dti); fld_nme = fld_nme(string_find(fld_nme,{'ltle'}));
scr_nme = { 'lm2.rci' };

%%
tbl_ref = mmil_readtext( [ out_dir '/' fld_nme{1} '/' scr_nme{1} '/' 'partial' '/' 'partial_correlation_' neu_var_nme{1} '_pvalues.csv' ] );
    num_row = size(tbl_ref,1) - 2; tot_num_row = num_row*numel(neu_var_nme); row_str_pnt = [ 1:num_row:tot_num_row tot_num_row+1];

for iS = 1:numel(scr_nme)
        
    out_tbl = cell( tot_num_row+1, numel(fld_nme)+1);
    for iF = 1:numel(fld_nme)
        out_tbl(1,iF) = fld_nme(iF);

        for iN = 1:numel(neu_var_nme)
                        
            tbl_rvl = mmil_readtext( [ out_dir '/' fld_nme{iF} '/' scr_nme{iS} '/' 'partial' '/' 'partial_correlation_' neu_var_nme{iN} '_rvalues.csv' ] );
            tbl_pvl = mmil_readtext( [ out_dir '/' fld_nme{iF} '/' scr_nme{iS} '/' 'partial' '/' 'partial_correlation_' neu_var_nme{iN} '_pvalues.csv' ] );
            
            col_ind = find(strcmpi(tbl_rvl(1,:),scr_nme{iS}));
            scr_ind = find(strcmpi(tbl_rvl(:,1),scr_nme{iS}));
            neu_ind = find(strcmpi(tbl_rvl(:,1),neu_var_nme{iN}));
            cov_ind = setxor(2:size(tbl_rvl,1),[ scr_ind neu_ind ]);
            
            out_tbl((row_str_pnt(iN):row_str_pnt(iN+1)-1)+1,1) = tbl_pvl( [ neu_ind  cov_ind] , 1);
            
            out_tbl( (row_str_pnt(iN):row_str_pnt(iN+1)-1)+1, iF+1 ) = ...
                strcat( strcat('p=',cellfun(@num2str,num2cell(cellfun(@(x) roundsd(x,2),tbl_pvl( [ neu_ind  cov_ind] , col_ind))),'uni',0)), ...
                        '; ', ...
                        strcat('rs=',cellfun(@num2str,num2cell(cellfun(@(x) roundsd(x,2),tbl_rvl( [ neu_ind  cov_ind] , col_ind))),'uni',0)));
            
        end
    end
       
    cell2csv([ prj_dir '/' prj_nme '/' 'FinalAnalysis' '/' 'tables' '/' 'Table4'  '/' 'PartialCorrelations_' scr_nme{iS} '.csv'], out_tbl)
    
end

