ovr_dir = [ prj_dir '/' prj_nme '/' 'InitialAnalysis_v2' '/' 'Cognitive' '/' 'Correlations' '/'];
out_dir = [ prj_dir '/' prj_nme '/' 'InitialAnalysis_v2' '/' 'Cognitive' '/' 'Correlations' '/' 'hld' '/'];

grp_nme = { 'total_pre_cog_tle_hc' 'surgery_ltle_atl' 'surgery_ltle_slah' 'surgery_rtle_atl' 'surgery_rtle_slah' };
tst_nme = { 'lm2.pre'              'lm2.chg'          'lm2.chg'           'lm2.chg'          'lm2.chg' };

%% Build tables
% Loop through neurobio
for iF = 1:numel(inp_fle)
    for iM = 1:numel(inp_mse{iF})
        for iR = 1:numel(inp_roi{iF}{iM})
                    
            % Identify 
            if inp_roi{iF}{iM}(iR)==0
                sve_fld     = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' inp_suf{iF}{iM} '_' 'laterality' ];
            else
                sve_fld     = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} ];
                sve_fld_lat = [ inp_mse{iF}{iM} '_' mmil_spec_char(roi_nme{inp_roi{iF}{iM}(iR)},{'.'}) '_' inp_suf{iF}{iM} '_' 'laterality' ];
            end
            
            % Load bilateral stats
            ref_tbl = mmil_readtext([ ovr_dir '/' grp_nme{1} '/' sve_fld '/' 'cross_correlation_rvalues.csv' ]);
            
            bil_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+1 );
            bil_tbl_hld(:,1) = ref_tbl(:,1);
            for iG = 1:numel(grp_nme)
                bil_tbl_hld{1,iG+1} = grp_nme{iG};
                
                cor_hld = mmil_readtext([ ovr_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_rvalues.csv' ]);
                    cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),tst_nme{iG}));
                    cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                pvl_hld = mmil_readtext([ ovr_dir '/' grp_nme{iG} '/' sve_fld '/' 'cross_correlation_pvalues.csv' ]);
                    pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),tst_nme{iG}));
                    pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                
                bil_tbl_hld(2:end,iG+1) = cellfun(@(x,y) [ 'p=' num2str(roundsd(y,2)) ' ; r=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,'uni',0);
                
            end
            
            cell2csv([out_dir '/' sve_fld '.csv' ],bil_tbl_hld)
            
            % Load lateral stats
            ref_tbl = mmil_readtext([ ovr_dir '/' grp_nme{1} '/' sve_fld_lat '/' 'cross_correlation_rvalues.csv' ]);
            
            lat_tbl_hld      = cell( size(ref_tbl,1), numel(grp_nme)+1 );
            lat_tbl_hld(:,1) = ref_tbl(:,1);
            for iG = 1:numel(grp_nme)
                lat_tbl_hld{1,iG+1} = grp_nme{iG};
                
                cor_hld = mmil_readtext([ ovr_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_rvalues.csv' ]);
                    cor_hld = cor_hld(2:end,strcmpi(cor_hld(1,:),tst_nme{iG}));
                    cor_hld(cellfun(@isstr,cor_hld)) = {NaN};
                pvl_hld = mmil_readtext([ ovr_dir '/' grp_nme{iG} '/' sve_fld_lat '/' 'cross_correlation_pvalues.csv' ]);
                    pvl_hld = pvl_hld(2:end,strcmpi(pvl_hld(1,:),tst_nme{iG}));
                    pvl_hld(cellfun(@isstr,pvl_hld)) = {NaN};
                
                lat_tbl_hld(2:end,iG+1) = cellfun(@(x,y) [ 'p=' num2str(roundsd(y,2)) ' ; r=' num2str(roundsd(x,2)) ],cor_hld,pvl_hld,'uni',0);
                
            end
            
            cell2csv([out_dir '/' sve_fld_lat '.csv' ],lat_tbl_hld)
            
        end
    end
end
