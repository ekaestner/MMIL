clc

load( [ prj_dir '/' prj_nme '/' 'groups.mat' ] )

out_dir =  '/home/ekaestner/Dropbox/McDonald Lab/Erik/Projects/McDLab/slh_atl_mem/explore/overall_cor';

%% Make Table for Cognitive
rvl_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/InitialAnalysis/Neurobio_compare/Correlation/Cognitive';

% Setup
grp_hld = dir(rvl_dir); 
    grp_hld = {grp_hld(:).name};
    grp_hld = grp_hld(3:end);
    
%{ 'TLE_Controls_pre_cog_3T' ...
%  'LTLE_post_cog_3T' }; % TLE_Controls_pre_cog_3T / LTLE_post_cog_3T

mse_hld = {    'cort_thick_ctx_aparc_annot_XREP' ...
               'fiber_FA_XREP' ...
               'subcort_vol_XREP_norm_IntracranialVolume' ...
               'wmparc_FA_wm_aparc_annot_XREP' ...
               'cort_thick_ctx_aparc_annot_XREP_LateralityIndex' ...
               'fiber_FA_XREP_LateralityIndex' ...
               'subcort_vol_XREP_LateralityIndex' ...
               'wmparc_FA_wm_aparc_annot_XREP_LateralityIndex' }; 
% cort_thick_ctx_aparc_annot_XREP                 / fiber_FA_XREP                 / subcort_vol_XREP_norm_IntracranialVolume / wmparc_FA_wm_aparc_annot_XREP
% cort_thick_ctx_aparc_annot_XREP_LateralityIndex / fiber_FA_XREP_LateralityIndex / subcort_vol_XREP_LateralityIndex         / wmparc_FA_wm_aparc_annot_XREP_LateralityIndex

for iG = 1:numel(grp_hld)
    
    if exist( [rvl_dir '/' grp_hld{iG}], 'dir' )
               
        for iM = 1:numel(mse_hld)
            
            clear tbl_dsp tbl_col
            
            grp = grp_hld{iG};
            mse = mse_hld{iM};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','238') '/' 'cross_correlation_rvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_238_rvl, mse_238_rvl_sbj, mse_238_rvl_col] = ejk_dta_frm( fcfg );
            mse_238_rvl(cellfun(@(x) strcmpi(x,'NA'),mse_238_rvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','238') '/' 'cross_correlation_pvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_238_pvl, mse_238_pvl_sbj, mse_238_pvl_col] = ejk_dta_frm( fcfg );
            mse_238_pvl(cellfun(@(x) strcmpi(x,'NA'),mse_238_pvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','dev') '/' 'cross_correlation_rvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_dev_rvl, mse_dev_rvl_sbj, mse_dev_rvl_col] = ejk_dta_frm( fcfg );
            mse_dev_rvl(cellfun(@(x) strcmpi(x,'NA'),mse_dev_rvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','dev') '/' 'cross_correlation_pvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_dev_pvl, mse_dev_pvl_sbj, mse_dev_pvl_col] = ejk_dta_frm( fcfg );
            mse_dev_pvl(cellfun(@(x) strcmpi(x,'NA'),mse_dev_pvl )) = {NaN};
            
            [ ~, ind_238, ind_dev] = intersect( mse_238_rvl_col, mse_dev_rvl_col);
            
            clear tbl_dsp tbl_col
            for iC = 1:size(mse_dev_pvl,1)
                
                tbl_dsp{iC} = cell(numel(ind_238),2);
                
                for iR = 1:numel(ind_238) %size(mse_dev_pvl,2)
                    
                    if mse_238_pvl{iC,ind_238(iR)} < .01
                        tbl_dsp{iC}{iR,1} = [ num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2)) '**'];
                    elseif mse_238_pvl{iC,ind_238(iR)} < .05
                        tbl_dsp{iC}{iR,1} = [ num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2)) '*'];
                    else
                        tbl_dsp{iC}{iR,1} = num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2));
                    end
                    
                    if mse_dev_pvl{iC,ind_dev(iR)} < .01
                        tbl_dsp{iC}{iR,2} = [ num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2)) '**'];
                    elseif mse_dev_pvl{iC,ind_dev(iR)} < .05
                        tbl_dsp{iC}{iR,2} = [ num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2)) '*'];
                    else
                        tbl_dsp{iC}{iR,2} = num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2));
                    end
                    
                end
                tbl_dsp{iC} = [ tbl_dsp{iC} repmat({' '},size(tbl_dsp{iC},1),1) ];
                tbl_col{iC} = { [mse_238_pvl_sbj{iC}, '_238'], [mse_238_pvl_sbj{iC}, '_dev'] '' };
            end
            
            tbl_col = [ {''} tbl_col{:} ];
            tbl_col = tbl_col(:,1:end-1);
            tbl_dsp = [ mse_dev_pvl_col(ind_dev)' tbl_dsp{:} ];
            tbl_dsp = tbl_dsp(:,1:end-1);
            %         cell2table( [tbl_col; tbl_dsp] );
            
            ejk_chk_dir([out_dir '/' grp '/' 'Neurobio_BY_Cognition' '/']);
            cell2csv( [out_dir '/' grp '/' 'Neurobio_BY_Cognition' '/' grp '_' mse '.csv'],[tbl_col; tbl_dsp]);
            
        end
    end
end

%% Make Table for Cognitive - SPEARMAN
rvl_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/InitialAnalysis/Neurobio_compare/Correlation/Cognitive_spearman';

% Setup
grp_hld = dir(rvl_dir); 
    grp_hld = {grp_hld(:).name};
    grp_hld = grp_hld(3:end);
    
%{ 'TLE_Controls_pre_cog_3T' ...
%  'LTLE_post_cog_3T' }; % TLE_Controls_pre_cog_3T / LTLE_post_cog_3T

mse_hld = {    'cort_thick_ctx_aparc_annot_XREP' ...
               'fiber_FA_XREP' ...
               'subcort_vol_XREP_norm_IntracranialVolume' ...
               'wmparc_FA_wm_aparc_annot_XREP' ...
               'cort_thick_ctx_aparc_annot_XREP_LateralityIndex' ...
               'fiber_FA_XREP_LateralityIndex' ...
               'subcort_vol_XREP_LateralityIndex' ...
               'wmparc_FA_wm_aparc_annot_XREP_LateralityIndex' }; 
% cort_thick_ctx_aparc_annot_XREP                 / fiber_FA_XREP                 / subcort_vol_XREP_norm_IntracranialVolume / wmparc_FA_wm_aparc_annot_XREP
% cort_thick_ctx_aparc_annot_XREP_LateralityIndex / fiber_FA_XREP_LateralityIndex / subcort_vol_XREP_LateralityIndex         / wmparc_FA_wm_aparc_annot_XREP_LateralityIndex


for iG = 1:numel(grp_hld)
    
    if exist( [rvl_dir '/' grp_hld{iG}], 'dir' )
               
        for iM = 1:numel(mse_hld)
            
            clear tbl_dsp tbl_col
            
            grp = grp_hld{iG};
            mse = mse_hld{iM};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','238') '/' 'cross_correlation_rvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_238_rvl, mse_238_rvl_sbj, mse_238_rvl_col] = ejk_dta_frm( fcfg );
            mse_238_rvl(cellfun(@(x) strcmpi(x,'NA'),mse_238_rvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','238') '/' 'cross_correlation_pvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_238_pvl, mse_238_pvl_sbj, mse_238_pvl_col] = ejk_dta_frm( fcfg );
            mse_238_pvl(cellfun(@(x) strcmpi(x,'NA'),mse_238_pvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','dev') '/' 'cross_correlation_rvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_dev_rvl, mse_dev_rvl_sbj, mse_dev_rvl_col] = ejk_dta_frm( fcfg );
            mse_dev_rvl(cellfun(@(x) strcmpi(x,'NA'),mse_dev_rvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','dev') '/' 'cross_correlation_pvalues.csv' ];
            fcfg.dta_col = 2;
            [ mse_dev_pvl, mse_dev_pvl_sbj, mse_dev_pvl_col] = ejk_dta_frm( fcfg );
            mse_dev_pvl(cellfun(@(x) strcmpi(x,'NA'),mse_dev_pvl )) = {NaN};
            
            [ ~, ind_238, ind_dev] = intersect( mse_238_rvl_col, mse_dev_rvl_col);
            
            clear tbl_dsp tbl_col
            for iC = 1:size(mse_dev_pvl,1)
                
                tbl_dsp{iC} = cell(numel(ind_238),2);
                
                for iR = 1:numel(ind_238) %size(mse_dev_pvl,2)
                    
                    if mse_238_pvl{iC,ind_238(iR)} < .01
                        tbl_dsp{iC}{iR,1} = [ num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2)) '**'];
                    elseif mse_238_pvl{iC,ind_238(iR)} < .05
                        tbl_dsp{iC}{iR,1} = [ num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2)) '*'];
                    else
                        tbl_dsp{iC}{iR,1} = num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2));
                    end
                    
                    if mse_dev_pvl{iC,ind_dev(iR)} < .01
                        tbl_dsp{iC}{iR,2} = [ num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2)) '**'];
                    elseif mse_dev_pvl{iC,ind_dev(iR)} < .05
                        tbl_dsp{iC}{iR,2} = [ num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2)) '*'];
                    else
                        tbl_dsp{iC}{iR,2} = num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2));
                    end
                    
                end
                tbl_dsp{iC} = [ tbl_dsp{iC} repmat({' '},size(tbl_dsp{iC},1),1) ];
                tbl_col{iC} = { [mse_238_pvl_sbj{iC}, '_238'], [mse_238_pvl_sbj{iC}, '_dev'] '' };
            end
            
            tbl_col = [ {''} tbl_col{:} ];
            tbl_col = tbl_col(:,1:end-1);
            tbl_dsp = [ mse_dev_pvl_col(ind_dev)' tbl_dsp{:} ];
            tbl_dsp = tbl_dsp(:,1:end-1);
            %         cell2table( [tbl_col; tbl_dsp] );
            
            ejk_chk_dir([out_dir '/' grp '/' 'Neurobio_BY_Cognition_spearman' '/']);
            cell2csv( [out_dir '/' grp '/' 'Neurobio_BY_Cognition_spearman' '/' grp '_' mse '.csv'],[tbl_col; tbl_dsp]);
            
        end
    end
end

%% Make Table for Clinical
rvl_dir = '/home/ekaestne/PROJECTS/OUTPUT/slh_atl_mem/InitialAnalysis/Neurobio_compare/Correlation/Clinical';

for iG = 1:numel(grp_hld)
    
    if exist( [rvl_dir '/' grp_hld{iG}], 'dir' )
               
        for iM = 1:numel(mse_hld)
            
            clear tbl_dsp tbl_col
            
            grp = grp_hld{iG};
            mse = mse_hld{iM};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','238') '/' 'cross_correlation_rvalues.csv' ];
            fcfg.dta_col = 2;
            fcfg.transpose = 1;
            [ mse_238_rvl, mse_238_rvl_sbj, mse_238_rvl_col] = ejk_dta_frm( fcfg );
            mse_238_rvl(cellfun(@(x) strcmpi(x,'NA'),mse_238_rvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','238') '/' 'cross_correlation_pvalues.csv' ];
            fcfg.dta_col = 2;
            fcfg.transpose = 1;
            [ mse_238_pvl, mse_238_pvl_sbj, mse_238_pvl_col] = ejk_dta_frm( fcfg );
            mse_238_pvl(cellfun(@(x) strcmpi(x,'NA'),mse_238_pvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','dev') '/' 'cross_correlation_rvalues.csv' ];
            fcfg.dta_col = 2;
            fcfg.transpose = 1;
            [ mse_dev_rvl, mse_dev_rvl_sbj, mse_dev_rvl_col] = ejk_dta_frm( fcfg );
            mse_dev_rvl(cellfun(@(x) strcmpi(x,'NA'),mse_dev_rvl )) = {NaN};
            
            fcfg = [];
            fcfg.dta_loc = [ rvl_dir '/' grp '/' strrep(mse,'XREP','dev') '/' 'cross_correlation_pvalues.csv' ];
            fcfg.dta_col = 2;
            fcfg.transpose = 1;
            [ mse_dev_pvl, mse_dev_pvl_sbj, mse_dev_pvl_col] = ejk_dta_frm( fcfg );
            mse_dev_pvl(cellfun(@(x) strcmpi(x,'NA'),mse_dev_pvl )) = {NaN};
            
            [ ~, ind_238, ind_dev] = intersect( mse_238_rvl_col, mse_dev_rvl_col);
            
            clear tbl_dsp tbl_col
            for iC = 1:size(mse_dev_pvl,1)
                
                tbl_dsp{iC} = cell(numel(ind_238),2);
                
                for iR = 1:numel(ind_238) %size(mse_dev_pvl,2)
                    
                    if mse_238_pvl{iC,ind_238(iR)} < .01
                        tbl_dsp{iC}{iR,1} = [ num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2)) '**'];
                    elseif mse_238_pvl{iC,ind_238(iR)} < .05
                        tbl_dsp{iC}{iR,1} = [ num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2)) '*'];
                    else
                        tbl_dsp{iC}{iR,1} = num2str(roundsd(mse_238_rvl{iC,ind_238(iR)},2));
                    end
                    
                    if mse_dev_pvl{iC,ind_dev(iR)} < .01
                        tbl_dsp{iC}{iR,2} = [ num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2)) '**'];
                    elseif mse_dev_pvl{iC,ind_dev(iR)} < .05
                        tbl_dsp{iC}{iR,2} = [ num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2)) '*'];
                    else
                        tbl_dsp{iC}{iR,2} = num2str(roundsd(mse_dev_rvl{iC,ind_dev(iR)},2));
                    end
                    
                end
                tbl_dsp{iC} = [ tbl_dsp{iC} repmat({' '},size(tbl_dsp{iC},1),1) ];
                tbl_col{iC} = { [mse_238_pvl_sbj{iC}, '_238'], [mse_238_pvl_sbj{iC}, '_dev'] '' };
            end
            
            tbl_col = [ {''} tbl_col{:} ];
            tbl_col = tbl_col(:,1:end-1);
            tbl_dsp = [ mse_dev_pvl_col(ind_dev)' tbl_dsp{:} ];
            tbl_dsp = tbl_dsp(:,1:end-1);
            %         cell2table( [tbl_col; tbl_dsp] );
            
            ejk_chk_dir([out_dir '/' grp '/' 'Neurobio_BY_Clinical' '/']);
            cell2csv( [out_dir '/' grp '/' 'Neurobio_BY_Clinical' '/' grp '_' mse '.csv'],[tbl_col; tbl_dsp]);
            
        end
    end
end
