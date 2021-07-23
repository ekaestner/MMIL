%%
fld_nme = { '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/DTI/wmparc_FA_wm/Raw' ...
            '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/DTI/fiber_FA/Raw' ...
            '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/MRI/subcort_vol_ICV_cor/Raw' };

grp_nme = { 'tle_controls_pre_3T_allSurg_all'   ...
            'tle_controls_pre_3T_allSurg_left'  ...
            'tle_controls_pre_3T_allSurg_right' ...
            'tle_post_3T_ATLonly_left' ...
            'tle_post_3T_ATLonly_right' };

for iF = 1:numel(fld_nme)
    for iG = 1:numel(grp_nme)
        
        rvl_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_rvalues.csv' ]);
        pvl_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_pvalues.csv' ]);
        num_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_n.csv' ]);
        
        stt_hld = rvl_hld;
        
        for iC=2:size(rvl_hld,2); for iR=2:size(rvl_hld,1)
            try num_str = num2str(num_hld{iR,iC}); catch; num_str = ''; end
            try rvl_str = num2str(roundsd(rvl_hld{iR,iC},2)); catch; num_str = ''; end
            try pvl_str = num2str(roundsd(pvl_hld{iR,iC},2)); catch; num_str = ''; end
            stt_hld{iR,iC} = [ 'r(' num_str ') = ' rvl_str(2:end)  '; p = ' pvl_str(2:end)];
        end; end
        cell2csv([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_stats.csv' ],stt_hld)
    end
end

%%
fld_nme = { '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Cognitive/Correlation' };

grp_nme = { 'TLE_pre_post'   ...
            'LTLE_pre_post'  ...
            'RTLE_pre_post' };

for iF = 1:numel(fld_nme)
    for iG = 1:numel(grp_nme)
        
        rvl_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_rvalues.csv' ]);
        pvl_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_pvalues.csv' ]);
        num_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_n.csv' ]);
        
        stt_hld = rvl_hld;
        
        for iC=2:size(rvl_hld,2); for iR=2:size(rvl_hld,1)
            try num_str = num2str(num_hld{iR,iC}); catch; num_str = ''; end
            try rvl_str = num2str(roundsd(rvl_hld{iR,iC},2)); catch; num_str = ''; end
            try pvl_str = num2str(roundsd(pvl_hld{iR,iC},2)); catch; num_str = ''; end
            stt_hld{iR,iC} = [ 'r(' num_str ') = ' rvl_str(2:end)  '; p = ' pvl_str(2:end)];
        end; end
        cell2csv([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_stats.csv' ],stt_hld)
    end
end

%%
fld_nme = { '/home/ekaestne/PROJECTS/OUTPUT/PostOperative/Naming_final_sample/SpecificCor/Clinical/Correlation' };

grp_nme = { 'TLE_Controls_pre_cln'   ...
            'LTLE_post_cln'  ...
            'RTLE_post_cln' };

for iF = 1:numel(fld_nme)
    for iG = 1:numel(grp_nme)
        
        rvl_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_rvalues.csv' ]);
        pvl_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_pvalues.csv' ]);
        num_hld = mmil_readtext([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_n.csv' ]);
        
        stt_hld = rvl_hld;
        
        for iC=2:size(rvl_hld,2); for iR=2:size(rvl_hld,1)
            try num_str = num2str(num_hld{iR,iC}); catch; num_str = ''; end
            try rvl_str = num2str(roundsd(rvl_hld{iR,iC},2)); catch; num_str = ''; end
            try pvl_str = num2str(roundsd(pvl_hld{iR,iC},2)); catch; num_str = ''; end
            stt_hld{iR,iC} = [ 'r(' num_str ') = ' rvl_str(2:end)  '; p = ' pvl_str(2:end)];
        end; end
        cell2csv([ fld_nme{iF} '/' grp_nme{iG} '/' 'cross_correlation_stats.csv' ],stt_hld)
    end
end



