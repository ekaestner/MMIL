function sbj_cog = ejk_load_redcap_06_fix_sbj_cog(sbj_cog,sbj_dem,sbj_srg)

%% Test Intervals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iS = 1:size(sbj_cog.sbj_nme,1)
    % Surgery -to- PostNeuroPsych
    if ~isempty(sbj_srg.srg_dte{iS,1}) && ~isempty(sbj_cog.neu_psy_tst_dte_pst{iS,1})
        srg_dte = cellfun(@(x) str2num(x),regexp(sbj_srg.srg_dte{iS,1},'-','split'));
        if numel(srg_dte)==1
            srg_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_srg_dte{iS,1},'/','split'));
            srg_dte = (srg_dte(3)*365) + (srg_dte(1)*30) + srg_dte(2);
        elseif numel(srg_dte)==3
            srg_dte = (srg_dte(1)*365) + (srg_dte(2)*30) + srg_dte(3);
        end
        
        cog_dte = cellfun(@(x) str2num(x),regexp(sbj_cog.neu_psy_tst_dte_pst{iS,1},'-','split'));
        if numel(cog_dte)==1;
            cog_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_age{iS,1},'/','split'));
            cog_dte = (cog_dte(3)*365) + (cog_dte(1)*30) + cog_dte(2);
        elseif numel(cog_dte)==3
            cog_dte = (cog_dte(1)*365) + (cog_dte(2)*30) + cog_dte(3);
        end
        
        sbj_cog.neu_psy_tst_dte_pst_gap(iS,1) = roundsd((srg_dte - cog_dte) / 365,3);
    end
    
    % MRI -to- PreNeuroPsych
    if ~isempty(sbj_dem.sbj_scn_dte{iS,1}) && ~isempty(sbj_cog.neu_psy_tst_dte{iS,1})
        scn_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_scn_dte{iS,1},'-','split'));
        if numel(scn_dte)==1
            scn_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_scn_dte{iS,1},'/','split'));
            scn_dte = (scn_dte(3)*365) + (scn_dte(1)*30) + scn_dte(2);
        elseif numel(scn_dte)==3
            scn_dte = (scn_dte(1)*365) + (scn_dte(2)*30) + scn_dte(3);
        end
        
        cog_dte = cellfun(@(x) str2num(x),regexp(sbj_cog.neu_psy_tst_dte{iS,1},'-','split'));
        if numel(cog_dte)==1;
            cog_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_age{iS,1},'/','split'));
            cog_dte = (cog_dte(3)*365) + (cog_dte(1)*30) + cog_dte(2);
        elseif numel(cog_dte)==3
            cog_dte = (cog_dte(1)*365) + (cog_dte(2)*30) + cog_dte(3);
        end
        
        sbj_cog.neu_psy_tst_dte_gap(iS,1) = roundsd((scn_dte - cog_dte) / 365,3);
    end
    
end

end