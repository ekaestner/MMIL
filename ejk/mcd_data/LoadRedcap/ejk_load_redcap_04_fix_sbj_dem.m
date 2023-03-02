function sbj_dem = ejk_load_redcap_04_fix_sbj_dem(sbj_dem)

%% Age %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iS = 1:size(sbj_dem.sbj_nme,1)
    if ( ~isempty(sbj_dem.sbj_scn_dte{iS,1}) && ~strcmpi(sbj_dem.sbj_scn_dte{iS,1},'') ) && ...
       ( ~isempty(sbj_dem.sbj_age{iS,1})     && ~strcmpi(sbj_dem.sbj_age{iS,1},'') )
        
        scn_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_scn_dte{iS,1},'-','split'));
        if numel(scn_dte)==1 
            scn_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_scn_dte{iS,1},'/','split')); 
            scn_dte = (scn_dte(3)*365) + (scn_dte(1)*30) + scn_dte(2);    
        elseif numel(scn_dte)==3 
            scn_dte = (scn_dte(1)*365) + (scn_dte(2)*30) + scn_dte(3);    
        end   
        
        
        brt_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_age{iS,1},'-','split'));
        if numel(brt_dte)==1; 
            brt_dte = cellfun(@(x) str2num(x),regexp(sbj_dem.sbj_age{iS,1},'/','split'));
            brt_dte = (brt_dte(3)*365) + (brt_dte(1)*30) + brt_dte(2);    
        elseif numel(brt_dte)==3
            brt_dte = (brt_dte(1)*365) + (brt_dte(2)*30) + brt_dte(3);
        end
        
        sbj_dem.sbj_age{iS,1} = roundsd((scn_dte - brt_dte) / 365,3);
    else
        sbj_dem.sbj_age{iS,1} = nan;
    end
end
    
sbj_dem.sbj_age = cell2mat(sbj_dem.sbj_age);

%% Sex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbj_sex_hld = {0 'F' ; 1 'M' };
sbj_hnd_hld = {7 'R' ; 8 'L' };

for iS = 1:size(sbj_dem.sbj_nme,1)
    try sbj_dem.sbj_sex{iS,1} = sbj_sex_hld{cell2mat(sbj_sex_hld(:,1))==sbj_dem.sbj_sex{iS,1},2}; catch sbj_dem.sbj_sex{iS,1} = ''; end % 0 = Woman, 1 = Man
    try sbj_dem.sbj_hnd{iS,1} = sbj_hnd_hld{cell2mat(sbj_hnd_hld(:,1))==sbj_dem.sbj_hnd{iS,1},2}; catch sbj_dem.sbj_hnd{iS,1} = ''; end % 7 = right, 8 = left
end

end