function sbj_sze = ejk_load_redcap_05_fix_sbj_sze(sbj_sze,sbj_dem);

%% Duration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iS = 1:size(sbj_sze.sbj_nme,1)
    sbj_sze.sbj_sze_dur(iS,1) = sbj_dem.sbj_age(iS,1) - sbj_sze.sbj_age_ons(iS,1);
end

%% MTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbj_sde_ons_hld = {1 'L'   ; 2 'R' };
sbj_mts_hld     = {1 'N/A' ; 2 'L' ; 3 'R'};

for iS = 1:size(sbj_sze.sbj_nme,1)
    try sbj_sze.sbj_sde_ons{iS,1} = sbj_sde_ons_hld{cell2mat(sbj_sde_ons_hld(:,1))== sbj_sze.sbj_sde_ons{iS,1},2}; catch sbj_sze.sbj_sde_ons{iS,1} = ''; end% 1 = left, 2 = right
    try sbj_sze.sbj_mts{iS,1}     = sbj_mts_hld{cell2mat(sbj_mts_hld(:,1))        == sbj_sze.sbj_mts{iS,1},2};     catch sbj_sze.sbj_mts{iS,1} = ''; end % 1 = no, 2 = left MTS, 3 = right MTS
end

end