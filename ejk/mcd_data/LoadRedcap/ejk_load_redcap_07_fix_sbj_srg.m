function sbj_srg = ejk_load_redcap_07_fix_sbj_srg(sbj_srg);

%% Surgery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
srg_typ_hld = {1 'diagnostic' ; 2 'ATL' ; 3 'ATL +' ; 4 'amygdalohippocampectomy' ; 5 'neuroablation (hippocampus & amygdala)' ; 6 'neuroablation (hippocampus & amygdala +)' ; 7 'lesionectomy' ; 8 'lesionectomy +' ; 9 'extratemporal resection' ; 10 'multi-lobar resection' ; 11 'hemispherectomy' ; 12 'neocortical resection only' ; 13 'corpus callosotomy' ; 14 'multiple subpial transections' ; 998 'other' ; 999 'unknown' };
eng_out_hld = {1 'I' ; 2 'II' ; 3 'III' ; 4 'IV' };
srg_sde_hld = {1 'L' ; 2 'R' }; 

for iS = 1:size(sbj_srg.sbj_nme,1)
    try sbj_srg.srg_typ{iS,1} = srg_typ_hld{cell2mat(srg_typ_hld(:,1)) == sbj_srg.srg_typ{iS,1},2}; catch sbj_srg.srg_typ{iS,1} = ''; end
    try sbj_srg.eng_out{iS,1} = eng_out_hld{cell2mat(eng_out_hld(:,1)) == sbj_srg.eng_out{iS,1},2}; catch sbj_srg.eng_out{iS,1} = ''; end
    try sbj_srg.srg_sde{iS,1} = srg_sde_hld{cell2mat(srg_sde_hld(:,1)) == sbj_srg.srg_sde{iS,1},2}; catch sbj_srg.srg_sde{iS,1} = ''; end
end

end