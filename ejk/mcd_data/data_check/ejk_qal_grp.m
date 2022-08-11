function grp_out = ejk_qal_grp(cfg)

%% Loop through and fix
grp_nme = fieldnames(cfg.grp);
for iG = 1:numel(grp_nme)
    grp_nme_sub = fieldnames(cfg.grp.(grp_nme{iG}));
    for iGS = 1:numel(grp_nme_sub)
        cat_nme = fieldnames(cfg.grp.(grp_nme{iG}).(grp_nme_sub{iGS}));
        for iCA = 1:numel(cat_nme)
        
        grp_out.([grp_nme{iG}]).(grp_nme_sub{iGS}).(cat_nme{iCA}) = cfg.grp.(grp_nme{iG}).(grp_nme_sub{iGS}).(cat_nme{iCA});        
        grp_out.([grp_nme{iG}]).(grp_nme_sub{iGS}).(cat_nme{iCA})(ismember(grp_out.([grp_nme{iG}]).(grp_nme_sub{iGS}).(cat_nme{iCA}),cfg.rmv_ind)) = [];
        
        end
    end
end

end

% ismember(tot_sbj{iS},sol_sbj)