
function [ grp_dta, grp_typ, grp_sbj] = ejk_group_create( cfg )

for iG = 1:numel( cfg.grp_inc )
    
    num_hld = [];
    grp_typ{iG} = cell(0);
    for iN = 1:numel( cfg.grp_inc{iG} )
        num_hld = [ num_hld ; cfg.grp.( cfg.grp_inc{iG}{iN} ) ];
        grp_typ{iG} = [ grp_typ{iG} ; repmat({cfg.grp_nme{iG}{iN}},numel(cfg.grp.( cfg.grp_inc{iG}{iN} )),1)];
    end
    
    grp_dta{iG} = cfg.dta( num_hld, : ); 
    grp_sbj{iG} = cfg.sbj( num_hld, : );
    
end


end