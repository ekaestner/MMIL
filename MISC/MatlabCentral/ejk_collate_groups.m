function grp_str = ejk_collate_groups(cfg)

grp_str.sbj_nme = cfg.grp_fle(2:end,1);
for iC = 2:size(cfg.grp_fle,2)
    
    grp_str.table.(cfg.grp_fle{1,iC}) = [ unique(cfg.grp_fle(2:end,iC)) num2cell(1:numel(unique(cfg.grp_fle(2:end,iC))))' ];
    grp_str.(cfg.grp_fle{1,iC}) = zeros(size(cfg.grp_fle,1)-1,1);
        
    for iCR = 1:size(grp_str.table.(cfg.grp_fle{1,iC}),1)
        grp_str.(cfg.grp_fle{1,iC})(strcmpi( cfg.grp_fle(2:end,iC) , grp_str.table.(cfg.grp_fle{1,iC}){iCR,1})) = grp_str.table.(cfg.grp_fle{1,iC}){iCR,2};
    end
    
end

end