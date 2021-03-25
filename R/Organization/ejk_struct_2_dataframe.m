function dta_out = ejk_struct_2_dataframe( cfg, dta )

dta_out = nan( numel(cfg.sbj_nme), numel(cfg.dta_fld));

[ ~, sbj_ind] = intersect( dta.sbj_nme, cfg.sbj_nme );

%%
for iC = 1:numel(cfg.dta_fld)
    dta_out(:,iC) = dta.(cfg.dta_fld{iC})(sbj_ind);
end

end