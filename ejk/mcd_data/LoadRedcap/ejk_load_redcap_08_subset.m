function dta_out = ejk_load_redcap_08_subset(dta_out,sbj_nme)

sbj_ind = ismember( dta_out.sbj_nme , sbj_nme );
dta_out_nme = fieldnames(dta_out);
for iF = 1:numel(dta_out_nme); dta_out.(dta_out_nme{iF}) = dta_out.(dta_out_nme{iF})(sbj_ind); end

end