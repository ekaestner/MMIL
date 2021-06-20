function [ dta_out, sbj_out, dta_col_out, msc_out ] = ejk_dta_frm( cfg )

if ~isfield( cfg, 'dta_col'); cfg.dta_col = 2; end

%%
dta_out     = mmil_readtext( cfg.dta_loc );
sbj_out     = dta_out(2:end,1);
dta_col_out = dta_out(1,cfg.dta_col:end);
dta_out     = dta_out(2:end,cfg.dta_col:end);

if cfg.dta_col > 2
    msc_out = dta_out(2:end,2:cfg.dta_col);
end

end