%%
function out_dta = ejk_clean_roi(cfg)

%%
if ~isfield(cfg,'dta');     error('No cfg.dta'); end
if ~isfield(cfg,'sbj_nme'); error('No cfg.sbj_nme'); end
if ~isfield(cfg,'sbj_col'); cfg.sbj_col = 1; end
if ~isfield(cfg,'dta_col'); cfg.dta_col = 2:size(cfg.dta,2); end
if ~isfield(cfg,'roi_nme'); cfg.roi_nme = repmat( {'all_roi'}, 1, numel(cfg.sbj_nme)); end

%%
ind_hld = 1:size(cfg.dta,2);
    ind_hld(cfg.dta_col) = [];
row_ind = 2:size(cfg.dta,1);

col_hld = cfg.dta(1,:);
dta_hld = cfg.dta(2:end, :);
    
col_wrk = col_hld(cfg.dta_col);
dta_wrk = dta_hld(:,cfg.dta_col);

col_hld = col_hld(ind_hld);
dta_hld = dta_hld(:,ind_hld);

%%
for iS = 1:numel( cfg.sbj_nme )
 
    sbj_ind = find( strcmpi( dta_hld(:,cfg.sbj_col), cfg.sbj_nme{iS}) );
    
    if numel(cfg.roi_nme{iS})~=1 && ~strcmpi(cfg.roi_nme{iS}{1},'all_roi')
        [ ~, roi_ind ] = ismember( col_wrk, cfg.roi_nme{iS});
        roi_ind = find(roi_ind);
    else
        roi_ind = 1:size(col_wrk,2);
    end
    
    dta_wrk( sbj_ind, roi_ind) = {NaN};
    
end

%%
out_dta(1,ind_hld)     = col_hld;
out_dta(1,cfg.dta_col) = col_wrk;

out_dta(row_ind,ind_hld)       = dta_hld;
out_dta(row_ind,cfg.dta_col)   = dta_wrk;

end