function [ dta_out ] = ejk_roi_summarize(cfg)

if ~isfield(cfg,'dta');     error('missing cfg.dta'); end
if ~isfield(cfg,'dta_col'); cfg.dta_col = 2:size(cfg.dta,2); end
if ~isfield(cfg,'typ');     cfg.typ = 'average'; end

ord_dta = cell2mat(cfg.dta(:,cfg.dta_col));

if strcmpi(cfg.typ,'average')

        [ out_add, out_ind] = sort(nanmean(ord_dta,2));
    
elseif strcmpi(cfg.typ,'subtract')
    if numel(cfg.dta_col)~=2; error('you made a mistake'); end
    
        
        [ out_add, out_ind] = sort( ord_dta(:,1) - ord_dta(:,2));
    
end

dta_out = [ cfg.dta(out_ind,:) num2cell(out_add)];

end