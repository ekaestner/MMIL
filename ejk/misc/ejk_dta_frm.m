function [ dta_out, sbj_out, dta_col_out, msc_out ] = ejk_dta_frm( cfg )

if ~isfield( cfg, 'dta_col');   cfg.dta_col = 2; end
if ~isfield( cfg, 'transpose'); cfg.transpose = 0; end
if ~isfield( cfg, 'all_num');   cfg.all_num = 0; end

%%
dta_out     = mmil_readtext( cfg.dta_loc );
if cfg.transpose; dta_out = dta_out'; end
sbj_out     = dta_out(2:end,1);
dta_col_out = dta_out(1,cfg.dta_col:end);
dta_out     = dta_out(2:end,cfg.dta_col:end);

if cfg.all_num
    dta_out = cell2mat(dta_out);
end

if cfg.dta_col > 2
    msc_out = dta_out(2:end,2:cfg.dta_col);
end

for iC = 1:numel(dta_col_out)
   if any(cellfun(@isempty,dta_out(:,iC)))
       if ~all(cellfun(@isempty,dta_out(:,iC))) && ~strcmpi( class(dta_out{find(cellfun(@isempty,dta_out(:,iC)),1),iC}), class(dta_out{find(~cellfun(@isempty,dta_out(:,iC)),1),iC}))
           if ischar( class(dta_out{find(~cellfun(@isempty,dta_out(:,iC)),1),iC}))
               dta_out(cellfun(@isempty,dta_out(:,iC)),iC) = {''};        
           end
       end
   end
end

end