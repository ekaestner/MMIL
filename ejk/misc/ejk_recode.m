function dta_out = ejk_recode(cfg)

if ~isfield(cfg,'swt_val'); cfg.swt_val=1; end

cfg.dta = ejk_fill_values(cfg.dta);
dta_out = cfg.dta;

%% Create recode data
for iR = 1:numel(cfg.rcd_nme)
    col_tgh{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))} = cfg.rcd_val{iR};
end

%% Recode
for iC = 1:numel(col_tgh)
    if ~isempty(col_tgh{iC})
        for iT = 1:size(col_tgh{iC},1)
            
            % Recode existing variables
            if isnumeric(col_tgh{iC}{iT,1})
                dta_out(cell2mat(cfg.dta(:,iC))==col_tgh{iC}{iT,1},iC) = {col_tgh{iC}{iT,2}};
            else
                dta_out(strcmpi(cfg.dta(:,iC),col_tgh{iC}{iT,1}),iC) = {col_tgh{iC}{iT,2}};
            end
            
            % Switch over other variables
            if cfg.swt_val
                if isnumeric(col_tgh{iC}{iT,1}) && ~isnumeric(col_tgh{iC}{iT,2})
                    dta_out(cellfun(@isnumeric,dta_out(:,iC)),iC) = {''};
                elseif ~isnumeric(col_tgh{iC}{iT,1}) && isnumeric(col_tgh{iC}{iT,2})
                    dta_out(cellfun(@ischar,dta_out(:,iC)),iC) = {NaN};
                end
            end
            
        end
    end
end

end