function [ dta_out, dta_out_col] = ejk_recode(cfg)

if ~isfield(cfg,'swt_val'); cfg.swt_val=1; end
if ~isfield(cfg,'frc_fll'); cfg.frc_fll={}; end
if ~isfield(cfg,'rcd_nme_new'); cfg.rcd_nme_new=repmat({''},1,numel(cfg.rcd_nme)); end

cfg.dta = ejk_fill_values(cfg.dta,ismember(cfg.dta_col,cfg.frc_fll));
dta_out = cfg.dta;
dta_out_col = cfg.dta_col;

%% Create recode data
col_tgh     = cell(1,size(dta_out,2));
col_tgh_new = cell(1,size(dta_out,2));
for iR = 1:numel(cfg.rcd_nme)
    if numel(col_tgh{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))}) == 0
        col_tgh{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))}{1}     = cfg.rcd_val{iR};
        if ~isempty(cfg.rcd_nme_new{iR})
            col_tgh_new{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))}{1} = [ cfg.rcd_nme{iR} '_' cfg.rcd_nme_new{iR} ];
        else
            col_tgh_new{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))}{1} = {};
        end
    else
        if ~isempty(cfg.rcd_nme_new{iR})
            col_tgh_new{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))}{numel(col_tgh{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))})+1} = [ cfg.rcd_nme{iR} '_' cfg.rcd_nme_new{iR} ];
        else
            col_tgh_new{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))}{numel(col_tgh{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))})+1} = {};
        end
        col_tgh{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))}{numel(col_tgh{find(strcmpi(cfg.dta_col,cfg.rcd_nme{iR}))})+1}     = cfg.rcd_val{iR};
    end
end

%% Recode
for iC = 1:numel(col_tgh)
    for iCC = 1:numel(col_tgh{iC})
        
        if ~isempty(col_tgh_new{iC}{iCC}); dta_out_hld = dta_out(:,iC); end
        new_col = size(dta_out_col,2)+1;
        
        for iT = 1:size(col_tgh{iC}{iCC},1)
            
            % Recode existing variables
            if isnumeric(col_tgh{iC}{iCC}{iT,1})
                num_hld = cellfun(@(x) x==col_tgh{iC}{iCC}{iT,1},dta_out(:,iC),'uni',0);
                num_hld(cellfun(@isempty,num_hld)) = {logical(0)};
                num_hld(cellfun(@numel,num_hld)>1) = {logical(0)};
                
                if isempty(col_tgh_new{iC}{iCC})
                    dta_out(cell2mat(num_hld),iC) = {col_tgh{iC}{iCC}{iT,2}};
                else
                    dta_out_hld(cell2mat(num_hld),1) = {col_tgh{iC}{iCC}{iT,2}};
                end
                
            else
                if isempty(col_tgh_new{iC}{iCC})
                    dta_out(strcmpi(dta_out(:,iC),col_tgh{iC}{iCC}{iT,1}),iC) = {col_tgh{iC}{iCC}{iT,2}};
                else
                    dta_out_hld(strcmpi(dta_out(:,iC),col_tgh{iC}{iCC}{iT,1}),1) = {col_tgh{iC}{iCC}{iT,2}};                    
                end
            end
                                   
        end
        
        % Switch over other variables
        if cfg.swt_val
            if isnumeric(col_tgh{iC}{iCC}{iT,1}) && ~isnumeric(col_tgh{iC}{iCC}{iT,2})
                dta_out(cellfun(@isnumeric,dta_out(:,iC)),iC) = {''};
            elseif ~isnumeric(col_tgh{iC}{iCC}{iT,1}) && isnumeric(col_tgh{iC}{iCC}{iT,2})
                dta_out(cellfun(@ischar,dta_out(:,iC)),iC) = {NaN};
            end
        end
        
        if ~isempty(col_tgh_new{iC}{iCC}); dta_out(:,new_col) = dta_out_hld; dta_out_col{new_col} = col_tgh_new{iC}{iCC}; end
        
    end
end

end