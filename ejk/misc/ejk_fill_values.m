function dta = ejk_fill_values(dta,varargin)

for iC = 1:size(dta,2)
    cls_tbl = tabulate(cellfun(@class,dta(~cellfun(@isempty,dta(:,iC)),iC),'uni',0));
    if sum(cell2mat(cls_tbl(:,3))>20)>1 && (isempty(varargin) || ~varargin{1}(iC))
        cls_tbl
        error(['inconsistent classes in column ' num2str(iC) ', please check']); 
    end
    [ ~, max_ind] = max(cell2mat(cls_tbl(:,2)));
    if strcmpi(cls_tbl(max_ind,1),'double')
        dta(cellfun(@isempty,dta(:,iC)),iC) = {NaN};
    elseif strcmpi(cls_tbl(max_ind,1),'char')
        dta(cellfun(@isempty,dta(:,iC)),iC) = {''};        
    elseif isempty(cls_tbl)
        sprintf(['Column ' num2str(iC) ' is empty\n'])
    else
        error(['unknown class in column ' num2str(iC) ', please add to function'])        
    end
end


end