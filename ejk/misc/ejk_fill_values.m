function dta = ejk_fill_values(dta)

for iC = 1:size(dta,2)
    cls_tbl = tabulate(cellfun(@class,dta(~cellfun(@isempty,dta(:,iC)),iC),'uni',0));
    if sum(cls_tbl(:,3)>20)>1; error('inconsistent classes, please check'); end
    [ ~, max_ind] = max(cell2mat(cls_tbl(:,2)));
    if strcmpi(cls_tbl(max_ind,1),'double')
        dta(cellfun(@isempty,dta(:,iC)),iC) = {NaN};
    elseif strcmpi(cls_tbl(max_ind,1),'char')
        dta(cellfun(@isempty,dta(:,iC)),iC) = {''};
    else
        error('unknown class, please add to function')        
    end
end


end