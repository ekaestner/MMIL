function col_nme = ejk_fix_column_names(col_nme)

for iC = 1:numel(col_nme)
    if ~isnan(num2str(col_nme{iC}(1)))
        col_nme{iC} = ['x' col_nme{iC}];
    end
end

end 