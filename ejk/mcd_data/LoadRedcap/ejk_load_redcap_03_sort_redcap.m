function sbj_out = ejk_load_redcap_03_sort_redcap(red_cap,var_nme)

for iN = 1:size(var_nme,1)
    
    if strcmpi(var_nme{iN,3},'cell')
        sbj_out.(var_nme{iN,1}) = cell(size(red_cap,1)-1,1);
    elseif strcmpi(var_nme{iN,3},'nan')
        sbj_out.(var_nme{iN,1}) = nan(size(red_cap,1)-1,1);
    end
    
    if ~strcmpi(var_nme{iN,2},'')
        
        if ischar(var_nme{iN,2})
            sbj_out_col = strcmpi(red_cap(1,:),var_nme{iN,2});
        else
            sbj_out_col = ismember(red_cap(1,:),var_nme{iN,2});
        end
        
        if sum(sbj_out_col)>0
            if strcmpi(var_nme{iN,3},'cell')
                for iS = 2:size(red_cap,1)
                    sbj_out.(var_nme{iN,1}){iS-1,1} = red_cap{iS,sbj_out_col};
                end
            end
            if strcmpi(var_nme{iN,3},'nan')
                for iS = 2:size(red_cap,1)
                    sbj_out.(var_nme{iN,1})(iS-1,1) = red_cap{iS,sbj_out_col};
                end
            end
        end
        
    end
end

end