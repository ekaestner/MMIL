function red_cap = ejk_load_redcap_01_load(cfg)

red_cap = mmil_readtext(cfg.red_fle,cfg.sep);
for iC = 1:size(red_cap,2)
    if isstr(red_cap{2,iC})
        red_cap(cellfun(@isempty,red_cap(:,iC)),iC) = repmat({''},sum(cellfun(@isempty,red_cap(:,iC))),1);
    else
        red_cap(cellfun(@isempty,red_cap(:,iC)),iC) = repmat({nan},sum(cellfun(@isempty,red_cap(:,iC))),1);
    end    
end


end