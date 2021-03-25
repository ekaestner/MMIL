function new_fle = order_sbj_inf(clr_fle)

bse_sbj_inf = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/bse_sbj_inf');
new_fle     = cell(numel(bse_sbj_inf),1);

fld_nme = cellfun(@(x) strfind(x,'|'),clr_fle,'uni',0);
fld_nme = cellfun(@(x,y) x(1:y-2),clr_fle,fld_nme,'uni',0);

for iR = 1:size(bse_sbj_inf,1)
    if ~isempty(bse_sbj_inf{iR});
        ind = find(strcmpi(bse_sbj_inf{iR},fld_nme));
        if ~isempty(ind)
            new_fle{iR,1} = clr_fle{ind,1};
        else
            new_fle{iR,1} = [bse_sbj_inf{iR} ' | []'];
        end
    end
end

end