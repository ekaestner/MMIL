function tbl = ejk_switch_names(tbl,varargin)

if isempty(varargin)
    tbl_ord = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/matlab/bespoke/ft/location/mmil_split_order.csv');
else
    tbl_ord = mmil_readtext([ varargin{1} '/' 'mmil_' mmil_spec_char(varargin{2},{'_'}) '_order.csv']);
end

for iR = 1:size(tbl,1)
    
    tbl_ind = find(strcmpi(tbl_ord(:,1),tbl{iR,1}));
    if ~isempty(tbl_ind)
        tbl{iR,1} = tbl_ord{tbl_ind,2};
    end
end

end