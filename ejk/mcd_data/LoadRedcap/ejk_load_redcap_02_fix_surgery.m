function red_cap = ejk_load_redcap_01_fix_surgery(red_cap)

srg_col = ismember( red_cap(1,:), [strcat( 'surgery_name___', cellfun(@num2str, num2cell(1:14),'uni',0)) ['surgery_name___998'] ['surgery_name___999'] ]);
for iS = 2:size(red_cap,1)
    srg_hld{iS,1} = find(cell2mat(red_cap(iS, srg_col))); 
    if isempty(srg_hld{iS,1}); srg_hld{iS,1}=NaN; end
end 
red_cap(:,srg_col) = [];
red_cap = [ red_cap [ 'surgery_name' ; srg_hld(2:end,:) ] ];

end