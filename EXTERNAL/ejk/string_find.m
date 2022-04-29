function found_numbers = string_find(string_cell,string_find,varargin)

if ischar(string_find); string_find = {string_find};; end

string_cell(cellfun(@isempty,string_cell)) = {'`'};
found_numbers = [];
if isempty(varargin)
    for iW = 1:numel(string_find)
        found_numbers = [ found_numbers find(~cellfun(@isempty,regexpi(string_cell,string_find{iW}))) ];
    end
else
    found_numbers = ~cellfun(@isempty,regexp(upper(string_cell),upper(string_find)));
end

end