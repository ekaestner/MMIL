function found_numbers = string_find(string_cell,string_find,varargin)

string_cell(cellfun(@isempty,string_cell)) = {'`'};
if isempty(varargin)
found_numbers = find(~cellfun(@isempty,regexp(upper(string_cell),upper(string_find))));
else
    found_numbers = ~cellfun(@isempty,regexp(upper(string_cell),upper(string_find)));
end

end