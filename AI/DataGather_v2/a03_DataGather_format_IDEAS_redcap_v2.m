age_col_nme = { 'Binned_Age_at_Scan' 'Binned_Onset_Age' 'Binned_Age_at_Surgery' };

for iAC = 1:numel(age_col_nme)
    age_col_ind = strcmpi(age_col_nme{iAC},red_cap_col);
    for iS = 1:numel(red_cap_sbj)
        % Get the binned age string
        age_str = red_cap_dta{iS,age_col_ind};
        if ischar(age_str) || (isstring(age_str) && ~isempty(age_str))
            tokens = regexp(age_str, '^(\d+)_to_(\d+)$', 'tokens');
            if ~isempty(tokens)
                nums = str2double(tokens{1});
                red_cap_dta{iS,age_col_ind} = mean(nums);
            else
                if strcmpi(age_str,'Over_55')
                    red_cap_dta{iS,age_col_ind} = 60;
                elseif strcmpi(age_str,'Less_than_20')
                    red_cap_dta{iS,age_col_ind} = 16;
                elseif strcmpi(age_str,'Less_than_1')
                    red_cap_dta{iS,age_col_ind} = 0.5;
                elseif strcmpi(age_str,'Over_40')
                    red_cap_dta{iS,age_col_ind} = 45;
                else
                    red_cap_dta{iS,age_col_ind} = NaN;
                end
            end
        else
            red_cap_dta{iS,age_col_ind} = NaN;
        end
    end
end