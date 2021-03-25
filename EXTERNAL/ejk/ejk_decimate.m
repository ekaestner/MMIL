% decimate returns a vector of category identifiers into the number of groups you specify 
% based on the splitting of a single continuous variable.
%
% OUTPUT
% category = index of category membership going from lowest (1) to highest
% (num_split)
% cat_check = choose to return this is you want to see the values and
% cutoffs for the category memberships
%
% INPUT
% data_vector = 
% datalist_start = the row index where the string-names begin
% 
% 
% 

function [category,cat_check] = ejk_decimate(data_vector,num_split)

% Test
% data_vector = [1 5 2 7 1 10 4 10 2 1 4 5 6 4 5 6];
% num_split = 4;

cut = 1/num_split;
cutoffs = cut:cut:1-cut;
cutoffs = quantile(data_vector,cutoffs);
cutoffs = [min(data_vector) cutoffs max(data_vector)];

cat_labels = strsplit(num2str(1:num_split),' ');
category = ordinal(data_vector,cat_labels,[],cutoffs);
category = str2num(char(category))'; %#ok<ST2NM>

[~,ind] = sort(data_vector);
category_check = category(ind); data_vector_check = data_vector(ind);
[~,ind] = sort(category_check);
cat_check = [data_vector_check(ind) category_check(ind)'];

return
end