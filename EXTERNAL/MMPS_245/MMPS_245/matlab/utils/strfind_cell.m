function resvec = strfind_cell(text,pattern)

resvec = zeros(1,length(text));
tmp = strfind(text,pattern);
for i = 1:length(text)
  resvec(i) = ~isempty(tmp{i});
end
