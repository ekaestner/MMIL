function [strlist] = rmquotes(strlist)

sz = size(strlist);

A = strlist(:);
for ii=1:length(A),
   tstr = A{ii};
   tstr(strfind(tstr,'"')) = '';
   A(ii) = {tstr};
end
strlist = reshape(A,sz);
