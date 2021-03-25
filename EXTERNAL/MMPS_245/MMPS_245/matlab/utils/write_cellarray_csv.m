function write_cellarray_csv(fname,data,separator,nanflag)
%function write_cellarray_csv(fname,data,separator,nanflag)
%
% optional arguments:
%   separator: separator between columns
%     default = ','
%   nanflag: if cell is empty, write 'NaN'
%     default = 1
%

if ~exist('separator','var') | isempty(separator), separator = ','; end
if ~exist('nanflag','var') | isempty(nanflag), nanflag = 1; end

fid = fopen(fname,'w');
for i = 1:size(data,1)
  for j = 1:size(data,2)
    if j>1, fprintf(fid,'%s',separator); end
    v = data{i,j};
    if isempty(v)
      if nanflag, fprintf(fid,'NaN'); end;
    elseif isnumeric(v)
      fprintf(fid,'%f',v(1));
    else
      fprintf(fid,'"%s"',v);
    end
  end
  fprintf(fid,'\n');
end
fclose(fid);
