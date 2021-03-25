function outstr = default_s(ask_str,fin,fout);
% function outstr = default_s(ask_str,fin,fout);
% Obtain string data from user, using defaults in fin, if non-empty, and
% writing response to fout, if non-empty.
% fin and fout must have been previously opened by the user.

% Copyright(c) 1993 John C. Mosher
% 12/2/93 Author

if(exist('fin') ~= 1),
  fin = [];
end
if(exist('fout') ~= 1),
  fout = [];
end

if(~isempty(fin)),		% a defaults file exists
  line = fgetl(fin);		% next line in file
else
  line = [];			% no default available
end

outstr = input([ask_str '[' line ']: '],'s');

if(isempty(outstr)),		% user wants the default
  outstr = line;
end

if(~isempty(fout)),		% user wants to write it out as well
  fprintf(fout,'%s\n',outstr);
end

return
