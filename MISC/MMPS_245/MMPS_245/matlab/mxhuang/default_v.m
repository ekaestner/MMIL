function outvec = default_v(ask_str,fin,fout);
% function outvec = default_v(ask_str,fin,fout);
% Obtain a single vector of data from user, using defaults in fin, if
% non-empty, and writing response to fout, if non-empty. fin and fout must
% have been previously opened by the user. To accept the default data in
% brackets, user simply hits return. If brackets are empty, or if user keys in
% ']', then user input is terminated.

% Copyright(c) 1993 John C. Mosher
% 12/2/93 Author

if(exist('fin') ~= 1),		% user gave no default file
  fin = [];
end
if(exist('fout') ~= 1),		% user gave no output file
  fout = [];
end

if(~isempty(fin)),		% a defaults file exists

  line = fgetl(fin);
  
else				% no default data exists
  
  line = [];			% no default available

end

outvec = input([ask_str '[' line ']: '],'s');

if(isempty(outvec)),		% user wants default
  outvec = eval(['[ ' line ' ]']); % convert to data
else
  outvec = eval(['[' outvec ']']);
end


if(~isempty(fout)),		% user wants to write it out as well
  fprintf(fout,' %g',outvec(:));
  fprintf(fout,'\n');		% end the line
end

return
