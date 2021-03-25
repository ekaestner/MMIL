function outmat = default_m(ask_str,fin,fout);
% function outmat = default_m(ask_str,fin,fout);
% Obtain matrix of data from user, using defaults in fin, if non-empty, and
% writing response to fout, if non-empty.
% fin and fout must have been previously opened by the user.
% To accept the default data in brackets, user simply hits return.
% If brackets are empty, or if user keys in ']', then user input is
% terminated.

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
  defaultmat = eval(['[ ' line ' ]']); % convert to data
  mdef = defaultmat(1);		% rows
  ndef = defaultmat(2);		% columns
  defaultmat(1:2) = [];		% leave data
  defaultmat = reshape(defaultmat,mdef,ndef);
  
else				% no default data exists
  
  mdef = 0;
  ndef = 0;
  defaultmat = [];			% no default available

end

outmat = [];			% accumulate user's input
outstr = [];			% initialize user's response
disp(['('']'' to terminate entries early, <return> to accept default):'])
disp(' ')
disp([ask_str ':'])

while(strcmp(outstr,']')==0),	% as long as any of the chars are not ']'

  [mout,nout] = size(outmat);	% present size of user's input

  if((mout+1) <= mdef),		% default data available

    outstr = input([sprintf('Row %.0f [',mout+1) ...
	sprintf(' %4g',defaultmat(mout+1,:)) ']: '],'s');

    if(isempty(outstr)),		% user wants the default
      outmat = [outmat;defaultmat(mout+1,:)];
    elseif(outstr(1)~=']'),	% user gave something else
      next_vec = eval(['[ ' outstr ']']); % interpret
      [mnext,nnext] = size(next_vec); % size of next user input
      if((nnext ~= nout) & (mout)), % wrong width and not first row
	disp('Inconsistent number of columns, try again.')
      else
	outmat = [outmat;next_vec]; % accumulate
      end
    end

  else				% default data not available

    outstr = input(sprintf('Row %.0f []:',mout+1),'s');

    if(isempty(outstr)),	% user gave nothing assume done
      outstr = ']';		% signal to while loop
    elseif(outstr(1)~=']'),	% user gave something else
      next_vec = eval(['[ ' outstr ']']); % interpret
      [mnext,nnext] = size(next_vec); % size of next user input
      if((nnext ~= nout) & (mout)), % wrong width and not first row
	disp('Inconsistent number of columns, try again.')
      else
	outmat = [outmat;next_vec]; % accumulate
      end
    end

  end				% ifelse default data

end				% while user inputs data

[mout,nout] = size(outmat);	% final size

if(~isempty(fout)),		% user wants to write it out as well
  fprintf(fout,' %g',mout,nout,outmat(:));
  fprintf(fout,'\n');		% end the line
end

return
