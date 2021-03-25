function allopt = ts_read_excel_setup(optfile,varargin)
% options = ts_read_excel_setup(file)
% options = ts_read_excel_setup(file, write_flag)
% options = ts_read_excel_setup(file, sheet)
% options = ts_read_excel_setup(file, sheet, write_flag)

% Created by Jason Sherfey on 19-Nov-2008

% opt = mmil_args2parms(varargin, { ...
%   'write_flag',0,[],...
%   'sheet_flag',0,[],...
% });

% look for spreadsheet
if ~exist('optfile','var') || ~exist(optfile,'file')
	[filename,pathname] = uigetfile('.xls');
	optfile = fullfile(pathname,filename);
end

% set parameters
write_flag  = 0;
sheet_flag  = 0;
switch nargin
    case 2
        if isstr(varargin{1}),  sheet_flag = 1; 
        else                    write_flag = all(varargin{1});
        end
    case 3
        if isstr(varargin{1}),  sheet_flag = 1; 
                                write_flag = all(varargin{2});
        else                    sheet_flag = 2;
                                write_flag = all(varargin{1});
        end
end

warning('off','MATLAB:xlsread:Mode')
if sheet_flag
    [num,txt,raw]=xlsread(optfile,varargin{sheet_flag});
else
    [num,txt,raw]=xlsread(optfile);			% read excel spreadsheet
end
warning('on','MATLAB:xlsread:Mode')

if write_flag
	[outpath,outname] = fileparts(optfile);
	outname = ['auto_' outname '.m'];
	outfile = fullfile(outpath,outname);
	fid = fopen(outfile,'wt');
else
    fid = 0;
end

n 		= size(raw,1); 							% # rows
allopt 	= [];

for i=1:n														% collect all options in one structure
	allopt = update(allopt,raw(i,:),fid);
end

if fid, fclose(fid); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opt = update(opt,col,fid)

if 			length(col) < 3, 	return;		% not enough info
elseif 	isnan(col{1}),		return;		% no option name
elseif 	isnan(col{3}), 		return;		% no  value
end

nam = col{1}; if isstr(nam) && any(regexp(nam,'\W')), return; end
fld = col{2}; if isstr(fld) && any(regexp(fld,'\W')), return; end
val = col{3}; if isnan(val), val=[]; end;  
valquotes_flag = 0;

% check if val is 'pwd'
if isstr(val) && strcmp(val,'pwd'), eval('val = pwd;'); end

if fid
	strval = val;
	if isnumeric(val), strval = num2str(val); end
end

% new field
if 			isstr(nam) && isstr(fld), str = [nam '.' fld];
elseif	isstr(nam) && isnan(fld), str = nam;
else		return;
end

% convert strings of cells and arrays into cells and arrays
if isstr(val)
% 	if length(regexp(val,'[{^}$]'))==2
if ~isempty(regexp(val,'^{.*}','match'))
		% val is a cell
		val = regexp(val,'[^,{}]+','match');
		% convert strings of arrays in cells to arrays in cells
		for k=1:length(val)
			if 	length(regexp(val{k},'[\[^\]$]'))==2
				% val{k} is an array
				val{k} = str2num(val{k});
			else
				% val{k} is a string in quotes
				val(k) = regexp(val{k},'[^'']+','match');
			end
		end
	elseif length(regexp(val,'[\[^\]$]'))==2
					% val is an array
		try val = str2num(val); end;
	else
		valquotes_flag = 1;
	end
end
eval(['opt.' str '=val;']);

if fid  && valquotes_flag
	fprintf(fid,'%s = ''%s'';\n',str,strval);
elseif fid
	fprintf(fid,'%s = %s;\n',str,strval);
end
	
