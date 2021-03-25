function info = mmil_read_HEAD_statsinfo(fname)
%function info = mmil_read_HEAD_statsinfo(fname)
%
% Required Input:
%   fname: file name of AFNI BRIK or HEAD file
%          should be an fbuc, e.g. output from 3dDeconvolve
%
% Output:
%   info: struct array containing labels and degrees of freedom for each
%         frame of fname
%
% 
% Created:  09/01/08 by Don Hagler
% Prev Mod: 10/27/17 by Dani Cornejo 
% Last Mod: 10/31/17 by Don Hagler
%

% based on code from BOLD_MMIL_Fourier_Analysis

if (~mmil_check_nargs(nargin,1)) return; end;
info = [];

[fpath,fstem,fext] = fileparts(fname);
if strcmp(fext,'.BRIK')
  fext = '.HEAD';
elseif ~strcmp(fext,'.HEAD')
  error('fname %s is neither BRIK nor HEAD file',fname);
end;
fname = sprintf('%s/%s%s',fpath,fstem,fext);
if ~exist(fname,'file')
  error('file %s not found',fname);
end;

fid=fopen(fname);
dofs = {};
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break; end;
  if strcmp(tline,'name = BRICK_LABS')
    tline = fgetl(fid);
    tline = fgetl(fid);
%    disp(tline);
    tline = tline(2:end); % get rid of initial '
    ind = findstr(tline,'~');
    k = 1;
    labels = {};
    for j=1:length(ind)
      labels{end+1} = tline(k:ind(j)-1);
      k=ind(j)+1;
    end;
  elseif strcmp(tline,'name = BRICK_STATSYM')
    tline = fgetl(fid);
    tline = fgetl(fid);
%    disp(tline);
    tline = tline(2:end-1); % get rid of initial ' and ending ~
    tline(end+1) = ';';
    ind = findstr(tline,';');
    k = 1;
    for j=1:length(ind)
      tmp = tline(k:ind(j)-1);
      k=ind(j)+1;
      n = regexp(tmp,'Ftest\((?<dof1>\d+),(?<dof2>\d+)\)','names');
      if ~isempty(n)
        dofs{end + 1} = [str2num(n.dof1),str2num(n.dof2)];
      else
        n = regexp(tmp,'Ttest\((?<dof>\d+)\)','names');
        if ~isempty(n)
          dofs{end + 1} = str2num(n.dof);
        else
          dofs{end + 1} = [];
        end;
      end;
    end;
  end;
end;
fclose(fid);

if length(labels) ~= length(dofs)
  error('number of labels did not match number of dofs');
end;

for i=1:length(labels)
  label = labels{i};
  info(i).name = label;
  if ~isempty(regexp(label,'Fstat'))
    info(i).type = 'fstat';
  elseif ~isempty(regexp(label,'Tstat'))
    info(i).type = 'tstat';
  elseif  ~isempty(regexp(label,'Coef'))
    info(i).type = 'coef';
  end;
  info(i).dofs = dofs{i};
end;


