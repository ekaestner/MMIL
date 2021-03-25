function [fname_T2, errcode] = MMIL_Choose_T2(ContainerPath,varargin)
%function [fname_T2, errcode] = MMIL_Choose_T2(ContainerPath,[options])
%
% Required Parameters:
%  ContainerPath: full path of MRIPROC Container
%
% Optional Parameters:
%  'fname_T2': full path file fname of T2 weighted image
%     if empty, will find appropriate T2 image in ContainerPath
%     {default = []}
%  'fstem_T2': file stem of T2 weighted image
%     {default = 'T2w_res'}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%     {default = '.mgz'}
%
% Created:  05/11/17 by Don Hagler
% Prev Mod: 05/11/17 by Don Hagler
% Last Mod: 10/26/17 by Feng Xue (To allow tar as input)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_T2 = []; errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin, { ...
  'fname_T2',[],[],...
  'fstem_T2','T2w_res',[],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if proc container exists
ContainerPath_tar = [];
switch exist(ContainerPath)
  case 2
    ContainerPath_tar = ContainerPath;
    ContainerPath = regexprep(ContainerPath,'\.tar$','');
  case 0
    fprintf('%s: ERROR: processed container %s not found\n',...
      mfilename,ContainerPath);
    errcode = 1;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose T2 file
if ~isempty(parms.fname_T2)
  fname_T2 = parms.fname_T2;
else
  fname_T2 = sprintf('%s/%s%s',ContainerPath,parms.fstem_T2,parms.ext);
end;
if isempty(ContainerPath_tar)
  if ~exist(fname_T2,'file')
    fprintf('%s: ERROR: T2 file %s not found\n',mfilename,fname_T2);
    errcode = 1;
    return;
  end;
else
  [~,fname_tar,ext] = fileparts(ContainerPath_tar);
  [~,fstem_T2,ext] = fileparts(fnames_T2);
  [err,result] = unix(sprintf('tar tf %s --occurrence %s/%s%s',ContainerPath_tar,fname_tar,fstem_T2,ext));
  if err, continue; end;
end
