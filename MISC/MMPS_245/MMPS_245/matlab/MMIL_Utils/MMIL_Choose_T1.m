function [fname_T1, errcode] = MMIL_Choose_T1(ContainerPath,varargin)
%function [fname_T1, errcode] = MMIL_Choose_T1(ContainerPath,[options])
%
% Required Parameters:
%  ContainerPath: full path of MRIPROC Container
%
% Optional Parameters:
%  'fname_T1': full path file fname of T1 weighted image
%     if empty, will find appropriate T1 image in ContainerPath
%     {default = []}
%  'T1type': which T1 series to use
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%  'FSContainerPath': full path of directory containing freesurfer recon
%     If supplied and fname_T1 is empty, will use nu.mgz
%     {default = []}
%  'bbregflag': [0|1] for FreeSurfer's bbregister, use orig.mgz
%     (requires FSContainerPath)
%     {default = 0}
%  'ext': file extension (i.e. '.mgh' or '.mgz')
%     {default = '.mgz'}
%
% Created:  11/02/11 by Don Hagler
% Prev Mod: 02/03/15 by Don Hagler
% Last Mod: 10/26/17 by Feng Xue (To allow tar as input)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = mmil_args2parms(varargin, { ...
  'fname_T1',[],[],...
  'T1type',2,[0,1,2,3],... 
  'FSContainerPath',[],[],...
  'bbregflag',false,[false true],...
  'ext','.mgz',{'.mgh','.mgz'},...
});

fname_T1 = [];
errcode = 0;

switch parms.T1type
  case 0, fnames_T1 = {['MPR_res' parms.ext]};
  case 1, fnames_T1 = {['hiFA_res' parms.ext]};
  case 2, fnames_T1 = {['MPR_res' parms.ext]; ['hiFA_res' parms.ext]};
  case 3, fnames_T1 = {['hiFA_res' parms.ext]; ['MPR_res' parms.ext]};
end

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
if ~isempty(parms.FSContainerPath) && ~exist(ContainerPath,'dir')
  fprintf('%s: ERROR: FreeSurfer container %s not found\n',...
    mfilename,parms.FSContainerPath);
  errcode = 1;
  return;
end;

if parms.bbregflag && isempty(parms.FSContainerPath)
  fprintf('%s: ERROR: bbregflag=1 and FSContainerPath is empty\n',mfilename);
  errcode = 1;
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose T1 file
if parms.bbregflag
  fname_T1 = [parms.FSContainerPath '/mri/orig.mgz'];
elseif ~isempty(parms.fname_T1)
  fname_T1 = parms.fname_T1;
elseif ~isempty(parms.FSContainerPath) % if recon path specified, reg to Freesurfer T1
  fname_T1 = sprintf('%s/mri/nu.mgz',parms.FSContainerPath);
else % if recon path not specified, use processed T1 (MPR or hiFA)
  foundflag = 0;
  for i=1:length(fnames_T1)
    fname_test = sprintf('%s/%s',ContainerPath,fnames_T1{i});
    if isempty(ContainerPath_tar)
      if ~exist(fname_test,'file'), continue; end;
    else
      [~,fname_tar,ext] = fileparts(ContainerPath_tar);
      [err,result] = unix(sprintf('tar tf %s --occurrence %s/%s',ContainerPath_tar,fname_tar,fnames_T1{i}));
      if err, continue; end;
    end
    foundflag = 1;
    fname_T1 = fname_test;
    break;
  end
  if ~foundflag
    fprintf('%s: ERROR: hiFA or MPR series not found in %s\n',...
      mfilename,ContainerPath);
    errcode = 1;
    return;
  end
end;
if isempty(ContainerPath_tar)
  if ~exist(fname_T1,'file')
    fprintf('%s: ERROR: T1 file %s not found\n',mfilename,fname_T1);
    errcode = 1;
    return;
  end;
end
