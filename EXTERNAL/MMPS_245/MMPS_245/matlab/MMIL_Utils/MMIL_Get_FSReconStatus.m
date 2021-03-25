function [status,message,errcode] = MMIL_Get_FSReconStatus(ContainerPath,FS_version,putflag);
%function [status,message,errcode] = MMIL_Get_FSReconStatus(ContainerPath,[FS_version],[putflag]);
%
% Required input:
%   ContainerPath: full path of freesurfer recon directory (for one subject)
%
% Optional parameters:
%  FS_version: freesurfer version number (e.g. 305, 450, etc.)
%    Only used if ContainerInfo.FS_version is not available
%    {default = 450}
%   putflag: [0|1] whether to write status to EVENT.xml file
%     {default: 0}
%
% status codes: messages
%  0: recon not started
%  1: recon started
%  2: recon complete
%  3: waiting for re-recon
%  4: re-recon started
%  5: re-recon complete
%  6: volume-only recon complete
%
% Created:  01/11/07 Don Hagler
% Rcnt mod: 11/04/11 Don Hagler
% Last Mod: 09/15/12 by Don Hagler
%

% check parameters
if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('putflag','var') | isempty(putflag), putflag = 0; end;
if ~exist('FS_version','var') | isempty(FS_version), FS_version = 450; end;

errcode = 0;
status = 0;
message = [];

if isempty(ContainerPath), return; end;
[ContainerRootDir,tmp_name,tmp_ext]=fileparts(ContainerPath);
ContainerDir = [tmp_name tmp_ext];
if isempty(ContainerDir)
  status = 0;
  message = 'recon not started';
  return;
end;  

% load ContainerInfo, get FS_version field
[FS_version,errcode] = MMIL_Get_FS_version(ContainerPath,FS_version);

[allexist,volexist] = fs_check_status(ContainerDir,ContainerRootDir,FS_version);

if allexist
  if ~isempty(dir(sprintf('%s/touch/remake*',ContainerPath)))
    status = 3;
    message = 'waiting for re-recon';
  elseif ~isempty(dir(sprintf('%s/touch/remade*',ContainerPath)))
    status = 5;
    message = 're-recon complete';
  else
    status = 2;
    message = 'recon complete';
  end;
else
  if ~isempty(dir(sprintf('%s/touch/remade*',ContainerPath)))
    status = 4;
    message = 're-recon started';
  elseif volexist
    status = 6;
    message = 'volume-only recon complete';
  else
    status = 1;
    message = 'recon started';
  end;
end;

if putflag
  errcode = MMIL_PutEvent_FSReconStatus(ContainerPath,status,message);
end;
