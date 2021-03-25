function mmil_mriRegister(fname_volA,fname_volB,varargin)
%function mmil_mriRegister(fname_volA,fname_volB,varargin)
%
% Purpose: perform rigid, affine, or nonlinear registration
%   between two MRI volumes
%   This is a wrapper for Dominic Holland's mriRegister
%
% Usage: mmil_mriRegister(fname_volA,fname_volB,'key1',val1,...)
%
% Required Input:
%   fname_volA: file name of volume A (baseline)
%   fname_volB: file name of volume B (followup)
%
% Optional Input:
%   'outdir' - output path
%     {default = pwd}
%   'fname_maskA' - full path of brain mask for volume A
%     If empty or does not exist, will generate one using mmil_dct_brainmask
%     {default = outdir/brain_maskA.mgh}
%   'fname_maskB' - full path of brain mask for volume B
%     If empty or does not exist, will generate one using mmil_dct_brainmask
%     {default = outdir/brain_maskB.mgh}
%   'fname_asegA' - full path of FreeSurfer aseg for volume A
%     If supplied, will create brain mask A from aseg (and save as fname_maskA)
%     using  mmil_dilate_mask
%     {default = []}
%   'fname_asegB' - full path of FreeSurfer aseg for volume B
%     If supplied, will create brain mask B from aseg (and save as fname_maskB)
%     using  mmil_dilate_mask
%     {default = []}
%   'fname_affine_mat' - full path of input affine registration matrix file
%     e.g. created by previous run of mriRegister
%     {default = []}
%   'fname_log' - file name of output log file (relative to outdir)
%     {default = 'reg.log'}
%   'paramfile' - full or relative path of parameter file
%     If relative, relative to $MMPS_PARMS/MRIREG
%     {default = 'inputParamsNonlin.txt'}
%   'binfile' - full or relative path of mriRegister binary file
%     {default = 'mriRegister'}
%   'forceflag' - [0|1] overwrite existing output files
%     {default = 0}
% 
% See Also:
%   mmil_dilate_mask: to create mask from directly specified aseg volume
%   mmil_dct_brainmask: to create mask from image volume using dct morph to atlas
%   mmil_aseg_brainmask: to create mask from FreeSurfer aseg
%
% Created:  08/17/09 by Don Hagler
% Last Mod: 04/23/10 by Don Hagler
%

%% todo: add parameters from paramfile or mriRegister command line as optional inputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,2)) return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',[],[],...
  'fname_maskA',[],[],...
  'fname_maskB',[],[],...
  'fname_asegA',[],[],...
  'fname_asegB',[],[],...
  'fname_affine_mat',[],[],...
  'fname_log','reg.log',[],...
  'paramfile','inputParamsNonlin.txt',[],...
  'binfile','mriRegister',[],...
  'forceflag',false,[false true],...
...
  'abm_smooth1',20,[0,100],... % smoothing for aseg brainmask (mm)
  'abm_thresh1',0.5,[0,1],... % threshold for aseg brainmask
  'abm_smooth2',40,[0,100],... % smoothing for aseg brainmask (mm)
  'abm_thresh2',0.2,[0,1],... % threshold for aseg brainmask
  'abm_smooth3',10,[0,100],... % smoothing for aseg brainmask (mm)
  'dbm_smooth',20,[0,100],... % smoothing for dct brainmask (mm)
  'roiflag',false,[false true],...
  'atlasname','pitu_afm.mat',[],...
  'atlasdir',[],[],...
  'bindir',[],[],...
  'parmsdir',[],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if isempty(parms.atlasdir)
  parms.atlasdir = [getenv('MMPS_DIR') '/atlases'];
end;
if isempty(regexp(parms.atlasname,'^/')) % relative path
  parms.atlasname = [parms.atlasdir '/' parms.atlasname];
end;

if isempty(parms.bindir)
  parms.bindir = [getenv('MMPS_DIR') '/bin'];
end;
if isempty(regexp(parms.binfile,'^/')) % relative path
  parms.binfile = [parms.bindir '/' parms.binfile];
end;

if isempty(parms.parmsdir)
  parms.parmsdir = [getenv('MMPS_PARMS') '/MRIREG'];
end;
if isempty(regexp(parms.paramfile,'^/')) % relative path
  parms.paramfile = [parms.parmsdir '/' parms.paramfile];
end;

% check mriRegister binary exists
if ~exist(parms.binfile,'file')
  error('binary file %s not found',parms.binfile);
end;

% check parameter file exists
if ~exist(parms.paramfile,'file')
  error('parameter file %s not found',parms.paramfile);
end;

% set full log file name
if isempty(parms.outdir)
  parms.outdir = pwd;
end;
parms.fname_log = [parms.outdir '/' parms.fname_log];

% check input files exist
if ~exist(fname_volA,'file'), error('file %s not found',fname_volA); end;
if ~exist(fname_volB,'file'), error('file %s not found',fname_volB); end;

if ~isempty(parms.fname_affine_mat)
  if ~exist(parms.fname_affine_mat,'file')
    error('file %s not found',parms.fname_affine_mat);
  end;
end;

% check if reg already done
parms.fname_touch = sprintf('%s/.completed_reg',parms.outdir);
if exist(parms.fname_touch,'file')
  if ~parms.forceflag
    fprintf('%s: registration already complete\n',mfilename);
    return;
  else
    cmd = ['rm ' parms.fname_touch];
    [status,result] = unix(cmd);
    if status, error('cmd %s failed:\n%s',cmd,result); end;
  end;
end;

% create output dir if needed
if ~exist(parms.outdir,'dir')
  [succ,msg] = mkdir(parms.outdir);
  if ~succ, error('failed to create outdir %s:\n%s',parms.outdir,msg); end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create masks if needed

if ~isempty(parms.fname_asegA)
  % create from aseg
  if ~exist(parms.fname_asegA,'file')
    error('file %s not found',parms.fname_asegA);
  end
  mmil_dilate_mask([],...
    'fname_in',parms.fname_asegA,...
    'fname_out',parms.fname_maskA,...
    'smooth1',parms.abm_smooth1,'thresh1',parms.abm_thresh1,...
    'smooth2',parms.abm_smooth2,'thresh2',parms.abm_thresh2,...
    'smooth3',parms.abm_smooth3,...
    'forceflag',parms.forceflag);
elseif ~isempty(parms.fname_maskA)
  % use specified mask
  if ~exist(parms.fname_maskA,'file')
    error('file %s not found',parms.fname_maskA);
  end
else
  % create using dct morph to atlas
  mmil_dct_brainmask([],...
    'fname_in',parms.fname_volA,...
    'fname_mask',parms.fname_maskA,...
    'smooth',parms.dbm_smooth,...
    'atlasname',parms.atlasname,...
    'forceflag',parms.forceflag);
end;

if ~isempty(parms.fname_asegB)
  % create from aseg
  if ~exist(parms.fname_asegB,'file')
    error('file %s not found',parms.fname_asegB);
  end
  mmil_dilate_mask([],...
    'fname_in',parms.fname_asegB,...
    'fname_out',parms.fname_maskB,...
    'smooth1',parms.abm_smooth1,'thresh1',parms.abm_thresh1,...
    'smooth2',parms.abm_smooth2,'thresh2',parms.abm_thresh2,...
    'smooth3',parms.abm_smooth3,...
    'forceflag',parms.forceflag);
elseif ~isempty(parms.fname_maskB)
  % use specified mask
  if ~exist(parms.fname_maskB,'file')
    error('file %s not found',parms.fname_maskB);
  end
else
  % create using dct morph to atlas
  mmil_dct_brainmask([],...
    'fname_in',parms.fname_volB,...
    'fname_mask',parms.fname_maskB,...
    'smooth',parms.dbm_smooth,...
    'forceflag',parms.forceflag);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run mriRegister

cmd = sprintf('%s -ip %s -s %s -sm %s -t %s -tm %s -od %s',...
  parms.binfile,parms.paramfile,...
  fname_volA,parms.fname_maskA,...
  fname_volB,parms.fname_maskB,...
  parms.outdir);
if ~isempty(parms.fname_asegA) && parms.roiflag
  cmd = sprintf('%s -sa %s',cmd,parms.fname_asegA);
end;
if ~isempty(parms.fname_affine_mat)
  cmd = sprintf('%s -arMatrix %s',cmd,parms.fname_affine_mat);
end;

% execute command
tic
[status,result] = unix(cmd);
seconds = toc;
hours   = floor(seconds/3600);
minutes = floor(((seconds-(3600*hours))/60));
seconds = seconds-((3600*hours)+(60*minutes));
if status
  error('cmd %s failed:\n%s',cmd,result);
else
  fprintf('%s: cmd %s succeeded:\n%s',mfilename,cmd,result);
end;
fprintf('%s: Time elapsed: %02.0f:%02.0f:%02.0f.\n',...
  mfilename,hours,minutes,seconds);

% write results to log file
fid = fopen(parms.fname_log,'w');
if fid==-1
  error('failed to open log file %s for writing\n',parms.fname_log);
end;
fprintf(fid,'%s',result);
fclose(fid);

% create touch file
cmd = ['touch ' parms.fname_touch];
[status,result] = unix(cmd);
if status
  error('failed to create file %s:\n%s',parms.fname_touch,result);
end;

%% todo: cleanup?
