function errcode = MMIL_Process_First(ContainerPath,varargin)
%function errcode = MMIL_Process_First(ContainerPath,[options])
%
% Purpose: run FSL's FIRST to segment structures including hippocampus
%
% Required Parameters:
%  ContainerPath: full path of MRIPROC Container
%   also output container for firsthippo subdirectory
%
% Optional Parameters:
%  'fname_T1': full path file fname of T1 weighted image
%     if empty, will find appropriate T1 image
%       in ContainerPath or FSContainerPath
%     {default = []}
%  'T1type': which T1 series to use
%     0=MPR; 1=hiFA; 2=Either (prefer MPR); 3=Either (prefer hiFA)
%     {default=2}
%  'FSContainerPath': full path of directory containing freesurfer recon
%     Use nu.mgz instead of MPR/hiFA in ContainerPath
%     '_nu' will appended to outdir (e.g. 'firsthippo_nu')
%     {default = []}
%  'outdir': output directory, relative to ContainerPath if not full path
%    {'firsthippo'}
%  'forceflag' - [0|1] whether to overwrite existing output
%    {default = 0}
%
% Notes: Currently outputs volumes for hippocampus only
%        Requires fsl-4.1.9_RH5_64 or later
%
% Created:  11/01/11 by Don Hagler
% Rcnt Mod: 11/17/11 by Don Hagler
% Last Mod: 09/09/12 by Don Hagler
%

%% todo: modify firsthippo.csh to output volumes for all structures
%% todo: do not require FS recon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errcode = 0;
if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_T1',[],[],...
  'outdir','firsthippo',[],...
  'T1type',2,[0,1,2,3],... 
  'FSContainerPath',[],[],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temporary, absolute dependence on FS recon
if isempty(parms.FSContainerPath)
  fprintf('%s: WARNING: FSContainerPath is empty, cannot continue!\n',mfilename);
  errcode = 1;
  return;
end;

fname_aseg = sprintf('%s/mri/aseg.mgz',parms.FSContainerPath);
fname_nu = sprintf('%s/mri/nu.mgz',parms.FSContainerPath);

if isempty(parms.fname_T1) && ~isempty(parms.FSContainerPath)
  parms.outdir = [parms.outdir '_nu'];
end;

% choose T1 file
tags = {'fname_T1','T1type','FSContainerPath'};
args = mmil_parms2args(parms,tags);
[parms.fname_T1,errcode] = MMIL_Choose_T1(ContainerPath,args{:});
if errcode, return; end;

if mmil_isrelative(parms.outdir)
  parms.outdir = [ContainerPath '/' parms.outdir];
end;

fname_out = sprintf('%s/subj.volume',parms.outdir);
if exist(fname_out,'file')
  if parms.forceflag
    cmd = sprintf('rm %s/grot.nii.gz; rm %s/subj*',parms.outdir,parms.outdir);
    [s,r] = mmil_unix(cmd);
    if s
      fprintf('%s: ERROR: failed to remove existing FIRST output:\n%s\n',...
        mfilename,r);
      errcode = 1;
      return;
    end;
  else
    return;
  end;
end;

% run FIRST, specifically for hippocampal volumes
if 0
  cmd = sprintf('firsthippo.csh %s %s',parms.fname_T1,parms.outdir);
else
  cmd = sprintf('firsthippo2.csh %s %s %s %s',...
    parms.fname_T1,fname_aseg,fname_nu,parms.outdir);
end;
fprintf('%s\n',cmd);
tic;
[s,r] = mmil_unix(cmd);

if s || ~exist(fname_out,'file')
  fprintf('%s: ERROR: firsthippo.csh failed:\n%s\n',mfilename,r);
  errcode = 1;
  return;
end;
toc;

