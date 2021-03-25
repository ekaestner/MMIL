function errcode = BOLD_MMIL_Export_Data_for_FSL(ContainerPath,FSContainerPath,varargin)
%function errcode = BOLD_MMIL_Export_Data_for_FSL(ContainerPath,FSContainerPath,[options])
%
% Purpose: Convert processed BOLD volumes to a nii
%   with customized steps to prepare for use with FSL
%
% Usage:
%  BOLD_MMIL_Export_Data_for_FSL(ContainerPath,FSContainerPath,'key1', value1,...);
%  e.g. BOLD_MMIL_Convert_for_FSL(ContainerPath,FSContainerPath,...
%         'snums_export',[3,4],'infix','corr');
%
% Required Input:
%  ContainerPath: Full path of BOOLDROC directory containing BOLD scans
%  FSContainerPath: full path of FSURF directory
%
% Optional Paramters:
%  'outdir': output directory
%    provide full path or will be relative to ContainerPath
%    {default = 'exportBOLDforFSL'}
%  'snums_export': vector of scan numbers to export
%    if empty, will export all valid BOLD scans in ContainerPath
%    {default = []}
%  'snums_valid': vector of scan numbers that were processed
%    if empty, will use snums_valid returned by BOLD_MMIL_Get_ScanInfo
%    {default = []}
%  'infix': BOLD file name infix (e.g. '', 'corr')
%    {default = []}
%  'skipTRs' - number of TRs at beginning of scan to be ignored in time shifting
%    {default = 0}
%  'smooth': FWHM blurring kernel applied to BOLD data
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  11/21/11 by Don Hagler
% Last Mod: 11/21/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir','exportBOLDforFSL',[],...
  'snums_export',[],[],...
  'snums_valid',[],[],...
  'infix',[],[],...
  'skipTRs',0,[0,100],...
  'smooth',0,[0 100],...
  'forceflag',false,[false true],...
...
  'fnamestem','BOLD',[],...
});
errcode = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(ContainerPath)
  fprintf('%s: ERROR: ContainerPath is empty\n',mfilename);
  errcode = 1;
  return;
elseif ~exist(ContainerPath,'file')
  fprintf('%s: ERROR: ContainerPath %s not found\n',mfilename,ContainerPath);
  errcode = 1;
  return;
end;
if isempty(FSContainerPath)
  fprintf('%s: ERROR: FSContainerPath is empty\n',mfilename);
  errcode = 1;
  return;
elseif ~exist(FSContainerPath,'file')
  fprintf('%s: ERROR: FSContainerPath %s not found\n',mfilename,FSContainerPath);
  errcode = 1;
  return;
end;

% set fname_mask and fname_T1    
[subjdir,subj,tmp] = fileparts(FSContainerPath);
subj = [subj tmp];
fname_T1 = sprintf('%s/mri/nu.mgz',FSContainerPath);
fname_mask = sprintf('%s/mri/aseg.mgz',FSContainerPath);
if ~exist(fname_T1,'file')
  fprintf('%s: ERROR: file %s not found\n',mfilename,fname_T1);
  errcode = 1;
  return;
end;
if ~exist(fname_mask,'file')
  fprintf('%s: ERROR: file %s not found\n',mfilename,fname_mask);
  errcode = 1;
  return;
end;

% load scan info
[ScanInfo,SessInfo,errcode] = BOLD_MMIL_Get_ScanInfo(ContainerPath,...
  'snums',parms.snums_valid,'fnamestem',parms.fnamestem);
if errcode || isempty(ScanInfo), return; end;
if isempty(parms.snums_export)
  parms.snums_export = SessInfo.snums_valid;
end;

if mmil_isrelative(parms.outdir)
  parms.outdir = [ContainerPath '/' parms.outdir];
end;
mmil_mkdir(parms.outdir);

% convert each volume
for i=1:length(parms.snums_export)
  s = parms.snums_export(i);
  nreps = ScanInfo(s).nreps;
  if nreps<=parms.skipTRs+2, continue; end;
  TR = ScanInfo(s).TR;  
  % check for bad scan num
  if s<1 | s>SessInfo.nscans
    fprintf('%s: ERROR: bad BOLD Scan Num (%d)\n',mfilename,s);
    errcode = 1;
    return;
  end;
  fstemlist = ScanInfo(s).fstemlist;
  for f=1:length(fstemlist)
    fstem = fstemlist{f};
    if ~isempty(parms.infix)
      fstem = [fstem '_' parms.infix];
    end;
    outdir = [parms.outdir '/' fstem];
    fname_out = sprintf('%s/summary.txt',outdir);
    if ~exist(fname_out,'file') || parms.forceflag
      % set fname_BOLD
      fname_BOLD = sprintf('%s/%s.mgz',ContainerPath,fstem);
      if ~exist(fname_BOLD,'file')
        fprintf('%s: ERROR: file %s not found\n',mfilename,fname_BOLD);
        errcode = 1;
        return;
      end;
      % set fname_motion
      fname_motion = sprintf('%s/%s_motion.1D',ContainerPath,fstem);
      if ~exist(fname_motion,'file')
        fprintf('%s: ERROR: motion 1D file %s not found\n',...
          mfilename,fname_motion);
        errcode = 1;
        return;
      end;
      cmd = sprintf('fsl_fmri_preproc.csh %s %s %s %s %s %0.3f %d %s %s %d',...
        fname_BOLD,fname_mask,fname_T1,outdir,fname_motion,...
        TR/1000,parms.skipTRs,subjdir,subj,parms.smooth);
      fprintf('%s: exporting %s for FSL...\n',mfilename,fstem);
      fprintf('%s\n',cmd);
      [s,r] = mmil_unix(cmd);
      if s
        fprintf('%s: cmd %s failed:\n%s\n',mfilename,cmd,r);
        errcode = 1;
        return;
      end;
    end;
  end;
end;

