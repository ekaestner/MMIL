function errcode = BOLD_MMIL_Export_Data(ContainerPath,varargin)
%function errcode = BOLD_MMIL_Export_Data(ContainerPath,[options])
%
% Purpose: Convert processed BOLD volumes to a different format
%
% Usage:
%  BOLD_MMIL_Export_Data(ContainerPath,'key1', value1,...);
%  e.g. BOLD_MMIL_Convert(ContainerPath,...
%         'snums_export',[3,4],'infix','corr');
%
% Required Input:
%  ContainerPath: Full path of BOLDPROC directory containing BOLD scans
%
% Optional Paramters:
%  'outdir': output directory
%    provide full path or relative to ContainerPath
%    {default = 'exportBOLD'}
%  'snums_export': vector of scan numbers to export
%    if empty, will export all valid BOLD scans in ContainerPath
%    {default = []}
%  'snums_valid': vector of scan numbers that were processed
%    if empty, will use snums_valid returned by BOLD_MMIL_Get_ScanInfo
%    {default = []}
%  'infix': BOLD file name infix (e.g. '', 'corr')
%    {default = []}
%  'out_type': output file type
%    supported types: 'nii','mgh','mgz'
%    'nii' format can be used by FSL and AFNI
%    'BRIK' format can be used by AFNI
%    {default = 'nii'}
%  'out_orient': output slice orientation
%    if empty or omitted, keep original orientation
%    e.g. 'LPS', 'RAS', etc.
%    For use with FSL, use 'LAS'
%    {default = []}
%  'forceflag': [0|1] whether to run calculations even if output files exist
%    {default = 0}
%
% Created:  04/22/07 by Don Hagler
% Last Mod: 03/09/13 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir','exportBOLD',[],...
  'out_type','nii',{'nii','BRIK','mgh','mgz'},...
  'snums_export',[],[],...
  'snums_valid',[],[],...
  'infix',[],[],...
  'out_orient',[],[],...
  'forceflag',false,[false true],...
...
  'fnamestem','BOLD',[],...
});
errcode = 1;

switch parms.out_type
  case 'nii'
    parms.out_ext = '.nii';
  case 'BRIK'
    parms.out_ext = '+orig.BRIK';
  case 'mgh'
    parms.out_ext = '.mgh';
  case 'mgz'
    parms.out_ext = '.mgz';
  otherwise
    error('invalid out_type %s',parms.out_type);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  % check for bad scan num
  if s<1 | s>SessInfo.nscans
    fprintf('%s: ERROR: bad BOLD Scan Num (%d)\n',mfilename,s);
    errcode = 1;
    return;
  end;
  fstemlist = ScanInfo(s).fstemlist;
  for f=1:length(fstemlist)
    fstem = fstemlist{f};
    if isempty(parms.infix)
      fname_in = sprintf('%s/%s.mgz',ContainerPath,fstem);
    else
      fname_in = sprintf('%s/%s_%s.mgz',ContainerPath,fstem,parms.infix);
    end;
    if ~exist(fname_in,'file')
      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
      continue;
    end;
    fname_out = sprintf('%s/%s%s',parms.outdir,fstem,parms.out_ext);
    if ~exist(fname_out,'file') || parms.forceflag
      fprintf('%s: converting mgh to %s for scan %d\n',...
        mfilename,parms.out_type,s);
      switch parms.out_type
        case {'nii','mgh','mgz'}
          fs_mri_convert(fname_in,fname_out,...
            'out_orient',parms.out_orient,'forceflag',parms.forceflag);
        case 'BRIK'
          mmil_mgh2BRIK(fname_in,fname_out,[],parms.out_orient,parms.forceflag);
      end;
    end;
  end;
end;

