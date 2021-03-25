function DTI_MMIL_Convert(ContainerPath,varargin)
%function DTI_MMIL_Convert(ContainerPath,[options])
%
% Purpose: Convert processed DTI volumes to a different format
%
% Usage:
%  DTI_MMIL_Convert(ContainerPath,'key1', value1,...);
%  e.g. DTI_MMIL_Convert(ContainerPath,...
%         'snums',[3,4],'infix','B0uw_gruw_mc');
%
% Required Input:
%  ContainerPath: Full path of DTIPROC directory containing DTI scans
%
% Optional Paramters:
%  'out_type' - output file type
%    supported types: 'nii', 'BRIK'
%    { default: 'nii' }
%  'snums' - vector of scan numbers to convert
%    if empty or omitted, will convert all DTI scans in ContainerPath
%    { default: [] }
%  'infix' - DTI file name infix (e.g. 'corr_resDTI')
%     if empty, will look for files like 'DTI1.mgz'
%     {default = []}
%  'out_orient' - output slice orientation
%    if empty or omitted, keep original orientation
%    e.g. 'LPS', 'RAS', etc.
%    For use with FSL, use 'LPS'
%    { default: []}
%  'forceflag' - [0|1] whether to run conversions even if output files exist
%    { default: 0 }
%
% Created:  05/17/07 by Don Hagler
% Rcnt Mod: 03/23/11 by Don Hagler
% Last Mod: 09/10/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'out_type','nii',{'nii','BRIK'},...
  'snums',[],[],...
  'infix',[],[],...
  'out_orient',[],[],...
  'forceflag',false,[false true],...
...
  'fnamestem','DTI',[],...
});

switch parms.out_type
  case 'nii'
    parms.out_ext = '.nii';
  case 'BRIK'
    parms.out_ext = '+orig.BRIK';
  otherwise
    error('invalid out_type %s',parms.out_type);
end;

% load scan info
[ScanInfo,SessInfo,errcode] = DTI_MMIL_Get_ScanInfo(ContainerPath);

% check for bad scan nums
if isempty(parms.snums), parms.snums=[1:length(ScanInfo)]; end;
for i=1:length(parms.snums)
  snum = parms.snums(i);
  if snum<1 | snum>length(ScanInfo)
    error('bad scan num: %d',snum);
  end;
end;

% convert each scan
for i=1:length(parms.snums)
  snum = parms.snums(i);
  if ~isempty(parms.infix)
    fstems = {...
      sprintf('%s%d_%s',parms.fnamestem,snum,parms.infix),...
      sprintf('%s%d_rev_%s',parms.fnamestem,snum,parms.infix)};
  else
    fstems = {...
      sprintf('%s%d',parms.fnamestem,snum),...
      sprintf('%s%d_rev',parms.fnamestem,snum)};
  end;
  for f=1:length(fstems)
    fname_in  = sprintf('%s/%s.mgz',ContainerPath,fstems{f});
    if ~exist(fname_in,'file')
%      fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
      continue;
    end;
    fname_out = sprintf('%s/%s%s',ContainerPath,fstems{f},parms.out_ext);
    if ~exist(fname_out,'file') || parms.forceflag
      fprintf('%s: converting mgz to %s for scan %d\n',...
        mfilename,parms.out_type,snum);
      switch parms.out_type
        case 'nii'
          fs_mri_convert(fname_in,fname_out,...
            'out_orient',parms.out_orient,'forceflag',parms.forceflag);
        case 'BRIK'
          mmil_mgh2BRIK(fname_in,fname_out,[],parms.out_orient,parms.forceflag);
      end;
    end;
  end;
end;

