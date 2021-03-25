function errcode = MMIL_Convert_Exam(ContainerPath,varargin)
%function errcode = MMIL_Convert_Exam(ContainerPath,[options])
%
% Purpose: Convert raw or processed data from one format to another
%
% Usage:
%  MMIL_Convert_Exam(ContainerPath,'key1', value1,...);
%
% Required Input:
%  ContainerPath: full path of MRIPROC, BOLDPROC, or DTIPROC directory
%
% Optional Paramters:
%  'outdir': output directory
%    provide full path or relative to ContainerPath
%    if empty, will write output to ContainerPath
%    {default = []}
%  'infix': file name infix (e.g. [], 'corr')
%    {default = []}
%  'in_type': input file type
%    supported types: 'nii','mgh','mgz'
%    {default = 'mgz'}
%  'out_type': output file type
%    supported types: 'nii','mgh','mgz'
%    'nii' format can be used by FSL and AFNI
%    {default = 'nii'}
%  'out_orient': output slice orientation
%    if empty or omitted, keep original orientation
%      e.g. 'LPS', 'RAS', etc.
%    for use with FSL, 'LAS' may be preferred
%    {default = []}
%  'verbose': [0|1] display status messages and warnings
%    {default: 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  02/19/16 by Don Hagler
% Last Mod: 02/20/16 by Don Hagler
%

%% todo: add support for BRIK? (need to use special functions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',[],[],...
  'in_type','mgz',{'nii','mgh','mgz'},...
  'out_type','nii',{'nii','mgh','mgz'},...
  'snums_export',[],[],...
  'infix',[],[],...
  'out_orient',[],[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
});
errcode = 0;

if ~isempty(parms.infix)
  parms.infix = ['_' parms.infix];
end;

switch parms.in_type
  case 'nii'
    parms.in_ext = '.nii';
  case 'mgh'
    parms.in_ext = '.mgh';
  case 'mgz'
    parms.in_ext = '.mgz';
  otherwise
    error('invalid in_type %s',parms.in_type);
end;

switch parms.out_type
  case 'nii'
    parms.out_ext = '.nii';
  case 'mgh'
    parms.out_ext = '.mgh';
  case 'mgz'
    parms.out_ext = '.mgz';
  otherwise
    error('invalid out_type %s',parms.out_type);
end;

if isempty(parms.outdir)
  parms.outdir = ContainerPath;
elseif mmil_isrelative(parms.outdir)
  parms.outdir = [ContainerPath '/' parms.outdir];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ContainerInfo,errcode] = MMIL_Load_ContainerInfo(ContainerPath);
if errcode || isempty(ContainerInfo), return; end;

mmil_mkdir(parms.outdir);

stypes = fieldnames(ContainerInfo.ScanInfo);
for i=1:length(stypes)
  stype = stypes{i};
  tmpinfo = ContainerInfo.ScanInfo.(stype);
  if isempty(tmpinfo), continue; end;
  for j=1:length(tmpinfo)
    fstem = sprintf('%s%d',stype,j);
    switch stype
      case 'BOLD'
        switch tmpinfo(j).pepolar
          case 0
            fstemlist = {[fstem '_for']};
          case 1
            fstemlist = {[fstem '_rev']};
          case 2
            fstemlist = {[fstem '_rev'],[fstem '_for']};
          case 3            
            fstemlist = {[fstem '_for'],[fstem '_rev']};
        end;
      case 'DTI'
        switch tmpinfo(j).pepolar
          case 0
            fstemlist = {fstem};
          case 1
            fstemlist = {[fstem '_rev']};
          case 2
            fstemlist = {[fstem '_rev'],fstem};
          case 3            
            fstemlist = {fstem,[fstem '_rev']};
        end;
      otherwise
        fstemlist = {fstem};
    end;
    for k=1:length(fstemlist)
      fstem = fstemlist{k};
      fname_in = sprintf('%s/%s%s%s',...
        ContainerPath,fstem,parms.infix,parms.in_ext);
      fname_out = sprintf('%s/%s%s%s',...
        parms.outdir,fstem,parms.infix,parms.out_ext);
      if ~exist(fname_in,'file')
        if parms.verbose
          fprintf('%s: WARNING: file %s not found\n',mfilename,fname_in);
        end;
        continue;
      end;
      if strcmp(fname_in,fname_out)
        if parms.verbose
          fprintf('%s: WARNING: input and output file names are identical: %s\n',...
            mfilename,fname_in);
        end;
        continue;
      end;
      if ~exist(fname_out,'file') || parms.forceflag      
        fs_mri_convert(fname_in,fname_out,...
          'out_orient',parms.out_orient,...
          'verbose',parms.verbose,...
          'forceflag',parms.forceflag);
      end;
    end;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

