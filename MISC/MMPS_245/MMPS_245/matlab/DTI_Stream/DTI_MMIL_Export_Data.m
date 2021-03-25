function [fname_data,fname_grad,fname_bval,max_bval,nb0]=...
  DTI_MMIL_Export_Data(ContainerPath,varargin)
%function [fname_data,fname_grad,fname_bval,max_bval,nb0]=...
%  DTI_MMIL_Export_Data(ContainerPath,[options])
%
% Purpose:
%   Concatenates raw data files (specified by snums)
%     and writes diffusion gradient table as text file (corrected for motion)
%
% Usage:
%  DTI_MMIL_Export_Data(ContainerPath,'key1', value1,...);
%
% Required Input Parameters:
%   ContainerPath: full path of directory containing processed diffusion data
%                  (mgh/mgz format)
%
% Optional Parameters specifying which data to load:
%   'snums': list of scan numbers to concatenate and analyze
%     if empty (or unspecified), use all DTI scans in container
%     {default = []}
%   'infix': if empty, will look for files like 'DTI1.mgz'
%     otherwise, input file will be sprintf('DTI%d_%s.mgz',snum,infix)
%     example infix = 'corr_resDTI'
%     {default = []}
%   'revflag': [0|1|2] specify whether to use non-rev or rev data
%     if revflag=0, use non-rev data
%     if revflag=1, use rev data
%       rev scans have names like 'DTI1_rev.mgz'
%     if revflag=2, use concatenated non-rev and rev data
%     {default = 0}
%   'min_ndirs': require at least this many diffusion directions to be valid
%     {default = 6}
%   'min_bval': minimum b value a scan must have to be included
%     {default = 1}
%   'flex_flag': [0|1] DTI_flex scans included in tensor fit
%     {default = 0}
%   'min_nb0': minimum number of b=0 image to be included
%     {default = 1}
%
% Optional Parameters specifying output
%  'outdir': output directory
%    {default = ContainerPath}
%  'outstem': output file stem (relative to outdir)
%    {default = 'data'}
%  'outtype': output file type
%    supported types: 'nii', 'mgh','mgz'
%    'nii' format can be used by FSL, AFNI, and ExploreDTI
%    'mgh' and 'mgz' formats can be used by FreeSurfer
%    { default = 'nii' }
%  'orient': output slice orientation
%    if empty or omitted, keep original orientation
%    e.g. 'LPS', 'RAS', etc.
%    For use with FSL, use 'LAS'
%    NOTE: no adjustments made to gradient table
%    { default = []}
%  'fname_data': full path of output concatenated data file
%    If supplied, will override outdir and outstem
%    {default = []}
%  'fname_grad': file name of gradient table text file
%    Relative to outdir unless full path is given
%    {default = 'gradtable.txt'};
%  'fname_bval': file name of b-value text file
%    Relative to dir unless full path is given
%    {default = 'bvals.txt'};
%  'scale_qmat_flag': [0|1] whether to scale gradient table by
%    sqrt(b-value)
%    { default = 1 }
%  'sort_b0vols_flag': [0|1] whether to sort b=0 images to beginning of volume
%    { default = 1 }
%  'nob0_qmat_flag': [0|1] do not include b=0 vectors in grad table
%    { default = 1 }
%  'swapxy_flag': [0|1] whether to swap x and y of qmat
%    { default = 1 }
%  'forceflag': [0|1] whether to run conversions even if output files exist
%    { default = 0 }
%
% Optional Output:
%   fname_data: full path of output data file
%   fname_grad: full path of output gradient table file
%   fname_bval: full path of output bvalue file
%   max_bval: maximum b value of exported data files
%
% see also: DTI_MMIL_Export_DTIstudio
%
% Created:  08/28/09 by Don Hagler
% Prev Mod: 02/04/13 by Don Hagler
% Last Mod: 02/19/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse input parameters
if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'snums',[],[],...
  'infix',[],[],...
  'revflag',0,[0,1,2],...
  'min_nb0',1,[],...
  'min_ndirs',6,[],...
  'min_bval',1,[],...
  'flex_flag',false,[false true],...
...
  'outdir',ContainerPath,[],...
  'outtype','nii',{'nii','mgh','mgz'},...
  'orient',[],[],...
  'outstem','data',[],...
  'fname_data',[],[],...
  'fname_grad','gradtable.txt',[],...
  'fname_bval','bvals.txt',[],...
  'scale_qmat_flag',true,[false true],...
  'sort_b0vols_flag',true,[false true],...
  'nob0_qmat_flag',true,[false true],...
  'swapxy_flag',true,[false true],...
  'forceflag',false,[false true],...
...
  'smf',10^-5,[10^-100,10^-1],...
  'orient_ref','LPS',[],...
});

fname_data = parms.fname_data;
fname_grad = parms.fname_grad;
fname_bval = parms.fname_bval;
max_bval = [];
nb0 = [];

if isempty(parms.fname_data)
  % construct output file name
  parms.fname_data = [parms.outdir '/' parms.outstem '.' parms.outtype];
else
  % add proper extension if missing from input fname_data
  [tmp_path,tmp_stem,tmp_ext] = fileparts(parms.fname_data);
  if ~strcmp(tmp_ext,['.' parms.outtype])
    parms.fname_data = [parms.fname_data '.' parms.outtype];
  end;  
end;

% make fname_grad a full path
if mmil_isrelative(parms.fname_grad)
  parms.fname_grad = [parms.outdir '/' parms.fname_grad];
end;
% make fname_bval a full path
if mmil_isrelative(parms.fname_bval)
  parms.fname_bval = [parms.outdir '/' parms.fname_bval];
end;

if nargout<=3 &...
   exist(parms.fname_data,'file') &...
   exist(parms.fname_grad,'file') & ~parms.forceflag
  return;
end;

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load concatenated volume and gradient directions
[vol,M,qmat,bvals,TEs]=DTI_MMIL_Load_Data(ContainerPath,...
  'snums',parms.snums,'infix',parms.infix,'revflag',parms.revflag,...
  'min_nb0',parms.min_nb0,'min_ndirs',parms.min_ndirs,...
  'min_bval',parms.min_bval,'flex_flag',parms.flex_flag);
if isempty(vol), return; end;

% set bval to 0 if length of qvec is zero
qlength = sqrt(sum(qmat.^2,2));
i_b0 = find(qlength<parms.smf);
nb0 = length(i_b0);
bvals(i_b0) = 0;

% sort b=0 images to beginning
if parms.sort_b0vols_flag
  [bvals,ind] = sort(bvals);
  qmat = qmat(ind,:);
  vol = vol(:,:,:,ind);
end;  

% scale qmat by sqrt(bvals/max(bvals))
max_bval = max(bvals);
if parms.scale_qmat_flag
  tmp_bvals = bvals/max_bval;
  tmp_bvals = sqrt(tmp_bvals * ones(1,3));
  qmat = tmp_bvals.*qmat;
end;

% exclude b=0 vectors from qmat
if parms.nob0_qmat_flag
  qmat = qmat(bvals~=0,:);
end;

% flip qmat components as needed to match orient_ref
if ~isempty(parms.orient)
  orient = parms.orient;
else
  orient = fs_read_orient([],M);
end;
[permvec,flipvec] = fs_compare_orient(orient,parms.orient_ref);
if any(flipvec==-1)
  qmat = qmat .* repmat(flipvec,size(qmat,1),1);
end;

% swap x and y of qmat
if parms.swapxy_flag
  tmp_qmat = qmat;
  tmp_qmat(:,1) = qmat(:,2);
  tmp_qmat(:,2) = qmat(:,1);
  qmat = tmp_qmat;
end;

% write bvals text file
fid = fopen(parms.fname_bval,'wt');
if fid==-1, error('failed to open file %s for writing',parms.fname_bval); end;
fprintf(fid,'%d\n',bvals);
fclose(fid);

% write gradtable text file
fid = fopen(parms.fname_grad,'wt');
if fid==-1, error('failed to open file %s for writing',parms.fname_grad); end;
for i=1:size(qmat,1)
  fprintf(fid,'%s\n',strtrim(sprintf('%0.6f ',qmat(i,:))));
end;
fclose(fid);

% save vol in correct format
if ~exist(parms.fname_data,'file') | parms.forceflag
  fprintf('%s: exporting data to %s...\n',mfilename,parms.fname_data);
  switch parms.outtype
    case 'nii'
      [tmp_path,tmp_stem] = fileparts(tempname);
      fname_tmp = [parms.outdir '/' tmp_stem '.mgh'];
      if ~isempty(parms.orient)
        [vol,M] = fs_reorient(vol,M,parms.orient);
      end;
      fs_save_mgh(vol,fname_tmp,M);
      fs_mri_convert(fname_tmp,parms.fname_data,...
        'forceflag',parms.forceflag);
      delete(fname_tmp);
    case {'mgh','mgz'}
      fs_save_mgh(vol,parms.fname_data,M);
    otherwise
      error('invalid outtype %s',parms.outtype);
  end;
end;

fname_data = parms.fname_data;
fname_grad = parms.fname_grad;
fname_bval = parms.fname_bval;
