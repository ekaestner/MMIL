function dti_CSD_tracto(vol,M,qmat,bvals,varargin)
%function dti_CSD_tracto(vol,M,qmat,bvals,[options])
%
% Purpose:
%  Create fiber streamlines using CSD tractography
%   from Alexander Leeman's ExploreDTI
%
% Usage:
%  dti_CSD_tracto(vol,M,qmat,bvals,'key1', value1,...);
%
% Required Parameters:
%   vol: 4D volume containing multiple diffusion weighted volumes
%   M: 4x4 vox2ras matrix for vol
%   qmat: matrix of diffusion direction vectors
%   bvals: vector of b values
%     one for all, or one for each diffusion direction
%
% Optional Parameters:
%   'outdir': full or relative path of output directory
%     {default = 'CSD_Tracto'}
%   'outstem': output file stem
%     {default = []}
%   'orient': reorient data to specified slice orientation
%     (e.g. 'LAS', 'LPI', etc.)
%     {default = []}
%   'forceflag': [0|1] whether to overwrite existing output
%     {default = 0}
%
% Optional Parameters for CSD tractography:
%   'see_point_sampling': seed point sampling in x, y, and z direction (in mm)
%     {default = [2 2 2]}
%   'step_size': step size (in mm)
%     {default = 1}
%   'FOD_thresh': FOD threshold
%     {default = 0.1}
%   'angle threshold': angle deviation threshold
%     {default = 45}
%   'fiber_length_range': min and max fiber length (in mm)
%     {default = [50 500]}
%   'max_order': highest order spherical harmonic (even integers)
%     {default = 4}
%   'thresh_FA': minimal FA to determine response function
%     {default = 0.8}
%
% Created:  05/11/12 by Don Hagler
% Last Mod: 10/27/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,4), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir','CSD_Tracto',[],...
  'outstem',[],[],...
  'orient',[],[],...
  'forceflag',false,[false true],...
...
  'seed_point_sampling',[2 2 2],[],...
  'step_size',1,[],...
  'FOD_thresh',0.1,[],...
  'angle_thresh',45,[],...
  'fiber_length_range',[50 500],[],...
  'max_order',4,[2:2:100],...
  'thresh_FA',0.8,[0 1],...
...
  'smf',10^-5,[10^-100,10^-1],...
  'orient_ref','LPS',[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if mmil_isrelative(parms.outdir)
  parms.outdir = [pwd '/' parms.outdir];
end;
if ~isempty(parms.outstem)
  parms.outstem = [parms.outstem '_'];
end;
parms.outstem = sprintf('%s/%s',parms.outdir,parms.outstem);

[nx,ny,nz,nf] = size(vol);
volsz = [nx,ny,nz];

if length(bvals)==1
  bvals = bvals*ones(nf,1);
end;
if length(bvals)~=nf
  error('number of bvals must match frames in vol');
end;

if size(qmat,1)~=nf
  error('number of directions in qmat must match frames in vol');
end;
if size(qmat,2)~=3
  error('qmat must have 3 columns');
end;

if ~isempty(parms.orient)
  orient = parms.orient;
else
  orient = fs_read_orient([],M);
end;

orient_swapxy = orient;
orient_swapxy(2) = orient(1);
orient_swapxy(1) = orient(2);

M_orig = M;

[permvec,flipvec] = fs_compare_orient(orient,parms.orient_ref);
flipflags = flipvec<0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data, qmat, and bvals

mmil_mkdir(parms.outdir);

% find b=0 images
qlength = max(sqrt(sum(qmat.^2,2)),eps);
i_b0 = find(qlength<parms.smf);
nb0 = length(i_b0);

% set bval to 0 if length of qvec is zero
bvals(i_b0) = 0;

% sort b=0 images to beginning
[bvals,ind_sort] = sort(bvals);

% find max b value
max_bval = max(bvals);

% write bvals text file
fname_bvals    = sprintf('%sbvals.txt',parms.outstem);
if ~exist(fname_bvals,'file') || parms.forceflag
  fid = fopen(fname_bvals,'wt');
  if fid==-1, error('failed to open file %s for writing',fname_bvals); end;
  fprintf(fid,'%d\n',bvals);
  fclose(fid);
end;

% write bvecs text file
fname_bvecs    = sprintf('%sbvecs.txt',parms.outstem);
if ~exist(fname_bvecs,'file') || parms.forceflag
  % ensure that qmat contains unit vectors
  qmat = qmat ./ repmat(qlength,1,3);

  % sort b=0 images to beginning
  qmat = qmat(ind_sort,:);

  % scale qmat by sqrt(bvals/max(bvals))
  qmat = qmat .* repmat(sqrt(bvals/max_bval),1,3);

  % exclude b=0 vectors from qmat
  qmat = qmat(bvals~=0,:);

  % flip qmat components as needed to match orient_ref
  if any(flipvec==-1)
    qmat = qmat .* repmat(flipvec,size(qmat,1),1);
  end;

  % write qmat to text file
  fid = fopen(fname_bvecs,'wt');
  if fid==-1, error('failed to open file %s for writing',fname_bvecs); end;
  for i=1:size(qmat,1)
    fprintf(fid,'%s\n',strtrim(sprintf('%0.6f ',qmat(i,:))));
  end;
  fclose(fid);
end;

fname_data_mgh = sprintf('%sdata.mgz',parms.outstem);
fname_data_swapxy = sprintf('%sdata_swapxy.mgz',parms.outstem);
if ~exist(fname_data_mgh,'file') ||...
   ~exist(fname_data_swapxy,'file') || parms.forceflag
  % sort b=0 images to beginning
  vol = vol(:,:,:,ind_sort);

  % reorient data
  if ~isempty(parms.orient)
    [vol,M] = fs_reorient(vol,M,parms.orient);
  end;

  % save to mgh file
  fs_save_mgh(vol,fname_data_mgh,M);

  % swap x and y dimensions (because csd_wbt reads volume like that)
  vol = fs_reorient(vol,M,orient_swapxy);

  % save to mgh file
  fs_save_mgh(vol,fname_data_swapxy,M);
elseif ~isempty(parms.orient)
  M = fs_read_header(fname_data_mgh);
end;

% convert data to nii
fname_data_nii = sprintf('%sdata_swapxy.nii',parms.outstem);
fs_mri_convert(fname_data_swapxy,fname_data_nii,'forceflag',parms.forceflag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run CSD tractography
fname_CSD = sprintf('%sdata_swapxy_Tracts_CSD.mat',parms.outstem);
if ~exist(fname_CSD,'file') || parms.forceflag
  % create text parameter file
  fname_parms = sprintf('%sparameters.txt',parms.outstem);
  fid = fopen(fname_parms,'wt');
  if fid==-1, error('failed to open file %s for writing',fname_parms); end;
  fprintf(fid,'%s\n',fname_data_nii);
  fprintf(fid,'%s\n',fname_bvecs);
  fprintf(fid,'%s\n',parms.outdir);
  fprintf(fid,'%d\n',nb0);
  fprintf(fid,'%d\n',max_bval);
  fprintf(fid,'%s\n',strtrim(sprintf('%d ',parms.seed_point_sampling)));
  fprintf(fid,'%0.2f\n',parms.step_size);
  fprintf(fid,'%0.2f\n',parms.FOD_thresh);
  fprintf(fid,'%0.2f\n',parms.angle_thresh);
  fprintf(fid,'%s\n',strtrim(sprintf('%d ',parms.fiber_length_range)));
  fprintf(fid,'%d\n',parms.max_order);
  fprintf(fid,'%0.2f\n',parms.thresh_FA);
  fclose(fid);
  % run csd_wbt
  fprintf('%s: running CSD tractography...\n',mfilename);
  csd_wbt(fname_parms);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save tract mask as mgz
fname_mask = sprintf('%stract_mask.mgh',parms.outstem);
if ~exist(fname_mask,'file') || parms.forceflag
  load(fname_CSD);
  if flipflags(3)
    TractMask = TractMask(:,:,end:-1:1);
  end;
  fs_save_mgh(TractMask,fname_mask,M);
end;

% convert fibers to DTI Studio format
fname_grp = sprintf('%sfiber_path.grp',parms.outstem);
dti_CSD_to_DTIStudio_fiber(fname_CSD,fname_grp,...
  'M',M,'forceflag',parms.forceflag);

fname_dat = sprintf('%sfiber_count.dat',parms.outstem);
if ~exist(fname_dat,'file') || parms.forceflag
  dti_fiber_path_to_mask(fname_grp,fname_dat,1);
end;

fname_mgh = sprintf('%sfiber_count.mgh',parms.outstem);
if ~exist(fname_mgh,'file') || parms.forceflag
  dti_dat2mgh(fname_dat,volsz,M,[],flipflags(3));
end;

