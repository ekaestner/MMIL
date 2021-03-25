function vals = fs_tal2surf(points,varargin)
%function vals = fs_tal2surf(points,[options])
%
% Purpose: sample Talairach coordinates to a subject's cortical surfaces
%
% Required Input:
%   points: npoints x 3 matrix of Talairach coordinates
%
% Optional Input:
%   'subj': FreeSurfer subject
%     {default = 'fsaverage'}
%   'subjdir': FreeSurfer subject root directory
%     {default = $FREESURFER_HOME/subjects}
%   'MNIflag': [0|1] whether input points are MNI Talaiarach
%       versus true Talaiarach coordinates
%     may be a vector, with one value for each point
%     {default = 1}
%   'surfname': FreeSurfer surface name
%     {default = 'white'}
%   'vol_sm': sigma for smoothing in volume
%     {default = 0}
%   'surf_sm': number of iterations for smoothing on surface
%     {default = 0}
%   'thresh': threshold applied after smoothing on surface
%     maximum value in output vals is 1
%     {default = 0}
%   'binarize_flag': create binarized mask instead of continuous values
%     {default = 0}
%   'hemilist': cell array of cortical hemispheres to search
%     {default = {'lh','rh'}}
%   'verbose': [0|1] display status messages
%     {default = 1}
%
% Output:
%   vals: cell array with two elements, one for each hemisphere (left,right)
%         inside each cell is matrix of values with size [nverts,npoints]
%      if hemilist has a single element, vals will be returned as a matrix
%         instead of a cell array
%
% Created:  11/02/11 by Don Hagler
% Last Mod: 01/05/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'subj','fsaverage',[],...
  'subjdir',[],[],...
  'MNIflag',true,[false true],...
  'surfname','white',[],...
  'vol_sm',0,[0,1000],...
  'surf_sm',0,[0,1000],...
  'thresh',0,[0,1],...
  'binarize_flag',false,[false true],...
  'hemilist',{'lh','rh'},{'lh','rh'},...
  'verbose',true,[false true],...
...
  'tmpdir',[pwd '/tmp'],[],...
});
vals = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(points,2)~=3, error('size of points must be npoints x 3'); end;
npoints = size(points,1);

if length(parms.MNIflag)==1 && npoints>1
  parms.MNIflag = boolean(ones(npoints,1)*parms.MNIflag);
end;
if length(parms.MNIflag)~=npoints
  error('MNIflag must have 1 or %d values',npoints);
end;

if ~iscell(parms.hemilist)
  parms.hemilist = {parms.hemilist};
end;

% check subjects dir
if isempty(parms.subjdir)
  fshome = getenv('FREESURFER_HOME');
  if isempty(fshome)
    error('FREESURFER_HOME not defined');
  end;
  parms.subjdir = [fshome '/subjects'];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% optionally transform from true Talairach to MNI Talairach using fs_tal2mni
if any(~parms.MNIflag)
  % transform points if MNIflag = 1
  ind_tal = find(~parms.MNIflag);
  points(ind_tal,:) = fs_tal2mni(points(ind_tal,:));
end;

% transform from MNI to RAS
if ~strcmp(parms.subj,'fsaverage') %  not necessary for fsverage
  % load transform using fs_read_talxfm
  fname = sprintf('%s/%s/mri/transforms/talairach.xfm',...
    parms.subjdir,parms.subj);
  if ~exist(fname,'file'), error('file %s not found',fname); end;
  T = fs_read_talxfm(fname);
  points = fs_mni2ras(points,T);
end;

% load surfaces
clear surfs;
for h=1:length(parms.hemilist)
  hemi = parms.hemilist{h};
  fname = sprintf('%s/%s/surf/%s.%s',...
    parms.subjdir,parms.subj,hemi,parms.surfname);
  surf = fs_read_surf(fname);
  if parms.surf_sm
    surf = fs_find_neighbors(surf,parms.verbose);
  end;
  surfs(h) = surf;
end;

if parms.vol_sm
  mmil_mkdir(parms.tmpdir);
  % read info from orig.mgz
  fname_orig = sprintf('%s/%s/mri/orig.mgz',...
    parms.subjdir,parms.subj);
  [M,volsz] = fs_read_header(fname_orig);
  % voxel indices from RAS coordinates
  XYZ = (inv(M)*[points';ones(1,size(points,1))])';
  XYZ = round(XYZ(:,1:3));
  for i=1:npoints
    % create volume for each point
    vol = zeros(volsz);
    vol(XYZ(i,1),XYZ(i,2),XYZ(i,3)) = 1;
    % smooth in volume
    vol_sm = mmil_smooth3d(vol,parms.vol_sm,parms.vol_sm,parms.vol_sm);
    % normalize to keep max value 1
    vol_sm = vol_sm * max(vol(:))/max(vol_sm(:));
    % save temporary file
    fname_tmp = sprintf('%s/tmp_vol_point%d.mgh',parms.tmpdir,i);
    fs_save_mgh(vol_sm,fname_tmp,M);
    % paint to surface
    fnames_out = fs_paint(parms.subj,fname_tmp,'subjdir',parms.subjdir);
    % load output of paint, put in vals
    for h=1:length(parms.hemilist)
      vals{h}(:,i) = fs_load_mgh(fnames_out{h});
    end;
  end;
  [s,r] = unix(sprintf('rm -r %s',parms.tmpdir));
  if s, error('failed to delete tmpdir %s:\n%s',parms.tmpdir,r); end;
else
  % find nearest vertex to RAS coords
  v = {};
  vdist = {};
  for h=1:length(parms.hemilist)
    [v{h},vdist{h}] = fs_nearest_verts(points,surfs(h));
    vals{h} = zeros(surfs(h).nverts,npoints);
  end;
  if length(parms.hemilist)>1
    % test whether lh or rh vertices are closer to each point
    [min_dist,min_hemi] = min([vdist{1};vdist{2}]);
  else
    min_hemi = ones(1,npoints);
  end;
  % set values for vertices to 1
  for i=1:npoints
    h = min_hemi(i);
    v_tmp = v{h}(i);
    vals{h}(v_tmp,i) = 1;
  end;
end;

% smooth on surface
if parms.surf_sm
  for h=1:length(parms.hemilist)
    tmp1 = vals{h};
    if nnz(tmp1)
      for i=1:npoints
        tmp2 = tmp1(:,i);
        if nnz(tmp2)
          % smooth values
          tmp3 = fs_smooth(surfs(h),tmp2,parms.surf_sm);
          % normalize to keep max value 1
          tmp3 = tmp3 * max(tmp2)/max(tmp3);
          tmp1(:,i) = tmp3;
        end;
      end;
      vals{h} = tmp1;
    end;
  end;
end;

% apply threshold
if parms.thresh
  for h=1:length(parms.hemilist)
    tmp = vals{h};
    tmp(tmp<parms.thresh) = 0;
    vals{h} = tmp;
  end;
end;

% binarize values
if parms.binarize_flag
  for h=1:length(parms.hemilist)
    tmp = vals{h};
    tmp(tmp>0) = 1;
    vals{h} = tmp;
  end;
end;

if length(parms.hemilist)==1
  vals = vals{1};
end;

return;

