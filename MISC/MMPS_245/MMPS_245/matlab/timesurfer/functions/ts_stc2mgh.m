function ts_stc2mgh(stcfile,mghfile,subj,hemi,sparsesmooth,postsmooth,t0,t1,...
  subjdir,mbmask_flag,forceflag)
% function ts_stc2mgh(stcfile,mghfile,subj,hemi,[sparsesmooth],[postsmooth],...
%   [t0],[t1],[subjdir],[mbmask_flag],[forceflag])
% 
% converts a stc (source time course) file to mgh format (readable by tksurfer)
% optionally applies smoothing on surface
%
%  Required parameters:
%    stcfile: full pathname of input stc file
%    mghfile: full pathname of output mgh file
%    subj:    subject name
%    hemi:    cortical hemisphere
%
%  Optional parameters:
%    sparsesmooth: number of sparse smoothing steps
%      {default = 0}
%      note: sparse smoothing is a fast way to fill
%            in gaps between sparse vertices
%    postsmooth: number of normal smoothing steps
%      {default = 0}
%      note: postsmoothing is additional nearest-neighbor average
%            smoothing applied after sparse smoothing
%    t0: first time sample to extract
%      {default = 1}
%    t1: last time sample to extract
%      {default = last}
%   subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%   mbmask_flag: [0|1] toggle mask out thalamus, corpus callosum
%    {default: 0}
%   forceflag : [0|1] toggle overwrite of existing mghfile
%    {default: 1}
%
% created:       <01/01/07 Don Hagler
% last modified:  10/15/10 Don Hagler
%

if (~mmil_check_nargs(nargin,4)) return; end;

if ~exist('sparsesmooth','var') | isempty(sparsesmooth), sparsesmooth=0; end;
if ~exist('postsmooth','var') | isempty(postsmooth), postsmooth=0; end;
if ~exist('t0','var') || isempty(t0), t0=1; end;
if ~exist('t1','var') || isempty(t1), t1=0; end;
if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;
if ~exist('mbmask_flag','var') | isempty(mbmask_flag), mbmask_flag=0; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag=1; end;

if(~exist(mghfile,'file') || forceflag)
  fprintf('%s: reading stcfile %s...\n',mfilename,stcfile);
  [starttime,sample_period,vertices,sol]=ts_read_stc(stcfile);

  fprintf('%s: start time = %0.2f ms\n',mfilename,starttime);
  fprintf('%s: sample_period = %0.2f ms\n',mfilename,sample_period);

  nverts = length(vertices);
  [ndips,tpoints] = size(sol);

  if t0<1, t0 = 1; end;
  if t1<1 | t1>tpoints, t1 = tpoints; end;
  if t0>t1, t0=t1; end;

  if (nverts ~= ndips)
    fprintf('%s: error: num vertices (%d) does not match num dips (%d)!\n',...
      mfilename,nverts,ndips);
    return;
  end;

  fprintf('%s: num vertices = %d\n',mfilename,nverts);
  fprintf('%s: num time points = %d\n',mfilename,tpoints);

  if sparsesmooth | postsmooth
    fprintf('%s: reading %s surface for subject %s...\n',...
      mfilename,hemi,subj);
    surf = fs_load_subj(subj,hemi,[],[],subjdir);
  else
    fprintf('%s: reading number of %s verts for subject %s...\n',...
      mfilename,hemi,subj);
    surf = fs_load_subj(subj,hemi,[],1,subjdir);
  end;

  tpoints = t1-t0+1;
  surfstats = zeros(surf.nverts,tpoints);
  surfstats(vertices+1,:) = sol(:,t0:t1);

  % mask midbrain
  if mbmask_flag
    surfstats = fs_mask_surfstats_with_aparc(surfstats,subj,hemi,subjdir);
  end;

  if sparsesmooth>0 | postsmooth>0
    fprintf('%s: smoothing timepoints %d-%d of %s...\n',mfilename,t0,t1,stcfile);
    for t=1:tpoints
      vals = surfstats(:,t);
      if sparsesmooth>0
        vals = fs_smooth_sparse(surf,vals,sparsesmooth);
      end;
      if postsmooth>0
        vals = fs_smooth(surf,vals,postsmooth);
      end;
      surfstats(:,t) = vals;
    end;
    % mask midbrain
    if mbmask_flag
      surfstats = fs_mask_surfstats_with_aparc(surfstats,subj,hemi,subjdir);
    end;
  end

  fprintf('%s: writing mghfile %s...\n',mfilename,mghfile);
  vol = reshape(surfstats,[surf.nverts,1,1,tpoints]);
  fs_save_mgh(vol,mghfile);

  if ~exist(mghfile,'file')
    warning('failed to write mghfile');
  else
    fprintf('%s: finished writing mghfile\n',mfilename);
  end;
end;
