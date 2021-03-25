function mmil_multiframe_subplots(fname,varargin)
% function mmil_multiframe_subplots(fname,[options])
%
% Purpose:
%   Creates output 'TIFF' figures of input images
%
% Required Parameters:
%   fname: file name of input volume file (mgz or mgh file)
%
% Optional Paramters: 
%  'outdir': output directory
%     {default = pwd}
%  'outstem': output file stem
%     {default = 'multiframe_subplot'}
%  'frames': vector of frame indices
%     use all if empty
%     {default = []}
%  'border_flag': [0|1] allow borders around each subplot
%     {default = 0}
%  'fig_size': figure size in inches
%     if not specified, will use default
%     {default = []}
%  'tif_dpi': resolution of tif files (dots per inch)
%     {default = 300}
%  'planestring' : plane of the output image
%     'HOR', 'SAG', or 'COR'
%     {default = 'HOR'}
%  'slice_fract': fractional slice position
%     {default = 0.5}
%  'frame_fracts': fractional frame index vector
%     if empty, use all frames
%     {default = []}
%  'clim': vector of lower and upper bounds for imagesc
%     if empty, will auto-scale
%     {default = []}
%  'plim_flag': [0|1] calculate clim from percentile values in plim
%     {default = 1}
%  'plim': percentile values (0-100) for intensity scaling range
%     {default = [1,99]}
%  'logtrans_flag': [0|1] calculate percentiles on log-transformed values
%     {default = 0}
%  'groups': vector of grouping variable
%     e.g. bvals for dMRI
%     length must be equal to number of frames
%     if supplied, will scale intensities differently for each group using plim
%     {default = []}
%  'forceflag': [0|1] whether to overwrite existing output file
%    {default: 0}
%    
% Created:  04/11/13 by Don Hagler
% Last Mod: 03/10/16 by Don Hagler
%

%% todo: option for cropping
%% todo: option for nrows/ncols

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'outdir',pwd,[],...
  'outstem','multiframe_subplot',[],...
  'frames',[],[],...
  'border_flag',false,[false,true],...
  'fig_size',[],[],...
  'tif_dpi',300,[10,10000],...
  'planestring','HOR',{'HOR','SAG','COR'},...
  'slice_fract',0.5,[0,1],...
  'frame_fracts',[],[0,1],...
  'clim',[],[],...
  'plim_flag',true,[false true],...
  'plim',[1,99],[],...
  'logtrans_flag',false,[false true],...
  'groups',[],[],...
  'forceflag',false,[false true],...
  ...
  'keepsingle',true,[false true],...
  'orient','PRS',[],...
});

% set the output filename
parms.fname_out = sprintf('%s/%s.tif',parms.outdir,parms.outstem);

if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

% check whether fname exists
if ~exist(fname,'file')
  error('file %s not found',fname);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmil_mkdir(parms.outdir);

figure;
set(gcf,'Visible','off');

% load fname
[vol,M] = fs_load_mgh(fname,[],parms.frames,[],parms.keepsingle);
nframes = size(vol,4);
%% todo: use mmil_load_mgh_info and wait to load data, one frame at a time?

if isempty(parms.frames)
  if ~isempty(parms.frame_fracts)
    parms.frames = unique(max(1,round(parms.frame_fracts*nframes)));
  else
    parms.frames = [1:nframes];
  end;
end;
nframes_plot = length(parms.frames);

% check groups
if ~isempty(parms.groups)
  if length(parms.groups) ~= nframes
    error('length of groups vector (%d) does not match nframes (%d)',...
      length(parms.groups),nframes);
  end;
  [uniq_groups,tmp,ind_groups] = unique(parms.groups);
  ngroups = length(uniq_groups);
else
  ngroups = 1;
  uniq_groups = 1;
  ind_groups = ones(nframes,1);
end;

% set intensity scaling limits
if isempty(parms.clim) && parms.plim_flag
  parms.clim = zeros(ngroups,2);
  for i=1:ngroups
    j = find(ind_groups==i);
    % calculate mean across frames within group (saves memory, esp. for prctile)
    vals = mmil_colvec(mean(vol(:,:,:,j),4));
    % apply log transform
    if parms.logtrans_flag
      vals = log(max(1,vals));
    end;
    % calulate intensity limits from percentile values
    clim = prctile(vals,parms.plim);
    % transform back to original scale
    if parms.logtrans_flag
      clim = exp(clim);
    end;
    if any(isnan(clim)) || length(unique(clim))==1
      clim = [0 1]; % avoid error in imagesc
    end;
    parms.clim(i,:) = clim;
  end;
end;

% determine number of rows and columns
nrows = ceil(sqrt(nframes_plot));
ncols = ceil(nframes_plot/nrows);

k = 1;
for f=1:nframes_plot
  frame = parms.frames(f);
  % select frame
  %% todo: load one frame from file at time to save memory?
  tvol = vol(:,:,:,frame);
  % reorient volume to standard orientation
  tvol = fs_reorient(tvol,M,parms.orient);
  % extract single slice
  img = extract_slice(tvol,parms.planestring,parms.slice_fract);
  if parms.border_flag
    subplot(nrows,ncols,k); 
  else
    mmil_subplot_tight(nrows,ncols,k); 
  end;
  k = k + 1;
  if ~isempty(parms.clim)
    clim = parms.clim(ind_groups(f),:);
    imagesc(img,clim);
  else
    imagesc(img);
  end;
  colormap(gray); axis off;
end;

save_plot(parms);
close(gcf);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function img = extract_slice(vol,plane,slice_fract)
  dim = size(vol);
  s = choose_slice(dim,plane,slice_fract);
  switch plane
    case 'COR'
      img = rot90(squeeze(vol(s,:,:)));
    case 'SAG'
      img = rot90(squeeze(vol(:,s,:)));
    case 'HOR'
      img = squeeze(vol(:,:,s));
  end;
return;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function slice_num = choose_slice(dim,plane,slice_fract)
  switch plane
    case 'COR'
      d = 1;
    case 'SAG'
      d = 2;
    case 'HOR'
      d = 3;
  end;
  slice_num = round(dim(d)*slice_fract);
  slice_num = min(max(slice_num,1),dim(d));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_plot(parms)
  if ~isempty(parms.fig_size)
    if length(parms.fig_size)==1
      parms.fig_size = [parms.fig_size parms.fig_size];
    end;
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [parms.fig_size(1) parms.fig_size(2)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 parms.fig_size(1) parms.fig_size(2)]);
  end;
  print(gcf,'-dtiff',parms.fname_out,sprintf('-r %d',parms.tif_dpi));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

