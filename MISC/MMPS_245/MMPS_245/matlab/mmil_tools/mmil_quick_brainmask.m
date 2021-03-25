function vol_out = mmil_quick_brainmask(vol_in,varargin)
%function vol_out = mmil_quick_brainmask(vol_in,[options])
%
% Purpose: create a filled mask from a brain image volume
%     does two rounds of filling and removes edge slices
%   this works well as an approximate brain mask from T2-weighted images
%
% Usage: vol_out = mmil_quick_brainmask(vol_in,'key1', value1,...);
%
% Required Input:
%   vol_in: input mask volume (ctx format)
%     Can be [] if 'fname_in' is specified
%
% Optional Input:
%  'fname_in': full path of input mask volume file
%    if supplied, vol_in is ignored and fname_in is loaded instead
%    {default = []}
%  'fname_out': full path of output mask volume file
%    If supplied, vol_out is saved as fname_out
%    {default = []}
%  'log_flag': calculate log of input volume
%     {default = 1}
%  'thresh': relative threshold applied to input volume
%     cumulative probability of intensity (max is 1)
%     {default = 0.9}
%  'fill1_smooth1': smoothing sigma for first step of first round of filling
%    {default = 10}
%  'fill1_thresh1': threshold applied to mask after first smoothing step
%    {default = 0.95}
%  'fill1_smooth2': smoothing sigma (voxels) for second dilation step
%    {default = 20}
%  'fill1_thresh2': threshold applied to mask after second smoothing step
%    {default = 0.1}
%  'fill1_smooth3': smoothing sigma (voxels) for third dilation step
%    {default = 0}
%  'fill1_erode_flag': [0|1] allow mask erosion by not resetting original values
%     after smoothing and thresholding
%    {default = 1}
%  'clip_edges_flag': [0|1] set edge slices to 0 after first round of filling
%    {default = 1}
%  'clip_edges_width': number of edge slices to set to 0
%    {default = 2}
%  'fill2_smooth1': smoothing sigma for first step of second round of filling
%    {default = 20}
%  'fill2_thresh1': threshold applied to mask after first smoothing step
%    {default = 0.95}
%  'fill2_smooth2': smoothing sigma (voxels) for second dilation step
%    {default = 20}
%  'fill2_thresh2': threshold applied to mask after second smoothing step
%    {default = 0.1}
%  'fill2_smooth3': smoothing sigma (voxels) for third dilation step
%    {default = 20}
%  'fill2_smooth3': smoothing sigma (voxels) for third dilation step
%    {default = 20}
%  'fill2_thresh3': threshold applied to mask after third smoothing step
%    sets small values to zero, but does not binarize mask
%    (if smooth3=0, this is ignored)
%    {default = 0.05}
%  'fill2_erode_flag': [0|1] allow mask erosion by not resetting original values
%     after smoothing and thresholding
%    {default = 1}
%  'binary_flag': make binary mask (otherwise smooth edges)
%     {default = 0}
%  'forceflag': [0|1] whether to overwrite fname_mask if it already exists
%    {default = 0}
%
% See Also:
%   mmil_dilate_mask: to create mask from mask volume
%   mmil_masksurf: to create surface mesh from mask
%   mmil_dct_brainmask: to create mask from image volume using dct morph to atlas
%
% Created:  05/07/12 by Don Hagler
% Last Mod: 12/29/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_in',[],[],...
  'fname_out',[],[],...
... % mask creation
  'log_flag',true,[false true],...
  'thresh',0.9,[],...
... % mask fill 1
  'fill1_smooth1',10,[],...
  'fill1_thresh1',0.95,[],...
  'fill1_smooth2',20,[],...
  'fill1_thresh2',0.1,[],...
  'fill1_smooth3',0,[],...
  'fill1_erode_flag',true,[false true],...
...
  'clip_edges_flag',true,[false true],...
  'clip_edges_width',2,[1 10],...
... % mask fill 2
  'fill2_smooth1',20,[],...
  'fill2_thresh1',0.95,[],...
  'fill2_smooth2',20,[],...
  'fill2_thresh2',0.1,[],...
  'fill2_smooth3',20,[],...
  'fill2_thresh3',0.05,[],...
  'fill2_erode_flag',true,[false true],...
... % final binarization
  'binary_flag',false,[false true],...
... % other
  'forceflag',false,[false true],...
});

if ~isempty(parms.fname_out) &&...
   exist(parms.fname_out,'file') && ~parms.forceflag
  if nargout>0
    vol_out = ctx_load_mgh(parms.fname_out);
  end;
  return;
end;

if ~isempty(parms.fname_in)
  if ~exist(parms.fname_in,'file')
    error('file %s not found',parms.fname_in);
  end;
  [vol_in,mr_parms] = ctx_load_mgh(parms.fname_in);
else
  mr_parms = [];
end;

vol_out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

volmask = vol_in;

% log transform intensities
if parms.log_flag
  volmask.imgs = log(max(1,volmask.imgs));
end;

% determine threshold based on cumulative probability
[hc,hv] = hist(mmil_colvec(volmask.imgs),1000);
chc = cumsum(hc)/sum(hc);
ind_thresh = max(1,min(find(chc>parms.thresh))-1);
thresh = hv(ind_thresh);

volmask.imgs = 1.0*(volmask.imgs>thresh);

% erode then dilate mask
if parms.fill1_smooth1>0
  volmask = mmil_dilate_mask(volmask,...
    'smooth1',parms.fill1_smooth1,'thresh1',parms.fill1_thresh1,...
    'smooth2',parms.fill1_smooth2,'thresh2',parms.fill1_thresh2,...
    'smooth3',parms.fill1_smooth3,...
    'erode_flag',parms.fill1_erode_flag);
end;

if parms.clip_edges_flag
  % set edge slices to 0
  volsz = size(volmask.imgs);
  if volsz(1) > parms.clip_edges_width*2
    volmask.imgs([1:parms.clip_edges_width,end-(parms.clip_edges_width-1):end],:,:) = 0;
  end;
  if volsz(2) > parms.clip_edges_width*2
    volmask.imgs(:,[1:parms.clip_edges_width,end-(parms.clip_edges_width-1):end],:) = 0;
  end;
  if volsz(3) > parms.clip_edges_width*2
    volmask.imgs(:,:,[1:parms.clip_edges_width,end-(parms.clip_edges_width-1):end]) = 0;
  end;
end;

% erode then dilate mask
if parms.fill2_smooth1>0
  volmask = mmil_dilate_mask(volmask,...
    'smooth1',parms.fill2_smooth1,'thresh1',parms.fill2_thresh1,...
    'smooth2',parms.fill2_smooth2,'thresh2',parms.fill2_thresh2,...
    'smooth3',parms.fill2_smooth3,'thresh3',parms.fill2_thresh3,...
    'erode_flag',parms.fill2_erode_flag);
end;

% binarize mask
if parms.binary_flag
  volmask.imgs = 1.0*(volmask.imgs>=eps);
end;

if ~isempty(parms.fname_out)
  ctx_save_mgh(volmask,parms.fname_out,mr_parms);
end;

if nargout>0
  vol_out = volmask;
end;

