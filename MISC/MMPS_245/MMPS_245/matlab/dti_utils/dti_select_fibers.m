function dti_select_fibers(fname_fiber,fname_ROI,varargin)
%function dti_select_fibers(fname_fiber,fname_ROI,[options])
%
% Purpose: select fiber streamlines that pass through specified ROIs
%
% Usage:
%  dti_select_fibers(fname_fiber,fname_ROI,'key1', value1,...)
%
% Required Input: 
%   fname_fiber: full or relative path of DTI Studio format fiber file
%   fname_ROI: name of mgh/mgz/mat format file containing mask for
%       region of interest specifying locations through which fibers must pass
%     may be cell array of multiple ROIs
%     files with ".mat" extension are assumed to be sparse
%       see mmil_save_sparse and mmil_load_sparse
%     {default = []}
%
% Optional Parameters:
%  'fiber': fiber struct; if supplied, will skip reading fname_fiber
%     {default = []}
%  'fname_out': output file name
%     {default = 'fibers.grp'}
%  'fname_NOT': name of mgh/mgz/mat format file containing binary mask for
%     exclusionary region of interest
%     May be cell array of multiple ROIs
%     {default = []}
%  'cutflag': [0|1] trim streamlines as they exit terminal ROIs
%     {default = 0}
%  'revsliceflag': [0|1] reverse slice order of input ROI volumes
%     {default = 1}
%  'thresh': threshold applied to ROI volumes
%     {default = 0.01}
%  'verbose': [0|1] display status messages
%     {default = 0}
%  'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  05/18/12 by Don Hagler
% Last Mod: 09/15/15 by Don Hagler
%

%% todo: min_fiberlen, fname_FA, FA_thresh
%% todo: make sure fibers pass through each ROI in correct order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_fiber',fname_fiber,[],...
  'fname_ROI',fname_ROI,[],...
...
  'fiber',[],[],...
  'fname_out','fibers.grp',[],...
  'fname_NOT',[],[],...
  'cutflag',false,[false true],...
  'revsliceflag',true,[false true],...
  'thresh',0.01,[],...
  'verbose',false,[false true],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

if exist(parms.fname_out,'file') && ~parms.forceflag
  return;
end;

if ~exist(parms.fname_fiber,'file') && isempty(parms.fiber)
  error('file %s not found',parms.fname_fiber);
end;

if ~isempty(parms.fname_ROI)
  if ~iscell(parms.fname_ROI), parms.fname_ROI = {parms.fname_ROI}; end;
  parms.nroi = length(parms.fname_ROI);
  for r=1:parms.nroi
    if ~exist(parms.fname_ROI{r},'file')
      error('file %s not found',parms.fname_ROI{r});
    end;
  end;
else
  parms.nroi = 0;
end;

if ~isempty(parms.fname_NOT)
  if ~iscell(parms.fname_NOT), parms.fname_NOT = {parms.fname_NOT}; end;
  parms.nnot = length(parms.fname_NOT);
  for r=1:parms.nnot
    if ~exist(parms.fname_NOT{r},'file')
      error('file %s not found',parms.fname_NOT{r});
    end;
  end;
else
  parms.nnot = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load fiber
if isempty(parms.fiber)
  if parms.verbose
    fprintf('%s: loading fiber file %s...\n',mfilename,parms.fname_fiber);
  end;
  fiber = dti_read_DTIStudio_fiber(parms.fname_fiber);
else
  fiber = parms.fiber;
end;

% select fibers for each ROI
for r=1:parms.nroi
  if parms.verbose
    fprintf('%s: loading ROI file %s...\n',mfilename,parms.fname_ROI{r});
  end;
  parms.vol_ROI{r} = load_ROI(parms.fname_ROI{r},parms.revsliceflag);
  if parms.verbose
    fprintf('%s: selecting fibers...\n',mfilename);
  end;
  fiber = select_fibers(parms,fiber,parms.vol_ROI{r});
  if isempty(fiber)
    fprintf('%s: WARNING: no fibers selected\n',mfilename);
    return;
  end;
end;

% cut ends of fibers
if parms.nroi && parms.cutflag
  if parms.verbose
    fprintf('%s: cutting fibers...\n',mfilename);
  end;
  fiber = cut_fibers(parms,fiber);
end;

% exclude fibers for each NOT
for r=1:parms.nnot
  if parms.verbose
    fprintf('%s: loading NOT file %s...\n',mfilename,parms.fname_NOT{r});
  end;
  vol_NOT = load_ROI(parms.fname_NOT{r},parms.revsliceflag);
  if parms.verbose
    fprintf('%s: excluding fibers...\n',mfilename);
  end;
  fiber = exclude_fibers(parms,fiber,vol_NOT);
  if isempty(fiber)
    fprintf('%s: WARNING: all fibers excluded\n',mfilename);
    return;
  end;
end;

if parms.verbose
  fprintf('%s: writing fiber file %s...\n',mfilename,parms.fname_out);
end;
dti_write_DTIStudio_fiber(fiber,parms.fname_out);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fiber = select_fibers(parms,fiber,vol_ROI)
  fnums = find_fibers(parms,fiber,vol_ROI);
  fiber = update_fibers(fiber,fnums);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fiber = exclude_fibers(parms,fiber,vol_NOT)
  excl_fnums = find_fibers(parms,fiber,vol_NOT);
  all_fnums = 1:fiber.header.nFiberNr;
  fnums = setdiff(all_fnums,excl_fnums);
  fiber = update_fibers(fiber,fnums);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fiber = update_fibers(fiber,fnums)
  if isempty(fnums)
    fiber = [];
    return;
  end;
  fiber.chain = fiber.chain(fnums);
  fiber_lengths = [fiber.chain.nLength];
  fiber.header.nFiberNr      = length(fiber_lengths);
  fiber.header.nFiberLenMax  = max(fiber_lengths);
  fiber.header.fFiberLenMean = mean(fiber_lengths);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnums = find_fibers(parms,fiber,vol_ROI)
  % get info about fibers and ROI
  nfibers = fiber.header.nFiberNr;
  volsz = size(vol_ROI);
  ROI_vox = find(vol_ROI>parms.thresh);
  nvox = length(ROI_vox);
  % find overlap between fibers and ROI
  selected = zeros(nfibers,1);
  for f=1:nfibers
    chain_xyz = ceil(fiber.chain(f).pxyzChain+1);
    for d=1:3
      chain_xyz(:,d) = min(chain_xyz(:,d),volsz(d));
    end;    
    chain_vox = sub2ind(volsz,chain_xyz(:,1),chain_xyz(:,2),chain_xyz(:,3));
    if ~isempty(intersect(chain_vox,ROI_vox)), selected(f) = 1; end;
  end;
  fnums = find(selected);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fiber = cut_fibers(parms,fiber)
  volsz = size(parms.vol_ROI{1});
  nfibers = fiber.header.nFiberNr;
  for f=1:nfibers
    % get voxel indices for each point along streamline
    chain_xyz = ceil(fiber.chain(f).pxyzChain+1);
    for d=1:3
      chain_xyz(:,d) = min(chain_xyz(:,d),volsz(d));
    end;    
    chain_vox = sub2ind(volsz,chain_xyz(:,1),chain_xyz(:,2),chain_xyz(:,3));
    chain_length = size(chain_xyz,1);
    % find ROI intersections for each point along streamline
    chain_roi_vec = zeros(chain_length,parms.nroi);
    % identify points along streamline for each roi
    for r=1:parms.nroi
      chain_roi_vec(:,r) = 1.0*(parms.vol_ROI{r}(chain_vox)>parms.thresh);
    end;
    % determine direction of streamline relative to order of rois
    chain_roi_first = chain_roi_vec(:,1);
    chain_roi_last = chain_roi_vec(:,parms.nroi);
    ind_roi_first = find(chain_roi_first);
    ind_roi_last = find(chain_roi_last);
    if min(ind_roi_last) < min(ind_roi_first)
      % swap
      tmp = ind_roi_last;
      ind_roi_last = ind_roi_first;
      ind_roi_first = tmp;
      tmp = chain_roi_last;
      chain_roi_last = chain_roi_first;
      chain_roi_first = tmp;
    end;
    % find beginning of chain (last time entering first roi)
    tmp_start = ind_roi_first(1);
    tmp_stop = ind_roi_first(end);
    roi_range = [tmp_start:tmp_stop];
    ind_zero = find(chain_roi_first(roi_range)==0);
    if isempty(ind_zero)
      chain_start = tmp_start;
    else
      chain_start = roi_range(ind_zero(end))+1;
    end;
    % find end of chain (first time leaving last roi)
    tmp_start = ind_roi_last(1);
    tmp_stop = ind_roi_last(end);
    roi_range = [tmp_start:tmp_stop];
    ind_zero = find(chain_roi_last(roi_range)==0);
    if isempty(ind_zero)
      chain_stop = tmp_stop;
    else
      chain_stop = roi_range(ind_zero(1))-1;
    end;
    % trim ends of streamline
    chain_xyz = fiber.chain(f).pxyzChain(chain_start:chain_stop,:);
    % update fiber
    fiber.chain(f).nLength = size(chain_xyz,1);
    fiber.chain(f).nSelEndIdx = fiber.chain(f).nLength - 1;
    fiber.chain(f).pxyzChain = chain_xyz;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vol_ROI = load_ROI(fname_ROI,revsliceflag)
  [tpath,tstem,text] = fileparts(fname_ROI);
  switch text
    case {'.mgh','.mgz'}
      vol_ROI = fs_load_mgh(fname_ROI);
    case '.mat'
      vol_ROI = mmil_load_sparse(fname_ROI);
    otherwise
      error('invalid ROI file type');
  end;
  if revsliceflag
    vol_ROI = vol_ROI(:,:,end:-1:1);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fiber_numbers = list_fiber_numbers(number,length)
  fiber_numbers = number*ones(length,1);
return;


