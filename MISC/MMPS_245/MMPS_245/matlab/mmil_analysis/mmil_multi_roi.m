function results = mmil_multi_roi(fname_data,fname_roi_list,varargin)
%function results = mmil_multi_roi(fname_data,fname_roi,[options]);
%
% Required Input:
%  fname_data: full path to data volume (mgh/mgz format)
%  fname_roi: full path to ROI volume (mgh/mgz format)
%    must be registered to and have same resolution as data volume
%    can be a cell array (curly bracketed list) of multiple ROI file names
%
% Optional Input:
%  'frames': for multi-frame data, data from only the specified frames
%     (e.g. time points) will be extracted
%     if empty or omitted, will use all available frames
%    {default = []}
%  'minval': minimum value; if 0, include all voxels except NaNs
%    {default = 1e-6}
%  'weighted_avg_flag': [0|1] toggle calculation of weighted average
%    if 1, use values in ROI volumes to weight the data volume
%    {default = 0}
%  'scalefact': scaling factor applied to fname_data
%    {default = 1}
%  'verbose': [0|1] display status meassages
%    {default = 0}
%
% Output:
%   results: struct array containing data for each ROI with fields:
%     vals: matrix of values with size = [nvals,nframes]
%     weights: vector of weights with size = [nvals,1]
%     nvals: number of values
%     nvals_valid: number of valid values (non NaN and abs(val)>= minval)
%     nvals_invalid: number of invalid values (i.e. NaN or abs(val)<minval)
%     avg: mean of valid values with size = [1,nframes]
%     median: median of valid values with size = [1,nframes]
%     stdv: standard deviation of valid values with size = [1,nframes]
%
% Created:  04/05/07 by Don Hagler
% Last mod: 02/05/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'frames',[],[],...
  'minval',1e-6,[0,1e6],...
  'weighted_avg_flag',false,[false true],...
  'scalefact',1,[1e-10,1e10],...
  'verbose',false,[false true],...
});

results = [];

if ~iscell(fname_roi_list), fname_roi_list = {fname_roi_list}; end;

if ~exist(fname_data,'file')
  error('data file %s not found',fname_data);
end;
data_volsz = read_volsz(fname_data,parms.frames);

for i=1:length(fname_roi_list)
  fname_roi = fname_roi_list{i};
  if ~isempty(fname_roi)
    if ~exist(fname_roi,'file')
      error('ROI file %s not found',fname_roi);
    end;
    roi_volsz = read_volsz(fname_roi);
    if any(data_volsz(1:3)~=roi_volsz(1:3))
      error('ROI vol size for file %s does not match data',fname_roi);
    end;
    if size(roi_volsz,2)>3 && roi_volsz(4)>1
      error('ROI volumes cannot be multi-frame (%s has %d frames)',...
          fname_roi,roi_volsz(4));
    end;
 end;
end;

% load data file
if parms.verbose
  fprintf('%s: loading data...\n',mfilename);
end;
vol_data = load_vol(fname_data,parms.frames);
if parms.scalefact~=1
  vol_data = parms.scalefact*vol_data;
end;
vol_data = reshape(vol_data,[prod(data_volsz(1:3)),data_volsz(4)]); % allow multiframe

if parms.verbose
  fprintf('%s: extracting values...\n',mfilename);
end;
for i=1:length(fname_roi_list)
  fname_roi = fname_roi_list{i};
  if isempty(fname_roi)
    roi = [];
  else
    vol_roi = load_vol(fname_roi);
    roi = find(vol_roi>0);
  end;
  results(i).vals = [];
  results(i).weights = [];
  results(i).nvals = 0;
  results(i).nvals_valid = 0;
  results(i).nvals_invalid = 0;
  results(i).avg = NaN;
  results(i).median = NaN;
  results(i).stdv = NaN;
  if ~isempty(roi)
    vals = vol_data(roi,:);
    nframes = size(vals,2); % allow multiframe
    nvals = size(vals,1);
    ind_valid = find(abs(vals(:,1))>=parms.minval & ~isnan(vals(:,1)));
    ind_invalid = setdiff([1:nvals],ind_valid);
    ind_nan = isnan(vals(:,1));
    vals(ind_nan,:) = 0;
    weights = zeros(nvals,1);
    weights(ind_invalid) = 0;
    results(i).vals = vals;
    results(i).weights = weights;
    results(i).nvals = nvals;
    results(i).nvals_valid = length(ind_valid);
    results(i).nvals_invalid = length(ind_invalid);
    if results(i).nvals_valid > 0
      if parms.weighted_avg_flag
        weights = vol_roi(roi);
        weights(ind_invalid) = 0;
        results(i).weights = weights;
        weights = repmat(weights,[1,nframes]);
        results(i).avg = mmil_wtd_mean(vals,weights,1);
        results(i).median = mmil_wtd_median(vals,weights,1);
        if results(i).nvals_valid>1
          results(i).stdv = mmil_wtd_std(vals,weights,1);
        end;
      else
        tmp_vals = vals(ind_valid,:);
        results(i).avg = mean(tmp_vals,1);
        results(i).median = median(tmp_vals,1);
        if results(i).nvals_valid>1
          results(i).stdv = std(tmp_vals,0,1);
        end;
      end;
    end;
  end;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,M,volsz] = load_vol(fname,frames)
  vol = []; M = []; volsz = [];
  if ~exist('frames','var'), frames = []; end;
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [vol,M,volsz] = mmil_load_sparse(fname);
      if ~isempty(frames)
        vol = vol(:,:,:,frames);
        volsz(4) = 1;
      end;
    case {'.mgh','.mgz'}
      [vol,M,tmp,volsz] = fs_load_mgh(fname,[],frames);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [volsz,M] = read_volsz(fname,frames)
  volsz = []; M = [];
  if ~exist('frames','var'), frames = []; end;
  [tmp,tmp,ext] = fileparts(fname);
  switch ext
    case '.mat'
      [tmp,M,volsz] = mmil_load_sparse(fname);
    case {'.mgh','.mgz'}
      [M,volsz] = fs_read_header(fname);
  end;
  if ~isempty(frames)
    volsz(4) = length(frames);
  end;
return;

