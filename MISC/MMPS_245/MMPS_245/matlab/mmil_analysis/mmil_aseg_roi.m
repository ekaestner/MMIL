function results = mmil_aseg_roi(fname_vals,fname_aseg,varargin)
%function results = mmil_aseg_roi(fname_vals,fname_aseg,[options])
%
% Required Input:
%  fname_vals: full path to volume with values to be extracted (mgh/mgz format)
%   (must be registered and resampled to T1)
%  fname_aseg: full path to segmentation volume (e.g. aseg.mgz, aparc+aseg.mgz)
%
% Optional Input:
%  'fname_weights': full path to weights volume (mgh/mgz format)
%    if supplied, summary statistics will be weighted by these values
%    i.e. avg will be a weighted mean, median will be a weighted median
%    {default = []}
%  'fname_mask': full path to mask volume (mgh/mgz format)
%    if supplied, values outside mask are set to zero
%    {default = []}
%  'roilist': list of acceptable ROI codes
%    The ROI codes are the values in the segmentation volume
%      and are also found in the FreeSurferColorLUT.txt
%    If this is empty, only the ROI codes in the segmentation volume will be used
%    If an ROI code is specified but not found in the segmentation volume,
%      zero values will be returned for that ROI
%    If an ROI code is specified but not found in the LUT, that ROI code will be ignored
%    {default = []}
%  'roigroups': struct array containing fields 'roicode' and 'roiname'
%    it will only process these ROIs, and will include regardless if found or
%    not in FreeSurferColorLUT.txt
%    {default = []}
%  'minval': minimum value for inclusion in summary statistics
%    if 0, include all voxels except NaNs
%    {default = 1e-6}
%  'mask_thresh': threshold value applied to mask volume
%    {default = 1e-6}
%  'scalefact': scaling factor applied to extracted values
%    {default = 1}
%  'frames': vector of frame numbers to extract
%    If empty, use all
%    {default = []}
%  'fname_colorlut': full path of color look up table file
%     if empty, use FreeSurferColorLUT.txt in $FREESURFER_HOME
%    {default = []}
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
% Created:  10/02/06 by Don Hagler
% Prev Mod: 09/19/17 by Don Hagler
% Last Mod: 10/31/17 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
results = [];

parms = check_input(fname_vals,fname_aseg,varargin);

[parms,data] = load_data(parms);

parms = get_roi_info(parms,data);

results = extract_values(parms,data);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_vals,fname_aseg,options)
  % parse input arguments
  parms = mmil_args2parms(options,{...
    'fname_vals',fname_vals,[],...
    'fname_aseg',fname_aseg,[],...
  ...
    'fname_weights',[],[],...
    'fname_mask',[],[],...
    'roilist',[],[],...
    'roigroups',[],[],...
    'minval',1e-6,[0,1e6],...
    'mask_thresh',1e-6,[],...
    'scalefact',1,[1e-10,1e10],...
    'frames',[],[],...
    'fname_colorlut',[],[],...
    'verbose',false,[false true],...
  });
  % check input files
  if ~exist(parms.fname_aseg,'file')
    error('segmentation volume file %s not found',parms.fname_aseg);
  end;
  if ~exist(parms.fname_vals,'file')
    error('functional volume file %s not found',parms.fname_vals);
  end;
  if ~isempty(parms.fname_weights) & ~exist(parms.fname_weights,'file')
    error('weights volume file %s not found',parms.fname_weights);
  end;
  if ~isempty(parms.fname_mask) & ~exist(parms.fname_mask,'file')
    error('mask volume file %s not found',parms.fname_mask);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,data] = load_data(parms)
  data = [];
  
  % load functional volume
  if parms.verbose
    fprintf('%s: loading functional volume %s...\n',mfilename,parms.fname_vals);
  end;
  [data.vol_vals,parms.M,parms.volsz,parms.nframes] = ...
    load_vol(parms.fname_vals,parms.frames);
  if parms.scalefact~=1
    data.vol_vals = parms.scalefact*data.vol_vals;
  end;

  % load segmentation volume
  if parms.verbose
    fprintf('%s: loading segmentation volume %s...\n',...
      mfilename,parms.fname_aseg);
  end;
  [data.vol_aseg,M_aseg,volsz_aseg] = load_vol(parms.fname_aseg);
  if any(parms.volsz ~= volsz_aseg) || any(parms.M(:) ~= M_aseg(:))
    error('mismatch between vals and aseg volumes');
  end;

  % load weights volume
  if ~isempty(parms.fname_weights)
    if parms.verbose
      fprintf('%s: loading weights volume %s...\n',...
        mfilename,parms.fname_weights);
    end;
    [data.vol_weights,M_weights,volsz_weights] = ...
      load_vol(parms.fname_weights);
    if any(parms.volsz ~= volsz_weights) || any(parms.M(:) ~= M_weights(:))
      error('mismatch between vals and weights volumes');
    end;
    % reshape to 2D matrix
    data.vol_weights = reshape(data.vol_weights,[prod(parms.volsz),1]);
    parms.weighted_avg_flag = 1;
  else
    parms.weighted_avg_flag = 0;
  end;

  % load mask volume
  if ~isempty(parms.fname_mask)
    if parms.verbose
      fprintf('%s: loading mask volume %s...\n',mfilename,parms.fname_mask);
    end;
    [data.vol_mask,M_mask,volsz_mask] = ...
      load_vol(parms.fname_mask);
    if any(parms.volsz ~= volsz_mask) || any(parms.M(:) ~= M_mask(:))
      error('mismatch between vals and mask volumes');
    end;
    % reshape to 2D matrix
    data.vol_mask = reshape(data.vol_mask,[prod(parms.volsz),1]);
    % threshold mask
    data.vol_mask = 1.0*(data.vol_mask>parms.mask_thresh);
    parms.mask_flag = 1;
  else
    parms.mask_flag = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = get_roi_info(parms,data)
  % load freesurfer LUT
  if parms.verbose
    fprintf('%s: loading freesurfer color lut...\n',mfilename);
  end;
  [roicodes,roinames] = fs_colorlut(parms.fname_colorlut);

  % find unique roicodes from color lut
  [roicodes,ind] = unique(roicodes);
  roinames = roinames(ind);
  nrois = length(roicodes);

  % find unique roi numbers in vol_aseg
  aseg_roicodes = unique(data.vol_aseg(:));
  parms.aseg_roicodes = aseg_roicodes(find(aseg_roicodes>0));
  parms.nsegrois = length(parms.aseg_roicodes);

  if isempty(parms.roilist), parms.roilist = parms.aseg_roicodes; end;
  [parms.roilist,ia,ib] = intersect(parms.roilist,roicodes);
  parms.roicodes=roicodes(ib);
  parms.roinames=roinames(ib);

  if isempty(parms.roigroups)
    parms.nresults = length(parms.roicodes);
  else
    parms.nresults = length(parms.roigroups);
  end;

  if ~isempty(parms.roigroups) && isfield(parms.roigroups,'roi')
    for i=1:parms.nresults
      parms.roigroups(i).roicode = parms.roigroups(i).roi;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = extract_values(parms,data)
  if parms.verbose
    fprintf('%s: extracting values...\n',mfilename);
  end;
  for i=1:parms.nresults
    invalid_flag = 0;
    if isempty(parms.roigroups)
      results(i).roicode = parms.roicodes(i);
      results(i).roiname = parms.roinames{i};
      i_roi = find(parms.aseg_roicodes==results(i).roicode);
      roi = find(data.vol_aseg==results(i).roicode);
      if isempty(i_roi), invalid_flag = 1; end;
    else
      results(i).roicode = parms.roigroups(i).roicode;
      results(i).roiname = deblank(parms.roigroups(i).roiname);
      i_roi = intersect(parms.aseg_roicodes,parms.roigroups(i).roicode);
      roi = find(ismember(data.vol_aseg(:),parms.roigroups(i).roicode));
      % if aseg file has 0 voxels for any ROI in parms.roigroups, do not exclude
      missing_rois = setdiff(parms.roigroups(i).roicode,parms.aseg_roicodes); 
      if ~isempty(missing_rois)
        invalid_flag = 1; 
      end;
    end;
    results(i).vals = [];
    results(i).weights = [];
    results(i).nvals = 0;
    results(i).nvals_valid = zeros(1,parms.nframes);
    results(i).nvals_invalid = zeros(1,parms.nframes);
    results(i).avg = nan(1,parms.nframes);
    results(i).median = nan(1,parms.nframes);
    results(i).stdv = nan(1,parms.nframes);
    if ~invalid_flag
      nvals = length(roi);
      results(i).nvals = nvals;
      results(i).vals = zeros(nvals,parms.nframes);
      results(i).weights = zeros(nvals,parms.nframes);
      for f=1:parms.nframes
        fvals = data.vol_vals(:,:,:,f);
        fvals = reshape(fvals,[prod(parms.volsz),1]);
        fvals = fvals(roi);
        results(i).vals(:,f) = fvals;
        ind_valid = find(abs(fvals)>=parms.minval & ~isnan(fvals));
        if parms.mask_flag
          ind_mask = find(data.vol_mask(roi));
          ind_valid = union(ind_valid,ind_mask);
        end;
        ind_invalid = setdiff([1:nvals],ind_valid);
        ind_nan = isnan(fvals);
        fvals(ind_nan,:) = 0;
        if parms.weighted_avg_flag
          weights = data.vol_weights(roi);
        else
          weights = ones(nvals,1);
        end;
        weights(ind_invalid) = 0;
        results(i).weights(:,f) = weights;
        results(i).nvals_valid(f) = length(ind_valid);
        results(i).nvals_invalid(f) = length(ind_invalid);
        if results(i).nvals_valid(f)>0
          if parms.weighted_avg_flag
            results(i).avg(f) = mmil_wtd_mean(fvals,weights);
            results(i).median(f) = mmil_wtd_median(fvals,weights);
            if results(i).nvals_valid(f)>1
              results(i).stdv(f) = mmil_wtd_std(fvals,weights);
            end;
          else
            tvals = fvals(ind_valid);
            results(i).avg(f) = mean(tvals);
            results(i).median(f) = median(tvals);
            if results(i).nvals_valid(f)>1
              results(i).stdv(f) = std(tvals);
            end;
          end;
        end;
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vol,M,volsz,nframes] = load_vol(fname,frames)
  if ~exist('frames','var'), frames = 1; end;
  vol = []; M = []; volsz = []; nframes = [];
  [vol,M,tmp,volsz] = fs_load_mgh(fname,[],frames);
  volsz = volsz(1:3);
  nframes = size(vol,4);
return;

