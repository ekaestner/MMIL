function results = mmil_surf_roi(fname_data,varargin);
%function results = mmil_surf_roi(fname_data,[options]);
%
% Purpose: extract average values or time courses for surface ROIs
%
% Usage:
%  mmil_surf_roi(fname_data,'key1', value1,...); 
%
% Required Parameters:
%  fname_data: full path to surface data (mgh format)
%
% Optional Input (entered as 'key',value pairs):
%  'fname_aparc': full path to aparc annotation file
%    (e.g. subjdir/subj/label/lh.aparc.annot)
%    {default = []}
%  'fname_label': full path name of FreeSurfer label file
%    can be a cell array (curly bracketed list) of multiple label file names
%    {default = []}
%  'fname_weights': full path name of mgh file containing weighted ROI
%     may be multi-frame for multiple ROIs
%    {default = []}
%  'frames': for multi-frame data, data from only the specified frames
%     (e.g. time points) will be extracted
%     if empty or omitted, will use all available frames
%    {default = []}
%  'minval': minimum value; if 0, include all voxels except NaNs
%    {default = 1e-6}
%  'scalefact': scaling factor applied to extracted values
%    {default = 1}
%  'fname_colorlut': full path of color look up table file
%    if empty, use FreeSurferColorLUT.txt in $FREESURFER_HOME
%    {default = []}
%  'hemi': cortical hemisphere
%     used to get ROI codes for aparc ROIs from color LUT
%     if empty, will infer from fname_aparc or fname_weights
%    {default = []}
%  'weights_thresh': threshold applied to weights file
%    {default = 0}
%  'verbose': [0|1] display status messages
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
% Note: must supply at least one of fname_aparc, fname_label, or fname_weights
%
% Created:  07/14/08 by Don Hagler
% Prev Mod: 09/19/17 by Don Hagler
% Last Mod: 10/24/17 by Don Hagler
%

results = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_aparc',[],[],...
  'fname_label',[],[],...
  'fname_weights',[],[],...
  'frames',[],[],...
  'minval',1e-6,[0,1e6],...
  'scalefact',1,[1e-10,1e10],...
  'fname_colorlut',[],[],...
  'hemi',[],{'lh','rh'},...
  'weights_thresh',0,[],...
  'verbose',false,[false true],...
...
  'annot_name',[],[],...
});

if isempty(parms.fname_aparc) &&...
   isempty(parms.fname_label) &&...
   isempty(parms.fname_weights)
  error('must supply fname_aparc, fname_label, or fname_weights as input ROIs');
end;

% check input file
if ~exist(fname_data,'file')
  error('error surface data file %s not found\n',fname_data);
end;

% get volsz and check frames
[tmp,volsz] = fs_read_header(fname_data);
if volsz(2)~=1 || volsz(3)~=1
  error('input data must be on surface, not volume');
end;
if min(parms.frames)<1
  error('frame numbers must be at least 1');
end;
if max(parms.frames)>volsz(4)
  error('frame numbers must be < number of frames (%d) in %s',volsz(4),fname_data);
end;
nverts = volsz(1);

nrois = 0;
roiverts = {};
roiweights = {};
roinames = {};
roicodes = [];

% load aparc ROIs
if ~isempty(parms.fname_aparc)
  if isempty(parms.hemi)
    % decide if fname_aparc is for lh or rh
    [tpath,tstem,text] = fileparts(parms.fname_aparc);
    n = regexp(tstem,'(?<hemi>[lr]h)\.(?<name>.+)','names');
    if isempty(n)
      fprintf('%s: WARNING: annot file name %s has unexpected pattern...\n',...
        mfilename,parms.fname_aparc);
    else
      parms.hemi = n.hemi;
      parms.annot_name = n.name;
    end;
  end;

  if ~isempty(parms.annot_name) && ~isempty(regexp(parms.annot_name,'aparc'))
    % load freesurfer LUT
    if parms.verbose
      fprintf('%s: loading freesurfer color lut...\n',mfilename);
    end;
    [all_fs_roicodes,all_fs_roinames] = fs_colorlut(parms.fname_colorlut);
    switch parms.hemi
      case 'lh' % lh roi codes and names
        [tmp,ind]=intersect(all_fs_roicodes,[1000:1300]);
      case 'rh' % rh roi codes and names
        [tmp,ind]=intersect(all_fs_roicodes,[2000:2300]);
    end;
    fs_roicodes = all_fs_roicodes(ind);
    fs_roinames = all_fs_roinames(ind);
  else
    fs_roicodes = [];
    fs_roinames = [];
  end;

  if parms.verbose
    fprintf('%s: loading aparc annotation from %s...\n',...
      mfilename,parms.fname_aparc);
  end;
  [aparc_roinums,aparc_names] = fs_read_annotation(parms.fname_aparc);
  if nverts~=length(aparc_roinums)
    error('number of data points in fname_data (%d) does not match aparc (%d)',...
      nverts,length(aparc_roinums));
  end;
  for i=1:length(aparc_names)
    nrois = nrois + 1;  
    roiverts{nrois} = find(aparc_roinums==i);
    roiweights{nrois} = [];
    if ~isempty(parms.hemi)
      roinames{nrois} = ['ctx-' parms.hemi '-' aparc_names{i}];
      switch parms.hemi
        case 'lh'
          roicode_offset = 21100;
        case 'rh'
          roicode_offset = 22100;
      end;
    else
      roinames{nrois} = aparc_names{i};
      roicode_offset = 23100;
    end;
    if isempty(fs_roicodes)
      roicodes{nrois} = roicode_offset + nrois - 1;
    else
      ind = find(strcmp(roinames{nrois},fs_roinames));
      if isempty(ind)
        roicodes{nrois} = roicode_offset + nrois - 1;
      else
        roicodes{nrois} = fs_roicodes(ind);
      end;
    end;
  end;
end;

% load label files
if ~isempty(parms.fname_label)
  if parms.verbose
    fprintf('%s: loading label files...\n',mfilename);
  end;
  if ~iscell(parms.fname_label), parms.fname_label = {parms.fname_label}; end;
  for i=1:length(parms.fname_label)
    fname_label = parms.fname_label{i};
    if ~exist(fname_label,'file')
      error('label file %s not found',fname_label);
    end;
    [tpath,tstem,text] = fileparts(fname_label);
    if ~strcmp(text,'.label')
      error('label file %s does not have .label extension',fname_label);
    end;
    nrois = nrois + 1;  
    roinames{nrois} = tstem;
    roicodes{nrois} = nrois;
    roiverts{nrois} = fs_read_label(fname_label);
    roiweights{nrois} = [];
    if max(roiverts{nrois})>nverts | min(roiverts{nrois})<1
      error('bad label vertex number');
    end;
  end;
end;

% load weights file
if ~isempty(parms.fname_weights)
  if parms.verbose
    fprintf('%s: loading weights files...\n',mfilename);
  end;
  [tpath,fstem,text] = fileparts(parms.fname_weights);
  if ~ismember(text,{'.mgh','.mgz'})
    error('weights file %s does not have .mgh or .mgz extension',...
      parms.fname_weights);
  end;
  n = regexp(fstem,'(?<stem>\w+)-(?<hemi>[lr]h$)','names');
  if ~isempty(n)
    fstem = ['ctx-' n.hemi '-' n.stem];
  elseif ~isempty(parms.hemi)
    fstem = ['ctx-' parms.hemi '-' fstem];
  else
    fstem = ['ctx-' fstem];
  end;
  weights = squeeze(fs_load_mgh(parms.fname_weights));
  if size(weights,1)~=nverts
    error('weights file %s has wrong number of vertices (%d instead of %d)',...
      fname_weights,size(weights,1),nverts);
  end;
  nwr = size(weights,2);
  for i=1:nwr
    nrois = nrois + 1;
    roinames{nrois} = sprintf('%s_%d',fstem,i);
    roicodes{nrois} = nrois;
    tmp_weights = weights(:,i);
    tmp_verts = find(tmp_weights>parms.weights_thresh);
    roiverts{nrois} = tmp_verts;
    roiweights{nrois} = tmp_weights(tmp_verts);
  end;
end;

% load surface data
if parms.verbose
  fprintf('%s: loading functional surface data from %s...\n',...
    mfilename,fname_data);
end;
surfdata = squeeze(double(fs_load_mgh(fname_data,[],parms.frames)));
parms.nframes = size(surfdata,2);

if parms.scalefact~=1
  surfdata = parms.scalefact*surfdata;
end;

if parms.verbose
  fprintf('%s: extracting values...\n',mfilename);
end;
for i=1:length(roiverts)
  results(i).roiname = roinames{i};
  results(i).roicode = roicodes{i};
  roi = roiverts{i};
  results(i).vals = [];
  results(i).weights = [];
  results(i).nvals = 0;
  results(i).nvals_valid = zeros(1,parms.nframes);
  results(i).nvals_invalid = zeros(1,parms.nframes);
  results(i).avg = nan(1,parms.nframes);
  results(i).median = nan(1,parms.nframes);
  results(i).stdv = nan(1,parms.nframes);
  if ~isempty(roi)
    vals = surfdata(roi,:);
    nvals = size(vals,1);
    results(i).vals = vals;
    results(i).nvals = nvals;
    results(i).weights = ones(size(vals));
    for f=1:parms.nframes
      fvals = vals(:,f);
      ind_valid = find(abs(fvals)>=parms.minval & ~isnan(fvals));
      ind_invalid = setdiff([1:nvals],ind_valid);
      ind_nan = isnan(fvals);
      fvals(ind_nan,:) = 0;
      if isempty(roiweights{i})
        weights = ones(nvals,1);
      else
        weights = roiweights{i};
      end;
      weights(ind_invalid) = 0;
      results(i).weights(:,f) = weights;
      results(i).nvals_valid(f) = length(ind_valid);
      results(i).nvals_invalid(f) = length(ind_invalid);
      if results(i).nvals_valid(f) > 0
        if ~isempty(roiweights{i})
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
