function mmil_aseg_weights(fname_aseg,fname_vals,dispersion,varargin)
%function mmil_aseg_weights(fname_aseg,fname_vals,dispersion,[options])
%
% Purpose: calculate weighting factors
%
% Required Input:
%  fname_aseg: full path of segmentation volume file (e.g. aseg.mgz)
%  fname_vals: full path of value file with same resolution
%  dispersion: vector of dispersion values
%     if roicodes is not supplied, disp_vals must match number of
%     unique values in aseg volume
%
% Optional Input:
%  'fname_out': output file name
%    {default = 'weights.mgz'}
%  'scalefact': scaling factor applied to values in fname_vals
%    {default = 1}
%  'dispfact': multiplicative factor applied to dispersion values
%    {default = 4}
%  'erode_flag': [0|1] whether to calculate weights only in shell of
%    ROIs, after eroding
%    {default = 0}
%  'erode_nvoxels': number of voxels to erode
%    ignored if erode_flag = 0 or if fname_aseg_erode exists
%  'fname_aseg_erode': full path to eroded aseg volume (mgh/mgz format)
%    if not supplied or does not exist, one will be generated from fname_aseg
%    {default = []}
%  'roicodes': vector of ROI codes corresponding to values in fname_aseg
%    must correspond to number of elements in disp_vals
%    {default = []}
%  'verbose': [0|1] display status meassages
%    {default = 0}
%  'forceflag': [0|1] overwrite existing output
%    {default = 0}
%
% Created:  10/23/12 by Don Hagler
% Last Mod: 10/23/12 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,3), return; end;

parms = check_input(fname_aseg,fname_vals,dispersion,varargin);

if exist(parms.fname_out,'file') && ~parms.forceflag, return; end;

[parms,data] = load_data(parms);

calc_weights(parms,data);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_aseg,fname_vals,dispersion,options)
  % parse input arguments
  parms = mmil_args2parms(options,{...
    'fname_aseg',fname_aseg,[],...
    'fname_vals',fname_vals,[],...
    'dispersion',dispersion,[],...
  ...
    'fname_out','weights.mgz',[],...
    'scalefact',1,[],...
    'dispfact',4,[1e-6,1e6],...
    'erode_flag',false,[false true],...
    'erode_nvoxels',1,[1,10],...
    'fname_aseg_erode',[],[],...
    'roicodes',[],[],...
    'verbose',false,[false true],...
    'forceflag',false,[false true],...
  });
  if ~exist(fname_aseg,'file')
    error('aseg file %s not found',fname_aseg);
  end;
  if ~exist(fname_vals,'file')
    error('vals file %s not found',fname);
  end;
  % created output directory
  parms.outdir = fileparts(parms.fname_out);
  if isempty(parms.outdir), parms.outdir = pwd; end;  
  mmil_mkdir(parms.outdir);
  % check that volume sizes and coordinate systems match
  [volsz_aseg,M_aseg] = read_volsz(parms.fname_aseg,parms);
  [volsz_vals,M_vals] = read_volsz(parms.fname_vals,parms);
  if any(volsz_aseg~=volsz_vals) || any(M_aseg(:)~=M_vals(:))
    error('mismatch in volumes between aseg and vals');
  end;
  parms.M = M_aseg;
  % check that roicodes and dispersion vectors match
  if ~isempty(parms.roicodes)
    if length(parms.roicodes) ~= length(parms.dispersion)
      error('length of dispersion does not match number of roicodes');
    end;
    parms.nrois = length(parms.roicodes);
  end;
  % scale dispersion values
  parms.dispersion = parms.dispfact * parms.dispersion;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,data] = load_data(parms)
  data = [];
  data.vol_aseg_erode = [];
  data.vol_aseg = fs_load_mgh(parms.fname_aseg);
  if isempty(parms.roicodes)
    parms.roicodes = setdiff(unique(data.vol_aseg(:)),0);
    if length(parms.roicodes) ~= length(parms.dispersion)
      error('length of dispersion does not match number of roicodes');
    end;
    parms.nrois = length(parms.roicodes);
  end;
  if parms.erode_flag
    if isempty(parms.fname_aseg_erode)
      [tmp,fstem] = fileparts(parms.fname_aseg);
      parms.fname_aseg_erode = [parms.outdir '/' fstem '_erode'];
      if parms.erode_nvoxels>1
        parms.fname_aseg_erode = ...
          sprintf('%s%d',parms.fname_aseg_erode,parms.erode_nvoxels);
      end;
      parms.fname_aseg_erode = [parms.fname_aseg_erode '.mgz'];  
    end;
    if ~exist(parms.fname_aseg_erode)
      fs_erode_aseg(parms.fname_aseg,parms.fname_aseg_erode,...
        'nvoxels',parms.erode_nvoxels,...
        'verbose',parms.verbose,'forceflag',parms.forceflag);
    end;
    data.vol_aseg_erode = fs_load_mgh(parms.fname_aseg_erode);
  end;
  data.vol_vals = fs_load_mgh(parms.fname_vals);
  if parms.scalefact ~= 1
    data.vol_vals = parms.scalefact * data.vol_vals;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_weights(parms,data)
  weightfun = @(x,mu,nu) max(0,1-((x-mu)/nu).^2); % Tukey's bisquare function
  vol_weights = 1.0 * (data.vol_aseg > 0);
  for r=1:parms.nrois
    % get values within ROI
    roicode = parms.roicodes(r);
    ind_roi = find(data.vol_aseg == roicode);
    x = data.vol_vals(ind_roi);
    mu = median(x);
    % get values in ROI shell
    if parms.erode_flag
      ind_roi_erode = find(data.vol_aseg_erode == roicode);
      ind_roi = setdiff(ind_roi,ind_roi_erode);
      x = data.vol_vals(ind_roi);
    end;
    % calulate weights
    nu = parms.dispersion(r);
    w = weightfun(x,mu,nu); 
    vol_weights(ind_roi) = w;
  end;
  fs_save_mgh(vol_weights,parms.fname_out,parms.M);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [volsz,M] = read_volsz(fname,parms)
  volsz = []; M = [];
  [M,volsz] = mmil_load_mgh_info(fname,parms.forceflag,parms.outdir);
  volsz = volsz(1:3);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

