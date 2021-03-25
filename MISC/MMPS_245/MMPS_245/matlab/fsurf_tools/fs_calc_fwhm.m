function [fwhm,nverts,surfarea] = fs_calc_fwhm(subj,lh_fnames,rh_fnames,varargin);
%function [fwhm,nverts,surfarea] = fs_calc_fwhm(subj,lh_fnames,rh_fnames,[options]);
%
% Usage:
%  [fwhm] = fs_calc_fwhm(subj,lh_fnames,rh_fnames,'key1', value1,...);
%
% Required input:
%  subj: subject ID (freesurfer recon directory name)
%  lh_fnames: full path name of input surface stats file(s) for left hemi
%                  (surface mgh format)
%  rh_fnames: full path name of input surface stats file(s) for right hemi
%    NOTE: If a cell array of multiple file names is supplied, 
%            the normalized residual error will be calculated for each,
%            and FWHM smoothness estimated from those
%          This feature is designed for surface-based group analyses.
%          Each stats file should be sampled onto icosahedral sphere and
%            the subject ID should be fsaverage or a custom average subject
%    NOTE: Either (but not both) can be supplied as empty ([]) to restrict
%            calculations to a single hemisphere. 
%
% Optional parameters:
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'surfname' - surface used for calculations
%    {default: white}
%  'surf_matfile' - matfile containing surfaces
%    if this file does not exist, will be created
%    if this file does exist, will be loaded, skipping some surface loading
%    {default: []}
%  'mask_midbrain_flag' - [0|1] toggle mask out thalamus, corpus collosum
%    {default: 1}
%
% Created:  07/05/07 Don Hagler
% Last Mod: 09/13/12 by Don Hagler
%

fwhm = NaN;
nverts = 0;
surfarea = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), g=struct(options{:}); 
  else g = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

if ~isfield(g,'subjdir'), g.subjdir = []; end;
if ~isfield(g,'surfname'), g.surfname = 'white'; end;
if ~isfield(g,'surf_matfile'), g.surf_matfile = []; end;
if ~isfield(g,'mask_midbrain_flag'), g.mask_midbrain_flag = 1; end;

gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'subjdir' 'surfname' 'surf_matfile'  'mask_midbrain_flag'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end;
end;

% get rid of options struct
subjdir = g.subjdir;
surfname = g.surfname;
surf_matfile = g.surf_matfile;
mask_midbrain_flag = g.mask_midbrain_flag;
clear g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

if nargin<3, help(mfilename); return; end;  

if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    fprintf('%s: ERROR: SUBJECTS_DIR not defined as an environment variable\n',mfilename);
    return;
  end;
end;

if iscell(subj) & length(subj)>1
  fprintf('%s: ERROR: can only specify one subj\n',mfilename);
  return;
end;

if ~iscell(lh_fnames), lh_fnames = {lh_fnames}; end;
if ~iscell(rh_fnames), rh_fnames = {rh_fnames}; end;

hemilist = {'lh','rh'};
if isempty(lh_fnames) & isempty(rh_fnames)
  fprintf('%s: ERROR: lh_fnames and rh_fnames cannot both be empty\n',...
    mfilename);
  return;
elseif isempty(lh_fnames)
  hemilist = {'rh'};
  hemi_data(1).fnames = rh_fnames;
  num_inputs = length(rh_fnames);
elseif isempty(rh_fnames)
  hemilist = {'lh'};
  hemi_data(1).fnames = lh_fnames;
  num_inputs = length(lh_fnames);
elseif length(lh_fnames) ~= length(rh_fnames)
  fprintf('%s: ERROR: lh_fnames and rh_fnames must have same length\n',...
    mfilename);
  return;
else
  hemi_data(1).fnames = lh_fnames;
  hemi_data(2).fnames = rh_fnames;
  num_inputs = length(lh_fnames);
end;

% check files exist and are mgh/mgz format
for h=1:length(hemilist)
  hemi = hemilist{h};
  surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subj,hemi,surfname);
  if ~exist(surffile,'file')
    fprintf('%s: ERROR: surface file %s not found\n',mfilename,surffile);
    return;
  end
  for i=1:num_inputs
    if ~exist(hemi_data(h).fnames{i},'file')
      fprintf('%s: ERROR: surface file %s not found\n',...
        mfilename,hemi_data(h).fnames{i});
      return;
    end
    [fpath,fstem,fext] = fileparts(hemi_data(h).fnames{i});
    if ~ismember(fext,{'.mgh','.mgz'})
      fprintf('%s: ERROR: surface files (e.g. %s) must have .mgh or .mgz extension\n',...
        mfilename,hemi_data(h).fnames{i});
      return;
    end;
  end;
end;
setenv('SUBJECTS_DIR',subjdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load subject's surfaces
fprintf('%s: loading surfaces\n',mfilename);
if isempty(surf_matfile) | ~exist(surf_matfile,'file')
  for h=1:length(hemilist)
    hemi = hemilist{h};
    surf = fs_load_subj(subj,hemi,surfname,[],subjdir);
    if isempty(surf)
      fprintf('%s: ERROR: failed to load %s %s surface for %s\n',...
        mfilename,hemi,surfname,subj);
      return;
    end;
    fprintf('%s: loaded %s surface with %d vertices\n',...
      mfilename,hemi,surf.nverts);
    fprintf('%s: computing surface area...\n',mfilename);
    surf = fs_calc_triarea(surf);
    fprintf('%s: average inter-vertex distance = %0.3f mm\n',...
      mfilename,surf.avgdist);
    hemi_surf(h).surf = surf;
  end;
  if ~isempty(surf_matfile)
    save(surf_matfile,'hemi_surf');
  end;
else
  load(surf_matfile);
end;

% get mask
for h=1:length(hemilist)
  hemi_data(h).mask = ones(hemi_surf(h).surf.nverts,1);
  if mask_midbrain_flag
    hemi_data(h).mask = fs_mask_surfstats_with_aparc(hemi_data(h).mask,...
      subj,hemilist{h},subjdir);
  end;
end;

% load data
fprintf('%s: loading data\n',mfilename);
for h=1:length(hemilist)
  hemi = hemilist{h};
  nverts = hemi_surf(h).surf.nverts;
  hemi_data(h).vals = zeros(nverts,num_inputs);
  for i=1:num_inputs
    % load first frame only
    vals = mmil_rowvec(fs_load_mgh(hemi_data(h).fnames{i},[],1));
    if isempty(vals) 
      fprintf('%s: ERROR: failed to load %s\n',...
        mfilename,hemi_data(h).fnames{i});
      return;
    elseif length(vals)~=nverts
      fprintf('%s: ERROR: %s has wrong number of vertices (%d)\n',...
        mfilename,hemi_data(h).fnames{i},length(vals));
      return;
    end;
    hemi_data(h).vals(:,i) = vals;
  end;
end;

% calculate normalized residuals
if num_inputs>1
  fprintf('%s: calculating normalized residuals...\n',mfilename);
  for h=1:length(hemilist)
    hemi_data(h).mean = 0;
    % calculate mean
    hemi_data(h).mean = mean(hemi_data(h).vals,2);
    % calculate residuals
    hemi_data(h).res = hemi_data(h).vals - ...
                       hemi_data(h).mean*ones(1,num_inputs);
    % calculate stdev of residuals
    hemi_data(h).res_stdev = std(hemi_data(h).res,0,2);
    % normalize residuals
    hemi_data(h).res_norm = hemi_data(h).res ./ ...
                            (hemi_data(h).res_stdev*ones(1,num_inputs));
  end;
else
  hemi_data(h).res_norm = hemi_data(h).vals;
end;

% calculate fwhm
fprintf('%s: calculating fwhm...\n',mfilename);
sum_dvals=0; sumsq_dvals=0;
sum_vals=0;  sumsq_vals=0;
total_dv=0;  total_v=0;
var_dv = 0;  var_v = 0;
avgdist = 0; nverts = 0;
surfarea = 0;
for h=1:length(hemilist)
  hemi = hemilist{h};
  surf = hemi_surf(h).surf;

  good_verts = find(hemi_data(h).mask);
  num_good_verts = length(good_verts);

  nverts = nverts + num_good_verts;
  surfarea = surfarea + sum(surf.vertex_area(good_verts));
  avgdist = avgdist + surf.avgdist*num_good_verts;

  nbrs = surf.nbrs;
  maxnbrs = size(nbrs,2);

  % exclude all vertices outside mask  
  bad_verts = find(hemi_data(h).mask==0);
  nbrs(bad_verts,:) = 0;
  bad_nbrs = find(ismember(nbrs(:),bad_verts));
  nbrs(bad_nbrs) = 0; 

  % change zeros to 1, add dummy vertex
  nbrs = [ones(1,maxnbrs);nbrs + 1];

  % check for nbr < vert, set to 1 so only count neighbor pair once
  vnums = (1:surf.nverts+1)';
  diff_nbrs = nbrs - vnums*ones(1,maxnbrs);
  nbrs(find(diff_nbrs<0))=1;

  num_nbrs = sum(nbrs>1,2)+1; % num neighbors for each vertex (including self)
  sum_num_nbrs = sum(num_nbrs); % num neighbors for all vertices

  for i=1:num_inputs
    % calculate interneighbor differences
    vals = hemi_data(h).res_norm(:,i);

    % catch bad values
    vals(bad_verts) = 0;
    nan_vals = find(isnan(vals));
    if ~isempty(nan_vals)
      fprintf('\n%s: WARNING: setting %d NaN valules to 0\n',...
        mfilename,length(nan_vals));
      vals(nan_vals)=0;
    end;
   
    vals = [0;vals]; % create dummy vertex with index 1, value 0
    nbr_vals = vals(nbrs(:,:));
    nbr_val_mask = zeros(size(nbr_vals));
    nbr_val_mask(find(nbr_vals))=1;
    dvals = (vals(nbrs(:,:)) - vals*ones(1,maxnbrs)).*nbr_val_mask;
    sum_dvals = sum_dvals + sum(dvals(:));
    sumsq_dvals = sumsq_dvals + sum(dvals(:).^2);
    total_dv = total_dv + sum_num_nbrs - num_good_verts - 1;
    sum_vals = sum_vals + sum(vals);
    sumsq_vals = sumsq_vals + sum(vals.^2);
    total_v = total_v + num_good_verts;
  end;
end;
avgdist = avgdist/nverts;
fprintf('\n%s: total vertices = %d\n',...
  mfilename,nverts);
fprintf('%s: total surface area = %0.3f mm^2\n',...
  mfilename,surfarea);
fprintf('%s: grand average inter-vertex distance = %0.3f mm\n',...
  mfilename,avgdist);

% calculate variance of interneighbor differences
var_dv = (sumsq_dvals - sum_dvals^2/total_dv)/(total_dv-1);
% calculate variance of all values
var_v = (sumsq_vals - sum_vals^2/total_v)/(total_v-1);
% calculate fwhm
varratio = -log(1.0 - 0.5*(var_dv/(var_v+eps)));
if (varratio<=0)
  varratio = 0.5*(var_dv/(var_v+eps));
end;
if (varratio<=0)
  fwhm=0;
else
 fwhm = sqrt(2.0*log(2.0)/varratio)*avgdist;
end;  

fprintf('\n%s: var_v = %0.4f\n',mfilename,var_dv);
fprintf('%s: var_dv = %0.4f\n',mfilename,var_v);
fprintf('%s: fwhm = %0.4f mm\n',mfilename,fwhm);

