function [kernelWidthMax,lambda2] = epi_optimzeB0(fname_for,fname_rev,varargin)
%function [kernelWidthMax,lambda2] = epi_optimizeB0(fname_for,fname_rev,[options])
%
% Purpse: search for optimal settings for estimating B0 distortion
%
% Usage: epi_optimizeB0(fname_for,fname_rev,'key','value'...)
%
% Required Parameters:
%   fname_for: input file name with "forward" phase encode polarity
%   fname_rev: input file name with "reverse" phase encode polarity
%
% Optional Parameters:
%  'swapxy_flag': [0|1] whether to swap x and y dimensions before estimating B0
%    {default = 0}
%  'kernelWidthMax_vec': vector of maximum smoothing kernel widths
%    reasonable values are between 23 and 35
%    {default = [25,31,35]}
%  'lambda2_vec': vector of gradient cost weighting (smoothness)
%    reasonable values are between 900 and 3300
%    {default = [1100,1500,1900]}
%  'lambda2_init': initial value of lambda2 used for optimization of
%     kernelWidthMax (ignored if multi_opt_flag = 1)
%    {default = 1100}
%  'multi_opt_flag': [0|1] optimize parameters simultaneously (i.e. 2D search)
%    Otherwise, first optimize kernelWidthMax, then lambda2
%    {default = 0}
%  'outstem': output file stem for optimization results
%    If empty, will create from fname_for file stem
%    {default = []}
%  'tmpdir': location for temporary directory
%    If empty, will create (and remove when finished) tmp_optimizeB0uw
%      directory in directory containing fname_for
%    {default = []}
%
% Created:  08/16/10 by Don Hagler
% Last Mod: 01/17/15 by Don Hagler
%

% based on code from Matt Erhart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms(varargin,{...
  'swapxy_flag',false,[false true],...
  'kernelWidthMax_vec',[25,31,35],[1:100],...
  'lambda2_vec',[1100,1500,1900],[1:10000],...
  'lambda2_init',1100,[1:10000],...
  'multi_opt_flag',false,[false true],...
  'outstem',[],[],...
  'tmpdir',[],[],...
...
  'cleanupflag',true,[false true],...
  'outext','.mgz',{'.mgz','.mgh'},...
});
kernelWidthMax = [];
lambda2 = [];

% check that input volumes exist
if ~exist(fname_for)
  error('file %s not found',fname_for);
end;
if ~exist(fname_rev)
  error('file %s not found',fname_rev);
end;

% check that input volumes share dimensions and orientations
[dimsmatch,orientmatch] = fs_vol_match(fname_for,fname_rev);
if ~dimsmatch
  error('mismatch in input volume dimensions for %s and %s',...
    fname_for,fname_rev);
end;
if ~orientmatch
  error('mismatch in volume orientation for %s and %s',...
    fname_for,fname_rev);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tmp_path,tmp_fstem,tmp_ext] = fileparts(fname_for);
if isempty(tmp_path), tmp_path = pwd; end;

% create output file stem
if isempty(parms.outstem)
  parms.outstem = [tmp_path '/' tmp_fstem];
end;

% create tmp directory
if isempty(parms.tmpdir)
  parms.tmpdir = tmp_path;
end;
parms.tmpdir = [parms.tmpdir '/tmp_optimizeB0uw'];
mmil_mkdir(parms.tmpdir);

% if swapxy_flag, swap x and y
if parms.swapxy_flag
  orig_orient = fs_read_orient(fname_for);
  pref_orient = orig_orient([2,1,3]);
  fprintf('%s: reorienting input volumes to %s...\n',mfilename,pref_orient);
  tmp_fname_for = [parms.tmpdir '/for' parms.outext];
  tmp_fname_rev = [parms.tmpdir '/rev' parms.outext];
  epi_reorient_vol(fname_for,tmp_fname_for,pref_orient);
  epi_reorient_vol(fname_rev,tmp_fname_rev,pref_orient);
else
  tmp_fname_for = fname_for;
  tmp_fname_rev = fname_rev;
end;

% set tmp output file names
tmp_fname_dx = [parms.tmpdir '/dx' parms.outext];
tmp_fname_for_out = [parms.tmpdir '/for_B0uw' parms.outext];
tmp_fname_rev_out = [parms.tmpdir '/rev_B0uw' parms.outext];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[orig_cost,orig_std,orig_mean] = ...
  uw_cost_function(tmp_fname_for,tmp_fname_rev);

if parms.multi_opt_flag
  % optimize kernelWidthMax and lambda2 simultaneously
  [kernelWidthMax,lambda2,uw_cost,results] = ...
    optimize_parameters(tmp_fname_for,tmp_fname_rev,...
      tmp_fname_dx,tmp_fname_for_out,tmp_fname_rev_out,...
      parms.kernelWidthMax_vec,parms.lambda2_vec);

  % save results
  fname_out = [parms.outstem '_B0uw_optimize_results.mat'];
  save(fname_out,'kernelWidthMax','lambda2',...
    'uw_cost','results',...
    'orig_cost','orig_std','orig_mean');
else
  % optimize kernelWidthMax, holding lambda2 constant
  [kernelWidthMax,lambda2,uw_cost,kresults] = ...
    optimize_parameters(tmp_fname_for,tmp_fname_rev,...
      tmp_fname_dx,tmp_fname_for_out,tmp_fname_rev_out,...
      parms.kernelWidthMax_vec,parms.lambda2_init);

  % optimize lambda2, holding kernelWidthMax constant
  [kernelWidthMax,lambda2,uw_cost,lresults] = ...
    optimize_parameters(tmp_fname_for,tmp_fname_rev,...
      tmp_fname_dx,tmp_fname_for_out,tmp_fname_rev_out,...
      kernelWidthMax,parms.lambda2_vec);

  % save results
  fname_out = [parms.outstem '_B0uw_optimize_results.mat'];
  save(fname_out,'kernelWidthMax','lambda2',...
    'uw_cost','kresults','lresults',...
    'orig_cost','orig_std','orig_mean');
end;

% save text file for easy reference
fname_out = sprintf('%s_B0uw_kWMax_%g_l2_%g_cost_%g',...
  parms.outstem,kernelWidthMax,lambda2,uw_cost);
cmd = ['touch ' fname_out];
[status,result]=unix(cmd);
if status
  fprintf('%s: WARNING: failed to create file %s:\n%s\n',...
    mfilename,fname_out,result);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove tmp dir
if parms.cleanupflag
  cmd = sprintf('rm -r %s',parms.tmpdir);
  [status,result]=unix(cmd);
  if status
    fprintf('%s: WARNING: failed to remove tmp dir:\n%s\n',mfilename,result);
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kernelWidthMax,lambda2,uw_cost,results] = optimize_parameters(...
    fname_for,fname_rev,fname_dx,fname_for_out,fname_rev_out,...
    kernelWidthMax_vec,lambda2_vec)

  nk = length(kernelWidthMax_vec);
  nl = length(lambda2_vec);
  results = zeros(nk*nl,5);
  i = 1;
  for k=1:length(kernelWidthMax_vec)
    kernelWidthMax = kernelWidthMax_vec(k);
    for l=1:length(lambda2_vec)
      lambda2 = lambda2_vec(l);
      epi_estimateB0(fname_for,fname_rev,...
        'fname_dx',fname_dx,...
        'fname_for_out',fname_for_out,...
        'fname_rev_out',fname_rev_out,...
        'kernelWidthMax',kernelWidthMax,...
        'lambda2',lambda2,...
        'forceflag',1);
      [uw_cost,uw_std,uw_mean]=uw_cost_function(fname_for_out,fname_rev_out);
      results(i,1) = kernelWidthMax;
      results(i,2) = lambda2;
      results(i,3) = uw_cost;
      results(i,4) = uw_std;
      results(i,5) = uw_mean;
      i = i+1;
    end;
  end

  % find minimum uw_cost
  [min_cost,ind_min] = min(results(:,3));

  % get best parameters
  kernelWidthMax = results(ind_min,1);
  lambda2 = results(ind_min,2);
  uw_cost = results(ind_min,3);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [uw_cost,uw_std,uw_mean] = uw_cost_function(fname_for,fname_rev)
  % load files (first frame only)
  [vol_for, M_for] = fs_load_mgh(fname_for,[],1);
  [vol_rev, M_rev] = fs_load_mgh(fname_rev,[],1);

  % calculate absolute difference between forward and reverse images
  vol_diff = abs(vol_for - vol_rev);

  % reshape to vector
  vec_diff = reshape(vol_diff,[numel(vol_diff),1]);

  % calculate std, mean, and std+mean as cost of quality of the unwarp
  uw_std = std(vec_diff);
  uw_mean = mean(vec_diff);
  uw_cost = uw_std + uw_mean;

  clear vol_for vol_rev vol_diff vec_diff
return


function rm_files(filelist)
  for f=1:length(filelist)
    if exist(filelist{f},'file')
      delete(filelist{f})
    end;
  end;
return
