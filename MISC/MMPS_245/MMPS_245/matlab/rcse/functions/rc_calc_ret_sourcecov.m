function [P,R]=rc_calc_ret_sourcecov(retmap,varargin);
%function [P,R]=rc_calc_ret_sourcecov(retmap,[options]);
%
% Purpose: create source covariance matrix for retinotopy inverse
%
% Usage:
%  rc_calc_ret_smooth_sourcecov(retmap,'key1', value1,...); 
%
% Depending on options, this may do the following:
%  1. assign priors to visual areas and additional dipoles
%  2. impose a smoothness constraint between nearest neighbor stim locs
%  3. impose a loose orientation constraint for visual area dipoles
%
% Required input:
%  retmap: output of rc_construct_ret_mapping
%
% Optional Input:
%  'indy_locs_flag' - [1|0] whether to calculate independent source estimates
%    for each stimulus location
%    {default: 0}
%  'theta_smfact' - Smoothness factor for independent location estimates
%     (neighbor covariance) related to proximity of polar angle
%     Only applicable when indy_locs_flag=1
%     {default: 0.999}
%  'ecc_smfact' - Smoothness factor for independent location estimates
%     (neighbor covariance) related to proximity of eccentricity
%     Only applicable when indy_locs_flag=1
%     {default: 0.999}
%  'upperlower_smfact' - Smoothness factor across upper and lower visual fields.
%     e.g. smfact=0 for independent sources
%          smfact=1 for identical sources
%    {default: 1}
%  'hemi_smfact' - Smoothness factor across left and right hemifields
%     e.g. smfact=0 for independent sources
%          smfact=1 for identical sources
%    {default: 1}
%  'loose_flag' - [1|0] whether to use loose orientation constraint
%    (allows additional tangential components to dipole)
%    indy_locs_flag must = 1
%    {default: 0}
%  'loose_tang_weight' - source covariance weight of tangential components
%     Only applicable when loose_flag=1
%     {default: 0.1}
%  'visual_area_weight' - source covariance matrix weighting for visual
%     area dipoles (prior estimate of relative strength)
%     {default: 1} (can range from 0-1)
%  'ret_dips_weight' - weighting for extra retinotopic dipoles
%     {default: 1} (can range from 0-1)
%  'nonret_dips_weight' - weighting for extra non-retinotopic dipoles
%     {default: 1} (can range from 0-1)
%
% created  01/04/06 by Don Hagler
% modified 02/19/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'indy_locs_flag',false,[false true],...
  'theta_smfact',0.999,[0 1],...
  'ecc_smfact',0.999,[0 1],...
  'upperlower_smfact',1,[0 1],...
  'hemi_smfact',1,[0 1],...
  'loose_flag',false,[false true],...
  'loose_tang_weight',0.1,[0 1],...
  'visual_area_weight',1,[0 1],...
  'ret_dips_weight',1,[0 1],...
  'nonret_dips_weight',1,[0 1],...
... % hidden parameters
  'theta_exp',8,[0 Inf],... % the larger this number, the more disimilar neighbors will be
  'ecc_exp',0.01,[0 Inf],... % the larger this number, the more disimilar neighbors will be
});

if parms.loose_flag && ~parms.indy_locs_flag,
  error('loose_flag=1 so indy_locs_flag must=1 too');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P=[];
R=[];

fprintf('%s: creating source covariance matrix...\n',mfilename);

% sources for areas
if parms.loose_flag
  area_sources = 3*retmap.num_locs; % 3 orthogonal orientations
elseif parms.indy_locs_flag
  area_sources = retmap.num_locs;
else
  area_sources = 1;
end;
total_sources = area_sources*retmap.num_areas;

% sources for extra dipoles
total_sources = total_sources + 3*retmap.num_ret_dips*retmap.num_ret_dip_locs;
total_sources = total_sources + 3*retmap.num_nonret_dips;

% initialize covariance matrix
P = eye(area_sources);
R = eye(total_sources);

%%% Smoothness Constraint %%%
% set values for P
if parms.indy_locs_flag
  if isfield(retmap,'cond_info') && ~isempty(retmap.cond_info)
    cond_thetas = cell2mat({retmap.cond_info.theta});
    cond_eccs = cell2mat({retmap.cond_info.ecc});
    upperfield_stims = find(cond_thetas>0 & cond_thetas<=180);
    lowerfield_stims = find(cond_thetas>180 | cond_thetas<=0);
    rightfield_stims = find(cond_thetas<90 & cond_thetas>=-90 |...
                            cond_thetas>=270);
    leftfield_stims = find(cond_thetas>=90 & cond_thetas<270 |...
                           cond_thetas<-90);
  end;
  for i=1:retmap.num_locs
    if parms.loose_flag
      s1 = 3*(i-1)+1;
    else
      s1 = i;
    end;

    if isfield(retmap,'cond_info') && ~isempty(retmap.cond_info)
      theta = cond_thetas(i);
      ecc = cond_eccs(i);

      theta_diffs = (theta - cond_thetas)/180;
      theta_diffs(theta_diffs>1) = theta_diffs(theta_diffs>1)-2;
      theta_diffs(theta_diffs<-1) = theta_diffs(theta_diffs<-1)+2;
      theta_weights = parms.theta_smfact.^(parms.theta_exp*abs(theta_diffs));

      ecc_diffs = (ecc - cond_eccs)/max(cond_eccs);
      ecc_weights = parms.ecc_smfact.^(parms.ecc_exp*abs(ecc_diffs));

      upperlower_weights = ones(size(theta_weights));
      if ismember(i,upperfield_stims)
        upperlower_weights(lowerfield_stims) = parms.upperlower_smfact;
      else
        upperlower_weights(upperfield_stims) = parms.upperlower_smfact;
      end;
      
      hemi_weights = ones(size(theta_weights));
      if ismember(i,rightfield_stims)
        hemi_weights(leftfield_stims) = parms.hemi_smfact;
      else
        hemi_weights(rightfield_stims) = parms.hemi_smfact;
      end;

      weights = theta_weights.*ecc_weights.*...
                upperlower_weights.*hemi_weights;
    else
      weights = zeros(1,retmap.num_locs);
      for j=1:retmap.num_locs
        if(j==i)
          weights(j) = 1;
        end
        k=abs(i-j);
        if k>retmap.num_locs/2
          k = retmap.num_locs - k;
        end;
        weights(j) = parms.theta_smfact^k;
      end;      
    end;  
    for j=1:retmap.num_locs
      if(j==i)
        continue;
      end
      if parms.loose_flag
        s2 = 3*(j-1)+1;
      else
        s2 = j;
      end;
      P(s1,s2) = weights(j);
      if parms.loose_flag
        P(s1+1,s1+1) = parms.loose_tang_weight;
        P(s1+2,s1+2) = parms.loose_tang_weight;
      end;
    end;
  end;
end;

% copy P into R multiple times
k=0;
for i=1:retmap.num_areas
  j = k + 1;
  k = j + area_sources - 1;
  R(j:k,j:k) = P*parms.visual_area_weight;
end

% set values for extra ret dipoles
for i=1:retmap.num_ret_dips
  j = k + 1;
  k = j + 3*retmap.num_ret_dip_locs - 1;
  R(j:k,j:k) = parms.ret_dips_weight*eye(3*retmap.num_ret_dip_locs);
end;

% set values for extra nonret dipoles
for i=1:retmap.num_nonret_dips
  j = k + 1;
  k = j + 2;
  R(j:k,j:k) = parms.nonret_dips_weight*eye(3);
end;

fprintf('%s: finished.\n',mfilename);

return;

