function [fiber_infix,fiber_ext] = dti_set_fiber_infix(varargin)
%function [fiber_infix,fiber_ext] = dti_set_fiber_infix(varargin)
%
% Purpose: set fiber infix and ext depending on parameters
%
% Optional Parameters:
%  'atlas_flag': whether to use atlas fibers and if so, what type of atlas fibers
%     0 - manually assisted fiber tracts generated with DTIStudio
%       (DTIStudio_fiber_masks directory must exist
%        -- imported with Import_DTIStudio_FiberMasks)
%     1 - location-only "count" atlas tracks
%     2 - location+direction "count" atlas tracks
%     3 - location-only "mask" atlas tracks
%     4 - location+direction "mask" atlas tracks
%       (fiber_masks_from_atlas directory must exist
%        -- created with AtlasTrack_Fibers)
%     {default = 2}
%   'count_flag': [0|1] for atlas_flag=0, whether to use fiber counts (or masks)
%     {default = 1}
%   'resT1flag': [0|1] whether to use fibers resampled to T1 resolution
%     {default = 0}
%   'xcg_flag': [0|1] exclude CSF and gray-mattter from fiber ROIs
%     {default = 0}
%   'xcg_suffix':  suffix attached to output fiber file names
%       after excluding CSF and gray matter
%     {default = 'xcg'}
%   'masksf_flag': [0|1] exclude voxels with multiple fibers
%     {default = 0}
%   'masksf_suffix': suffix attached to output fiber file names
%     after excluding voxels with multiple fibers
%     {default = 'masksf'}
%   'disp_flag': [0|1] generate fiber ROIs excluding CSF and gray matter
%     as defined by FreeSurfer aseg
%     {default = 0}
%   'disp_suffix':  suffix attached to output fiber file names
%       after calculating dispersion weighting
%     {default = 'dwtd'}
%   'dispfact': multiplicative factor applied to dispersion values
%     {default = 4}
%   'thresh_prob': fiber probability threshold applied to fiber ROIs
%     {default = 0}
%   'thresh_FA': fractional anisotropy threshold applied to fiber ROIs
%     {default = 0}
%
% Created:   11/10/09 by Don Hagler
% Last Mod:  02/23/13 by Don Hagler

fiber_infix = []; fiber_ext = [];

parms = mmil_args2parms(varargin,{...
  'atlas_flag',2,[0:4],...
  'count_flag',true,[false true],...
  'resT1flag',false,[false true],...
  'xcg_flag',false,[false true],...
  'xcg_suffix','xcg',[],...
  'masksf_flag',false,[false true],...
  'masksf_suffix','masksf',[],...
  'disp_flag',false,[false true],...
  'disp_suffix','dwtd',[],...
  'dispfact',4,[1e-6,1e6],...
  'thresh_prob',0,[0 1],...
  'thresh_FA',0,[0 1],...
...
  'suffix_tags',{'xcg_flag','xcg_suffix',...
                  'masksf_flag','masksf_suffix',...
                  'disp_flag','disp_suffix','dispfact',...
                  'thresh_prob','thresh_FA'},[],...
});

% set fiber_infix and fiber_ext
switch parms.atlas_flag
  case 0 % manual
    if parms.count_flag
      fiber_infix = 'count';
    else
      fiber_infix = 'mask';
    end;
    fiber_ext = '.mat';
  case 1 % loc only, count atlas
    fiber_infix = 'countatlas';
    fiber_ext = '.mat';
  case 2 % loc+dir, count atlas
    fiber_infix = 'prob_countatlas';
    fiber_ext = '.mat';
  case 3 % loc only, mask atlas
    fiber_infix = 'maskatlas';
    fiber_ext = '.mat';
  case 4 % loc+dir, mask atlas
    fiber_infix = 'prob_maskatlas';
    fiber_ext = '.mat';
end;

if parms.resT1flag
  fiber_infix = [fiber_infix '_resT1'];
end;

args = mmil_parms2args(parms,parms.suffix_tags);
suffix = dti_set_fiber_suffix(args{:});
if ~isempty(suffix)
  fiber_infix = [fiber_infix '_' suffix];
end;

