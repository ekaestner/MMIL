function suffix = dti_set_fiber_suffix(varargin)
%function suffix = dti_set_fiber_uffix(varargin)
%
% Purpose: set fiber suffix depending on parameters
%
% Optional Parameters:
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
% Created:   02/23/13 by Don Hagler
% Last Mod:  02/23/13 by Don Hagler

suffix = [];

parms = mmil_args2parms(varargin,{...
  'xcg_flag',false,[false true],...
  'xcg_suffix','xcg',[],...
  'masksf_flag',false,[false true],...
  'masksf_suffix','masksf',[],...
  'disp_flag',false,[false true],...
  'disp_suffix','dwtd',[],...
  'dispfact',4,[1e-6,1e6],...
  'thresh_prob',0,[0 1],...
  'thresh_FA',0,[0 1],...
});

% build suffix according to parameters
if parms.xcg_flag
  suffix = concat_suffix(suffix,parms.xcg_suffix);
end;
if parms.masksf_flag
  suffix = concat_suffix(suffix,parms.masksf_suffix);
end;
if parms.disp_flag
  suffix = concat_suffix(suffix,sprintf('%s%0.1f',...
    parms.disp_suffix,parms.dispfact));
end;
if parms.thresh_prob>0
  suffix = concat_suffix(suffix,sprintf('pthresh%0.2f',parms.thresh_prob));
end;
if parms.thresh_FA>0
  suffix = concat_suffix(suffix,sprintf('FAthresh%0.2f',parms.thresh_FA));
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function suffix = concat_suffix(prefix,suffix)
  if ~isempty(prefix)
    suffix = [prefix '_' suffix];
  end;
return;

