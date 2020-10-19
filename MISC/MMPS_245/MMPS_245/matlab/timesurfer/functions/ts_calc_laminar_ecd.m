function data_ecd = ts_calc_laminar_ecd(data_csd,varargin)
%function data_ecd = ts_calc_laminar_ecd(data_csd,[options])
%
% Purpose: calculate equivalent current dipole from
%                    current source density
%
% Required Input:
%   data_csd: current source density
%     matrix with size = [nchans+1,ntpoints]
%     units = uA/mm^3
%
% Optional Input: ('key',value pair)
%   'dz': intrachannel spacing (mm)
%     {default = 0.15}
%
% Output:
%   data_ecd: equivalent current dipole
%     vector with size = [1,ntpoints]
%     units = uA*mm/mm^2
%
% Created:  08/17/16 by Don Hagler
% Last Mod: 08/17/16 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_ecd = [];
if ~mmil_check_nargs(nargin,1), return; end;

% parse input parameters
[parms,data_csd] = check_input(data_csd,varargin);

data_ecd = calc_ecd(data_csd,parms);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [parms,data_csd] = check_input(data_csd,options)
  parms = mmil_args2parms(options,{...
    'dz',0.15,[],...
  });

  % check data matrix
  if ndims(data_csd)~=2
    error('data matrix must be 2D (nchans x ntpoints)');
  end;
  [parms.nchans,parms.ntpoints] = size(data_csd);
  parms.chans = [1:parms.nchans];
  % distance along the electrode in mm
  parms.z = (parms.chans-1)*parms.dz; 
  ind_mid = round(parms.nchans/2);
  parms.z_cent = parms.z - parms.z(ind_mid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data_ecd = calc_ecd(data_csd,parms)
  % initialize output
  data_ecd = zeros(1,parms.ntpoints);
  % scale CSD by distance from center
  z_data_csd = bsxfun(@times,parms.z_cent',data_csd);
  % evaluate ECD using trapezoidal integration
  data_ecd = trapz(parms.z_cent,z_data_csd,1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


