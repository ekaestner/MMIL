function motion_data = mmil_load_motion_par(fname_motion,varargin)
%function motion_data = mmil_load_motion_par(fname_motion,varargin)
%
% Purpose: load mcf.par file from FSL mcflirt
%
% Usage:
%  motion_data = mmil_load_motion(fname_motion,'key1', value1,...); 
%
% Required Parameters:
%   fname_motion: text file containing 6 columns of motion estimates
%
% Optional parameters:
%  'skipTRs': number of initial repetitions to remove
%    {default: 0}
%  'nframes': number of frames in corresponding data file
%    If supplied, will check that number of time points in fname_motion matches
%    {default: []}
%
% Output:
%   motion_data: matrix of motion parameter estimates with size [nframes,npar]
%     order of parameters: dx, dy, dz, rx, ry, rz
%     with units of mm and degrees
%
% Created:  03/23/16 Don Hagler
% Last Mod: 03/23/16 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

motion_data = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'skipTRs',0,[0 Inf],...
  'nframes',[],[],...
...
  'reorder_vec',[4,5,6,1,2,3],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  load fname_motion
if ~exist(fname_motion,'file')
  error('file %s not found',fname_motion);
else
  motion_data = load(fname_motion);
  if isempty(parms.nframes)
    parms.nframes = size(motion_data,1);
  end;
  if size(motion_data,1) ~= parms.nframes
    error('number of time points (%d) in fname_motion %s does not match number of frames (%d)',...
      size(motion_data,1),parms.fname_motion,parms.nframes);
  end;
  if parms.skipTRs>0 && parms.skipTRs<parms.nframes
    motion_data = motion_data(parms.skipTRs+1:end,:);
  end;
end;

% reorder motion parameter estimates
%   from rz,rx,ry,dz,dx,dy to dx,dy,dz,rx,ry,rz
motion_data = motion_data(:,parms.reorder_vec);

% convert from radians to degrees    
motion_data(:,4:6) = 180 * motion_data(:,4:6) / pi;

