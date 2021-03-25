function M = mmil_construct_M(varargin)
%function M = mmil_construct_M(varargin)
%
% Purpose: create a 4x4 transformation matrix
%   with translation, rotation, and scaling
%
% Usage: mmil_construct_M('key1', value1,...);
%
% Optional Input:
%  'scale': three value vector of x, y, z scaling factors
%    {default: [1,1,1]}
%  'rot': three value vector of x, y, z rotation (in degrees clockwise)
%    {default: [0,0,0]}
%  'trans': three value vector of x, y, z translation
%    {default: [0,0,0]}
%  'order': order of operations
%    {default: 'scale','rot','trans'}
%  'rot_order': order of rotations
%    {default: {'x','y','z'}}
%  'orient': 3 character orientation string e.g. 'RAS','LPI', etc.
%    overrides 'rot' and 'rot_order' options and flips scale(1) as necessary
%    {default: []}
%
% Output:
%   4x4 transformation matrix
%
% created:  01/12/10 by Don Hagler
% last mod: 0/25/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

%if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'scale',[1,1,1],[],...
  'rot',[0 0 0],[],...
  'trans',[0 0 0],[],...
  'order',{'scale','rot','trans'},{'scale','rot','trans'},...
  'rot_order',{'x','y','z'},{'x','y','z'},...
  'orient',[],[],...
...
  'smf',10*eps,[],...
});

% make sure they are three-value vectors
parms.scale(length(parms.scale)+1:3) = 0;
parms.scale = parms.scale(1:3);
parms.rot(length(parms.rot)+1:3) = 0;
parms.rot = parms.rot(1:3);
parms.trans(length(parms.trans)+1:3) = 0;
parms.trans = parms.trans(1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(parms.orient)
  parms.rot_order = {'x','y','z'};
  switch upper(parms.orient)
    case 'AIL'
      parms.rot = [ 90   0 -90 ];
      f = 1;
    case 'AIR'
      parms.rot = [ 90   0  90 ];
      f = -1;
    case 'ALI'
      parms.rot = [  0 180 -90 ];
      f = -1;
    case 'ALS'
      parms.rot = [  0   0 -90 ];
      f = 1;
    case 'ARI'
      parms.rot = [  0 180  90 ];
      f = 1;
    case 'ARS'
      parms.rot = [  0   0  90 ];
      f = -1;
    case 'ASL'
      parms.rot = [ 90 180 -90 ];
      f = -1;
    case 'ASR'
      parms.rot = [ 90 180  90 ];
      f = 1;
    case 'IAL'
      parms.rot = [  0  90   0 ];
      f = -1;
    case 'IAR'
      parms.rot = [  0 -90   0 ];
      f = 1;
    case 'ILA'
      parms.rot = [  0 -90 -90 ];
      f = 1;
    case 'ILP'
      parms.rot = [  0  90 -90 ];
      f = -1;
    case 'IPL'
      parms.rot = [  0 -90 180 ];
      f = 1;
    case 'IPR'
      parms.rot = [  0  90 180 ];
      f = -1;
    case 'IRA'
      parms.rot = [  0  90  90 ];
      f = -1;
    case 'IRP'
      parms.rot = [  0 -90  90 ];
      f = 1;
    case 'LAI'
      parms.rot = [  0 180   0 ];
      f = 1;
    case 'LAS'
      parms.rot = [  0   0   0 ];
      f = -1;
    case 'LIA'
      parms.rot = [ 90   0   0 ];
      f = -1;
    case 'LIP'
      parms.rot = [ 90   0 180 ];
      f = 1;
    case 'LPI'
      parms.rot = [  0 180 180 ];
      f = -1;
    case 'LPS'
      parms.rot = [  0   0 180 ];
      f = 1;
    case 'LSA'
      parms.rot = [ 90 180   0 ];
      f = 1;
    case 'LSP'
      parms.rot = [ 90 180 180 ];
      f = -1;
    case 'PIL'
      parms.rot = [ 90   0 -90 ];
      f = -1;
    case 'PIR'
      parms.rot = [ 90   0  90 ];
      f = 1;
    case 'PLI'
      parms.rot = [  0 180 -90 ];
      f = 1;
    case 'PLS'
      parms.rot = [  0   0 -90 ];
      f = -1;
    case 'PRI'
      parms.rot = [  0 180  90 ];
      f = -1;
    case 'PRS'
      parms.rot = [  0   0  90 ];
      f = 1;
    case 'PSL'
      parms.rot = [ 90 180 -90 ];
      f = 1;
    case 'PSR'
      parms.rot = [ 90 180  90 ];
      f = -1;
    case 'RAI'
      parms.rot = [  0 180   0 ];
      f = -1;
    case 'RAS'
      parms.rot = [  0   0   0 ];
      f = 1;
    case 'RIA'
      parms.rot = [ 90   0   0 ];
      f = 1;
    case 'RIP'
      parms.rot = [ 90   0 180 ];
      f = -1;
    case 'RPI'
      parms.rot = [  0 180 180 ];
      f = 1;
    case 'RPS'
      parms.rot = [  0   0 180 ];
      f = -1;
    case 'RSA'
      parms.rot = [ 90 180   0 ];
      f = -1;
    case 'RSP'
      parms.rot = [ 90 180 180 ];
      f = 1;
    case 'SAL'
      parms.rot = [  0   90  0 ];
      f = 1;
    case 'SAR'
      parms.rot = [  0  -90  0 ];
      f = -1;
    case 'SLA'
      parms.rot = [  0  -90 -90 ];
      f = -1;
    case 'SLP'
      parms.rot = [  0  90 -90 ];
      f = 1;
    case 'SPL'
      parms.rot = [  0 -90 180 ];
      f = -1;
    case 'SPR'
      parms.rot = [  0  90 180 ];
      f = 1;
    case 'SRA'
      parms.rot = [  0  90  90 ];
      f = 1;
    case 'SRP'
      parms.rot = [  0 -90  90 ];
      f = -1;
    otherwise
      error('invalid orientation')
  end;
  parms.scale = parms.scale.*[f 1 1];
end;

parms.rot = parms.rot*pi/180; % convert to radians

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = eye(4);
for i=1:length(parms.order)
  switch parms.order{i}
    case 'scale'
      M = M_scale(parms.scale)*M;
    case 'rot'
      for j=1:length(parms.rot_order)
        switch parms.rot_order{j}
          case 'x'
            M = M_rotx(parms.rot(1))*M;
          case 'y'
            M = M_roty(parms.rot(2))*M;
          case 'z'
            M = M_rotz(parms.rot(3))*M;
        end;
      end;
    case 'trans'
      M = M_trans(parms.trans)*M;
  end;
end;

if ~isempty(parms.orient)
  M(abs(M)<parms.smf) = 0;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = M_rotx(a)
  sa = sin(a); ca = cos(a);
  M = [1 0 0 0; 0 ca sa 0; 0 -sa ca 0; 0 0 0 1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = M_roty(a)
  sa = sin(a); ca = cos(a);
  M = [ca 0 -sa 0; 0 1 0 0; sa 0 ca 0; 0 0 0 1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = M_rotz(a)
  sa = sin(a); ca = cos(a);
  M = [ca sa 0 0; -sa ca 0 0; 0 0 1 0; 0 0 0 1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = M_scale(s)
  M = [s(1) 0 0 0; 0 s(2) 0 0; 0 0 s(3) 0; 0 0 0 1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = M_trans(t)
  M = [1 0 0 t(1); 0 1 0 t(2); 0 0 1 t(3); 0 0 0 1];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
