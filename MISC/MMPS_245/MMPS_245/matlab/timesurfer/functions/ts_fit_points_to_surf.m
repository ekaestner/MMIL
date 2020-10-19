function T = ts_fit_points_to_surf(hptsfile,surffile,varargin)
%function T = ts_fit_points_to_surf(hptsfile,surffile,[options])
%
% Purpose: find a rigid body transformation (rotation+translation)
%   for a set of points to fit them as closely as possible to a surface
%
% Required Input:
%  hptsfile: text file containing head point coordinates
%  surffile: scalp surface file in "tri" format
%
% Optional Parameters:
%  'intransfile': file name for input 4x4 transformation matrix
%    if empty, use T_init
%    {default = []}
%  'T_init': initial 4x4 transformation matrix
%    if empty, use identity matrix
%    {default = []}
%  'outtransfile': file name for output 4x4 transformation matrix
%     if empty, no output file created
%    {default = []}
%  'badpoints': vector of head point indices to ignore
%      (including cardinal points)
%    {default = []}
%  'trans_type': ['head2mri' or 'mri2head']
%     if 'head2mri': input and output transformation matrices are head2mri
%       i.e. head  points registered to surface (derived from MRI)
%     if 'mri2head': input and output transformation matrices are mri2head
%    {default = 'mri2head'}
%  'oldtrans_flag': [0|1] whether intransfile is old format
%    with mm units for translations
%    {default = 0}
%  'verbose': [0|1] whether to output messages to terminal
%    {default = 0} 
%
% Output:
%   T: 4x4 transformation matrix
%      head2mri or mri2head depending on 'trans_type'
%
% Created:  11/06/11 by Eran Mukamel  emukamel@ucsd.edu
% Last Mod: 11/07/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
parms = mmil_args2parms( varargin, {...
  'intransfile',[],[],...
  'T_init',[],[],...
  'outtransfile',[],[],...
  'badpoints',[],[],...
  'trans_type','mri2head',{'head2mri','mri2head'},...
  'oldtrans_flag',false,[false true],...
  'verbose',false,[false true],...
...
  'niters',0,[0,Inf],...
});
T = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(hptsfile)
  error('hptsfile %s not found',hptsfile);
end;
if ~exist(surffile)
  error('surffile %s not found',surffile);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize transformation
T0 = init_trans(parms);

% load head points
points = importdata(hptsfile);
points_coords = points.data(:,2:4);

% exclude outliers
if ~isempty(parms.badpoints)
  points_coords = points_coords(setdiff([1:size(points_coords,1)],parms.badpoints),:);
end

% read surface
surf = fs_read_trisurf(surffile);
surf_verts = surf.vertices/1000; % convert to meters

% initialize transformation
X0 = T0(1:3,4)'; % displacement
Ang0 = orth2ang(T0(1:3,1:3)); % Euler angles of rotation
XAng0 = [X0, Ang0];
lb = [-0.2,-0.2,-0.2, -pi,-pi,-pi]; % lower bounds
ub = -lb; % upper bounds
opts = optimset('Display','off','Algorithm','active-set');

dist0 = pointdist(points_coords, surf_verts, XAng0);
if parms.verbose
  fprintf('%s: initial RMSd = %3.3g.\n',mfilename,dist0);
  fprintf('%s: initial transformation matrix:\n',mfilename);
  disp(T0);
end;

% find best fitting transformation
fprintf('%s: finding best fit...\n',mfilename);
[XAng,finaldist,exitflag] =...
  fmincon(@(XAng)pointdist(points_coords,surf_verts,XAng),XAng0,[],[],[],[],lb,ub,[],opts);

% return the final transformation matrix
T_head2mri = eye(4);
T_head2mri(1:3,4) = XAng(1:3);
T_head2mri(1:3,1:3) = ang2orth(XAng(4:6));
T = T_head2mri;
if strcmp(parms.trans_type,'mri2head')
  T = inv(T);
end;
if ~isempty(parms.outtransfile)
  ts_write_transfile(parms.outtransfile,T);
end;

if parms.verbose
  fprintf('%s: final RMSd = %3.3g.\n',mfilename,finaldist);
  fprintf('%s: final transformation matrix:\n',mfilename);
  disp(T);
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T0 = init_trans(parms)
  T0 = [];
  % load trans file
  if parms.verbose
    fprintf('%s: initializing transformation matrix...\n',mfilename);
  end;
  if ~isempty(parms.intransfile)
    if ~exist(parms.intransfile,'file')
      error('intransfile %s not found',parms.intransfile);
    end;
    if parms.verbose
      fprintf('%s: reading transfile %s...\n',mfilename,parms.intransfile);
    end;
    if parms.oldtrans_flag
      T0 = ts_read_transfile(parms.intransfile,0.001);
    else
      T0 = ts_read_transfile(parms.intransfile);
    end;
  elseif ~isempty(parms.T_init)
    if parms.verbose
      fprintf('%s: using T_init...\n',mfilename);
    end;
    T0 = parms.T_init;
  else
    if parms.verbose
      fprintf('%s: using identity matrix...\n',mfilename);
    end;
    T0 = eye(4);
  end;
  % internal trans_type is always head2mri (applied to head points)
  if strcmp(parms.trans_type,'mri2head')
    T0 = inv(T0);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang = orth2ang(orthm)
  ang(1) = asin(orthm(1,3));
  ang(2) = angle( orthm(1,1:2)*[1 ;i] );
  yz = orthm* ...
      [orthm(1,:)',...
      [-sin(ang(2)); cos(ang(2)); 0],...
      [-sin(ang(1))*cos(ang(2)); -sin(ang(1)*sin(ang(2)));
      cos(ang(1))] ];
  ang(3) = angle(yz(2,2:3)* [1; i]);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function orthm = ang2orth(ang)
  sa = sin(ang(2)); ca = cos(ang(2));
  sb = sin(ang(1)); cb = cos(ang(1));
  sc = sin(ang(3)); cc = cos(ang(3));
  ra = [  ca,  sa,  0; ...
         -sa,  ca,  0; ...
           0,   0,  1];
  rb = [  cb,  0,  sb; ...
           0,  1,  0; ...
         -sb,  0,  cb];
  rc = [  1,   0,   0; ...
          0,   cc, sc;...
          0,  -sc, cc];
  orthm = rc*rb*ra;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dist = pointdist(points,surf_verts,XAng)
  T = eye(4);
  T(1:3,4) = XAng(1:3);
  T(1:3,1:3) = ang2orth(XAng(4:6));

  npoints = size(points,1);
  nvert = size(surf_verts,1);

  % Transform the fiducial points
  points = [points, ones(npoints,1)] * T';

  % Find the minimum distance to the surface
  dist = 0;
  for jxyz = 1:3
    dist = dist + (points(:,jxyz)*ones(1,nvert) - ones(npoints,1)*surf_verts(:,jxyz)').^2;
  end
  % dist = sum(min(dist, [], 2));
  dist = sum(sqrt(min(dist, [], 2)));
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
