function [trans,rot] = ts_calc_fif_motion(fnames,varargin)
%function [trans,rot] = ts_calc_fif_motion(fnames,[options])
%
% Purpose: calculate head motion between multiple fif files
%   using device2head transformation in each
%
% Required Input:
%   fnames: cell array of fif file names
%
% Optional Parameters ('key',value pairs):
%  'plot_flag': [0|1] whether to plot results
%    {default = 1}
%  'plot_outstem': output file stem for plots
%    supply full path, or else relative to current directory
%    {default = 'motion'}
%  'visible_flag': [0|1] whether plots should be visible
%    {default = 1}
%
% Output:
%   trans: matrix of head translations (mm)
%     size is N x 3, with N = number of fif files
%     and one column each for x, y, and z
%   rot: matrix of head rotations (degrees)
%
% Created:  08/10/11 by Don Hagler
% Last Mod: 11/14/14 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'plot_flag',true,[false true],...
  'visible_flag',true,[false true],...
  'plot_outstem','motion',[],...
...
  'colstrs',{'ro-','gx-','b*-'},[],...
  'legstrs',{'x','y','z'},[],...
});
trans = [];
rot = [];

N_axes = 3;

if ~iscell(fnames), fnames = {fnames}; end;
N = length(fnames);

if mmil_isrelative(parms.plot_outstem)
  parms.plot_outstem = [pwd '/' parms.plot_outstem];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_device2head_array = cell(N,1);
for f=1:N
  if ~exist(fnames{f},'file')
    error('file %s not found',fnames{f});
  end;
  M_device2head_array{f} = loadtrans(fnames{f});
end;

M0 = M_device2head_array{1};
M_motion_array = cell(N,1);
trans = zeros(N,N_axes);
rot = zeros(N,N_axes);
for f=2:N
  M = M_device2head_array{f};
  M_motion_array{f} = M*inv(M0);
  motion_vec = mmil_M_mat2vec(M_motion_array{f});
  trans(f,:) = motion_vec(1:3)*1000;
  rot(f,:) = motion_vec(4:6);
end;

if parms.plot_flag
  figure(1); clf; hold on;
  for i=1:N_axes
    plot(1:N,trans(:,i),parms.colstrs{i});
  end;
  legend(parms.legstrs);
  title('Subject Head Motion (translation)');
  xlabel('time point');
  ylabel('translation (mm)');
  xlim([1,N]);
  ylim([-10,10]);
  if ~parms.visible_flag
    set(gcf,'Visible','Off');
  end;
  fname_out = [parms.plot_outstem '_trans.tif'];
  print(gcf,'-dtiff',fname_out);
  if ~parms.visible_flag
    close(gcf);
  end;

  figure(2); clf; hold on;
  for i=1:N_axes
    plot(1:N,rot(:,i),parms.colstrs{i});
  end;
  legend(parms.legstrs);
  title('Subject Head Motion (rotation)');
  xlabel('time point');
  ylabel('rotation (degrees)');
  xlim([1,N]);
  ylim([-10,10]);
  if ~parms.visible_flag
    set(gcf,'Visible','Off');
  end;
  fname_out = [parms.plot_outstem '_rot.tif'];
  print(gcf,'-dtiff',fname_out);
  if ~parms.visible_flag
    close(gcf);
  end;
end;

fprintf('\n');
fprintf('time point    trans x    trans y    trans z    rot x     rot y    rot z\n');
for f=1:N
  fprintf('%4d          %6.3f     %6.3f     %6.3f     %6.3f   %6.3f   %6.3f\n',...
    f,trans(f,1),trans(f,2),trans(f,3),rot(f,1),rot(f,2),rot(f,3));
end;
fprintf('\n');

