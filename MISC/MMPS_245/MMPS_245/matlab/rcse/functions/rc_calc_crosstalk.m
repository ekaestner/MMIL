function CT = rc_calc_crosstalk(varargin)
%function CT = rc_calc_crosstalk([options])
%
% Purpose: calculate crosstalk between RCSE sources
%
% Usage:
%  CT = rc_calc_crosstalk('key1', value1,...);
%
% Optional Inupt:
%  'prefix': prefix of RCSE output files
%    {default = 'RCSE'}
%  'rootdir': directory containing matfiles dir
%    {default = pwd}
%  'plot_flag': [0|1] whether to plot color image of matrix
%    {default = 1}
%  'outstem': output file stem for plot
%    if empty, use prefix
%    full path or relative to rootdir
%    {default = []}
%  'outfix': string attached to outstem
%    {default = 'crosstalk'}
%  'visible_flag': [0|1] whether to make plot visible
%    {default = 0}
%  'clim': 2-element vector with lower and upper value limits for imagesc
%    {default = [-1,1]}
%
% Output:
%   CT: crosstalk matrix
%
% Created:  08/19/12 by Don Hagler
% Last Mod: 01/23/13 by Don Hagler
%

%% todo: for indy solutions, sum across dipoles for each area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if ~mmil_check_nargs(nargin,1), return; end;
CT = [];
parms = mmil_args2parms(varargin, { ...
  'prefix','RCSE',[],...
  'rootdir',pwd,[],...
  'plot_flag',true,[false true],...
  'outstem',[],[],...
  'outfix','crosstalk',[],...
  'visible_flag',false,[false true],...
  'clim',[0,1],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load RCSE files
matname=sprintf('%s/matfiles/%s_inverse.mat',parms.rootdir,parms.prefix);
load(matname);
matname=sprintf('%s/matfiles/%s_ret_forward.mat',parms.rootdir,parms.prefix);
load(matname);

% calculate crosstalk

CT = (inverse.W*retforward.F).^2;
for i=1:size(CT,1)
  CT(i,:) = CT(i,:) / CT(i,i);
end;

% plot crosstalk
if parms.plot_flag
  if isempty(parms.outstem)
    parms.outstem = parms.prefix;
  end;
  if mmil_isrelative(parms.outstem)
    parms.outstem = [rootdir '/' parms.outstem];
  end;
  fname = sprintf('%s-%s.tif',parms.outstem,parms.outfix);
  figure;
  imagesc(CT,parms.clim);
  axis off;
  colorbar;
  if ~parms.visible_flag, set(gcf,'visible','off'); end;
  print('-dtiff',fname);
  if ~parms.visible_flag, close(gcf); end;
end;

