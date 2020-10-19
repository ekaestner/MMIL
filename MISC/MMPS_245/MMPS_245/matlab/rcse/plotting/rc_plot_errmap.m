function output = rc_plot_errmap(prefix,varargin)
%function output = rc_plot_errmap(prefix,[options])
%
% Optional Parameters:
%   'rootdir': MEG session analysis directory
%     {default = pwd}
%   'sign_flag': [0|1|2] whether to use abs, signed, or squared error
%     0: absolute error
%     1: signed error
%     2: squared error
%     {default = 0}
%   'std_flag': [0|1] calculate standard deviation of error across sensors
%     otherwise calculate mean
%     {default = 0} 
%   'irls_err_flag': [0|1] whether to plot normalized error calculatd
%      during iteratively reweighted least-squares (IRLS)
%      otherwise, calculate error from results.E
%     {default = 0}
%   'irls_weights_flag': [0|1] whether to plot weights from IRLS instead of err
%     {default = 0}
%   'irls_init_flag': [0|1] whether IRLS err or weights should be from initial
%     iteration (1) or last (0)
%     {default = 0}
%   'tif_flag': [0|1] save plots as tifs
%     {default = 0}
%   'eps_flag': [0|1] save plots as eps
%     {default = 0}
%   't0': start time of range used to calculate error
%     {default = -Inf}
%   't1': end time of range used to calculate error
%     {default = Inf}
%
% Output:
%   output: struct with fields:
%     Z: image of normalized error
%     W: image of unnormalized error
%     M: mask of stimulus locations
%     err_norm: vector of normalized error for each stimulus location
%     err_raw: vector of unnormalized error
%     weight: vector of IRLS weights (empty if irls_weights_flag = 0)
%     max_std: greatest standard deviation time across channels and locations
%     max_std: greatest variance deviation time across channels and locations
%
% Created:  06/21/10 by Don Hagler
% Last Mod: 08/01/13 by Don Hagler
%

output = [];
if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'rootdir',pwd,[],...
  'sign_flag',0,[0 1 2],...
  'std_flag',false,[false true],...
  'irls_err_flag',false,[false true],...
  'irls_weights_flag',false,[false true],...
  'irls_init_flag',false,[false true],...
  'tif_flag',false,[false true],...
  'eps_flag',false,[false true],...
  't0',-Inf,[],...
  't1',Inf,[],...
...
  'visible_flag',true,[false true],...
  'outdir',pwd,[],...
  'outstem',prefix,[],...
  'cmap','jet',[],...
  'label_flag',true,[false true],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
...
  'grid_unit',0.05,[],...
  'clim',[0,1],[],...
  'z_offset',0,[],...
  'err_xlim',[],[],...
  'cbar_flag',true,[false true],...
...
  'fg_color',0.4,[],...
  'bg_color',0,[],...
...
  'plot_flag',true,[false true],...  
});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(parms.fg_color)==1
  parms.fg_color = parms.fg_color*ones(1,3);
end;
if length(parms.bg_color)==1
  parms.bg_color = parms.bg_color*ones(1,3);
end;

parms.outstem = [parms.outdir '/' parms.outstem];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fprintf('%s: plotting fit variance for %s\n',mfilename,prefix);
matname=sprintf('%s/matfiles/%s_parms.mat',parms.rootdir,prefix);
ri = load(matname);
matname=sprintf('%s/matfiles/%s_results.mat',parms.rootdir,prefix);
load(matname);

time = results.time*1000;
[~,parms.s0] = min(abs(time - parms.t0));
if parms.t1 > time(end)
  parms.s1 = length(time);
else
  [~,parms.s1] = min(abs(time - parms.t1));
end;

cond_info = ri.parms.cond_info;
unique_location_conds = ri.parms.unique_location_conds;

num_sensors = length(ri.parms.goodchans);

max_std = max(std(results.Y,0,2));
max_var = max(var(results.Y,0,2));

r_max = max([cond_info.ecc]) + max([cond_info.ecc_width])/2;

x = [-r_max:parms.grid_unit:r_max];
y = [r_max:-parms.grid_unit:-r_max];

[X,Y] = meshgrid(x,y);
Z = zeros(size(X)); % norm err image
W = zeros(size(X)); % raw err image
M = zeros(size(X)); % mask of stimulus locations

% calculate polar coordinates for grid points
R = sqrt(X.^2 + Y.^2);
T = atan2(Y,X);

c = 1;
errvec = [];
errvec_raw = [];
weightvec = [];
for j=1:length(unique_location_conds);
  i = unique_location_conds(j);
  if cond_info(i).contrast==0, continue; end;

  % calculate average error for this condition
  p = 1+(c-1)*num_sensors;
  q = p+num_sensors-1;
  err = results.E(parms.s0:parms.s1,p:q);
  data = results.Y(parms.s0:parms.s1,p:q);
  switch parms.sign_flag
    case 0 % abs
      err = abs(err);
    case 1 % signed
      if ~parms.std_flag
        pol = ones(size(err,1),1)*sign(mean(data,1));
        err = err .* pol;
      end;
    case 2 % squared
      err = err.^2;
  end;  
  if parms.std_flag
    err_raw = mean(std(err,0,2));
  else
    err_raw = mean(err(:));
  end;

  if parms.irls_weights_flag || parms.irls_err_flag
    if parms.irls_init_flag
      weight = results.init_cond_weights(c);
    else
      weight = results.cond_weights(c);
    end;
    weightvec = [weightvec weight];
  end;

  % calculate average normalized error
  if parms.irls_weights_flag
    err_norm = 1-weight;
  elseif parms.irls_err_flag
    if parms.irls_init_flag
      err_norm = results.init_cond_err_norm(c);
    else
      err_norm = results.cond_err_norm(c);
    end;
  elseif parms.sign_flag == 2
    err_norm = err_raw/max_var;
  else
    err_norm = err_raw/max_std;
  end;

  r = cond_info(i).ecc;
  t = cond_info(i).theta*pi/180;
  if t>pi, t = t - 2*pi; end;
  
  dr = cond_info(i).ecc_width/2;
  dt = cond_info(i).theta_width*pi/360;

  ind = find(...
    (R>r-dr & R<r+dr) & ...
    (T>t-dt & T<t+dt));
    
  Z(ind) = err_norm;
  W(ind) = err_raw;
  M(ind) = 1;

  errvec = [errvec err_norm];
  errvec_raw = [errvec_raw err_raw];

  c = c + 1;
end

if nargout
  output.Z = Z;
  output.W = W;
  output.M = M;
  output.err_norm = errvec;
  output.err_raw = errvec_raw;
  output.weight = weightvec;
  output.max_std = max_std;
  output.max_var = max_var;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~parms.plot_flag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% err map

if parms.visible_flag, figure(1); end;

tags = {'outdir','outstem','tif_flag','eps_flag','visible_flag',...
  'cbar_flag','clim','cmap','title','fontname','fontsize','offset'};
tmp_parms = parms;
if parms.z_offset
  Z(M>0) = Z(M>0) + parms.z_offset;
end;
tmp_parms.offset = parms.z_offset;
if parms.label_flag
  if parms.irls_weights_flag
    tmp_parms.title = '1 - weights from IRLS';
  elseif parms.irls_err_flag
    tmp_parms.title = 'normalized error from IRLS';
  else
    tmp_parms.title = 'mean sq err / max var';
  end;
else
  tmp_parms.title = [];
end;
tmp_parms.outstem = sprintf('%s-errmap',parms.outstem);
args = mmil_parms2args(tmp_parms,tags);
mmil_imagesc(Z,args{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% histogram

if ~parms.visible_flag
  close(gcf);
end;

if ~parms.visible_flag
  figure();
  set(gcf,'Visible','off');
else
  figure(2); clf;
end;

hist(errvec,10);
tmp = get(gca,'XLim');
if isempty(parms.err_xlim)
  set(gca,'XLim',[0,tmp(2)]);
else
  set(gca,'XLim',parms.err_xlim);
end;

h = findobj(gca,'Type','patch');
set(h,'FaceColor',parms.fg_color,'EdgeColor',parms.bg_color)

if parms.label_flag
  if parms.irls_weights_flag
    title_str = 'histogram of IRLS weights';
  elseif parms.irls_err_flag
    title_str = 'histogram of IRLS normalized error';
  elseif parms.std_flag
    switch parms.sign_flag
      case 0
        title_str = 'histogram of std abs(err) / max std';
      case 1
        title_str = 'histogram of std err / max std';
      case 2
        title_str = 'histogram of std err^2 / max var';
    end;  
  else
    switch parms.sign_flag
      case 0
        title_str = 'histogram of mean abs(err) / max std';
      case 1
        title_str = 'histogram of mean err / max std';
      case 2
        title_str = 'histogram of mean err^2 / max var';
    end;  
  end;
  title(title_str,'FontSize',parms.fontsize,'FontName',parms.fontname);

  if parms.irls_weights_flag
    xlabel('IRLS weight',...
      'FontSize',parms.fontsize,'FontName',parms.fontname);
  else
    xlabel('normalized error',...
      'FontSize',parms.fontsize,'FontName',parms.fontname);
  end;
  ylabel('number of stimulus locations',...
    'FontSize',parms.fontsize,'FontName',parms.fontname);
end;

if parms.tif_flag
  print('-dtiff',sprintf('%s-err_hist.tif',parms.outstem));
end;
if parms.eps_flag
  mmil_printeps(gcf,sprintf('%s-err_hist',parms.outstem));
end;

if ~parms.visible_flag
  close(gcf);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weights vs. error

if ~parms.irls_weights_flag, return; end;

if ~parms.visible_flag
  close(gcf);
end;

if ~parms.visible_flag
  figure();
  set(gcf,'Visible','off');
else
  figure(3); clf;
end;

hold on;
plot(errvec,weightvec,'*','color',parms.fg_color);
plot(errvec,weightvec,'o','color',parms.bg_color);
set(gca,'YLim',[0,1]);

if parms.label_flag
  xlabel('normalized error',...
      'FontSize',parms.fontsize,'FontName',parms.fontname);
  ylabel('IRLS weights',...
      'FontSize',parms.fontsize,'FontName',parms.fontname);
  title('weights vs. error',...
      'FontSize',parms.fontsize,'FontName',parms.fontname);
end;

if parms.tif_flag
  print('-dtiff',sprintf('%s-irls_weights.tif',parms.outstem));
end;
if parms.eps_flag
  mmil_printeps(gcf,sprintf('%s-irls_weights',parms.outstem));
end;

if ~parms.visible_flag
  close(gcf);
end;

