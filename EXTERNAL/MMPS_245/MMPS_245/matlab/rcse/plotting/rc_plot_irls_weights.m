function rc_plot_irls_weights(prefix,varargin)
%function rc_plot_irls_weights(prefix,[options])
%
% Purpose: plot weights vs. error from reweighted least-squares
%
% Optional ('key',value,...)
%   'norm_flag': [0|1] use error, minimum error subtracted,
%     normalized by median absolute deviation
%     {default = 0}
%   'tif_flag': [0|1] save plots as tifs
%     {default = 0}
%
% Created:  07/13/10 by Don Hagler
% Last Mod: 02/20/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
orig_parms = mmil_args2parms(varargin, { ...
  'norm_flag',false,[false true],...
  'tif_flag',false,[false true],...
... % todo: use these
  'label_flag',true,[false true],...
  'fontname','Arial',[],...
  'fontsize',12,[],...
});

fprintf('%s: plotting IRLS weights and errors for %s\n',mfilename,prefix);
matname=sprintf('matfiles/%s_parms.mat',prefix);
load(matname);
matname=sprintf('matfiles/%s_results.mat',prefix);
load(matname);

if orig_parms.norm_flag
  err = results.cond_err_norm;
else
  err = results.cond_err;
end;
weights = results.cond_weights;

% error historgram
if orig_parms.tif_flag
  figure(1); clf;
  set(gcf,'Visible','off');
else
  figure;
end;
hist(err)
if orig_parms.tif_flag
  fname = sprintf('%s-irls_weights_histo.tif',prefix);
  print('-dtiff',fname);
end;  
  
% weights vs. error
if orig_parms.tif_flag
  figure(1); clf;
  set(gcf,'Visible','off');
else
  figure;
end;
plot(err,weights,'bo');
set(gca,'YLim',[0,1]);
title('weights vs. error');
if orig_parms.tif_flag
  fname = sprintf('%s-irls_weights_vs_error.tif',prefix);
  print('-dtiff',fname);
  close(gcf);
end;  




