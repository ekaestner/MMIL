function fname = rc_plot_normvarerr(prefix,types_flag)
%function fname = rc_plot_normvarerr(prefix,[types_flag])
%
% Purpose: plot normalized variance of RCSE error
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%   types_flag: [0|1] whether to keep grad, mag, and EEG separate
%     {default = 0}
%
% Early Mod: 05/03/09 by Don Hagler
% Last Mod:  08/01/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('types_flag','var') || isempty(types_flag), types_flag = 0; end;
if ~exist('fontname','var') || isempty(fontname), fontname = 'Arial'; end;

fprintf('%s: plotting normalized variance of error for %s\n',mfilename,prefix);

matname=sprintf('matfiles/%s_results.mat',prefix);
load(matname);

clf;
time = results.retmap.areas(1).time*1000;
plot(time,results.norm_var_E,'k','LineWidth',2);
set(gca,'FontName',fontname)
if types_flag
  hold on;
  plot(time,results.norm_var_Eeeg,'r');
  plot(time,results.norm_var_Egrad,'g');
  plot(time,results.norm_var_Emag,'b');
  legend('mean','EEG','grad','mag','Location','SouthEast');
end;

title('normalized variance of error');
axis('tight');
set(gca,'YLim',[0 1]); % normalized residual variance
set(gca,'XGrid','on','XTick',[-100 0 100 200 300 400])

outstem = [prefix '-normvarerr'];
fname = sprintf('%s.tif',outstem);
%print('-dtiff',fname);
%mmil_printeps(gcf,[outstem '.eps']);
