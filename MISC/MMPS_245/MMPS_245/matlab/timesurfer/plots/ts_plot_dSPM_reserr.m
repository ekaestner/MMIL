function ts_plot_dSPM_reserr(prefix,cond,types_flag)
%function ts_plot_dSPM_reserr(prefix,cond,[types_flag])
%
% last mod: 08/05/09 by Don Hagler
%

if (~mmil_check_nargs(nargin,2)) return; end;
if ~exist('types_flag','var') | isempty(types_flag), types_flag = 0; end;

matfile=sprintf('matfiles/%s_results_cond%02d.mat',prefix,cond);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

matfile=sprintf('matfiles/%s_avg_data.mat',prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

fprintf('%s: plotting residual error for %s cond %d\n',mfilename,prefix,cond);

clf;
time = avg_data.averages(1).time*1000;
plot(time,norm_var_E,'k','LineWidth',2);
if types_flag
  hold on;
  plot(time,norm_var_Eeeg,'r');
  plot(time,norm_var_Egrad,'g');
  plot(time,norm_var_Emag,'b');
  legend('mean','EEG','grad','mag','Location','SouthEast');
end;

title('normalized variance of error');
axis('tight');
set(gca,'YLim',[0 1]); % normalized residual variance

print('-dtiff',sprintf('%s-reserr.tif',prefix));

