function ts_plot_dSPM_reserr(prefix,cond)
%function ts_plot_dSPM_reserr(prefix,cond)
%
% last mod: 08/05/09 by Don Hagler
%

if (~mmil_check_nargs(nargin,2)) return; end;

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

fprintf('%s: plotting data, fit, and residual error for %s cond %d\n',...
  mfilename,prefix,cond);

time = avg_data.averages(1).time*1000;

titles = {'Data' 'Fit' 'Error' 'norm var err'};

for i=1:4
  h(i) = subplot(2,2,i);
end

plot(h(1),time,Y);
plot(h(2),time,Yfit);
plot(h(3),time,E);
plot(h(4),time,norm_var_E);

for i=1:4
  title(h(i),titles{i});
  axis(h(i),'tight');
end
set(h(4),'YLim',[0 1]); % normalized residual variance

print('-dtiff',sprintf('%s-fit.tif',prefix));

