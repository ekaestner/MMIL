function ts_plot_dSPM_sources(prefix,cond)
%function ts_plot_dSPM_sources(prefix,cond)
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

fprintf('%s: plotting source waveforms for %s cond %d\n',...
  mfilename,prefix,cond);

time = avg_data.averages(1).time*1000;
plot(time,S);
axis tight;
title('Source Waveforms');

print('-dtiff',sprintf('%s-sources.tif',prefix));

