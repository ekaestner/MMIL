function ts_plot_dSPM_crosstalk(prefix)
%function ts_plot_dSPM_crosstalk([prefix])

if ~exist('prefix','var'), prefix = 'dSPM'; end;

ds = 10;

matfile=sprintf('matfiles/%s_inverse.mat',prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
load(matfile);

fprintf('%s: calculating crosstalk for %s...\n',mfilename,prefix);

keyboard

CT = M*G_xyz;

tmp = CT(1:ds:end,1:ds:end);
imagesc(tmp);
colorbar;

print('-dtiff',sprintf('%s-crosstalk.tif',prefix));

