function fname = rc_plot_fitvar(prefix,fontname,linewidth,label_flag)
%function fname = rc_plot_fitvar(prefix,[fontname],[linewidth],[label_flag])
%
% Purpose: plot variance of data, fit, and residual error of RCSE
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%   fontname
%     {default = 'Arial'}
%   linewidth
%     {default = 1.5 }
%   label_flag
%     {default = 1}
%
% Early Mod: 05/03/09 by Don Hagler
% Last Mod:  08/01/13 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('fontname','var') || isempty(fontname), fontname = 'Arial'; end;
if ~exist('linewidth','var') || isempty(linewidth), linewidth = 1.5; end;
if ~exist('label_flag','var') || isempty(label_flag), label_flag = 1; end;

fprintf('%s: plotting fit variance for %s\n',mfilename,prefix);

matname=sprintf('matfiles/%s_parms.mat',prefix);
load(matname);

matname=sprintf('matfiles/%s_results.mat',prefix);
load(matname);

rc_plot_fitvar_results(parms,results,label_flag,linewidth,[],fontname);

set(gca,'XGrid','on','XTick',[-100 0 100 200 300 400]);

outstem = [prefix '-fitvar'];
fname = sprintf('%s.tif',outstem);
%print('-dtiff',fname);
%mmil_printeps(gcf,[outstem '.eps']);

