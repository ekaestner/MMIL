function mmil_save_fig(fig_num,fname_tif,fig_size,tif_dpi,eps_flag,visible_flag)
%function mmil_save_fig(fig_num,fname_tif,fig_size,tif_dpi,eps_flag,visible_flag)
%
% created:  04/20/16 by Don Hagler
% last mod: 04/20/16 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('fig_size','var') || isempty(fig_size), fig_size = [4 4]; end;
if ~exist('tif_dpi','var') || isempty(tif_dpi), tif_dpi = 300; end;
if ~exist('eps_flag','var') || isempty(eps_flag), eps_flag = 0; end;
if ~exist('visible_flag','var') || isempty(visible_flag), visible_flag = 0; end;

if ~visible_flag, set(fig_num,'visible','off'); end;
set(fig_num, 'PaperUnits', 'inches');
set(fig_num, 'PaperSize', [fig_size(1) fig_size(2)]);
set(fig_num, 'PaperPositionMode', 'manual');
set(fig_num, 'PaperPosition', [0 0 fig_size(1) fig_size(2)]);
print(fig_num,'-dtiff',fname_tif,sprintf('-r %d',tif_dpi));
if eps_flag
  fname_eps = regexprep(fname_tif,'tif','eps');
  mmil_printeps(fig_num,fname_eps);
end;
if ~visible_flag, close(fig_num); end;

return;

