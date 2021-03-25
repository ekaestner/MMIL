function rc_plot_condinfo_stimlocs(fname_conds)
%function rc_plot_condinfo_stimlocs(fname_conds)
%
% Purpose: plot 2D locations of stimuli from cond_info.csv file
%
% Required Input:
%   fname_conds: full or relative path file name of
%     cond_info csv file (e.g. cond_info.csv)
%
% Early Mod: 10/01/10 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%


if ~mmil_check_nargs(nargin,1), return; end;

cond_info = rc_read_cond_info(fname_conds);
r_max = 12.5;

text_scalef = 1;
text_xoff = 0;
%text_scalef = 1.2;
%text_xoff = -0.5;

clf;
hold on;
for i=1:length(cond_info)
  if cond_info(i).contrast==0, continue; end;
  x = cond_info(i).ecc*cos(cond_info(i).theta*pi/180);
  y = cond_info(i).ecc*sin(cond_info(i).theta*pi/180);
  x2 = text_scalef*x + text_xoff;
  y2 = text_scalef*y;

  plot(x,y,'b.');
  text(x2,y2,num2str(cond_info(i).event_code));
end;

% plot axes
line([-r_max r_max],[0 0],[0 0],...
       'Color','k','LineWidth',1);
xlabel('x');

line([0 0],[-r_max r_max],[0 0],...
       'Color','k','LineWidth',1);
ylabel('y');

