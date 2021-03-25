function rc_plot_fitvar_results(parms,results,label_flag,linewidth,fontsize,fontname,baselineflag,avgflag)
%function rc_plot_fitvar_results(parms,results,[label_flag],[linewidth],[fontsize],[fontname],[baselineflag],[avgflag])
%
% Early Mod: 09/02/09 by Don Hagler
% Last Mod:  02/20/11 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

if ~exist('label_flag','var') | isempty(label_flag), label_flag = 1; end;
if ~exist('fontname','var') | isempty(fontname), fontname = 'Arial'; end;
if ~exist('fontsize','var') | isempty(fontsize), fontsize = 12; end;
if ~exist('linewidth','var') | isempty(linewidth), linewidth = 1.5; end;
if ~exist('baselineflag','var') | isempty(baselineflag), baselineflag = 0; end;
if ~exist('avgflag','var') | isempty(avgflag), avgflag = 1; end;
linewidth_mult = 1.5;

time = results.time*1000;

Ygrad = results.Ygrad;
Yfitgrad = results.Yfitgrad;
Egrad = results.Egrad;

Ymag = results.Ymag;
Yfitmag = results.Yfitmag;
Emag = results.Emag;

Yeeg = results.Yeeg;
Yfiteeg = results.Yfiteeg;
Eeeg = results.Eeeg;

if avgflag
  var_Yavg = 0;
  var_Yfitavg = 0;
  var_Eavg = 0;
  ntype = 0;
end;

clf;
hold on;

if ~isempty(Ygrad)
  var_Ygrad = var(Ygrad,0,2);
  var_Yfitgrad = var(Yfitgrad,0,2);
  var_Egrad = var(Egrad,0,2);
  if baselineflag
    baseline_var_Egrad = ...
      mean(var_Egrad(parms.baseline_start_samp:parms.baseline_end_samp));
    var_Ygrad = var_Ygrad - baseline_var_Egrad;
    var_Egrad = var_Egrad - baseline_var_Egrad;
  end;
  max_var_Ygrad = max(var_Ygrad);
  var_Ygrad = var_Ygrad/max_var_Ygrad;
  var_Yfitgrad = var_Yfitgrad/max_var_Ygrad;
  var_Egrad = var_Egrad/max_var_Ygrad;
  if ~avgflag
    plot(time,var_Ygrad,'k','LineWidth',linewidth);
    plot(time,var_Yfitgrad,'b','LineWidth',linewidth);
    plot(time,var_Egrad,'r','LineWidth',linewidth);
    linewidth = linewidth*linewidth_mult;
  else
    var_Yavg = var_Yavg + var_Ygrad;
    var_Yfitavg = var_Yfitavg + var_Yfitgrad;
    var_Eavg = var_Eavg + var_Egrad;
    ntype = ntype + 1;
  end;
end;

if ~isempty(Ymag)
  var_Ymag = var(Ymag,0,2);
  var_Yfitmag = var(Yfitmag,0,2);
  var_Emag = var(Emag,0,2);
  if baselineflag
    baseline_var_Emag = ...
      mean(var_Emag(parms.baseline_start_samp:parms.baseline_end_samp));
    var_Ymag = var_Ymag - baseline_var_Emag;
    var_Emag = var_Emag - baseline_var_Emag;
  end;
  max_var_Ymag = max(var_Ymag);
  var_Ymag = var_Ymag/max_var_Ymag;
  var_Yfitmag = var_Yfitmag/max_var_Ymag;
  var_Emag = var_Emag/max_var_Ymag;
  if ~avgflag
    plot(time,var_Ymag,'k','LineWidth',linewidth);
    plot(time,var_Yfitmag,'b','LineWidth',linewidth);
    plot(time,var_Emag,'r','LineWidth',linewidth);
    linewidth = linewidth*linewidth_mult;
  else
    var_Yavg = var_Yavg + var_Ymag;
    var_Yfitavg = var_Yfitavg + var_Yfitmag;
    var_Eavg = var_Eavg + var_Emag;
    ntype = ntype + 1;
  end;
end;

if ~isempty(Yeeg)
  var_Yeeg = var(Yeeg,0,2);
  var_Yfiteeg = var(Yfiteeg,0,2);
  var_Eeeg = var(Eeeg,0,2);
  if baselineflag
    baseline_var_Eeeg = ...
      mean(var_Eeeg(parms.baseline_start_samp:parms.baseline_end_samp));
    var_Yeeg = var_Yeeg - baseline_var_Eeeg;
    var_Eeeg = var_Eeeg - baseline_var_Eeeg;
  end;
  max_var_Yeeg = max(var_Yeeg);
  var_Yeeg = var_Yeeg/max_var_Yeeg;
  var_Yfiteeg = var_Yfiteeg/max_var_Yeeg;
  var_Eeeg = var_Eeeg/max_var_Yeeg;
  if ~avgflag
    plot(time,var_Yeeg,'k','LineWidth',linewidth);
    plot(time,var_Yfiteeg,'b','LineWidth',linewidth);
    plot(time,var_Eeeg,'r','LineWidth',linewidth);
  else
    var_Yavg = var_Yavg + var_Yeeg;
    var_Yfitavg = var_Yfitavg + var_Yfiteeg;
    var_Eavg = var_Eavg + var_Eeeg;
    ntype = ntype + 1;
  end;
end;

if avgflag
  var_Yavg = var_Yavg/ntype;
  var_Yfitavg = var_Yfitavg/ntype;
  var_Eavg = var_Eavg/ntype;
  plot(time,var_Yavg,'k','LineWidth',linewidth);
  plot(time,var_Yfitavg,'b','LineWidth',linewidth);
  plot(time,var_Eavg,'r','LineWidth',linewidth);
end;

axis('tight');
set(gca,'YLim',[0 1]);

if label_flag
  set(gca,'FontSize',fontsize,'FontName',fontname)
  legend({'data','fit','residual'},'FontSize',fontsize,'FontName',fontname)
  title('variance of data, fit, and residual','FontSize',fontsize,'FontName',fontname)
  xlabel('Time (msec)','FontSize',fontsize,'FontName',fontname)
  ylabel('normalized variance','FontSize',fontsize,'FontName',fontname)
end;

