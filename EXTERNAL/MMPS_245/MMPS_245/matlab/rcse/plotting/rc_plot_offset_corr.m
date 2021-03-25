function rc_plot_offset_corr(prefix,varargin)
%function rc_plot_offset_corr(prefix,[options])
%
% Required Input:
%   prefix: RCSE prefix
%
% Optional Input:
%  'self_corr_flag': calculate correlation between offset grid locations
%     within a quarterfield
%    {default = 1}
%  'ref_corr_flag': calculate correlation with reference waveform
%    and use as weighting or mask for selection of dipoles
%    reference waveform assumes zero baseline, 1 between fit_time0 and fit_time1
%    {default = 1}
%  'ref_corr_pol': [-1,0,1] polarity of reference waveform
%    if -1, assumes negative polarity
%    if 0, allows either negative or positive polarity
%    if 1, asssumes positive polarity
%    {default = 0}
%  'ref_corr_thresh': correlation with reference threshold
%    0 keeps values, but must be positive
%    otherise, apply threshold and binarize
%    {default = 0}
%  'expvar_thresh': explained variance threshold
%    0 keeps values, but must be positive
%    otherise, apply threshold and binarize
%    {default = 0}
%  'expvar_fact': exponential weighting applied to explained variance in
%    determining which offsets to choose
%    {default = 0}
%  'plot_outdir': output directory for plot tif files
%    {default = 'plots-offsets'}
%
% Created:  06/09/09 by Don Hagler
% Last Mod: 02/20/11 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin, { ...
  'prefix',prefix,[],...
  'expvar_thresh',0,[0,1],...
  'ref_corr_thresh',0,[0,1],...
  'expvar_fact',5,[0,1000],...
  'self_corr_flag',true,[false true],...
  'ref_corr_flag',false,[false true],...
  'ref_corr_pol',0,[-1 0 1],...
  'plot_outdir','plots',[],...
...
  'expvar_clim',[],[],...
  'color_order',{'b' 'g' 'r' 'c' 'm' 'y' 'k'},[],...
  'hemistrs',{'right','left'},{'right','left'},...
  'uplowstrs',{'upper','lower'},{'upper','lower'},...
  'areanames',{'V1','V2','V3'},{'V1','V2','V3'},...
  'wf_ylim',[-15,15],[],...
  'norm_wf_flag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create reference waveform

% load parameters from retinv
matfile = sprintf('matfiles/%s_parms.mat',parms.prefix);
if ~exist(matfile,'file')
  error('file %s not found',matfile);
end;
tmp_parms = parms;
load(matfile);
ri_parms = parms;
parms = tmp_parms;

baseline_start = ri_parms.baseline_start_samp;
baseline_end = ri_parms.baseline_end_samp;
fit_start = ri_parms.fit_t0;
fit_end = ri_parms.fit_t1;

% load time vector
matfile = sprintf('matfiles/%s_avg_data.mat',parms.prefix);
load(matfile);
time = avg_data.averages(1).time;

if isempty(ri_parms.use_areas)
  parms.areas = [1:length(parms.areanames)];
else
  parms.areas = ri_parms.use_areas;
end;

% load offset waveforms
wforms = [];
err_wforms = [];
matfile = sprintf('matfiles/%s_ret_forward.mat',parms.prefix);
if ~exist(matfile,'file')
  error('file %s not found\n',mfilename,matfile);
end;
clear retforward;
load(matfile);
wforms = retforward.retmap.offset_wforms;
err_wforms = retforward.retmap.offset_err_wforms;

matsize = size(wforms);
gridsize = matsize(1:2);
ngrid = prod(gridsize);
nareas = matsize(3);
ntpoints = matsize(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlations with reference

if parms.ref_corr_flag
%  fprintf('%s: calculating correlations to reference...\n',mfilename);
  % create reference waveform
  baseline_range = [baseline_start:baseline_end];
  fit_range = [fit_start:fit_end];
  nbaseline = length(baseline_range);
  nfit = length(fit_range);
  wform_ref = zeros(nbaseline+nfit,1);
  if parms.ref_corr_pol==0
    wform_ref(nbaseline+1:end) = 1;
  else
    wform_ref(nbaseline+1:end) = parms.ref_corr_pol;
  end;
  R_ref = zeros(nareas,ngrid);
  for a=parms.areas
    tmp_wforms = squeeze(wforms(:,:,a,:));
    tmp_wforms = reshape(tmp_wforms,[ngrid,ntpoints])';
    tmp_wforms = tmp_wforms([baseline_range,fit_range],:);
    tmp_R = corr(tmp_wforms,wform_ref);
    R_ref(a,:) = tmp_R;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create probability maps based on explained variance and ref corr
probmap = zeros([nareas,gridsize]);
err = mean(err_wforms(:,:,fit_start:fit_end),3);
expvar = 1-err;
if parms.expvar_thresh
  expvar(expvar<parms.expvar_thresh)=0;
  expvar(expvar>=parms.expvar_thresh)=1;
else
  expvar(expvar<0)=0;
end;
expvar = expvar.^parms.expvar_fact;
for a=parms.areas
  if parms.ref_corr_flag
    tmp_R_ref = reshape(R_ref(a,:),gridsize);
    if parms.ref_corr_pol==0
      tmp_R_ref = abs(tmp_R_ref);
    end;
    if parms.ref_corr_thresh
      tmp_R_ref(tmp_R_ref<parms.ref_corr_thresh) = 0;
      tmp_R_ref(tmp_R_ref>=parms.ref_corr_thresh) = 1;
    else
      tmp_R_ref(tmp_R_ref<0) = 0;
    end;
    probmap(a,:,:) = expvar.*tmp_R_ref;
  else
    probmap(a,:,:) = expvar;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find offsets

offsets_1D = zeros(nareas,1);
offsets_2D = zeros(nareas,2);
for a=parms.areas
  maskvec = reshape(probmap(a,:,:),[ngrid,1]);
  if all(maskvec==1)
    ind = round(length(maskvec)/2);
  else
    [max_r,ind] = max(maskvec);
  end;
  [ind_r,ind_th] = ind2sub(gridsize,ind);
  offsets_1D(a) = ind;
  offsets_2D(a,1) = ind_r;
  offsets_2D(a,2) = ind_th;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate correlation to self for best offsets

R_self = zeros([nareas,gridsize]);
if parms.self_corr_flag
%  fprintf('%s: calculating correlations to self...\n',mfilename);
  for a=parms.areas
    ind_r = offsets_2D(a,1);
    ind_th = offsets_2D(a,2);
    best_wform = squeeze(wforms(ind_r,ind_th,a,:));
    best_wform = reshape(best_wform,[1,ntpoints])';

    tmp_wforms = squeeze(wforms(:,:,a,:));
    tmp_wforms = reshape(tmp_wforms,[ngrid,ntpoints])';
    
    tmp_R_self = corr(tmp_wforms,best_wform);
    tmp_R_self = reshape(tmp_R_self,gridsize);
    R_self(a,:,:) = tmp_R_self;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots

%fprintf('%s: plotting results...\n',mfilename);
f = 1;
if ~strncmp(parms.plot_outdir,'/',1) % relative
  parms.plot_outdir = [pwd '/' parms.plot_outdir];
end;
if ~exist(parms.plot_outdir,'dir')
  [s,m]=mkdir(parms.plot_outdir);
  if ~s
    error('failed to create plot_outdir %s:\n%s',parms.plot_outdir,m);
  end;
end;

% plot explained variance on offset grid
figure(f); clf;
plotname = 'explained variance';
if ~isempty(parms.expvar_clim)
  imagesc(expvar,parms.expvar_clim);
else
  imagesc(expvar);
end;
colorbar;
title(plotname);
fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
set(gcf,'Visible','off'); print('-dtiff',fname);

for a=parms.areas
  % plot correlation with reference
  if parms.ref_corr_flag
    figure(f); clf;
    plotname = 'correlation with reference';
    plotname = [parms.areanames{a} ' ' plotname];
    im = reshape(R_ref(a,:),gridsize);
    imagesc(im,[-1,1]); colorbar;
    title(plotname);
    fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
    set(gcf,'Visible','off'); print('-dtiff',fname);
  end;

  % plot dipole location probability
  figure(f); clf;
  plotname = 'dipole location probability';
  plotname = [parms.areanames{a} ' ' plotname];
  im = squeeze(probmap(a,:,:));
  im = im/(eps+max(im(:)));
  imagesc(im); colorbar;
  title(plotname);
  fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
  set(gcf,'Visible','off'); print('-dtiff',fname);

  % plot self-correlations with dipole offsets
  ylim = [-1,1];
  figure(f); clf; p = 1;
  plotname = 'dipole with self-correlation';
  plotname = [parms.areanames{a} ' ' plotname];
  subplot(1,2,p); p=p+1;
  im = zeros(gridsize);
  ind = offsets_1D(a);
  im(ind) = 1;
  imagesc(im,ylim); colorbar;
  title('grid loc');
  subplot(1,2,p); p=p+1;
  im = squeeze(R_self(a,:,:));
  imagesc(im,ylim); colorbar;
  title('self corr');
  suptitle(plotname);
  fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
  set(gcf,'Visible','off'); print('-dtiff',fname);
end;

% plot source waveforms
figure(f); clf;
plotname = 'source waveforms';
ylim = parms.wf_ylim;
hold on;
for a=parms.areas
  ind_r = offsets_2D(a,1);
  ind_th = offsets_2D(a,2);
  wform = squeeze(wforms(ind_r,ind_th,a,:));
  if parms.norm_wf_flag
    wform = wform/abs(mean(wform(fit_start:fit_end)));
  end;
  plot(1000*time,wform,parms.color_order{a});
end;
legend(parms.areanames(parms.areas),'Location','EastOutside');
xlabel('time (ms)')
axis tight;
set(gca,'YLim',ylim);
title(plotname);
fname = sprintf('%s/%s.tif',parms.plot_outdir,regexprep(plotname,' ','_'));
set(gcf,'Visible','off'); print('-dtiff',fname);

close(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

