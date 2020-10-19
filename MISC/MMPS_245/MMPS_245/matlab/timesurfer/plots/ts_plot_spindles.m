function ts_plot_spindles(results,varargin)
%function ts_plot_spindles(results,[options])
%
% Purpose: plot spindles, spindle power, and fft power in subplots
%
% Required input:
%   results: output struct from ts_find_spindles
%
% Optional parameters:
%   'outdir': output directory for plots
%     {default = pwd}
%   'outstem_spindle': output file stem for spindle plots
%     {default = 'spindle'}
%   'outstem_reject': output file stem for reject plots
%     {default = 'reject'}
%   'image_flag': [0|1] plot chan x time matrices as colorscale images
%     {default = 0}
%   'reject_flag': plot rejected spindles
%     {default = 0}
%   'filt_broad_flag': [0|1] plot broad bandpass-filtered data time courses
%     {default = 1}
%   'filt_narrow_flag': [0|1] plot narrow bandpass-filtered data time courses
%     {default = 1}
%   'peak_flag': [0|1] plot spindle peak power time courses
%     {default = 1}
%   'edge_flag': [0|1] plot spindle edge power time courses
%     {default = 1}
%   'fft_flag': [0|1] plot Fourier power spectra
%     {default = 1}
%   'spindle_freq_low': lower bound of spindle frequency range (in Hz)
%     {default = 9}
%   'spindle_freq_high': upper bound of spindle frequency range (in Hz)
%     {default = 17}
%   'nonspindle_freq_low': lower bound of non-spindle frequency range (in Hz)
%     {default = 5}
%   'nonspindle_freq_high': upper bound of non-spindle frequency range (in Hz)
%     {default = 8}
%   'min_peak_deflection': fraction of channel-wise median abs deviation
%      if the difference between a local maximum and surrounding minima
%      is smaller than min_peak_deflection * med_abs
%      it is not considered a likely peak in the spindle oscillation
%      for annotation of npeaks
%     {default = 1}
%   'min_peak_ratio': exclude peaks from npeaks count if deflection
%      is smaller than min_peak_ratio * largest peak deflection for the epoch
%     {default = 0.25}
%   'time_max': maximum time for time course plots
%     {default = 1000}
%   'fft_max': maximum frequency for FFT plots
%     {default = 50}
%   'vert_flag': [0|1] arrange subplots vertically
%     otherwise, use a grid
%     {default = 0}
%   'chan_colors': cell array of color codes for each channel
%     {default = {'b','g','m','y','k','r'}}
%   'annot_flag': [0|1] include annotations on plots
%     {default = 1}
%   'fixed_ylim_flag': [0|1] used fixed limits for vertical axes
%     otherwise, fit to max vals for each spindle
%     {default = 0}
%   'visible_flag': [0|1] make plots visible
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  06/10/14 by Don Hagler
% Last Mod: 04/15/15 by Don Hagler
%

%% todo: option to arrange channels vertically with some spacing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'outdir',pwd,[],...
  'outstem_spindle','spindle',[],...
  'outstem_reject','reject',[],...
  'image_flag',false,[false true],...
  'reject_flag',false,[false true],...
  'filt_broad_flag',true,[false true],...
  'filt_narrow_flag',true,[false true],...
  'peak_flag',true,[false true],...
  'edge_flag',true,[false true],...
  'fft_flag',true,[false true],...
  'spindle_freq_low',9,[0.5,100],...
  'spindle_freq_high',17,[0.5,100],...
  'nonspindle_freq_low',5,[0,100],...
  'nonspindle_freq_high',8,[0,100],...
  'min_peak_deflection',1,[0,Inf],...
  'min_peak_ratio',0.25,[0,Inf],...
  'time_max',1000,[],...
  'fft_max',50,[],...
  'vert_flag',false,[false true],...
  'chan_colors',{'b','g','m','c','y','k','r'},[],...
  'annot_flag',true,[false true],...
  'fixed_ylim_flag',false,[false true],...
  'visible_flag',false,[false true],...
  'fontsize',8,[],...
  'forceflag',false,[false true],...
  ...
  'ylim_data',[],[],...
  'ylim_filt_broad',[],[],...
  'ylim_filt_narrow',[],[],...
  'ylim_peak',[],[],...
  'ylim_edge',[],[],...
  'ylim_fft',[],[],...
  'fft_freq_norm_flag',false,[false true],...
  'fft_freq_pow',1.2,[1,10],...
});

if parms.vert_flag
  parms.nplots = 2;
else
  parms.nplots = 1;
end;
if parms.filt_broad_flag
  parms.nplots = parms.nplots + 1;
end;
if parms.filt_narrow_flag
  parms.nplots = parms.nplots + 1;
end;
if parms.peak_flag
  parms.nplots = parms.nplots + 1;
end;
if parms.edge_flag
  parms.nplots = parms.nplots + 1;
end;
if parms.fft_flag
  parms.nplots = parms.nplots + 1;
end;

if parms.reject_flag
  parms.reject_flags = [0,1];
else
  parms.reject_flags = 0;
end;

if ~isempty(results.data_filt_broad_mad)
  parms.data_mad = results.data_filt_broad_mad;
elseif ~isempty(results.data_filt_narrow_mad)
  parms.data_mad = results.data_filt_narrow_mad;
else
  parms.data_mad = results.data_mad;
end;

mmil_mkdir(parms.outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parms.fixed_ylim_flag
  % find ylim for data
  if isempty(parms.ylim_data)
    parms.ylim_data = 1.05*prctile(results.spindle_data(:),[0.1,99.9]);
  end;
  % find ylim for broad-band filtered data
  if parms.filt_broad_flag && isempty(parms.ylim_filt_broad)
    parms.ylim_filt_broad = 1.05*prctile(results.spindle_data_filt_broad(:),[0.1,99.9]);
  end;
  % find ylim for narrow-band filtered data
  if parms.filt_narrow_flag && isempty(parms.ylim_filt_narrow)
    parms.ylim_filt_narrow = 1.05*prctile(results.spindle_data_filt_narrow(:),[0.1,99.9]);
  end;
  % find ylim for peak wave
  if parms.peak_flag && isempty(parms.ylim_peak)
    parms.ylim_peak = 1.05*prctile(results.spindle_power_peak(:),[0.1,99.9]);
  end;
  % find ylim for edge wave
  if parms.edge_flag && isempty(parms.ylim_edge)
    parms.ylim_edge = 1.05*prctile(results.spindle_power_edge(:),[0.1,99.9]);
  end;
  % set ylim for fft
  if parms.fft_flag && isempty(parms.ylim_fft)
    parms.ylim_fft = [0,0.25*parms.ylim_data(2)];
    if parms.fft_freq_norm_flag
      parms.ylim_fft = ...
        parms.ylim_fft .* (parms.spindle_freq_low.^parms.fft_freq_pow);
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r=1:length(parms.reject_flags)
  reject_flag = parms.reject_flags(r);
  if reject_flag
    if ~isfield(results,'reject_data'), continue; end;
    outstem = parms.outstem_reject;
    nspindles = size(results.ind_reject,1);
    if nspindles==0, continue; end;
    spindle_data = results.reject_data;
    spindle_data_filt_broad = results.reject_data_filt_broad;
    spindle_data_filt_narrow = results.reject_data_filt_narrow;
    spindle_power_peak = results.reject_power_peak;
    spindle_power_edge = results.reject_power_edge;
    ind_spindle_rel = results.ind_reject_rel;
  else
    outstem = parms.outstem_spindle;
    nspindles = size(results.ind_spindle,1);
    if nspindles==0, continue; end;
    spindle_data = results.spindle_data;
    spindle_data_filt_broad = results.spindle_data_filt_broad;
    spindle_data_filt_narrow = results.spindle_data_filt_narrow;
    spindle_power_peak = results.spindle_power_peak;
    spindle_power_edge = results.spindle_power_edge;
    ind_spindle_rel = results.ind_spindle_rel;
  end;

  for i=1:nspindles
    fname_tif = sprintf('%s/%s_%04d.tif',parms.outdir,outstem,i);
    if ~exist(fname_tif,'file') || parms.forceflag

      figure(1); clf;
      if ~parms.visible_flag
        set(gcf,'visible','off');
      end;
      k = 1;

      % plot spindle time course
      k = plot_spindle_data(spindle_data,...
                            parms,results,ind_spindle_rel,i,k,2,...
                            'data','time (sec)',1,parms.ylim_data);
      if isempty(spindle_data_filt_narrow) &&...
         isempty(spindle_data_filt_broad)
        % display number of sufficiently large peaks
        if parms.annot_flag
          annot_npeaks(spindle_data,parms,results,ind_spindle_rel,i);
        end;
      end;

      % plot broad-bandpass-filtered data
      if parms.filt_broad_flag && ~isempty(spindle_data_filt_broad)
        k = plot_spindle_data(spindle_data_filt_broad,...
                              parms,results,ind_spindle_rel,i,k,1,...
                              'broad-band filtered data','time (sec)',1,...
                              parms.ylim_filt_broad);
        % display number of sufficiently large peaks
        if parms.annot_flag
          annot_npeaks(spindle_data_filt_broad,parms,results,ind_spindle_rel,i);
        end;
      end;

      % plot narrow-bandpass-filtered data
      if parms.filt_narrow_flag && ~isempty(spindle_data_filt_narrow)
        k = plot_spindle_data(spindle_data_filt_narrow,...
                              parms,results,ind_spindle_rel,i,k,1,...
                              'narrow-band filtered data','time (sec)',1,...
                              parms.ylim_filt_narrow);
        % display number of sufficiently large peaks
        if isempty(spindle_data_filt_broad)
          if parms.annot_flag
            annot_npeaks(spindle_data_filt_narrow,parms,results,ind_spindle_rel,i);
          end;
        end;
      end;

      % plot spindle peak power (wavelet convolved)
      if parms.peak_flag
        k = plot_spindle_data(spindle_power_peak,...
                            parms,results,ind_spindle_rel,i,k,1,...
                              'peak wavelet','time (sec)',0,parms.ylim_peak);
        % display peak amplitude
        if parms.annot_flag
          annot_amp(spindle_power_peak,parms,results,ind_spindle_rel,i);
        end;
      end;

      % plot spindle edge power (wavelet convolved)
      if parms.edge_flag
        k = plot_spindle_data(spindle_power_edge,...
                              parms,results,ind_spindle_rel,i,k,1,...
                              'edge wavelet','time (sec)',0,parms.ylim_edge);
        % display duration
        if parms.annot_flag
          annot_dur(spindle_power_edge,parms,results,ind_spindle_rel,i);
        end;
      end;

      % plot Fourier power spectrum
      if parms.fft_flag
        k = plot_fft(spindle_data,parms,results,ind_spindle_rel,i,k,...
                     parms.ylim_fft);
      end;

      % save plot as tif
      print(gcf,'-dtiff',fname_tif);
      if ~parms.visible_flag
        close(gcf);
      end;
    end;
  end;
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = plot_spindle_data(data,parms,results,ind_spindle_rel,...
                               i,k,np,title_str,xlabel_str,signed_flag,vlim)
  if ~exist('vlim','var'), vlim = []; end;

  if parms.vert_flag
    subplot(parms.nplots,1,k:k+np-1);
    k = k + np;
  else
    nc = ceil(sqrt(parms.nplots));
    nr = ceil(parms.nplots/nc);
    subplot(nr,nc,k);
    k = k + 1;
  end;

  hold on;
  ntpoints = size(data,3);
  time = [0:ntpoints-1]*1000/results.sfreq;
  max_time = time(end)/2;
  time = time - max_time;
  tmp_data = reshape(data(i,:,:),...
    [results.nchans,ntpoints]);
  max_val = 1.01*max(abs(tmp_data(:)));
  if isempty(vlim)
    if ~signed_flag
      vlim = [0,max_val];
    else
      vlim = [-max_val,max_val];
    end;
  end;
  if parms.image_flag
    imagesc(time,[1:results.nchans],tmp_data,vlim);
    axis ij;
    colorbar
    ylabel('channel','FontSize',parms.fontsize);
    ylim([1,results.nchans]);
  else
    j = 1;
    for c=1:results.nchans
      chan_col = parms.chan_colors{j};
      j = j + 1;
      if j>length(parms.chan_colors), j = 1; end;
      plot(time,tmp_data(c,:),chan_col);
    end;
    t0 = ind_spindle_rel(i,1)*1000/results.sfreq - max_time;
    t1 = ind_spindle_rel(i,2)*1000/results.sfreq - max_time;
    line([t0,t0],vlim,'color','r');
    line([t1,t1],vlim,'color','r');
    ylim(vlim);
  end;
  xlim([-parms.time_max,parms.time_max]);
  if ~isempty(title_str)
    title(title_str,'FontSize',parms.fontsize);
  end;
  if ~isempty(xlabel_str)
    xlabel(xlabel_str,'FontSize',parms.fontsize);
  end;
  set(gca,'FontSize',parms.fontsize);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = plot_fft(data,parms,results,ind_spindle_rel,i,k,vlim)
  if ~exist('vlim','var'), vlim = []; end;
  title_str = 'FFT power spectrum';
  xlabel_str = 'frequency (Hz)';

  if parms.vert_flag
    subplot(parms.nplots,1,k);
  else
    nc = ceil(sqrt(parms.nplots));
    nr = ceil(parms.nplots/nc);
    subplot(nr,nc,k);
  end;
  ind_first = ind_spindle_rel(i,1);
  ind_last  = ind_spindle_rel(i,2);
  k = k + 1;
  if parms.image_flag
    fft_data = [];
    for c=1:results.nchans
      % extract data for epoch
      epoch_data = squeeze(data(i,c,ind_first:ind_last))';
      [fft_freqs,fft_power] = calc_fft(epoch_data,parms,results);
      if isempty(fft_data)
        fft_data = nan(results.nchans,length(fft_freqs));
      end;
      if parms.fft_freq_norm_flag
        fft_power = fft_power .* (fft_freqs.^parms.fft_freq_pow);
      end;
      fft_data(c,:) = fft_power;
    end;
    max_val = 1.01*max(mmil_colvec(fft_data(:,1:parms.fft_max)));
    if isempty(vlim)
      vlim = [0,max_val];
    end;
    imagesc(fft_freqs,[1:results.nchans],fft_data,vlim);
    axis ij;
    colorbar
    ylabel('channel','FontSize',parms.fontsize);
    ylim([1,results.nchans]);
  else
    hold on;
    j = 1;
    for c=1:results.nchans
      chan_col = parms.chan_colors{j};
      j = j + 1;
      if j>length(parms.chan_colors), j = 1; end;
      % extract data for epoch
      epoch_data = squeeze(data(i,c,ind_first:ind_last))';
      [fft_freqs,fft_power,fft_ratio,spindle_power,nonspindle_power] =...
           calc_fft(epoch_data,parms,results);
      if parms.fft_freq_norm_flag
        fft_power = fft_power .* (fft_freqs.^parms.fft_freq_pow);
      end;
      plot(fft_freqs,fft_power,chan_col);
      if parms.annot_flag
        x = 0.5*parms.fft_max;
        y = 0.9*spindle_power;
        text(x,y,sprintf('F ratio = %0.1f',fft_ratio),'color',chan_col);
      end;
    end;
    if isempty(vlim)
      axis tight;
    else
      ylim(vlim);
    end;
  end;
  xlim([0,parms.fft_max]);
  title(title_str,'FontSize',parms.fontsize);
  xlabel(xlabel_str,'FontSize',parms.fontsize);
  set(gca,'FontSize',parms.fontsize);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fft_freqs,fft_power,fft_ratio,spindle_power,nonspindle_power] =...
         calc_fft(epoch_data,parms,results)
  % remove linear trend
  epoch_data = detrend(epoch_data);
  % calculate Fourier power
  ntpoints = length(epoch_data);
  nfft = 2^nextpow2(ntpoints);
  fft_power = fft(epoch_data,nfft)/ntpoints;
  fft_power = 2*abs(fft_power(1:nfft/2+1));
  % calculate frequencies
  fft_freqs = results.sfreq/2*linspace(0,1,nfft/2+1);
  % calculate indices for frequencies of interest
  [tmp,ind_freq_low] = min(abs(fft_freqs - parms.spindle_freq_low));
  [tmp,ind_freq_high] = min(abs(fft_freqs - parms.spindle_freq_high));
  ind_spindle_freqs = ind_freq_low:ind_freq_high;
  [tmp,ind_freq_low] = min(abs(fft_freqs - parms.nonspindle_freq_low));
  [tmp,ind_freq_high] = min(abs(fft_freqs - parms.nonspindle_freq_high));
  ind_nonspindle_freqs = ind_freq_low:ind_freq_high;
  % average power across frequency ranges
  spindle_power = max(fft_power(ind_spindle_freqs));
  nonspindle_power = max(fft_power(ind_nonspindle_freqs));
  % calculate ratio of spindle power to non-spindle power
  fft_ratio = spindle_power / nonspindle_power;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function annot_npeaks(data,parms,results,ind_spindle_rel,i)
  ntpoints = size(data,3);
  tmp_data = reshape(data(i,:,:),...
    [results.nchans,ntpoints]);
  j = 1;
  for c=1:results.nchans
    chan_col = parms.chan_colors{j};
    j = j + 1;
    if j>length(parms.chan_colors), j = 1; end;
    dmad = max(parms.data_mad(c),eps);
    % extract data for epoch
    ind_first = ind_spindle_rel(i,1);
    ind_last  = ind_spindle_rel(i,2);
    epoch_data = squeeze(data(i,c,ind_first:ind_last))';
    % remove linear trend
    epoch_data = detrend(epoch_data);
    % find local maxima with large enough peak to peak amplitudes
    [maxtab,mintab] = mmil_peakdet(epoch_data,...
                                   parms.min_peak_deflection*dmad);
    tmp_minima = get_minmax(mintab);
    tmp_maxima = get_minmax(maxtab);
    tmp_minima = get_deflection(tmp_minima,tmp_maxima);
    tmp_maxima = get_deflection(tmp_maxima,tmp_minima);
    % exclude peaks with deflection much smaller than max peak
    peak_vals = tmp_maxima.deflection;
    max_peak_val = max(peak_vals);
    peak_vals = peak_vals(peak_vals >= parms.min_peak_ratio*max_peak_val);
    npeaks = length(peak_vals);
    % calculate duration and amplitude
    dur = (1+diff(ind_spindle_rel(i,:)))/results.sfreq;
    amp = max(tmp_data(c,ind_spindle_rel(i,1):ind_spindle_rel(i,2)));
    x = -0.2*parms.time_max;
    y = 0.95*amp;
    text(x,y,sprintf('%d peaks, %0.1f Hz',npeaks,npeaks/dur),...
        'color',chan_col);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax = get_minmax(mtab)
  minmax = struct('latency',[],'amplitude',[],'deflection',[]);
  j = 1;
  for i=1:size(mtab,1)
    latency = mtab(i,1);
    amplitude = mtab(i,2);
    minmax.latency(j) = latency;
    minmax.amplitude(j) = amplitude;
    j = j + 1;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minmax1 = get_deflection(minmax1,minmax2)
  % calculate deflection for minima or maxima
  for i=1:length(minmax1.latency)
    latency = minmax1.latency(i);
    % find preceding minimum
    tmp = latency - minmax2.latency;
    tmp(tmp<0)=Inf;
    [tmp,ind1] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp1 = 0;
    else
      amp1 = minmax2.amplitude(ind1);
    end;

    % find subsequent maxima or minima
    tmp = minmax2.latency - latency;
    tmp(tmp<0)=Inf;
    [tmp,ind2] = min(tmp);
    if isempty(tmp) || isinf(tmp)
      amp2 = 0;
    else
      amp2 = minmax2.amplitude(ind2);
    end;

    minmax1.deflection(i) =...
      abs(minmax1.amplitude(i) - mean([amp1 amp2]));
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function annot_amp(data,parms,results,ind_spindle_rel,i)
  ntpoints = size(data,3);
  tmp_data = reshape(data(i,:,:),...
    [results.nchans,ntpoints]);
  j = 1;
  for c=1:results.nchans
    chan_col = parms.chan_colors{j};
    j = j + 1;
    if j>length(parms.chan_colors), j = 1; end;
    amp = max(tmp_data(c,ind_spindle_rel(i,1):ind_spindle_rel(i,2)));
    x = -0.2*parms.time_max;
    y = 0.8*amp;
    text(x,y,sprintf('amp = %0.2f',amp),'color',chan_col);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function annot_dur(data,parms,results,ind_spindle_rel,i)
  ntpoints = size(data,3);
  tmp_data = reshape(data(i,:,:),...
    [results.nchans,ntpoints]);
  max_val = 1.01*max(abs(tmp_data(:)));
  dur = (1+diff(ind_spindle_rel(i,:)))/results.sfreq;
  x = -0.2*parms.time_max;
  y = 0.8*max_val;
  text(x,y,sprintf('dur = %0.3f sec',dur),'color','k');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



