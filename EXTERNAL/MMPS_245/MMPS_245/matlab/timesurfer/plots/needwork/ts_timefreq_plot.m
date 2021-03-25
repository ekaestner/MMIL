function ts_timefreq_plot(timefreq_data, varargin)
%
%
% This is essentially a massive wrapper around a large set of FIELDTRIP
% plotting functions.
%
% This code basically converts parameters to fieldtrip CFG objects, then
% converts the timefreq object to a FIELDTRIP object, and finally does
% the necessary plotting.

  % single filename input
%  if (iscell(timefreq))
%    timefreq  = { timefreq };
%  end;
  parms.chans=[];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % plot results with multiplotTFR
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  cfg = [];
  cfg.baseline = [-80 -5]/1000;
  cfg.baselinetype = 'relative';
  %cfg.baselinetype = 'absolute';
  cfg.showlabels = 'yes';
  cfg.colorbar = 'yes';
  %cfg.zlim = [0 10^-24];
  cfg.zlim = [0 1.75];

  if (isempty(parms.chans))
    multiplotTFR(cfg, timefreq_data);
  else
    singleplotTFR(cfg,timefreq_data);
  end;
