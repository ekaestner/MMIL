function ts_summarize_resamp_results(results,varargin)
%function ts_summarize_resamp_results(results,[options])
%
% Purpose: summarize waveform resampling results as csv text file
%
% Required Input:
%   results: struct containing output from ts_resamp_wforms
%     may be struct array
%
% Optional Parameters:
%   'fstem_out': output file stem
%     multiple files may be created, e.g. for peaks, auc, onset, ranges
%     {default = 'resamp_results'}
%   'forceflag': [0|1] overwrite existing output
%     {default = 1}
%
% Created:  05/01/12 by Don Hagler
% Last Mod: 05/02/12 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;

parms = check_input(results,varargin);

if parms.peak_flag
  write_peak_amplitude(results,parms);
  write_peak_latency(results,parms);
end;

if parms.onset_flag
  write_onset_latency(results,parms);
end;

if parms.auc_flag
  write_auc(results,parms);
end;

if parms.ranges_flag
  write_ranges(results,parms,'ranges_sig');
  write_ranges(results,parms,'ranges_data');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(results,options)
  parms_filter = {...
    'fstem_out','resamp_results',[],...
    'forceflag',true,[false true],...
...
    'peak_flag',false,[false true],...
    'onset_flag',false,[false true],...
    'auc_flag',false,[false true],...
    'ranges_flag',false,[false true],...
    'diff_flag',false,[false true],...
  };
  parms = mmil_args2parms(options,parms_filter);
  
  % get info about how many rois and conds in results struct array
  [parms.nroi,parms.ncond] = size(results);

  % create row labels
  parms.nrow = parms.nroi*parms.ncond;
  parms.row_labels = cell(1,parms.nrow);
  k = 1;
  for r=1:parms.nroi
    for c=1:parms.ncond
      tmp_results = results(r,c);
      if isfield(tmp_results,'roi') & isfield(tmp_results,'cond')
        parms.row_labels{k} = sprintf('%s %s',tmp_results.roi,tmp_results.cond);
      elseif isfield(tmp_results,'roi') 
        parms.row_labels{k} = sprintf('%s cond%d',tmp_results.roi,c);
      elseif isfield(tmp_results,'cond') 
        parms.row_labels{k} = sprintf('roi%d %s',r,tmp_results.cond);
      else
        parms.row_labels{k} = sprintf('roi%d cond%d',r,c);
      end;      
      k = k + 1;
    end;
  end;

  % check for fields present in results
  if isfield(results,'peak')
    parms.peak_flag = 1;
    if isfield(results,'peak1'), parms.diff_flag = 1; end;
  end;
  if isfield(results,'auc')
    parms.auc_flag = 1;
    if isfield(results,'auc1'), parms.diff_flag = 1; end;
  end;
  if isfield(results,'onset')
    parms.onset_flag = 1;
    if isfield(results,'onset1'), parms.diff_flag = 1; end;
  end;
  if isfield(results,'ranges_sig'), parms.ranges_flag = 1; end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_peak_amplitude(results,parms)
  fname_out = [parms.fstem_out '_peak_amplitude.csv'];
  if ~exist(fname_out,'file') || parms.forceflag
    col_labels = {'mean','lo','hi'};
    if parms.diff_flag
      col_labels = cat(2,col_labels,...
        {'pval','mean1','lo1','hi1','mean2','lo2','hi2'});
    end;
    ncol = length(col_labels);
    data = zeros(parms.nrow,ncol);
    k = 1;
    for r=1:parms.nroi
      for c=1:parms.ncond
        tmp_results = results(r,c);
        j = 1;
        data(k,j) = tmp_results.peak.amplitude; j=j+1;
        data(k,j) = tmp_results.peak.amplitude_ci(1); j=j+1;
        data(k,j) = tmp_results.peak.amplitude_ci(2); j=j+1;
        if parms.diff_flag
          data(k,j) = tmp_results.peak.amplitude_pval; j=j+1;
          data(k,j) = tmp_results.peak1.amplitude; j=j+1;
          data(k,j) = tmp_results.peak1.amplitude_ci(1); j=j+1;
          data(k,j) = tmp_results.peak1.amplitude_ci(2); j=j+1;
          data(k,j) = tmp_results.peak2.amplitude; j=j+1;
          data(k,j) = tmp_results.peak2.amplitude_ci(1); j=j+1;
          data(k,j) = tmp_results.peak2.amplitude_ci(2); j=j+1;
        end;
        k = k + 1;
      end;
    end;
    mmil_write_csv(fname_out,data,...
      'col_labels',col_labels,...
      'row_labels',parms.row_labels,...
      'firstcol_label','roi x cond');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_peak_latency(results,parms)
  fname_out = [parms.fstem_out '_peak_latency.csv'];
  if ~exist(fname_out,'file') || parms.forceflag
    col_labels = {'mean','lo','hi'};
    if parms.diff_flag
      col_labels = cat(2,col_labels,...
        {'pval','mean1','lo1','hi1','mean2','lo2','hi2'});
    end;
    ncol = length(col_labels);
    data = zeros(parms.nrow,ncol);
    k = 1;
    for r=1:parms.nroi
      for c=1:parms.ncond
        tmp_results = results(r,c);
        j = 1;
        data(k,j) = tmp_results.peak.latency; j=j+1;
        data(k,j) = tmp_results.peak.latency_ci(1); j=j+1;
        data(k,j) = tmp_results.peak.latency_ci(2); j=j+1;
        if parms.diff_flag
          data(k,j) = tmp_results.peak.latency_pval; j=j+1;
          data(k,j) = tmp_results.peak1.latency; j=j+1;
          data(k,j) = tmp_results.peak1.latency_ci(1); j=j+1;
          data(k,j) = tmp_results.peak1.latency_ci(2); j=j+1;
          data(k,j) = tmp_results.peak2.latency; j=j+1;
          data(k,j) = tmp_results.peak2.latency_ci(1); j=j+1;
          data(k,j) = tmp_results.peak2.latency_ci(2); j=j+1;
        end;
        k = k + 1;
      end;
    end;
    mmil_write_csv(fname_out,data,...
      'col_labels',col_labels,...
      'row_labels',parms.row_labels,...
      'firstcol_label','roi x cond');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_onset_latency(results,parms)
  fname_out = [parms.fstem_out '_onset_latency.csv'];
  if ~exist(fname_out,'file') || parms.forceflag
    col_labels = {'mean','lo','hi'};
    if parms.diff_flag
      col_labels = cat(2,col_labels,...
        {'pval','mean1','lo1','hi1','mean2','lo2','hi2'});
    end;
    ncol = length(col_labels);
    data = zeros(parms.nrow,ncol);
    k = 1;
    for r=1:parms.nroi
      for c=1:parms.ncond
        tmp_results = results(r,c);
        j = 1;
        data(k,j) = tmp_results.onset.latency; j=j+1;
        data(k,j) = tmp_results.onset.latency_ci(1); j=j+1;
        data(k,j) = tmp_results.onset.latency_ci(2); j=j+1;
        if parms.diff_flag
          data(k,j) = tmp_results.onset.latency_pval; j=j+1;
          data(k,j) = tmp_results.onset1.latency; j=j+1;
          data(k,j) = tmp_results.onset1.latency_ci(1); j=j+1;
          data(k,j) = tmp_results.onset1.latency_ci(2); j=j+1;
          data(k,j) = tmp_results.onset2.latency; j=j+1;
          data(k,j) = tmp_results.onset2.latency_ci(1); j=j+1;
          data(k,j) = tmp_results.onset2.latency_ci(2); j=j+1;
        end;
        k = k + 1;
      end;
    end;
    mmil_write_csv(fname_out,data,...
      'col_labels',col_labels,...
      'row_labels',parms.row_labels,...
      'firstcol_label','roi x cond');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_auc(results,parms)
  fname_out = [parms.fstem_out '_auc.csv'];
  if ~exist(fname_out,'file') || parms.forceflag
    col_labels = {'mean','lo','hi'};
    if parms.diff_flag
      col_labels = cat(2,col_labels,...
        {'pval','mean1','lo1','hi1','mean2','lo2','hi2'});
    end;
    ncol = length(col_labels);
    data = zeros(parms.nrow,ncol);
    k = 1;
    for r=1:parms.nroi
      for c=1:parms.ncond
        tmp_results = results(r,c);
        j = 1;
        data(k,j) = tmp_results.auc.mean; j=j+1;
        data(k,j) = tmp_results.auc.ci(1); j=j+1;
        data(k,j) = tmp_results.auc.ci(2); j=j+1;
        if parms.diff_flag
          data(k,j) = tmp_results.auc.pval; j=j+1;
          data(k,j) = tmp_results.auc1.mean; j=j+1;
          data(k,j) = tmp_results.auc1.ci(1); j=j+1;
          data(k,j) = tmp_results.auc1.ci(2); j=j+1;
          data(k,j) = tmp_results.auc2.mean; j=j+1;
          data(k,j) = tmp_results.auc2.ci(1); j=j+1;
          data(k,j) = tmp_results.auc2.ci(2); j=j+1;
        end;
        k = k + 1;
      end;
    end;
    mmil_write_csv(fname_out,data,...
      'col_labels',col_labels,...
      'row_labels',parms.row_labels,...
      'firstcol_label','roi x cond');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_ranges(results,parms,ranges_type)
  if ~exist('ranges_type','var') || isempty(ranges_type)
    ranges_type = 'ranges_sig';
  end;
  fname_out = [parms.fstem_out '_' ranges_type '.csv'];
  if ~exist(fname_out,'file') || parms.forceflag
    max_nranges = 0;
    for r=1:parms.nroi
      for c=1:parms.ncond
        if isempty(results(r,c).(ranges_type))
          nranges = 0;
        else
          nranges = results(r,c).(ranges_type).nranges;
        end;
        if nranges > max_nranges
          max_nranges = nranges;
        end;
      end;
    end;

    if max_nranges == 0, return; end;
    col_labels = cell(1,max_nranges*2);
    j = 1;
    for n=1:max_nranges
      col_labels{j} = sprintf('range %d onset',n); j=j+1;
      col_labels{j} = sprintf('range %d offset',n); j=j+1;
      col_labels{j} = sprintf('range %d pval',n); j=j+1;
    end;
    ncol = length(col_labels);
    data = cell(parms.nrow,ncol);
    k = 1;
    for r=1:parms.nroi
      for c=1:parms.ncond
        tmp_results = results(r,c);
        j = 1;
        if isempty(tmp_results.(ranges_type))
          nranges = 0;
        else
          nranges = tmp_results.(ranges_type).nranges;
        end;
        for n=1:nranges
          data{k,j} = tmp_results.(ranges_type).t_onset(n); j=j+1;
          data{k,j} = tmp_results.(ranges_type).t_offset(n); j=j+1;
          data{k,j} = tmp_results.(ranges_type).pval(n); j=j+1;
        end;
        k = k + 1;
      end;
    end;
    mmil_write_csv(fname_out,data,...
      'col_labels',col_labels,...
      'row_labels',parms.row_labels,...
      'firstcol_label','roi x cond');
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

