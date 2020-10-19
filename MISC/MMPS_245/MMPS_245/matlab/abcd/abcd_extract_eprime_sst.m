function abcd_extract_eprime_sst(fname,varargin)
%function abcd_extract_eprime_sst(fname,[options])
%
% Purpose: extract condition time courses
%   from eprime data files for SST task
%
% Required input:
%   fname: name of input eprime file
%
% Optional input:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     if empty, will use filestem of fname
%     {default = []}
%   'TR': scan repetition time (s)
%     {default = 0.8}
%   'numTRs': number of repetitions
%     {default = 500}
%   'minfrac': minimum fraction of a TR to register an event
%     {default = 0.1}
%   'nskipTRs': number of TRs to remove from beginning of time course
%     {default = 0}
%   'switch_thresh': percent correct below which to swap button assignements
%     {default = 0.3}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  10/07/16 by Don Hagler
% Prev Mod: 09/13/17 by Don Hagler
% Last Mod: 09/26/17 by Dani Cornejo 
% 

%% todo: use forceflag for all output files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: Based on make-SST-refs.py provided by Michael Michael Riedel, Ph.D.
%       miriedel@fiu.edu. 
%       Using style from abcd_extract_eprime_mid.m provided by Donald
%       Hagler, Ph.D. dhagler@mail.ucsd.edu.  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check arguments 
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname,varargin);   

% create output directory
mmil_mkdir(parms.outdir);

% create structure with info for each event
% write field _events from original; no processing  
[event_info,start_time,all_types] = get_event_info(parms); 
if isempty(all_types), exit; end;  

% switch buttons if necesary 
[event_info,parms.switch_flag] = sst_switch(event_info,parms); 

% find types such that we have no-resp and go-late cases 
all_types_new = get_sst_types(event_info,parms);

% write with the correct types 
write_event_revised(parms,all_types_new);

% get behav data and write it to a csv
%% todo: use revised event_info instead of only orig or switched?
behav = get_behavioral_data(event_info,parms,all_types_new); 

% create files for each condition
for i=1:parms.nconds
  type  = parms.typenames{i};
  cond  = parms.condnames{i};
  pat = '\w_(?<Event>\w+)'; 
  get_event = regexp(cond,pat,'names');
  event = get_event.Event;
  %% todo: use revised event_info instead of only orig or switched?
  ind_type = find(strcmp(type,all_types_new));  
  switch event  
    case 'go'
      % get times of stim onset and offset
      parms.mindur = 0.08;  
      onset  = [event_info(ind_type).(['go_onset_time'])];    
      offset = [event_info(ind_type).(['go_offset_time'])];          
    case 'stop'
      % get times of stim onset and offset
      parms.mindur = 0;  
      onset  = [event_info(ind_type).(['ssd_onset_time'])];  
      offset = onset+2;         
  end
  % find most recent start time for each event
  [ind_start,event_ref_time] = set_ref(onset,start_time);  
  % create files for each scan
  write_files(cond,onset,offset,ind_start,event_ref_time,parms); 
end 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms = mmil_args2parms(options,{...
    'fname',fname,[],...
    ...
    'outdir',pwd,[],...
    'outstem',[],[],...
    'TR',0.8,[],...
    'numTRs',500,[],...
    'minfrac',0.1,[0,1],...
    'nskipTRs',0,[],...
    'switch_thresh',0.3,[0,0.5],...
    'forceflag',false,[false true],...
    ...
    'colnames',  {'NARGUID','SessionDate','SessionTime','ExperimentName','ExperimentVersion',...
                  'Subject','Procedure[Block]','Block','Trial','GetReady.RTTime', 'BeginFix.StartTime',...
                  'Go.RT', 'Fix.RT', 'StopSignal.RT','Go.Duration', 'SSDDur',...
                  'StopSignal.Duration','SSD.RT', 'Go.OnsetTime', 'Go.OffsetTime',  'Go.ACC' , 'Go.RESP',...
                  'Fix.RESP', 'SSD.OnsetTime','SSD.OffsetTime' , 'StopSignal.ACC', 'Go.CRESP',...
                  'StopSignal.RESP', 'SSD.ACC','TrialCode'},[],...
    'fieldnames',{'narguid','date','time','experiment','version',...
                  'subject','procedure_block','block','trial','getready_rttime','beginfix_st',...
                  'go_rt','fix_rt','stop_rt','go_dur','ssd_dur',...
                  'stop_dur','ssd_rt','go_onset_time','go_offset_time' , 'go_acc','go_resp',...
                  'fix_resp', 'ssd_onset_time','ssd_offset_time','stop_acc','go_cresp',...
                  'stop_resp', 'ssd_acc', 'type',},[],...
    'typenames',{'CorrectGo', 'CorrectLateGo', 'IncorrectGo', 'IncorrectLateGo','NoRespGo',... 
                 'CorrectStop', 'IncorrectStop','SSDStop'},[],...
    'condnames',{'correct_go','correctlate_go', 'incorrect_go', 'incorrectlate_go' ...
                  'noresp_go', 'correct_stop', 'incorrect_stop','ssd_stop'},[],...
  });

  % check conditions
  parms.nconds = length(parms.condnames);
  parms.ntypes = length(parms.typenames);
  parms.mindur = parms.minfrac*parms.TR;
  if parms.nconds ~= parms.ntypes
    error('condnames and typenames length mismatch');
  end;
  
  % check input file
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  [fdir,fstem,fext] = fileparts(parms.fname);
  if isempty(parms.outstem)
    parms.outstem = fstem;
  end;

  % compute time at the end of each TR
  parms.TR_offset = linspace(parms.TR,parms.numTRs * parms.TR,parms.numTRs);
  parms.TR_onset  = parms.TR_offset - parms.TR;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,start_time,all_types] = get_event_info(parms)
  event_info=[]; start_time=[]; all_types=[]; 

  % write event info to file
  fname_csv = sprintf('%s/%s_events.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv,'file') || parms.forceflag
    %% todo: use new function for robust reading of eprime files
    % check for UTF-16, convert to ASCII if necessary
    fname = abcd_check_eprime_encoding(parms.fname,parms.outdir,parms.forceflag); 
    % read input file, allowing either comma or tab delimited
    A = mmil_readtext(fname,'[,\t]');
    % remove usual first row unless it is missing
    if ~strcmp(char(A(1,1)),'ExperimentName')
      A = A(2:end,:);
    end;
    % select columns of interest
    colnames = A(1,:);
    vals = A(2:end,:);
    [~,i_all,i_sel] = intersect(colnames,parms.colnames); 
    [i_sel,i_sort] = sort(i_sel);
    i_all = i_all(i_sort);
    % create new matrix with replacement column labels
    vals = vals(:,i_all);
    tags = parms.fieldnames(i_sel);
    B = cat(1,tags,vals);
    % write to csv
    mmil_write_csv(fname_csv,B);
  end;

  % read event info into struct array
  event_info = mmil_csv2struct(fname_csv);    
  % get start times
  ind_start = find(~cellfun(@isempty,{event_info.beginfix_st})); 
  start_time = [event_info(ind_start).beginfix_st]; 
  % remove non-events
  all_types = {event_info.type}; 
  ind_events = find(~cellfun(@isempty,all_types)); 
  event_info = event_info(ind_events); 
  all_types = {event_info.type};    
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [event_info,switch_flag] = sst_switch(event_info,parms)
  switch_flag = 0;
  all_types = {event_info.type};
  correct_total = size(cell2mat({event_info.go_cresp}),2); 
  correct = 0; 
  for i=1:length(all_types)
    if strcmp(all_types{i},'CorrectGo')
      correct = correct+1; 
    end;
  end;
  accuracy = 100*correct/correct_total;
  switch_pass = 0;
  fprintf('%s: accuracy = %0.1f%%\n',mfilename,accuracy);
  if accuracy < 100*parms.switch_thresh
    fprintf('%s: accuracy < %0.1f%%, switching button responses\n',...
      mfilename,100*parms.switch_thresh);
    switch_pass = 1;
  end
  if switch_pass 
    for i=1:length({event_info.go_cresp})
      if event_info(i).go_resp == 1
        event_info(i).go_resp=2; 
      elseif event_info(i).go_resp == 2
        event_info(i).go_resp=1;
      end 
      if event_info(i).fix_resp == 1
        event_info(i).fix_resp=2; 
      elseif event_info(i).fix_resp == 2
        event_info(i).fix_resp=1;
      end 
      if event_info(i).stop_resp == 1
        event_info(i).stop_resp=2; 
      elseif event_info(i).stop_resp == 2
        event_info(i).stop_resp=1;
      end     
      if event_info(i).go_acc == 1
        event_info(i).go_acc=0; 
      end 
      if event_info(i).stop_acc == 1
        event_info(i).stop_acc=0; 
      end 
      if event_info(i).ssd_acc == 1
        event_info(i).ssd_acc=0; 
      end 
    end 
    % write to csv
    fname_csv_out  = sprintf('%s/%s_events_switched.csv',...
      parms.outdir,parms.outstem); 
    mmil_struct2csv(event_info,fname_csv_out) 
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% todo: what is the purpose of this?
function write_event_revised(parms,all_types_new)
  % output file name 
  if parms.switch_flag
    fname_csv_in  = sprintf('%s/%s_events_switched.csv',parms.outdir,parms.outstem);
  else
    fname_csv_in  = sprintf('%s/%s_events.csv',parms.outdir,parms.outstem);
  end;
  fname_csv_out = sprintf('%s/%s_events_revised.csv',parms.outdir,parms.outstem); 
  % output file writing 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    % event info 
    event_info = mmil_csv2struct(fname_csv_in);    
    % find events
    all_types = {event_info.type}; 
    ind_events = find(~cellfun(@isempty,all_types));  
    for i=1:length(ind_events)
      event_info(ind_events(i)).type = char(all_types_new(i));     
    end 
    % write to csv
    mmil_struct2csv(event_info,fname_csv_out);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all_types_new = get_sst_types(event_info,parms)
  all_types_new = {event_info.type};  
  % correct_go
  go_cresp = {event_info.go_cresp};
  go_resp = {event_info.go_resp};  
  correct_go_ind = [];
  for i=1:size(go_cresp,2)  
    if go_resp{i} == go_cresp{i}
      correct_go_ind = [correct_go_ind i]; 
      all_types_new{i} = 'CorrectGo';
    end 
  end
  % correct_go_late 
  fix_resp = {event_info.fix_resp}; 
  correct_go_late_ind = [];
  for i=1:size(go_cresp,2)  
    if fix_resp{i} == go_cresp{i} 
      if go_resp{i} == go_cresp{i}
      else 
        correct_go_late_ind = [correct_go_late_ind i]; 
        all_types_new{i} = 'CorrectLateGo';
     end 
    end 
  end
  % no_resp_go 
  go_event_ind = find(~cellfun(@isempty,go_cresp)); 
  noresp_ind = find(cellfun(@isempty,go_resp)); 
  nofix_ind = find(cellfun(@isempty,fix_resp));
  noresp_nofix_ind = [intersect(noresp_ind,nofix_ind)];  
  no_resp_go_ind = [intersect(noresp_nofix_ind,go_event_ind)];   
  for i=no_resp_go_ind
    all_types_new{i} = 'NoRespGo';
  end
  % incorrect_go
  incorrect_go_ind = [];
  incorrect_go_late_ind = [];
  for i=1:size(go_cresp,2)  
    if go_resp{i} ~= go_cresp{i}
      incorrect_go_ind = [incorrect_go_ind i]; 
      all_types_new{i} = 'IncorrectGo'; %legit 1~=2
    elseif go_resp{i} == go_cresp{i}
    else 
      if fix_resp{i} ~= go_cresp{i}
        incorrect_go_late_ind = [incorrect_go_late_ind i]; 
        all_types_new{i} = 'IncorrectLateGo';
      end 
    end 
  end
  % correct & incorrect _stop 
  ssd_rt = {event_info.ssd_rt};
  ssd_dur = {event_info.ssd_dur};
  stop_resp = {event_info.stop_resp};  
  fix_resp = {event_info.fix_resp}; 
  incorrect_stop_ind = [];
  correct_stop_ind = []; 
  ssd_stop_ind = []; 
  for i=1:size(ssd_dur,2) 
    if ssd_dur{i} > 0  % only stop trials  
      if stop_resp{i} > 0
        incorrect_stop_ind = [incorrect_stop_ind i]; 
        all_types_new{i} = 'IncorrectStop';
      elseif fix_resp{i} > 0
        incorrect_stop_ind = [incorrect_stop_ind i]; 
        all_types_new{i} = 'IncorrectStop';     
      else 
        correct_stop_ind = [correct_stop_ind i];
        all_types_new{i} = 'CorrectStop';
      end
      if ssd_rt{i} > 0
        % stop too early 
        ssd_stop_ind = [ssd_stop_ind i];   
        all_types_new{i} = 'SSDStop';
      end
    end
  end
  % write to txt file %% todo: what this the purpose of this?
  fname_txt = sprintf('%s/%s_events_revised.txt',parms.outdir,parms.outstem);
  if ~exist(fname_txt,'file') || parms.forceflag
    fid = fopen(fname_txt,'wt');
    for i=1:length(all_types_new)
      fprintf(fid,'%s\n',all_types_new{i});
    end
    fclose(fid);
  end 
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_start,event_ref_time] = set_ref(onset,start_time)
  if ~isempty(onset)
    d = bsxfun(@minus,onset,start_time');
    d(d<0) = Inf;
    [d,ind_start] = min(d,[],1); 
    event_ref_time  = start_time(ind_start);
  else
    ind_start = start_time;
    event_ref_time = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_files(eventname,onset,offset,ind_start,event_ref_time,parms)
  nscans = length(unique(ind_start));   
  for s=1:2
    ind_scan = find(ind_start == s); 
    fname_1D = sprintf('%s/%s_scan%d_%s.1D',parms.outdir,parms.outstem,s,eventname);
    fname_txt = sprintf('%s/%s_scan%d_%s.txt',parms.outdir,parms.outstem,s,eventname);
    if ~exist(fname_1D,'file') || ~exist(fname_txt,'file') || parms.forceflag
      if ~isempty(onset(ind_scan))
        % subtract start time and convert to seconds
        rel_onset  = (onset(ind_scan) - event_ref_time(ind_scan))/1000;   
        rel_offset = (offset(ind_scan) - event_ref_time(ind_scan))/1000; 
        % loop over events
        TR_mask = zeros(parms.numTRs,1); 
        for k=1:length(rel_onset)  
          ind_rep = find(rel_onset(k)<(parms.TR_offset-parms.mindur) &...
                        rel_offset(k)>(parms.TR_onset+parms.mindur));  
          TR_mask(ind_rep) = 1;
        end; 
        % calculate event times from TR_onset
        event_times = parms.TR_onset(find(TR_mask)); 
        % convert time units to multiples of TRs
        event_TRs = event_times/parms.TR; 
        % remove  nskipTRs
        if parms.nskipTRs>0
          TR_mask = TR_mask(parms.nskipTRs+1:end);
          event_times = event_times(parms.nskipTRs+1:end) - parms.nskipTRs*parms.TR;
        end;
        % write to 1D file
        fid = fopen(fname_1D,'wt');
        if fid<0, error('failed to open %s for writing',fname_1D); end;
        for k=1:length(TR_mask)
          fprintf(fid,'%d\n',TR_mask(k));
        end;
        fclose(fid);
        % write to txt file with stim times in seconds
        fid = fopen(fname_txt,'wt');
        if fid<0, error('failed to open %s for writing',fname_txt); end;
        fprintf(fid,'%s\n',sprintf('%0.2f ',event_TRs*parms.TR));
        fclose(fid);  
      else
        TR_mask = zeros(parms.numTRs,1);
        % write to 0s 1D file
        fid = fopen(fname_1D,'wt');
        if fid<0, error('failed to open %s for writing',fname_1D); end;
        for k=1:length(TR_mask)
          fprintf(fid,'%d\n',TR_mask(k));
        end;
        fclose(fid);
        % write * txt file
        fid = fopen(fname_txt,'wt');
        if fid<0, error('failed to open %s for writing',fname_txt); end;
        fprintf(fid,'*\n');
        fclose(fid);   
      end;
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function behav = get_behavioral_data(event_info,parms,all_types_new) 

  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.switch_flag = parms.switch_flag;

  blocks = [event_info.block]; 
  trials = [event_info.trial]; 
  % determine version of task
  n_blocks = size(unique(blocks),2);
  n_trials = size(unique(trials),2);
  if n_trials > 2, 
    blocks = [event_info.block]; 
  else 
    blocks = [event_info.trial];
  end
  
  nruns = length(unique(blocks));  
  
  % event durations
  go_dur = {event_info.go_dur}; 
  time_go_dur = [event_info.go_dur];
  time_go_dur = time_go_dur(1); 
  stop_dur = {event_info.stop_dur}; 
  go_event_ind = find(~cellfun(@isempty,go_dur)); 
  stop_event_ind = find(~cellfun(@isempty,stop_dur));
  % total numbers of trials    
  behav.('Go_total') = size(go_event_ind,2);
  behav.('Stop_total') = size(stop_event_ind,2);
  behav.('Total') = size(go_event_ind,2)+size(stop_event_ind,2);
  % number of trials per run    
  for j=1:2 %runs 
    % go
    type_counts_run_go = sprintf('Go_total_run_%d',j);
    ind_type_run = [intersect(find(blocks==j),go_event_ind)];   
    behav.(type_counts_run_go) = size(ind_type_run,2); 
    % stop
    type_counts_run_stop = sprintf('Stop_total_run_%d',j);
    ind_type_run = [intersect(find(blocks==j),stop_event_ind)];   
    behav.(type_counts_run_stop) = size(ind_type_run,2); 
  end
  if nruns < 2
    behav.(type_counts_run_go) = NaN; 
    behav.(type_counts_run_stop) = NaN;
  end 
  
  % results for each condition
  for i=1:parms.nconds
    % total counts per type
    type  = parms.typenames{i};         
    type_counts_total = sprintf('%s_counts_total',type);  
    ind_type = find(strcmp(type,all_types_new)); 
    behav.(type_counts_total) = size(ind_type,2); 
    % total counts per type per run  
    for j=1:2 % runs        
      type_counts_run = sprintf('%s_counts_run_%d',type,j); 
      ind_type_run = [intersect(find(blocks==j),ind_type)];  
      behav.(type_counts_run) = size(ind_type_run,2); 
    end 
    if nruns < 2
      behav.(type_counts_run) = NaN; 
    end 
    
    % selecting go or stop event 
    cond  = parms.condnames{i};  
    pat = '\w_(?<Event>\w+)'; 
    get_event = regexp(cond,pat,'names');
    event = get_event.Event;    
    switch event  
      case 'go'
        event_ind = go_event_ind; 
        event_rt = {event_info.('go_rt')}; 
        if strcmp(cond,'correctlate_go')
           event_rt_fix = {event_info.('fix_rt')};
           event_rt = cellfun(@(x) x+time_go_dur, event_rt_fix,'UniformOutput',false);
        end 
        if strcmp(cond,'incorrectlate_go')
           event_rt_fix = {event_info.('fix_rt')};
           event_rt = cellfun(@(x) x+time_go_dur, event_rt_fix,'UniformOutput',false);
        end 
      case 'stop'
        event_ind = stop_event_ind;
        event_rt = {event_info.('stop_rt')};
    end 
    % total percent per type 
    type_percent_total = sprintf('%s_percent_total',type);  
    behav.(type_percent_total) = size(ind_type,2)./size(event_ind,2);  
    % percent per type per run  
    for j=1:2 % runs 
      type_counts_run = sprintf('%s_percent_run_%d',type,j); 
      ind_type_run = [intersect(find(blocks==j),ind_type)];
      ind_type_percent_run = [intersect(find(blocks==j),event_ind)]; 
      behav.(type_counts_run) = size(ind_type_run,2)./size(ind_type_percent_run,2);  
    end
    % total rt per type 
    type_rt_total = sprintf('%s_rt_total',type); 
    all_rt = cell2mat(event_rt(ind_type)); 
    behav.(type_rt_total) = mean(all_rt(all_rt~=0));  
    type_rt_std_total = sprintf('%s_rt_std_total',type);
    behav.(type_rt_std_total) = std(all_rt(all_rt~=0)); 
    
    % total rt per type per run 
    for j=1:2 % runs 
      type_rt_runs = sprintf('%s_rt_run_%d',type,j);
      ind_type_run = [intersect(find(blocks==j),ind_type)];
      all_rt = cell2mat(event_rt(ind_type_run));  
      behav.(type_rt_runs) = mean(all_rt(all_rt~=0)); 
      type_rt_std_runs = sprintf('%s_rt_std_run_%d',type,j);
      behav.(type_rt_std_runs) = std(all_rt(all_rt~=0));
    end     
  end 
  % ssd and ssrt 
  ssd = {event_info.ssd_dur}; 
  behav.('SSD_mean_total')= mean(cell2mat(ssd));   
  behav.('SSRT_mean_total') = behav.CorrectGo_rt_total-behav.SSD_mean_total; 
  % rm meaningless fields 
  rmfields_names = {'NoRespGo_rt_total','NoRespGo_rt_run_1','NoRespGo_rt_run_2', ... 
                    'NoRespGo_rt_std_total','NoRespGo_rt_std_run_1','NoRespGo_rt_std_run_2', ... 
                    'CorrectStop_rt_total','CorrectStop_rt_run_1', ...
                    'CorrectStop_rt_std_total','CorrectStop_rt_std_run_1', ...
                    'CorrectStop_rt_run_2','SSDStop_rt_total','SSDStop_rt_run_1',...
                    'CorrectStop_rt_std_run_2','SSDStop_rt_std_total','SSDStop_rt_std_run_1',...
                    'SSDStop_rt_run_2','SSDStop_rt_std_run_2'}; 
                  
  behav = rmfield(behav,rmfields_names); 
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

