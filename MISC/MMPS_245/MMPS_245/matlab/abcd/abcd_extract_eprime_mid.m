function abcd_extract_eprime_mid(fname,varargin)
%function abcd_extract_eprime_mid(fname,[options])
%
% Purpose: extract condition time courses
%   from eprime data files for MID task
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
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
% Created:  10/07/16 by Don Hagler
% Prev Mod: 09/13/17 by Don Hagler
% Last Mod: 09/26/17 by Dani Cornejo 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: based on abcd_mid_extract.m
%       provided by Mary Soules (mfield@med.umich.edu) 10/03/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname,varargin);

% create output directory
mmil_mkdir(parms.outdir);

% create struct with info for each event
[event_info,start_time,all_types] = get_event_info(parms);

%get behav data and write it to a csv
behav = get_behavioral_data(event_info,parms,all_types);  

% separate file for each combination of condnames and stimnames (events)
for i=1:parms.nconds
  type = parms.typenames{i};
  cond = parms.condnames{i};
  ind_type = find(strcmp(type,all_types));
  for j=1:parms.nstims
    stim = parms.stimnames{j};
    eventname = sprintf('%s_%s',cond,stim);
    % get times of stim onset and offset
    switch stim
      case 'pos_feedback'
        tmp_stim = 'feedback';
        % find events with positive feedback (acc = 1)
        acc = [event_info(ind_type).acc];
        ind_keep = find(acc==1);
        onset = [event_info(ind_type(ind_keep)).([tmp_stim '_onset'])];
        offset = [event_info(ind_type(ind_keep)).([tmp_stim '_offset'])];
        % find most recent start time for each event
        [ind_start,event_ref_time] = set_ref(onset,start_time);
        % create files for each scan
        write_files(eventname,onset,offset,ind_start,event_ref_time,parms)
      case 'neg_feedback'
        tmp_stim = 'feedback';
        % find events with negative feedback (acc = 0)
        acc = [event_info(ind_type).acc];
        ind_keep = find(acc==0);
        onset = [event_info(ind_type(ind_keep)).([tmp_stim '_onset'])];
        offset = [event_info(ind_type(ind_keep)).([tmp_stim '_offset'])];
        % find most recent start time for each event
        [ind_start,event_ref_time] = set_ref(onset,start_time);
        % create files for each scan
        write_files(eventname,onset,offset,ind_start,event_ref_time,parms)    
      otherwise
        tmp_stim = 'cue';
        ind_keep = [1:length(ind_type)];
        onset = [event_info(ind_type(ind_keep)).([tmp_stim '_onset'])]; 
        offset = [event_info(ind_type(ind_keep)).(['probe_offset'])]; 
        % find most recent start time for each event
        [ind_start,event_ref_time] = set_ref(onset,start_time); 
        % create files for each scan
        write_files(eventname,onset,offset,ind_start,event_ref_time,parms)    
    end; 
  end;
end; 

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
    'forceflag',false,[false true],...
    ...
    'colnames',  {'NARGUID','SessionDate','SessionTime','ExperimentName','MIDVERSION',...
                  'Block','SubTrial','Condition','PrepTime.OnsetTime','PrepTime.OffsetTime',...
                  'Cue.OnsetTime','Anticipation.OnsetTime','Probe.OnsetTime','Feedback.OnsetTime',...
                  'Cue.OffsetTime','Anticipation.OffsetTime','Probe.OffsetTime','Feedback.OffsetTime',...
                  'prbacc','prbrt','RunMoney',... 
                  'Anticipation.Duration','Cue.Duration','Probe.Duration','Probe.OffsetTime'},[],...
    'fieldnames',{'narguid','date','time','experiment','version',...
                  'run','trial','type','prep_onset','prep_offset',...
                  'cue_onset','antic_onset','probe_onset','feedback_onset',...
                  'cue_offset','antic_offset','probe_offset','feedback_offset',...
                  'acc','rt','money',...
                  'antic_dur','cue_dur','probe_dur','probe_offset'},[],...
    'typenames',{'SmallReward','LgReward','SmallPun','LgPun','Triangle'},[],...
    'condnames',{'small_reward','large_reward','small_loss','large_loss','neutral'},[],...
    'stimnames',{'antic','pos_feedback','neg_feedback'},[],...
  });

  parms.nconds = length(parms.condnames);
  parms.ntypes = length(parms.typenames);
  parms.nstims = length(parms.stimnames);
  parms.mindur = parms.minfrac*parms.TR;
  if parms.nconds ~= parms.ntypes
    error('condnames and typenames length mismatch');
  end;
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  [fdir,fstem,fext] = fileparts(parms.fname);
  if isempty(parms.outstem)
    parms.outstem = fstem;
  end;
  % calculate time at the end of each TR
  parms.TR_offset = linspace(parms.TR,parms.numTRs * parms.TR,parms.numTRs);
  parms.TR_onset = parms.TR_offset - parms.TR;

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
    [tmp,i_all,i_sel] = intersect(colnames,parms.colnames);
    [i_sel,i_sort] = sort(i_sel);
    i_all = i_all(i_sort);
    % create new matrix with replacement column labels
    vals = vals(:,i_all);
    tags = parms.fieldnames(i_sel);
    B = cat(1,tags,vals);
    % write to csv
    mmil_write_csv(fname_csv,B);
  end;
  event_info = mmil_csv2struct(fname_csv);
  % get start times
  ind_start = find(~cellfun(@isempty,{event_info.prep_onset}));
  start_time = [event_info(ind_start).prep_onset];
  % remove non-events
  all_types = {event_info.type};
  ind_events = find(~cellfun(@isempty,all_types));
  event_info = event_info(ind_events);
  all_types = {event_info.type};
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function behav = get_behavioral_data(event_info,parms,all_types) 

  behav = [];
  behav.('SubjID') = []; behav.('VisitID') = []; 
  
  behav.switch_flag = 0; 
  
  runs = [event_info.run]; 
  nruns = length(unique(runs));  
  
  for i=1:parms.nconds
  type = parms.typenames{i}; 
  cond = parms.condnames{i}; 
  ind_type = find(strcmp(type,all_types)); 
  
  type_counts_total = sprintf('%s_total',cond); 
  behav.(type_counts_total) = size(ind_type,2); 
  for j=1:2 %runs
    type_counts_run = sprintf('%s_run_%d_total',cond,j); 
    ind_type_run = [intersect(find(runs==j),ind_type)];    
    behav.(type_counts_run) = size(ind_type_run,2); 
  end 
  if nruns < 2
    behav.(type_counts_run) = NaN; 
  end 
  
  for j=1:parms.nstims
    stim = parms.stimnames{j}; 
    eventname = sprintf('%s_%s',cond,stim);
    
    acc = [event_info(ind_type).acc]; 
    rt = [event_info(ind_type).rt]; 
     
    switch stim
      case 'pos_feedback'
             
        ind_keep = find(acc==1); 
        type_counts_stim = sprintf('%s_counts',eventname);
        behav.(type_counts_stim) = size(ind_keep,2); 
        type_acc_stim = sprintf('%s_rate',eventname);
        behav.(type_acc_stim) = size(ind_keep,2)./size(ind_type,2); 
        type_rt_stim = sprintf('%s_rt',eventname);
        rt_keep = [event_info(ind_type(ind_keep)).(['rt'])];   
        behav.(type_rt_stim) = mean(rt_keep);  
        type_rt_stim = sprintf('%s_rt_std',eventname);
        behav.(type_rt_stim) = std(rt_keep);
        
        for j=1:2 %runs
          type_counts_run = sprintf('%s_run_%d_counts',eventname,j); 
          ind_type_run = [intersect(find(runs==j),ind_type)];  
          ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
          behav.(type_counts_run) = size(ind_keep_run,2); 
          type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
          behav.(type_acc_run) = size(ind_keep_run,2)./size(ind_type_run,2); 
          type_rt_run = sprintf('%s_run_%d_rt',eventname,j);
          rt_keep = [event_info(ind_keep_run).(['rt'])];    
          behav.(type_rt_run) = mean(rt_keep);
          type_rt_run = sprintf('%s_run_%d_rt_std',eventname,j);
          behav.(type_rt_run) = std(rt_keep);
        end
        if nruns < 2
          behav.(type_counts_run) = NaN; 
        end 
        
      case 'neg_feedback'
        ind_keep = find(acc==0); 
        type_counts_stim = sprintf('%s_counts',eventname);
        behav.(type_counts_stim) = size(ind_keep,2); 
        type_acc_stim = sprintf('%s_rate',eventname);
        behav.(type_acc_stim) = size(ind_keep,2)./size(ind_type,2); 
          
        for j=1:2 %runs
          type_counts_run = sprintf('%s_run_%d_counts',eventname,j); 
          ind_type_run = [intersect(find(runs==j),ind_type)];  
          ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
          behav.(type_counts_run) = size(ind_keep_run,2); 
          type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
          behav.(type_acc_run) = size(ind_keep_run,2)./size(ind_type_run,2);          
        end
        if nruns < 2
          behav.(type_counts_run) = NaN; 
        end 
    end
  end 
  end

  
  %collapse across small or large 
  collapse_cond = {'reward', 'loss'}; 
  collapse_ind_type.reward = [];
  collapse_ind_type.loss = [];
  
  for i=1:parms.nconds
    type = parms.typenames{i}; 
    cond = parms.condnames{i};  
    ind_type = find(strcmp(type,all_types)); 
    if findstr('reward', cond) > 0
      collapse_ind_type.reward = [collapse_ind_type.reward ind_type]; 
    elseif findstr('loss', cond) > 0
      collapse_ind_type.loss = [collapse_ind_type.loss ind_type]; 
    end
  end  
  
  for i=1:size(collapse_cond,2)
    cond = collapse_cond{i}; 
    type_counts_total = sprintf('%s_total',cond); 
    behav.(type_counts_total) = size(collapse_ind_type.(cond),2); 
    for j=1:2 %runs
     type_counts_run = sprintf('%s_run_%d_total',cond,j); 
     ind_type_run = [intersect(find(runs==j),collapse_ind_type.(cond))];    
     behav.(type_counts_run) = size(ind_type_run,2); 
    end 
    if nruns < 2
      behav.(type_counts_run) = NaN; 
    end 
   
    for j=1:parms.nstims
      stim = parms.stimnames{j}; 
      eventname = sprintf('%s_%s',cond,stim);   
    
      acc = [event_info(collapse_ind_type.(cond)).acc];  
      rt = [event_info(collapse_ind_type.(cond)).rt]; 
    
      switch stim
        case 'pos_feedback'
    
          ind_keep = find(acc==1); 
          ind_type = collapse_ind_type.(cond); 
          type_counts_stim = sprintf('%s_counts',eventname);
          behav.(type_counts_stim) = size(ind_keep,2); 
          type_acc_stim = sprintf('%s_rate',eventname);
          behav.(type_acc_stim) = size(ind_keep,2)./size(ind_type,2); 
          type_rt_stim = sprintf('%s_rt',eventname);
          rt_keep = [event_info(ind_type(ind_keep)).(['rt'])];   
          behav.(type_rt_stim) = mean(rt_keep); 
          type_rt_stim = sprintf('%s_rt_std',eventname);
          behav.(type_rt_stim) = std(rt_keep); 
          
          for j=1:2 %runs
            type_counts_run = sprintf('%s_run_%d_counts',eventname,j);  
            ind_type_run = [intersect(find(runs==j),ind_type)];  
            ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
            behav.(type_counts_run) = size(ind_keep_run,2); 
            type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = size(ind_keep_run,2)./size(ind_type_run,2); 
            type_rt_run = sprintf('%s_run_%d_rt',eventname,j);
            rt_keep = [event_info(ind_keep_run).(['rt'])];    
            behav.(type_rt_run) = mean(rt_keep);    
            type_rt_run = sprintf('%s_run_%d_rt_std',eventname,j);
            behav.(type_rt_run) = std(rt_keep);
          end 
          if nruns < 2
            behav.(type_counts_run) = NaN; 
          end 
        
        case 'neg_feedback'
    
          ind_keep = find(acc==0); 
          ind_type = collapse_ind_type.(cond); 
          type_counts_stim = sprintf('%s_counts',eventname);
          behav.(type_counts_stim) = size(ind_keep,2); 
          type_acc_stim = sprintf('%s_rate',eventname);
          behav.(type_acc_stim) = size(ind_keep,2)./size(ind_type,2);  
          
          for j=1:2 %runs
            type_counts_run = sprintf('%s_run_%d_counts',eventname,j);  
            ind_type_run = [intersect(find(runs==j),ind_type)];  
            ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
            behav.(type_counts_run) = size(ind_keep_run,2); 
            type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = size(ind_keep_run,2)./size(ind_type_run,2);         
          end
          if nruns < 2
            behav.(type_counts_run) = NaN; 
          end 
          
      end        
    end
  end
  
  %money 
  money = [event_info.money];
  behav.('money_total') = sum(money); 
  for j=1:2 %runs
    type_counts_run = sprintf('money_run_%d',j);  
    behav.(type_counts_run) = sum(money(find(runs==j)));   
  end   
  if nruns < 2
    behav.(type_counts_run) = NaN; 
  end 
  
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end; 

  
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



