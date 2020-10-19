function [data] = ts_session(parmfile,varargin)
% Last Mod: 09/15/12 by Don Hagler
%

if nargin == 0
  fprintf('\nERROR: %s requires at least one input.\n',mfilename);
  disp('Function help and wiki documentation are in progress.');
  disp('Call Sequence: [data] = ts_session(parmfile,varargin);');
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STARTUP 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart = tic;
inparms = mmil_args2parms(varargin,...
						{'verbose',0,{0,1},...
             'gui_study',[],[],...
						},false);
          
% load common function parameters
if ischar(parmfile) && strfind(parmfile,'.csv')
  % parse csv file
  if inparms.verbose
    fprintf('%s: reading parameter file: %s\n',mfilename,parmfile);
  end
  base = load_parms(parmfile);
elseif isstruct(parmfile)
  % user supplied parameter structure
  base = parmfile;
elseif isstruct(inparms.gui_study)
  % parse STUDY from GUI (timesurfer)
  if inparms.verbose
    fprintf('%s: running timesurfer study\n',mfilename);
  end
  funs = inparms.gui_study.funlist;
  glob = find(ismember(funs,'global'));
  base = [];
  for k = 1:length(funs)
    if strcmpi(funs{k},'global'), continue; end
    tmp = inparms.gui_study;
    tmp.protocol = tmp.protocol([glob k]);
    tmp.funlist  = tmp.funlist([glob k]);
    tmp.funindex = tmp.funindex([glob k]);
    out = load_parms('gui_study',tmp);
    if isempty(base)
      base = out;
    else
      base.function(end+1) = out.function;
    end
  end
  if isempty(base), return; end
  clear funs glob tmp out k
%   base = load_parms('gui_study',inparms.gui_study);
end

% return if there are no functions to call
if ~isfield(base,'function')
  fprintf('WARNING: Exiting session - no functions specified.\n');
  return;
end

% load subject-specific parameters
if issubfield(base,'global.subjectfile') && ~isempty(base.global.subjectfile)
  % process multiple subjects
  if inparms.verbose
    fprintf('%s: reading subject file: %s\n',mfilename,base.global.subjectfile);
  end  
  subj = load_subject_parameters(base);
else
  % process one subject
  subj = base;
end

data = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over subjects
for s = 1:length(subj)
  allparms = subj(s);
  subj_tstart = tic;
  % initialize subject-specific session
  session = []; orig_session = [];
  if issubfield(allparms,'global.sessionfile')
    sessionfile = allparms.global.sessionfile;
  else
    sessionfile = '';
  end
  if ~issubfield(allparms,'global.rootoutdir') || ~ischar(allparms.global.rootoutdir)
    allparms.global.rootoutdir = pwd;
  end
  if ~issubfield(allparms,'global.rootindir') || ~ischar(allparms.global.rootindir)
    allparms.global.rootindir = allparms.global.rootoutdir;
  end
  if ~exist(sessionfile,'file'), sessionfile = sprintf('%s/session.mat',allparms.global.rootoutdir); end
  if ~exist(sessionfile,'file'), sessionfile = sprintf('%s/session.mat',allparms.global.rootindir);  end
  if ~exist(sessionfile,'file'), sessionfile = sprintf('%s/session.mat',pwd); end  
  if  exist(sessionfile,'file'), load(sessionfile); orig_session=session; else sessionfile = '';   end
  logfid = 1;
  if ~isempty(sessionfile)
    [fpath,fname] = fileparts(sessionfile);
    logfile = [fpath '/session.log'];
  elseif issubfield(allparms,'global.rootoutdir')
    logfile = [allparms.global.rootoutdir '/session.log'];
  else
    logfile = [pwd '/session.log'];
  end
  if issubfield(allparms,'global.logfid') && ~isempty(allparms.global.logfid) && ~isequal(allparms.global.logfid,1)
    logfile = [];
    logfid  = allparms.global.logfid;
  end
    
  % loop over functions
  for f = 1:length(allparms.function)
    % initialize function
    fun            = allparms.function(f).name;
    funcall        = allparms.function(f).funcall;
    itype          = allparms.function(f).itype;
    otype          = allparms.function(f).otype;
    run_flag       = allparms.function(f).process_flags.run_flag;
    cluster_flag   = allparms.function(f).process_flags.cluster_flag;
    script_flag    = allparms.function(f).process_flags.script_flag; 
    load_flag      = allparms.function(f).spec_flags.load_flag;
    save_flag      = allparms.function(f).spec_flags.save_flag;
    parms_orig     = checkparms(allparms.function(f).parms);
    parms_orig.function = fun;
    parms_orig.itype    = itype; % parms.indata   = itype;
    parms_orig.otype    = otype; % parms.outdata  = otype;
    
    % rootoutdir
    if ~isfield(parms_orig,'rootoutdir') || isempty(parms_orig.rootoutdir)
      parms_orig.rootoutdir = allparms.global.rootoutdir;
    end
    % rootindir
    if ~isfield(parms_orig,'rootindir') || isempty(parms_orig.rootindir)
      parms_orig.rootindir = allparms.global.rootindir;
    end    
    % set the number of function calls w.r.t. the loop parameter
    if isempty(parms_orig.loop)
      ncalls = 1; 
			parms_orig.loop = [];
    else
      ncalls = length(parms_orig.(parms_orig.loop{1}));
    end

    % loop over function calls
    for i = 1:ncalls
      parms = parms_orig;
      if ~isfield(parms,'logfile') || ~ischar(parms.logfile)
        parms.logfile = logfile;
      end
      if ~isfield(parms,'logfid') || isempty(parms.logfid)
        parms.logfid = logfid;
      end
      if i==1 && f==1 && s==1
        % mmil_logstr(parms,'\nData analysis started at %s\n',datestr(now));
      end
      % set the iteration-specific parameter value
        for n = 1:length(parms.loop)
          if ~iscell(parms.(parms.loop{n})), continue; end
          parms.(parms.loop{n}) = parms.(parms.loop{n}){i};
        end
        
      % don't load data if in it's already in memory or if cluster computing
      if (f > 1 && ischar(itype) && exist(itype,'var')) % || cluster_flag
        load_flag  = 0;
        found_flag = 1;
      else
        found_flag = 0;
      end
%% determine the input datafiles     
      datafile = {};
      % check session for datafiles
      if ~found_flag && isfield(session,'function_id') && isfield(parms,'input_id') && isfield(session,'filename')
        tempfile = {};
        for k = 1:length(parms.input_id)
          this_id = parms.input_id(k);
          ind = find([session.function_id] == this_id);
          if ~isempty(ind)
            temptemp = {};
            for ix = 1:length(ind)
              sessionfiles  = {};              
              sessionfiles_ = session(ind(ix)).filename;
              if isempty(sessionfiles_), continue; end
              for ixs = 1:length(sessionfiles_)
                if exist(sessionfiles_{ixs},'file')
                  sessionfiles = {sessionfiles{:} sessionfiles_{ixs}};
                end
              end
              temptemp = {temptemp{:},sessionfiles{:}};
              clear sessionfiles
              clear sessionfiles_
            end
            if isempty(temptemp) || ~exist(temptemp{1},'file'), continue; end
            tempfile = {tempfile{:},temptemp{:}};
            clear temptemp
          end
          clear ind this_id
        end
        if ~isempty(tempfile)
          datafile = unique(tempfile); found_flag = 1; 
        end
        clear tempfile
      end      
      
      % check parms for datafiles 
      if ~found_flag && isfield(parms,'datafile')
        if ~iscell(parms.datafile), parms.datafile = {parms.datafile}; end
        nfiles = length(parms.datafile);
        found  = zeros(1,nfiles);
        for k  = 1:nfiles
          datafile   = parms.datafile{k}; if exist(datafile,'file')==2, found(k)=1; end
          if ~found(k) && isfield(parms,'datapath')
            datafile = [parms.datapath '/' parms.datafile{k}];   
            if exist(datafile,'file')==2, parms.datafile{k}=datafile;   found(k)=1; end
          end
          if ~found(k) && isfield(parms,'rootindir')
            datafile = [parms.rootindir '/' parms.datafile{k}];  
            if exist(datafile,'file')==2, parms.datafile{k}=datafile;   found(k)=1; end
          end
          if ~found(k) && isfield(parms,'rootindir') && isfield(parms,'inpath')
            datafile = [parms.rootindir '/' parms.inpath '/' parms.datafile{k}];  
            if exist(datafile,'file')==2, parms.datafile{k}=datafile;   found(k)=1; end
          end          
          if ~found(k) && isfield(parms,'rootoutdir')
            datafile = [parms.rootoutdir '/' parms.datafile{k}]; 
            if exist(datafile,'file')==2, parms.datafile{k}=datafile;   found(k)=1; end
          end      
          if ~found(k)
            datafile = [pwd '/' parms.datafile{k}]; 
            if exist(datafile,'file')==2, parms.datafile{k}=datafile;   found(k)=1; end
          end                 
        end
        if any(found)
          found_flag = 1;
          datafile   = parms.datafile(found > 0);
        else
          datafile   = {};
        end
        clear nfiles        
      end
      if ~found_flag && (~isfield(parms,'datapath') || ~exist(parms.datapath,'dir'))
        if all(isfield(parms,{'rootindir','inpath'}))
          parms.datapath = fullfile(parms.rootindir,parms.inpath);
        else
          parms.datapath = parms.rootindir;
        end
      end
      if ~found_flag && isfield(parms,'datapath') && exist(parms.datapath,'dir') % datapath
        DIR      = dir(parms.datapath);
        files    = {DIR.name};
        files    = files(~[DIR.isdir]);
        for k = 1:length(files)
          datafile{end+1} = fullfile(parms.datapath,files{k});
        end
        clear DIR files
        found_flag = 1;
%       elseif ~found_flag && isfield(parms,'rootindir') && exist(parms.rootindir) % rootindir
%         DIR      = dir(parms.rootindir);
%         files    = {DIR.name};
%         files    = files(~[DIR.isdir]);
%         for k = 1:length(files)
%           datafile{end+1} = fullfile(parms.rootindir,files{k});
%         end
%         clear DIR files
%         found_flag = 1;        
      end
      
      % check hdr for datafiles
      if ~found_flag % hdr...
      end
      
      if ~iscell(datafile), datafile = {datafile}; end
      
      % apply datastring constraints (this is ugly)
      if length(datafile) > 1 && isfield(parms,'datastring') && ~isempty(parms.datastring)
        if ~iscell(parms.datastring), parms.datastring = {parms.datastring}; end
        for pp = 1:length(parms.datastring)
          strloc   = strfind(datafile,parms.datastring{pp});
          keep     = ~cellfun('isempty',strloc);
          datafile = datafile(keep);
        end
        clear strloc keep
      end        
      if length(datafile) > 1 && isfield(parms,'dataskip') && ~isempty(parms.dataskip)
        if ~iscell(parms.dataskip), parms.dataskip = {parms.dataskip}; end
        for pp = 1:length(parms.dataskip)
          strloc   = strfind(datafile,parms.dataskip{pp});
          keep     = cellfun('isempty',strloc);
          datafile = datafile(keep);
        end
        clear strloc keep
      end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      % load_flag
        % 0 - no (do nothing)      
        % 1 - find datafiles
        % 2 - find and load datafiles   
      if     load_flag == 0
        % do nothing
      elseif load_flag == 1 || load_flag == 2
        % find datafiles
        if ~isempty(datafile)
          parms.datafile = datafile;
        else
          % find_datafiles()
        end
      end
      if load_flag == 2 && run_flag
        % load datafiles
        hdr.parms.filename = {};
        hdr.parms.filename = datafile;
%% MOVE THIS TO ts_checkdata_header (e.g., the stat_data type handling)
% %         for f = 1:length(datafile)
% %           tmpdata = load(datafile{f});
% %           fdnames = fieldnames(tmpdata);
% %           if any(strcmp(itype,fdnames))
% %             % this matfile contains itype
% %             hdr.parms.filename = {hdr.parms.filename{:},datafile{f}};
% %           elseif isfield(parms,fdnames{1}) && ~isstruct(parms.(fdnames{1}))
% %             % keep if this is a function parameter that needs to be defined
% %             parms.(fdnames{1}) = tmpdata.(fdnames{1});
% %           end
% %           clear tmpdata fdnames
% %         end
%%
        eval(sprintf('%s = hdr;',itype));
      end
      
      % get output filename & check whether it already exists
      if save_flag == 1 && (~isfield(parms,'filename') || isempty(parms.filename))
        args = mmil_parms2args(parms);
        parms.filename = ts_create_filename(fun,args{:});
        clear args;
%         parms   = ts_make_output_filename(parms);
        if ~iscell(parms.filename)
          parms.filename = {parms.filename}; 
        end
        overwrite = 0; 
        if isfield(parms,'overwrite'), overwrite = parms.overwrite; end
        if exist(parms.filename{1},'file') && ~overwrite
          fprintf('%s: skipping %s, output file already exists: %s\n',mfilename,fun,parms.filename{1});
          % mmil_logstr(parms,'%s: skipping %s, output file already exists: %s\n',mfilename,fun,parms.filename{1});
          continue;
        end
        pathstr = fileparts(parms.filename{1});
        if ~exist(pathstr,'dir')
%           fprintf('making output directory: %s\n',pathstr);
          % mmil_logstr(parms,'making output directory: %s\n',pathstr);
          unix(['mkdir -p ' pathstr]);
        end        
      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      % process function
        % cluster_flag (cluster overrides local machine)
        % run_flag
        % script_flag      
      if cluster_flag
        % write script to batchdirs
        if ~exist('batchdir','var')
          funscripts = {}; 
          [jnk user] = unix('whoami');
          user       = user(1:end-1);
          if isfield(parms,'rootindir') && ischar(parms.rootindir)
%             strdir = strrep(parms.rootindir,'/','_');
%             strdir = strrep(strdir,' ','_');
            [jnk strdir] = fileparts(parms.rootindir);
            batchdir = sprintf('%s_%s_%s',fun,datestr(now,30),strdir);
          else
            batchdir = sprintf('%s_%s',fun,datestr(now,30));
          end
        end
        outpath  = sprintf('/home/%s/batchdirs/%s',user,batchdir);
        thisscript = sprintf('run_%s_%g',fun,i);
        outfile    = sprintf('%s/%s.m',outpath,thisscript);
        funscripts{end+1} = thisscript;
        if isfield(parms,'filename')
          ts_write_mscript(parms,'function',fun,'cmdstr',funcall,'outfile',outfile,'hdr_flag',1,'filename',parms.filename{1})
        else
          ts_write_mscript(parms,'function',fun,'cmdstr',funcall,'outfile',outfile,'hdr_flag',1)
        end
        clear thisscript outfile;        
      end
      if run_flag
        % execute function call
        try
          eval(funcall);
        catch
          % save subject-specific session
            if ~isempty(session)
              session = close_session(parms,inparms,session,orig_session,sessionfile,allparms);
            end
            fprintf('%s: function failed: %s\n   Subject %g/%g, Function %g/%g, Call %g/%g\n   %s\n',mfilename,fun,s,length(subj),f,length(allparms.function),i,ncalls,datestr(now));
            % mmil_logstr(parms,'%s: function failed: %s\n   Subject %g/%g, Function %g/%g, Call %g/%g\n   %s\n',mfilename,fun,s,length(subj),f,length(allparms.function),i,ncalls,datestr(now));
            rethrow(lasterror);
        end
      elseif script_flag
        % write script to rootoutdir/scripts
      end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      % whether to return output
      keep_flag = (nargout > 0) && (i==ncalls) && (f==length(allparms.function));
      if keep_flag && ischar(otype)
        eval(sprintf('data = %s;',otype));
      end
      if (i==ncalls) && (f<length(allparms.function)) && ...
        ((ischar(otype) && strcmp(otype,allparms.function(f+1).itype)) || ...
        (ischar(allparms.function(f+1).itype) && exist(allparms.function(f+1).itype,'var'))) && ...
        ~(isfield(allparms.function(f+1).parms,'load_flag') && allparms.function(f+1).parms.load_flag)
        keep_flag = 1;
        % hold the output data in memory for the next function
      elseif ischar(otype) && ~save_flag
        % clear the output structure from memory
        eval(sprintf('clear %s',otype));
      end      
      % save_flag
        % 0 - no
        % 1 - save data structure
        % 2 - save figure: save_figure()
        % 3 - hold data structure in memory
      if     save_flag == 0 && ~exist('filename','var')
%         if     ischar(otype) && exist(otype,'var') && ...
%                (i==ncalls) && (f < length(allparms.function)) && ...
%                strcmp(otype,allparms.function(f+1).itype)
               % hold the output data in memory for the next function
%         elseif ischar(otype) && exist(otype,'var') && ...
%                (i==ncalls) && (f==length(allparms.function))
%                % return the output data
%                eval(sprintf('data = %s;',otype));
%         elseif ischar(otype)
%           % clear the output structure from memory
%           eval(sprintf('clear %s',otype));
%         end
      elseif save_flag == 1 || exist('filename','var')
        % save output data and add info to the session structure
        % note: filename contains all files saved during this fun call
        if ~isfield(parms,'filename') || isempty(parms.filename)
          % make empty cell
          parms.filename = {};
        elseif ~iscell(parms.filename)
          % make filename a cell array of strings
          parms.filename = {parms.filename};
        end
        if exist('filename','var') && ~isempty(filename)
          % use independent filename variable returned by this fun
          if ~iscell(filename)  , filename = {filename};  end
          if iscell(filename{1}), filename = filename{1}; end
          if allparms.function(f).spec_flags.oflag == 4
              parms.filename = filename;
          else
              parms.filename = {parms.filename{:} filename{:}}; 
          end
          clear filename;
        end
        if ischar(otype) && run_flag
          if eval(sprintf('isempty(%s)',otype))
            % nothing was returned => do not save or add to session struct
            continue;
          end
          % look for filename field in output structure
          if eval(sprintf('issubfield(%s,''parms.filename'')',otype))
             eval(sprintf('filename = %s.parms.filename;',otype));
             if eval(sprintf('length(fieldnames(%s))==1',otype))
               % the only field in otype is parms; the data is not present
               % so do not save the structure and do not keep it's filename
               parms.filename = {};
               if ~keep_flag, eval(sprintf('clear %s',otype)); end
               otype          = [];
             end
             if ~iscell(filename), filename = {filename}; end
             parms.filename = {parms.filename{:} filename{:}}; 
             try eval(sprintf('%s.parms.filename = parms.filename;',otype)); end
          end
        end
        if ~isfield(parms,'function_id') || isempty(parms.function_id)
          if length(session) >= 1
            parms.function_id = max([session.function_id]) + 1;
          else
            parms.function_id = 1;
          end
        end
        parms.filename = unique(parms.filename);
%           if ~isfield(parms,'filename'), parms.filename = ''; end
%           if ~iscell (parms.filename)  , parms.filename = {parms.filename}; end
        session(end+1).function_id    = parms.function_id;
        session(end)  .name           = fun;
        session(end)  .filename       = parms.filename;
        session(end)  .parms          = parms;
        session(end)  .date           = datestr(now);
        if load_flag == 1 || load_flag == 2
          if ~isfield(parms,'input_id') || isempty(parms.input_id)
            parms.input_id = 0;
          end
          session(end)  .parms.input_id = parms.input_id;
          session(end)  .parms.datafile = datafile;
          if ischar(itype)
            try eval(sprintf('try parms.previous = %s.parms; end',itype)); end
          end
        end
        if ischar(otype) && run_flag
          fprintf('%s: saving data: %s\n',mfilename,parms.filename{1});     
          % mmil_logstr(parms,'%s: saving data: %s\n',mfilename,parms.filename{1});     
          if eval(sprintf('~isfield(%s,''parms'')',otype))
            eval(sprintf('%s.parms = parms;',otype));
          elseif eval(sprintf('~issubfield(%s,''parms.filename'')',otype))
            eval(sprintf('%s.parms.filename = parms.filename;',otype));
          end
          eval(sprintf('save(%s.parms.filename{1},''%s'',''-v7.3'');',otype,otype));
          if ~keep_flag, eval(sprintf('clear %s',otype)); end
        end
%           write_session_log(session(end));
        clear filename;         
      elseif save_flag == 2 && ~cluster_flag
        % save figures
        save_figure(args{:}); 
        if ischar(otype) && ~keep_flag, eval(sprintf('clear %s',otype)); end
      elseif save_flag == 3
        % do nothing
      end

      % clear the input structure from memory
      if ischar(itype) && ~keep_flag && i==ncalls, eval(sprintf('clear %s',itype)); end
      
      % save new session info and update for external sessions
      session = close_session(parms,inparms,session,orig_session,sessionfile,allparms);
      orig_session = session;
      clear filename
    end
    % end loop over ncalls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % submit this function's jobs to the cluster
    if cluster_flag
      write_scriptlist(outpath,funscripts);
      % which cluster
      if isfield(parms,'cluster') && ischar(parms.cluster)
        mmilcluster = parms.cluster;
      else
        switch mod(f,3)
          case 0
            mmilcluster = 'mmilcluster.ucsd.edu';
          case 1
            mmilcluster = 'mmilcluster2.ucsd.edu';
          case 2
            mmilcluster = 'mmilcluster3.ucsd.edu';
        end
      end
      if isfield(parms,'clusterscript') && ischar(parms.clusterscript)
          clusterscript = parms.clusterscript;
      else
          clusterscript = 'qmatjobs';
      end
%       fprintf('submitting job to cluster ''%s'':\n%s\n',mmilcluster,batchdir);
      fprintf('submitting jobs to cluster ''%s'' using %s\n',mmilcluster,clusterscript);
      % mmil_logstr(parms,'submitting jobs to cluster ''%s'' using %s\n',mmilcluster,clusterscript);
      % job statement
      jobstring = sprintf('%s %s',clusterscript,batchdir);
      eval(sprintf('!ssh -X %s "%s"',mmilcluster,jobstring));
      clear jobstring batchdir mmilcluster funscripts outpath
    end
  end
  % end loop over functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  % save subject-specific session
    if ~isempty(session)
%       session = close_session(parms,inparms,session,orig_session,sessionfile,allparms);
      clear session orig_session
    end
    % mmil_logstr(parms,'Total processing time: %g minutes',toc(subj_tstart)/60);
    % mmil_logstr(parms,'Finished at %s\n',datestr(now));
end
% end loop over subjects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEANUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Total processing time: %g minutes\n',toc(tstart)/60);
% nargout

% toc
end 
%% subfunctions
function session = close_session(parms,inparms,session,orig_session,sessionfile,allparms)
if isempty(session), return; end
  % write the session log
  if inparms.verbose
    fprintf('saving parameter logs.\n');
  end
  
  % query sessionfile
  if isempty(sessionfile)
    if ~exist(sessionfile,'file'), sessionfile = sprintf('%s/session.mat',allparms.global.rootoutdir); end
    if ~exist(sessionfile,'file'), sessionfile = sprintf('%s/session.mat',allparms.global.rootindir);  end
    if ~exist(sessionfile,'file'), sessionfile = sprintf('%s/session.mat',pwd); end  
    if ~exist(sessionfile,'file') % no session file exists
      sessionfile = ''; 
    elseif isempty(orig_session)  % session created since start
      orig_session.function_id = [];  
    end
  end
  if ~isempty(sessionfile)
    % check for updates
    % TODO: determine whether updates have occurred by looking at
    % session.date and comparing it to the present time:
    %  if datenum(curr_session(end).date) > datenum(orig_session(end).date)
    %     add additions to curr_session, not orig_session
    %  end
    curr_session = getfield(load(sessionfile),'session');
    if ~isequal([curr_session.function_id],[orig_session.function_id])
      if ~(isfield(curr_session,'date') && isfield(orig_session,'date')) || ...
          (datenum(curr_session(end).date) > datenum(orig_session(end).date))
        % session file has been updated by another session since the start of this one
  %       fprintf('Updating session for subject-specific info added by other processes since this session began...\n');
        try
          % try to add updates to this session
          funids = unique([curr_session.function_id]);
          for k  = 1:length(funids)
            n    = sum([curr_session.function_id] == funids(k)) - sum([orig_session.function_id] == funids(k));
            if n == 0, continue; end
            ix   = find([curr_session.function_id] == funids(k));
            session(end+1:end+n) = curr_session(ix(end-n+1:end));
          end
        catch
          fprintf('\nWARNING: failed to update session.  Some info may have been lost!\n');
          % mmil_logstr(parms,'\nWARNING: failed to update session.  Some info may have been lost!\n');
        end
      end
    end
  end
  if isempty(sessionfile)
    if isfield(parms,'rootoutdir')
      sessionfile = [parms.rootoutdir '/session.mat'];
    else
      sessionfile = [pwd '/session.mat'];
    end
  end
  try
    % write log listing files created and function ids
    write_session_log(session,sessionfile,parms.logfile,parms.logfid); 
  catch
    fprintf('%s: failed to write file log\n',mfilename);
    % mmil_logstr(parms,'%s: failed to write file log\n',mfilename);
  end
%   try
%     % write human-readable list of parameters
%     write_session_parms(session,sessionfile,parms.logfile,parms.logfid);
%   catch
%     fprintf('%s: Failed to write parameter log.\n',mfilename);
%     % mmil_logstr(parms,'%s: Failed to write parameter log.\n',mfilename);
%   end
  try
    % write session MAT file
    save(sessionfile,'session');
    % mmil_logstr(parms,'Saving session MAT file: %s',sessionfile);
  catch
    fprintf('%s: Failed to save session MAT file.\n',mfilename);
    % mmil_logstr(parms,'%s: Failed to save session MAT file.\n',mfilename);
  end
end
function parms = checkparms(parms)
  % check function-specific loop parameter
  % note: loop parameter is the fieldname of a cell array that varies per iteration
  % note: loop parameter may be a cell array of strings
  if (~isfield(parms,'loop') || isempty(parms.loop)) && (isfield(parms,'loop_param') && ~isempty(parms.loop_param))
    parms.loop = parms.loop_param;
    parms      = rmfield(parms,'loop_param');
  elseif ~isfield(parms,'loop') || isempty(parms.loop)
    parms.loop = [];
  end
  % force the loop parameter to be a cell array
  if ~isempty(parms.loop) 
    if ~iscell(parms.loop), parms.loop = {parms.loop}; end
    parms.loop = parms.loop(isfield(parms,parms.loop));
    tmp  = cellfun(@(x) length(parms.(x)),parms.loop);
    parms.loop = parms.loop(tmp>0);
    loop = parms.loop;
    for i = 1:length(loop)
      if ~isfield(parms,loop{i})
        parms.loop = setdiff(parms.loop,loop{i});
        continue;
      end
      if ~iscell(parms.(loop{i}))
        parms.(loop{i}) = {parms.(loop{i})};
      end
    end
  end
end
function write_session_parms(session,sessionfile,logfile,logfid)
[fpath,fname] = fileparts(sessionfile);
outfile = sprintf('%s/parms.log',fpath);  
fid     = fopen(outfile,'wt');
fprintf(fid,'Parameter log: %s\n',date);
fprintf(fid,'Session file: %s\n\n',sessionfile);
fprintf(fid,'---------------------');
for n = 1:length(session);
    fun = session(n).name;
    id = session(n).function_id;
    fprintf(fid,'\n---------------------\nfunction:%s\n',fun);
    if ~isnan(id)
        fprintf(fid,'User assigned ID %g to this function call.\n',id);
    end
    fields = fieldnames(session(n).parms);
    for f=1:length(fields)
        thisparm = fields{f};
        thisval = session(n).parms.(thisparm);
        %str=''; %this will be for the string version of this val
        if ischar(thisval)
            str = thisval;
        elseif isnumeric(thisval)
            str = num2str(thisval);  
        elseif iscellstr(thisval)    %this may not be the correct function
            c  = cell2mat(thisval);
            str = mat2str(c);
        %elseif %continue checking for other datatypes
        end
        fprintf(fid, '%-30s: %s\n',thisparm,str);
    end %end loop over parameters
end %end loop over elements of the "session" structure array
%write the footer
fclose(fid);
% fprintf('%s: parameter list compiled on %s\n',mfilename,date);
parms.logfile = logfile; 
parms.logfid  = logfid;
% mmil_logstr(parms,'Saving parameter log: %s',outfile);
end
function write_session_log(session,sessionfile,logfile,logfid)
  % saves "outfiles.log" => { function_id   function_name   [input_id] }
  % note: append if it already exists
  [fpath,fname] = fileparts(sessionfile);
  outfile = sprintf('%s/outfiles.log',fpath);
  fid = fopen(outfile,'wt');
  fprintf(fid,'%-15s %-20s [%s] %-20s\n','function id','function name','input id','outfile');
  for i = 1:length(session)
    nfile = 0;
    if isnumeric(session(i).function_id)    && ~isempty(session(i).function_id)
      function_id = session(i).function_id;
    else
      function_id = 0;
    end
    if issubfield(session(i),'parms.input_id') && isnumeric(session(i).parms.input_id) && ~isempty(session(i).parms.input_id)
      input_id = session(i).parms.input_id;
    else
      input_id = 0;
    end  
    if isfield(session(i),'filename') && ~isempty(session(i).filename)
      if iscell(session(i).filename)
        filename = session(i).filename{1};
        nfile   = length(session(i).filename);
      else
        filename = session(i).filename;
      end
    else
      filename = '';
    end
    fprintf(fid,'%-15g %-20s [%g] %-20s\n',function_id,session(i).name,input_id,filename);
    clear function_id input_id
    if nfile > 1
      for f = 1:nfile-1
        fprintf(fid,'%45s %s\n',' ',session(i).filename{f+1});
      end
    end
  end
  fprintf(fid,'\n\n');
  fclose(fid);
  parms.logfile = logfile; 
  parms.logfid  = logfid;
  % mmil_logstr(parms,'Saving file log: %s',outfile);
end
function write_scriptlist(outpath,funscripts)
  outfile = fullfile(outpath,'scriptlist.txt');
  % fprintf('%s: writing script list for cluster computing: %s\n',mfilename,outfile);
  fid = fopen(outfile,'wt');
  for i = 1:length(funscripts), fprintf(fid,'%s\n',funscripts{i}); end
  fclose(fid);
end
function save_figure(varargin)
  parms = mmil_args2parms(varargin,...
              {'format','eps',[],...
               'rootoutdir',pwd,[],...
               'savepath',[],[],...
                           'filename',[],[],...
               'save',1,{1,0},...
               'close',0,{1,0},...
                           'prefix','proc',[],...
                           'overwrite',0,{1,0}...
              },false);

  if isempty(parms.filename)
      parms.filename{1} = sprintf('%s/images/%s_%s',parms.rootoutdir,parms.prefix,datestr(now,30));
  else
      [pathstr,fname] = fileparts(parms.filename{1});
      parms.filename{1} = sprintf('%s/images/%s',parms.rootoutdir,fname);
  end
  [pathstr,fname,ext] = fileparts(parms.filename{1});
  if ~exist(pathstr,'dir'),
      fprintf('making output directory: %s\n',pathstr);
      unix(['mkdir -p ' pathstr]);
  end
  if strcmp(parms.format,'jpg'), parms.printcmd = {'-djpeg'}; end
  if strcmp(parms.format,'eps'), parms.printcmd = {'-depsc','-tiff','-r150'}; end
  if strcmp(parms.format,'tif'), parms.printcmd = {'-dtiff'}; end
  fprintf('%s: saving figure: %s\n',mfilename,parms.filename{1});
  if parms.save, print(gcf,parms.printcmd{:},parms.filename{1}); end
  if parms.close, close; end
end

% % function writesummary(session)
% % fprintf('%s: writing session summary and saving session parameters.\n',mfilename);
% % fid = fopen('session_parameters.log','wt');
% % fprintf(fid,'parameter record for successful functions (cells and structs are not shown)\n');
% % for i = 1:length(session)         % subjects
% %   subj = session(i);
% %   funs = fieldnames(subj);
% %   for j = 1:length(funs)          % functions
% %     fun = funs{i};
% %     funparm = subj.(fun);
% %     for k = 1:length(funparm)     % function iterations
% %       fprintf(fid,'\n---------------------------------------------------------------------------\n');
% %       fprintf(fid,'subject %d, iteration %d of function ''%s''\n',i,k,fun);
% %       fprintf(fid,'---------------------------------------------------------------------------\n\n');
% %       opt   = funparm(k);
% %       parms = fieldnames(opt);
% %       for p = 1:length(parms)
% %         parm = parms{p};
% % 				pclass = class(opt.(parm));
% %         switch pclass
% %           case 'logical'
% %             value = opt.(parm);
% %             cchar = '%d';
% %           case 'char'
% %             value = opt.(parm);
% %             cchar = '%s';
% %           case 'struct'
% %             value = pclass;
% %             cchar = '%s';
% % 					case 'cell'
% % 						if length(opt.(parm))==1 && ischar(opt.(parm){1})
% % 							value = opt.(parm){1};
% % 						else
% % 							value = pclass;
% % 						end
% % 						cchar = '%s';
% %           otherwise
% %             if isnumeric(opt.(parm))
% %               value = num2str(opt.(parm));
% %               cchar = '[%s]';
% %             else
% %               continue;
% %             end
% %         end
% %         cmdstring = sprintf('%s = %s;\\n',parm,cchar);
% %         eval(sprintf('fprintf(fid,''%s'',value);',cmdstring));
% %       end
% %     end
% %   end
% % end    
% % fclose(fid);
