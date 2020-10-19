function varargout = ts_timefreq_rejection(datafile,varargin)
%% script for TF wave-based visual trial rejection

parms = mmil_args2parms(varargin,...
						{'events',[],[],...
						 'conditions',[],[],...
             'badchanfile',[],[],...
             'badchans',[],[],...
             'rejectfile',[],[],...
             'overwrite',0,{0,1},...
						 'toilim',[],[],...
             'toi',[],[],...
             'foi',[],[],...
             'prefix','proc',[],...
             'fpath',[],[],...
             'fname',[],[],...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...                
						},false);


if iscell(datafile) && ischar(datafile{1})
  datafile = datafile{1};
end
if ischar(datafile) && exist(datafile,'file')
  load(datafile) % epoch_data
elseif isstruct(datafile)
  epoch_data = datafile; 
  epoch_data = ts_checkdata_header(epoch_data,'events',parms.events);
  if issubfield(datafile,'parms.filename')
    datafile = datafile.parms.filename;
  end
    if iscell(datafile), datafile = datafile{1}; end
else
  mmil_logstr(parms,'%s: TF waves not found!\n',mfilename);
%   fprintf('%s: TF waves not found!\n',mfilename);
  return;
end

if exist(parms.badchanfile,'file')
  epoch_data = ts_data_selection(epoch_data,'badchanfile',parms.badchanfile,'badchans',parms.badchans,'removebadchans',1);
end

[data badchans badtrials] = ts_visual_reject(epoch_data,varargin{:});

% store reject info in matlab structure
reject_data = [];
reject_data.badchans  = badchans;
reject_data.badtrials = [];
for i = 1:length(epoch_data.epochs)
  reject_data.badtrials{i}  = badtrials{i};
  reject_data.event_code(i) = epoch_data.epochs(i).event_code;
end

if isempty(parms.fpath)
  fpath = fileparts(datafile);
else
  fpath = parms.fpath;
end
if isempty(parms.fname)
  fname = sprintf('%s_tf_reject_data.mat',parms.prefix);
else
  fname = parms.fname;
end

outfile = fullfile(fpath,fname);
mmil_logstr(parms,'%s: saving TF reject data: %s\n',mfilename,outfile);
% fprintf('%s: saving TF reject data: %s\n',mfilename,outfile);
save(outfile,'reject_data')

if nargout == 1
  varargout{1} = outfile;
else
  varargout{1} = data;
  varargout{2} = badchans;
  varargout{3} = badtrials;
end

% write reject info to text file
[datatype datafield dataparam] = ts_object_info(data);
outfile = strrep(outfile,'.mat','.txt');
fid = fopen(outfile,'w+');
badchans = find([data.sensor_info.badchan]);  
fprintf(fid,'Rejected Channels\n\n');
if ~isempty(badchans)
  for j = 1:length(badchans)
    fprintf(fid,'%s\n',data.sensor_info(badchans(j)).label);
  end
else
  fprintf(fid,'No bad channels.\n');
end
fprintf(fid,'\n');
fprintf(fid,'Trial Info\n\n');
fprintf(fid,'         \t Good\t    Rejected     \t Bad\n');
fprintf(fid,'Condition\tTrials\tEEG File\tManual\tTrials\n\n');
for j = 1:length([data.epochs.event_code])
  if ~ismember(data.(datafield)(j).event_code,[reject_data.event_code])
    continue;
  end
  if j > length(reject_data.badtrials), break; end
 fprintf(fid,'%-9s\t%-6s\t%-8s\t%-6s\t%s\n',num2str(j),num2str(data.(datafield)(j).num_trials),...
                                                   num2str(data.(datafield)(j).num_rejects.eeg),...
                                                   num2str(data.(datafield)(j).num_rejects.manual),...
                                                   num2str(reject_data.badtrials{j}));
end
fclose(fid);  
