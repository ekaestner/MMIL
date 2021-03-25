function guiplot

% Will plot the average over trials for one sensor

hMainFig = figure('Visible','on','Position',[400,500,600,750]);%,'Color','w');
% [distance from left, distance from bottom, width, height]

% Construct components
hplot     = axes('Units','normalized','Parent',hMainFig,...
    'Position',[.10 .50 .80 .37]);
hxlim1    = uicontrol('Style','edit','String','-.5','units','normalized',...
    'Position',[.75 .88 .07 .03],'Callback',{@axislimits_Callback});
hxlim2    = uicontrol('Style','edit','String','1.5','units','normalized',...
    'Position',[.83 .88 .07 .03],'Callback',{@axislimits_Callback});
hylim1    = uicontrol('Style','edit','String','-100','units','normalized',...
    'Position',[.91 .80 .07 .03],'Callback',{@axislimits_Callback});
hylim2    = uicontrol('Style','edit','String','100','units','normalized',...
    'Position',[.91 .84 .07 .03],'Callback',{@axislimits_Callback});
hxleft    = uicontrol('Style','pushbutton','String','<<','units','normalized',...
    'Position',[.75 .91 .07 .03],'Callback',{@axis_xleft_Callback});
hxright   = uicontrol('Style','pushbutton','String','>>','units','normalized',...
    'Position',[.83 .91 .07 .03],'Callback',{@axis_xright_Callback});
hxlimset  = uicontrol('Style','pushbutton','String','XLim','units','normalized',...
    'Position',[.68 .88 .06 .03],'Callback',{@set_XLim_Callback});
hxreset   = uicontrol('Style','pushbutton','String','reset','units','normalized',...
    'Position',[.60 .88 .07 .03],'Callback',{@reset_XLim_Callback});
hylimset  = uicontrol('Style','pushbutton','String','YLim','units','normalized',...
    'Position',[.91 .76 .07 .03],'Callback',{@set_YLim_Callback});
hyreset   = uicontrol('Style','pushbutton','String','reset','units','normalized',...
    'Position',[.91 .72 .07 .03],'Callback',{@reset_YLim_Callback});
hopenfile = uicontrol('Style','pushbutton','String','Load file','units','normalized',...
    'Position',[.05 .09 .15 .03],'Callback',{@open_file_Callback});
hclose    = uicontrol('Style','pushbutton','String','Close','units','normalized',...
    'Position',[.05 .03 .15 .03],'Callback',{@close_Callback});  
hrefchan  = uicontrol('Style','edit','String','STI101','units','normalized',...
    'Position',[.10 .07 .10 .02]);    
hreflbl   = uicontrol('Style','text','String','trig','units','normalized',...
    'Position',[.05 .07 .05 .02]);    
hfiletext = uicontrol('Style','text','String','No file has been loaded.','units','normalized',...
    'Position',[.05 .00 .90 .02]);
hica      = uicontrol('Style','pushbutton','String','Manual ICA','units','normalized',...
    'Position',[.25 .02 .15 .03],'Callback',{@ICAanalysis_Callback});
hpreproc  = uicontrol('Style','pushbutton','String','Preprocess','units','normalized',...
    'Position',[.25 .05 .15 .03],'Callback',{@preprocess_Callback});
hprocall  = uicontrol('Style','checkbox','String','All sensors','Value',0,'units','normalized',...
    'Position',[.25 .08 .15 .02]);
hprocFT   = uicontrol('Style','checkbox','String','FT','Value',0,'units','normalized',...
    'Position',[.25 .17 .15 .02]);
hprocdet  = uicontrol('Style','checkbox','String','detrend','Value',1,'units','normalized',...
    'Position',[.25 .15 .15 .02]); 
hprocblc  = uicontrol('Style','checkbox','String','blc','Value',1,'units','normalized',...
    'Position',[.25 .13 .07 .02]);
hprocbc1   = uicontrol('Style','edit','String','-inf','units','normalized',...
    'Position',[.32 .13 .05 .02]);
hprocbc2   = uicontrol('Style','edit','String','0','units','normalized',...
    'Position',[.37 .13 .03 .02]);
hprocbp   = uicontrol('Style','checkbox','String','bp','Value',1,'units','normalized',...
    'Position',[.25 .11 .07 .02]);  
hprocf1   = uicontrol('Style','edit','String','.2','units','normalized',...
    'Position',[.32 .11 .04 .02]);
hprocf2   = uicontrol('Style','edit','String','30','units','normalized',...
    'Position',[.36 .11 .04 .02]);
hstats    = uicontrol('Style','pushbutton','String','MC Statistics','units','normalized',...
    'Position',[.45 .05 .15 .03],'Callback',{@stats_Callback});
hstatall  = uicontrol('Style','checkbox','String','All sensors','Value',0,'units','normalized',...
    'Position',[.45 .08 .15 .02]);
hstatn    = uicontrol('Style','edit','String','200','units','normalized',...
    'Position',[.45 .11 .07 .02]);
hstatp1    = uicontrol('Style','edit','String','.05','units','normalized',...
    'Position',[.45 .13 .07 .02]);
hstatp2    = uicontrol('Style','edit','String','.05','units','normalized',...
    'Position',[.53 .13 .07 .02]);  
hstatp1lbl = uicontrol('Style','text','String','samp','units','normalized',...
    'Position',[.45 .15 .07 .02]);
hstatp2lbl = uicontrol('Style','text','String','clust','units','normalized',...
    'Position',[.53 .15 .07 .02]);    
hstatc    = uicontrol('Style','edit','String','b','units','normalized',...
    'Position',[.52 .11 .04 .02]);
hstatt    = uicontrol('Style','edit','String','.1','units','normalized',...
    'Position',[.56 .11 .04 .02]);  
htfstat   = uicontrol('Style','pushbutton','String','TF Stats','units','normalized',...
    'Position',[.65 .02 .15 .03],'Callback',{@TFstats_Callback});
hpower    = uicontrol('Style','pushbutton','String','TF Analysis','units','normalized',...
    'Position',[.65 .05 .15 .03],'Callback',{@TFanalysis_Callback});
hband     = uicontrol('Style','pushbutton','String','Band Plot','units','normalized',...
    'Position',[.83 .05 .15 .03],'Callback',{@BandPlot_Callback});  
hbandlo   = uicontrol('Style','edit','String','0','units','normalized',...
    'Position',[.83 .03 .07 .02]);
hbandhi   = uicontrol('Style','edit','String','0','units','normalized',...
    'Position',[.91 .03 .07 .02]);  
hfoibg    = uibuttongroup('Parent',hMainFig,'Title','foi',...
    'Position',[.65 .08 .12 .08]);
hfoimeg   = uicontrol(hfoibg,'Style','radiobutton','String','meeg','units','normalized',...
    'Position',[.10 .00 .90 .50]);
hfoieeg   = uicontrol(hfoibg,'Style','radiobutton','String','ieeg','units','normalized',...
    'Position',[.10 .50 .90 .50]);
hprevious = uicontrol('Style','pushbutton','String','Previous Sensor','units','normalized',...
    'Position',[.10 .43 .20 .03],'Callback',{@previous_sensor_Callback});
hnext     = uicontrol('Style','pushbutton','String','Next Sensor','units','normalized',...
    'Position',[.70 .43 .20 .03],'Callback',{@next_sensor_Callback});
hsensor   = uicontrol('Style','popupmenu','String','Select a Sensor','units','normalized',...
    'Position',[.70 .38 .20 .05],'Callback',{@select_sensor_Callback});
label1    = uicontrol('Style','text','String','Conditions','units','normalized',...
    'Position',[.10 .35 .12 .02]);
hconds    = uicontrol('Style','listbox','String','Select a Condition','units','normalized','Max',10,'Min',0,...
    'Position',[.10 .25 .25 .10],'Callback',{@select_condition_Callback});  
label2    = uicontrol('Style','text','String','Layouts','units','normalized',...
    'Position',[.39 .35 .08 .02]);
hlay      = uicontrol('Style','listbox','String','Select a Layout','units','normalized',...
    'Position',[.39 .25 .30 .10],'Callback',{@select_layout_Callback});    
happlay   = uicontrol('Style','checkbox','String','Apply','Value',0,'units','normalized',...
    'Position',[.47 .35 .10 .02],'Callback',{@select_layout_Callback});  
hchtoggle = uicontrol('Style','pushbutton','String','toggle','Value',0,'units','normalized',...
    'Position',[.57 .35 .10 .02],'Callback',{@toggle_channel_Callback});  
hplottype = uibuttongroup('Parent',hMainFig,'Title','Plot Type',...
    'Position',[.71 .25 .17 .11]);
hploterp  = uicontrol(hplottype,'Style','radiobutton','String','waveform','units','normalized',...
    'Position',[.10 .66 .90 .33]);
hplotpow  = uicontrol(hplottype,'Style','radiobutton','String','power','units','normalized',...
    'Position',[.10 .33 .90 .33]);
hplotplv  = uicontrol(hplottype,'Style','radiobutton','String','plv','units','normalized',...
    'Position',[.10 .00 .90 .33]);  
hgptrials = uibuttongroup('Parent',hMainFig,...
    'Position',[.88 .25 .11 .10]);
hplotavg  = uicontrol(hgptrials,'Style','radiobutton','String','avg','units','normalized',...
    'Position',[.10 .50 .90 .45]);
hplottrl  = uicontrol(hgptrials,'Style','radiobutton','String','trials','units','normalized',...
    'Position',[.10 .00 .90 .45]);
hreject   = uicontrol('Style','pushbutton','String','reject','units','normalized',...
    'Position',[.88 .22 .11 .03],'Callback',{@reject_Callback});  

% Assign name to appear in the window title
set(hMainFig,'Name','MEG, EEG, and iEEG visualization');

% Load data for plot
function open_file_Callback(hObject,eventdata,handles)
[filename,pathname] = uigetfile({'*.mat;*.fif;*.eeg;*.avg;*.cnt;*.vhdr;*.set'},'Pick a file.','MultiSelect','on');
handles = guidata(gcbo);
if isequal(filename,0) || isequal(pathname,0)
  return;
end
tic
if ~iscell(filename)
  filename = {filename};
end
datatype  = 'epoch_data';    
datafield = 'epochs';   
zparam    = 'data';
for k = 1:length(filename)
  fname    = filename{k};
  datafile = [pathname filename{k}];
  [jnk1 jnk2 ftype] = fileparts(fname);
  switch ftype
    case {'.eeg','.vhdr','.set'}
      alldata(k) = ts_iEEG_eeg2epoch(datafile);
    case '.cnt'
      alldata(k) = ts_loadcnt(datafile);
    case '.avg'
      tmp         = ts_matrix2avg(loadavg(datafile));
      tmp.epochs  = rmfield(tmp.averages,'std_dev');
      alldata(k)  = rmfield(tmp,'averages');
      clear tmp
    case '.fif'
      hdr = ts_read_fif_header(datafile,0);
      try
        % assume data should be epoched
        ts_process_fif_data(datafile,...
          'saveepochs_flag',1,...
          'trigchan',get(hrefchan,'String'),...
          'valid_event_codes',[],...
          'stim_delay',0,...
          'prestim_dur',100,...
          'poststim_dur',400,...
          'reject_mag',10000,...
          'reject_grad',6000,...
          'reject_eeg',0,...
          'reject_eog',200,...
          'detrend_flag',0,...
          'baseline_flag',0,...
          'forceflag',0,...
          'saveperevent',0);
        infile      = sprintf('%s/matfiles/proc_epoch_data_1.mat',pwd);
        alldata(k)  = getfield(load(infile),'epoch_data');
        clear infile
      catch
        fprintf('failed to read epochs from file\nreading continuous data...\n');
        try
          % try loading as continuous data
          alldata(k)  = ts_loadfif(datafile,'epochs',1);
          if k == length(filename)
            try
              matfile = 'raw_cont_data.mat';
              cont_data = alldata;
              fprintf('saving cont_data to mat file: %s\n',matfile);
              save(matfile,'cont_data');
              clear cont_data
            end
%             try
%               fprintf('saving cont_data to fif file\n');
%               ts_avg2fif(ts_trials2avg(alldata),datafile,'raw_cont_data',0,1);
%             end
          end
        catch
          % report sensors in the file that failed to load
          fprintf('failed to load the following channels:\n');
          hdr.sensors.label'
          fprintf('try changing the trigger channel to one of the following:\n');
          for n = 1:length(hdr.sensors.label)
            if findstr('STI',hdr.sensors.label{n})
              fprintf('%s\n',hdr.sensors.label{n});
            end
          end
        end
      end
    case '.mat'
      s = load(datafile);
      if      isfield(s,'cont_data')
        datatype  = 'cont_data';    
      elseif  isfield(s,'epoch_data')
        datatype  = 'epoch_data';    
        datafield = 'epochs';   
        zparam    = 'data';           
      elseif  isfield(s,'avg_data')
        datatype  = 'avg_data';
        datafield = 'averages';
        zparam    = 'data';     
      elseif  isfield(s,'timefreq_data')
        datatype  = 'timefreq_data'; 
        datafield = 'timefreq';
        zparam    = 'power';
        set(hplotpow,'Value',1);
      elseif  isfield(s,'stat_data')
%         if isfield(handles,'current_datatype'),  origtype  = handles.current_datatype;  end
%         if isfield(handles,'current_datafield'), datafield = handles.current_datafield; end
%         if isfield(handles,'current_zparam'),    zparam    = handles.current_zparam;    end
        if ndims(s.stat_data.stats(1).mask) > 2
          datatype = 'tfstat_data';
        else
          datatype = 'stat_data';
        end
      else
        disp('mat file data type not recognized.');
        return;
      end
      alldata(k) = s.(datatype);
      clear s;
      % identify cont_data as epoch_data
      if strcmp(datatype,'cont_data')
        datatype = 'epoch_data';
      end
  end
end
if findstr(datatype,'stat')
  % only keep the between condition stats
  for k = 1:length(alldata)
    if k == 1
      temp = alldata(1);
      temp.stats = temp.stats(1);
    elseif isequal({alldata(1).sensor_info.label},{alldata(k).sensor_info.label})
      temp.stats(end+1) = alldata(k).stats(1);
    else
      continue;
    end
    try
      temp.stats(end).event_code = [alldata(k).stats(2).event_code alldata(k).stats(3).event_code]; 
    catch
      temp.stats(end).event_code = handles.current_condition;
    end
  end
  if isfield(handles,datatype) && isequal({handles.(datatype)(1).sensor_info.label},{temp.sensor_info.label})
    handles.(datatype).stats(end+1:end+length(temp.stats)) = temp.stats;
  else
    handles.(datatype) = temp;
  end
  clear temp
  guidata(gcbo,handles);
  axislimits_Callback(hObject,eventdata,handles);
  return;
end

if k > 1
  fprintf('%s: combining data from multiple files\n',mfilename);
  alldata = ts_combine_data(alldata);
end
[nullvector{1:length(alldata.sensor_info)}] = deal([]);
for k = 1:length(alldata.(datafield))
  handles.rejects(k).badtrials = nullvector;
  handles.rejects(k).powtrials = nullvector;
  handles.rejects(k).badchans  = find([alldata.sensor_info.badchan]);
end
fprintf('%s: %g file(s) loaded successfully.\n',mfilename,k);
toc
set(hfiletext,'String',[pathname fname]);
% sensors
set(hsensor,'String',{alldata.sensor_info.label});
set(hsensor,'Value',1);
% event codes
evcodes = [alldata.(datafield).event_code];
% update datatype if timefreq_data has trials
if strcmp(datafield,'timefreq')
  if isfield(alldata.(datafield),'cmplx')
    zparam = 'cmplx';
  end
  if ndims(alldata.(datafield)(1).(zparam)) > 3
    datatype = 'timefreq_trials';
  end
end
% determine whether the data is iEEG or MEEG
if (ismember('eeg',unique({alldata.sensor_info.typestring})) ...
             && isempty(intersect({'grad1','grad2','mag','meg'},{alldata.sensor_info.typestring})))
	handles.iEEG  = 1;
  handles.plots.ylabel = 'amplitude (uV)';
else
  handles.iEEG = 0;
  handles.plots.ylabel = 'amplitude (T)';
  label = {alldata.sensor_info.label};
  ix = find(cellfun('isempty',regexp(label,'MEG\s+\d+')));
  for j = 1:length(ix)
    alldata.sensor_info(ix(j)).label = strrep(alldata.sensor_info(ix(j)).label,'MEG','MEG ');
  end
end
% save current state info
handles.evcodes           = evcodes;
handles.(datatype)        = alldata;
clear alldata;
handles.plots.maxmin      = 1;
handles.plots.graphcolor  = 'brgkyrgbkyrgbkyrgbkyrgbkyrgkyrgbkyrgbkyrgbkyrgbkyrgkyrgbkygkyrgbky';
handles.plots.caxis       = [-10 10];
handles.current_sensor    = 1;
handles.current_condition = 1;
handles.preproc.hprocbc1 = str2num(get(hprocbc1,'String'));
handles.preproc.hprocbc2 = str2num(get(hprocbc2,'String'));
handles.tempfiles = {};
handles = update_datatype(handles,datatype);
handles = set_layouts(handles);
tmp = {};
if ~isempty(handles.cond_labels)
  for k = 1:length(evcodes)
    tmp = {tmp{:} sprintf('%g - %s (%s)',evcodes(k),handles.cond_labels{[handles.cond_evcode]==evcodes(k)},handles.plots.graphcolor(k))}; 
  end
else
  for k = 1:length(evcodes)
    tmp = {tmp{:} sprintf('%g (%s)',evcodes(k),handles.plots.graphcolor(k))}; 
  end
end
set(hconds,'String',tmp);
set(hconds,'Value',1);
guidata(gcbo,handles);
axislimits_Callback(hObject,eventdata,handles);
end

function previous_sensor_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
ch      = handles.current_sensor - 1;
if ch < 1
  ch = handles.(handles.current_datatype).num_sensors;
end
handles.current_sensor = ch;
guidata(gcbo,handles);
update_figure(handles);
end

function next_sensor_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
ch      = handles.current_sensor + 1;
if ch > handles.(handles.current_datatype).num_sensors
  ch = 1;
end
handles.current_sensor = ch;
guidata(gcbo,handles);
update_figure(handles);
end

function select_sensor_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
ch      = get(hObject,'Value');
handles.current_sensor = ch;
guidata(gcbo,handles);
update_figure(handles);
end  

function select_condition_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
handles.current_condition = get(hObject,'Value');
guidata(gcbo,handles);
update_figure(handles);
end

function select_layout_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
handles.current_layout = get(hObject,'Value');
guidata(gcbo,handles);
update_figure(handles);
end

function preprocess_Callback(hObject,eventdata,handles)
tic
handles = guidata(gcbo);
cf1 = str2num(get(hprocf1,'String'));
cf2 = str2num(get(hprocf2,'String'));
bc1 = str2num(get(hprocbc1,'String'));
bc2 = str2num(get(hprocbc2,'String'));
datatype = handles.current_datatype;
datafield = handles.current_datafield;
if ~get(hprocall,'Value')
  epochs  = ts_data_selection(handles.(datatype),'channels',handles.current_sensor);
  if get(hprocFT,'Value')    
    if get(hprocdet,'Value'), det='yes'; else det='no'; end
    if get(hprocblc,'Value'), blc='yes'; else blc='no'; end
    if get(hprocbp,'Value'),  lpf='yes'; else lpf='no'; end
    epochs = ts_preproc(epochs,'lpfilter',lpf,'lpfreq',cf2,...
      'detrend',det,'blc',blc,'blcwindow',[bc1 bc2],'fieldtrip_flag',1);    
  else
    epochs = ts_preproc(epochs,'bandpass_flag',get(hprocbp,'Value'),'bpfreq',[cf1 cf2],...
      'detrend_flag',get(hprocdet,'Value'),'baseline_flag',get(hprocblc,'Value'),...
      'baseline_start',bc1,'baseline_end',bc2);
  end    
  for k = 1:length(epochs.(datafield))
    handles.(datatype).(datafield)(k).data(handles.current_sensor,:,:) = epochs.(datafield)(k).data;
  end
else
  if get(hprocFT,'Value')    
    if get(hprocdet,'Value'), det='yes'; else det='no'; end
    if get(hprocblc,'Value'), blc='yes'; else blc='no'; end
    if get(hprocbp,'Value'),  lpf='yes'; else lpf='no'; end
    handles.(datatype) = ts_preproc(handles.(datatype),'lpfilter',lpf,'lpfreq',cf2,...
      'detrend',det,'blc',blc,'blcwindow',[bc1 bc2],'fieldtrip_flag',1);      
  else
    handles.(datatype) = ts_preproc(handles.(datatype),'bandpass_flag',get(hprocbp,'Value'),'bpfreq',[cf1 cf2],...
      'detrend_flag',get(hprocdet,'Value'),'baseline_flag',get(hprocblc,'Value'),...
      'baseline_start',bc1,'baseline_end',bc2);
  end
end
toc
guidata(gcbo,handles);
update_figure(handles);
end

function stats_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
if length(handles.current_condition) ~= 2, return; end
if ~strcmp(handles.current_datatype,'epoch_data'), return; end
if ~get(hstatall,'Value')
  chans = handles.current_sensor;
else
  chans = 1:handles.epoch_data.num_sensors;
end
tic
cs = handles.current_condition;
for k = 1:length(chans)
  if handles.epoch_data.sensor_info(chans(k)).badchan
    fprintf('skipping bad channel: %s\n',handles.epoch_data.sensor_info(chans(k)).label);
    continue; 
  end
  if length(chans)>1
    fprintf('running stats on channel %g of %g\n',k,length(chans)); 
  end
  ch = chans(k);
  epochs = ts_data_selection(handles.epoch_data,'channels',ch);
  if ~isstruct(epochs)
    fprintf('skipping empty channel: %s\n',handles.epoch_data.sensor_info(ch).label);
    continue; 
  end
  for c = 1:length(cs)
    if isempty(handles.rejects(cs(c)).badtrials{ch}), continue; end
    trl = setdiff(1:handles.epoch_data.epochs(cs(c)).num_trials,handles.rejects(cs(c)).badtrials{ch});
    epochs.epochs(cs(c)).data       = epochs.epochs(cs(c)).data(:,:,trl);
    epochs.epochs(cs(c)).num_trials = length(trl);
  end 
  stats = ts_statistics(epochs,'numrandomization',str2num(get(hstatn,'String')),'type','between',...
    'events',[handles.epoch_data.epochs(handles.current_condition).event_code],...
    'clusteralpha',str2num(get(hstatp1,'String')),'alpha',str2num(get(hstatp2,'String')));
  stats.stats(1).event_code = [handles.epoch_data.epochs(handles.current_condition).event_code];
  statix = 0;
  if ~isfield(handles,'stat_data')
    n  = handles.epoch_data.num_sensors;
    nt = length(handles.epoch_data.epochs(handles.current_condition(1)).time);
    handles.stat_data.sensor_info            = handles.epoch_data.sensor_info;
%     handles.stat_data.sensor_info            = init_sensor_info(n);
%     [handles.stat_data.stats(1).mask]        = false(n,nt);
%     [handles.stat_data.stats(1).prob]        = zeros(n,nt);
    [handles.stat_data.stats(1).label{1:n}]  = deal('');
    statix = 1;  
  else
    tgt = stats.stats(1).event_code;
    for k = 1:length(handles.stat_data.stats)
      if isequal(tgt,[handles.stat_data.stats(k).event_code])
        statix = k;
        break;
      end
    end    
    if statix == 0
      statix = length(handles.stat_data.stats)+1;
      [handles.stat_data.stats(statix).label{1:handles.epoch_data.num_sensors}]  = deal('');
    end
  end
  [ch sel2] = match_str({handles.stat_data.sensor_info.label},{stats.sensor_info.label});
  try
    handles.stat_data.sensor_info(ch)  = stats.sensor_info;
  catch
    stats.sensor_info = orderfields(stats.sensor_info,handles.stat_data.sensor_info);
    handles.stat_data.sensor_info(ch)  = stats.sensor_info;
  end
  handles.stat_data.stats(statix).mask(ch,:) = stats.stats(1).mask;
  handles.stat_data.stats(statix).prob(ch,:) = stats.stats(1).prob;
  handles.stat_data.stats(statix).label(ch)  = stats.stats(1).label;
  handles.stat_data.stats(statix).event_code = stats.stats(1).event_code;
  handles.stat_data.stats(statix).time       = stats.stats(1).time;
  clear epochs;
end
toc
guidata(gcbo,handles);
update_figure(handles);
end

function TFstats_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
if ~isfield(handles,'timefreq_trials')
  fprintf('no trial power found. aborting stats.\n');
  return;
end
if length(handles.current_condition) ~= 2, return; end
if ~get(hstatall,'Value')
  chans = handles.current_sensor;
else
  chans = 1:handles.epoch_data.num_sensors;
end
tic
cs = handles.current_condition;
for k = 1:length(chans)
  if handles.timefreq_trials.sensor_info(chans(k)).badchan
    fprintf('skipping bad channel: %s\n',handles.timefreq_trials.sensor_info(chans(k)).label);
    continue; 
  end
  if length(chans)>1
    fprintf('running stats on channel %g of %g\n',k,length(chans)); 
  end  
  ch = chans(k);
  tfdata = ts_data_selection(handles.timefreq_trials,'channels',ch);  
  if ~isstruct(tfdata)
    fprintf('skipping empty channel: %s\n',handles.timefreq_trials.sensor_info(ch).label);
    continue; 
  end
  for c = 1:length(cs)
    if isempty(handles.rejects(cs(c)).powtrials{ch}), continue; end
    trl = setdiff(1:handles.timefreq_trials.timefreq(cs(c)).num_trials,handles.rejects(cs(c)).powtrials{ch});
    tfdata.timefreq(cs(c)).power       = tfdata.timefreq(cs(c)).power(:,:,trl);
    tfdata.timefreq(cs(c)).num_trials = length(trl);
  end 
  stats = ts_freqstatistics(tfdata,'numrandomization',str2num(get(hstatn,'String')),'type','between',...
    'events',[handles.timefreq_trials.timefreq(handles.current_condition).event_code],...
    'clusteralpha',str2num(get(hstatp1,'String')),'alpha',str2num(get(hstatp2,'String')));
  if ~isfield(handles,'tfstat_data')
    n  = handles.timefreq_trials.num_sensors;
    nt = length(handles.timefreq_trials.timefreq(handles.current_condition(1)).time);
    nf = length(handles.timefreq_trials.timefreq(handles.current_condition(1)).frequencies);
    handles.tfstat_data.sensor_info            = init_sensor_info(n);
    [handles.tfstat_data.stats(1).mask]        = false(n,nt,nf);
    [handles.tfstat_data.stats(1).prob]        = zeros(n,nt,nf);
    [handles.tfstat_data.stats(1).label{1:n}]  = deal('');
  end
  try
    handles.tfstat_data.sensor_info(ch)  = stats.sensor_info;
  catch
    stats.sensor_info = orderfields(stats.sensor_info,handles.tfstat_data.sensor_info);
    handles.tfstat_data.sensor_info(ch)  = stats.sensor_info;
  end
  handles.tfstat_data.stats(1).mask(ch,:,:) = stats.stats(1).mask;
  handles.tfstat_data.stats(1).prob(ch,:,:) = stats.stats(1).prob;
  handles.tfstat_data.stats(1).label(ch)    = stats.stats(1).label;
  handles.tfstat_data.stats(1).time         = stats.stats(1).time;
  handles.tfstat_data.stats(1).frequencies  = stats.stats(1).frequencies;
  clear tfdata;
end
toc
guidata(gcbo,handles);
update_figure(handles);
end
  
function TFanalysis_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
if length(handles.current_condition) > 1
  fprintf('Select only one condition to run TF analysis.\n');
  return;
else
  c = handles.current_condition;
end
tic
fc = str2num(get(hprocf2,'String'));
if get(hprocbp,'Value'),  lpf='yes'; else lpf='no'; end
% if all(get(hprocf2,'String')~=0) && all(get(hprocf2,'String')~=0)
%   fc = str2num(get(hprocf2,'String'));
%   lpfilt = 'yes';
% else
%   lpfilt = 'no';
% end
if get(hfoimeg,'Value')
  foi = 'meg';
elseif get(hfoieeg,'Value')
  foi = 'ieeg';
else
  foi = 'ieeg';
end
toilim = [str2num(get(hxlim1,'String')) str2num(get(hxlim2,'String'))];
nc = handles.epoch_data.num_sensors;
for k = 1:length(handles.current_sensor)  
  ch = handles.current_sensor(k);
  fc = str2num(get(hprocf2,'String'));
  tfdata = ts_freqanalysis_fieldtrip(...
    ts_data_selection(handles.epoch_data,'channels',ch,'events',handles.evcodes(c)),...
    'savecomplex_flag',0,'lnfilter','yes','lpfilter',lpf,'lpfreq',fc,'foi',foi,...
    'toilim',toilim);
  nf = length(tfdata.timefreq.frequencies);  
  if ~isfield(handles,'timefreq_data') && isfield(handles,'epoch_data')
    handles.timefreq_data = rmfield(handles.epoch_data,{'epochs','sensor_info'});   
    handles.timefreq_data.sensor_info = init_sensor_info(handles.epoch_data.num_sensors);
    for j = 1:length(handles.epoch_data.epochs)
      nt = length(handles.epoch_data.epochs(j).time);
      handles.timefreq_data.timefreq(j).event_code  = handles.epoch_data.epochs(j).event_code;
      handles.timefreq_data.timefreq(j).num_trials  = handles.epoch_data.epochs(j).num_trials;
      handles.timefreq_data.timefreq(j).num_rejects = handles.epoch_data.epochs(j).num_rejects;
      handles.timefreq_data.timefreq(j).time        = handles.epoch_data.epochs(j).time;
      handles.timefreq_data.timefreq(j).frequencies = tfdata.timefreq.frequencies;  
      handles.timefreq_data.timefreq(j).power       = [];%nan(nc,nt,nf);
      handles.timefreq_data.timefreq(j).cmplx       = [];%nan(nc,nt,nf);
    end
  end
  nt = length(handles.timefreq_data.timefreq(c).time);
  try
    handles.timefreq_data.sensor_info(ch)         = tfdata.sensor_info;
  catch
    handles.timefreq_data.sensor_info(ch)         = orderfields(tfdata.sensor_info,handles.timefreq_data.sensor_info);
  end
  handles.timefreq_data.num_sensors             = length(handles.timefreq_data.sensor_info);
  tidx = nearest(handles.timefreq_data.timefreq(c).time,tfdata.timefreq.time(1)):...
         nearest(handles.timefreq_data.timefreq(c).time,tfdata.timefreq.time(end));
  fidx = nearest(handles.timefreq_data.timefreq(c).frequencies,tfdata.timefreq.frequencies(1)):...
         nearest(handles.timefreq_data.timefreq(c).frequencies,tfdata.timefreq.frequencies(end));  
  otix = 1:length(tfdata.timefreq.time);
  ofix = 1:length(tfdata.timefreq.frequencies);
  if length(tidx) ~= length(otix)
    [otoi tidx otix] = intersect(handles.timefreq_data.timefreq(c).time,tfdata.timefreq.time);
  end
  if length(fidx) ~= length(ofix)
    [ofoi fidx ofix] = intersect(handles.timefreq_data.timefreq(c).frequencies,tfdata.timefreq.frequencies);
  end
  if isfield(tfdata.timefreq,'cmplx')
    if size(handles.timefreq_data.timefreq(c).cmplx,2) < nt
      cmplx              = nan(1,nt,nf);
      cmplx(1,otix,ofix) = tfdata.timefreq.cmplx(1,otix,ofix);
      if isempty(handles.timefreq_data.timefreq(c).cmplx)
        handles.timefreq_data.timefreq(c).cmplx               = cmplx;
      else
        handles.timefreq_data.timefreq(c).cmplx(ch,tidx,fidx) = cmplx;
      end
    else
      handles.timefreq_data.timefreq(c).cmplx(ch,tidx,fidx) = tfdata.timefreq.cmplx(1,otix,ofix);
    end
  end
  if isfield(tfdata.timefreq,'power')
    if size(handles.timefreq_data.timefreq(c).power,2) < nt
      power              = nan(1,nt,nf);
      power(1,otix,ofix) = tfdata.timefreq.power(1,otix,ofix);
      if isempty(handles.timefreq_data.timefreq(c).power)
        handles.timefreq_data.timefreq(c).power               = power;      
      else
        handles.timefreq_data.timefreq(c).power(ch,tidx,fidx) = power;      
      end
    else
      handles.timefreq_data.timefreq(c).power(ch,tidx,fidx) = tfdata.timefreq.power(1,otix,ofix);
    end
  elseif isfield(tfdata.timefreq,'cmplx')
    handles.timefreq_data.timefreq(c).power = nan(size(handles.timefreq_data.timefreq(c).cmplx));
    handles.timefreq_data.timefreq(c).power(ch,tidx,fidx) = 2*abs(tfdata.timefreq.cmplx(1,otix,ofix)).^2;
  end  
end
toc
ymin = tfdata.timefreq.frequencies(1);
ymax = tfdata.timefreq.frequencies(end);
set(hplotpow,'Value',1);
set(happlay,'Value',0);
set(hylim1,'String',num2str(ymin));
set(hylim2,'String',num2str(ymax));
handles.plots.ymin = ymin;
handles.plots.ymax = ymax;
handles.current_datatype = 'timefreq_data';
clear tfdata
guidata(gcbo,handles);
update_figure(handles);
end

function ICAanalysis_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
tic 
% plottype: {'alltrials','activations'}
handles.epoch_data = ts_manualICA(handles.epoch_data,'plottype','activations',...
  'maxsteps',5,'ntrial',1,'ncomponents',80,'rescale',1); %'event_codes',handles.evcodes(handles.current_condition));
toc
update_figure(handles);
guidata(gcbo,handles);
end

function axis_xleft_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
tmin  = handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).time(1);
tmax  = handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).time(end);
x1    = str2num(get(hxlim1,'String'));  
x2    = str2num(get(hxlim2,'String'));  
if x1 == tmin, return; end
xlow  = x1 - (x2 - x1);
xhigh = x1;
if xlow < tmin
  xlow  = tmin;
  xhigh = tmin + (x2 - x1);
end
set(hxlim1,'String',num2str(xlow));
set(hxlim2,'String',num2str(xhigh));
axislimits_Callback(hObject,eventdata,handles);
end

function axis_xright_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
tmin  = handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).time(1);
tmax  = handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).time(end);
x1      = str2num(get(hxlim1,'String'));  
x2      = str2num(get(hxlim2,'String'));  
if x2 == tmax, return; end
xlow    = x2;
xhigh   = x2 + (x2 - x1);
if xhigh > tmax
  xlow  = tmax - (x2 - x1);
  xhigh = tmax;  
end
set(hxlim1,'String',num2str(xlow));
set(hxlim2,'String',num2str(xhigh));
axislimits_Callback(hObject,eventdata,handles);
end

function axislimits_Callback(hObject,eventdata,handles)
handles = guidata(gcbo);
handles.preproc.hprocbc1 = str2num(get(hprocbc1,'String'));
handles.preproc.hprocbc2 = str2num(get(hprocbc2,'String'));
handles.plots.xlow   = str2num(get(hxlim1,'String'));  
handles.plots.xhigh  = str2num(get(hxlim2,'String'));  
handles.plots.ylow   = str2num(get(hylim1,'String'));  
handles.plots.yhigh  = str2num(get(hylim2,'String'));  
handles.plots.maxmin = 0;
guidata(gcbo,handles);
handles = update_figure(handles);
end

function set_XLim_Callback(hObject,eventdata,handles)
[x y] = ginput(2);
if x(1) > x(2)
  tmp = x(1);
  x(1) = x(2);
  x(2) = tmp;
  clear tmp;
end
set(gca,'XLim',[x(1) x(2)]);
end
function reset_XLim_Callback(hObject,eventdata,handles)
set(gca,'XLim',[str2num(get(hxlim1,'String')) str2num(get(hxlim2,'String'))]);
end
function set_YLim_Callback(hObject,eventdata,handles)
[x y] = ginput(2);
if y(1) > y(2)
  tmp = y(1);
  y(1) = y(2);
  y(2) = tmp;
  clear tmp;
end
set(gca,'YLim',[y(1) y(2)]);  
end
function reset_YLim_Callback(hObject,eventdata,handles)
set(gca,'YLim',[str2num(get(hylim1,'String')) str2num(get(hylim2,'String'))]);  
end

function toggle_channel_Callback(hObject,eventdata,handles)
if ~(get(happlay,'Value') && get(hchtoggle,'Value'))
  set(hchtoggle,'Value',0); 
  return; 
end  
% hMeanFig = figure('Visible','on','Position',[1010,500,500,500]);
% hMeanAx  = axes('Parent',hMeanFig);
% axes(hplot)
% set(hMeanAx,'NextPlot','replacechildren')
% plot(hMeanAx,t,dat); grid on; xlabel('time (s)'); ylabel(handles.plots.ylabel); title('average');
% ylim = [min(dat(:)) max(dat(:))];
handles     = guidata(gcbo);
type        = handles.current_datatype;
field       = handles.current_datafield;
c           = handles.current_condition(1);
cfg.layout  = handles.layout_list{handles.current_layout};
lay         = prepare_layout(cfg);
badchans    = find([handles.(type).sensor_info.badchan]);
[x y]       = ginput(1);
while ~isempty(x)
  idx       = find(lay.pos(:,1) < x & lay.pos(:,2) < y);
  tmp       = lay.pos(idx,1) + lay.pos(idx,2); 
  idx       = idx(find(tmp == max(tmp)));
  if ismember(idx,badchans)
    badchans = setdiff(badchans,idx);
  else
    badchans = [badchans idx];
  end
  fprintf('bad channels: %s\n',num2str(badchans));
%   % plot mean
%   trl = setdiff(1:handles.(type).(field)(c).num_trials,badtrials);
%   dat = squeeze(handles.epoch_data.epochs(c).data(ch,:,trl));
%   if length(trl) > 1, dat = mean(dat,2); end
%   set(hMeanAx,'NextPlot','replacechildren')
%   plot(hMeanAx,t,dat); grid on;
%   xlabel('time (s)'); ylabel(handles.plots.ylabel); title('new average');
%   set(hMeanAx,'XLim',[handles.plots.xlow handles.plots.xhigh],'YLim',ylim);
% %   axis(hMeanAx,[handles.plots.xlow handles.plots.xhigh handles.plots.ylow handles.plots.yhigh]);
%   axes(hplot)
  [x y]     = ginput(1);  
end
[badchans sel2] = match_str({handles.(type).sensor_info.label},lay.label(badchans));

if ~isempty(badchans)
  [handles.(type).sensor_info(badchans).badchan] = deal(1);
  handles.rejects(c).badchans = sort(badchans);
else
  [handles.(type).sensor_info(1:length(handles.(type).sensor_info)).badchan] = deal(0);
  handles.rejects(c).badchans = [];  
end
guidata(gcbo,handles);
set(hchtoggle,'Value',0)
axislimits_Callback(hObject,eventdata,handles);
% close(hMeanFig)
end

function reject_Callback(hObject,eventdata,handles)
if get(happlay,'Value') || ~get(hplottrl,'Value')
  return;
end
hMeanFig = figure('Visible','on','Position',[1010,500,500,500]);
hMeanAx  = axes('Parent',hMeanFig);
axes(hplot)
handles = guidata(gcbo);
ch          = handles.current_sensor;
lay         = handles.trial_layout;
c           = handles.current_condition(1);
t           = handles.epoch_data.epochs(c).time;
badtrials   = handles.rejects(c).badtrials{handles.current_sensor};
trl         = setdiff(1:handles.epoch_data.epochs(c).num_trials,badtrials);
dat         = squeeze(handles.epoch_data.epochs(c).data(ch,:,trl));
if length(trl) > 1
  dat = mean(dat,2);
end
ylim = [min(dat(:)) max(dat(:))];
set(hMeanAx,'NextPlot','replacechildren')
plot(hMeanAx,t,dat); grid on;
xlabel('time (s)'); ylabel(handles.plots.ylabel); title('average');
[x y]       = ginput(1);
while ~isempty(x)
%   fprintf('%g trials marked BAD: %s\n',length(badtrials),num2str(badtrials));
  idx       = find(lay.pos(:,1) < x & lay.pos(:,2) < y);
  tmp       = lay.pos(idx,1) + lay.pos(idx,2); 
  idx       = idx(find(tmp == max(tmp)));
  if ismember(idx,badtrials)
    badtrials = setdiff(badtrials,idx);
    text(lay.pos(idx,1),lay.pos(idx,2)+lay.width(idx),...
      sprintf(' %0s\n ',['trial ' num2str(idx) ' (-)']),'Fontsize',8);
  else
    badtrials = [badtrials idx];
    text(lay.pos(idx,1),lay.pos(idx,2)+lay.width(idx),...
      sprintf(' %0s\n ',['trial ' num2str(idx) ' (x)']),'Fontsize',8);
  end
  fprintf('bad trials: %s\n',num2str(badtrials));
  % plot mean
  trl = setdiff(1:handles.epoch_data.epochs(c).num_trials,badtrials);
  dat = squeeze(handles.epoch_data.epochs(c).data(ch,:,trl));
  if length(trl) > 1, dat = mean(dat,2); end
  set(hMeanAx,'NextPlot','replacechildren')
  plot(hMeanAx,t,dat); grid on;
  xlabel('time (s)'); ylabel(handles.plots.ylabel); title('new average');
  set(hMeanAx,'XLim',[handles.plots.xlow handles.plots.xhigh],'YLim',ylim);
%   axis(hMeanAx,[handles.plots.xlow handles.plots.xhigh handles.plots.ylow handles.plots.yhigh]);
  axes(hplot)
  [x y]     = ginput(1);  
end
handles.rejects(c).badtrials{handles.current_sensor} = sort(badtrials);
guidata(gcbo,handles);
close(hMeanFig)
end

function BandPlot_Callback(hObject,eventdata,handles)
if ~get(hplotpow,'Value')
  return;
end
axes(hplot)
set(hplot,'NextPlot','replacechildren')
handles = guidata(gcbo);
ch      = handles.current_sensor;
if ~isfield(handles,'timefreq_data')
  fprintf('error: could not find power averages\n');
  return; 
elseif handles.timefreq_data.sensor_info(ch).badchan
  fprintf('structure does not contain power for this channel\n');
  return;
end  
flo     = str2num(get(hbandlo,'String'));
fhi     = str2num(get(hbandhi,'String'));
cs      = handles.current_condition;
t       = handles.timefreq_data.timefreq(cs(1)).time;
f       = handles.timefreq_data.timefreq(cs(1)).frequencies;
xmin    = str2num(get(hxlim1,'String'));  
xmax    = str2num(get(hxlim2,'String'));  
ymin    = str2num(get(hylim1,'String'));  
ymax    = str2num(get(hylim2,'String'));  

% tidx    = find(t >= xmin & t <= xmax);
% tsel    = t(tidx);
% xticks  = round(linspace(1,length(tsel),5));
% yticks  = round(linspace(1,length(fsel),5));
% xlabels = tsel(round(linspace(1,length(tsel),5)));
% ylabels = fliplr(fsel(round(linspace(1,length(fsel),5))));

if (flo==0 && fhi==0) || (flo > fhi)
  [x y] = ginput(2);
  flo = f(nearest(1:length(f),y(1)));
  fhi = f(nearest(1:length(f),y(2)));
  if flo > fhi
    tmp = flo;
    flo = fhi;
    fhi = tmp;
  end
end
fidx = find(f >= flo & f <= fhi);
tidx = find(t <= handles.preproc.hprocbc2 & t >= handles.preproc.hprocbc1);
if isempty(tidx)
  tidx = 1:length(t);
end
figure(200); set(gcf,'Position',[1010,375,600,675]);
subplot(2,1,1)
for k = 1:length(cs)
  c = cs(k);
  x = squeeze(handles.timefreq_data.timefreq(c).power(ch,:,:));
  x = (x - repmat(nanmean(x(tidx,:),1),[length(t) 1])) ./ repmat(nanstd(x(tidx,:),0,1),[length(t) 1]);
  x = mean(x(:,fidx),2);
  plot(t,x,handles.plots.graphcolor(k));
  if k == 1, hold on; end
end
ylabel(sprintf('mean power zscore (%g-%gHz)',flo,fhi));
title(handles.timefreq_data.sensor_info(ch).label);
set(gca,'Xlim',[xmin xmax]);
hold off
subplot(2,1,2)
for k = 1:length(cs)
  c = cs(k);
  x = nanmean(squeeze(handles.timefreq_data.timefreq(c).power(ch,:,fidx)),2);
  plot(t,x,handles.plots.graphcolor(k));
  if k == 1, hold on; end
end
ylabel(sprintf('mean power (%g-%gHz)',flo,fhi));
xlabel('time (s)')
set(gca,'Xlim',[xmin xmax]);
hold off
end

function close_Callback(hObject,eventdata,handles)
  handles = guidata(gcbo);
  if isfield(handles,'tempfiles')
    rmfiles = unique(handles.tempfiles);    
    for k = 1:length(rmfiles)
      if exist(rmfiles{k})
        unix(['rm ' rmfiles{k}]);
      end
    end
  end
  close all;
end

function handles = update_figure(handles)
set(hsensor,'Value',handles.current_sensor);
multi_flag = get(happlay,'Value');
if get(hploterp,'Value'),     type = 'wave';
elseif get(hplotpow,'Value'), type = 'power';
elseif get(hplotplv,'Value'), type = 'plv';
else
  disp('plot type not recognized.');
  return;
end
if get(hplotavg,'Value'),     trials_flag = 0;
elseif get(hplottrl,'Value'), trials_flag = 1;
else
  disp('plot type not recognized.');
  return;
end

if      multi_flag   && strcmp(type,'wave')  && trials_flag
  multiplot_wave_trials(handles);
elseif  multi_flag   && strcmp(type,'wave')  && ~trials_flag
  multiplot_wave_averages(handles);
elseif  multi_flag   && strcmp(type,'power') && ~trials_flag
  multiplot_pow_averages(handles);
elseif  multi_flag   && strcmp(type,'all')   && ~trials_flag
  multiplot_all_averages(handles);
elseif  ~multi_flag  && strcmp(type,'wave')  && trials_flag
  singleplot_wave_trials;
%   handles = singleplot_wave_trials(handles);
elseif  ~multi_flag  && strcmp(type,'wave')  && ~trials_flag
  singleplot_wave_averages(handles);
elseif  ~multi_flag  && strcmp(type,'power') && ~trials_flag
  singleplot_pow_averages(handles);
elseif  ~multi_flag  && strcmp(type,'power') && trials_flag
  singleplot_pow_trials(handles);
elseif  ~multi_flag  && strcmp(type,'all')   && ~trials_flag
  singleplot_all_averages(handles);
else
  disp('plot type not recognized.');
  return;
end
end    

function singleplot_wave_trials
handles = guidata(gcbo);
handles = update_datatype(handles,'epoch_data');
if ~isfield(handles,'epoch_data')
  fprintf('error: could not find waveform trials\n');
  return; 
end  
set(hplot,'NextPlot','replacechildren')
hold off

ch = handles.current_sensor;
cs = handles.current_condition;
% get maxmin over all trials
xmin = str2num(get(hxlim1,'String'));  
xmax = str2num(get(hxlim2,'String'));  
t    = handles.epoch_data.epochs(cs(1)).time;
xidc = nearest(t,xmin):nearest(t,xmax);
ymin = min(min(min(handles.epoch_data.epochs(cs).data(ch,xidc,:))));
ymax = max(max(max(handles.epoch_data.epochs(cs).data(ch,xidc,:))));
set(hylim1,'String',num2str(ymin));
set(hylim2,'String',num2str(ymax));

% TODO: add loop over conditions
ntrl = handles.epoch_data.epochs(cs(1)).num_trials;
ncol = ceil(sqrt(ntrl))+1;
nrow = ceil(sqrt(ntrl))+1;
sens = repmat(handles.epoch_data.sensor_info(ch),1,ntrl);
k = 0;
for i = 1:nrow
  for j = 1:ncol
    k = k + 1;
    if k > ntrl, break; end
    xpos = (j-1)/ncol;
    ypos = (nrow-i-1)/nrow;
    w    = 0.8 * 1/ncol;
    h    = 0.8 * 1/nrow;    
    x    = xpos + w*(t-xmin)/(xmax-xmin);
    y    = squeeze(handles.epoch_data.epochs(cs(1)).data(ch,:,k)); 
    y(y > ymax) = ymax;
    y(y < ymin) = ymin;
    y = ypos + h*(y-ymin)/(ymax-ymin);
    plot(x,y);
    % draw y-axis
    xs2 = xpos + w*([0 0]-xmin)/(xmax-xmin);
    ys2 = ypos + h*([ymin ymax]-ymin)/(ymax-ymin);
    plot(xs2,ys2,'k');
    % draw x-axis
    xs2 = xpos + w*([xmin xmax]-xmin)/(xmax-xmin);
    ys2 = ypos + h*([0 0]-ymin)/(ymax-ymin);   
    plot(xs2,ys2,'k');
    % add toggle label
    if ismember(k,handles.rejects(cs(1)).badtrials{ch})
      text(xpos,ypos+h,sprintf(' %0s\n ',['trial ' num2str(k) ' (x)']),'Fontsize',8);
    else
      text(xpos,ypos+h,sprintf(' %0s\n ',['trial ' num2str(k)]),'Fontsize',8);
    end    
    if k == 1, hold on; end
    sens(k).label = num2str(k);
  end
end
title([handles.epoch_data.sensor_info(ch).label])
hold off
axis off
handles.trial_layout = get_ordered_layout(sens,0);
guidata(gcbo,handles);
end

function singleplot_wave_averages(handles)
handles = guidata(gcbo);
if isfield(handles,'avg_data')
  handles = update_datatype(handles,'avg_data');
elseif isfield(handles,'epoch_data')
  handles = update_datatype(handles,'epoch_data');
else
  fprintf('error: could not find waveform averages\n');
  return;   
end
type  = handles.current_datatype;
field = handles.current_datafield;

set(hplot,'NextPlot','replacechildren')
ch = handles.current_sensor;
cs = handles.current_condition;
t  = handles.(type).(field)(cs(1)).time;
xmin = str2num(get(hxlim1,'String'));  
xmax = str2num(get(hxlim2,'String'));  
plotstat = 0;
if isfield(handles,'stat_data') && length(cs) == 2
  evcodes = handles.evcodes(cs);
  for k = 1:length(handles.stat_data.stats)
    if isequal(evcodes,[handles.stat_data.stats(k).event_code])
      plotstat = 1;
      statix   = k;
      break;
    end
  end
end
if plotstat && ismember(handles.(type).sensor_info(ch).label,handles.stat_data.stats(statix).label)
  [sel sel2] = match_str(handles.stat_data.stats(statix).label,{handles.(type).sensor_info(ch).label});
  st = handles.stat_data.stats(statix).time;
  xmin = max([xmin min(t) min(st)]);
  xmax = min([xmax max(t) max(st)]);
  set(hxlim1,'String',num2str(xmin));  
  set(hxlim2,'String',num2str(xmax));  
  handles.plots.xlow  = xmin;
  handles.plots.xhigh = xmax;
  xdix = nearest(t,xmin):nearest(t,xmax);
  xsix = nearest(st,xmin):nearest(st,xmax);
  t    = t(xdix);
  st   = st(xsix);  
  if ndims(handles.(type).(field)(cs(1)).data) > 2    
    trl1 = setdiff(1:handles.(type).(field)(cs(1)).num_trials,handles.rejects(cs(1)).badtrials{ch});
    trl2 = setdiff(1:handles.(type).(field)(cs(2)).num_trials,handles.rejects(cs(2)).badtrials{ch});
  else
    trl1 = 1;
    trl2 = 1;
  end
  x1   = squeeze(handles.(type).(field)(cs(1)).data(ch,xdix,trl1));
  x2   = squeeze(handles.(type).(field)(cs(2)).data(ch,xdix,trl2));
%   if handles.(type).(field)(cs(1)).num_trials > 1
  if length(trl1) > 1
    x1 = mean(x1,2);
  end
%   if handles.(type).(field)(cs(2)).num_trials > 1
  if length(trl2) > 1
    x2 = mean(x2,2);
  end
  hold off
  plot(t,x1,handles.plots.graphcolor(cs(1)),t,x2,handles.plots.graphcolor(cs(2)));
  hold on
  tt = t(squeeze(handles.stat_data.stats(statix).mask(sel,xsix)) > 0);
  if ~isempty(tt) %length(tt) > 0
%       % thicken x-axis where there are significant differences
%       scatter(tt,zeros(1,length(tt)),'k','filled');
    % shade regions of significant differences
    mask = squeeze(handles.stat_data.stats(statix).mask(sel,xsix)) > 0;
    n = length(mask);
    mask  = [0 mask(2:end-1).*(((mask(1:n-2)==0) + (mask(3:n)==0)) > 0) 0];
    edges = find(mask);
    for k = 1:length(edges)/2
      ix1 = t(edges(2*k-1));
      ix2 = t(edges(2*k));
      iy1 = handles.plots.ylow;
      iy2 = handles.plots.yhigh;
      fill([ix1 ix1 ix2 ix2],[iy1 iy2 iy2 iy1],get(hstatc,'String'),'EdgeColor','k',...
        'FaceAlpha',str2num(get(hstatt,'String')),'EdgeAlpha',str2num(get(hstatt,'String')));
    end
  end
else
  hold off  
  for k = 1:length(cs)
    if strcmp(type,'avg_data')
      trl = 1;
    else
      trl = setdiff(1:handles.(type).(field)(cs(k)).num_trials,handles.rejects(cs(k)).badtrials{ch});
    end
    x   = squeeze(handles.(type).(field)(cs(k)).data(ch,:,trl));
  %   if handles.(type).(field)(cs(1)).num_trials > 1
    if length(trl) > 1
      x = mean(x,2);
    end  
    plot(t,x,handles.plots.graphcolor(cs(k)));
    hold on
  end
end
axis([handles.plots.xlow handles.plots.xhigh handles.plots.ylow handles.plots.yhigh]);
xlabel('time (s)');
ylabel(handles.plots.ylabel); 
title([handles.(type).sensor_info(ch).label]);
grid on
hold off    
guidata(gcbo,handles);
end

function singleplot_pow_averages(handles)
handles = guidata(gcbo);
handles = update_datatype(handles,'timefreq_data');
ch      = handles.current_sensor;
chtf    = ch;
if ~isfield(handles,'timefreq_data')
  fprintf('error: could not find power averages\n');
  return; 
elseif handles.timefreq_data.sensor_info(ch).badchan
  fprintf('structure does not contain power for this channel\n');
  return;
end  
if isfield(handles,'epoch_data') && ismember(handles.epoch_data.sensor_info(ch).label,{handles.timefreq_data.sensor_info.label})
  singleplot_all_averages(handles);
%   return;
end
% [ jk] = match_str({handles.timefreq_data.sensor_info.label},{handles.epoch_data.sensor_info(ch).label});
c       = handles.current_condition(1);
ctf     = find([handles.timefreq_data.timefreq.event_code] == handles.evcodes(c));
t       = handles.timefreq_data.timefreq(c).time;
f       = handles.timefreq_data.timefreq(c).frequencies;
xmin    = str2num(get(hxlim1,'String'));  
xmax    = str2num(get(hxlim2,'String'));  
ymin    = str2num(get(hylim1,'String'));  
ymax    = str2num(get(hylim2,'String'));  

fidx    = find(f >= ymin & f <= ymax);
fsel    = f(fidx);
tidx    = find(t >= xmin & t <= xmax);
tsel    = t(tidx);
xbc1 = handles.preproc.hprocbc1;
xbc2 = handles.preproc.hprocbc2;
tidz = find(t <= xbc2 & t >= xbc1);
if isempty(tidz) 
  tidz = 1:length(t); 
end
data    = squeeze(handles.timefreq_data.timefreq(ctf).power(chtf,:,:));

xticks  = round(linspace(1,length(tsel),5));
yticks  = round(linspace(1,length(fsel),5));
xlabels = tsel(round(linspace(1,length(tsel),5)));
ylabels = fliplr(fsel(round(linspace(1,length(fsel),5))));

% calc power zscore
tfr     = data(tidx,fidx);

if all(all(isnan(tfr)))
  TFanalysis_Callback;
  return;
end

tfz     = (tfr - repmat(nanmean(data(tidz,fidx),1),  [length(tsel) 1])) ./ ...
                 repmat(nanstd (data(tidz,fidx),0,1),[length(tsel) 1]);
% % calc freq correction factor
% frq     = repmat(fliplr(fsel)',1,size(tfr,1));

% plot power zscore
set(hplot,'NextPlot','replacechildren')
hold off
imagesc(flipud(tfz'));
title('power zscore');         caxis(handles.plots.caxis); 
set(gca,'XTick',xticks);       set(gca,'YTick',yticks);  
set(gca,'XTickLabel',xlabels); set(gca,'YTickLabel',ylabels);
grid on;                       ylabel('frequency (Hz)');
xlabel('time (s)');
title([handles.timefreq_data.sensor_info(chtf).label ' (power zscore)']);
% axis tight

% todo: add time-freq stats
if isfield(handles,'tfstat_data')
  mask = handles.tfstat_data.stats(1).mask;
  %% HOW TO SHOW?? LOOK AT ERP STATS....
end

hold off
guidata(gcbo,handles);
end

function multiplot_wave_averages(handles)
handles = guidata(gcbo);
if isfield(handles,'avg_data')
  handles = update_datatype(handles,'avg_data');
elseif isfield(handles,'epoch_data')
  handles = update_datatype(handles,'epoch_data');
else
  fprintf('error: could not find waveform averages\n');
  return;   
end
type  = handles.current_datatype;
field = handles.current_datafield;

set(hplot,'NextPlot','replacechildren')
ch = handles.current_sensor;
cs = handles.current_condition;
t  = handles.(type).(field)(cs(1)).time;
xmin = str2num(get(hxlim1,'String'));  
xmax = str2num(get(hxlim2,'String'));  
title ''
avg = handles.(type);
avg.averages = avg.(field);
if ~strcmp(field,'averages')
  avg = rmfield(avg,field);
end

plotstat = 0;
if isfield(handles,'stat_data') && length(cs) == 2
  evcodes = handles.evcodes(cs);
  for k = 1:length(handles.stat_data.stats)
    if isequal(evcodes,[handles.stat_data.stats(k).event_code])
      plotstat = 1;
      statix   = k;
      break;
    end
  end
end
if plotstat
  st = handles.stat_data.stats(statix).time;
  xmin = max([xmin min(t) min(st)]);
  xmax = min([xmax max(t) max(st)]);
  set(hxlim1,'String',num2str(xmin));  
  set(hxlim2,'String',num2str(xmax));  
  handles.plots.xlow  = xmin;
  handles.plots.xhigh = xmax;
  xdix = nearest(t,xmin):nearest(t,xmax);
  xsix = nearest(st,xmin):nearest(st,xmax);
  t    = t(xdix);
  st   = st(xsix);  
  mask = handles.stat_data.stats(statix).mask(:,xsix);
  
  for k = 1:length(avg.averages)
    if ~ismember(avg.averages(k).event_code,handles.evcodes(cs)), continue; end
    avg.averages(k).data = avg.averages(k).data(:,xdix,:);
    avg.averages(k).time = t;
  end    
%   cfg.layout = handles.layout_list{handles.current_layout};
%   lay        = prepare_layout(cfg);
%   [dch sch] = match_str({avg.sensor_info.label},{handles.stat_data.sensor_info.label});
%   [lch sch] = match_str(lay.label,{handles.stat_data.sensor_info.label});  
  % update time limits in mask
  ts_ezmultiplot(avg,'newfig',0,'baseline','no','axes','no','zerolines','no',...
    'xlim',[handles.plots.xlow handles.plots.xhigh],'ylim','maxmin',...
    'layout',handles.layout_list{handles.current_layout},...
    'events',handles.evcodes(handles.current_condition),...
    'statcolor',get(hstatc,'String'),'transparency',str2num(get(hstatt,'String')),...
    'mask',mask,'graphcolor',handles.plots.graphcolor(cs));
else
  ts_ezmultiplot(ts_data_selection(avg,'toilim',[xmin xmax]),'newfig',0,'baseline','no','axes','yes','zerolines','no',...
    'xlim',[handles.plots.xlow handles.plots.xhigh],'ylim','maxmin',...
    'layout',handles.layout_list{handles.current_layout},...
    'events',handles.evcodes(handles.current_condition),'graphcolor',handles.plots.graphcolor(cs));
end    
hold off
guidata(gcbo,handles);
end

function multiplot_pow_averages(handles)
handles = guidata(gcbo); 
handles = update_datatype(handles,'timefreq_data');
if ~isfield(handles,'timefreq_data')
  fprintf('error: could not find power averages\n');
  return; 
end  
set(hplot,'NextPlot','replacechildren')
cs = handles.current_condition;
t  = handles.timefreq_data.timefreq(cs(1)).time;
xmin = str2num(get(hxlim1,'String'));  
xmax = str2num(get(hxlim2,'String'));  
ymin = str2num(get(hylim1,'String'));  
ymax = str2num(get(hylim2,'String')); 
avg  = ts_data_selection(handles.timefreq_data,'toilim',[xmin xmax],'foilim',[ymin ymax]);
title ''
% prepare stats
if isfield(handles,'tfstat_data')
  mask = handles.tfstat_data.stats(1).mask;
  % update time limits in mask
  mask = mask(:,find(t >= handles.plots.xlow & t <= handles.plots.xhigh));
  [channels jnk] = match_str({avg.sensor_info.label},{handles.stat_data.sensor_info.label});
  ts_ezmultiplot(avg,'channels',channels,'newfig',0,'baseline','no','axes','no','zerolines','yes',...
    'xlim',[handles.plots.xlow handles.plots.xhigh],'zlim','maxmin','ylim','maxmin',...
    'layout',handles.layout_list{handles.current_layout},...
    'events',handles.evcodes(handles.current_condition),...
    'statcolor',get(hstatc,'String'),'transparency',str2num(get(hstatt,'String')),...
    'mask',mask,'graphcolor',handles.plots.graphcolor(cs));
else
  ts_ezmultiplot(avg,'newfig',0,'baseline','zscore','axes','no','zerolines','yes',...
    'xlim',[xmin xmax],'zlim',handles.plots.caxis,'ylim','maxmin',...
    'layout',handles.layout_list{handles.current_layout},...
    'events',handles.evcodes(handles.current_condition),'graphcolor',handles.plots.graphcolor(cs));
end    
hold off
guidata(gcbo,handles);
end

function singleplot_all_averages(handles)
set(hplot,'NextPlot','replacechildren')
handles = guidata(gcbo);  
handles = update_datatype(handles,'timefreq_data');
if isfield(handles,'epoch_data')
  type  = 'epoch_data';
  field = 'epochs';
elseif isfield(handles,'avg_data')
  type  = 'avg_data';
  field = 'averages';
else
  fprintf('error: could not find waveforms\n');
end

ch        = handles.current_sensor;
[chtf jk] = match_str({handles.timefreq_data.sensor_info.label},{handles.(type).sensor_info(ch).label});
c       = handles.current_condition(1);
ctf     = find([handles.timefreq_data.timefreq.event_code] == handles.evcodes(c));
t       = handles.timefreq_data.timefreq(c).time;
f       = handles.timefreq_data.timefreq(c).frequencies;
xmin    = str2num(get(hxlim1,'String'));  
xmax    = str2num(get(hxlim2,'String'));  
ymin    = str2num(get(hylim1,'String'));  
ymax    = str2num(get(hylim2,'String'));  

fidx    = find(f >= ymin & f <= ymax);
fsel    = f(fidx);
tidx    = find(t >= xmin & t <= xmax);
tsel    = t(tidx);
tidz    = find(t <= handles.preproc.hprocbc2 & t >= handles.preproc.hprocbc1);
% tidz    = find(t <=0 & t >= str2num(get(hxlim1,'String'))); 
if isempty(tidz)
  tidz = 1:length(t); 
end
data    = squeeze(handles.timefreq_data.timefreq(ctf).power(chtf,:,:));

xticks  = round(linspace(1,length(tsel),5));
yticks  = round(linspace(1,length(fsel),5));
xlabels = tsel(round(linspace(1,length(tsel),5)));
ylabels = fliplr(fsel(round(linspace(1,length(fsel),5))));

% calc power plot
tfr     = data(tidx,fidx);

% calc baseline corrected plot
tfz     = (tfr - repmat(nanmean(data(tidz,fidx),1),  [length(tsel) 1])) ./ ...
                 repmat(nanstd (data(tidz,fidx),0,1),[length(tsel) 1]);
  
% calc freq correction factor
frq     = repmat(fliplr(fsel)',1,size(tfr,1));

% calc waveform plot
x       = squeeze(handles.(type).(field)(c).data(ch,tidx,:));
if size(x,2) > 1
  x = mean(x,2);
end

figure(100); set(gcf,'Visible','on','Position',[1010,500,300,725]);
subplot(4,1,1),imagesc(flipud(tfr')); title('power');
set(gca,'XTick',xticks);       set(gca,'YTick',yticks);  
set(gca,'XTickLabel',xlabels); set(gca,'YTickLabel',ylabels);
grid on;                       ylabel('frequency (Hz)');        colorbar;
% subplot(4,1,2),imagesc(flipud(tfr'.*(frq.^3))); title('power * freq^3');
subplot(4,1,2),imagesc(flipud(tfr'.*(frq))); title('power * freq');
set(gca,'XTick',xticks);       set(gca,'YTick',yticks);  
set(gca,'XTickLabel',xlabels); set(gca,'YTickLabel',ylabels);
grid on;                       ylabel('frequency (Hz)');        colorbar;
subplot(4,1,3),imagesc(flipud(tfz')); title('power zscore'); 
caxis(handles.plots.caxis); 
set(gca,'XTick',xticks);       set(gca,'YTick',yticks);  
set(gca,'XTickLabel',xlabels); set(gca,'YTickLabel',ylabels);
grid on;                       ylabel('frequency (Hz)');
subplot(4,1,4),plot(tsel,x);
set(gca,'Xlim',[handles.plots.xlow handles.plots.xhigh]);
grid on;
xlabel('time (s)'); ylabel(handles.plots.ylabel); 
title([handles.timefreq_data.sensor_info(chtf).label ', condition ' num2str(handles.timefreq_data.timefreq(ctf).event_code)]);
hold off    
guidata(gcbo,handles);
axes(hplot)
end

function singleplot_pow_trials(handles)
handles = guidata(gcbo);  
handles = update_datatype(handles,'timefreq_data');
if length(handles.current_condition) > 1
  fprintf('Select only one condition to run TF trial analysis.\n');
  return;
else
  c = handles.current_condition;
end
if length(handles.current_sensor) > 1
  fprintf('Select only one channel to run TF trial analysis.\n');
  return;
else
  ch = handles.current_sensor;
end
set(hplot,'NextPlot','replacechildren')
hold off
fc = str2num(get(hprocf2,'String'));
if get(hfoimeg,'Value')
  foi = 'meg';
elseif get(hfoieeg,'Value')
  foi = 'ieeg';
else
  foi = 'ieeg';
end
if ~isfield(handles,'timefreq_trials')  || ...
  ~ismember(handles.timefreq_trials.sensor_info(ch).label,{handles.timefreq_trials.sensor_info.label}) || ...
  ~ismember(handles.timefreq_trials.timefreq(c).event_code,[handles.timefreq_trials.timefreq.event_code])
  if get(hprocbp,'Value'),  lpf='yes'; else lpf='no'; end
  tfdata = ts_freqanalysis_fieldtrip(...
  ts_data_selection(handles.epoch_data,'channels',ch,'events',handles.evcodes(c)),...
  'save_flag',0,'trials_flag',1,'lnfilter',lpf,'lpfilter','yes','lpfreq',fc,'foi',foi);  
  if ~isfield(handles,'timefreq_trials')
    handles.timefreq_trials = tfdata;
  return;
  else
    try handles.timefreq_trials.sensor_info(ch) = tfdata.sensor_info(1); end
    try handles.timefreq_trials.timefreq(end+1) = tfdata.timefreq; end
  end
else
  tfdata = handles.timefreq_trials;
  tfdata.timefreq = tfdata.timefreq(c);
  tfdata.timefreq.power = tfdata.timefreq.power(ch,:,:,:);
end
xmin = str2num(get(hxlim1,'String'));  
xmax = str2num(get(hxlim2,'String'));  
ymin = tfdata.timefreq.frequencies(1);   
ymax = tfdata.timefreq.frequencies(end); 
t    = tfdata.timefreq.time;
f    = tfdata.timefreq.frequencies;
ntrl = tfdata.timefreq.num_trials;
ncol = ceil(sqrt(ntrl))+1;
nrow = ceil(sqrt(ntrl))+1;
sens = repmat(tfdata.sensor_info,1,ntrl);
if isfield(tfdata.timefreq,'cmplx')
  data = squeeze(2*abs(tfdata.timefreq.cmplx).^2);  % time x freq x trials
else
  data = squeeze(tfdata.timefreq.power);
end
clear tfdata;
zmin = min(min(min(data(:))));
zmax = max(max(max(data(:))));

k = 0;
for i = 1:nrow
  for j = 1:ncol
    k = k + 1;
    if k > ntrl, break; end
    xpos = (j-1)/ncol;
    ypos = (nrow-i-1)/nrow;
    w    = 0.8 * 1/ncol;
    h    = 0.8 * 1/nrow;    
    x    = xpos + w*(t-xmin)/(xmax-xmin);
    y    = ypos + h*(f-ymin)/(ymax-ymin);
    imagesc(x,y,flipud(squeeze(data(:,:,k))'),[zmin zmax]);
    if k == 1, hold on; end
    % add trial labels
    if issubfield(handles,'rejects.powtrials') && ismember(k,handles.rejects(c).powtrials{ch})
      text(xpos,ypos,sprintf(' %0s\n ',['trial ' num2str(k) ' (x)']),'Fontsize',8);
    else
      text(xpos,ypos,sprintf(' %0s\n ',['trial ' num2str(k)]),'Fontsize',8);
    end    
    % add zerolines
    x2 = xpos + w*([0 0]-xmin)/(xmax-xmin);
    y2 = ypos + h*([ymin ymax]-ymin)/(ymax-ymin);
    plot(x2,y2,'k');
    sens(k).label = num2str(k);
  end
end
title([handles.timefreq_trials.sensor_info(ch).label])
axis tight
axis off
hold off
handles.trial_layout = get_ordered_layout(sens,0);
handles.trial_power  = data;
clear data;
guidata(gcbo,handles);  
end

function multiplot_wave_trials(handles)
% handles = guidata(gcbo);  
% handles = update_datatype(handles,'epoch_data');
% guidata(gcbo,handles);
end

function multiplot_all_averages(handles)
% handles = guidata(gcbo); 
% handles = update_datatype(handles,'both');
% guidata(gcbo,handles);
end

function sens = init_sensor_info(nchan)
[sens(1:nchan).label]      = deal('');
[sens(1:nchan).typestring] = deal('');
[sens(1:nchan).type]       = deal(1);
[sens(1:nchan).kind]       = deal(2);
% [sens(1:nchan).badchan]    = deal(0);
[sens(1:nchan).badchan]    = deal(1);
[sens(1:nchan).lognum]     = deal([]);
[sens(1:nchan).loc]        = deal([]);  
end

function handles = set_layouts(handles)
% handles           = guidata(gcbo); 
[layoutlist tempflag]        = get_ordered_layout(handles.(handles.current_datatype).sensor_info);
layoutlist = {layoutlist};
if tempflag
  handles.tempfiles{end+1} = layoutlist{end};
end
% try
  if handles.iEEG
    % ieeg data
      [jnk user] = unix('whoami');
      user       = user(1:end-1);
      [lays tempflag] = write_layout_ieeg(sprintf('/home/%s',user),handles.(handles.current_datatype).sensor_info,10,10);
      if ~iscell(lays), lays = {lays}; end
      if tempflag
        handles.tempfiles = {handles.tempfiles{:} lays{:}};
      end      
      layoutlist = {layoutlist{:} lays{:}};
  else
    % meeg data
      files = dir('/neuro/lout');
      files = files(~[files.isdir]);
      for k = 1:length({files.name})
        if findstr(files(k).name,'.lout')
          layoutlist = {layoutlist{:} files(k).name};
        end
      end
  end
% end
set(hlay,'String',layoutlist);
set(hlay,'Value',1);
handles.current_layout = 1;
handles.layout_list = layoutlist;
% look for condition key
handles.cond_mat_file = [];
handles.cond_csv_file = [];
recentfile = get(hfiletext,'String');
recentpath = fileparts(recentfile);
if exist(fullfile(recentpath,'cond_key.mat'),'file')
  handles.cond_mat_file = fullfile(recentpath,'cond_key.mat');
  cond_key = get_cond_key(handles.cond_mat_file);
  idx = find(ismember([cond_key.event_code],handles.evcodes));
  handles.cond_evcode = [cond_key.event_code];
  handles.cond_labels = [cond_key.name];
  handles.cond_labels = handles.cond_labels(idx);  
  handles.cond_evcode = handles.cond_evcode(idx);  
elseif exist(fullfile(recentpath,'conditionkey.csv'),'file')
  handles.cond_csv_file = fullfile(recentpath,'conditionkey.csv');
  cond_key = get_cond_key(handles.cond_mat_file);
  idx = find(ismember([cond_key.event_code],handles.evcodes));
  handles.cond_evcode = [cond_key.event_code];
  handles.cond_labels = [cond_key.name];
  handles.cond_labels = handles.cond_labels(idx);  
  handles.cond_evcode = handles.cond_evcode(idx);
else
  handles.cond_labels = [];
  handles.cond_evcode = [];
end
% guidata(gcbo,handles);
end

function [layfile tempflag] = get_ordered_layout(sens,outtype)
  tempflag = 0;
  if ~exist('outtype','var')
    outtype = 1; % write layout to file
  end
  nchan = length({sens.label});
  ncol = ceil(sqrt(nchan))+1;
  nrow = ceil(sqrt(nchan))+1;
  k = 0;
  for i=1:nrow
    for j=1:ncol
      k = k+1;
      if k<=nchan
        x = (j-1)/ncol;
        y = (nrow-i-1)/nrow;
        lay.pos(k,:) = [x y];
        lay.width(k,1)  = 0.8 * 1/ncol;
        lay.height(k,1) = 0.8 * 1/nrow;
      end
    end
  end
  lay.label = {sens.label};

  lay.label{end+1}  = 'SCALE';
  lay.width(end+1)  = mean(lay.width);
  lay.height(end+1) = mean(lay.height);
  x = (ncol-2)/ncol;
  y = 0/nrow;
  lay.pos(end+1,:) = [x y];

  lay.label{end+1}  = 'COMNT';
  lay.width(end+1)  = mean(lay.width);
  lay.height(end+1) = mean(lay.height);
  x = (ncol-1)/ncol;
  y = 0/nrow;
  lay.pos(end+1,:) = [x y];  
  
  if outtype
    % write temporary layout file
    [jnk user] = unix('whoami');
    user       = user(1:end-1);
    layfile = sprintf('/home/%s/ordered.lay',user);
    if ~exist(layfile,'file')
      tempflag = 1;
    else
      layfile = sprintf('/home/%s/ordered_%s.lay',user,datestr(now,30)); 
      tempflag = 1;
    end
    fid = fopen(layfile,'wt');
    for k = 1:length(lay.label)
      fprintf(fid,'%g %f %f %f %f %s\n',k,lay.pos(k,1),lay.pos(k,2),lay.width(k),lay.height(k),lay.label{k});
    end
    fclose(fid);
  else
    layfile = lay;
  end
end

function handles = update_datatype(handles,datatype)
if isfield(handles,'current_datatype') && strcmp(handles.current_datatype,datatype)
  if str2num(get(hylim1,'String'))==0 && str2num(get(hylim2,'String'))==0  
    reset_axis_limits = 2;
  else
    reset_axis_limits = 0;
  end
else
  reset_axis_limits = 1;
end
switch datatype
  case 'avg_data'
    handles.current_datatype  = 'avg_data';
    handles.current_datafield = 'averages';
    handles.current_zparam    = 'data';
    set(hploterp,'Value',1);
    set(hplotavg,'Value',1)
    if isfield(handles,'epoch_data')
      handles = rmfield(handles,'epoch_data');
    end
  case 'epoch_data'
    handles.current_datatype  = 'epoch_data';
    handles.current_datafield = 'epochs';
    handles.current_zparam    = 'data';
    set(hploterp,'Value',1);
    if isfield(handles,'avg_data')
      handles = rmfield(handles,'avg_data');
    end    
  case 'timefreq_data'
    handles.current_datatype  = 'timefreq_data';
    handles.current_datafield = 'timefreq';
    handles.current_zparam    = 'power';    
    set(hplotpow,'Value',1);
    if isfield(handles.timefreq_data.timefreq,'cmplx') && ~isfield(handles.timefreq_data.timefreq,'power')
      for k = 1:length(handles.timefreq_data.timefreq)
        handles.timefreq_data.timefreq(k).power = 2*abs(handles.timefreq_data.timefreq(k).cmplx).^2;    
      end
    end
  case 'timefreq_trials'
    handles.current_datatype  = 'timefreq_trials';
    handles.current_datafield = 'timefreq';
    handles.current_zparam    = 'power';      
    set(hplotpow,'Value',1);
    if isfield(handles.timefreq_trials.timefreq,'cmplx') && ...
      (~isfield(handles.timefreq_trials.timefreq,'power') || isempty(handles.timefreq_trials.timefreq.power))
      for k = 1:length(handles.timefreq_trials.timefreq)
        handles.timefreq_trials.timefreq(k).power = 2*abs(handles.timefreq_trials.timefreq(k).cmplx).^2;    
      end
    end
    if ~isfield(handles,'timefreq_data')
      handles.timefreq_data = handles.timefreq_trials;
      for k = 1:length(handles.timefreq_data.timefreq)
        try
          handles.timefreq_data.timefreq(k).power = nanmean(handles.timefreq_trials.timefreq(k).power,4);
        end
        try
          handles.timefreq_data.timefreq(k).cmplx = nanmean(handles.timefreq_trials.timefreq(k).cmplx,4);
        end
      end
    end
  case 'both'
%     handles.current_datatype  = {'epoch_data' 'timefreq_data'};
%     handles.current_datafield = {'epochs'     'timefreq'     };
%     handles.current_zparam    = {'data'       'power'        };
  otherwise
    disp('datatype not recognized.');
end
if handles.plots.maxmin || reset_axis_limits==1
  if get(happlay,'Value')
    if strcmp(handles.current_datafield,'timefreq')
      ymin = min(handles.(handles.current_datatype).(handles.current_datafield)(1).frequencies);
      ymax = max(handles.(handles.current_datatype).(handles.current_datafield)(1).frequencies);
    else  
      ymin = min(min(squeeze(mean(handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).(handles.current_zparam),3))));
      ymax = max(max(squeeze(mean(handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).(handles.current_zparam),3))));  
    end
  else
    if strcmp(handles.current_datafield,'timefreq')
      ymin = min(handles.(handles.current_datatype).(handles.current_datafield)(1).frequencies);
      ymax = max(handles.(handles.current_datatype).(handles.current_datafield)(1).frequencies);
    else  
      ymin = min(min(squeeze(mean(handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).(handles.current_zparam)(handles.current_sensor,:,:),3))));
      ymax = max(max(squeeze(mean(handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).(handles.current_zparam)(handles.current_sensor,:,:),3))));  
    end    
  end
  set(hylim1,'String',num2str(ymin));
  set(hylim2,'String',num2str(ymax));
  handles.plots.ylow  = str2num(get(hylim1,'String'));  
  handles.plots.yhigh = str2num(get(hylim2,'String'));  
  if handles.plots.maxmin
    xmin = min(handles.(handles.current_datatype).(handles.current_datafield)(1).time);
    xmax = max(handles.(handles.current_datatype).(handles.current_datafield)(1).time);  
    set(hxlim1,'String',num2str(xmin));
    set(hxlim2,'String',num2str(xmax));
    handles.plots.xlow  = str2num(get(hxlim1,'String'));  
    handles.plots.xhigh = str2num(get(hxlim2,'String')); 
  end 
elseif reset_axis_limits==2
  if strcmp(handles.current_datafield,'timefreq')
    ymin = min(handles.(handles.current_datatype).(handles.current_datafield)(1).frequencies);
    ymax = max(handles.(handles.current_datatype).(handles.current_datafield)(1).frequencies);
  else  
    ymin = min(min(squeeze(mean(handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).(handles.current_zparam)(handles.current_sensor,:,:),3))));
    ymax = max(max(squeeze(mean(handles.(handles.current_datatype).(handles.current_datafield)(handles.current_condition(1)).(handles.current_zparam)(handles.current_sensor,:,:),3))));  
  end
  set(hylim1,'String',num2str(ymin));
  set(hylim2,'String',num2str(ymax));
  handles.plots.ylow  = str2num(get(hylim1,'String'));  
  handles.plots.yhigh = str2num(get(hylim2,'String'));    
end
end

end