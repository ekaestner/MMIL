function visualizer(varargin)
% cont_data
% UserData.Raw.data/time = []; % orig data
% UserData.Proc.map = []; % [ProcSamples(:) OrigSamples(:)]
% UserData.Hilb.map = []; % Hilbert transform-based analytic signal
% UserData.TFR.map  = [];
% UserData.FFT.map  = [];
% UserData.Events   = []; % could be triggers or feature detection times
InitializeGUI;
% parse inputs
if nargin >= 1
  iscont = any(strmatch('cont',varargin(strmatch('char',cellfun(@(x)class(x),varargin,'uniformoutput',false)))));
  for k = 1:length(varargin)
    if isstruct(varargin{k})
      if iscont
        Load_Data([],[],varargin{k},'continuous');
      else
        Load_Data([],[],varargin{k});
      end
    elseif isnumeric(varargin{k})
      if iscont
        Load_Data([],[],varargin{k},'continuous');
      else
        Load_Data([],[],varargin{k});
      end
    elseif ischar(varargin{k}) && ~strmatch(varargin{k},'continuous')
      % convert remaining key/value pairs => parms structure
      % ...
      break
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALLBACKS & SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeGUI
% need tags: bp & bpfreq for CalcProc's call to ts_preproc

% create figure + menu
FigTag = 'MainFig1'; k = 1;
while any(findobj('tag',FigTag))
  k = k + 1;
  FigTag = sprintf('MainFig%g',k);
end
MainFig = figure('tag',FigTag,'Visible','off','WindowButtonMotionFcn',@ShowValue,'WindowScrollWheelFcn',@ZoomFunction,'WindowKeyPressFcn',@arrowpress);
set(MainFig,'MenuBar','none');
file_m    = uimenu(MainFig  ,'Label','File');
view_m    = uimenu(MainFig  ,'Label','View');
uimenu(file_m,'Label','Load data','Callback',@Load_Data);
uimenu(file_m,'Label','Exit','Callback',['global options; clear options; close(findobj(''tag'',''' FigTag '''));']);
% uimenu(view_m,'Label','imagesc','Callback',

% initialize panels
BorderWidth = .2;
BorderType  = 'line';%none,etchedin,etchedout,beveledin,beveledout,line
p1  = uipanel('Title','','FontSize',10,'BackgroundColor','white','Position',[0 .05 1 1]   ,'tag','plots','BorderWidth',BorderWidth,'BorderType',BorderType); % plots
p2  = uipanel('Title',''   ,'FontSize',10,'BackgroundColor','white','Position',[0 0 .1 .1],'BorderWidth',BorderWidth,'BorderType',BorderType  ); % time
% p2  = uipanel('Title','T'    ,'FontSize',10,'BackgroundColor','white','Position',[0 .1 .1 .1]  );
% p3  = uipanel('Title','F'    ,'FontSize',10,'BackgroundColor','white','Position',[0 0 .1 .1]   );
p4  = uipanel('Title','','FontSize',10,'BackgroundColor','white','Position',[.1 .05 .1 .05],'BorderWidth',BorderWidth,'BorderType',BorderType ); % scale
p5  = uipanel('Title','','FontSize',10,'BackgroundColor','white','Position',[.1 0 .1 .05],'BorderWidth',BorderWidth,'BorderType',BorderType  ); % events
p6  = uipanel('Title',''  ,'FontSize',10,'BackgroundColor','white','Position',[.2 .05 .1 .05],'BorderWidth',BorderWidth,'BorderType',BorderType ); % lay
p7  = uipanel('Title',''  ,'FontSize',10,'BackgroundColor','white','Position',[.2 0 .1 .05],'BorderWidth',BorderWidth,'BorderType',BorderType  ); % loc
p8  = uipanel('Title',''  ,'FontSize',10,'BackgroundColor','white','Position',[.3 0 .1 .1],'BorderWidth',BorderWidth,'BorderType',BorderType  ); % mri
p9  = uipanel('Title','','FontSize',10,'BackgroundColor','white','Position',[.6 0 .4 .1],'BorderWidth',BorderWidth,'BorderType',BorderType  ); % params
p10 = uipanel('Title','','FontSize',10,'BackgroundColor','white','Position',[.4 0 .2 .1],'BorderWidth',BorderWidth,'BorderType',BorderType  ); % channels

% panel 1 (plots)
ax(1) = subplot(1,1,1,'Parent',p1,'tag','tplot'); axis on; box on; title(''); set(ax(1),'color','red','XTickLabel','','YTickLabel','');
% ax(2) = subplot(2,1,2,'Parent',p1,'tag','tplot'); axis on; box on; title('plot 2'); set(ax(2),'color','blue','Position',[.6 .1 .2 .2],'XTickLabel','','YTickLabel','');

% panel 2 (time-control)
geometry = {[1 1] [1 1] 1};
geomvert = [.5 1 .5];
listui   = { ...
           {'parent',p2,'style','pushbutton','tag','tleft' ,'string','<==','Callback',{@ShiftTime,'left' }} ...
           {'parent',p2,'style','pushbutton','tag','tright','string','==>','Callback',{@ShiftTime,'right'}} ...
           {'parent',p2,'style','edit'  ,'tag','t0','Callback',{@UpdateManager}} ...
           {'parent',p2,'style','edit'  ,'tag','tf','Callback',{@UpdateManager}} ...
           {'parent',p2,'style','slider','tag','tslide','Callback',{@ShiftTime,'slide'}} ...           
%            {'parent',p2,'style','text','string','Time','BackgroundColor','k','ForegroundColor','k'} ...           
           };
supergui(MainFig,geometry,geomvert,listui{:});

% panel 3 (freq-control)
% panel 4 (y-control)
geometry = {[1] [1]};
geomvert = [1 1];
listui = { ...
         {'parent',p4,'style','checkbox','tag','autoscale','backgroundcolor','w','value',0,'string','autoscale','Callback',{@UpdateManager}} ...
         {'parent',p4,'style','checkbox','tag','eventchk','backgroundcolor','w','value',0,'string','events','Callback',{@UpdateManager}} ...
         };        
supergui(MainFig,geometry,geomvert,listui{:});

% panel 5 (events)
geometry = {[1 1],[1 1 1]};
geomvert = [1 1];
listui   = { ...
%            {'parent',p5,'style','pushbutton','string','Load events','Callback',@Load_Events} ...
           {'parent',p5,'style','pushbutton','string','Load','Callback',@Load_Events} ...
           {'parent',p5,'style','pushbutton','string','Extract','Callback',@Extract_Events} ...
           {'parent',p5,'style','pushbutton','tag','saveevents' ,'string','Save','Callback',@Save_Events} ...
           {'parent',p5,'style','popupmenu' ,'tag','evcodes','string',{'1'},'value',1} ...
           {'parent',p5,'style','pushbutton','tag','addevcode','string','+','Callback',@Add_EventCode} ...
           };
supergui(MainFig,geometry,geomvert,listui{:});
%            {'parent',p5,'style','pushbutton','tag','evright','string','==>','Callback',{@ShiftEvent,'right'}} ...

% panel 6 (misc)
geometry = {[1 1] [1 1]};
geomvert = [1 1];   
listui = { ...
         {'parent',p6,'style','checkbox','tag','diffchk','value',0,'string','diff','backgroundcolor','w','visible','on','Callback',{@UpdateManager}} ...
         {'parent',p6,'style','checkbox','tag','invertchk','value',0,'string','0)invert','backgroundcolor','w','Callback',{@UpdateManager}    } ...
         {'parent',p6,'style','checkbox','tag','phasechk','value',0,'string','phase','backgroundcolor','w','Callback',{@UpdateManager}} ...
         {'parent',p6,'style','checkbox','tag','envelope','value',0,'string','3)env','backgroundcolor','w','visible','on','Callback',{@UpdateManager}} ...
         };        
supergui(MainFig,geometry,geomvert,listui{:});

% panel 7
geometry = {[1 1] [1 1]};
geomvert = [1 1];   
listui = { ...
         {'parent',p7,'style','checkbox','tag','hilbfreqchk','value',0,'string','Hfreq','backgroundcolor','w','Callback',{@UpdateManager}} ...
         {'parent',p7,'style','checkbox',  'tag','imagesc','string','img','Callback',{@UpdateManager},'backgroundcolor','w'} ...           
         {'parent',p7,'style','checkbox','tag','topoplotchk','value',0,'string','topo','backgroundcolor','w','Callback',{@UpdateManager}} ...
         {'parent',p7,'style','edit','tag','nframes','string','1','Callback',{@UpdateManager}    } ...
         };        
supergui(MainFig,geometry,geomvert,listui{:});

% panel ? (lay/loc)?
% panel 9 (parms)
geometry = {[1 .3 .3 .3 .3]};
geomvert = [1];
listui = { ...
         {'parent',p9,'style','listbox'   ,'tag','chanlabels','value',[]  ,'string','','backgroundcolor','w','Max',350,'Min',0,'Callback',{@UpdateManager}} ...
         {'parent',p9,'style','checkbox','Visible','on','tag','avgbtn'  ,'string','avg','Callback',@UpdateManager} ...
         {'parent',p9,'style','checkbox'  ,'tag','stdchk'  ,'string','2)std','value',0,'Callback',{@UpdateManager},'backgroundcolor','w'} ...
         {'parent',p9,'style','checkbox'  ,'tag','histchk' ,'value',0,'string','hist','backgroundcolor','w','Callback',{@UpdateManager}} ...
         {'parent',p9,'style','checkbox'  ,'tag','xcorrchk','value',0,'string','xcorr','backgroundcolor','w','Callback',{@UpdateManager}} ... % ,'Visible','on'
         };
supergui(MainFig,geometry,geomvert,listui{:});
% {'parent',p9,'style','checkbox','tag','chanoverlaychk'  ,'string','overlay','value',0,'Callback',{@UpdateManager},'backgroundcolor','w'} ...

geometry = {[1 1 1] [1 1 1] [1 1 1] [1 1 1]};
geomvert = [1 1 1 1];
listui   = { ...
           {'parent',p10,'style','checkbox'  ,'tag','selgrad' ,'string','grad'   ,'Callback',{@ShiftChannel,'redraw'},'backgroundcolor','w' } ...
           {'parent',p10,'style','edit'      ,'tag','rejgrad' ,'string','3000'   ,'Callback',{@ShiftChannel,'redraw'},'backgroundcolor','w' } ...
           {'parent',p10,'style','pushbutton','tag','chanup'  ,'string','Chan UP','Callback',{@ShiftChannel,'up'}} ...
           {'parent',p10,'style','checkbox'  ,'tag','selmag'  ,'string','mag'    ,'Callback',{@ShiftChannel,'redraw'} ,'backgroundcolor','w' } ...
           {'parent',p10,'style','edit'      ,'tag','rejmag'  ,'string','6000'    ,'Callback',{@ShiftChannel,'redraw'},'backgroundcolor','w' } ...
           {'parent',p10,'style','edit'      ,'tag','chanstep','string','10'     ,'Callback',{@ShiftChannel,'redraw'}} ...  
           {'parent',p10,'style','checkbox'  ,'tag','seleeg'  ,'string','eeg'    ,'Callback',{@ShiftChannel,'redraw'} ,'backgroundcolor','w'} ...
           {'parent',p10,'style','edit'      ,'tag','rejeeg'  ,'string','200'    ,'Callback',{@ShiftChannel,'redraw'} ,'backgroundcolor','w' } ...
           {'parent',p10,'style','pushbutton','tag','chandown','string','DOWN'   ,'Callback',{@ShiftChannel,'down'}} ...  
           {'parent',p10,'style','checkbox'  ,'tag','selsti'  ,'string','sti'    ,'Callback',{@ShiftChannel,'redraw'} ,'backgroundcolor','w'} ...
           {'parent',p10,'style','edit'      ,'tag','rejsti'  ,'string','10'     ,'Callback',{@ShiftChannel,'redraw'} ,'backgroundcolor','w' } ...
           {'parent',p10,'style','checkbox'  ,'tag','lnfilter','string','line'   ,'Callback',{@UpdateManager} ,'backgroundcolor','w'} ...
           };
supergui(MainFig,geometry,geomvert,listui{:});

geometry = {[1 1] [1 1] [1 1]};
geomvert = [1 1 1];
listui = { ...
         {'parent',p8,'style','checkbox','tag','bpfilter','value',0,'string','1)bp','backgroundcolor','w','Callback',{@UpdateManager}} ...
         {'parent',p8,'style','edit','tag','bpfreq','string','[.1 30]','Callback',{@UpdateManager}    } ...
         {'parent',p8,'style','checkbox','tag','rejcheck','string','rej','value',0,'Callback',{@UpdateManager},'backgroundcolor','w'} ...
         {'parent',p8,'style','edit'    ,'tag','zlim','string','[]','Callback',{@UpdateManager}    } ...
         {'parent',p8,'style','checkbox','tag','clusterchk','value',0,'string','cluster','backgroundcolor','w','Callback',{@UpdateManager}} ...
         {'parent',p8,'style','edit','tag','clustercnt','string','3' } ...
         };        
supergui(MainFig,geometry,geomvert,listui{:});
%          {'parent',p8,'style','checkbox','tag','smoothchk','string','4)smooth','value',0,'Callback',{@UpdateManager},'backgroundcolor','w'} ...

% initialize structures
UserData.Raw.data = {}; % orig data (TimeSurfer structure)
UserData.Raw.map  = [];
UserData.Proc.map = []; % [ProcSamples(:) OrigSamples(:)]
UserData.Hilb.map = []; % Hilbert transform-based analytic signal
UserData.TFR.map  = [];
UserData.TFR.data = {};
UserData.FFT.map  = [];
UserData.time     = [];
UserData.frequencies = [];
UserData.sensor_info = [];
pos = get(0,'screensize');
set(MainFig,'position',[pos(3)/10 pos(4)/10 .8*pos(3) .8*pos(4)],'UserData',UserData);

set(MainFig,'Visible','on');

global options
options = [];
options.FigTag              = FigTag;
h = findobj('tag','chanlabels');
if length(h) > 1
  options.FigIndex          = cellfun(@(x)x==gcf,cellfun(@(x)get(x,'parent'),get(h,'parent'),'uniformoutput',false));
else
  options.FigIndex          = 1;
end
options.FigHandle           = MainFig;
options.procparms.blc       = 0;
options.procparms.blcwindow = [-inf inf];
options.procparms.bpfilter  = 0;
options.procparms.bpfreq    = [.1 30];
options.procparms.hilbert   = 0;
options.waveplots.xlim      = [];
options.waveplots.ylim      = [];
options.waveplots.labels    = {};
options.waveplots.datatypes = {}; % Raw, Proc, Hilb
options.TFRplots.xlim       = [];
options.TFRplots.ylim       = [];
options.TFRplots.zlim       = [];
options.TFRplots.band       = [];
options.TFRplots.labels     = {};
options.TFRplots.datatypes  = {}; % Amp, power, cmplx, Phase, Band, Zscore
options.events.label        = '';
options.events.time         = []; % events specified by channel
options.events.type         = [];
options.mode                = '';
options.waveplots.waitbar   = 0;
options.markers.time        = []; % events that are common to all channels
options.markers.type        = []; % event codes
options.markers.handle      = []; % handles to present vlines
options.markers.position    = []; % current position when marker is drawn
options.rejection.parms.prestim       = -.5;
options.rejection.parms.poststim      = .5;
options.rejection.parms.reject_mag    = 6000;
options.rejection.parms.reject_grad   = 3000;
options.rejection.parms.reject_eeg    = 200;
options.rejection.parms.reject_eog    = 200;
options.rejection.parms.reject_ieeg   = [];
options.rejection.parms.prescale_mag  = 10^15;  % to fT
options.rejection.parms.prescale_grad = 10^13;  % to fT/cm
options.rejection.parms.prescale_eeg  = 10^6;   % to uV
options.rejection.parms.prescale_eog  = 10^6;   % to uV
options.rejection.rejects.badtrials   = [];
options.rejection.rejects.badchans    = {};
options.rejection.rejects.time        = [];
options.rejection.update_flag         = 1;
options.subplots.handles              = [];
set(gcf,'CurrentObject',ax(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_Data(src,evnt,varargin)
global options
MainFig  = findobj('tag',options.FigTag);
UserData = get(MainFig,'UserData');
load_flag = 1;
if (nargin>2) && ischar(varargin{1})
  datafile = varargin{1};
elseif (nargin>2) && isstruct(varargin{1})
  type = gettype(varargin{1});
  UserData.(type).data{end+1} = varargin{1};
  load_flag = 0;
  filename  = 'tmp.mat';
elseif (nargin>2) && isnumeric(varargin{1})
  if nargin>3 && ischar(varargin{1}) && strcmp(varargin{1},'continuous')
    tmp  = ts_matrix2data(varargin{1},'continuous');
  else
    tmp  = ts_matrix2data(varargin{1});
  end
  type = gettype(tmp);
  UserData.(type).data{end+1} = tmp;
  load_flag = 0;
  filename  = 'tmp.mat';  
  clear tmp
else
  [filename,pathname] = uigetfile({'*.mat;*.fif;*.eeg;*.avg;*.cnt;*.edf;*.vhdr'},'Pick a file.','MultiSelect','on');% ;*.set
  if isequal(filename,0) || isequal(pathname,0), return; end
  if iscell(filename)
    datafile = cellfun(@(x)fullfile(pathname,x),filename,'uniformoutput',false);
    filename = filename{1};
  else
    datafile = [pathname filename];
    fprintf('Loading data: %s\n',datafile);
  end
  % re-initialize UserData
  UserData.Raw.data = {}; % orig data (TimeSurfer structure)
  UserData.Raw.map  = [];
  UserData.Proc.map = []; % [ProcSamples(:) OrigSamples(:)]
  UserData.Hilb.map = []; % Hilbert transform-based analytic signal
  UserData.TFR.map  = [];
  UserData.TFR.data = {};
  UserData.FFT.map  = [];
  UserData.time     = [];
  UserData.frequencies = [];
  UserData.sensor_info = [];  
end
options.flags(1) = 1;
[jnk1 jnk2 ftype] = fileparts(filename);
if ~strcmp(filename,'tmp.mat')
  UserData.Raw.data = {};
end
switch ftype
  case '.mat'
    if load_flag
      res = load(datafile);
      fld = fieldnames(res);
      for k = 1:length(fld)
        type = gettype(res.(fld{k}));
        if ~isempty(type)
          break
        elseif k==length(fld)
          error('Data type not recognized.');
        end
      end
      UserData.(type).data{end+1} = res.(fld{k});
    end
  case '.fif'
		type = 'Raw';
		tmp = ts_MNE_loadfif(datafile);
    UserData.(type).data{end+1} = tmp;
    clear tmp
  case '.eeg'
    type = 'Raw';
    tmp = ts_iEEG_eeg2epoch(datafile);
    UserData.(type).data{end+1} = tmp;
    clear tmp
  case '.avg'
    type = 'Raw';
    tmp = ts_iEEG_neuroscan2avg(datafile);
    UserData.(type).data{end+1} = tmp;
    clear tmp
  case '.cnt'
    type = 'Raw';
    tmp = ts_loadcnt(datafile);
    UserData.(type).data{end+1} = tmp;
    clear tmp
  case '.vhdr'
    type = 'Raw';
    tmp = ts_iEEG_eeg2epoch(datafile);
    UserData.(type).data{end+1} = tmp;
    clear tmp
  case '.edf'
    type = 'Raw';
    tmp = ts_loadedf(datafile);
    UserData.(type).data{end+1} = tmp;
    clear tmp
  case '.dat'
    type = 'Raw';
    tmp = ts_load_data(datafile);
    UserData.(type).data{end+1} = tmp;
    clear tmp    
end
    [jnk,field,param]           = ts_object_info(UserData.(type).data{end});
    UserData.(type).hdr         = rmfield(UserData.(type).data{end},field);
    UserData.(type).hdr.(field) = rmfield(UserData.(type).data{end}.(field),param);
    UserData.(type).time        = UserData.(type).data{end}.(field)(1).time;
    if isfield(UserData.(type).data{end}.(field),'frequencies')
      UserData.(type).frequencies = UserData.(type).data{end}.(field)(1).frequencies;
    end
    datastruct = UserData.(type).data{end}.(field);
    UserData.(type).data{end} = [];
    for k = 1:length(datastruct)
      if k == 1
        UserData.(type).data{end}   = datastruct(k).(param{1}); 
      else
        UserData.(type).data{end+1} = datastruct(k).(param{1}); 
      end
    end
    clear datastruct
    if ismember(field,{'epoch','epochs','cont','continuous'})
      if ndims(UserData.(type).data{end}) == 2
        options.mode = 'continuous';
      else
        options.mode = 'epoched';
      end
    elseif ismember(field,{'timefreq','sync','plv','coh','mscoh'})
      if ndims(UserData.(type).data{end}) == 3
        options.mode = 'continuous';
      else
        options.mode = 'epoched';
      end
    elseif ismember(field,{'average','averages'})
      UserData.(type).hdr.epochs  = UserData.(type).hdr.(field);
      UserData.(type).hdr         = rmfield(UserData.(type).hdr,field);
      field                       = 'epochs';
      options.mode                = 'continuous';
%       error('Averaged data is not yet supported.');
    else
      error('Data type not recognized.');
    end
% end
N                        = length(UserData.(type).time);
UserData.(type).map      = [1:N;1:N]';
if strcmp(options.mode,'continuous')
  channels               = 1:min(UserData.(type).hdr.num_sensors,10);
elseif strcmp(options.mode,'epoched')
  channels               = 1;
end
if isempty(UserData.sensor_info)
  UserData.sensor_info = UserData.(type).hdr.sensor_info;
end
if strcmp(type,'TFR')
  options.TFRplots.xlim   = [UserData.(type).time(1) min(UserData.(type).time(1)+3,UserData.(type).time(end))];
  options.TFRplots.ylim   = [UserData.(type).frequencies(1) UserData.(type).frequencies(end)];
  tmpdata = [UserData.(type).data{:}];
	options.TFRplots.zlim   = mean([tmpdata(:)]) + std([tmpdata(:)])*[-1 1];
  set(gethandle('zlim'),'string',num2str(options.TFRplots.zlim));
  clear tmpdata
  options.TFRplots.zlim(1)=0;
%   options.TFRplots.zlim   = [min([UserData.(type).data(:)]) max([UserData.(type).data(:)])];
  options.TFRplots.band   = [];
  options.TFRplots.labels = {UserData.(type).hdr.sensor_info(channels).label};
  options.TFRplots.datatypes = {param{1}};
  if isempty(get(gethandle('t0'),'string'))
    set(gethandle('t0'),'string',num2str(options.TFRplots.xlim(1)));
    set(gethandle('tf'),'string',num2str(options.TFRplots.xlim(end)));
    set(gethandle('tslide'),'value',0);
    UserData.time = UserData.(type).time;
  end
  if isempty(UserData.frequencies)
    UserData.frequencies = UserData.(type).frequencies;
  end
else
  options.waveplots.xlim   = [UserData.(type).time(1) min(UserData.(type).time(1)+3,UserData.(type).time(end))];
  tmpdata = cat(ndims(UserData.(type).data{end}),UserData.(type).data{:});
  options.waveplots.ylim   = [min([tmpdata(:)]) max([tmpdata(:)])];
  clear tmpdata
  options.waveplots.labels = {UserData.(type).hdr.sensor_info(channels).label};
  options.waveplots.datatypes = {'Raw'};
  set(gethandle('t0'),'string',num2str(options.waveplots.xlim(1)));
  set(gethandle('tf'),'string',num2str(options.waveplots.xlim(end)));  
  set(gethandle('tslide'),'value',0);
  if isempty(UserData.time), UserData.time = UserData.(type).time; end
end
set(MainFig,'UserData',UserData);
if isempty(get(gethandle('chanlabels'),'string'))
  set(gethandle('chanlabels'),'string',{UserData.(type).hdr.sensor_info.label});
  set(gethandle('chanlabels'),'value',channels);
end
UpdateManager(src,evnt);
disp('finished')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShiftTime(src,evnt,action)
global options
UserData = get(findobj('tag',options.FigTag),'UserData');
if isempty(UserData.Raw.data), return; end
ht0 = gethandle('t0');%findobj('tag','t0'); ht0 = ht0(options.FigIndex);
htf = gethandle('tf');%findobj('tag','tf'); htf = htf(options.FigIndex);
t0  = str2num(get(ht0,'string'));
tf  = str2num(get(htf,'string'));
period  = diff(options.waveplots.xlim);
offset  = -period;
switch action
  case 'right'
    offset = -offset;
  case 'slide'
    offset = 0;
    h      = gethandle('tslide');%findobj('tag','tslide'); h = h(options.FigIndex);
    t0     = min(UserData.time) + get(h,'value')*(UserData.time(end)-UserData.time(1));
    tf     = t0 + period;
end
if (t0+offset) < min(UserData.time)
  t0 = min(UserData.time)-offset;
  tf = t0 + period;
end
if t0 >= tf
  tf = t0 + period;
end
if (tf+offset) > max(UserData.time)
  t0 = max(UserData.time)-period-offset;
  tf = t0 + period;
end
set(ht0,'string',num2str(t0+offset));
set(htf,'string',num2str(tf+offset));
UpdateManager(src,evnt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShiftChannel(src,evnt,action)
global options
UserData = get(findobj('tag',options.FigTag),'UserData');
if isempty(UserData) || isempty(UserData.Raw.data), return; end
sens = UserData.sensor_info;
typestring = {sens.typestring};
mag_i   = find(strcmp ('mag' ,typestring));
grad_i  = find(strncmp('grad',typestring,4));
eeg_i   = find(strcmp ('eeg' ,typestring));
sti_i   = find(strcmp ('sti' ,typestring));
selchan = [];
if get(gethandle('selgrad'),'value'), selchan = [selchan grad_i]; end
if get(gethandle('selmag') ,'value'), selchan = [selchan mag_i];  end
if get(gethandle('seleeg') ,'value'), selchan = [selchan eeg_i];  end
if get(gethandle('selsti') ,'value'), selchan = [selchan sti_i];  end
if isempty(selchan), selchan = 1:length(sens); else selchan = sort(selchan); end
sens = sens(selchan);

step = str2num(get(gethandle('chanstep'),'string'));
list = gethandle('chanlabels');
curr = get(list,'string');
curr = curr(get(list,'value'));
ch0  = strmatch(curr{1}  ,{sens.label},'exact'); % ch0  = strmatch(curr{1}  ,{sens(selchan).label},'exact');
chf  = strmatch(curr{end},{sens.label},'exact'); % chf  = strmatch(curr{end},{sens(selchan).label},'exact');
n    = length(sens);                             % n    = length(sens(selchan));
if isempty(ch0), ch0 = 1;    end
if isempty(chf), chf = step; end
switch action
  case 'up'
    ch0 = max(ch0-step,1);
    chf = min(ch0+step-1,n);
  case 'down'
    chf = min(chf+step,n);
    ch0 = max(chf-step+1,1);
  case 'redraw'
    chf = min(ch0+step-1,n);
    ch0 = max(chf-step+1,1);
end
if length(ch0:chf) ~= step
  set(gethandle('chanstep'),'string',num2str(length(ch0:chf)));
end
set(gethandle('chanlabels'),'value',selchan([ch0:chf]));
UpdateManager(src,evnt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_Events(src,evnt)
[filename,pathname] = uigetfile({'*.mat;*.ev2;*.log;*.txt;*.dio;*.evt'},'Pick a file.','MultiSelect','on');
if isequal(filename,0) || isequal(pathname,0), return; end
if iscell(filename)
  datafile = cellfun(@(x)fullfile(pathname,x),filename,'uniformoutput',false);
  filename = filename{1};
else
  datafile = [pathname filename];
  fprintf('Loading events: %s\n',datafile);
end
global options
hMain         = findobj('tag',options.FigTag);
UserData      = get(hMain,'UserData');
[a,b,fext] = fileparts(datafile);
switch fext
  case '.mat'
    res = load(datafile);
    if isfield(res,'events')
      if numel(res.events) == 1
        options.markers.time = res.events.time;
        options.markers.type = res.events.type;
      else
        options.events = res.events;
      end
    elseif isfield(res,'peaks')
      options.events  = res.peaks;
    end
  case '.ev2'
    [a,types,c,d,e,markers] = textread(datafile,'%f %f %f %f %f %f');
    UserData = get(findobj('tag',options.FigTag),'UserData');
    options.markers.time = markers/UserData.Raw.hdr.sfreq;
    options.markers.type = types;
  case '.txt'
    [content,res] = mmil_readtext(datafile,'\t');
    evnt = content(2:end,:);
    options.markers.time = [evnt{:,2}]/UserData.Raw.hdr.sfreq;
    options.markers.type = [evnt{:,3}];
  case '.evt' % BESA
    [content,res] = mmil_readtext(datafile,'\t');
    % BESA event file
    Tunits = content{1,1};
    if strcmp(Tunits,'Tms')
      c  = 1E3;
    elseif strcmp(Tunits,'Tmu')
      c  = 1E6;
    else
      c  = 1;
    end
    evnt = content(2:end,:);
    options.markers.time = ([evnt{:,1}]')/c;
    options.markers.type = [evnt{:,2}]';
    clear Tunits c
end
if iscell(options.events(1).label)
  for k = 1:length(options.events)
    options.events(k).label = options.events(k).label{1};
  end
end
h         = gethandle('evcodes');
lastcodes = get(h,'string');
lastcodes = cellfun(@str2num,lastcodes);
currcodes = unique([lastcodes' options.markers.type]);
set(h,'string',cellfun(@num2str,num2cell(currcodes),'uniformoutput',false));
% turn on event vlines
set(gethandle('eventchk'),'value',1);
UpdateManager(src,evnt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Extract_Events(src,evnt)
global options
labels        = options.waveplots.labels; 
hMain         = findobj('tag',options.FigTag);
UserData      = get(hMain,'UserData');
if isempty(UserData) || isempty(UserData.Raw.data), return; end
[sel1,chans]  = match_str(labels,{UserData.Raw.hdr.sensor_info.label});

if (length(labels) < 1) || (length(labels) > 12)
  error('Select 1-12 channels before extracting events.');
end
% parse trigger channels
allsamp = [];
allcond = [];
for k = 1:length(UserData.Raw.data)
  x   = UserData.Raw.data{k}(chans,:);
  t   = UserData.time;
  Fs  = UserData.Raw.hdr.sfreq;
  trigdata    = ts_matrix2data(x,'time',t,'sfreq',Fs,'continuous',1);
  [samp,cond] = extract_triggers(trigdata);
  allsamp = [allsamp samp];
  allcond = [allcond cond];
end
options.markers.time = t(allsamp);
options.markers.type = allcond;
% add to event code list
h         = gethandle('evcodes');
lastcodes = get(h,'string');
lastcodes = cellfun(@str2num,lastcodes);
currcodes = unique([lastcodes' allcond]);
set(h,'string',cellfun(@num2str,num2cell(currcodes),'uniformoutput',false));
% turn on event vlines
set(gethandle('eventchk'),'value',1);
UpdateManager(src,evnt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Add_EventCode(src,evnt)
h   = gethandle('evcodes');
str = get(h,'string');
num = str2num(str{end});
str{end+1} = num2str(num+1);
set(h,'string',str);
set(h,'value',length(str));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Save_Events(src,evnt)
global options
% ... save marked events to text file
[filename,pathname] = uiputfile({'*.log;'},'Save as','marked_events.log');
if isequal(filename,0) || isequal(pathname,0)
  return;
end
outfile = fullfile(pathname,filename);
[fpath,fname,fext] = fileparts(outfile);

UserData        = get(findobj('tag',options.FigTag),'UserData');
events.latency  = round(options.markers.time * UserData.Raw.hdr.sfreq);
events.time     = options.markers.time;
events.type     = options.markers.type;
save(strrep(outfile,fext,'.mat'),'events');

% write reject info to text file
fid = fopen(outfile,'w+');
fprintf(fid,'Marked events (%s)\n\n',datestr(now));
fprintf(fid,'   Latency       Time       Type\n');
for k = 1:length(events.time)
  fprintf(fid,'%10g %10g %3g\n',events.latency(k),events.time(k),events.type(k));
end
fclose(fid);
fprintf('Saving log files:\n%s\n%s\n',outfile,strrep(outfile,fext,'.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShiftEvent(src,evnt)
global options
UserData = get(findobj('tag',options.FigTag),'UserData');
if isempty(UserData.Raw.data) || isempty(options.events.label), return; end
ht0 = gethandle('t0');%findobj('tag','t0'); ht0 = ht0(options.FigIndex);
htf = gethandle('tf');%findobj('tag','tf'); htf = htf(options.FigIndex);
period  = diff(options.waveplots.xlim);
[sel1,sel2] = match_str({options.events.label},options.waveplots.labels);
if isempty(sel1), return; end
evtimes = [options.events.time];
switch action
  case 'right'
    evtime = evtimes(find(evtimes>options.waveplots.xlim(2),1,'first'));
  case 'left'
    evtime = evtimes(find(evtimes<options.waveplots.xlim(1),1,'first'));
end
if isempty(evtime), return; end
t0 = evtime - period/2;
tf = t0 + period;
set(ht0,'string',num2str(t0));
set(htf,'string',num2str(tf));
UpdateManager(src,evnt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateManager(src,evnt)
UpdateOptions(src);
global options
types = options.waveplots.datatypes;
if ~isempty(options.TFRplots.datatypes)
  types = {types{:} 'TFR'};
end
for k = 1:length(types)
  UpdateData(types{k});
end
UpdateRejection % get rejection info for the next plot
UpdatePlots;    % update plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowValue(src,evnt)
% global options
% UserData = get(findobj('tag',options.FigTag),'UserData');
% if isempty(UserData) || isempty(UserData.Raw.data), return; end
% % figpos  = get(gcf,'position');
% % figh    = figpos(4);
% nplots  = length(options.subplots.handles);
% % ploth   = figh / nplots;
% % plotnum = find(ismember(options.subplots.handles,gca));
% % if isempty(plotnum), return; end
% % ypos    = plotnum*(ploth) + ploth/2;
% % pt=get(gca,'currentpoint');
% 
% figpt = get(gcf,'currentpoint');
% figptx = figpt(1);
% figpty = figpt(2);
% fight = get(gcf,'position');
% fight = fight(4);
% gca = options.subplots.handles(round((figpty/fight)*nplots));
% 
% % get height of figure
% % divide by # of subplots from options.handles.subplots
% % determine which subplot the mouse is over (use ypos,#subplots,figheight)
% % get sensor label from the subplots ylabel string
% % get time and amplitude from the subplots XData and YData
% % display sensor label, time, and amplitude in the figure name
% 
% pt  = get(gca,'currentpoint');
% sen = get(gca,'UserData');%get(get(gca,'ylabel'),'string');
% str = sprintf('%s: %g at %gs',sen,pt(1,2),pt(1,1));
% set(src,'Name',str);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateOptions(src)
global options
ht0 = gethandle('t0');%findobj('tag','t0'); ht0 = ht0(options.FigIndex);
htf = gethandle('tf');%findobj('tag','tf'); htf = htf(options.FigIndex);
t0  = str2double(get(ht0,'string'));
tf  = str2double(get(htf,'string'));
h   = gethandle('chanlabels');%findobj('tag','chanlabels'); h = h(options.FigIndex);
str = get(h,'string');
val = get(h,'value' );
options.waveplots.xlim   = [t0 tf];
options.waveplots.labels = str(val);
options.TFRplots.xlim    = [t0 tf];
options.TFRplots.labels  = str(val);
UserData = get(findobj('tag',options.FigTag),'UserData');
if isempty(UserData.time), return; end
trange = [UserData.time(1) UserData.time(end)];
tmin   = min(UserData.time);
newval = (t0-tmin)/(trange(2)-trange(1));
if newval < 0, newval = 0; end
if newval > 1, newval = 1; end
h = gethandle('tslide');%findobj('tag','tslide'); h = h(options.FigIndex);
% Determine whether to update rejection info or not (do not if only displayed
% channels has changed)
tags     = {'chanlabels','chanstep','chanup','chandown','selgrad','selmag',...
            'seleeg','selsti','autoscale','eventchk','invertchk','phasechk',...
            'diffchk','hilbfreqchk','imagesc','topoplotchk','nframes',...
            'avgbtn','histchk','xcorrchk','zlim','clusterchk','clustercnt'};
options.rejection.update_flag = 1;
for k = 1:length(tags)
  if isequal(src,findobj('tag',tags{k}))
    options.rejection.update_flag = 0;
    break;
  end
end
set(h,'value',newval);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function options = UpdateRejection(data,options)
function UpdateRejection
% Note: this function is always called before updating the plots panel
% Purpose: use thresholding to reject trials in selected non-rejected chans
% (given markers, rejcheck) -> calc trial rejects & add to rejects list
% (else) -> do nothing

% TODO: increase efficiency; this function reprocesses data every time the
% plots are updated to recalculate rejected trials. modify to avoid recalc.

global options
if ~options.rejection.update_flag, return; end
if isempty(options.markers.time) || ~get(gethandle('rejcheck'),'value'), return; end
UserData   = get(findobj('tag',options.FigTag),'UserData');
type       = 'Raw';
if isempty(UserData.(type).data)   , return; end
parms      = options.rejection.parms;
sens       = UserData.(type).hdr.sensor_info;
[rej,jk]   = match_str({sens.label},options.rejection.rejects.badchans);
chans      = setdiff(1:length(sens),rej);
sens       = sens(chans);
typestring = {sens.typestring};
mag_i   = strcmp ('mag' ,typestring);
grad_i  = strncmp('grad',typestring,4);
eeg_i   = strcmp ('eeg' ,typestring);
eog_i   = strcmp ('eog' ,typestring);
% generate rejection thresholds for each channel
parms.reject_grad = str2num(get(gethandle('rejgrad'),'string'));
parms.reject_mag  = str2num(get(gethandle('rejmag') ,'string'));
parms.reject_eeg  = str2num(get(gethandle('rejeeg') ,'string'));
parms.reject_sti  = str2num(get(gethandle('rejsti') ,'string'));
parms.reject_thresh         = zeros(length(chans),1);
parms.reject_thresh(mag_i)  = parms.reject_mag/parms.prescale_mag;
parms.reject_thresh(grad_i) = parms.reject_grad/parms.prescale_grad;
parms.reject_thresh(eeg_i)  = parms.reject_eeg/parms.prescale_eeg;
parms.reject_thresh(eog_i)  = parms.reject_eog/parms.prescale_eog;
parms.reject_thresh(parms.reject_thresh<=0) = inf;
options.rejection.parms = parms;
xlim = options.waveplots.xlim;
tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2)); 
data = UserData.(type).data{1}(chans,tix); % (chans,tix)
if get(gethandle('bpfilter'),'value')
  data = blc(data);
  data = ts_freq_filt(data',UserData.(type).hdr.sfreq,str2num(get(gethandle('bpfreq'),'string')),[0 0],'bandpass')';
end
if get(gethandle('lnfilter'),'value')
  data  = blc(data);
  data  = ts_freq_filt(data',UserData.(type).hdr.sfreq,60,5,'notch')';
end          
if get(gethandle('stdchk'),'value')
  StdDevWindow = min(50,floor(length(data)/10));
  StdDevWindow = StdDevWindow + ~mod(StdDevWindow,2);
  npoints      = length(data) - StdDevWindow + 1;
  offset       = floor(StdDevWindow/2);
  tmp          = zeros(size(data));
  for l = 1:npoints
    tmp(l+offset) = std(data(l:l+StdDevWindow-1));
  end
  data = tmp;
end          
if get(gethandle('envelope'),'value')
  data = abs(hilbert(data));
end

poststim = parms.poststim;
prestim  = parms.prestim;
evnts = options.markers.time;
% evnts = evnts(evnts>=xlim(1)&evnts<=xlim(2));
t     = UserData.time(tix);
level = options.rejection.parms.reject_thresh * ones(1,length(t)); % chan x time
tmp   = abs(data) > level;
tmp   = sum(tmp,1);
ind   = tmp > 0;
tover = t(ind);
rejix = find(cellfun(@(x)any(find((tover>(x+prestim))&(tover<(x+poststim)))),num2cell(evnts)));
options.rejection.rejects.badtrials = rejix;
options.rejection.rejects.time      = evnts(rejix);
% UpdateManager(src,evnt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateFreq
  % UpdateFcontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateChan
  % UpdateLay
  % UpdateLoc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateTcontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateFcontrol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateLay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateLoc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateData(type)
global options
% get values from uicontrols
h           = gethandle('t0');%findobj('tag','t0'); h = h(options.FigIndex);
t0          = str2num(get(h,'string'));
h           = gethandle('tf');%findobj('tag','tf'); h = h(options.FigIndex);
tf          = str2num(get(h,'string'));
% user data
hMain       = findobj('tag',options.FigTag);
UserData    = get(hMain,'UserData');
T           = UserData.time; % complete, original time vector
map         = UserData.(type).map;
ActiveIndex = [nearest(T,t0):nearest(T,tf)]';
ProcIndex   = setdiff(ActiveIndex,map(:,2));
% is it necessary to process anything?
if ~isempty(ProcIndex) % yes
  for k = 1:length(UserData.Raw.data)    
    ProcData  = UserData.Raw.data{k}(:,ProcIndex,:,:);
    ProcData  = eval(sprintf('Calc%s(ProcData,UserData.Raw.time(ProcIndex),UserData.Raw.sfreq);',type));
    % where to insert this new data:
    k = find(map(:,2) > ProcIndex(1),1,'first');
    % add new results to UserData
    tmp = UserData.(type).data{k};
    UserData.(type).data{k} = cat(2,tmp(:,1:k-1,:,:),ProcData,tmp(:,k:end));
    if k == 1
      tmp = UserData.(type).map;
      UserData.(type).map  = [ones(size(tmp,1)+numel(ProcIndex),1) cat(1,tmp(1:k-1,2),ProcIndex,tmp(k:end,2))];
      tmp = UserData.(type).time;
      UserData.(type).time = [tmp(1:k-1) T(ProcIndex) tmp(k:end)];
    end
    clear tmp
    if k == length(UserData.Raw.data)
      set(hMain,'UserData',UserData);
    end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outdat = CalcProc(seldat,t,Fs)
global options
seldat = ts_matrix2epoch(seldat,'time',t,'sfreq',Fs);
h1     = gethandle('bpfilter');%findobj('tag','bpfilter'); h1 = h1(options.FigIndex);
h2     = gethandle('bpfreq');%findobj('tag','bpfreq');   h2 = h2(options.FigIndex);
outdat = ts_preproc(seldat,'bpfilter',get(h1,'string'),'bpfreq',str2num(get(h2,'string')));
outdat = outdat.epochs.data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outdat = CalcHilbert(seldat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outdat = CalcTFR(seldat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outdat = CalcFFT(seldat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdatePlots
global options
if options.waveplots.xlim(1) >= options.waveplots.xlim(2), return; end
if isempty(options.waveplots.labels) && isempty(options.TFRplots.labels), return; end
colors   = 'kbrgmy';
lntype   = {'-',':','-.','--'};
hMain    = findobj('tag',options.FigTag);
UserData = get(hMain,'UserData');
wavtypes = options.waveplots.datatypes;
wavchans = options.waveplots.labels;
tfrtypes = options.TFRplots.datatypes;
tfrchans = options.TFRplots.labels;
sens     = {wavchans{:} tfrchans{:}};
uniqsens = unique(sens);
uniqcnt  = 0;
% [chans jnk] = match_str({UserData.Raw.hdr.sensor_info.label},sens);
labels = {UserData.sensor_info.label};
children = get(gethandle('plots'),'Children');
delete(children);
tfrlims  = str2num(get(gethandle('zlim'),'string'));
hwavplots= [];
% plot data
if strcmp(options.mode,'continuous')
  % plot continuous data
  nrow   = ~isempty(wavtypes)*length(wavchans) + ~isempty(tfrtypes)*length(tfrchans);
  ncol   = 1;
  cnt    = 0;
  maxy   = [];
  miny   = [];
  xlim   = [];
  if get(gethandle('imagesc'),'value')
    h = subplot(2,1,1,'parent',gethandle('plots'));
    set(h,'ButtonDownFcn',@onclick);
    type      = wavtypes{1};
    sens      = setdiff(sens,options.rejection.rejects.badchans);
    if isempty(sens), return; end
    [sel,jnk] = match_str({UserData.(type).hdr.sensor_info.label},sens);
    sel  = unique(sel);
    xlim = options.waveplots.xlim;
    tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2)); 
    j    = 1;
    dat  = UserData.(type).data{j}(sel,tix);
    if get(gethandle('autoscale'),'value')
      for ch = 1:length(sel)
        dat(ch,:) = dat(ch,:) / max(abs(dat(ch,:)));
      end
    end
    if get(gethandle('bpfilter'),'value')
      dat  = blc(dat);
      dat  = ts_freq_filt(dat',UserData.(type).hdr.sfreq,str2num(get(gethandle('bpfreq'),'string')),[0 0],'bandpass')';
    end
    if get(gethandle('lnfilter'),'value')
      dat  = blc(dat);
      dat  = ts_freq_filt(dat',UserData.(type).hdr.sfreq,60,5,'notch')';
    end
    h = imagesc(UserData.(type).time(tix),1:length(sel),dat);
    set(h,'ButtonDownFcn',@onclick);
    axis tight;
%     for ch = 1:size(alldata{1},1)
%       dat(ch,:) = alldata{j}(ch,:);%allchannels(ch).allsets{j};
%     end
    h = subplot(2,1,2,'parent',gethandle('plots'));
    set(h,'ButtonDownFcn',@onclick);
    plot(UserData.(type).time(tix),dat);
    set(gca,'ButtonDownFcn',@onclick);
    axis tight   
    maxvals  = max(dat,[],2);
    maxval   = max(maxvals);
    maxchan  = find(maxvals == maxval,1,'first');
    maxlabel = UserData.(type).hdr.sensor_info(sel(maxchan)).label;
    maxtind  = find(maxval==dat(maxchan,:),1,'first');
    maxtime  = UserData.(type).time(tix(maxtind));
%     htxt = text(maxtime,maxval,maxlabel);
    htxt = text(double(maxtime),double(maxval),maxlabel);
    set(gca,'UserData',maxlabel);
    set(htxt,'fontsize',8,'ButtonDownFcn',@onclick);    
%     options.rejection.rejects.badchans = {options.rejection.rejects.badchans{:} this};
%     UpdateManager;
    clear dat 
    options = drawmarkers(options,xlim);  
    return
  end
  % loop over channels
  for n  = 1:length(sens)
%     if options.waveplots.waitbar
%       if n==1
%         waitbar(0,'Updating plots...');
%       end
%       waitbar(n/length(sens));
%     end    
    this = sens{n};
    if ismember(this,options.rejection.rejects.badchans)
      badchan_flag = 1;
    else
      badchan_flag = 0;
    end
    if n > 1 && ismember(this,sens(1:n-1)), continue; end
    uniqcnt       = uniqcnt + 1;
    plot_wave     = ~isempty(wavtypes) && ismember(this,wavchans);
    plot_timefreq = ~isempty(tfrtypes) && ismember(this,tfrchans);  
%     chan = strmatch(this,labels,'exact');
    if plot_wave
      % plot t-domain signals
      cnt = cnt + 1;
      allplots(cnt) = subplot(nrow,ncol,cnt,'parent',gethandle('plots'));
      hwavplots(end+1) = allplots(cnt);
      pos = get(gca,'position');
      set(gca,'UserData',this,'position',[.05 pos(2) .9 pos(4)],'clipping','off');      
      if badchan_flag
        hylabel = ylabel(this);
        set(hylabel,'ButtonDownFcn',@onclick);
        set(gca,'Xtick',[],'Ytick',[]);
        continue;
      end
      % overlay different wave types
      cntoverlay = 0;
      for k  = 1:length(wavtypes)
        type = wavtypes{k};
        chan = strmatch(this,{UserData.(type).hdr.sensor_info.label},'exact');
        xlim = options.waveplots.xlim;
        tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2)); 
        nset = numel(UserData.(type).data);
        ymintype = [];
        ymaxtype = [];
        for j = 1:nset
          cntoverlay = cntoverlay + 1;
          dat  = UserData.(type).data{j}(chan,tix);
          if isempty(dat), continue; end
          if get(gethandle('invertchk'),'value')
            dat = -dat;
          end
          if get(gethandle('bpfilter'),'value')
            dat = blc(dat);
            dat = ts_freq_filt(dat',UserData.(type).hdr.sfreq,str2num(get(gethandle('bpfreq'),'string')),[0 0],'bandpass')';
          end
          if get(gethandle('lnfilter'),'value')
            dat  = blc(dat);
            dat  = ts_freq_filt(dat',UserData.(type).hdr.sfreq,60,5,'notch')';
          end          
          if get(gethandle('stdchk'),'value')
            StdDevWindow = min(50,floor(length(dat)/10));
            StdDevWindow = StdDevWindow + ~mod(StdDevWindow,2);
            npoints      = length(dat) - StdDevWindow + 1;
            offset       = floor(StdDevWindow/2);
            tmp          = zeros(size(dat));
            for l = 1:npoints
              tmp(l+offset) = std(dat(l:l+StdDevWindow-1));
            end
            dat = tmp;
          end          
          if get(gethandle('envelope'),'value')
            dat = abs(hilbert(dat));
          end
%           if get(gethandle('smoothchk'),'value')
%             tmp = str2num(get(gethandle('zlim'),'string'));
%             if ~isempty(tmp) && numel(tmp)==1 && isint(tmp)
%               npoint = tmp;
%             else
%               npoint = 10;
%             end
%             dat = smooth(dat,npoint,'lowess');
%           end
          ylim = [min(dat) max(dat)];
          if isempty(ylim) || ylim(1)==ylim(2), continue; end
          if isempty(miny)
            miny = ylim(1);
            maxy = ylim(2);
          else
            miny = min(miny,ylim(1));
            maxy = max(maxy,ylim(2));
          end
          if isempty(ymintype)
            ymintype = ylim(1);
            ymaxtype = ylim(2);
          else
            ymintype = min(ymintype,ylim(1));
            ymaxtype = max(ymaxtype,ylim(2));
          end
          hold on
          if get(gethandle('histchk'),'value')
            try
              hist(dat,round(min(100,max(20,length(dat)/100))));
            end
            axis tight
            xlim = get(gca,'xlim');
            ylim = get(gca,'ylim');            
            miny = ylim(1); ymintype = miny;
            maxy = ylim(2); ymaxtype = maxy;            
            if get(gethandle('clusterchk'),'value')
              nclusters = str2num(get(gethandle('clustercnt'),'string'));
              try
                [cidx,ctrs]  = kmeans(dat,nclusters);
              catch
                fprintf('Failed to find %g clusters in %s; change number and try again\n',nclusters,this);
                continue
              end
              [jnk,ctrord] = sort(ctrs);
              for cnum = 1:nclusters
                vline(ctrs(cnum),'Color',colors(ctrord==cnum),'LineWidth',1);
                if ctrs(cnum) < max(ctrs)
                  vline(max(dat(cidx==cnum)),'Color','k','LineWidth',5);
                end
              end
            end  
          elseif get(gethandle('clusterchk'),'value')
            nclusters = str2num(get(gethandle('clustercnt'),'string'));
            if isempty(nclusters), nclusters = 3; end
            t = UserData.(type).time(tix)';
            try
              [cidx,ctrs]  = kmeans(dat,nclusters);
            catch
              fprintf('Failed to find %g clusters in %s; change number and try again\n',nclusters,this);
              continue
            end              
            [jnk,ctrord] = sort(ctrs);
            for cnum = 1:nclusters
              plot(t(cidx==cnum),dat(cidx==cnum),[colors(ctrord==cnum) '.'])
              hline(ctrs(cnum),'Color','k','LineWidth',1)
            end
            clear t
          else
            lnstyle = [colors(mod(cntoverlay,length(colors))) lntype{mod(k,length(colors))}];
            if k==length(wavtypes)
              h=plot(UserData.(type).time(tix),dat,lnstyle);
              set(h,'ButtonDownFcn',@onclick);
            else
              plot(UserData.(type).time(tix),dat,lnstyle);
            end
          end
          hold on
          if get(gethandle('diffchk'),'value')
            tmp = gradient(dat);
            lnstyle = [colors(mod(cntoverlay,length(colors))) lntype{mod(k,length(colors))}];
            plot(UserData.(type).time(tix),.5*ymintype*(tmp/max(tmp(:))),lnstyle);
            set(gca,'ButtonDownFcn',@onclick);
          end
          if get(gethandle('phasechk'),'value')
            tmp = angle(hilbert(dat));
            lnstyle = [colors(mod(cntoverlay,length(colors))) lntype{mod(k,length(colors))}];
            plot(UserData.(type).time(tix),.5*ymintype*(tmp/max(tmp(:))),lnstyle);
            set(gca,'ButtonDownFcn',@onclick);
            clear tmp
          end
          if get(gethandle('hilbfreqchk'),'value')
            tmp = angle(hilbert(dat))';
            tmp = (1/(2*pi)) * diff(unwrap(tmp))./diff(UserData.(type).time(tix));
            lnstyle = [colors(mod(cntoverlay,length(colors))) lntype{mod(k,length(colors))}];
            plot(UserData.(type).time(tix(1:end-1)),.5*ymintype*(tmp/max(tmp(:))),lnstyle);
            set(gca,'ButtonDownFcn',@onclick);
            clear tmp
          end
          if get(gethandle('avgbtn'),'value') || get(gethandle('rejcheck'),'value') || get(gethandle('xcorrchk'),'value') || get(gethandle('topoplotchk'),'value') %|| get(gethandle('chanoverlaychk'),'value') || get(gethandle('imagesc'),'value')
            % collect data for calculation later
            alldata{j}(n,:) = dat;
%             allchannels(n).allsets{j} = dat;
          end
        end
      end
      set(gca,'xlim',xlim);
      if ~isempty(ymintype) && ~isempty(ymaxtype)
        set(gca,'ylim',[ymintype ymaxtype]);
      end
      set(gca,'ButtonDownFcn',@onclick);
      % mark events
      if get(gethandle('eventchk'),'value') && ~isempty(options.events(1).label) && ismember(this,{options.events.label})
        chix = strmatch(this,{options.events.label},'exact');
        tix  = options.events(chix).time>xlim(1) & options.events(chix).time<xlim(2);
        evt  = options.events(chix).time(tix);
        typ  = unique(options.events(chix).type(tix));
        for k = 1:length(typ)
          tmp = evt([options.events(chix).type(tix)]==typ(k));
          for j = 1:length(tmp)
            vline(tmp(j),'Color',colors(k),'ButtonDownFcn',@onclick);
          end
        end
        clear tmp
      end
      hylabel = ylabel(this);
      set(hylabel,'ButtonDownFcn',@onclick);
      hold off
    end
    set(gca,'ytick',[]);
    if plot_timefreq
%       if (ylim(1)<=0 && ylim(2)>=0), hline(0,'k'); end
      axis off
    elseif uniqcnt<length(uniqsens)
      set(gca,'xtick',[]);
    end
    box off
    grid on
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % maybe remove this later
        if length(uniqsens)>20 && uniqcnt<length(uniqsens), axis off; end
%         hline(0,'k')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % plot tf-domain images
    if plot_timefreq
      if uniqcnt==length(uniqsens)
        % turn off x-axis on wave plot since this will be below it
        set(gca,'xtick',[]);
      end
      cnt = cnt + 1;
      allplots(cnt) = subplot(nrow,ncol,cnt,'parent',gethandle('plots'));
      set(gca,'UserData','TFR');
      pos = get(gca,'position');
      set(gca,'UserData',this,'position',[.05 pos(2) .9 pos(4)]);         
      % for each tfrtype
        type = 'TFR';
        chan = strmatch(this,{UserData.(type).hdr.sensor_info.label},'exact');
        if isempty(chan), continue; end
        xlim = options.TFRplots.xlim;
        ylim = options.TFRplots.ylim;
        if isempty(tfrlims)
          zlim = options.TFRplots.zlim;
        else
          zlim = tfrlims;
        end
        tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2));      
        fix  = nearest(UserData.(type).frequencies,ylim(1)):nearest(UserData.(type).frequencies,ylim(2));
        if length(UserData.(type).data) > 1
          warning('Only showing the first TFR data set.  Cannot overlay TFRs.');
        end
        dat  = UserData.(type).data{1}(chan,tix,fix,:);
        if ~isreal(dat),   dat = abs(double(dat)).^2; end
        if ndims(dat) > 3, dat = nanmean(dat,4);      end
        dat  = squeeze(dat)';
        himage = imagesc(UserData.(type).time(tix),UserData.(type).frequencies(fix),dat,zlim); 
        set(gca,'YDir','normal')
        set(himage,'ButtonDownFcn',@onclick);
        foi = UserData.(type).frequencies(fix); nn = length(get(gca,'ytick'));
        dy=floor(length(foi)/nn); yticklabels = foi(dy:dy:length(foi));
        nn=length(yticklabels); dy=(max(foi)-min(foi))/nn; yticks=min(foi)+dy:dy:max(foi); 
        set(gca,'ytick',yticks,'yticklabel',yticklabels);
        set(gca,'xlim',xlim);
        set(gca,'ylim',ylim);
        set(gca,'clim',zlim);      
        hold off
    end
    if uniqcnt<length(uniqsens)
      set(gca,'xtick',[]);
      set(gca,'ytick',[]);
      hylabel = ylabel(this);
      set(hylabel,'ButtonDownFcn',@onclick);
    end
    box off
%     grid on
  end % end loop over channels
  % scale all wave plots to max/min over all channels being displayed
  if ~isempty(miny) && ~get(gethandle('autoscale'),'value')
    set(hwavplots,'ylim',[miny maxy]);
  end
%   set(hwavplots,'clipping','off');
%   % adjust positions to use as much space as possible;
%   % start at bottom and move up to avoid overlap deletion
% %   if n > 1
% %     for k = length(allplots):-1:1
% %       axes(allplots(k));
% %       pos = get(gca,'position');
% %       set(gca,'position',[.05 .95*pos(2) .9 1.05*pos(4)]);
% %     end
% %   end
%   if get(gethandle('rejcheck'),'value')
%     options = UpdateRejection(alldata,options);
%   end
  if isempty(xlim), return; end
  options = drawmarkers(options,xlim);
  if get(gethandle('xcorrchk'),'value')
    % calculate cross-correlations and plot them in a new figure
    for j = 1:length(alldata)
      % open a new figure for each set of cross-correlations
      hfig  = findobj('tag','xcorrFig');
      if isempty(hfig)
        figure('Name',['Normalized cross-correlations: set ' num2str(j)],'tag','xcorrFig');
      else
        figure(hfig);
      end
      nchan = size(alldata{1},1);
      for n = 1:nchan
        for m = 1:nchan
          if m < n, continue; end
          subplot(nchan,nchan,(n-1)*nchan+m)
          y1 = alldata{j}(n,:);%allchannels(n).allsets{j};
          y2 = alldata{j}(m,:);%allchannels(m).allsets{j};
          [xc,lags] = xcorr(y1,y2,'coeff');
          plot(lags/UserData.Raw.hdr.sfreq,xc); axis tight; vline(lags(ceil(length(lags)/2))/UserData.Raw.hdr.sfreq,'k'); hline(0,'k');
          if n==1, title(sens{m});  end
          if m==n, hylabel=ylabel(sens{m}); set(hylabel,'ButtonDownFcn',@onclick); end
          clear y1 y2
        end
      end
    end
  end
%   if get(gethandle('chanoverlaychk'),'value')
%     try
%       j    = 1;
%       tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2)); 
%       hfig = findobj('tag','OverlayFig');
%       if isempty(hfig)
%         figure('tag','OverlayFig');
%       else
%         figure(hfig);
%       end
%       for ch = 1:size(alldata{1},1)
%         tmp(ch,:) = alldata{j}(ch,:);%allchannels(ch).allsets{j};
%       end
%       plot(UserData.(type).time(tix),tmp);
%       axis tight
%       clear tmp
%     catch
%       close(gcf)
%     end
%     options = drawmarkers(options,xlim);
%   end    
  if get(gethandle('imagesc'),'value')    
    j    = 1;
    tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2)); 
    hfig = findobj('tag','imagescFig');
    if isempty(hfig)
      figure('tag','imagescFig');
    else
      figure(hfig);
    end
    for ch = 1:size(alldata{j},1)
      tmp(ch,:) = alldata{j}(ch,:);
    end    
    imagesc(UserData.(type).time(tix),1:size(alldata{j},1),tmp);
    axis tight;
    clear tmp    
    options = drawmarkers(options,xlim);
  end
  if get(gethandle('topoplotchk'),'value')
    % prepare topoplot plotting parameters
    cfg = [];
    cfg.emarker     = 'o';
    cfg.ecolor      = [0 0 0];
    cfg.emarkersize = 2;
    cfg.hlmarker    = 'o';
    cfg.hlcolor     = [0 0 0];
    cfg.hlmarkersize= 6;
    cfg.hllinewidth = 3;
    cfg.hcolor      = [0 0 0];
    cfg.hlinewidth  = 2;
    cfg.fontsize    = 6;
    cfg.efsize      = get(0,'DefaultAxesFontSize');    
    % prepare position Info
    nchan = UserData.(type).hdr.num_sensors;
    T     = UserData.(type).hdr.coor_trans.device2head;
    pos   = zeros(nchan,3); % (x,y,z) for each channel
    for k = 1:nchan
      loc         = UserData.(type).hdr.sensor_info(k).loc;
      if any(strmatch('grad',UserData.(type).hdr.sensor_info(k).typestring)) || ...
         any(strmatch('mag' ,UserData.(type).hdr.sensor_info(k).typestring))
        loc       = T*loc;
      end
      pos(k,1:3)  = loc(1:3,4);
    end
    method = 'stereographic'; % gnomic, stereographic, ortographic, inverse, polar
    prj    = elproj(pos, method); % * [0 1; -1 0];
              % ELPROJ makes a azimuthal projection of a 3D electrode cloud
              %  on a plane tangent to the sphere fitted through the electrodes
              %  the projection is along the z-axis
    X = prj(:,1);   % x-coordinates
    Y = prj(:,2);   % y-coordinates
    cfg.highlight   = [];
    cfg.normal      = setdiff(1:length(X),cfg.highlight);   
    % Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45
    y = 0.9*((X-min(X))/(max(X)-min(X))-0.5); % NOTE: x becomes y and y becomes x, in griddata is also reversed which makes x x and y y again
    x = 0.9*((Y-min(Y))/(max(Y)-min(Y))-0.5);
    interplimits = 'head';
      % 'electrodes' to furthest electrode
      % 'head' to edge of head
    % Find limits for interpolation:
    if strcmp(interplimits,'head')
      xmin = min(-.5,min(x)); xmax = max(0.5,max(x));
      ymin = min(-.5,min(y)); ymax = max(0.5,max(y));
    else
      xmin = max(-.5,min(x)); xmax = min(0.5,max(x));
      ymin = max(-.5,min(y)); ymax = min(0.5,max(y));
    end
    gridscale  = 100;                             % resolution
    interp     = 'v4';                            % 'linear','cubic','nearest','v4'
    xi         = linspace(xmin,xmax,gridscale);   % x-axis description (row vector)
    yi         = linspace(ymin,ymax,gridscale);   % y-axis description (row vector)
    delta      = xi(2)-xi(1);   
    cfg.labels = {UserData.(type).hdr.sensor_info.label};
    % Prepare data and create topoplot series
    for j = 1:length(alldata)
      % open a new figure for each topoplot series
      FigTag    = sprintf('TopoplotSet%g',j);
      hTopoplot = findobj('tag',FigTag);
      if isempty(hTopoplot)
        hTopoplot = figure('Name',['Topoplot: set ' num2str(j)],'tag',FigTag);
      else
        figure(hTopoplot);
      end
      clf
      % divide xlim into N subintervals
      % loop over subintervals
      N   = str2num(get(gethandle('nframes'),'string'));
      if isempty(N) || ~isint(N), N = 1; end
      tix = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2));
      t   = UserData.(type).time(tix); % time vector
      % alldata & t only contains tix times
      Nix = floor(length(tix)/N);
      for interval = 1:N
        subplot(1,N,interval), hold on
        % average each channel over this subinterval
        ix  = [1:Nix] + (interval-1)*Nix; % index to this subinterval
        if size(alldata{j},1) < nchan
          % get data from UserData
          tmptix = tix(ix); 
          for chan = 1:nchan
            dat = UserData.(type).data{j}(chan,tmptix);
            % preprocessing
            if get(gethandle('invertchk'),'value'), dat = -dat; end
            if get(gethandle('bpfilter'),'value')
              dat = blc(dat);
              dat = ts_freq_filt(dat',UserData.(type).hdr.sfreq,str2num(get(gethandle('bpfreq'),'string')),[0 0],'bandpass')';
            end
            if get(gethandle('lnfilter'),'value')
              dat = blc(dat);
              dat = ts_freq_filt(dat',UserData.(type).hdr.sfreq,60,5,'notch')';
            end            
            if get(gethandle('stdchk'),'value')
              StdDevWindow = min(50,floor(length(dat)/10));
              StdDevWindow = StdDevWindow + ~mod(StdDevWindow,2);
              npoints      = length(dat) - StdDevWindow + 1;
              offset       = floor(StdDevWindow/2);
              tmp          = zeros(size(dat));
              for l = 1:npoints
                tmp(l+offset) = std(dat(l:l+StdDevWindow-1));
              end
              dat = tmp;
              clear tmp
            end          
            if get(gethandle('envelope'),'value'), dat = abs(hilbert(dat)); end
%             if get(gethandle('smoothchk'),'value')
%               tmp = str2num(get(gethandle('zlim'),'string'));
%               if ~isempty(tmp) && numel(tmp)==1 && isint(tmp), npoint = tmp; else npoint = 10; end
%               dat = smooth(dat,npoint,'lowess');
%               clear tmp
%             end
            tmpdat(chan,:) = dat;
          end
          dat = mean(tmpdat,2); clear tmpdat
        else
          dat = mean(alldata{j}(:,ix),2);
        end
        tt  = t(ix);
        tlim= [min(tt) max(tt)];
        cfg.comment = sprintf('t=[%g %g]',tlim);
        % Interpolate data; NOTE: undo the reversal of x & y
        [Xi,Yi,Zi] = griddata(y,x,dat,yi',xi,interp);
          % griddata uses meshgrid to create evenly-spaced XI & YI
        % draw the topoplot for this interval
        draw_topoplot(Xi,Yi,Zi,x,y,delta,cfg);
        hold off
      end % end loop over subintervals
    end
  end
  if get(gethandle('avgbtn'),'value') && ~isempty(options.markers.type)
    codes = unique(options.markers.type);
    conds = options.markers.type;
    t     = UserData.Raw.time;
    Fs    = UserData.Raw.hdr.sfreq;
    label = get(allplots,'UserData');
    [chix,jnk] = match_str({UserData.Raw.hdr.sensor_info.label},label);
    tix   = nearest(t,xlim(1)):nearest(t,xlim(2));
    t     = t(tix);
    nchan = length(chix);
    nrow  = ceil(sqrt(nchan));
    ncol  = nrow;
    figure('Name',sprintf('average of events from %g - %g sec',min(t),max(t)),'color','w','units','normalized','position',[0 0 1 1]);
    legstr={};
    h = gethandle('evcodes');
    v = get(h,'value');
    for k = v%1:length(codes)
      color   = colors(mod(find(conds(k)==codes)-1,length(colors))+1);
      evntind = find(conds==codes(k));
      evtimes = setdiff(options.markers.time(evntind),options.rejection.rejects.time);
      evtimes = evtimes(evtimes>=xlim(1)&evtimes<=xlim(2));
      evtimes = evtimes - t(1);
      if isempty(evtimes), continue; end
      evsamps = round(evtimes*Fs);
      avg     = ts_matrix2data(alldata{1},'Fs',Fs,'sens',UserData.Raw.hdr.sensor_info(chix),'time',t,'continuous',1);
      avg     = ts_epoch_data(avg,'samples',evsamps,'prestim',1,'poststim',1);
      avg     = ts_trials2avg(avg);
      legstr{end+1} = sprintf('%g (n=%g)',codes(k),length(evsamps));
      pcnt    = 1;
      for jr = 1:nrow
        for jc = 1:ncol
          if pcnt > nchan, break; end
          subplot(nrow,ncol,pcnt);
          if ~ismember(avg.sensor_info(pcnt).label,options.rejection.rejects.badchans)
            plot(avg.averages.time,avg.averages.data(pcnt,:),color);
            hline(0,'k'); vline(0,'k'); box off;
          end
          title(avg.sensor_info(pcnt).label);
          axis tight; hold on; pcnt=pcnt+1;
        end
      end
    end
    legend(legstr{:});
  end
  clear alldata
elseif strcmp(options.mode,'epoched')
  % plot epoched data
  if length(uniqsens) > 1
    fprintf('Select only one channel before plotting epoched data.\n');
    return;
  end
  label = sens{1}; this = label;
  chan  = strmatch(label,labels,'exact');
  plot_wave     = ~isempty(wavtypes) && ismember(label,wavchans);
  plot_timefreq = ~isempty(tfrtypes) && ismember(label,tfrchans);  
  cntoverlay = 0;
  if plot_timefreq
%     if n==length(sens)
%       % turn off x-axis on wave plot since this will be below it
%       set(gca,'xtick',[]);
%     end
%     cnt = cnt + 1;
    % for each tfrtype
      type = 'TFR';
%       chan = strmatch(this,{UserData.(type).hdr.sensor_info.label},'exact');
      xlim = options.TFRplots.xlim;
      ylim = options.TFRplots.ylim;
      if isempty(tfrlims) || tfrlims(1)==tfrlims(2)
        zlim = options.TFRplots.zlim;
      else
        zlim = tfrlims;
      end
      tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2));      
      fix  = nearest(UserData.(type).frequencies,ylim(1)):nearest(UserData.(type).frequencies,ylim(2));
      if length(UserData.(type).data) > 1
        warning('Only showing the first TFR data set.  Cannot overlay TFRs.');
      end
      data = UserData.(type).data{1}(chan,tix,fix,:);
      if ~isreal(data),   data = abs(double(data)).^2; end
      
      ntrl  = size(data,4);
      maxrow = 10;
      nrow  = min([ceil(sqrt(ntrl)) maxrow]);
      ncol  = nrow;
      ntrl  = min(ntrl,nrow*ncol);
      trial_range = str2num(get(gethandle('nframes'),'string'));
      if length(trial_range)==2 && isint(trial_range)
        ntrl    = min(ntrl,diff(trial_range)+1);
        nrow    = min([nrow ceil(sqrt(ntrl))]);
        ncol    = nrow;
        trials  = trial_range(1):trial_range(2);
      else
        trial_range = [1 ntrl];
        trials      = 1:ntrl;
      end
      set(gethandle('nframes'),'string',sprintf('[%g %g]',trial_range));
%       % adjust num rows and num cols so that 
%       if nrow~=ncol
%         ncol = ceil((nrow*ncol)/ceil(.8*nrow));
%         nrow = ceil(.8*nrow);
%       end
      % set offsets to zero to turn them off
      yoffset = 2;                                    % correction for axes extending outside of frame vertically
%       nrow    = nrow + yoffset;
      
      xoffset = 1;                                    % correction for showing y-axis labels on in first column
%       ncol    = ncol + xoffset;
      
      xstepsize = 1 / (ncol+xoffset);
      ystepsize = 1 / (nrow+yoffset);

      xstep = mod((1:ncol+xoffset)-1,ncol);
      ystep = mod((1:nrow+yoffset)-1,nrow)+1;
      xpos  = xstep*xstepsize;
      ypos  = 1-ystep*ystepsize;
      xpos  = xpos(1:end-xoffset) + (xoffset/2)*xstepsize;  % correction; does nothing if offset=0
      ypos  = ypos(1:end-yoffset) - (yoffset/2)*ystepsize;  % correction; does nothing if offset=0
%       nrow  = nrow - yoffset;
%       ncol  = ncol - xoffset;

% screensize = get(0,'screensize');
% figure('position',screensize);
      for n = 1:ntrl
        if n > nrow*ncol, break; end
        type = 'TFR';
        tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2));      
        xi  = mod(n-1,ncol)+1;%xi  = ncols-mod(i,ncols);
        yi  = floor((n-1)./ncol)+1;
%         allplots(n) = subplot('Position',[xpos(xi) ypos(yi) 1/ncol 1/nrow],'parent',gethandle('plots')); set(gca,'units','normalized');
        allplots(n) = subplot('Position',[xpos(xi) ypos(yi) xstepsize ystepsize],'parent',gethandle('plots')); set(gca,'units','normalized');%findobj('tag','MainFig1')
%         set(gca,'parent',findobj('tag','MainFig1'),'Units','normalized','Position',[xpos(xi) ypos(yi) 1/ncol 1/nrow]);
%         allplots(n) = subplot(nrow,ncol,n,'parent',gethandle('plots'));
        set(gca,'UserData','TFR');
        if n == 1, hold off; end
        dat = squeeze(data(:,:,:,trials(n)))';
        himage = imagesc(UserData.(type).time(tix),UserData.(type).frequencies(fix),dat,zlim);
        set(himage,'ButtonDownFcn',@onclick); hold on
        set(gca,'YDir','normal')
        foi = UserData.(type).frequencies(fix); 
        nn  = 4;%length(get(gca,'ytick'));
        yticklabels = foi(round(linspace(1,length(foi),nn)));
%         dy  = floor(length(foi)/nn); 
%         yticklabels = foi(1:dy:length(foi));%foi(dy:dy:length(foi));
        nn  = length(yticklabels); 
        yticks = linspace(ylim(1),ylim(2),nn);
%         dy  = (max(foi)-min(foi))/nn; 
%         dy  = diff(ylim)/nn;
%         yticks      = min(foi):dy:max(foi);%min(foi)+dy:dy:max(foi); 
        set(gca,'ytick',yticks,'yticklabel',yticklabels);        
        set(gca,'xlim',xlim);
        set(gca,'ylim',ylim);
        set(gca,'clim',zlim);      
        % add vertical line at t=0
        if xlim(1)<0 && xlim(2)>0, vline(0,'color','k','linewidth',.5); end
        % add lines & '+' at x and y axis tick labels in all axes
        % y-axis
        vline(xlim(1),'color','k','linewidth',3);
        scatter(xlim(1)*ones(1,length(yticks)),yticks,'+k','linewidth',.5,'Clipping','on','sizedata',30);
        % x-axis
        hline(ylim(2),'color','k','linewidth',3);
        nticks = 5;
        xticks = unique(round(linspace(xlim(1),xlim(2),nticks)*1000))/1000;
        scatter(xticks,ylim(2)*ones(1,length(xticks)),'+k','linewidth',.5,'Clipping','on','sizedata',30);
        box off
        % set tick labels
        if xi==1 && yi==nrow
          set(gca,'xtick',xticks);
          set(gca,'xticklabel',xticks);          
        else
          set(gca,'xticklabel',[]); 
        end
        if n > 1, set(gca,'yticklabel',[]); end        
        % add text indicating the trial number
        ax   = axis;
        xtxt = ax(1)+.1*(ax(2)-ax(1));
        ytxt = ax(3)+.8*(ax(4)-ax(3));
        text(xtxt,ytxt,num2str(trials(n)),'fontsize',14);

        % add waveform if given
        if plot_wave
          cntoverlay   = 0;
          for k  = 1:length(wavtypes)
            type = wavtypes{k};
            tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2));
            nset = length(UserData.(type).data);
            for j = 1:nset
              if n == 1
                maxval = max(max(max(abs(UserData.(type).data{j}(chan,tix,trials)))));
              end
              cntoverlay = cntoverlay + 1;
              lnstyle    = [colors(mod(cntoverlay,length(colors))) lntype{mod(k,length(colors))}];
              dat        = UserData.(type).data{j}(chan,tix,trials(n));
              if get(gethandle('invertchk'),'value'), dat = -dat; end
              if get(gethandle('bpfilter'),'value')
                dat = blc(dat);
                dat = ts_freq_filt(dat',UserData.(type).hdr.sfreq,str2num(get(gethandle('bpfreq'),'string')),[0 0],'bandpass')';
              end
              if get(gethandle('lnfilter'),'value')
                dat = blc(dat);
                dat = ts_freq_filt(dat',UserData.(type).hdr.sfreq,60,5,'notch')';
              end            
              if get(gethandle('stdchk'),'value')
                StdDevWindow = min(50,floor(length(dat)/10));
                StdDevWindow = StdDevWindow + ~mod(StdDevWindow,2);
                npoints      = length(dat) - StdDevWindow + 1;
                offset       = floor(StdDevWindow/2);
                tmp          = zeros(size(dat));
                for l = 1:npoints
                  tmp(l+offset) = std(dat(l:l+StdDevWindow-1));
                end
                dat = tmp;
                clear tmp
              end          
              if get(gethandle('envelope'),'value'), dat = abs(hilbert(dat)); end
              % get factors for waveform range
              if ~get(gethandle('autoscale'),'value')
                scale_factor = max(max(max(abs(UserData.(type).data{j}(chan,tix,trials(n)))))) / maxval;
                wave_ylim    = [-maxval maxval];
              else
                scale_factor = 1;
                wave_ylim = [min(dat) max(dat)];
              end              
              % add y-axis labels for wave plot
              if (xi==ncol && yi==1) && ~get(gethandle('autoscale'),'value')
                vline(xlim(2),'color','k','linewidth',2);
                yticks = linspace(ylim(1),ylim(2),length(yticks));
                scatter(xlim(2)*ones(1,length(yticks)),yticks,'+k','linewidth',.5,'Clipping','on','sizedata',30);
                clear wave_yticklabel
                wave_yticks = linspace(wave_ylim(1),wave_ylim(2),length(yticks));
                for tmpind=1:length(wave_yticks)
                  wave_yticklabels{tmpind} = sprintf('%0.2i',wave_yticks(tmpind));
                end
                set(gca,'yaxislocation','right','ytick',yticks,'yticklabel',wave_yticklabels);
              elseif get(gethandle('autoscale'),'value')
                % todo: add labels indicating subplot-specific y-limits
              end
              % rescale and shift data so that it plots on top of TF data
              dat = dat - mean(dat);      % center on y=0
              dat = dat / max(abs(dat));  % scale to within [-1 1]
              dat = dat * diff(ylim)/2;   % scale to within +/- half subplot height
              dat = dat * scale_factor;   % rescale if all trials scaled the same
              dat = dat + mean(ylim);     % center y on subplot
              % plot rescaled and shifted waveform
              plot(UserData.(type).time(tix),dat,lnstyle);
            end
          end
          options = drawmarkers(options,xlim);
        end
      end
elseif plot_wave
    for k  = 1:length(wavtypes)
      type = wavtypes{k};
      xlim = options.waveplots.xlim;
      trials = str2num(get(gethandle('zlim'),'string'));
      if length(trials>=2) && isint(trials(1)) && isint(trials(2))
        if length(trials)==2
          trials = trials(1):trials(2);
        end
      else
        trials = [];
      end
      tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2));
      nset = length(UserData.(type).data);
      for j = 1:nset
        cntoverlay = cntoverlay + 1;
        data = UserData.(type).data{j}(chan,tix,:);
        ylim = [min([data(:)]),max([data(:)])];
        if all(ylim==0), ylim=[-1 1]; end
        if ylim(1)==ylim(2),ylim(1)=.5*ylim(1); ylim(2)=1.5*ylim(2); end
        % what trials should we plot?
        if isempty(trials)
          trials = 1:min(100,size(data,3));
        end
        trials  = trials(trials<=size(data,3));
        ntrl    = length(trials);
        nrow    = ceil(sqrt(ntrl));
        ncol    = nrow;
        lnstyle = [colors(mod(cntoverlay,length(colors))) lntype{mod(k,length(colors))}];        
        for trl = 1:ntrl
          n = trials(trl);
          allplots(trl) = subplot(nrow,ncol,trl,'parent',gethandle('plots'));
          dat  = data(:,:,n);
            if get(gethandle('invertchk'),'value'), dat = -dat; end
            if get(gethandle('bpfilter'),'value')
              dat = blc(dat);
              dat = ts_freq_filt(dat',UserData.(type).hdr.sfreq,str2num(get(gethandle('bpfreq'),'string')),[0 0],'bandpass')';
            end
            if get(gethandle('lnfilter'),'value')
              dat = blc(dat);
              dat = ts_freq_filt(dat',UserData.(type).hdr.sfreq,60,5,'notch')';
            end            
            if get(gethandle('stdchk'),'value')
              StdDevWindow = min(50,floor(length(dat)/10));
              StdDevWindow = StdDevWindow + ~mod(StdDevWindow,2);
              npoints      = length(dat) - StdDevWindow + 1;
              offset       = floor(StdDevWindow/2);
              tmp          = zeros(size(dat));
              for l = 1:npoints
                tmp(l+offset) = std(dat(l:l+StdDevWindow-1));
              end
              dat = tmp;
              clear tmp
            end          
            if get(gethandle('envelope'),'value'), dat = abs(hilbert(dat)); end
%             if get(gethandle('smoothchk'),'value')
%               tmp = str2num(get(gethandle('zlim'),'string'));
%               if ~isempty(tmp) && numel(tmp)==1 && isint(tmp), npoint = tmp; else npoint = 25; end
%               dat = smooth(dat,npoint,'lowess');
%               clear tmp
%             end
          if get(gethandle('histchk'),'value')
            try
              hist(dat,round(min(100,length(dat)/100)));
            end
            axis tight
            xlim = get(gca,'xlim');
            ylim = get(gca,'ylim');            
          else
            plot(UserData.(type).time(tix),dat,lnstyle);
      %       title(sprintf('trial %g',n));
            hold on            
            set(gca,'xlim',xlim);
            if get(gethandle('autoscale'),'value')
              axis tight
            else
              set(gca,'ylim',ylim);
            end
            if trl < ntrl
              set(gca,'xtick',[]);
            end
%             set(gca,'ytick',[]);
            if xlim(1)<0 && xlim(2)>0, vline(0,'k'); end
            box off
            grid on
            if get(gethandle('eventchk'),'value') && ~isempty(options.events(1).label) && ismember(this,{options.events.label})
              chix = strmatch(this,{options.events.label},'exact');
              typs = unique(options.events(chix).type);
              if trl==1
                for typ = 1:length(typs)
                  if size(dat,3)==sum(options.events(chix).type==typ)
                    break;
                  elseif typ==length(typs)
                    typ = nan;
                  end
                end
              end
              if ~isnan(typ)
                title(num2str(options.events(chix).time(trl)));
              end
            end
          end
          if j==nset && k==length(wavtypes), hold off; end
          if get(gethandle('xcorrchk'),'value') || get(gethandle('topoplotchk'),'value') %|| get(gethandle('chanoverlaychk'),'value') || get(gethandle('imagesc'),'value')
            % collect data for later calculation
            alldata{j}(trl,:) = dat;
          end   
          options = drawmarkers(options,xlim);
        end
      end
    end
  end
  if get(gethandle('xcorrchk'),'value')
    % calculate cross-correlations and plot them in a new figure
    for j = 1:length(alldata)
      % open a new figure for each set of cross-correlations
      hfig  = findobj('tag','xcorrFig');
      if isempty(hfig)
        figure('Name',['Normalized cross-correlations: set ' num2str(j)],'tag','xcorrFig');
      else
        figure(hfig);
      end
      for n = 1:ntrl
        for m = 1:ntrl
          if m < n, continue; end
          subplot(ntrl,ntrl,(n-1)*ntrl+m)
          y1 = alldata{j}(n,:);
          y2 = alldata{j}(m,:);
          [xc,lags] = xcorr(y1,y2,'coeff');
          plot(lags/UserData.Raw.hdr.sfreq,xc); axis tight; vline(lags(ceil(length(lags)/2))/UserData.Raw.hdr.sfreq,'k'); hline(0,'k');
          if n==1, title(sprintf('trial %g',trials(m)));  end
          if m==n, ylabel(sprintf('trial %g',trials(m))); end
          clear y1 y2
        end
      end
    end
  end
%   if get(gethandle('chanoverlaychk'),'value')
%     try
%       j    = 1;
%       tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2)); 
%       hfig = findobj('tag','OverlayFig');
%       if isempty(hfig)
%         figure('tag','OverlayFig');
%       else
%         figure(hfig);
%       end
%       for trl = 1:size(alldata{1},1)
%         tmp(trl,:) = alldata{j}(trl,:);
%       end
%       plot(UserData.(type).time(tix),tmp);
%       axis tight
%       clear tmp
%     catch
%       close(gcf)
%     end
%   end     
  if get(gethandle('imagesc'),'value')    
    j    = 1;
    tix  = nearest(UserData.(type).time,xlim(1)):nearest(UserData.(type).time,xlim(2)); 
    hfig = findobj('tag','imagescFig');
    if isempty(hfig)
      figure('tag','imagescFig');
    else
      figure(hfig);
    end
    for trl = 1:size(alldata{j},1)
      tmp(trl,:) = alldata{j}(trl,:);
    end    
    imagesc(UserData.(type).time(tix),1:size(alldata{j},1),tmp);
    axis tight;
    clear tmp    
  end  
else
  error('Mode not recognized; aborting plot.');
end
options.subplots.handles = allplots;
% if options.waveplots.waitbar
%   close(hwait)
% end
function draw_topoplot(Xi,Yi,Zi,x,y,delta,cfg)
% Take data within head
rmax       = .5;
mask       = (sqrt(Xi.^2+Yi.^2) <= rmax);
Zi(mask==0)= NaN;

% Draw topoplot on head
numcontours = 6;
shading     = 'flat';                         % 'flat' or 'interp'
contour(Xi,Yi,Zi,numcontours,'k');
h = surface(Xi-delta/2,Yi-delta/2,zeros(size(Zi)),Zi,'EdgeColor','none', 'FaceColor',shading);

% calculate colormap limits
zmin = min(Zi(:));
zmax = max(Zi(:));
caxis([zmin zmax])

% draw head
  % Define the outline of the head, ears and nose:
  l     = 0:2*pi/100:2*pi;
  tip   = rmax*1.15; base = rmax-.004;
  EarX  = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
  EarY  = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
  % Plot head, ears, and nose:
  plot(cos(l).*rmax, sin(l).*rmax, 'color', cfg.hcolor, 'Linestyle', '-', 'LineWidth', cfg.hlinewidth);
  plot([0.18*rmax;0;-0.18*rmax], [base;tip;base], 'Color', cfg.hcolor, 'LineWidth', cfg.hlinewidth);
  plot( EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth)
  plot(-EarX, EarY, 'color', cfg.hcolor, 'LineWidth', cfg.hlinewidth)

  % draw electrodes
  hp2 = plot(y(cfg.normal),    x(cfg.normal),    cfg.emarker,  'Color', cfg.ecolor,  'markersize', cfg.emarkersize);
  hp2 = plot(y(cfg.highlight), x(cfg.highlight), cfg.hlmarker, 'Color', cfg.hlcolor, 'markersize', cfg.hlmarkersize,'linewidth', cfg.hllinewidth);
  % add labels
%   for ch = 1:length(cfg.labels)
%     text(y(ch), x(ch), cfg.labels{ch}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
%                                   'Color', cfg.ecolor, 'FontSize', cfg.fontsize);%cfg.efsize);
%   end
  x_COMNT =  0.6; 
  y_COMNT = -0.6;
  HorAlign = 'right';
  VerAlign = 'bottom';
% Write comment:
% if isfield(cfg, 'comment') 
%   if strcmp(cfg.commentpos, 'title')
%     title(cfg.comment, 'Fontsize', cfg.fontsize);
%   else
    text(x_COMNT, y_COMNT, cfg.comment, 'Fontsize', cfg.fontsize, 'HorizontalAlignment', HorAlign, 'VerticalAlignment', VerAlign);
%   end
% end

% colorbar
% hold off
axis off
xlim([-.6 .6]);
ylim([-.6 .6]);
  
function type = gettype(data)
datatype = ts_object_info(data,'verbose',0);
if strcmp(datatype,'epoch_data') || strcmp(datatype,'cont_data') || strcmp(datatype,'avg_data')
  type = 'Raw';
elseif strcmp(datatype,'timefreq_data')
  type = 'TFR';
else
  type = '';
end

function ZoomFunction(src,evnt)
global options
tfrtypes = options.TFRplots.datatypes;
% tfrchans = options.TFRplots.labels;
% wavtypes = options.waveplots.datatypes;
% wavchans = options.waveplots.labels;
% sens     = unique({wavchans{:} tfrchans{:}});
% plot_timefreq = ~isempty(tfrtypes);
% tfr_plot_ind  = find(ismember(sens,tfrchans));
% if ~isempty(tfrtypes) || strcmp(get(gco,'Tag'),'chanlabels'), return; end % TODO: add check for TF data & scale clim appropriately
if strcmp(get(gco,'Tag'),'chanlabels'), return; end % TODO: add check for TF data & scale clim appropriately
hpanel = gethandle('plots');
hplots = get(hpanel,'Children');
if get(gethandle('imagesc'),'value')        % adjust image & wave plots separately
  CLIM = get(hplots(2),'clim'); if isempty(CLIM), return; end
  YLIM = get(hplots(1),'ylim'); if isempty(YLIM), return; end
  if evnt.VerticalScrollCount < 0           % zoom in
    set(hplots(2),'clim',CLIM/1.5);
    set(hplots(1),'ylim',YLIM/1.5);
  else                                      % zoom out
    set(hplots(2),'clim',CLIM*1.5);
    set(hplots(1),'ylim',YLIM*1.5);    
  end
elseif get(gethandle('autoscale'),'value') && isempty(tfrtypes) % adjust zoom in selected channel
  YLIM = get(gca,'ylim'); if isempty(YLIM), return; end
  if evnt.VerticalScrollCount < 0           % zoom in
    set(gca,'ylim',YLIM/1.5);
  else                                      % zoom out
    set(gca,'ylim',YLIM*1.5);
  end
elseif ~isempty(tfrtypes) && strcmp(options.mode,'epoched')          % adjust zoom in all channels
  YLIM = get(hplots(1),'ylim'); if isempty(YLIM), return; end
  CLIM = get(hplots(1),'clim'); if isempty(CLIM), return; end
  if evnt.VerticalScrollCount < 0           % zoom in
    set(hplots,'clim',CLIM/1.5);
  else                                      % zoom out
    set(hplots,'clim',CLIM*1.5);
  end
  set(findobj('tag','zlim'),'string',sprintf('[%g,%g]',CLIM));
else                                                                 % adjust zoom in all channels
  YLIM = get(hplots(1),'ylim'); if isempty(YLIM), return; end
  if evnt.VerticalScrollCount < 0           % zoom in
    set(hplots,'ylim',YLIM/1.5);
  else                                      % zoom out
    set(hplots,'ylim',YLIM*1.5);
  end
end

function onclick(src,evnt)
% Note: this function is called every time you click on a plot or line
% Purpose: interactive plotting element
% (~eventchk, ~rejcheck) -> select current axis
% ( eventchk, ~rejcheck) -> add/remove markers (channel independent events)
%                           (vlines for waveplots; hlines for tfrplots)
% (~eventchk,  rejcheck) -> add/remove channel from reject list (then replot)
% ( eventchk,  rejcheck) -> add/remove markers (TODO: instead, add/remove marker from reject list only)

global options
hpanel   = gethandle('plots');
children = get(hpanel,'Children');
childpos = get(children,'position');
xlim     = get(gca,'xlim');
ylim     = get(gca,'ylim');
if iscell(childpos), x = childpos{1}(1); else x = childpos(1); end
y  = 0;
if iscell(childpos), w = childpos{1}(3); else w = childpos(3); end
h  = 1;
P  = get(gca,'Currentpoint');
p  = P(1,1);
% xL = (p-xlim(1))/diff(xlim);
this = get(gca,'UserData');
if ~get(gethandle('eventchk'),'value') && ischar(this)
  if ismember(this,options.rejection.rejects.badchans)
    % restore channel
    chix = strmatch(this,options.rejection.rejects.badchans,'exact');
    options.rejection.rejects.badchans = options.rejection.rejects.badchans(setdiff(1:length(options.rejection.rejects.badchans),chix));
    UpdateManager(src,evnt);
  elseif get(gethandle('rejcheck'),'value')
    % reject channel
    options.rejection.rejects.badchans = {options.rejection.rejects.badchans{:} this};
    UpdateManager(src,evnt);
%     return
  end
  return
elseif get(gethandle('eventchk'),'value') && ...
    ~isempty(options.events(1).label) && ...
    ischar(this) && ismember(this,{options.events.label})
  chix = strmatch(this,{options.events.label},'exact');
  typs = unique(options.events(chix).type);
%   chans = options.waveplots.labels; if ~iscell(chans), chan={chans}; end
%   [s1 s2] = match_str({options.events.label},chans);
%   tmpevents = options.events(s1);
  if length(typs) >= 4
    fprintf('Sorting channels by event times...');
    % find 3rd type closest to point p
    T3       = options.events(chix).time(options.events(chix).type==typs(3));
    evnum3   = nearest(T3,p);
    T3       = T3(evnum3);
    % find 4th type closest to point p
    T4       = options.events(chix).time(options.events(chix).type==typs(4));
    evnum4   = nearest(T4,p);
    T4       = T4(evnum4);
    wavchans = options.waveplots.labels;
    if (T3-p) < (T4-p)
      % clicked start of epoch => sort by type 1 times b/w T3 & T4
      sorttype = typs(1);
      T4 = options.events(chix).time(options.events(chix).type==typs(4));
      T4 = T4(T4>T3);
      T4 = T4(nearest(T4,p));
    else
      % clicked end of epoch => sort by type 2 times b/w T3 & T4
      sorttype = typs(2);
      T3 = options.events(chix).time(options.events(chix).type==typs(3));
      T3 = T3(T3<T4);
      T3 = T3(nearest(T3,p));
    end             
    % sort by times b/w T3 & T4
    pktimes = inf(1,length(wavchans));
    for k = 1:length(wavchans)
      if ismember(wavchans{k},{options.events.label})
        tmpchix = strmatch(wavchans{k},{options.events.label},'exact');
        tmptime = options.events(tmpchix).time(options.events(chix).type==sorttype);
        tmptime = tmptime(tmptime>=T3 & tmptime<=T4);
        if isempty(tmptime)
          pktimes(k) = inf;
        else
          % keep first event of appropriate type after T3
          pktimes(k) = tmptime(1);
        end
      end
    end % end loop over wavchans
    [jnk,sortidx] = sort(pktimes);
    options.waveplots.labels = wavchans(sortidx);
    options.TFRplots.labels  = wavchans(sortidx);
    allsensors = get(gethandle('chanlabels'),'string');
    hiddensens = setdiff(allsensors,wavchans);
    allsensors = {options.waveplots.labels{:} hiddensens{:}};
    set(gethandle('chanlabels'),'string',allsensors);
    set(gethandle('chanlabels'),'value',1:length(wavchans));
    UpdateManager(src,evnt);
    % TODO: add delay map topoplot and streamlines for sorted channels
    fprintf('done.\n');
    return
  end
end
colors = 'rkbgymrkbgymrkbgymrkbgymrkbgymrkbgymrkbgymrkbgymrkbgym';
MinDistance = diff(xlim/200); % eps
if ~any(abs(options.markers.time - p) < MinDistance)
  h = gethandle('evcodes');
  color = colors(get(h,'value'));
  % add this marked event to the markers vector
  options.markers.time(end+1)   = p;
  % use default event code of 1 for now
  options.markers.type(end+1)   = get(h,'value');%1;
  % draw vline
  ax = axes('Position',[x y w h],'parent',hpanel);
  set(gca,'xlim',xlim,'ButtonDownFcn',@onclick);
  options.markers.handle(end+1)   = vline(p,'Color',color,'ButtonDownFcn',@onclick,'LineWidth',.1);
  options.markers.position(end+1) = p; % this may not be necessary
  set(ax,'Xtick',[],'Ytick',[],'Visible','off');  
  fprintf('Added marker at time = %g sec\n',p);
  if isequal(get(gca,'ylim'),options.TFRplots.ylim)
    currpos = get(gca,'position');
    ax = axes('Position',currpos,'parent',hpanel);
    set(gca,'xlim',xlim,'ButtonDownFcn',@onclick);
    set(gca,'ylim',ylim);
    hline(options.TFRplots.ylim(2)-P(1,2),'Color',color,'ButtonDownFcn',@onclick,'LineWidth',.1);
    set(ax,'Xtick',[],'Ytick',[],'Visible','off');
  end  
else
  % remove this event from the markers
  ind = find(abs(options.markers.time - p) < MinDistance);
  tmptime = options.markers.time(ind);
  options.markers.time(ind) = [];
  options.markers.type(ind) = [];
  ind = find(abs(options.markers.position - p) < MinDistance);
  options.markers.position(ind) = [];
  for k = 1:length(ind)
    try
      delete(options.markers.handle(ind(k)))
      if k==1, fprintf('Removed marker at time = %g sec\n',p); end
    catch
%       options.markers.handle(ind(k)) = [];
    end
  end
  options.markers.handle(ind)   = [];
  if get(gethandle('rejcheck'),'value') && ismember(tmptime,options.rejection.rejects.time)
    % remove marker from reject list & re-draw markers
    tmpind = find(tmptime==options.rejection.rejects.time);
    options.rejection.rejects.time(tmpind) = [];
    options.rejection.rejects.badtrials(tmpind) = [];
  end
end
if get(gethandle('rejcheck'),'value')
  UpdateManager(src,evnt);
end

function options = drawmarkers(options,xlim)
% Note: this function is called every time the plot panel is updated
% Purpose: draw the desired markers & patches around trials & rejects
% (~eventchk, ~rejcheck) -> do nothing
% ( eventchk, ~rejcheck) -> draw vlines at all markers
% (~eventchk,  rejcheck) -> draw vlines at rejected markers
% ( eventchk,  rejcheck) -> draw vlines at all markers & patches around rejected markers

colors = 'rkbgymrkbgymrkbgymrkbgymrkbgymrkbgymrkbgymrkbgymrkbgym';
if isempty(options.markers.time) || ~(get(gethandle('eventchk'),'value') || get(gethandle('rejcheck'),'value'))
  return
end
ind = find(options.markers.time>=xlim(1) & options.markers.time<=xlim(2));
if isempty(ind)
  return
end
hpanel   = gethandle('plots');
children = get(hpanel,'Children');
childpos = get(children,'position');
if iscell(childpos), x = childpos{1}(1); else x = childpos(1); end
y   = 0;
if iscell(childpos), w = childpos{1}(3); else w = childpos(3); end
h   = 1;
pos = options.markers.time(ind);
if get(gethandle('rejcheck'),'value') && ~get(gethandle('eventchk'),'value')
  ind2 = ismember(pos,options.rejection.rejects.time);
  pos  = pos(ind2);
  ind  = ind(ind2);
end
conds = options.markers.type(ind);
codes = unique(options.markers.type);
ax  = axes('Position',[x y w h],'parent',hpanel); hold on
set(gca,'xlim',xlim,'ylim',[0 .93],'ButtonDownFcn',@onclick);
for k = 1:length(pos)
  htxt  = [];
  color = colors(mod(find(conds(k)==codes)-1,length(colors))+1);
  if get(gethandle('eventchk'),'value')
    options.markers.handle(k) = vline(pos(k),'Color',color,'ButtonDownFcn',@onclick,'LineWidth',.1);
    htxt = text(pos(k),.94,num2str(conds(k)));
    set(htxt,'fontsize',8);
  end
  if ismember(pos(k),options.rejection.rejects.time) && get(gethandle('rejcheck'),'value')
    if get(gethandle('eventchk'),'value')
      x1 = pos(k) + options.rejection.parms.prestim;
      x2 = pos(k) + options.rejection.parms.poststim;
      hh = patch([x1 x1 x2 x2],[y y+h y+h y],color);
      set(hh,'EdgeColor',get(hh,'FaceColor'),'FaceAlpha',.2,'EdgeAlpha',.2,'ButtonDownFcn',@onclick);
    else
      options.markers.handle(k) = vline(pos(k),'Color',color,'ButtonDownFcn',@onclick,'LineWidth',.1);
    end
    if isempty(htxt)
      htxt = text(pos(k),.94,num2str(conds(k)));
      set(htxt,'fontsize',8);
    end      
  end
end
set(gca,'ylim',[0 1]);
options.markers.position = pos;
set(ax,'Xtick',[],'Ytick',[],'Visible','off');  

function h = gethandle(tag)
global options
h = findobj('tag',tag);
h = h(options.FigIndex);

function arrowpress(src,evnt)     
if strcmpi('uparrow',evnt.Key)
  ShiftChannel(src,evnt,'up');
elseif strcmpi('downarrow',evnt.Key)
	ShiftChannel(src,evnt,'down');
elseif strcmpi('leftarrow',evnt.Key)
  ShiftTime(src,evnt,'left');
elseif strcmpi('rightarrow',evnt.Key)
  ShiftTime(src,evnt,'right');
end

%% NOTES
% 
% function onclick(src,evnt)
% % Note: this function is called every time you click on a plot or line
% % Purpose: interactive plotting element
% % (~eventchk, ~rejcheck) -> select current axis
% % ( eventchk, ~rejcheck) -> add/remove markers (channel independent events)
% %                           (vlines for waveplots; hlines for tfrplots)
% % (~eventchk,  rejcheck) -> add/remove channel from reject list
% % ( eventchk,  rejcheck) -> add/remove markers (TODO: instead, add/remove marker from reject list)
% 
% function options = drawmarkers(options,xlim)
% % Note: this function is called every time the plot panel is updated
% % Purpose: draw the desired markers & patches around trials & rejects
% % (~eventchk, ~rejcheck) -> do nothing
% % ( eventchk, ~rejcheck) -> draw vlines at all markers
% % (~eventchk,  rejcheck) -> draw vlines at rejected markers
% % ( eventchk,  rejcheck) -> draw vlines at all markers & patches around rejected markers
% 
% function UpdateRejection
% % Note: this function is always called before updating the plots panel
% % Purpose: use thresholding to reject trials in selected non-rejected chans
% % (given markers, rejcheck) -> calc trial rejects & add to rejects list
% % (else) -> do nothing
