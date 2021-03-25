function [epoch_data, reject_data] = ts_manual_reject(epoch_data,varargin)
% Purpose: To perform visual rejection of neuroimaging data
% 
% Usage: [data results] = ts_manual_reject(data,'key1',val1, ... )
%
% Example: [epoch_data, reject_data] = ts_manual_reject(epoch_data, ...
%     'chantype', 'grad', ...
%     'events', [], ...
%     'adv_flag',1,...
%     'numstd', 6, ...
%     'reject_flag', 1, ...    %default to 1, this cause the data from the specified rejectfile to be applied to the current data set
%     'rootoutdir', rootoutdir, ...
%     'prefix','Vis_Rej', ...
%     'prescale_grad',10^24,...
%     'rawscale_grad',10^15,...
%     'reject_file', rejectfile, ...
%     'trial_info_flag', 0, ...
%     'save_epoch_flag', 0, ...
%     'trial_info_rej_flag', 0, ...
%     'toilim', [-.50 .50],...
%     'toilim_raw', [-.50 .50], ...
%     'raw_file', fname_raw, ...
%     'raw_flag', 1, ...
%     'ylim', [0 30], ...
%     'ylim_raw', [-800 800]););
% %
% Inputs: 
%     epoch_data: timesurfer data struct of epoch data 
%
% Outputs:
%     epoch_data: timesurfer data struct of epoch data
%       -AND-
%     reject_data: timesurfer data struct of rejected data
%
% Parameters: default : options :  description 

%     events :[] :: Define in order to only select a few events not the
%     entire data set

%     adv_flag :0 :0,1: Set to use inspect mode

%     chantype : 'all' : 'all', 'grad', 'mag', 'eeg', 'eog' : Define to only view a
%     specific channel type

%     numstd : 6 : 1-Inf : Define number of standard deviations before a trial crosses channel threshold 

%     prescale_mag : 10^15 :: default converts to fT, change to set scale
%     value for mags

%     prescale_grad : 10^13 :: converts to fT/cm change for grad scaling

%     prescale_eeg  : 10^6  :: converts to uV by default, change for eeg
%     scaling

%     prescale_eog  : 10^6  ::  converts to uV by default, change for eog

%    rawscale_mag  : prescale_mag :: defaults prescale value: use for
%   scaling raw data differently

%    rawscale_grad : prescale_grad :: defaults prescale value: use for
%   scaling raw data differently

%    rawscale_eeg : prescale_eeg :: defaults prescale value: use for
%   scaling raw data differently

%    rawscale_eog : prescale_eog :: defaults prescale value: use for
%   scaling raw data differently

%    keepbadchans : 0 : 0,1 : Keep channels marked as bad in sensor info
%    when set to 1

%    keepbadtrials : 1 : 0,1 : Keep Trials marked as bad in the trial info
%    section when set to 1

%    reject_flag : 1 : 0,1 :    %default to 1, this cause the data from the
%    specified rejectfile to be applied to the current data set

%    rootoutdir : pwd :: Changing this will set where the reject file is saved 

%    prefix : 'Untitled1':: Changing this will set a particular prefix on
%    your saved reject file

%    reject_file :[] :: specifies the location of a reject file for this data set, the same file will be written to for saving 

%    trial_info_flag : 1 : 0,1 : Set to 1 by default, set to 0 if no trial
%    info section or wish it to be ignored

%   save_epoch_flag : 1 :0,1 : By default will update epoch_data's to mark
%   rejected trials as bad

%   trial_info_rej_flag : 1 : 0,1 : set to consider badtrials from
%   trial_info as rejected trials

%   usr_file : 0 : 0,1 : Set to 1 by default will load in screen positions
%   and save screen positions

%   toilim : [] : -Inf,Inf :  %Set value to define time window seen for given epoch, by default is set to view entire epoch

%   toilim_raw : [] :-Inf Inf : %Set value to define time window seen for given raw epoch, by default is set to view entire raw epoch

%   raw_file : [] :: defines location of raw data file 

%   raw_flag : 1 : 0,1 : To view raw data set to 1 , if no raw file
%   specified with use epoch data as raw data

%   ylim : [0 180] : [-Inf,Inf] : Set to define ylimits of data set view,
%   default is for band based data

%   ylim_raw : [-400 400] : [-Inf,Inf] : Set to define ylimits for viewing
%   raw data


% Created by Andrew Schulman on Jul-2010
% Modified last by ADS on 15-Apr-201

if (get(0,'ScreenPixelsPerInch') < 90) ... % if resolution from matlab is incorrect
    set(0,'ScreenPixelsPerInch', 90); ...
    sprintf('done');
end

parms = mmil_args2parms(varargin, { ...
  'events',[],[],... %need to add in event code support
  'adv_flag',[0],[0,1],...
  'chantype','all',[],...
  'numstd',6,[1,Inf],...
  'prescale_mag'       ,10^15,[],...  % to fT
  'prescale_grad'      ,10^13,[],... % to fT/cm
  'prescale_eeg'       ,10^6,[],...   % to uV
  'prescale_eog'       ,10^6,[],...   % to uV
  'rawscale_mag'       ,[],[],...
  'rawscale_grad'      ,[],[],...
  'rawscale_eeg'       ,[],[],...
  'rawscale_eog'       ,[],[],...
  'keepbadchans',0,{0,1},...   %Default to remove badchans
  'keepbadtrials',1,{0,1},...  % Default is to keep bad trials
  'reject_flag',1,[0,1],...    %default to 1, this cause the data from the specified rejectfile to be applied to the current data set
  'rootoutdir',[pwd],[],...
  'prefix','Untitled1',[],...
  'reject_file',[],[],...
  'trial_info_flag',1,[0,1],...
  'usr_file',1,[0,1],...
  'save_epoch_flag',1,[0,1],...
  'trial_info_rej_flag',1,[0,1],....
  'toilim',[],[-Inf,Inf],...  %Set value to define time window seen for given epoch, by default is set to view entire epoch
  'toilim_raw',[],[-Inf Inf],... %Set value to define time window seen for given raw epoch, by default is set to view entire raw epoch
  'raw_file',[],[],...
  'raw_flag',1,[0,1],...
  'ylim',[0 180],[-Inf,Inf],...
  'ylim_raw',[-400 400],[-Inf,Inf],...
},...
false);                                        % or calculate noise covar for each trial and then avg

if ~isfield(epoch_data.epochs,'trial_info')
    parms.trial_info_flag=0;
end;

if ~parms.trial_info_flag || ~parms.trial_info_rej_flag %if no trial info section or do not use previous reject info from trial info than modfying bad markers in epoch_data is disabled.
    parms.save_epoch_flag=0;
end;

%% DEFAULTS
global CHANNELDATA
global TRIALDATA
eeg               =  0;       % 1 for eeg data rejection, 0 for meg data rejection
rej_input = 'No';
if exist('reject_data','var') % && ~parms.reject_flag
    rej_input = questdlg('Do you want to use the previous reject data','Reject Menu');
    parms.reject_flag=1;
end;
data = epoch_data;
if ~isempty(parms.raw_file)  %if raw file is specified load it
    raw_data = load(parms.raw_file);
    raw_data = raw_data.epoch_data;
elseif(parms.raw_flag) %if not specified set epoch_data to raw_data unless flag is set to false
    raw_data=epoch_data;
end;

data         = ts_data_selection(data,'chantype',parms.chantype,'events',parms.events,'keepbadchans',parms.keepbadchans,'keepbadtrials',parms.keepbadtrials);
raw_data     = ts_data_selection(raw_data,'chantype',parms.chantype,'events',parms.events,'keepbadchans',parms.keepbadchans,'keepbadtrials',parms.keepbadtrials);

typestring = {data.sensor_info.typestring};
mag_i      = strcmp ('mag' ,typestring);
grad_i     = strncmp('grad',typestring,4);
eeg_i      = strcmp ('eeg' ,typestring);
eog_i      = strcmp ('eog' ,typestring);

if isempty(parms.rawscale_mag); parms.rawscale_mag=parms.prescale_mag; end;
if isempty(parms.rawscale_eeg); parms.rawscale_eeg=parms.prescale_eeg; end;
if isempty(parms.rawscale_eog); parms.rawscale_eog=parms.prescale_eog; end;
if isempty(parms.rawscale_grad); parms.rawscale_grad=parms.prescale_grad; end;
%% Convert to human friendly units
for cc = 1:length(data.epochs)
  data.epochs(cc).data(mag_i,:,:)  = data.epochs(cc).data(mag_i,:,:) * parms.prescale_mag;
  data.epochs(cc).data(grad_i,:,:) = data.epochs(cc).data(grad_i,:,:)* parms.prescale_grad;
  data.epochs(cc).data(eeg_i,:,:)  = data.epochs(cc).data(eeg_i,:,:) * parms.prescale_eeg;
  data.epochs(cc).data(eog_i,:,:)  = data.epochs(cc).data(eog_i,:,:) * parms.prescale_eog;
  raw_data.epochs(cc).data(mag_i,:,:)  = raw_data.epochs(cc).data(mag_i,:,:) * parms.rawscale_mag;
  raw_data.epochs(cc).data(grad_i,:,:) = raw_data.epochs(cc).data(grad_i,:,:)* parms.rawscale_grad;
  raw_data.epochs(cc).data(eeg_i,:,:)  = raw_data.epochs(cc).data(eeg_i,:,:) * parms.rawscale_eeg;
  raw_data.epochs(cc).data(eog_i,:,:)  = raw_data.epochs(cc).data(eog_i,:,:) * parms.rawscale_eog;  
end;

    if  ~isempty(parms.reject_file)
        if exist(parms.reject_file, 'file')
            load(parms.reject_file)
        else
            error('Specified Reject File does not exist');
        end;
    end;


%% Establish Data Sets
ScreenSize = get(0,'ScreenSize');
if exist('userinfo','var') && parms.usr_file
    sum_size  = userinfo{1};
    plot_size = userinfo{2};
    raw_size  = userinfo{3};
else
    sum_size   = [.5 .5 (2*(ScreenSize(4)/3)) (ScreenSize(4)/2)];
    plot_size  = [0 0 (4*(ScreenSize(4)/3)) ScreenSize(4)];
    raw_size   = [0 0 (4*(ScreenSize(4)/3)) ScreenSize(4)];
end;
sumplot=figure('Position',sum_size,'NumberTitle','off','Name','Summary su','Visible','off','WindowScrollWheelFcn',{@VR_zoom_in_out,2});


% CHANNEL DATA SET
CHANNELDATA=[];
CHANNELDATA.chan_num    = zeros(data.num_sensors,1); %Channel's Ref #
CHANNELDATA.label       =cell(data.num_sensors); %Channel label
CHANNELDATA.max(1).data = zeros(data.num_sensors,3); %Max(1) stores a absmax value for a channel based on the entire data set Val,Trial,CC
CHANNELDATA.max(2).data = zeros(data.num_sensors,3); %Max(2) stores a current absmax based on what's been rejected
CHANNELDATA.patch       = zeros(data.num_sensors,1);
CHANNELDATA.toggle      = true(data.num_sensors,1);
CHANNELDATA.panel       = uipanel('Title','Channel','FontSize',10,'BackgroundColor','white','Position',[0 .5 .75 .5],'Tag','chan_plots','BorderWidth',.2,'BorderType','line','Parent',sumplot);
CHANNELDATA.badpanel    = uitable('Units','normalized','Position',[.75 .5 .25 .5],'tag','bad_chans','Parent',sumplot);
CHANNELDATA.plot_axes   = axes('Parent',CHANNELDATA.panel);
CHANNELDATA.hg          = hggroup('Parent',CHANNELDATA.plot_axes);
CHANNELDATA.eoline      = []; %extreme outlier line
CHANNELDATA.oline       = []; %outlier line.
CHANNELDATA.plotdata    = [];
CHANNELDATA.plotfig     = figure('Position',plot_size,'NumberTitle','off','WindowScrollWheelFcn',{@VR_zoom_in_out,0},'WindowKeyPressFcn',@arrowpress);
CHANNELDATA.view_mode   = 'Good';
CHANNELDATA.sumplot     = sumplot;

for i=1:length(data.epochs)
    CHANNELDATA.plotdata(i).plot_lines      = zeros(data.num_sensors,data.epochs(i).num_trials);
    CHANNELDATA.plotdata(i).toggle          = true(data.num_sensors,data.epochs(i).num_trials);
    CHANNELDATA.plotdata(i).threshold_lines = zeros(data.num_sensors,2);
    CHANNELDATA.plotdata(i).axes            = zeros(data.num_sensors,1);
    CHANNELDATA.plotdata(i).panels          = [];
    CHANNELDATA.plotdata(i).raw_pan         = [];
    CHANNELDATA.plotdata(i).raw_axes        = zeros(data.num_sensors,1);
    CHANNELDATA.plotdata(i).raw_lines       = zeros(data.num_sensors,data.epochs(i).num_trials);
end;
CHANNELDATA.raw_plot  = figure('Position',raw_size,'NumberTitle','off','Name','Raw Plot','Visible','off','WindowScrollWheelFcn',{@VR_zoom_in_out,1});
CHANNELDATA.c_cond    = 1;
CHANNELDATA.num_conds = length(data.epochs);
CHANNELDATA.outliers  = false(data.num_sensors,2); %first dimension is extreme 2nd is regular

set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);

CHANNELDATA.text      = zeros(length(data.epochs),1);
CHANNELDATA.etext     = uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 70 500 20],'Visible','on','HorizontalAlignment','left');
CHANNELDATA.otext     = uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 50 500 20],'Visible','on','HorizontalAlignment','left');
CHANNELDATA.print_txt = uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[400 50 500 20],'Visible','on','HorizontalAlignment','left');
CHANNELDATA.sumtext   = uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 90 500 20],'Visible','on','HorizontalAlignment','left','String','Summary of Outlier Data');
CHANNELDATA.quit      = 0;
%Trial Data Set
TRIALDATA=[];
TRIALDATA.trial_num          = zeros(sum([data.epochs(:).num_trials]),1);
TRIALDATA.trial_cc           = zeros(sum([data.epochs(:).num_trials]),1);
TRIALDATA.max(1).data        = zeros(sum([data.epochs(:).num_trials]),2);
TRIALDATA.max(2).data        = zeros(sum([data.epochs(:).num_trials]),2);
TRIALDATA.patch              = zeros(sum([data.epochs(:).num_trials]),1);
TRIALDATA.toggle             = true(sum([data.epochs(:).num_trials]),1);
TRIALDATA.panel              = uipanel('Title','Trial','FontSize',10,'BackgroundColor','white','Position',[0 0 .75 .5]   ,'tag','trial_plots','BorderWidth',.5,'BorderType','line','Parent',sumplot);
TRIALDATA.badpanel           = uitable('Units','normalized','Position',[.75 0 .25 .5],'Tag','bad_trials','Parent',sumplot);
TRIALDATA.plot_axes          = axes('Parent',TRIALDATA.panel);
TRIALDATA.hg                 = hggroup('Parent',TRIALDATA.plot_axes);
TRIALDATA.eoline             = [];
TRIALDATA.oline              = [];
TRIALDATA.outliers           = false(sum([data.epochs(:).num_trials]),2); %first dimension is extreme 2nd is regular

if ~isempty(parms.toilim)
tidx = nearest(data.epochs(1).time,parms.toilim(1)):nearest(data.epochs(1).time,parms.toilim(end));
else
    tidx = 1:length(data.epochs(1).time);
end;

if parms.raw_flag
    if ~isempty(parms.toilim_raw)
        tidx_raw = nearest(raw_data.epochs(1).time,parms.toilim_raw(1)):nearest(raw_data.epochs(1).time,parms.toilim_raw(end));
    else
        tidx_raw = 1:length(raw_data.epochs(1).time);
    end;
end;

CHANNELDATA.Ylim = parms.ylim;

if parms.raw_flag
 
    CHANNELDATA.Ylim_raw = parms.ylim_raw;

end;


tic
for cc = 1:length(data.epochs)
    pages = ceil(data.num_sensors/35);
    page_count=1;
    for pindex=1:pages
        CHANNELDATA.plotdata(cc).panels(pindex)=uipanel('Parent',CHANNELDATA.plotfig,'Visible','off');
        if parms.raw_flag
            CHANNELDATA.plotdata(cc).raw_pan(pindex)=uipanel('Parent',CHANNELDATA.raw_plot,'Visible','off');
        end;
    end;
    row=0;
    column=0;
    g_count=1;
    clear i 
    i  = 1:data.num_sensors;  %Define Channel Vector
    ti = 1:data.epochs(cc).num_trials; %Defines Trial Vector
    if cc == 1
        CHANNELDATA.chan_num                                        = i ;
        CHANNELDATA.label                                           = [arrayfun(@(x) data.sensor_info(1,x).label,i,'UniformOutput',0)]';
        [CHANNELDATA.max(1).data(:,1) CHANNELDATA.max(1).data(:,2)] = arrayfun(@(x) max(max(abs(data.epochs(1).data(x,tidx,(1:data.epochs(1).num_trials))),[],2)),i);
        [CHANNELDATA.max(1).data(:,3)]                              = deal(1);
        CHANNELDATA.patch                                           = arrayfun(@(x) patch(CHANNELDATA.max(1).data(x,1),CHANNELDATA.chan_num(x),[0 0 0],'Parent',CHANNELDATA.hg,'Marker','.','EdgeColor',[0 0 0],'ButtonDownFcn',{@VR_Sum_points},'tag','chan'),i);    %Need move to outside of this loop, assign patches once do not reassign each time max exceeds
    else
        [c_max c_max_tri] = arrayfun(@(x) max(max(abs(data.epochs(cc).data(x,tidx,(1:data.epochs(cc).num_trials))),[],2)),i);
        c_max             = c_max';
        c_max_tri         = c_max_tri';
        temp              = find(c_max > CHANNELDATA.max(1).data(:,1));
        CHANNELDATA.max(1).data(temp,1) = c_max(temp);
        CHANNELDATA.max(1).data(temp,2) = c_max_tri(temp);
        CHANNELDATA.max(1).data(temp,3) = deal(cc);
    end;

    page_num = arrayfun(@(x) ceil(x/35),1:data.num_sensors); %Vector List of a channels assigned page
    col        = [.075 .2 .325 .4500 .575 .700 .825];
    row(1:7)   = .85;
    row(8:14)  = .675;
    row(15:21)  = .500;
    row(22:28) = .325;
    row(29:35) = .150;
    k          = [0 0 0];
    while length(row)<data.num_sensors || length(col)<data.num_sensors
        col = [col col];
        row = [row row];
    end;
    
        
    left(1:data.num_sensors)   = col(1:data.num_sensors);
    bottom(1:data.num_sensors) = row(1:data.num_sensors);
%     
%     
    
    CHANNELDATA.plotdata(cc).axes = arrayfun(@(x) axes('Parent',CHANNELDATA.plotdata(cc).panels(page_num(x)),'Position',[left(x) bottom(x) .1 .1],'Xlim', [data.epochs(cc).time(tidx(1)) data.epochs(cc).time(tidx(end))], 'Ylim',parms.ylim),i);
    arrayfun(@(chan) set(get(CHANNELDATA.plotdata(cc).axes(chan),'Title'),'String',data.sensor_info(1,chan).label,'ButtonDownFcn',{@VR_Channel_disable}),i);
    
    if parms.raw_flag
        CHANNELDATA.plotdata(cc).raw_axes=arrayfun(@(x) axes('Parent',CHANNELDATA.plotdata(cc).raw_pan(page_num(x)),'Position',[left(x) bottom(x) .1 .1],'Xlim', [raw_data.epochs(cc).time(tidx_raw(1)) raw_data.epochs(cc).time(tidx_raw(end))],'Ylim',parms.ylim_raw),i);
        arrayfun(@(chan) set(get(CHANNELDATA.plotdata(cc).raw_axes(chan),'Title'),'String',data.sensor_info(1,chan).label),i);
    end;
    
%     CHANNELDATA.plotdata(cc).plot_lines = arrayfun(@(i) gen_plotline(i),i,'UniformOutput',0);
%     if parms.raw_flag
%         CHANNELDATA.plotdata(cc).raw_lines  = arrayfun(@(i) gen_rawline(i),i,'UniformOutput',0);
%     end;
    
    CHANNELDATA.plotdata(cc).threshold_lines(:,1)=arrayfun(@(x) line(data.epochs(cc).time(tidx),abs(mean(data.epochs(cc).data(x,tidx,:),3))+(parms.numstd)*std(data.epochs(cc).data(x,tidx,:),0,3),'Parent',CHANNELDATA.plotdata(cc).axes(x),'LineStyle','-','LineWidth',.2,'Color',[0 1 0]),i);
    CHANNELDATA.plotdata(cc).threshold_lines(:,2)=arrayfun(@(x) line(data.epochs(cc).time(tidx),-(abs(mean(data.epochs(cc).data(x,tidx,:),3))+(parms.numstd)*std(data.epochs(cc).data(x,tidx,:),0,3)),'Parent',CHANNELDATA.plotdata(cc).axes(x),'LineStyle','-','LineWidth',.2,'Color',[0 1 0]),i);
    
    
    
    for i = 1:data.num_sensors

        
        ti = 1:data.epochs(cc).num_trials;
        
            CHANNELDATA.plotdata(cc).plot_lines(i,:) = arrayfun(@(x) line(data.epochs(cc).time(tidx),(data.epochs(cc).data(i,tidx,x)),'Parent',CHANNELDATA.plotdata(cc).axes(i),'Color',k,'ButtonDownFcn',{@VR_Plot_lines}),ti);
          
            
            if parms.raw_flag
                CHANNELDATA.plotdata(cc).raw_lines(i,:)  = arrayfun(@(x) line(raw_data.epochs(cc).time(tidx_raw),(raw_data.epochs(cc).data(i,tidx_raw,x)),'Parent',CHANNELDATA.plotdata(cc).raw_axes(i),'Color',k,'Visible','off'),ti);
            end;
        clear ti

    end;    
    
    
    
    if cc==1
        CHANNELDATA.text(cc)=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 30 500 20],'Visible','on','String','hold','HorizontalAlignment','left');
    else
        CHANNELDATA.text(cc)=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 30 500 20],'Visible','off','String','hold','HorizontalAlignment','left');
    end;
end;
toc
CHANNELDATA.max(2)=CHANNELDATA.max(1);
set(CHANNELDATA.plotdata(1).panels(1),'Visible','on');

%ADD UIcontrol buttons
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[5 5 40 18],'String','<','Callback',@prev_page,'TooltipString','Prev page');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[55 5 40 18],'String','>','Callback',@next_page,'TooltipString','Next page');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[115 5 40 18],'String','<<','Callback',@prev_cc,'TooltipString','Prev Condition');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[165 5 40 18],'String','>>','Callback',@next_cc,'TooltipString','Next Condition');
if parms.raw_flag
    uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[215 5 80 18],'String','View Raw','Callback',@view_raw,'TooltipString','See Raw Data');
else
     uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[215 5 80 18],'String','View Raw','Callback',@view_raw,'TooltipString','See Raw Data','Enable','off');
end;
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[300 5 80 18],'String','View all','Callback',@view_all,'TooltipString','View all data');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[400 5 80 18],'String','Selection Mode','Callback',@selection,'TooltipString','Enter Data Selection');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[500 5 80 18],'String','Good only','Callback',{@good_only,data,tidx},'TooltipString','View only good data');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[700 5 80 18],'String','Plot Options','CallBack',{@change_numstd,parms.numstd,data,tidx},'TooltipString','Change number of std for threshold');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[800 5 80 18],'String','Update','CallBack',{@VR_update,data,tidx},'TooltipString','Change number of std for threshold');

if parms.adv_flag
    uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[600 5 80 18],'String','Inspect','Callback',@Inspection,'TooltipString','Inspection mode');
else
    uicontrol(CHANNELDATA.plotfig,'Enable','off','units','pixels','position',[600 5 80 18],'String','Inspect','Callback',@Inspection,'TooltipString','Inspection mode');
end;

if parms.raw_flag
    uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[900 5 80 18],'String','View Channel','CallBack',{@view_channel,raw_data},'TooltipString','Change number of std for threshold');
else
    uicontrol(CHANNELDATA.plotfig,'Enable','off','units','pixels','position',[900 5 80 18],'String','View Channel','CallBack',{@view_channel,raw_data},'TooltipString','Change number of std for threshold');
end;

uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[1000 5 80 18],'String','Restore','Callback',{@restore_all},'TooltipString','Restore All');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[1100 5 80 18],'String','Exit','TooltipString','Quit and Save','Callback',{@quitting});
CHANNELDATA.save_menu   = uimenu(CHANNELDATA.plotfig,'Label','Save Menu');
                          uimenu(CHANNELDATA.save_menu,'Label','Save Rejects','Callback',@save_rejects); 
%Raw Plot Ylim buttons
if parms.raw_flag
uicontrol(CHANNELDATA.raw_plot,'Style','text','units','pixels','position',[5 5 40 18],'String','Min');                                                  %min text
CHANNELDATA.raw_min=uicontrol(CHANNELDATA.raw_plot,'Style','edit','units','pixels','position',[55 5 40 18],'Callback',{@VR_zoom_in_out,3},'String',num2str(CHANNELDATA.Ylim_raw(1)));  %min box
uicontrol(CHANNELDATA.raw_plot,'Style','text','units','pixels','position',[115 5 40 18],'String','Max');                                                  %max text
CHANNELDATA.raw_max=uicontrol(CHANNELDATA.raw_plot,'Style','edit','units','pixels','position',[165 5 40 18],'Callback',{@VR_zoom_in_out,3},'String',num2str(CHANNELDATA.Ylim_raw(2)));  %max box
end;

%toggle bad trials
mean_all_chan = mean(CHANNELDATA.max(1).data(:,1));
std_all_chan = std(CHANNELDATA.max(1).data(:,1));
CHANNELDATA.oline = line([mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan],[1 length(CHANNELDATA.chan_num)],'LineStyle','--','Parent',CHANNELDATA.plot_axes,'Color',[0 1 0]);
CHANNELDATA.eoline = line([mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan],[1 length(CHANNELDATA.chan_num)],'LineStyle','-','Parent',CHANNELDATA.plot_axes,'Color',[0 1 0]);


%Trial Data
temp = find(CHANNELDATA.toggle);
for cc = 1:length(data.epochs)
    for i=1:data.epochs(cc).num_trials
        if cc == 1
            ind=i;
            TRIALDATA.trial_num(i) = i;
            TRIALDATA.trial_cc(i,1) = i;
            TRIALDATA.trial_cc(i,2) = cc;
            [TRIALDATA.max(1).data(i,1) TRIALDATA.max(1).data(i,2)] = max(max(abs(data.epochs(1).data(temp,tidx,i)),[],2)); %max val, channel
            TRIALDATA.patch(i) = patch(TRIALDATA.trial_num(i),TRIALDATA.max(1).data(i,1),[0 0 0],'Tag' ,'trial','Parent',TRIALDATA.hg,'Marker','.','EdgeColor',[0 0 0],'ButtonDownFcn',{@VR_Sum_points});
        else
            ind=sum([data.epochs(1:(cc-1)).num_trials])+i;
            TRIALDATA.trial_num(ind) = ind;
            TRIALDATA.trial_cc(ind,1)=i;
            TRIALDATA.trial_cc(ind,2)=cc;
            [TRIALDATA.max(1).data(ind,1) TRIALDATA.max(1).data(ind,2)] = max(max(abs(data.epochs(cc).data(temp,tidx,i)),[],2)); %max val, channel
            TRIALDATA.patch(ind) = patch(TRIALDATA.trial_num(ind),TRIALDATA.max(1).data(ind,1),[0 0 0],'Tag' ,'trial','Parent',TRIALDATA.hg,'Marker','.','EdgeColor',[0 0 0],'ButtonDownFcn',{@VR_Sum_points});
        end;
    end;
end;
TRIALDATA.max(2)=TRIALDATA.max(1);
clear temp
mean_all_trial = mean(TRIALDATA.max(1).data(:,1));
std_all_trial = std(TRIALDATA.max(1).data(:,1));
TRIALDATA.oline = line([1 length(TRIALDATA.trial_num)],[mean_all_trial+2.5*std_all_trial mean_all_chan+2.5*std_all_trial],'LineStyle','--','Parent',TRIALDATA.plot_axes,'Color',[0 1 0]);
TRIALDATA.eoline = line([1 length(TRIALDATA.trial_num)],[mean_all_trial+5*std_all_trial mean_all_chan+5*std_all_trial],'LineStyle','-','Parent',TRIALDATA.plot_axes,'Color',[0 1 0]);
set(CHANNELDATA.plot_axes, 'Xlim',[0 max(CHANNELDATA.max(1).data(:,1))],'Ylim',[0 data.num_sensors]);
set(TRIALDATA.plot_axes, 'Ylim',[0 max(TRIALDATA.max(1).data(:,1))],'Xlim',[0 sum([data.epochs(:).num_trials])]);



%Setup Condition lines
for cc =  1:(length(data.epochs)-1)
    x=sum([data.epochs(1:cc).num_trials]);
    y=get(TRIALDATA.plot_axes,'Ylim');
    line([x x],[0 y(2)],'Parent',TRIALDATA.plot_axes,'HitTest','off');
end;
chan_data = find(~CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(chan_data),'ColumnName','Channel','CellSelectionCallback',{@VR_Restore,0});
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(trial_data,:),'ColumnName',{'Trial','Cond'},'CellSelectionCallback',{@VR_Restore,1},'ColumnWidth',{70 70});

disp('Accounting for all trials and channels');
VR_update([],[],data,tidx);
if exist('reject_data','var') && parms.reject_flag 
    for i=1:length(reject_data.badchans)
        x=reject_data.badchans(i);
        CHANNELDATA.toggle(x)=false;
        tmp=get(CHANNELDATA.plotdata(1).axes(x),'Title');
        VR_Channel_disable(tmp,[]);
    end;
    
    for cc=1:length(reject_data.badtrials)
        for i=1:length(reject_data.badtrials{cc})
            x = reject_data.badtrials{cc}(i); %x is relative trial number
            index=find((x==TRIALDATA.trial_cc(:,1)) & (cc==TRIALDATA.trial_cc(:,2))); %index is absolute trial number
            set(TRIALDATA.patch(index),'Visible','off','HitTest','off','Selected','off');
            TRIALDATA.toggle(index)=false;
            set(CHANNELDATA.plotdata(cc).plot_lines(:,x),'Visible','off','Selected','off','HitTest','off')
        end;
    end;
end;


if isfield(data.epochs,'trial_info') && parms.trial_info_rej_flag
    badchans = find([data.sensor_info.badchan]);
    CHANNELDATA.toggle(badchans) = false;
    tmp=get(CHANNELDATA.plotdata(1).axes(badchans),'Title');
    
    for i=1:length(tmp)
        VR_Channel_disable(tmp(i),[]);
    end;
    
    trial_info = [data.epochs.trial_info];
    trial_list = vertcat([trial_info.badtrial]);
    badtrials  = find(trial_list);
    if ~isempty(badtrials)
        set(TRIALDATA.patch(badtrials),'Visible','off','HitTest','off','Selected','off');
        TRIALDATA.toggle(badtrials)=false;
        indexed_trials = TRIALDATA.trial_cc(badtrials,:);
        ccind = unique(indexed_trials(:,2));
        for bt_ind=1:length(ccind)
            trials = find(indexed_trials(:,2)==ccind(bt_ind));
            set(CHANNELDATA.plotdata(ccind(bt_ind)).plot_lines(:,trials),'Visible','off','Selected','off','HitTest','off');
        end;
    end
end;
set(sumplot,'Visible','on');
if parms.raw_flag
set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan(1),'Visible','on');
end;
disp('Accounting for only the good data');
VR_update([],[],data,tidx);



repeatrej = 1;
while repeatrej && ishandle(CHANNELDATA.plotfig)
    refresh(CHANNELDATA.plotfig)
    if CHANNELDATA.quit == 0,
        uiwait(CHANNELDATA.plotfig);
        continue;
    else
%         exit_input = questdlg('Are you sure you want to exit','Exit Menu');
%         if strcmpi(exit_input,'Yes')
            save_rejects;
%         else
            CHANNELDATA.quit=0;
%         end;
    end;
end;
close all
        
function save_rejects(varargin)
        if ~exist('reject_data','var') || ~isfield(reject_data,'event_code')
            reject_data.event_code = [data.epochs.event_code];
        end;
        if ~isfield(reject_data,'badtrials'), reject_data.badtrials = []; end;
        
        if ~isfield(reject_data,'chanindex'), reject_data.chanindex = parms.chantype; end;
        
        if ~isfield(reject_data,'badchans'), reject_data.badchans=[]; end;
        badtrials={};
        for ind=1:CHANNELDATA.num_conds
            badtrials{ind}=[];
        end;
        reject_data.badchans = find(~CHANNELDATA.toggle);
        reject_data.badchanlabels = {data.sensor_info(reject_data.badchans).label};
        
        if parms.save_epoch_flag
            if parms.keepbadchans
                indexs = [];
                for i=1:length(reject_data.badchanlabels)
                    tmp = find(strcmp(reject_data.badchanlabels(i),{epoch_data.sensor_info.label}));
                    indexs = [indexs tmp];
                end;
                
            else
                [epoch_data.sensor_info(:).badchan]        = deal(0);
                [epoch_data.sensor_info(badchans).badchan] = deal(1);
            end;
        end;
        
        disp('Current reject data structure:');
        reject_data
        
            % bad trials
            abs_bad_trials = find(~TRIALDATA.toggle); %Listing of absolute bad trials by absolute number
            for ti = 1:length(abs_bad_trials)
                tmp=abs_bad_trials(ti);
                trial_cc=TRIALDATA.trial_cc(tmp,1);
                cc=TRIALDATA.trial_cc(tmp,2);
                badtrials{cc}=[badtrials{cc} trial_cc];
            end;
            reject_data.badtrials=badtrials;
            
           
            
            if parms.trial_info_flag
                for D_ind=1:length(data.epochs)
                    for ED_ind=1:length(epoch_data.epochs)
                        if epoch_data.epochs(ED_ind).event_code==data.epochs(D_ind).event_code
                            epoch_index(D_ind)=ED_ind;
                            break;
                        end
                    end;    
                end;
                trial_info=[];
                for bi=1:length(badtrials)
                    if ~isempty(badtrials{bi})

                        
                        trial_info(epoch_index(bi)).number        = epoch_data.epochs(epoch_index(bi)).trial_info.number(badtrials{bi});
                        trial_info(epoch_index(bi)).latency       = epoch_data.epochs(epoch_index(bi)).trial_info.latency(badtrials{bi});
                        trial_info(epoch_index(bi)).badtrial(1:length(badtrials{bi})) = 1;%epoch_data.epochs(epoch_index(bi)).trial_info.badtrial(badtrials{bi});
                        trial_info(epoch_index(bi)).event_code    = epoch_data.epochs(epoch_index(bi)).trial_info.event_code(badtrials{bi});
                        trial_info(epoch_index(bi)).duration      = epoch_data.epochs(epoch_index(bi)).trial_info.duration(badtrials{bi});
                        trial_info(epoch_index(bi)).datafile      = epoch_data.epochs(epoch_index(bi)).trial_info.datafile(badtrials{bi});
                        trial_info(epoch_index(bi)).events_fnames = epoch_data.epochs(epoch_index(bi)).trial_info.events_fnames(badtrials{bi});
                    else
                        
                        trial_info(epoch_index(bi)).number        = [];
                        trial_info(epoch_index(bi)).latency       = [];
                        trial_info(epoch_index(bi)).badtrial      = [];
                        trial_info(epoch_index(bi)).event_code    = [];
                        trial_info(epoch_index(bi)).duration      = [];
                        trial_info(epoch_index(bi)).datafile      = [];
                        trial_info(epoch_index(bi)).events_fnames = [];
                    end;
                    if parms.save_epoch_flag
                        if parms.keepbadtrials
                            indexs = arrayfun(@(x) find(trial_info(epoch_index(bi)).latency(x) == epoch_data.epochs(epoch_index(bi)).trial_info.latency),1:length(trial_info(epoch_index(bi)).latency));
                            epoch_data.epochs(epoch_index(bi)).trial_info.badtrial(indexs) = 1;
                        else
                            epoch_data.epochs(epoch_index(bi)).trial_info.badtrial(:)             = 0;
                            epoch_data.epochs(epoch_index(bi)).trial_info.badtrial(badtrials{bi}) = 1;
                        end;
                       epoch_data.epochs(epoch_index(bi)).num_rejects.manual = epoch_data.epochs(epoch_index(bi)).num_rejects.manual + length(badtrials{bi});
                    end;
                    
                end
            end;
 
            reject_data.badchans = find(~CHANNELDATA.toggle);
            reject_data.badchanlabels = {data.sensor_info(reject_data.badchans).label};
            if exist('trial_info','var')
                reject_data.badtrial_info = trial_info;
            end;
            disp('Corrected reject data structure: ')
            for i = 1: length(reject_data.event_code)
                reject_data.badtrials{i}=sort(reject_data.badtrials{i});
                Num_badtrials(i) = length(reject_data.badtrials{i});
            end;
            reject_data
            clear summary_trials
            Num_trials_orig = [data.epochs(:).num_trials];
            Num_trials_new = Num_trials_orig - Num_badtrials;
            print_out={'Please review the rejection below'; 'click yes to save'};
            print_out{end+1}=['Rejected Channels:'];
            if ~isempty(reject_data.badchanlabels)
                tmpstr = '';
                for indo=1:length(reject_data.badchanlabels) % added to put badchans in row instead of column in output 6.6.11 -BQR
                    tmpstr = [tmpstr ' ' reject_data.badchanlabels{indo}]; 
                end;
                print_out{end+1}=tmpstr;
                clear tmpstr
            else
                print_out{end+1}='No bad channels.';
            end
            print_out{end+1}=['Original number of trials: ',num2str(Num_trials_orig)];
            print_out{end+1}=['Number of good trials: ',num2str(Num_trials_new)];
            summary_trials {1} = [' ',num2str(reject_data.badtrials{1,1}),'; '];
            clear sumbadtrid*;
            sumbadtridc{1} = reject_data.badtrials{1,1};
            print_out{end+1}='Bad trials in each condition:';
            for i = 1: length(reject_data.event_code)
                print_out{end+1} = ['Condition' num2str(i) '(Event code' num2str(reject_data.event_code(i)),') ' num2str(reject_data.badtrials{1,i})];
                %print_out{end+1} = [num2str(reject_data.badtrials{1,i})];
                if i >= 2
                    summary_trials {end+1} = [' ', num2str(sum(Num_trials_orig(1:i-1))+reject_data.badtrials{1,i}),'; '];
                    sumbadtridc{i} = [sum(Num_trials_orig(1:i-1))+reject_data.badtrials{1,i}];
                end
            end
            print_out{end+1}=['Summary method bad trial index:'];
            print_out{end+1}=[summary_trials{:}];
            happy = questdlg(print_out,'Save Reject File');
            for i=3:length(print_out);
                disp(print_out{i});
            end;
            if strcmpi(happy,'Yes')
                filename = [parms.prefix , '_VR.mat'];
                rejectfilenew = fullfile(parms.rootoutdir,filename);
                disp(['Saving reject data to ',rejectfilenew])
                save(rejectfilenew,'reject_data')
                repeatrej = 0;
                if parms.usr_file
                    userinfo = {};
                    userinfo{1} = get(sumplot,'Position');
                    userinfo{2} = get(CHANNELDATA.plotfig,'Position');
                    userinfo{3} = get(CHANNELDATA.raw_plot,'Position');
                    save('~/matlab/ts_man_userinfo','userinfo');
                end;
            elseif strcmpi(happy,'No')
                disp(['Not saving Data']);
                repeatrej = 0;
            elseif strcmpi(happy,'Cancel')
               repeatrej=1;
               CHANNELDATA.quit=0;
            end;      
 end

function arrowpress(~,evnt)

if strcmpi('control',evnt.Modifier)
    if strcmpi(evnt.Key,'R');
        save_rejects;
    elseif strcmpi(evnt.Key,'A') %Autoscaling
        autoscale;
    end;
else         
    if strcmpi('uparrow',evnt.Key)
        prev_cc;
    elseif strcmpi('downarrow',evnt.Key)
        next_cc;
    elseif strcmpi('leftarrow',evnt.Key)
    	prev_page;
    elseif strcmpi('rightarrow',evnt.Key)
    	next_page;
    end;
end;
end

function autoscale(varargin)
    axes     = {};
    axes_raw = {};
    for ind=1:length(CHANNELDATA.plotdata)
        axes     = vertcat(axes, get(CHANNELDATA.plotdata(1).panels,'Children'));
        axes_raw = vertcat(axes_raw, get(CHANNELDATA.plotdata(1).panels,'Children'));
    end;
        axes     = vertcat(axes{:});
        axes_raw = vertcat(axes_raw{:});
        
        set(axes,'Ylim',parms.ylim)
        set(axes_raw,'Ylim',parms.ylim_raw);
    
end

function next_page(varargin)

h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on'); %change from 1
current=find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels==h);

if ((current+1)<=length(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels))
    new_h = CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels(current+1);
    axes = get(h,'Children');
    Ylim = get(axes(1),'Ylim');
    set(axes,'Ylim',CHANNELDATA.Ylim);
    n_axes= get(new_h,'Children');
    set(n_axes,'Ylim',Ylim);
    set(h,'Visible','off');
    set(new_h,'Visible','on');
end;

%Raw Plot
h_raw=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on'); %change from 1
current=find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan==h_raw);
if ((current+1)<=length(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan))
    new_h_raw = CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan(current+1);
    axes_raw = get(h_raw,'Children');

    set(h_raw,'Visible','off');

    set(new_h_raw,'Visible','on');
end;
uiresume(CHANNELDATA.plotfig);
end

function prev_page(varargin)
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on'); %change from 1
current=find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels==h);
if ((current-1)>=1)
    new_h = CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels(current-1);
    set(h,'Visible','off');
    set(new_h,'Visible','on');
     axes = get(h,'Children');
     Ylim_n = get(axes(1),'Ylim');
     set(axes,'Ylim',CHANNELDATA.Ylim);
     new_axes = get(new_h,'Children');
     set(new_axes,'Ylim',Ylim_n);
end;
h_raw=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on'); %change from 1
current=find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan==h_raw);
if ((current-1)>=1)
    new_h_raw = CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan(current-1);
    set(h_raw,'Visible','off');
    set(new_h_raw,'Visible','on');
    axes_raw = get(h_raw,'Children');
    Ylim_raw =  get(axes_raw(1),'Ylim');
    set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
    n_axes_raw = get(new_h_raw,'Children');
    set(n_axes_raw,'Ylim',Ylim_raw);   
end;
uiresume(CHANNELDATA.plotfig);
end

function next_cc(varargin)
if (CHANNELDATA.c_cond+1<=length(CHANNELDATA.plotdata))
    h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
    h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
    set(CHANNELDATA.text(CHANNELDATA.c_cond),'Visible','off');
    CHANNELDATA.c_cond=CHANNELDATA.c_cond+1;
    set(CHANNELDATA.text(CHANNELDATA.c_cond),'Visible','on');
    new_h=CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels(1);
    new_h_raw = CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan(1);
    set(h,'Visible','off');
    axes = get(h,'Children');

    set(new_h,'Visible','on');
    set(h_raw,'Visible','off');
    axes_raw = get(h_raw,'Children');

    set(new_h_raw,'Visible','on');
    set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
    tmp_raw_lines = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines,'Visible','on');
    if ~isempty(tmp_raw_lines)
        [raw_info(1) raw_info(2)] = find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines == tmp_raw_lines(1));
        set(CHANNELDATA.raw_plot,'Name',['Raw Plot showing Trial: ',num2str(raw_info(2)), 'in Condition Number: ',num2str(CHANNELDATA.c_cond)]);
    else
        set(CHANNELDATA.raw_plot,'Name',['Raw Plot']);
    end;
end;
uiresume(CHANNELDATA.plotfig);
end

function prev_cc(varargin)

if (CHANNELDATA.c_cond-1>=1)
    h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
    h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
    set(CHANNELDATA.text(CHANNELDATA.c_cond),'Visible','off');
    CHANNELDATA.c_cond=CHANNELDATA.c_cond-1;
    set(CHANNELDATA.text(CHANNELDATA.c_cond),'Visible','on');
    new_h=CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels(1);
    new_h_raw = CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan(1);
    %Restore zoom on panels
    axes = get(h,'Children');

    
    
    set(h,'Visible','off');
    set(new_h,'Visible','on');
    set(h_raw,'Visible','off');
    set(new_h_raw,'Visible','on');
    set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
    tmp_raw_lines = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines,'Visible','on');
    if ~isempty(tmp_raw_lines)
        [raw_info(1) raw_info(2)] = find(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_lines == tmp_raw_lines(1));
        set(CHANNELDATA.raw_plot,'Name',['Raw Plot showing Trial: ',num2str(raw_info(2)), 'in Condition Number: ',num2str(CHANNELDATA.c_cond)]);
    else
        set(CHANNELDATA.raw_plot,'Name',['Raw Plot']);
    end;
end;
uiresume(CHANNELDATA.plotfig);
end

function view_all(varargin)

CHANNELDATA.view_mode='All'; %Change view mode
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);

for cc=1:CHANNELDATA.num_conds
    tmp_lines=findobj(CHANNELDATA.plotdata(cc).plot_lines,'flat','Visible','off');
    tmp_sel_lines=findobj(CHANNELDATA.plotdata(cc).plot_lines,'flat','Color',[0 0 1],'-or','Color',[1 0 0]);
    tmp_threslines=findobj(CHANNELDATA.plotdata(cc).threshold_lines,'flat','Visible','off');
    tmp_chanpatch=findobj(CHANNELDATA.patch,'flat','Visible','off');
    tmp_trialpatch=findobj(TRIALDATA.patch,'flat','Visible','off');
    tmp_sel_chan_patch=findobj(CHANNELDATA.patch,'flat','Selected','on');
    tmp_sel_trial_patch=findobj(TRIALDATA.patch,'flat','Selected','on');
    tmp_titles=findobj('Type','text','Color',[0 0 1]);
    set(tmp_sel_chan_patch,'Visible','on','Selected','off','HitTest','on','ButtonDownFcn',{@VR_Sum_points});
    %Set Trial Patches to on
    set(tmp_sel_trial_patch,'Visible','on','Selected','off','HitTest','on','ButtonDownFcn',{@VR_Sum_points});
    %Set Trial Patches to on
    set(tmp_trialpatch,'Visible','on','Selected','off','HitTest','on','ButtonDownFcn',{@VR_Sum_points});
    set(tmp_chanpatch,'Visible','on','Selected','off','HitTest','on','ButtonDownFcn',{@VR_Sum_points});
    set(tmp_lines,'Visible','on','Selected','off','HitTest','on','ButtonDownFcn',{@VR_Plot_lines});
    set(tmp_sel_lines,'Color',[0 0 0],'Selected','off','HitTest','on','LineWidth',.5,'ButtonDownFcn',{@VR_Plot_lines});
    set(tmp_threslines,'Visible','on');
    set(tmp_titles,'Color',[0 0 0]);
end;
for i=1:length(CHANNELDATA.patch)
    set(CHANNELDATA.patch(i),'Xdata',CHANNELDATA.max(1).data(i,1));
end;
for i=1:length(TRIALDATA.patch)
    set(TRIALDATA.patch(i),'Ydata',TRIALDATA.max(1).data(i,1));
end;
%E_line and O_line CHANNELDATA.oline = line([mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan],[1 length(CHANNELDATA.chan_num)]
mean_all_trial = mean(TRIALDATA.max(1).data(:,1));
std_all_trial = std(TRIALDATA.max(1).data(:,1));
mean_all_chan = mean(CHANNELDATA.max(1).data(:,1));
std_all_chan = std(CHANNELDATA.max(1).data(:,1));
set(TRIALDATA.oline,'Xdata',[1 length(TRIALDATA.trial_num)],'Ydata',[mean_all_trial+2.5*std_all_trial mean_all_trial+2.5*std_all_trial]);
set(TRIALDATA.eoline,'Xdata',[1 length(TRIALDATA.trial_num)],'Ydata',[mean_all_trial+5*std_all_trial mean_all_trial+5*std_all_trial]);
set(CHANNELDATA.oline,'Xdata',[mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan],'Ydata', [1 length(CHANNELDATA.chan_num)]);
set(CHANNELDATA.eoline,'Xdata',[mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan],'Ydata', [1 length(CHANNELDATA.chan_num)]);
set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
uiresume(CHANNELDATA.plotfig);
end

function selection(varargin)

chan_bad=find(~CHANNELDATA.toggle(:));
trial_bad=find(~TRIALDATA.toggle(:));
view_all();
CHANNELDATA.view_mode='Selection Mode';
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);

for i=1:length(CHANNELDATA.patch)
    set(CHANNELDATA.patch(i),'Xdata',CHANNELDATA.max(1).data(i,1),'EdgeColor',[0 0 0]);
end;

for i=1:length(TRIALDATA.patch)
    set(TRIALDATA.patch(i),'Ydata',TRIALDATA.max(1).data(i,1),'EdgeColor',[0 0 0]);
end;

for i=1:length(chan_bad)
    tmp_index=chan_bad(i); %index of bad channel
    for cc=1:CHANNELDATA.num_conds
        tmp = get(CHANNELDATA.plotdata(cc).axes(tmp_index),'Title');
        tmp_patch = CHANNELDATA.patch(tmp_index);
        set(tmp,'Color',[0 0 1]);
        set(tmp_patch,'EdgeColor',[0 0 1],'HitTest','on');
    end;
    
    if length(varargin)<1
        VR_Sum_points(CHANNELDATA.patch(tmp_index),[]);
    end;

end;

for id=1:length(trial_bad)
    tmp_index=trial_bad(id); %absolute trial number
    tmp_trial=TRIALDATA.trial_cc(tmp_index,1);
    cc=TRIALDATA.trial_cc(tmp_index,2);
    set(CHANNELDATA.plotdata(cc).plot_lines(:,tmp_trial),'Color',[0 0 1],'HitTest','on','LineWidth',2,'ButtonDownFcn',{@VR_Plot_lines});
    VR_Sum_points(TRIALDATA.patch(tmp_index),[]);
    set(TRIALDATA.patch(tmp_index),'EdgeColor',[0 0 1],'HitTest','on','Selected','off');
    set(findobj(CHANNELDATA.plotdata(cc).plot_lines,'Color',[1 0 0]),'Color',[0 0 1]);
end;



set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
uiresume(CHANNELDATA.plotfig);
end

function Inspection(varargin)
selection(1);
CHANNELDATA.view_mode='Inspection';
set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
set(CHANNELDATA.patch,'ButtonDownFcn',{@VR_Sum_points});
set(TRIALDATA.patch,'ButtonDownFcn',{@VR_Sum_points});
uiresume;
end

function good_only(varargin)

CHANNELDATA.view_mode='Good';
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines,'Color',[0 0 0]);
chan_bad=find(~CHANNELDATA.toggle(:));
trial_bad=find(~TRIALDATA.toggle(:));
for i=1:length(chan_bad)
    tmp_index=chan_bad(i);
    tmp=get(CHANNELDATA.plotdata(CHANNELDATA.c_cond).axes(tmp_index),'Title');
    set(tmp,'Color',[0 0 1]);
    set(CHANNELDATA.patch(tmp_index),'Selected','off','HitTest','off','Visible','off');
    for ti=1:length(CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines(tmp_index,:))
        x=CHANNELDATA.plotdata(CHANNELDATA.c_cond).plot_lines(tmp_index,ti);
        set(x,'HitTest','off','Visible','off');
    end;
    set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).threshold_lines(tmp_index,1),'Visible','off');
    set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).threshold_lines(tmp_index,2),'Visible','off');
end;
for i=1:length(trial_bad)
    x=trial_bad(i);
    tmp_index=TRIALDATA.trial_cc(x,1);
    cc=TRIALDATA.trial_cc(x,2);
    set(TRIALDATA.patch(x),'Selected','off','HitTest','off','Visible','off','EdgeColor',[0 0 0]);
    h=CHANNELDATA.plotdata(cc).plot_lines(:,tmp_index);
    set(h,'HitTest','off','Visible','off');
end;
for i=1:length(CHANNELDATA.patch)
    set(CHANNELDATA.patch(i),'Xdata',CHANNELDATA.max(2).data(i,1),'EdgeColor',[0 0 0]);
end;
for i=1:length(TRIALDATA.patch)
    set(TRIALDATA.patch(i),'Ydata',TRIALDATA.max(2).data(i,1),'EdgeColor',[0 0 0]);
end;
set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
VR_update(0,[],varargin{3},varargin{4});
uiresume(CHANNELDATA.plotfig);
end

function view_channel(h,x,raw_data)

    
prompt = {'Channel to view','Event Code'};
def={data.sensor_info(1).label,num2str(data.epochs(1).event_code)};
dlg_title = 'Plot channel';
num_lines = 1;
chan = inputdlg(prompt,dlg_title,num_lines,def);

cond_ind=0; %Determine which epoch contains that event
for ci=1:length(raw_data.epochs)
    if str2double(chan{2})==raw_data.epochs(ci).event_code
        cond_ind = ci;
        break;
    end;
end;
        
if cond_ind==0
    error('Event was not found')
end;
    
e_code = chan(2);
chan = chan(1);
chan_num = find(strcmpi(CHANNELDATA.label,chan));
trials = find(TRIALDATA.toggle & TRIALDATA.trial_cc(:,2)==cond_ind);
trials = TRIALDATA.trial_cc(trials);
dat=raw_data.epochs(cond_ind).data(chan_num,:,trials);
dat=mean(dat,3);
fig = figure; %('NumberTitle','off','Name',['Showing ' chan, 'on condition : ',num2str(cc)]);
name_txt=['Average Raw Data for ',chan{1} ,' (Event code ' num2str(e_code{1}),')'];
set(fig,'NumberTitle','off','Name',name_txt);
ax = axes('Parent',fig);
line(raw_data.epochs(cond_ind).time,dat,'Parent',ax);
end

    




function change_numstd(h,x,numstd,data,tidx)

prompt = {'Number of Std for Threshold:','New Lower Y-limit','New Upper Y-limit'};
dlg_title = 'Change Plot Properties';
num_lines = 1;
def = {num2str(numstd),num2str(CHANNELDATA.Ylim(1)),num2str(CHANNELDATA.Ylim(2))};
values = cell2num(inputdlg(prompt,dlg_title,num_lines,def));
if isempty(values)
  return
end
userstd=values(1);
Ylim = [values(2) values(3)];
for cc=1:length(CHANNELDATA.plotdata)
    for ci=1:length(CHANNELDATA.plotdata(cc).threshold_lines(:,1))
        set(CHANNELDATA.plotdata(cc).threshold_lines(ci,1),'Ydata',abs(mean(data.epochs(cc).data(ci,tidx,:),3))+userstd*std(data.epochs(cc).data(ci,tidx,:),0,3));
        set(CHANNELDATA.plotdata(cc).threshold_lines(ci,2),'Ydata',-(abs(mean(data.epochs(cc).data(ci,tidx,:),3))+userstd*std(data.epochs(cc).data(ci,tidx,:),0,3)));
    end;
end;
    CHANNELDATA.Ylim = Ylim;
    for cc=1:CHANNELDATA.num_conds
    set(CHANNELDATA.plotdata(cc).axes,'Ylim',CHANNELDATA.Ylim);
    end;
   
    
set(h,'CallBack',{@change_numstd,userstd,data,tidx});
end

function quitting(varargin)

h=varargin{1};
CHANNELDATA.quit=1;
set(h,'UserData','done');
uiresume(CHANNELDATA.plotfig);
end

function restore_all(varargin)

c_cond=CHANNELDATA.c_cond;
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);

CHANNELDATA.toggle(:)=true;
TRIALDATA.toggle(:)=true;
for cc=1:CHANNELDATA.num_conds
    CHANNELDATA.c_cond=cc;
    view_all();
    CHANNELDATA.view_mode='Good';
end;
CHANNELDATA.c_cond=c_cond;
chan_data = find(~CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(chan_data),'ColumnName','Channel','CellSelectionCallback',{@VR_Restore,0});
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(trial_data,:),'ColumnName',['Trial' 'Condition'],'CellSelectionCallback',{@VR_Restore,1});
set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);

%if keepbadtrials restore matrix
if parms.keepbadtrials
    if ~isempty(parms.events)
        eve_index=[];
        eve_index=arrayfun(@(x) find(parms.events(x)==[epoch_data.epochs.event_code]),1:length(parms.events));
    else
        eve_index=1:length(epoch_data.epochs)
    end;
    for idx = eve_index
        if isfield(epoch_data.epochs,'trial_info');
            epoch_data.epochs(idx).trialinfo.badtrial(find(epoch_data.epochs(idx).trial_info.badtrial))=0;
        end;
        
        epoch_data.epochs(idx).num_rejects.mag    = 0; 
        epoch_data.epochs(idx).num_rejects.grad   = 0;
        epoch_data.epochs(idx).num_rejects.eeg    = 0;
        epoch_data.epochs(idx).num_rejects.eog    = 0;
        epoch_data.epochs(idx).num_rejects.manual = 0;
        epoch_data.epochs(idx).num_rejects.skip   = 0;
    end;
end;
    
uiresume(CHANNELDATA.plotfig);
end

function view_raw(h,x)
set(CHANNELDATA.raw_plot,'Visible','on');
set(h,'Callback',@hide_raw,'String','Hide Raw');
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
end

function hide_raw(h,x)
%Restore zoom on panels
h_pan=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h_pan,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
set(CHANNELDATA.raw_plot,'Visible','off');
set(h,'Callback',@view_raw,'String','Show Raw');
end


end


