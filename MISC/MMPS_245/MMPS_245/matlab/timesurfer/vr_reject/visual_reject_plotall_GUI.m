%Last modified July 28th 2010 5 PM by ADS Version 1.0.0 
% for major support issues contract at aschulma@ucsd.edu
% This version adds in the ability to also print out a text file for
% behavioral analysis in the save data dialog. It also adds in more
% information about whats selected, and selects the chan/trial related to
% the trial/chan absmax.
function visual_reject_plotall_GUI
close all; clear all

%% SET YOUR SUBJECT AND BAND:

prompt={'Please Enter a Subject name example: ''S28a''','Enter 1 if ICA is complete','Please Enter a Band: 4-8Hz, 8-12Hz, 12-25Hz,25-55Hz'};
def={'AM_1a_test','0','4-8Hz'};
if exist(['~/matlab/VR_gui_userinfo.mat'], 'file')
    load ('~/matlab/VR_gui_userinfo.mat');
    def = userinfo{4};
end


subjinfo=inputdlg(prompt,'Load Info',1,def);
subj        = subjinfo{1};  % example: 'S28a'
ICA         =  str2double(subjinfo{2});     % 1 if ICA was done or 0 if ICA was not done for this subject
band        = subjinfo{3}; %'4-8Hz', '8-12Hz', '12-25Hz' or '25-55Hz


% %
%  subj        = 'S53j';  % example: 'S28a'
%  ICA         =  1;      % 1 if ICA was done or 0 if ICA was not done for this subject
%  band        = '8-12Hz'; %'4-8Hz', '8-12Hz', '12-25Hz' or '25-55Hz'


%% DEFAULTS
global CHANNELDATA
global TRIALDATA
eeg               =  0;       % 1 for eeg data rejection, 0 for meg data rejection
bandavg       =  1;       % 1 for displaying the band averages, or 0 for displaying raw data
numstd         =  6;       % controls the threshold (dotted black line) for bad trial identification, 5 standard deviations has worked the best so far for the old SL SNG wavelet results (1sec long wavelet); 6 works better for megshort results
ylimuser       =  [];      % changes the axis scaling, [] to get limits from the data, or enter your scale here [-8e-23 8e-23] for theta and alpha, [-2e-23 2e-23] for beta
threshuser    =  [];      % changes threshold to set values, [] to get it from the data using numstd, or enter your choice here [0 3e-23]. These numbers will be added to the mean
numstduser  =  [];
evnts            =  [1 2 3 5 6 7 11 12];      % condition numbers to plot (exmpl [1 2] or [] to plot all)
rej_input = 'No';
UNITS     = 10^24;      % convert data from [T/m]^2/Hz to 100*[fT/cm]^2/Hz
UNITS_RAW  = 10^13;     % convert data from  T/m to fT/cm
% Data root directory
%root = '/data3/kmdev/SNg'; % ip101
%root = '/local_mount/space/svarog/3/users/sarah/SNg'; %ip50 SL
%root = '/data/kmdev/SNg'; % ip82
%root = '/mnt/101data3/kmdev/SNg'; % ip32
%root = '/mnt/82data/kmdev/SNg'; % ip32 or ip50
%root = '/local_mount/space/svarog/3/users/aschulman/dev/EF' %AS Development folder
%root = '/local_mount/space/svarog/3/users/sarah/EF'; % ip50 SL EF
%root = '/mnt/82data/kmdev/Cog';  % ip50 SL Simon tasks
root = '/home/aschulman/aschulman/dummy_data/VR'; %ADS dummy data folder
% root = '/space/mdkm1/9/kmdev/projects/AM/test'; %data from AM test


if ICA == 1, icastr = '_ICA'; else icastr = ''; end;
if bandavg == 1
    toilim  = [-.5 .5];  %time window of interest (in sec)
    toilim_raw = [-.75 .75];
    disp(['Loading data for subject ',subj, icastr,' ',band,' band...'])
    if eeg == 0
        
        %         %SL ip50
        %         toilim  = [-.45 .75];
        %         load ([root,'/',subj,'/TF/matfiles/timefreq',icastr,'/power_',band,'/',subj,'_SNg',icastr,'_timefreq_band_average_freq',band,'_conds1_5_conds_1_2_3_5_avg.mat'])
        %         rejectfile = [root,'/',subj,'/TF/matfiles/timefreq',icastr,'/power_',band,'/',subj,'_SNg',icastr,'_',band,'_tf_reject_data.mat'];
        %         rejectfilenew = [root,'/',subj,'/TF/matfiles/timefreq',icastr,'/power_',band,'/',subj,'_SNg',icastr,'_',band,'_tf_reject_data_new.mat'];
        %         %SL new (only have done this for S29, S34, S31, and S33)
        %         toilim  = [-.45 .75];
        %         load ([root,'/',subj,'/TF/matfiles/timefreq',icastr,'_megshort/power_',band,'/',subj,'_SNg',icastr,'_average_freqband_',band,'_toi-0.45-0.75sec_events1_2_3_5.mat'])
        %         rejectfile = [root,'/',subj,'/TF/matfiles/',subj,'_SNg_tf_reject_data.mat'];
        %         rejectfilenew = [root,'/',subj,'/TF/matfiles/',subj,'_SNg_tf_reject_data_new.mat'];

                %SL EF
                toilim  = [-.45 .6]; toilim_raw = [-.7 .85]; %toilim_raw = [-.45 .6];  
                load ([root,'/',subj,'/TF/matfiles/timefreq',icastr,'/power_',band,'/',subj,'_EF',icastr,'_average_freqband_',band,'_events2_3_4.mat'])
                if ICA == 1
                    raw=load ([root,'/',subj,'/TF/matfiles/',subj,'_EF_ICA_epoch_data_ICA.mat']);
                    raw_data=raw.epoch_data;
                else
                    raw=load ([root,'/',subj,'/TF/matfiles/',subj,'_EF_epoch_data_1.mat']);
                    raw_data=raw.epoch_data;
                end;
                clear raw
                rejectfile = [root,'/',subj,'/TF/matfiles/',subj,'_EF',icastr,'_tf_reject_data.mat'];
                rejectfilenew = [root,'/',subj,'/TF/matfiles/',subj,'_EF',icastr,'_tf_reject_data_new.mat'];
%         
%         % SL Simon tasks
%         task = 'SC';
%         toilim  = [-.25 .75]; toilim_raw = [-.25 .75];  %toilim_raw = [-.7 .85];
%         load ([root,'/',subj,'/',subj, '_',task,'/TF/matfiles/timefreq/power_',band,'/',subj,'_',task,'_average_freqband_',band,'_toi-0.25-0.75sec_events21_22.mat']);
%         if ICA == 1
%             raw=load ([root,'/',subj,'/',subj,'_',task,'/TF/matfiles/',subj,'_',task,'_epoch_data_ICA.mat']);
%             raw_data=raw.epoch_data;
%         else
%             raw=load ([root,'/',subj,'/',subj,'_',task,'/TF/matfiles/',subj,'_',task,'_epoch_data_1.mat']);
%             raw_data=raw.epoch_data;
%         end;
%         clear raw
%         rejectfile = [root,'/',subj,'/',subj,'_',task,'/TF/matfiles/',subj,'_',task,'_tf_reject_data.mat'];
%         rejectfilenew =[root,'/',subj,'/',subj,'_',task,'/TF/matfiles/',subj,'_',task,'_tf_reject_data_new.mat'];
%         


        if exist(rejectfilenew, 'file')
            load(rejectfilenew)
        elseif exist(rejectfile, 'file')
            load(rejectfile)
        end;
    else
        fpath  = [root,'/', subj, '/TF/matfiles/timefreq_eeg_megshort/power_',band];
        fname  = [subj, '_SNg_megshort_eeg_average_freqband_', band, '_toi-0.45-0.75sec_events1_2_3_5.mat'];
        load (fullfile(fpath,fname))
        epoch_data = ts_data_selection(epoch_data,'chanlabel',{'EEG 001' 'EEG 002' 'EEG 003' 'EEG 004'},...
            'removebadchans',1);
        rejectfilenew = fullfile(fpath,[subj, '_SNg_eeg_', band,'_tf_reject_data.mat']);
        if exist(rejectfilenew, 'file')
            load(rejectfilenew)
        end;
    end;
    if exist('reject_data','var')
        disp('Most current reject data:')
        reject_data
    end
else
  % EF TASK
%     toilim = [-.5 .5];
%     toilim_raw  = [-.75 .75];
%     if  ICA == 1
%         disp(['Loading data for subject ',subj,': ', subj,'_EF_epoch_data_ICA.mat ...'])
%         load ([root,'/',subj,'/TF/matfiles/',subj,'_EF_epoch_data_ICA.mat'])
%     else
%         disp(['Loading data for subject ',subj,': ', subj,'_EF_epoch_data_1.mat ...'])
%         load ([root,'/',subj,'/TF/matfiles/',subj,'_EF_epoch_data_1.mat'])
%     end
%     if eeg == 1
%         epoch_data = ts_data_selection(epoch_data,'chanlabel',{'EEG 001' 'EEG 002' 'EEG 003' 'EEG 004'},...
%             'removebadchans',1);
%         %         epoch_data = ts_preproc (epoch_data, 'bandpass_flag',0,'bandpass_detrend_flag',1, ...
%         %             'bandpass_low_cf',4,'bandpass_low_tb',0,'bandpass_high_cf',12,'bandpass_high_tb',0, ...
%         %             'detrend_flag', 1);
%         rejectfilenew = [root,'/',subj,'/TF/matfiles/',subj,'_SNg_eeg_tf_reject_data.mat'];
%     else
%         rejectfilenew = [root,'/',subj,'/TF/matfiles/',subj,'_SNg',icastr,'_tf_reject_data.mat'];
%     end;
    
    
%         % AM task %SK 10/06/2010
%         UNITS     = 10^13;      % convert data from [T/m]^2/Hz to 100*[fT/cm]^2/Hzvisual_reject_plotall_GUI.m
%         UNITS_RAW  = 10^13;     % convert data from  T/m to fT/cm
%         task = 'test';
%         toilim  = [-.3 .8]; toilim_raw = [-.3 .8];  % ONLY TD data available
%         load ([root,'/',subj,'_test/TD/matfiles/',subj,'_',task,'_epoch_data_1.mat']);
%         raw_data = epoch_data;
%         
%         rejectfile = [root,'/',subj,'_test/TD/matfiles/',subj,'_',task,'_tf_reject_data.mat'];
%         rejectfilenew =[root,'/',subj,'_test/TD/matfiles/',subj,'_',task,'_tf_reject_data_new.mat'];
%         eventsfnames = [root,'/',subj,'_test/', subj,'_all_events_fnames.txt'];
%         
%     if exist(rejectfilenew, 'file')
%             load(rejectfilenew)
%     elseif exist(rejectfile, 'file')
%             load(rejectfile)
%     end;
%         
    if exist('reject_data','var')
        disp('Most current reject data:')
        reject_data
        
    end
end
if exist('reject_data','var')
    rej_input = questdlg('Do you want to use the previous reject data','Reject Menu');
end;
%% Convert to human friendly units
for cc = 1:length(epoch_data.epochs)
    epoch_data.epochs(cc).data = epoch_data.epochs(cc).data* UNITS;
    raw_data.epochs(cc).data = raw_data.epochs(cc).data * UNITS_RAW;
end;


%% Establish Data Sets
ScreenSize = get(0,'ScreenSize');
if exist('userinfo','var')
    sum_size = userinfo{1};
    plot_size = userinfo{2};
    raw_size = userinfo{3};
else
    sum_size = [.5 .5 (2*(ScreenSize(4)/3)) (ScreenSize(4)/2)];
    plot_size = [0 0 (4*(ScreenSize(4)/3)) ScreenSize(4)];
    raw_size = [0 0 (4*(ScreenSize(4)/3)) ScreenSize(4)];
end;
sumplot=figure('Position',sum_size,'NumberTitle','off','Name','Summary Plot','Visible','off','WindowScrollWheelFcn',{@VR_zoom_in_out,2});
if eeg==0
epoch_data = ts_data_selection(epoch_data,'chantype','grad','removebadchans',1,'events',evnts);
raw_data     = ts_data_selection(raw_data,'chantype','grad','removebadchans',1,'events',evnts);
else
epoch_data = ts_data_selection(epoch_data,'events',evnts);
raw_data     = ts_data_selection(raw_data,'events',evnts);
end;

% CHANNEL DATA SET
CHANNELDATA=[];
CHANNELDATA.chan_num = zeros(epoch_data.num_sensors,1); %Channel's Ref #
CHANNELDATA.label={}; %Channel label
CHANNELDATA.max(1).data = zeros(epoch_data.num_sensors,3); %Max(1) stores a absmax value for a channel based on the entire data set Val,Trial,CC
CHANNELDATA.max(2).data = zeros(epoch_data.num_sensors,3); %Max(2) stores a current absmax based on what's been rejected
CHANNELDATA.patch = zeros(epoch_data.num_sensors,1);
CHANNELDATA.toggle = logical(ones(epoch_data.num_sensors,1));
CHANNELDATA.panel = uipanel('Title','Channel','FontSize',10,'BackgroundColor','white','Position',[0 .5 .75 .5],'Tag','chan_plots','BorderWidth',.2,'BorderType','line','Parent',sumplot);
CHANNELDATA.badpanel = uitable('Units','normalized','Position',[.75 .5 .25 .5],'tag','bad_chans','Parent',sumplot);
CHANNELDATA.plot_axes = axes('Parent',CHANNELDATA.panel);
CHANNELDATA.hg = hggroup('Parent',CHANNELDATA.plot_axes);
CHANNELDATA.eoline = []; %extreme outlier line
CHANNELDATA.oline = []; %outlier line.
CHANNELDATA.plotdata = [];
CHANNELDATA.plotfig=figure('Position',plot_size,'NumberTitle','off','WindowScrollWheelFcn',{@VR_zoom_in_out,0});
CHANNELDATA.view_mode = 'Good';
CHANNELDATA.sumplot=sumplot;
for i=1:length(epoch_data.epochs)
    CHANNELDATA.plotdata(i).plot_lines = zeros(epoch_data.num_sensors,epoch_data.epochs(i).num_trials);
    CHANNELDATA.plotdata(i).toggle = logical(ones(epoch_data.num_sensors,epoch_data.epochs(i).num_trials));
    CHANNELDATA.plotdata(i).threshold_lines = zeros(epoch_data.num_sensors,2);
    CHANNELDATA.plotdata(i).axes = zeros(epoch_data.num_sensors,1);
    CHANNELDATA.plotdata(i).panels = [];
    CHANNELDATA.plotdata(i).raw_pan = [];
    CHANNELDATA.plotdata(i).raw_axes = zeros(epoch_data.num_sensors,1);
    CHANNELDATA.plotdata(i).raw_lines = zeros(epoch_data.num_sensors,epoch_data.epochs(i).num_trials);
end;
CHANNELDATA.raw_plot = figure('Position',raw_size,'NumberTitle','off','Name','Raw Plot','Visible','off','WindowScrollWheelFcn',{@VR_zoom_in_out,1});
CHANNELDATA.c_cond = 1;
CHANNELDATA.num_conds = length(epoch_data.epochs);
CHANNELDATA.outliers = logical(zeros(epoch_data.num_sensors,2)); %first dimension is extreme 2nd is regular
set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);

CHANNELDATA.text = zeros(length(epoch_data.epochs),1);
CHANNELDATA.etext=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 70 500 20],'Visible','on','HorizontalAlignment','left');
CHANNELDATA.otext=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 50 500 20],'Visible','on','HorizontalAlignment','left');
CHANNELDATA.print_txt=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[400 50 500 20],'Visible','on','HorizontalAlignment','left');
CHANNELDATA.sumtext=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 90 500 20],'Visible','on','HorizontalAlignment','left','String','Summary of Outlier Data');
CHANNELDATA.quit=0;
%Trial Data Set
TRIALDATA=[];
TRIALDATA.trial_num = zeros(sum([epoch_data.epochs(:).num_trials]),1);
TRIALDATA.trial_trial_datacc = zeros(sum([epoch_data.epochs(:).num_trials]),1);
TRIALDATA.max(1).data = zeros(sum([epoch_data.epochs(:).num_trials]),2);
TRIALDATA.max(2).data = zeros(sum([epoch_data.epochs(:).num_trials]),2);
TRIALDATA.patch = zeros(sum([epoch_data.epochs(:).num_trials]),1);
TRIALDATA.toggle = logical(ones(sum([epoch_data.epochs(:).num_trials]),1));
TRIALDATA.panel = uipanel('Title','Trial','FontSize',10,'BackgroundColor','white','Position',[0 0 .75 .5]   ,'tag','trial_plots','BorderWidth',.5,'BorderType','line','Parent',sumplot);
TRIALDATA.badpanel = uitable('Units','normalized','Position',[.75 0 .25 .5],'Tag','bad_trials','Parent',sumplot);
TRIALDATA.plot_axes = axes('Parent',TRIALDATA.panel);
TRIALDATA.hg=hggroup('Parent',TRIALDATA.plot_axes);
TRIALDATA.eoline = [];
TRIALDATA.oline = [];
TRIALDATA.outliers=logical(zeros(sum([epoch_data.epochs(:).num_trials]),2)); %first dimension is extreme 2nd is regular




tidx = nearest(epoch_data.epochs(1).time,toilim(1)):nearest(epoch_data.epochs(1).time,toilim(end));
tidx_raw = nearest(raw_data.epochs(1).time,toilim_raw(1)):nearest(raw_data.epochs(1).time,toilim_raw(end));
% %remove from here after test
% goodchans=find(CHANNELDATA.toggle);
% for cc = 1:length(epoch_data.epochs)
%     if min(min(min(epoch_data.epochs(cc).data(goodchans,tidx,:)))) < 0
%         ylim(cc,:) = [min(min(min(epoch_data.epochs(cc).data(goodchans,tidx,:)))) max(max(max(epoch_data.epochs(cc).data(goodchans,tidx,:))))];
%     else
%         ylim(cc,:) = [0 80*mean(mean(mean(epoch_data.epochs(cc).data(goodchans,tidx,:))))];
%     end
% end
% ylimdata = [min(ylim(:,1)) max(ylim(:,2))];

%Change here to change Ylimits and bands
if strcmpi(band,'4-8Hz')
    ylimdata = [0 180];
elseif strcmpi(band,'8-12Hz')
    ylimdata = [0 180];
elseif strcmpi(bantrial_datad,'12-25Hz')
    ylimdata = [0 60];
elseif strcmpi(band,'25-55Hz')
    ylimdata = [0 30];
end;
ylimdata_raw = [-400 400];


CHANNELDATA.Ylim = ylimdata;
CHANNELDATA.Ylim_raw = ylimdata_raw;

for cc = 1:length(epoch_data.epochs)
    pages = ceil(epoch_data.num_sensors/35);
    page_count=1;
    for pindex=1:pages
        CHANNELDATA.plotdata(cc).panels(pindex)=uipanel('Parent',CHANNELDATA.plotfig,'Visible','off');
        CHANNELDATA.plotdata(cc).raw_pan(pindex)=uipanel('Parent',CHANNELDATA.raw_plot,'Visible','off');
    end;
    row=0;
    column=0;
    g_count=1;
    for i = 1:epoch_data.num_sensors
        k=[0 0 0];
        if cc == 1
            CHANNELDATA.chan_num(i) = i;
            CHANNELDATA.label{i,1} = epoch_data.sensor_info(1,i).label;
            [CHANNELDATA.max(1).data(i,1) CHANNELDATA.max(1).data(i,2)] = max(max(abs(epoch_data.epochs(1).data(i,tidx,(1:epoch_data.epochs(1).num_trials))),[],2));
            CHANNELDATA.max(1).data(i,3) = 1;
            CHANNELDATA.patch(i) = patch(CHANNELDATA.max(1).data(i,1),CHANNELDATA.chan_num(i),[0 0 0],'Parent',CHANNELDATA.hg,'Marker','.','EdgeColor',[0 0 0],'ButtonDownFcn',{@VR_Sum_points},'tag','chan');
        else
            [x y] = max(max(abs(epoch_data.epochs(cc).data(i,tidx,(1:epoch_data.epochs(cc).num_trials))),[],2));
            if x > CHANNELDATA.max(1).data(i,1)
                CHANNELDATA.max(1).data(i,1) = x;
                CHANNELDATA.max(1).data(i,2) = y;
                CHANNELDATA.max(1).data(i,3) = cc;
                h=CHANNELDATA.patch(i);
                set(h,'XData',CHANNELDATA.max(1).data(i,1));
                set(h,'ButtonDownFcn',{@VR_Sum_points});
            end;
        end;
        
        if(column==7) %change here for number of columns on plot page
            column = 0;
            row = row+1;
        end;
        
        left = .075+(column*.125);
        bottom = .85-(.175*row);
        %Main Plot
        CHANNELDATA.plotdata(cc).axes(i)=axes('Parent',CHANNELDATA.plotdata(cc).panels(page_count),'Position',[left bottom .1 .1],'Xlim', [epoch_data.epochs(cc).time(tidx(1)) epoch_data.epochs(cc).time(tidx(end))], 'Ylim',ylimdata);
        set(get(CHANNELDATA.plotdata(cc).axes(i),'Title'),'String',epoch_data.sensor_info(1,i).label,'ButtonDownFcn',{@VR_Channel_disable});
        %Raw Plot
        CHANNELDATA.plotdata(cc).raw_axes(i)=axes('Parent',CHANNELDATA.plotdata(cc).raw_pan(page_count),'Position',[left bottom .1 .1],'Xlim', [raw_data.epochs(cc).time(tidx_raw(1)) raw_data.epochs(cc).time(tidx_raw(end))],'Ylim',ylimdata_raw);
        set(get(CHANNELDATA.plotdata(cc).raw_axes(i),'Title'),'String',epoch_data.sensor_info(1,i).label);
        
        
        for ti = 1:epoch_data.epochs(cc).num_trials
            CHANNELDATA.plotdata(cc).plot_lines(i,ti)=line(epoch_data.epochs(cc).time(tidx),(epoch_data.epochs(cc).data(i,tidx,ti)),'Parent',CHANNELDATA.plotdata(cc).axes(i),'Color',k,'ButtonDownFcn',{@VR_Plot_lines});
            CHANNELDATA.plotdata(cc).raw_lines(i,ti)=line(raw_data.epochs(cc).time(tidx_raw),(raw_data.epochs(cc).data(i,tidx_raw,ti)),'Parent',CHANNELDATA.plotdata(cc).raw_axes(i),'Color',k,'Visible','off');
        end;
        CHANNELDATA.plotdata(cc).threshold_lines(i,1)=line(epoch_data.epochs(cc).time(tidx),abs(mean(epoch_data.epochs(cc).data(i,tidx,:),3))+numstd*std(epoch_data.epochs(cc).data(i,tidx,:),0,3),'Parent',CHANNELDATA.plotdata(cc).axes(i),'LineStyle','-','LineWidth',.2,'Color',[0 1 0]);
        CHANNELDATA.plotdata(cc).threshold_lines(i,2)=line(epoch_data.epochs(cc).time(tidx),-(abs(mean(epoch_data.epochs(cc).data(i,tidx,:),3))+numstd*std(epoch_data.epochs(cc).data(i,tidx,:),0,3)),'Parent',CHANNELDATA.plotdata(cc).axes(i),'LineStyle','-','LineWidth',.2,'Color',[0 1 0]);
        column=column+1;
        
        g_count=g_count+1;
        if(g_count==36)
            g_count=1;
            page_count=page_count+1;
            row=0;
            column=0;
        end
    end;
    if cc==1
        CHANNELDATA.text(cc)=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 30 500 20],'Visible','on','String','hold','HorizontalAlignment','left');
    else
        CHANNELDATA.text(cc)=uicontrol(CHANNELDATA.plotfig,'units','pixels','Style','text','position',[100 30 500 20],'Visible','off','String','hold','HorizontalAlignment','left');
    end;
end;
CHANNELDATA.max(2)=CHANNELDATA.max(1);
set(CHANNELDATA.plotdata(1).panels(1),'Visible','on');

%ADD UIcontrol buttons
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[5 5 40 18],'String','<','Callback',@prev_page,'TooltipString','Prev page');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[55 5 40 18],'String','>','Callback',@next_page,'TooltipString','Next page');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[115 5 40 18],'String','<<','Callback',@prev_cc,'TooltipString','Prev Condition');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[165 5 40 18],'String','>>','Callback',@next_cc,'TooltipString','Next Condition');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[215 5 80 18],'String','View Raw','Callback',@view_raw,'TooltipString','See Raw Data');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[300 5 80 18],'String','View all','Callback',@view_all,'TooltipString','View all data');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[400 5 80 18],'String','Selection Mode','Callback',@selection,'TooltipString','Enter Data Selection');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[500 5 80 18],'String','Good only','Callback',{@good_only,epoch_data,tidx},'TooltipString','View only good data');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[600 5 80 18],'String','Inspect','Callback',@Inspection,'TooltipString','Inspection mode');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[700 5 80 18],'String','Plot Options','CallBack',{@change_numstd,6,epoch_data,tidx},'TooltipString','Change number of std for threshold');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[800 5 80 18],'String','Update','CallBack',{@VR_update,epoch_data,tidx},'TooltipString','Change number of std for threshold');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[900 5 80 18],'String','View Channel','CallBack',{@view_channel,raw_data},'TooltipString','Change number of std for threshold');
uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[1000 5 80 18],'String','Restore','Callback',{@restore_all},'TooltipString','Restore All');
quit_but=uicontrol(CHANNELDATA.plotfig,'units','pixels','position',[1100 5 80 18],'String','Exit','TooltipString','Quit and Save','Callback',{@quitting});
%Raw Plot Ylim buttons
uicontrol(CHANNELDATA.raw_plot,'Style','text','units','pixels','position',[5 5 40 18],'String','Min');                                                  %min text
CHANNELDATA.raw_min=uicontrol(CHANNELDATA.raw_plot,'Style','edit','units','pixels','position',[55 5 40 18],'Callback',{@VR_zoom_in_out,3},'String',num2str(CHANNELDATA.Ylim_raw(1)));  %min box
uicontrol(CHANNELDATA.raw_plot,'Style','text','units','pixels','position',[115 5 40 18],'String','Max');                                                  %max text
CHANNELDATA.raw_max=uicontrol(CHANNELDATA.raw_plot,'Style','edit','units','pixels','position',[165 5 40 18],'Callback',{@VR_zoom_in_out,3},'String',num2str(CHANNELDATA.Ylim_raw(2)));  %max box
%toggle bad trials
mean_all_chan = mean(CHANNELDATA.max(1).data(:,1));
std_all_chan = std(CHANNELDATA.max(1).data(:,1));
CHANNELDATA.oline = line([mean_all_chan+2.5*std_all_chan mean_all_chan+2.5*std_all_chan],[1 length(CHANNELDATA.chan_num)],'LineStyle','--','Parent',CHANNELDATA.plot_axes,'Color',[0 1 0]);
CHANNELDATA.eoline = line([mean_all_chan+5*std_all_chan mean_all_chan+5*std_all_chan],[1 length(CHANNELDATA.chan_num)],'LineStyle','-','Parent',CHANNELDATA.plot_axes,'Color',[0 1 0]);


%Trial Data
temp = find(CHANNELDATA.toggle);
for cc = 1:length(epoch_data.epochs)
    for i=1:epoch_data.epochs(cc).num_trials
        if cc == 1
            ind=i;
            TRIALDATA.trial_num(i) = i;
            TRIALDATA.trial_cc(i,1) = i;
            TRIALDATA.trial_cc(i,2) = cc;
            [TRIALDATA.max(1).data(i,1) TRIALDATA.max(1).data(i,2)] = max(max(abs(epoch_data.epochs(1).data(temp,tidx,i)),[],2)); %max val, channel
            TRIALDATA.patch(i) = patch(TRIALDATA.trial_num(i),TRIALDATA.max(1).data(i,1),[0 0 0],'Tag' ,'trial','Parent',TRIALDATA.hg,'Marker','.','EdgeColor',[0 0 0],'ButtonDownFcn',{@VR_Sum_points});
        else
            ind=sum([epoch_data.epochs(1:(cc-1)).num_trials])+i;
            TRIALDATA.trial_num(ind) = ind;
            TRIALDATA.trial_cc(ind,1)=i;
            TRIALDATA.trial_cc(ind,2)=cc;
            [TRIALDATA.max(1).data(ind,1) TRIALDATA.max(1).data(ind,2)] = max(max(abs(epoch_data.epochs(cc).data(temp,tidx,i)),[],2)); %max val, channel
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
set(CHANNELDATA.plot_axes, 'Xlim',[0 max(CHANNELDATA.max(1).data(:,1))],'Ylim',[0 epoch_data.num_sensors]);
set(TRIALDATA.plot_axes, 'Ylim',[0 max(TRIALDATA.max(1).data(:,1))],'Xlim',[0 sum([epoch_data.epochs(:).num_trials])]);



%Setup Condition lines
for cc =  1:(length(epoch_data.epochs)-1)
    x=sum([epoch_data.epochs(1:cc).num_trials]);
    y=get(TRIALDATA.plot_axes,'Ylim');
    line([x x],[0 y(2)],'Parent',TRIALDATA.plot_axes,'HitTest','off');
end;
chan_data = find(~CHANNELDATA.toggle);
trial_data = find(~TRIALDATA.toggle);
set(CHANNELDATA.badpanel,'Data',CHANNELDATA.label(chan_data),'ColumnName','Channel','CellSelectionCallback',{@VR_Restore,0});
set(TRIALDATA.badpanel,'Data',TRIALDATA.trial_cc(trial_data,:),'ColumnName',{'Trial','Cond'},'CellSelectionCallback',{@VR_Restore,1},'ColumnWidth',{70 70});

disp('Accounting for all trials and channels');
VR_update([],[],epoch_data,tidx);
if strcmp(rej_input,'Yes')
    for i=1:length(reject_data.badchans)
        x=reject_data.badchans(i);
        CHANNELDATA.toggle(x)=logical(0);
        tmp=get(CHANNELDATA.plotdata(1).axes(x),'Title');
        VR_Channel_disable(tmp,[]);
    end;
    
    for cc=1:length(reject_data.badtrials)
        for i=1:length(reject_data.badtrials{cc})
            x = reject_data.badtrials{cc}(i); %x is relative trial number
            index=find((x==TRIALDATA.trial_cc(:,1)) & (cc==TRIALDATA.trial_cc(:,2))); %index is absolute trial number
            set(TRIALDATA.patch(index),'Visible','off','HitTest','off','Selected','off');
            TRIALDATA.toggle(index)=logical(0);
            set(CHANNELDATA.plotdata(cc).plot_lines(:,x),'Visible','off','Selected','off','HitTest','off')
        end;
    end;
end;
set(sumplot,'Visible','on');
set(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan(1),'Visible','on');
disp('Accounting for only the good data');
VR_update([],[],epoch_data,tidx);
repeatrej = 1;
while repeatrej && ishandle(CHANNELDATA.plotfig)
    refresh(CHANNELDATA.plotfig)
    if CHANNELDATA.quit == 0,
        uiwait;
    else
        %     waitfor(quit_but,'UserData','done');
        % saving reject data structure
        if ~exist('reject_data','var') || ~isfield(reject_data,'event_code')
            reject_data.event_code = [epoch_data.epochs.event_code];
        end;
        if ~isfield(reject_data,'badtrials'), reject_data.badtrials = []; end;
        if eeg == 0
            if ~isfield(reject_data,'chanindex'), reject_data.chanindex = 'grads'; end;
        else
            % eeg analysis
            if ~isfield(reject_data,'chanindex'), reject_data.chanindex = 'eeg'; end;
        end;
        if ~isfield(reject_data,'badchans'), reject_data.badchans=[]; end;
        badtrials={};
        for ind=1:CHANNELDATA.num_conds
            badtrials{ind}=[];
        end;
        reject_data.badchans = find(~CHANNELDATA.toggle);
        reject_data.badchanlabels = {epoch_data.sensor_info(reject_data.badchans).label};
        disp('Current reject data structure:');
        reject_data
        savereject = 'n';
        exit_input = questdlg('Are you sure you want to exit','Exit Menu');
        if strcmpi(exit_input,'Yes')
            % bad trials
            abs_bad_trials = find(~TRIALDATA.toggle); %Listing of absolute bad trials by absolute number
            for ti=1:length(abs_bad_trials)
                tmp=abs_bad_trials(ti);
                trial_cc=TRIALDATA.trial_cc(tmp,1);
                cc=TRIALDATA.trial_cc(tmp,2);
                badtrials{cc}=[badtrials{cc} trial_cc];
            end;
            reject_data.badtrials=badtrials;
            reject_data.badchans = find(~CHANNELDATA.toggle);
            reject_data.badchanlabels = {epoch_data.sensor_info(reject_data.badchans).label};
            
            disp('Corrected reject data structure: ')
            for i = 1: length(reject_data.event_code)
                reject_data.badtrials{i}=sort(reject_data.badtrials{i});
                Num_badtrials(i) = length(reject_data.badtrials{i});
            end;
            reject_data
            clear summary_trials
            Num_trials_orig = [epoch_data.epochs(:).num_trials];
            Num_trials_new = Num_trials_orig - Num_badtrials;
            print_out={'Please review the rejection below'; 'click yes to save'};
            print_out{end+1}=['Rejected Channels:'];
            if ~isempty(reject_data.badchanlabels)
                for indo=1:length(reject_data.badchanlabels)
                    print_out{end+1}=reject_data.badchanlabels{indo};
                end;
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
            happy = questdlg(print_out,'Save Dialog','Yes','Yes with log','No','Yes');
            for i=3:length(print_out);
                disp(print_out{i});
            end;
            if strcmpi(happy,'Yes')
                disp(['Saving reject data to ',rejectfilenew])
                save(rejectfilenew,'reject_data')
                repeatrej = 0;
            elseif strcmpi(happy,'Yes with log')
                disp(['Saving reject data to ',rejectfilenew])
                save(rejectfilenew,'reject_data')
                repeatrej = 0;
                log_out(raw_data,epoch_data,root,subj,task);
            elseif strcmpi(happy,'No')
                disp(['Not saving Data']);
                repeatrej = 0;
            end;
        else
            set(quit_but,'UserData','idle');
        end;
    end;
end;
    userinfo = {};
    userinfo{1} = get(sumplot,'Position');
    userinfo{2} = get(CHANNELDATA.plotfig,'Position');
    userinfo{3} = get(CHANNELDATA.raw_plot,'Position');
    userinfo{4} = {subj num2str(ICA) band};
    save('~/matlab/VR_gui_userinfo','userinfo');


close all
clear all


function next_page(varargin)
global CHANNELDATA
global TRIALDATA
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
    Ylim_raw=get(axes_raw(1),'Ylim');
    set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
    set(h_raw,'Visible','off');
    new_raw_axes = get(new_h_raw,'Children');
    set(new_raw_axes,'Ylim',Ylim_raw);
    set(new_h_raw,'Visible','on');
end;
uiresume;


function prev_page(varargin)
global CHANNELDATA
global TRIALDATA
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
uiresume;


function next_cc(varargin)
global CHANNELDATA
global TRIALDATA

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
    set(axes,'Ylim',CHANNELDATA.Ylim);
    set(new_h,'Visible','on');
    set(h_raw,'Visible','off');
    axes_raw = get(h_raw,'Children');
    set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
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
uiresume;

function prev_cc(varargin)
global CHANNELDATA
global TRIALDATA
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
    set(axes,'Ylim',CHANNELDATA.Ylim);
    axes_raw = get(h_raw,'Children');
    set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
    
    
    
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
uiresume;

function view_all(varargin)
global CHANNELDATA
global TRIALDATA
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
uiresume;

function selection(varargin)
global CHANNELDATA
global TRIALDATA
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
uiresume;

function good_only(varargin)
global CHANNELDATA
global TRIALDATA
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
uiresume;

function Inspection(varargin)
global CHANNELDATA
global TRIALDATA
selection(1);
CHANNELDATA.view_mode='Inspection';
set(CHANNELDATA.plotfig,'Name',['Condition: ',num2str(CHANNELDATA.c_cond), ' View Mode: ',CHANNELDATA.view_mode]);
set(CHANNELDATA.patch,'ButtonDownFcn',{@VR_Sum_points});
set(TRIALDATA.patch,'ButtonDownFcn',{@VR_Sum_points});
uiresume;

function view_channel(h,x,raw_data)
global CHANNELDATA
global TRIALDATA
    
prompt = {'Which Channel would you like to view','Which Condition'};
def={'MEG 0000','1'};
dlg_title = 'Plot channel';
num_lines = 1;
chan = inputdlg(prompt,dlg_title,num_lines,def);
if (chan{2}>length(raw_data.epochs))
    cc=length(raw_data.epochs);
else
    cc = chan{2};
end;
chan = chan(1);
chan_num = find(strcmpi(CHANNELDATA.label,chan));
trials = find(TRIALDATA.toggle & TRIALDATA.trial_cc(:,2)==3);
trials = TRIALDATA.trial_cc(trials);
data=raw_data.epochs(cc).data(chan_num,:,trials);
data=mean(data,3);
fig = figure; %('NumberTitle','off','Name',['Showing ' chan, 'on condition : ',num2str(cc)]);
ax = axes('Parent',fig);
line(raw_data.epochs(cc).time,data,'Parent',ax);


    




function change_numstd(h,x,numstd,epoch_data,tidx)
global CHANNELDATA
global TRIALDATA
prompt = {'Enter Number of Std for Threshold:','Enter New Y-limit'};
dlg_title = 'Change Number Plot Properties';
num_lines = 1;
def = {'6',num2str(CHANNELDATA.Ylim(2))};
values = cell2num(inputdlg(prompt,dlg_title,num_lines,def));
userstd=values(1);
Ylim = values(2);
for cc=1:length(CHANNELDATA.plotdata)
    for ci=1:length(CHANNELDATA.plotdata(cc).threshold_lines(:,1))
        set(CHANNELDATA.plotdata(cc).threshold_lines(ci,1),'Ydata',abs(mean(epoch_data.epochs(cc).data(ci,tidx,:),3))+userstd*std(epoch_data.epochs(cc).data(ci,tidx,:),0,3));
        set(CHANNELDATA.plotdata(cc).threshold_lines(ci,2),'Ydata',-(abs(mean(epoch_data.epochs(cc).data(ci,tidx,:),3))+userstd*std(epoch_data.epochs(cc).data(ci,tidx,:),0,3)));
    end;
end;
    CHANNELDATA.Ylim = [0 Ylim];
    for cc=1:CHANNELDATA.num_conds
    set(CHANNELDATA.plotdata(cc).axes,'Ylim',CHANNELDATA.Ylim);
    end;
   
    
set(h,'CallBack',{@change_numstd,userstd,epoch_data,tidx});

function quitting(varargin)
global CHANNELDATA
h=varargin{1};
CHANNELDATA.quit=1;
set(h,'UserData','done');
uiresume;

function restore_all(varargin)
global CHANNELDATA
global TRIALDATA
c_cond=CHANNELDATA.c_cond;
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);

CHANNELDATA.toggle(:)=logical(1);
TRIALDATA.toggle(:)=logical(1);
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
uiresume;

function view_raw(h,x)
global CHANNELDATA
set(CHANNELDATA.raw_plot,'Visible','on');
set(h,'Callback',@hide_raw,'String','Hide Raw');
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);


function hide_raw(h,x)
global CHANNELDATA
%Restore zoom on panels
h=findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).panels,'flat','Visible','on');
h_raw = findobj(CHANNELDATA.plotdata(CHANNELDATA.c_cond).raw_pan,'flat','Visible','on');
axes = get(h,'Children');
set(axes,'Ylim',CHANNELDATA.Ylim);
axes_raw = get(h_raw,'Children');
set(axes_raw,'Ylim',CHANNELDATA.Ylim_raw);
set(CHANNELDATA.raw_plot,'Visible','off');
set(h,'Callback',@view_raw,'String','Show Raw');

function log_out(raw_data,epoch_data,root,subj,task)
global CHANNELDATA
global TRIALDATA
disp(['Saving log file to ' root,'/',subj,'/',subj,'_',task,'_edited_events_rejected.txt']);
epoch_info = [];
reject_info = [];
for cc=1:CHANNELDATA.num_conds
    latency_info=[];
    for ti = 1:length(raw_data.epochs(cc).trial_info)
        if raw_data.epochs(cc).trial_info(ti).badtrial == 0
            latency_info(end+1) = raw_data.epochs(cc).trial_info(ti).latency;
        else
            reject_info(end+1) = raw_data.epochs(cc).trial_info(ti).latency;
        end;
    end;
    clear ti
    
    epoch_info(cc).latency_info = latency_info;
end;
bad_trial=find(~TRIALDATA.toggle);
for  ti = 1:length(bad_trial)
    absindex = bad_trial(ti);
    ccind = TRIALDATA.trial_cc(absindex,1);
    cond = TRIALDATA.trial_cc(absindex,2);
    reject_info(end+1) = epoch_info(cond).latency_info(ccind);
end;

events=ts_import_events([root,'/',subj,'/',subj,'_',task,'_edited_events_rejected.txt']);
%Edit event to have rejected trials
for ev=1:length(events)
    if(any(events(ev).latency==reject_info))
        events(ev).type='reject';
        events(ev).condition=0;
    end;
end;
ts_export_events(events,[root,'/',subj,'/',subj,'_',task,'_edited_events_rejected.txt']);

