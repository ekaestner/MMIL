function [] = ts_iEEG_movie_maker(datain,start_time,stop_time,time_window,labels,filename,varargin)
%
% [] = ts_iEEG_movie_maker(datain,start_time,stop_time,time_window,labels,filename,event_ind,cfg)
%
% This program makes avi files from epoch data. 
%
% ARGUMENTS
% datain -   structure in the form of epoch_data or avg_data
% start_time -  start time in seconds to begin movie
% stop_time -   stop time in seconds for end of movie
% time_window - size of the time window to take mean from for each frame
% labels -      a list of the grid names, this might be helpful:
%                labels=strcat(repmat('G',64,1),num2str([1:64]','%d'));   
% filename -    filename to save avi file to.
% event_ind -    the event number (NOT REQUIRED FOR epoch_data, IS REQUIRED FOR avg_data)
%
% OPTIONAL
% cfg - optional structure for customizing the look of the movies.
%
%   PARAMETER       DEFAULT VALUE                   DESCRIPTION
%   colorbar        cfg.colorbar = true;            show colorbar
%   showcontour     cfg.showcontour = true;         show contourlines
%   numcontourlines cfg.numcontourlines = 6;        number of contourlines
%   gridnames       cfg.gridnames = true;           show grid names
%   showgrid        cfg.showgrid = true;            show grid node locations
%   gridorient      cfg.gridorient = 'SWNW';        what corners G1 and G8 are, respectively.
%                                                   other values:'NESE','NWSW','SENE','SESW' or 'NWNE'
%   encode_offline  cfg.encode_offline = false;     writes out .png files instead of an .avi
%   clims           cfg.clims = 'minmax';           color and plot thresholds,
%                                                   values:'absmax',and
%                                                   user specified [min max]
%
% For example, if we wanted to make a movie of a sliding window of 10ms
% from -500ms to 1000ms of the first condition:
%
%   Lbl=strcat(repmat('G',64,1),num2str([1:64]','%d'));
%   load('matfiles\TT2_S100.rej.epoch.event1.mat')
%   ts_iEEG_movie_maker(epoch_data,-0.5,1.0,.010,Lbl,'TT2_S100_C1_from_-500_to_1000ms_wind_10ms.avi')
%
% OR, using the avg_data:
%
%   load('matfiles\TT2_S100.reg.avg.mat')
%   ts_iEEG_movie_maker(avg_data,-0.5,1.0,.010,Lbl,'TT2_S100_C1_from_-500_to_1000ms_wind_10ms.avi',1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: Brian T. Quinn 06/23/2008
%
%TODO:
%1[X]get cfg args working - DONE
%2[X]arg for which channels one wants - DONE labels
%3[\ ]get offline rendering commands
%4[\ ]orient grid when given specification
%5[X]grid dimensions! ^SEE #2
%6[ ]ah1 ah2 stuff
%7[ ]highlight certain channels


%PRINT HELP IF NO ARGS PASSED
if nargin==0,
    help ts_iEEG_movie_maker
    return
end

%PARSE_ARGS
cfg=[];
event_ind={};
if ~isempty(varargin),
  for i=1:length(varargin),
      if isstruct(varargin{i}),
          cfg=varargin{i};                            
      else event_ind=varargin{i};
      end
  end
end

% cfg.plotTCxlim = 'roi'
if ~isfield(cfg, 'colorbar'),       cfg.colorbar = true;                    end
if ~isfield(cfg, 'showcontour'),    cfg.showcontour = true;                 end
if ~isfield(cfg, 'numcontourlines'),cfg.numcontourlines = 6;                end
if ~isfield(cfg, 'gridnames'),      cfg.gridnames = true;                   end
if ~isfield(cfg, 'showgrid'),       cfg.showgrid = true;                    end
if ~isfield(cfg, 'gridorient'),     cfg.gridorient = 'SWNW';                end
if ~isfield(cfg, 'clims'),          cfg.clims = 'minmax';                   end

%below are redundant variables, commented out for now....
%if ~isfield(cfg, 'event_code'),     cfg.event_code = 1;                     end
%if ~isfield(cfg, 'filename'),       cfg.filename = 'temp.avi';              end 
%if ~isfield(cfg, 'start_time'),     cfg.start_time = 0;                     end
%if ~isfield(cfg, 'stop_time'),      cfg.stop_time = 1;                      end
%if ~isfield(cfg, 'time_window'),    cfg.time_window = .01;                  end

%for future development...
if ~isfield(cfg, 'encode_offline'), cfg.encode_offline = false;             end
% cfg.plottitle = ''

%datain
Lbl = {datain.sensor_info.label};
[COI,ia]=intersect(Lbl,labels);
[grid_inds,grid_sort_ind]=sort(ia);
COI=COI(grid_sort_ind);

%load data into array we will manipulate later
if isfield(datain,'averages'), %if avg_data structure passed
    if isempty(event_ind), error('Must specify an event code when passing avg_data'); end
    event_ind=find([datain.averages(:).event_code]==event_ind); %determine which event cond is desired
    mCond=datain.averages(event_ind).data(grid_inds,:,:);                   %grad necessary data
    timevect=datain.averages(event_ind).time;
    event_code=datain.averages(event_ind).event_code;
else                                                       %otherwise epoch_data was passed
    mCond=squeeze(mean(datain.epochs.data(grid_inds,:,:),3));              %take avg
    timevect=datain.epochs.time;
    event_code=datain.epochs.event_code;
end

%get indecies for time
[dummy_val,start_ind]=min(abs(timevect - start_time));
[dummy_val,temp_wind]=min(abs(timevect - (start_time+time_window)));
ind_window=temp_wind-start_ind;ind_window=ind_window+mod(ind_window,2); %make it even so ind_window/2 is an integer
[dummy_val,stop_ind]=min(abs(timevect - stop_time));

%time_size=size(mCond,2);

% Color Limits
if isvector(cfg.clims), amin=cfg.clims(1);amax=cfg.clims(2);    end  %user defined
if strcmp(cfg.clims,'minmax'),  %find min and max of time segment and add 10%
    amin=min(min(mCond(:,start_ind:stop_ind)))*1.1;
    amax=max(max(mCond(:,start_ind:stop_ind)))*1.1;
end
if strcmp(cfg.clims,'absmax'),  %use absolute max +/-
    amax=max(max(abs(mCond(:,start_ind:stop_ind))));
    amin=-amax;
end


%Grid Dimensions
%find the smallest square grid to use for length of data
%this should be somewhere in data structure in future 
if ~isfield(cfg,'griddims'),
    if strncmp(COI(1),'G',1), cfg.griddims=[ceil(sqrt(length(COI))) ceil(sqrt(length(COI)))];
    else cfg.griddims=[length(COI) 1];
    end
end
griddims=cfg.griddims;

%create grids for interpolation later - ASSUMES SQUARE GRID CURRENTLY
[x,y]=meshgrid(1:griddims(1));
[xi,yi]=meshgrid(linspace(0,griddims(1)+1,160));

%Open a new figure and start frame/image time count
figure
k=0;

%Proceed from start to stop, with 50% overlapping windows
for t=start_ind:ind_window/2:stop_ind-ind_window,
    clf
    k=k+1;
    smallframe=nan(griddims); %blank the smallframe
    smallframe(1:length(COI))=mean(mCond(:,t:t+ind_window),2); %default gridorient is SWNW
    gridx=x;gridy=y;
    if strcmp(cfg.gridorient,'NESE'),
        smallframe=fliplr(flipud(smallframe));
        gridx=fliplr(x);gridy=flipud(y);
    end
    if strcmp(cfg.gridorient,'NWSW'),
        smallframe=flipud(smallframe);
        gridx=x;gridy=flipud(y);
    end
    if strcmp(cfg.gridorient,'SENE'),
        smallframe=fliplr(smallframe);
        gridx=fliplr(x);gridy=y;
    end
    if strcmp(cfg.gridorient,'SESW'),
        smallframe=rot90(smallframe,-1);
        gridx=rot90(gridx);gridy=rot90(gridy);
    end
    if strcmp(cfg.gridorient,'NWNE'),
        smallframe=rot90(smallframe);
        gridx=rot90(gridx,-1);gridy=rot90(gridy,-1);
    end
    if (sum(sum(isnan(smallframe)))>0),                 %if badchans or nonexisting chans
        smallframe=interpolate_missing(smallframe);        
    end
    %TODO: move this out of loop and set position as well.do same for ah2. use
    %axes(ah1) to set active. will allow ah1 to fill figure window when the
    %timeseriesplots are not wanted.
    ah1=axes('OuterPosition',[0 .1 1.0 .85]);
    zi=griddata(x,y,double(smallframe),xi,yi,'v4');%consider other interp methods
    if (cfg.showcontour),   contour(xi,yi,zi,cfg.numcontourlines,'k');      end
    surface(xi,yi,zeros(size(zi)),zi,'EdgeColor','none','FaceColor','flat');
    alim([0 1]);
    caxis([amin amax]);
    hold on
    %plot gridpoints and grid labels
    if (cfg.showgrid), plot(gridx,gridy,'ko','MarkerFaceColor','g');        end
    if (cfg.gridnames),
        for i=1:length(COI),
            text(gridx(i),gridy(i),COI{i},'HorizontalAlignment','center','VerticalAlignment','middle','Color','black','FontSize',10,'FontWeight','bold')
            if (datain.sensor_info(i).badchan==1), plot(gridx(i),gridy(i),'ko','MarkerFaceColor','r'); end
        end
    end
    xlim([.5 griddims(1)+.5]);ylim([.5 griddims(2)+.5]); %set axis limits to .5 border around data
    axis equal;axis off
    titletext=sprintf('Event Code %i, time=%-4.0f to %-4.0fms', event_code,timevect(floor(t))*1000,timevect(floor(t+ind_window))*1000);
    title(titletext,'FontSize',16);%,'HorizontalAlignment','left');
    if (cfg.colorbar),            colorbar;                                 end
    
    %begin timeseries plot
    ah2=axes('OuterPosition',[0 0 1.0 .25]);
    plot(timevect(start_ind:stop_ind).*1000,mCond(:,start_ind:stop_ind),'k')
    vline(timevect(t:t+ind_window).*1000,'y')
    hold on
    plot(timevect(start_ind:stop_ind).*1000,mCond(:,start_ind:stop_ind),'k')
    vline(0,'g')
    xlabel('Time (ms)')
    ylabel('mV')
    axis tight;axis fill
    ylim([amin amax]);
    if (cfg.encode_offline),                %if encode offline, create image filenames
        %fileparts=strsplit(filename,'.');
        %imagename=sprintf('im%04d_%s.png',k,fileparts{1});
        imagename=sprintf('im%04d_event%i_%-0.0fto%-0.0f.png',k,event_code,timevect(floor(t))*1000,timevect(floor(t+ind_window))*1000);
        print('-dpng', imagename);          %save as .png file
    else
        F(k) = getframe(gcf);
    end
end
%write out avi file
if ~(cfg.encode_offline),   movie2avi(F,filename,'FPS',5);      end
clear F
%close all %close figure windows


function [outdat]=interpolate_missing(indat)
% [outdat]=interpolate_missing(indat)
% date:07/03/2008
% by:BTQ
% takes in a matrix and uses linear interpolation of the adjacent points to
% reconstruct the missing point
%
dims=size(indat);
[xnan,ynan]=find(isnan(indat)==1);
indat(isnan(indat))=0; %turn nans to zeros, we'll ignore later
wn=1;
wd=1/sqrt(2);
kern=[wd wn wd;wn 0 wn;wd wn wd]; %kernel of weights we will use for convolution
if ~isempty(xnan),
    for i=1:length(xnan),
        xn=xnan(i);yn=ynan(i);
        blank=zeros(size(indat));blank(xn,yn)=1; %put a 1 at location to interpolate
        cblank=conv2(blank,kern);cblank=cblank(2:dims(1)+1,2:dims(2)+1); %cblank is blank conv'd w/ kernel
        [nbrx,nbry]=find(cblank~=0);
        wtot=0;val=0;
        for j=1:length(nbrx),
            val=val+cblank(nbrx(j),nbry(j))*indat(nbrx(j),nbry(j));
            wtot=wtot+cblank(nbrx(j),nbry(j))*logical(indat(nbrx(j),nbry(j))); %adds nothing if indat==0;
        end
        indat(xn,yn)=val/wtot;
    end
end
outdat=indat;

function [] = vline(x_pos,varargin)
% usage: [] = vline(x_pos)
% date: 6/11/2007
% by: BTQ
% purpose: put a vertical dotted on an axis
%         without changing it
if nargin==2,
    color=cell2mat(varargin);
else
    color=[1 0 0];
end
%obtain the y range from the current axis
yrange = get(gca,'YLim');
%plot the line
hold on
if length(x_pos)>1
    line([x_pos;x_pos],repmat(yrange,length(x_pos),1)','LineStyle','-','Color',color)
else
    plot([x_pos x_pos],yrange,'LineStyle','-','Color',color)
end
hold off
