function ts_iEEG_Plot_TFR (data,varargin)

% Use:  ts_iEEG_Plot (data,'conditions',[1 4 5],...)
%
% Note: If you plot lots of conditions together, MATLAB can get bogged down 
%       considerably.  In addition more than 16 may cause anomalies with
%       the legends or cause the program to crash.
%
% Required Input:
%    
%   data - a TimeSurfer data structure that can be plotted: either should
%          be avg_data or timefreq_data
%
% Optional Parameters: 
%
%    conditions     - the conditions (not event numbers) to plot  
%                     {default is the smaller of: all conditions or first 4 conditions}
%    channels       - cell array of channel names to plot.  {default = all}
%    badchans       - list of channels that should not be plotted {default = []} 
%    baseline       - 'yes','relative','relchange','absolute', 'no' 
%                      {default = 'no'} {'yes' = 'absolute'}
%    baseline_start - start time for calculation of baseline {default = first time point}
%    baseline_end   - end time for calculation of baseline   {default = 0}
%    plotzscores    - 'yes' or 'no' {default = 'yes'} - note that the
%                      baseline_start and _end options may be used to
%                      specifed the time period over which to calculate the
%                      mean and standard deviation
%    freqcorr       - when plotting power, perform a 1/f^3 correction {default = 'no'}
%    xlim           - [xmin xmax] - x-axis limits {default based on max and min of time}
%    ylim           - [ymin ymax] - y-axis limits {default based on max and min of data}
%    zlim           - [zmin zmax] - z-axis limits {default based on max and min of data}
%                     or 'absmax'
%    showaxis       - 'yes' or 'no' - show MATLAB axis {default = 'yes'}
%    showzeroaxis   - 'yes' or 'no' - show axis where x and y = zero {default = 'yes'} 
%    conditionnames - cell array of strings describing each condition that
%                     is plotted {default = condition numbers}
%    footnotes      - string describing data {default = ''}
%    figname        - template file name for the figure
%    title          - title placed in figure
%    layoutfile     - a valid neuroscan montage or fieldtrip layout file
%    savefigs       - 'yes' or 'no' - save figure to file {default = 'no'}
%    clearfigs      - 'yes' or 'no' - clear figures off screen {default ='no'}
%    stat_data       - a valid stat_data data structure
%    laminar        - 'yes' or 'no'
%
% Advanced options:
%
%    max_rows, max_columns, max_group_size, plot_colors, figbackground,
%    zeroaxiswidth, zeroaxisstyle, axisfontsize, papertype, leftmargin,
%    rightmargin, topmargin, bottommargin, paperorient, shading, whichfigs,
%    rowpad, columnpad.
%
% Created:  09/25/07 by Rajan Patel
% Rcnt Mod: 04/09/08 by Rajan Patel
% Last Mod: 09/15/12 by Don Hagler
%
% See also: ts_iEEG_PlotAvg, ts_iEEG_PlotWavelets

%% Check for errors in calling function

if nargin < 1
   help(mfilename);
   return;
end

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), opt=struct(options{:}); 
  else opt= []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

%% Data Options
                                                   
if isfield (data,'timefreq')
   datafield = 'timefreq';
else
  fprintf('ERROR: Supplied data is not supported by this function.\n');
  return;
end

if isfield(opt,'conditions')
  conds = opt.conditions;
  for i = 1:length(conds)
   if conds(i) > length(data.(datafield))
    fprintf('ERROR: Condition %s is outside of the range of data supplied.\n',num2str(conds(i)));
    return
   end
  end
   if length(conds)>16
    fprintf('WARNING: More than 16 conditions may cause problems.\n');
   end
else
  num_conditions  = length(data.(datafield));
  num_condstoplot = min(4, num_conditions); 
  conds           = 1:num_condstoplot;                                     % default plot all or first 16 conditions
end

if isfield(opt, 'channels') && ~isempty(opt.channels)
    if ~iscell(opt.channels), opt.channels = {opt.channels}; end
    channels = opt.channels;  
    channels = intersect({data.sensor_info.label},channels);               % find the channels that match those in the data
    if isempty(channels)
        fprintf('ERROR: %s not found as channels in given data.\n',channels{:});
        return;
    elseif length(channels) ~= length(opt.channels)
        fprintf('WARNING: One or more channels not found in supplied data:\n%s\n',cell2mat(setdiff(opt.channels,channels)));
    end
else
    channels = {};
end

if isfield(opt, 'badchans')                                                % channels not to plot
  ignore_channels = opt.badchans;
else
  ignore_channels = [];
end

if isfield(opt,'baseline_data') && ~isempty(opt.baseline_data)
    baseline_data = opt.baseline_data;
    [data.sensor_info(find([baseline_data.sensor_info.badchan])).badchan] = deal(1); % add bad channels of baseline
else
    baseline_data = [];
end

if isfield(opt,'showbadchans') && strcmpi(opt.showbadchans,'no')
 badchan_idx     = find([data.sensor_info.badchan]);                        % don't plot channels marked bad
 ignore_channels = [{data.sensor_info(badchan_idx).label} ignore_channels];
end

baseline = 4;
if isfield(opt,'baseline')
  if strcmpi(opt.baseline,'yes') || ...
     strcmpi(opt.baseline,'absolute'), baseline = 1; end
  if strcmpi(opt.baseline,'relative'), baseline = 2; end
  if strcmpi(opt.baseline,'relchange'),baseline = 3; end
  if strcmpi(opt.baseline,'no'),       baseline = 0; end
end

if isfield(opt,'plotzscores') && strcmpi(opt.plotzscores,'yes')
    baseline = 4; 
end

if isfield(opt,'baseline_start')
  baseline_start = opt.baseline_start;
else
  baseline_start = [];
end

if isfield(opt,'baseline_end')
 baseline_end = opt.baseline_end;
else
 baseline_end = [];
end

plot_stats = 0;
if isfield(opt,'stat_data') && ~isempty(opt.stat_data)
  plot_stats = 1;
end

freqcorr = 0;
if isfield(opt,'freqcorr') && strcmpi(opt.freqcorr,'yes')
  freqcorr = 1;
end

timedsfact = 0;
if isfield(opt,'timedsfact') && opt.timedsfact > 0
    timedsfact = opt.timedsfact;
end

freqdsfact = 0;
if isfield(opt,'freqdsfact') && opt.freqdsfact > 0
    freqdsfact = opt.freqdsfact;
end

%% Figure Options 

if isfield(opt,'xlim') && ~isempty(opt.xlim)                                                    
  xmin = opt.xlim(1);
  xmax = opt.xlim(2);
else
  xmin = [];
  xmax = [];
  for i=1:length(conds)                                                     
    xmin = min([xmin data.(datafield)(conds(i)).time]);                     % look for x min across all selected conditions
    xmax = max([xmax data.(datafield)(conds(i)).time]);                     % look for x max across all selected conditions
  end
end

if isfield(opt,'ylim') && ~isempty(opt.ylim)
  ymin = opt.ylim(1);
  ymax = opt.ylim(2);
  calc_ylim = 0;
else
  calc_ylim = 1;
  ymin = [];
  ymax = [];
end

abs_zlim=0;
calc_zlim = 1;
zmin = [];
zmax = [];
if isfield(opt,'zlim')
 if strcmpi(opt.zlim,'absmax')
  abs_zlim=1; 
 elseif ~isempty(opt.zlim)
  zmin = opt.zlim(1);
  zmax = opt.zlim(2);
  calc_zlim = 0;
 end
end

show_axis = 1 ;
if isfield(opt,'showaxis')
    if strcmpi(opt.showaxis,'no'), show_axis=0; end
end

center_axis=1;
if isfield(opt,'showzeroaxis')
    if strcmpi(opt.showzeroaxis,'no'), center_axis=0; end
end

center_axis_x_ticks = 1;
center_axis_y_ticks = 0;
if center_axis
  if isfield(opt,'showzeroaxisxticks')
    if strcmpi(opt.showzeroaxisxticks,'no'), center_axis_x_ticks = 0; end
  end
  if isfield(opt,'showzeroaxisyticks')
    if strcmpi(opt.showzeroaxisyticks,'yes'), center_axis_y_ticks = 1; end
  end
end

show_colorbar = 1;
if isfield(opt,'showcolorbar')
    if strcmpi(opt.showcolorbar,'no'), show_colorbar = 0; end
end

if isfield(opt,'conditionnames')
 condition_names = opt.conditionnames;
else
 condition_names = {};
end

if isfield(opt,'footnotes')
  foot_notes = opt.footnotes;
else
  foot_notes = '';
end

if isfield(opt,'figname')
 fig_name = opt.figname;
 [fig_path, fig_name, fig_ext] = fileparts(fig_name);
 fig_name = [fig_path fig_name];
else
 fig_name = {};
end

if isfield(opt,'title')
 title_string = opt.title;
else
 title_string = '';
end

if isfield(opt,'resolution') && ~isempty(opt.resolution)
  resolution = num2str(opt.resolution);
else
  resolution = '300';
end

save_figs = 0;
if isfield(opt,'savefigs');
  opt.savefigs = lower(opt.savefigs);
  if ~strcmp(opt.savefigs,'no')
    save_figs = 1;
    switch opt.savefigs
      case 'bmp'
        print_command{1} = '-dbmp';
      case 'emf'
        print_command{1} = '-dmeta';
        print_command{2} = sprintf('-r%s',resolution);
      case 'eps'
        print_command{1} = '-depsc';
        print_command{2} = '-tiff';
        print_command{3} = sprintf('-r%s',resolution);
      case 'hdf'
        print_command{1} = '-dhdf';
        print_command{2} = sprintf('-r%s',resolution);
      case 'ill'
        print_command{1} = '-dill';
        print_command{2} = sprintf('-r%s',resolution);
        opt.paperorient  = 'portrait';
      case 'jpg'
        print_command{1} = '-djpeg';
        print_command{2} = sprintf('-r%s',resolution);
      case 'pbm'
        print_command{1} = '-dpbm';
      case 'pcx'
        print_command{1} = '-dpcx24b';
      case 'pdf'
        print_command{1} = '-dpdf';
      case 'png'
        print_command{1} = '-dpng';
        print_command{2} = sprintf('-r%s',resolution);
      case 'ppm'
        print_command{1} = '-dppm';
      case 'tif'
        print_command{1} = '-dtiff';
        print_command{2} = sprintf('-r%s',resolution);
      otherwise
        fprintf('Using default EPS format to save figures.\n');
        opt.savefigs = 'eps';
        print_command{1} = '-depsc';
        print_command{2} = sprintf('-tiff');
        print_command{3} = sprintf('-r%s',resolution);
    end
  end
end

clear_figs = 0;
if isfield(opt,'clearfigs')
 if strcmpi(opt.clearfigs,'yes'), clear_figs = 1; end
end

autolayout = 1;
if isfield(opt,'layoutfile') && ~isempty(opt.layoutfile)
  if exist(opt.layoutfile,'file')
    autolayout = 0;
  else
    fprintf('Could not find specified layout file: %s.\n',opt.layoutfile);
  end
end

laminar = 0;
channamepos   = [.5 1];
channamealign = 'center';
% if isfield(opt,'laminar') && strcmpi(opt.laminar,'yes')
%     laminar = 1;
%     channamepos = [0 .5];
%     channamealign = 'right';
% end

%% Hidden Options

if isfield(opt,'max_rows')                                                 % maximum rows of subplots
 max_rows = opt.max_rows;
else
 max_rows = 10;
end
 
if isfield(opt,'max_columns')                                              % maximum columns of subplots
 max_columns = opt.max_columns;
else
 max_columns = 10;
end

if isfield(opt, 'max_group_size')                                          % if a chan group is larger than this it gets its own page
 max_group_size = opt.max_group_size;
else
 max_group_size = 30;
end

if isfield(opt,'plot_colors')                                             % list of colors to toggle through for condition averages
 plot_colors = opt.plot_colors;
else
 plot_colors = 'brcmykwg';
end

if isfield(opt,'plotlinewidth')
  plot_line_width = opt.plotlinewidth;
else
  plot_line_width = 0.5;
end

if isfield(opt,'statslinewidth')
  stats_line_width = opt.statslinewidth;
else
 if plot_line_width < 1
  stats_line_width = 0.1;
 else
  stats_line_width = 0.1;
 end
end

if isfield(opt,'statsmarker')
  stats_marker = opt.statsmarker;
else
  stats_marker = '*';
end

if isfield(opt,'markersize');
  marker_size = opt.markersize;
else
  marker_size = 3;
end

if isfield(opt,'figbackground')                                           % background color of figure
 fig_background = opt.figbackground;
else
 fig_background = 'w';
end

if isfield(opt,'zeroaxiswidth')                                           % width of zero axis line if used
 zero_line_width = opt.zeroaxiswidth;
else
 zero_line_width = 0.25;                                                          
end

if isfield(opt,'zeroaxisstyle')
 line_style = opt.zeroaxisstyle;                                          % style of zero axis line
else
 line_style = '-';                    
end

if isfield(opt,'axisfontsize')                                            % size of axis font in points
 axis_font = opt.axisfontsize;
else
 axis_font  = 4;      
end

if isfield(opt,'papertype')                                               % paper size when saving figure
 papertype = opt.papertype;
else
 papertype = 'usletter';
end

if isfield(opt,'leftmargin')                                              % margins for paper
 leftmargin = opt.leftmargin;
else
 leftmargin = 0.5;
end

if isfield(opt,'rightmargin')
 rightmargin = opt.rightmargin;
else
 rightmargin = 0.5;
end

if isfield(opt,'topmargin')
 topmargin = opt.topmargin;
else
 topmargin   = 0.5;
end

if isfield(opt,'bottommargin')
 bottommargin = opt.bottommargin;
else
 bottommargin = 0.5;
end

if isfield(opt,'paperorient')                                              % printing landscape
 paperorient = opt.paperorient;
else
 paperorient = 'landscape';
end

if isfield(opt,'shading')                                                  % type of shading for time frequency plots
 shading_style = opt.shading;                                              % is a pcolor option
else
 shading_style = 'interp';
end

if isfield(opt,'whichfigs')                                                % if there is more than one page of figures which one to plot
    whichfigs = opt.whichfigs;
else
    whichfigs = [];
end

if isfield(opt,'rowpad')                                                   % these control how close plots are to each other           
  rowpad = opt.rowpad;
else
  rowpad = 0.02;
end

if isfield(opt,'columnpad')
  columnpad = opt.columnpad;
else
  columnpad = 0.01;
end

if isfield(opt,'showonlyfirstaxis');
  showonlyfirstaxis = opt.showonlyfirstaxis;
  if strcmpi(showonlyfirstaxis,'no')
    rowpad = rowpad + 0.02;
    columnpad = columnpad  + 0.02;
  end
else
  showonlyfirstaxis = 'yes';
end

invalid_stats = 0;

%% Set up for plotting.

% Load a layout or montage File

if ~autolayout
  [dum1,dum2,ext] = fileparts(opt.layoutfile);
  clear dum1 dum2
  ext = lower(ext);
  switch (ext)
    case '.lay'
      figs = ts_iEEG_readFTlayout(opt.layoutfile,data.sensor_info);
    case '.asc'
      figs = ts_iEEG_readNSmontage(opt.layoutfile,data.sensor_info);
    otherwise
      fprintf('The layout file must be a valid NeuroScan .asc montage file or a FieldTrip .lay file.\n');
      autolayout = 1;
  end  
  
end

% Use sensor info to determine number and layout of figures.
if autolayout
  if length(channels) == 1                                                   % set up for single plot
    figs.plots.row   = 1;
    figs.plots.column= 1;
    [figs.plots.name,figs.plots.index] = intersect({data.sensor_info.label},channels);
    figs.name        = channels{1};                                        % labels figure with channel name as opposed to page #
  else                                                                       % set up plotting for multiple channels
    if ~isempty(channels)
      ignore_channels = [setdiff({data.sensor_info.label},channels)...     % if given list of channels set up ignore_channels to those not given
        ignore_channels];                                 % and preserve any bad channels
    end
    figs = ts_iEEG_makefigs (data.sensor_info,max_rows,max_columns,max_group_size,ignore_channels,laminar);
    if ~isempty(channels), figs.name = cell2mat(channels); end
    if ~isempty(whichfigs), figs = figs(whichfigs); end
  end

  % Set up some constants for subplot command

  for i = 1:length(figs)
    figs(i).max_rows    = 0;
    figs(i).max_columns = 0;
    for j=1:length(figs(i).plots)
      figs(i).max_rows    = max(figs(i).max_rows,figs(i).plots(j).row);       % adjust size based on figure
      figs(i).max_columns = max(figs(i).max_columns,figs(i).plots(j).column);
    end

    figs(i).width = ((.85-(figs(i).max_columns*columnpad))/figs(i).max_columns);  % allowable width for each plot
    % this accounts for padding between plots
    figs(i).height = ((.9-(figs(i).max_rows*rowpad))/figs(i).max_rows);           % allowable height for each plot with a
    % pad between rows

    % Prevent wierd aspect ratios.
    if figs(i).width > (1.75*figs(i).height)
      figs(i).width  = figs(i).width * .75;
    end
    if figs(i).height > (1.75*figs(i).width)
      figs(i).height = figs(i).height *.75;
    end

    figs(i).ypos = [];
      for j = 1:figs(i).max_rows
        figs(i).ypos = cat(1,figs(i).ypos,(0.9-(j*figs(i).height))-(rowpad*(j-1)));
      end                                                                           
    figs(i).ypos = repmat(figs(i).ypos,[1 figs(i).max_columns]);

    figs(i).xpos = [];
    for j = 1:figs(i).max_columns
      figs(i).xpos=cat(2,figs(i).xpos,((j-1)*figs(i).width)+(columnpad*(j-1))+.05); % where to place each plot on the row
    end
    figs(i).xpos = repmat(figs(i).xpos,[figs(i).max_rows 1]);

    figs(i).width  = repmat(figs(i).width,[figs(i).max_rows figs(i).max_columns]);
    figs(i).height = repmat(figs(i).height,[figs(i).max_rows figs(i).max_columns]);

  end
end
%% Plotting
scrsz = get(0,'ScreenSize');                                               % initial figure size based on screen size
plot_colors = strrep(plot_colors,fig_background,'');                       % remove background color from plot colors
bctag = '';
zscorescaling = 1;
if baseline
    fprintf('Subtracting baseline...\n');
    for j=1:length(conds)
        if isempty(baseline_start), baseline_start = data.(datafield)(conds(j)).time(1); end
           if isempty(baseline_end),   baseline_end   = 0; end
        if isempty(baseline_data)
           fprintf('%s: Calculating baseline for individual event %d.\n',mfilename,data.(datafield)(conds(j)).event_code)
           samples = find(data.(datafield)(conds(j)).time >= baseline_start & data.(datafield)(conds(j)).time <= baseline_end);
           mean_baseline = repmat(nanmean(data.(datafield)(conds(j)).power(:,samples,:),2),...
                                  [1 size(data.(datafield)(conds(j)).data,2) 1]);
           std_baseline  = repmat(nanstd (data.(datafield)(conds(j)).power(:,samples,:),[],2),...
                                 [1 size(data.(datafield)(conds(j)).data,2) 1]);
           blevents = data.(datafield)(conds(j)).event_code;
           zscorescaling = 1;
        else
           fprintf('%s: Using previously calculated baseline.\n',mfilename);
           blsamples      = find(baseline_data.timefreq.time >= baseline_start & baseline_data.timefreq.time <= baseline_end); 
           mean_baseline  = repmat(nanmean(baseline_data.timefreq.power(:,blsamples,:),2),   [1 size(data.(datafield)(conds(j)).data,2) 1]);
           std_baseline   = repmat(nanstd (baseline_data.timefreq.power(:,blsamples,:),[],2),[1 size(data.(datafield)(conds(j)).data,2) 1]);
           baseline_start = baseline_data.timefreq.time(blsamples(1));
           baseline_end   = baseline_data.timefreq.time(blsamples(end));     
           if length(data.(datafield)(conds(j)).num_trials)==1
               if baseline == 4
                   zscorescaling  = (1/sqrt(data.(datafield)(conds(j)).num_trials/baseline_data.timefreq.num_trials));
               else
                   zscorescaling  = 1;
               end
           else 
            for ch = 1:length(baseline_data.timefreq.num_trials)
              if baseline ==4
                   zscorescaling(ch) = (1/sqrt(data.(datafield)(conds(j)).num_trials(ch)/baseline_data.timefreq.num_trials(ch)));
              else
                   zscorescaling(ch) = 1;
              end
            end
           end          
        end
        switch baseline
            case 1 % absolute
                foot_notes = sprintf('%s\nBaseline: Absolute\nStart: %s\nBaeline End: %s',foot_notes,num2str(baseline_start),num2str(baseline_end));
                bctag = sprintf('_bc-abs%s-%sVSallevents',num2str(baseline_start),num2str(baseline_end));
                data.(datafield)(conds(j)).power = data.(datafield)(conds(j)).power - mean_baseline;
            case 2 % relative
                foot_notes = sprintf('%s\nBaseline: Relative\nStart: %s\nBaeline End: %s',foot_notes,num2str(baseline_start),num2str(baseline_end));
                bctag = sprintf('_bc-rel%s-%sVSallevents',num2str(baseline_start),num2str(baseline_end));
                data.(datafield)(conds(j)).power = data.(datafield)(conds(j)).power./mean_baseline;
            case 3 % relative change
                foot_notes = sprintf('%s\nBaseline: Relative Change\nStart: %s\nBaeline End: %s',foot_notes,num2str(baseline_start),num2str(baseline_end));
                bctag = sprintf('_bc-relchng%s-%sVSallevents',num2str(baseline_start),num2str(baseline_end));
                data.(datafield)(conds(j)).power = (data.(datafield)(conds(j)).power - mean_baseline)./mean_baseline;
            case 4 % z-scores
                foot_notes = sprintf('%s\nZ-scores\nBaseline Start: %s\nBaeline End: %s',foot_notes,num2str(baseline_start),num2str(baseline_end));
                bctag = sprintf('_zscores%s-%sVSallevents',num2str(baseline_start),num2str(baseline_end));
                data.(datafield)(conds(j)).power = (data.(datafield)(conds(j)).power - mean_baseline)./std_baseline;                
        end
        clear mean_baseline std_baseline
    end
else
    foot_notes = sprintf('%s\nBaseline: None',foot_notes);
end
if freqcorr
    for j = 1:length(conds)
        norm_fact = repmat(permute(data.(datafield)(conds(j)).frequencies.^3,[1 3 2]),...
            [size(data.(datafield)(conds(j)).power,1)...
            size(data.(datafield)(conds(j)).power,2)...
            1]...
            );
        data.(datafield)(conds(j)).power = data.(datafield)(conds(j)).power.*norm_fact;
        clear norm_fact
    end
end
if calc_ylim
    for k=1:length(conds)
        ymin = min([ymin data.(datafield)(conds(k)).frequencies]);
        ymax = max([ymax data.(datafield)(conds(k)).frequencies]);
    end
end
if timedsfact ~=0
    timedstag = sprintf('_tds%s',num2str(timedsfact));
    foot_notes = sprintf('%s\n Time Downsample: %s',foot_notes,num2str(timedsfact));
else
    timedstag = '';
end
if freqdsfact ~=0
    freqdstag = sprintf('_fds%s',num2str(freqdsfact));
    foot_notes = sprintf('%s\n Freq Downsample: %s',foot_notes,num2str(freqdsfact));
else
    freqdstag = '';
end
for i = 1:length(figs)
   if calc_zlim || abs_zlim
       for j = 1:length(figs(i).plots)                                                                                % channels on this figure
           for k=1:length(conds)                                                                                      % conditions - always 1 for timefreq
               time_index = find(data.(datafield)(conds(k)).time >= xmin & data.(datafield)(conds(k)).time <= xmax); % should be same but this is more flexible
               freq_index = find(data.(datafield)(conds(k)).frequencies >= ymin & data.(datafield)(conds(k)).frequencies <= ymax);       
               for l = 1:length(freq_index)                                                                            % selected frequency range
                    zmin = min([zmin data.(datafield)(conds(k)).power(figs(i).plots(j).index,time_index,freq_index(l))]);% look for z min/max across channels from this page
                    zmax = max([zmax data.(datafield)(conds(k)).power(figs(i).plots(j).index,time_index,freq_index(l))]);% selected frequency range and all conditions
               end
           end
       end     
       if abs_zlim
           abs_max = max([abs(zmin) abs(zmax)]);
           zmax = abs_max;
           zmin = -abs_max;
       end       
       if baseline == 4
           zmin = -3;
           zmax = 3;
       end
   end
   
   limits_tag = sprintf('_freq%s-%s_sec%s-%s',num2str(ymin),num2str(ymax),num2str(xmin),num2str(xmax));
   if ~isempty(fig_name)
       this_figures_name = sprintf('%s%s%s%s%s_%s',fig_name,limits_tag,timedstag,freqdstag,bctag,figs(i).name);
   else
       this_figures_name = 'wavelet_plot';
       for j=i:length(conds)
           this_figures_name = [this_figures_name '_cond' num2str(conds(j))];
       end
       this_figures_name = sprintf('%s%s%s%s%s_%s',this_figures_name,limits_tag,timedstag,freqdstag,bctag,figs(i).name);
   end

   figure('Name',this_figures_name,...                                    % set up figure window
       'NumberTitle','off',...
       'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2],...
       'Color',fig_background);                                           % background color
   for j = 1:length(figs(i).plots)
       subplot('Position',...
           [figs(i).xpos(figs(i).plots(j).row,figs(i).plots(j).column) ...
           figs(i).ypos(figs(i).plots(j).row,figs(i).plots(j).column)...
           figs(i).width(figs(i).plots(j).row,figs(i).plots(j).column) ...
           figs(i).height(figs(i).plots(j).row,figs(i).plots(j).column)]);
       hold on                                                            % only add to plot, make sure nothing is replaced
       for k=1:length(conds)
           if data.sensor_info(figs(i).plots(j).index).badchan ~= 1             
               TFP = squeeze(data.(datafield)(conds(k)).power(figs(i).plots(j).index,:,:));
               %TFP(find(isnan(TFP)))  = 0;
               timevec = data.(datafield)(conds(k)).time;
               freqvec = data.(datafield)(conds(k)).frequencies;
               if timedsfact ~= 0
                   TFP = resample(TFP,1,timedsfact,0);
                   timevec = resample(timevec,1,timedsfact,0);
               end
               if freqdsfact ~= 0
                   TFP = permute(TFP,[2 1]); % work on freq dim
                   TFP = resample(TFP,1,freqdsfact,0);
                   TFP = permute(TFP,[2 1]);
                   freqvec = resample(freqvec,1,freqdsfact,0);
               end 
               if ~clear_figs 
                 % takes less resources so better when staying on screen to use imagesc
                 imagesc(timevec,freqvec,TFP');
               else
                 % looks better so when saving files with clear screen use pcolor
                 pcolor(timevec,freqvec,TFP');                 
                 shading(shading_style);
               end
           end
           
           box off;                                                           % remove MATLAB box from plot
           xlim([xmin xmax]);                                                 % set up time axis limits
           ylim([ymin ymax]);                                                 % set up freq axis limits
           if ~isempty(zmin) && ~isempty(zmax)                                % set limits for color axis or auto
              if length(zscorescaling) == 1   
               zminch = zmin * zscorescaling;
               zmaxch = zmax * zscorescaling;
               caxis([zminch zmaxch]);
               title(figs(i).plots(j).name,'FontSize',8,...
               'VerticalAlignment','Middle','Clipping','off');                % show channel name
              else
               zminch = zmin * zscorescaling(figs(i).plots(j).index);
               zmaxch = zmax * zscorescaling(figs(i).plots(j).index);
               caxis([zmin zmax]);
               title(sprintf('%s - [%0.1f %0.1f]\nTrls: %.0f BLTrls: %.0f',...
                     char(figs(i).plots(j).name),zminch,zmaxch,...
                     data.(datafield)(conds(k)).num_trials(figs(i).plots(j).index),...
                     baseline_data.timefreq.num_trials(figs(i).plots(j).index)),...
                   'FontSize',6,...
                   'VerticalAlignment','Middle','Clipping','off');                % show channel name with ranges
               show_colorbar = false;
              end
           end      
           
           
           if show_axis  && data.sensor_info(figs(i).plots(j).index).badchan ~= 1 % display MATLAB axis
               if strcmpi(showonlyfirstaxis,'no')
                   axis on;
                   set(gca,'FontSize',axis_font);
                   x_tick_markers = get(gca,'xtick');
                   y_tick_markers = get(gca,'ytick');
               else
                   if j==1
                       axis on;
                       set(gca,'FontSize',axis_font);
                       x_tick_markers = get(gca,'xtick');
                       y_tick_markers = get(gca,'ytick');
                   else
                       axis off;
                   end
               end
           else
               axis off;
               x_tick_markers = get(gca,'xtick');
               y_tick_markers = get(gca,'ytick');
           end
           if center_axis && data.sensor_info(figs(i).plots(j).index).badchan ~= 1                                         % display manually created zero axis
               if center_axis_x_ticks, plot(repmat(x_tick_markers',2)',repmat(ylim,length(x_tick_markers)*2,1)',...        % tick lines on zero axis
                       'color','white','linewidth',zero_line_width); end
               plot(zeros(2),ylim,line_style,'color','black','linewidth',zero_line_width);
               if center_axis_y_ticks, plot(repmat(y_tick_markers',2)',repmat(xlim,length(y_tick_markers)*2,1)',...        % tick lines on zero axis
                       'color','white','linewidth',zero_line_width); end
           end
           new_foot_notes = sprintf('%s\n TOI: %s to %s',foot_notes,...
               num2str(data.(datafield)(conds(k)).time(1)),...
               num2str(data.(datafield)(conds(k)).time(end)));
           new_foot_notes = sprintf('%s\n FOI: %s to %s',new_foot_notes,...
               num2str(data.(datafield)(conds(k)).frequencies(1)),...
               num2str(data.(datafield)(conds(k)).frequencies(end)));
       end
   end
   if show_colorbar, colorbar('Position',[.01 .1 .025 .75],'FontSize',6); end;
   annotation('textbox','Position',[0 .9 1 .1],...                        % set up figure's title box
       'VerticalAlignment','middle',...
       'HorizontalAlignment','center',...
       'Color','k',...
       'FontSize',10,...
       'FontWeight','bold',...
       'FitHeightToText','on',...
       'LineStyle','none',...
       'String',title_string);
   annotation('textbox','Position',[.9 .7 .1 .1],...                     % set up figure's footnotes
       'VerticalAlignment','top',...
       'HorizontalAlignment','center',...
       'Color','k',...
       'FontSize',6,...
       'FitHeightToText','on',...
       'LineStyle','none',...
       'String',new_foot_notes);
   hold off
   if save_figs
       set(gcf,'PaperOrientation',paperorient,...                         % set paper properties first in order to get
           'PaperType',papertype,...                                      % proper paper size for figure
           'PaperUnits','inches');
       papersize = get(gcf,'PaperSize');
       fig_width = papersize(1) - leftmargin - rightmargin;               % fill page
       fig_height= papersize(2) - topmargin - bottommargin;
       fig_left  = (papersize(1)-fig_width)/2;                            % center figure
       fig_bottom= (papersize(2)-fig_height)/2;
       set(gcf,'PaperPosition', [fig_left, fig_bottom, fig_width, fig_height]);
       this_figures_name = [this_figures_name '.' opt.savefigs];
       if isfield(opt,'figpath'), this_figures_dir = opt.figpath;
       else this_figures_dir  = fullfile(pwd,'images');
       end
       if ~exist(this_figures_dir,'dir'), mkdir(this_figures_dir); end
       this_figures_name = fullfile(this_figures_dir,this_figures_name);
       fprintf('Saving figure: %s\n',this_figures_name);
       print(print_command{:}, this_figures_name);              % print to an eps file
   end
   if clear_figs, close; end                                              % don't keep figures on screen option
end
