 function fig = ts_iEEG_makefigs (sens,max_rows,max_columns,max_group_size,ignore_channels,laminar)

%  Use:  ts_iEEG_makefigs (sensor_info,max_rows,max_columns,max_group_size,ignore_channels);
%
%    This will allow the automatic creation of figures for use with
%    iEEG data.  Large grids will also be set to have their own figure.
%    Only accepts TimeSurfer data structure now.
%
%    Used by ts_iEEG_Plot.
%
%  Require input:
%
%    sensor_info    - sensor_info field from TimeSurfer data structure
%    max_rows       - how many rows per figure
%    max_columns    - how many columns per figure
%    max_group_size - max size of channel group requiring its own page
%    ignore_channels- cell array of channels not to be plotted
%
%  Ouput:
%    
%    figs - an array the size of number of fiures to be created 
%           containing the following field
%     plots - an array the size of the number of plots for this figure
%             containing the following fields:
%         row    - which row of the figure
%         column - which column of the figure
%         index  - data index of the channel
%         name   - the channel name for this plot
%
%  Created:       09/25/07 - Rajan Patel
%  Last Modified: 09/27/07 - Rajan Patel
%
%  See also: ts_iEEG_plot

% remove periods
orig_labels = {sens.label};
if any(~cellfun('isempty',regexp({sens.label},'\.')))
  ix = find(~cellfun('isempty',regexp({sens.label},'\.')));
  for i = 1:length(ix)
    sens(ix(i)).label = strrep(sens(ix(i)).label,'.','');
  end
end

% add typestring prefix if label contains digits only
if any(cellfun('isempty',regexp({sens.label},'\D+')))
  ix = find(cellfun('isempty',regexp({sens.label},'\D+')));
  for i = 1:length(ix)
    sens(ix(i)).label = [sens(ix(i)).typestring sens(ix(i)).label];
  end
end
proc_labels = {sens.label};

exlabels = ignore_channels;
[jnk ix] = setdiff({sens.label},exlabels);
sens     = sens(sort(ix));
labels   = {sens.label};

found    = zeros(1,length(labels));
groups   = unique(regexp([labels{:}],'[A-Za-z]*(?=[\d*])','match'));

sz = [];
for k = 1:length(groups), sz(k) = length(groups{k}); end
[jnk ix] = sort(sz,'descend');
groups   = groups(ix);    % sort strings from longest to shortest
grouped  = [];

for k = 1:length(groups)
  grouped(k).name    = groups{k};
  grouped(k).indices = [];  
  labnum = [];
  str    = groups{k};
  pat1   = sprintf('%s(?=\\d)',str);
  pat2   = sprintf('(?<=%s)\\d*',str);
  n      = 1;
  for ch = 1:length(labels)
    % label must not belong to a prior charset
    if found(ch) ~= 0, continue; end    
    lab       = labels{ch};
    [ix1 ix2] = regexp(lab,pat1,'start','end');  % indices of 1st and last character of str in lab
    [jx1 jx2] = regexp(lab,pat2,'start','end');  % indices of 1st and last number following str in lab
    if length(ix1) ~= 1,           continue; end  % label must contain exactly one copy of str
    if isempty(jx1),               continue; end  % str must be followed by a number
    if (ix2-ix1+1) ~= length(str), continue; end  % first & last character of str must be separated by length(str) in lab
    labnum(n)            = str2double(regexp(lab,pat2,'match'));
    grouped(k).indices(n) = ch;
    found(ch)            = 1;
    n = n + 1;
  end
  % sort indices based on the label number that follows this str
  [jnk ix]           = sort(labnum);
  grouped(k).indices = grouped(k).indices(ix);
end
% group remaining channels
ix = find(~found);
if ~isempty(ix)
  grouped(k+1).name    = 'Others';
  grouped(k+1).indices = ix;
end

curr_fig        = 0;
curr_row        = 0;
curr_column     = 1;
count           = 1;
for g = 1:length(grouped)
  if isempty(grouped(g).indices), continue; end
  switch(grouped(g).name)
      case {'G' 'GA' 'GR' 'GB' 'GC'}
          curr_fig           = curr_fig + 1;            
          fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
          % find the number of the last sensor of the group
          % assume it is a square grid
          grid_size          = ceil(sqrt(str2double(...
                            regexp(sens(grouped(g).indices(end)).label,...
                            sprintf('(?<=%s)\\d*',grouped(g).name),...
                            'match'))));
          for ch = 1:length(grouped(g).indices)
            curr_sensor = str2double( ...
                          regexp(sens(grouped(g).indices(ch)).label,...
                          sprintf('(?<=%s)\\d*',grouped(g).name),...
                          'match'));
            fig(curr_fig).plots(ch).name  = sens(grouped(g).indices(ch)).label;
            fig(curr_fig).plots(ch).index = grouped(g).indices(ch); 
            fig(curr_fig).plots(ch).row   = ceil(curr_sensor/grid_size);
            fig(curr_fig).plots(ch).column= rem(curr_sensor,grid_size);
            if fig(curr_fig).plots(ch).column == 0, fig(curr_fig).plots(ch).column = grid_size; end
          end
          curr_fig = curr_fig + 1;
          fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
          curr_row = 0;
          count    = 1;
          continue;
      case 'Others'
          if curr_fig == 0
              curr_fig = 1;
              curr_row = 0;
              count    = 1;
              fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
          else
            curr_fig = curr_fig + 1;
            curr_row = 0;
            count    = 1;
            fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
          end
          curr_row = curr_row + 1;
          curr_column = 1;
          for ch = 1:length(grouped(g).indices)
            if curr_column > max_columns
                curr_column = 1;
                curr_row    = curr_row + 1;
            end
            if curr_row > max_rows
                curr_fig    = curr_fig + 1;
                curr_column = 1;
                curr_row    = 1;
                count       = 1;
                fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
            end  
            fig(curr_fig).plots(count).name   = sens(grouped(g).indices(ch)).label;
            fig(curr_fig).plots(count).index  = grouped(g).indices(ch); 
            fig(curr_fig).plots(count).row    = curr_row;
            fig(curr_fig).plots(count).column = curr_column;
            curr_column = curr_column + 1;
            count = count + 1;
          end         
      otherwise
          if curr_fig == 0
              curr_fig = 1;
              curr_row = 0;
              count    = 1;
              fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
          end
          curr_row = curr_row + 1;
          if curr_row > max_rows
              curr_fig = curr_fig + 1;
              curr_row = 1;
              count    = 1;
              fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
          end
          % find the number of the last sensor of the group
          strip_size = str2double(...
                            regexp(sens(grouped(g).indices(end)).label,...
                            sprintf('(?<=%s)\\d*',grouped(g).name),...
                            'match'));
          no_spaces  = 0;
          if strip_size > max_columns
              no_spaces   = 1;        
              curr_column = 1;
          end
          for ch = 1:length(grouped(g).indices)
            if ~no_spaces
               curr_sensor = str2double(...
                               regexp(sens(grouped(g).indices(ch)).label,...
                               sprintf('(?<=%s)\\d*',grouped(g).name),...
                               'match'));                        
               fig(curr_fig).plots(count).name   = sens(grouped(g).indices(ch)).label;
               fig(curr_fig).plots(count).index  = grouped(g).indices(ch); 
               fig(curr_fig).plots(count).row    = curr_row;
               fig(curr_fig).plots(count).column = rem(curr_sensor,strip_size);
               if fig(curr_fig).plots(count).column == 0 , fig(curr_fig).plots(count).column = strip_size; end
               count = count + 1;
            else  %% NOT IDEAL BUT WORKS.        
               fig(curr_fig).plots(count).name   = sens(grouped(g).indices(ch)).label;
               fig(curr_fig).plots(count).index  = grouped(g).indices(ch); 
               fig(curr_fig).plots(count).row    = curr_row;
               fig(curr_fig).plots(count).column = curr_column; 
               count = count + 1;
               curr_column = curr_column + 1;
               if curr_column > max_columns
                   curr_row = curr_row + 1;
                   if curr_row > max_rows
                      curr_fig = curr_fig + 1;
                      curr_row = 1;
                      count    = 1;
                      fig(curr_fig).name = sprintf('page%s',num2str(curr_fig));
                   end
                   curr_column = 1;
               end
            end
          end          
  end
  if curr_fig~=0 && isempty(fig(curr_fig).plots)
    curr_fig = curr_fig - 1;
  end
end
if ~isequal(orig_labels,proc_labels)
  for i = 1:length(fig)
    for j = 1:length(fig(i).plots)
      ix = find(ismember(proc_labels,fig(i).plots(j).name));
      if numel(ix)==1
        fig(i).plots(j).name = orig_labels{ix};
      end
    end
  end
end


