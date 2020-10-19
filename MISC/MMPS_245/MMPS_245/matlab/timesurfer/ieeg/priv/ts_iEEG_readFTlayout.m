function fig = ts_iEEG_readFTlayout(filename,sensor_info)
% Use: fig = ts_iEEG_readFTlayout (filename,sensor_info)
%
% Reads a layout file created for FieldTrip and makes the appropriate fig
% structure
%
% Required Input
%
% filename - valid FieldTrip Layout File
% sensor_info - the sensor_info field of a valid TimeSurfer data structure
%
% Ouput:
%
%    fig - an array the size of number of fiures to be created 
%           containing the following field
%     plots - an array the size of the number of plots for this figure
%             containing the following fields:
%         row    - which row of the figure
%         column - which column of the figure
%         index  - data index of the channel
%         name   - the channel name for this plot
%
% Created by Rajan Patel 01/17/2008
% Last Modified by Rajan Patel 01/17/2008
%
% See FieldTrip for more info on layout files:
%
% prepare_layout

%% Correct MEG labels if they lack spaces
label = {sensor_info.label};
ix = find(cellfun('isempty',regexp(label,'MEG\s+\d+')));
for k = 1:length(ix)
  sensor_info(ix(k)).label = strrep(sensor_info(ix(k)).label,'MEG','MEG ');
end
clear label ix

%% Read in Layout File (taken from FieldTrip code)

if ~exist(filename, 'file')
  error(sprintf('could not open layout file: %s', filename));
end
[chNum,X,Y,Width,Height,Lbl,Rem] = textread(filename,'%f %f %f %f %f %q %q');
for i=1:length(Lbl)
  if ~isempty(Rem{i})
    % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
    Lbl{i} = [Lbl{i} ' ' Rem{i}];
  end
end
lay.pos    = [X Y];
lay.width  = Width;
lay.height = Height;
lay.label  = Lbl;

%% Convert to Subplot Reference Units [0 1]


newwidth  =fixpt_interp1(linspace(0,max(abs(lay.pos(:,1))),size(lay.pos,1)),...
                         linspace(0,1,size(lay.pos,1)),...
                         lay.width,...
                         float('single'),1^-14,...
                         float('single'),1^-14,...
                         'Floor');
newheight =fixpt_interp1(linspace(0,max(abs(lay.pos(:,2))),size(lay.pos,1)),...
                         linspace(0,1,size(lay.pos,1)),...
                         lay.height,...
                         float('single'),1^-14,...
                         float('single'),1^-14,...
                         'Floor');
                       
max_right = 1-(max(newwidth)/2);
max_top   = 1-(max(newheight)/2);

newpos(:,1) = fixpt_interp1(linspace(min(lay.pos(:,1)),max(lay.pos(:,1)),size(lay.pos,1)),...
                            linspace(0,max_right,size(lay.pos,1)),...
                            lay.pos(:,1),...
                            float('single'),1^-14,...
                            float('single'),1^-14,...
                            'Floor');
newpos(:,2) = fixpt_interp1(linspace(min(lay.pos(:,2)),max(lay.pos(:,2)),size(lay.pos,1)),...
                            linspace(0,max_top,size(lay.pos,1)),...
                            lay.pos(:,2),...
                            float('single'),1^-14,...
                            float('single'),1^-14,...
                            'Floor');
%% Test

% for i = 1:length(lay.label)
%   sensor_info(i).label = lay.label(i);
% end                                  
%% Make it Into fig Structure

% fig.name = 'page1';
[idxlay jnk] = match_str(lay.label,{sensor_info.label});
lay.label = lay.label(idxlay);

fig = [];
fig.plots = [];
fig.name  = 'page1';
plotnum = 1;
for i = 1:length(lay.label)
 fig.plots(plotnum).index = [];
 for c = 1:length(sensor_info)
   if strcmp(lay.label(i),sensor_info(c).label)
     fig.plots(plotnum).index = c;
     fig.plots(plotnum).name  = sensor_info(c).label;
   end
 end
 if ~isempty(fig.plots(plotnum).index)
   fig.plots(plotnum).column = plotnum;
   fig.plots(plotnum).row    = 1;
   fig.xpos  (fig.plots(plotnum).row,fig.plots(plotnum).column) = newpos(i,1);
   fig.ypos  (fig.plots(plotnum).row,fig.plots(plotnum).column) = newpos(i,2);
   fig.width (fig.plots(plotnum).row,fig.plots(plotnum).column) = newwidth(i)/2;
   fig.height(fig.plots(plotnum).row,fig.plots(plotnum).column) = newheight(i)/2;
   plotnum = plotnum + 1;
 else
   fprintf('WARNING: Channel %s was not found in the data set.\n',num2str(i),lay.label(i));
 end
end
  
%% Test 

% for i = 1:length(fig.plots)
%       hold on;
%       subplot('Position',[fig.xpos(fig.plots(i).row,fig.plots(i).column) fig.ypos(fig.plots(i).row,fig.plots(i).column) fig.width(fig.plots(i).row,fig.plots(i).column) fig.height(fig.plots(i).row,fig.plots(i).column)]);
%       title(fig.plots(i).name,'fontsize',4);
%       axis off;
% end