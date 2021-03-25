function fig = ts_iEEG_readNSmontage (filename,sensor_info)
% Use: fig = ts_iEEG_readNSmontage (filename,sensor_info)
%
% Reads montage file (channel layout) exported by Neuroscan 
%
% Required Input:
%
% filename - valid Neuroscan asc montage file
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
% Based on function written by Chunmao Wang, 5/7/2004
% Modified for iEEG Stream by Rajan Patel 01/16/2008
%
% See Neuroscan for more information on the montage file.

%% Mao's Read Montage Scripts

fid=fopen(filename);
i=0;j=0;k=0;
while 1
   line = fgetl(fid);
   if ~isstr(line), break, end
   
   % read channel names
   if line(1) == '#'
      line2=line(2:size(line,2));
      j=sscanf(line2,'%d');
      line3=line2(6:size(line2,2));
      s=char(zeros(1,10));
      t=sscanf(line3,'%s');
      s(1:size(t,2))=t;
      ChanName(j(1),:)=s;
   end
   
	% read channel coordinates   
   
   if and( line(1) ~= ';', line(1) ~= '#' )
      k=k+1;
      t=sscanf(line,'%f');
      ChanLayout(k,:,:,:,:,:,:)=t;
   end
   
   %disp(line)
   i=i+1;
end
fclose(fid);

%% Reconfigure Mao's Structures

for i = 1:max(ChanLayout(:,1)+1)
  layouts(i).layout  = ChanLayout(find(ChanLayout(:,1)==i-1),3:end);
  layouts(i).channel = ChanLayout(find(ChanLayout(:,1)==i-1),2);
end
clear ChanLayout

ChanName = cellstr(ChanName);

%% Make it into a fig structure

fig = [];
for i = 1:length(layouts);
 fig(i).plots  = [];
 fig(i).name = sprintf('page%s',num2str(i));
 plotnum = 1;
 for j = 1:length(layouts(i).channel)
   fig(i).plots(plotnum).index = [];
   for c = 1:length(sensor_info)
     if strcmp(ChanName{layouts(i).channel(j)},sensor_info(c).label)
        fig(i).plots(plotnum).index = c;
        fig(i).plots(plotnum).name  = sensor_info(c).label;
     end
   end
   if ~isempty(fig(i).plots(plotnum).index)
     fig(i).plots(plotnum).column = plotnum;
     fig(i).plots(plotnum).row    = 1;
     fig(i).xpos  (fig(i).plots(plotnum).row,fig(i).plots(plotnum).column) = layouts(i).layout(j,1);
     fig(i).ypos  (fig(i).plots(plotnum).row,fig(i).plots(plotnum).column) = (1-layouts(i).layout(j,2))-(layouts(i).layout(j,4)-layouts(i).layout(j,2));
     fig(i).width (fig(i).plots(plotnum).row,fig(i).plots(plotnum).column) = layouts(i).layout(j,3) - layouts(i).layout(j,1);
     fig(i).height(fig(i).plots(plotnum).row,fig(i).plots(plotnum).column) = layouts(i).layout(j,4) - layouts(i).layout(j,2);
     plotnum = plotnum+1;
   else
     fprintf('WARNING: For page %s, channel %s was not found in the data set.\n',num2str(i),ChanName{layouts(i).channel(j)});
   end
 end
end