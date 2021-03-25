function mmil_barh(data,varargin)
%function mmil_barh(data,[options])
%
% Purpose: 
%   creates horizontal bar plot with optional labels and error bars
%
% Usage:
%   mmil_barh(data,'key1', value1,...);
%
% Required input:
%  data: data matrix
%     rows should be different measurements
%     columns should be different groups
%
% Optional parameters:
%  'err': error matrix for error bars
%    should have same size as data matrix
%    {default = []}
%  'data_names': cell array of labels for different measurements
%    length should match number of rows of data matrix
%    {default = []}
%  'group_names': cell array of labels for different groups
%    length should match number of columns of data matrix
%    if ommitted, no legend will be displayed
%    {default = []}
%  'barwidth': relative width of bars
%    {default = 1}
%  'title': plot title
%    {default = []}
%  'xlabel': x axis label
%    {default = []}
%  'xlim': x axis limits
%    {default = []}
%  'fontname'
%    {default = 'Arial'}
%  'fontsize'
%    {default = 12}
%  'linewidth': line width for error bars
%    {default = 1.5}
%  'color_order': cell array containing matlab color characters
%    {default = {'b','g','r','c','m','y','k'}}
%  'legend_loc': legend location
%    {default = 'NorthEast'}
%
%  Created:  11/29/07 by Don Hagler
%  Last Mod: 05/18/15 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~mmil_check_nargs(nargin,1), return; end;
parms = mmil_args2parms(varargin,{...
  'err',[],[],...
  'data_names',[],[],...
  'group_names',[],[],...
  'barwidth',1,[0.1,2],...
  'title',[],[],...
  'xlabel',[],[],...
  'xlim',[],[],...
  'fontname','Arial',[],...
  'fontsize',12,[1,100],...
  'linewidth',1.5,[0.1,100],...
  'color_order',{'b','g','r','c','m','y','k'},[],...
  'legend_loc','NorthEast',[],...
...
  'text_v_offset',0,[-1,1],...
  'text_h_offset',0.01,[-1,1],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

num_meas = size(data,1);
num_groups = size(data,2);

% empirically determined bar spacings
if num_groups==1
  bar_spacing = 1;
elseif num_groups==2
  bar_spacing = 0.26;
elseif num_groups==3
  bar_spacing = 0.22;
elseif num_groups==4
  bar_spacing = 0.17;
elseif num_groups==5
  bar_spacing = 0.15;
else
  bar_spacing = 0.13;
end;
errwidth = bar_spacing/2;

%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

if ~isempty(parms.data_names)
  if ~iscell(parms.data_names), parms.data_names = {parms.data_names}; end;
  if length(parms.data_names) ~= num_meas
    error('length of data_names must match number of rows of data');
  end;
end;

if ~isempty(parms.group_names)
  if ~iscell(parms.group_names), parms.group_names = {parms.group_names}; end;
  if length(parms.group_names) ~= num_groups
    error('length of group_names must match number of columns of data');
  end;
end;

if ~isempty(parms.err) && any(size(data)~=size(parms.err))
  error('size of err must match size of data');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%

clf;
if num_meas==1
  tmp_data = cat(1,data,zeros(size(data)));
else
  tmp_data = data;
end;
barh(tmp_data,parms.barwidth);

h_ax_data = gca;
axis ij;
set(h_ax_data,'Color','none')
set(h_ax_data,'FontSize',parms.fontsize,'FontName',parms.fontname);
if ~isempty(parms.title), title(parms.title); end;
if ~isempty(parms.xlabel), xlabel(parms.xlabel); end;
if ~isempty(parms.xlim), set(h_ax_data,'XLim',parms.xlim); end;
set(h_ax_data,'YLim',[0,num_meas+1]);

if ~isempty(parms.color_order)
  if ~iscell(parms.color_order), parms.color_order = {parms.color_order}; end;
  while length(parms.color_order)<num_groups
    parms.color_order = {parms.color_order;parms.color_order};
  end;
  if length(parms.color_order)>num_groups
    parms.color_order = {parms.color_order{1:num_groups}};
  end;
  cmap = [];
  for i=1:length(parms.color_order)
    switch parms.color_order{i}
      case 'y'
        rgbvec = [1 1 0];
      case 'm'
        rgbvec = [1 0 1];
      case 'c'
        rgbvec = [0 1 1];
      case 'r'
        rgbvec = [1 0 0];
      case 'g'
        rgbvec = [0 1 0];
      case 'b'
        rgbvec = [0 0 1];
      case 'w'
        rgbvec = [1 1 1];
      case 'k'
        rgbvec = [0 0 0];
      otherwise
        error('bad color specifier %s',parms.color_order{i});
    end;
    cmap(i,:) = rgbvec;
  end;
  colormap(cmap);
end;

if ~isempty(parms.group_names)
  legend(parms.group_names,'Location',parms.legend_loc);
end;

max_text_width = 0;
total_text_height = 0;
pos = get(h_ax_data,'Position');
if isempty(parms.data_names)
  set(h_ax_data,'YTick',[0.5:1:num_meas+0.5]);
  set(h_ax_data,'YTickLabel',[]);
else
  set(h_ax_data,'YTick',[0.5:1:num_meas+0.5]);
  set(h_ax_data,'YTickLabel',[]);
  h_ax_text = axes('Position',[0 0 1 1]);
  axis ij;
  set(h_ax_text,'Visible','off')
  x = pos(1) - parms.text_h_offset;
  h_text = zeros(num_meas,1);
  y0 = 1 - (pos(2) + pos(4));
  yh = pos(4);
  ys = yh / (num_meas+1);
  y1 = parms.text_v_offset + y0 + ys;

  % draw text labels
  for i=1:num_meas
    if num_meas>1
      y = y1 + ((i-1)/(num_meas-1))*(yh-2*ys);
    else
      y = y1;
    end;
    h_text(i) = text(x,y,parms.data_names{i},...
      'HorizontalAlignment','right',...
      'FontSize',parms.fontsize,'FontName',parms.fontname);

    % get width of text labels
    tmp = get(h_text(i));
    text_width = tmp.Extent(3);
    if text_width > max_text_width
      max_text_width = text_width;
    end;
    text_height = tmp.Extent(4);
    total_text_height = total_text_height + text_height;
  end;

  % adjust width of plot to compensate for width of text labels
  text_h_space = max_text_width + 2*parms.text_h_offset;
  pos(3) = max([pos(3) + pos(1) - text_h_space,0.1]);
  pos(1) = text_h_space;
  set(h_ax_data,'Position',pos);
  x = pos(1) - parms.text_h_offset;

  % reposition text labels
  for i=1:num_meas
    if num_meas>1
      y = y1 + ((i-1)/(num_meas-1))*(yh-2*ys);
    else
      y = y1;
    end;
    set(h_text(i),'Position',[x,y,0]);
  end;
end;

if ~isempty(parms.err)
  h_ax_err = axes('Position',pos);
  axis ij;
  set(h_ax_err,'Visible','off','Color','none')
  hold(h_ax_err,'on');
  for m=1:num_meas
    for g=1:num_groups
      colstr = parms.color_order{g};
      tmp_data = data(m,g);
      tmp_err = parms.err(m,g);
      if tmp_data>=0
        x = [tmp_data,tmp_data+tmp_err];
      else
        x = [tmp_data,tmp_data-tmp_err];
      end;
      if mod(num_groups,2) % odd
        tmp_y = m - ((num_groups-1)/2)*bar_spacing +...
                    (g-1)*bar_spacing;
      else % even
        tmp_y = m - num_groups/2*bar_spacing +...
                    (g-1)*bar_spacing + bar_spacing/2;
      end;
      y = [tmp_y,tmp_y];
      line(x,y,'LineWidth',parms.linewidth,'Color',colstr);
      if tmp_data>=0
        x = [tmp_data+tmp_err,tmp_data+tmp_err];
      else
        x = [tmp_data-tmp_err,tmp_data-tmp_err];
      end;
      y = [tmp_y-errwidth,tmp_y+errwidth];
      line(x,y,'LineWidth',parms.linewidth,'Color',colstr);
    end;
  end;
  set(h_ax_err,'XLim',get(h_ax_data,'XLim'),'Layer','top')
  set(h_ax_err,'YLim',get(h_ax_data,'YLim'),'Layer','top')
  set(h_ax_err,'Visible','off','Color','none')
end;

axes(h_ax_data);

