function [layfile tempflag] = write_grid_layout(filepath,sens,maxrows,maxcols,padfactor,infile)
if ~exist('padfactor','var'), padfactor = .1; end
if ~exist('infile','var'), infile = ''; end
%	write_grid_layout(layout1,data.sensor_info,chans,8,8)
tempflag = 1;
if ~exist(filepath,'file')
  fprintf('%s: making output directory: %s\n',mfilename,filepath);
  mkdir(filepath);
end

% create fig structure
n = 1;
if isempty(infile) || ~exist(infile,'file')
  ignore_channels = [];
  max_group_size  = length(sens) - length(ignore_channels);
  fig = ts_iEEG_makefigs (sens,maxrows,maxcols,max_group_size,ignore_channels,0);

  nr = []; nc = [];
  for i = 1:length(fig)
    if isempty(fig(i).plots), continue; end
    nr = max([fig(i).plots.row]);
    nc = max([fig(i).plots.column]);
    if isempty(nr), continue; end
    if isempty(nc), continue; end

    nr = max(nr,nc);
    nc = max(nr,nc);

    nr = nr + 1;
    nc = nc + 1;

    w = 0.8 * 1/nc;
    h = 0.8 * 1/nr;

    rpad = padfactor/sqrt(nr); % .1/sqrt(nr)
    cpad = padfactor/sqrt(nc); % .1/sqrt(nc)

    h = ((1-padfactor) - nr*rpad) / nr; % (.9 - nr*rpad)/nr
    w = ((1-1.5*padfactor) - nc*cpad) / nc; % (.85 - nc*cpad)/nc
    h = max(w,h);
    w = h;

    k = 1:nr;  y = [(1-padfactor) - (k-1)*rpad - k*h]; % [.9 - (k-1)...			% y pos
    k = 1:nc;  x = [(.5*padfactor) + (k-1)*cpad + (k-1)*w]; % .05 + (k-1)...	% x pos
    w = w*ones(1,nc);
    h = h*ones(1,nr);

    layfile{n} = fullfile(filepath,sprintf('layout_page%g.lay',n));%fig(i).name)); % sprintf('%s/layout_%s.lay',filepath,fig(i).name);
    if exist(layfile{n},'file'), tempflag = 0; end
    fid = fopen(layfile{n},'wt');
    for j = 1:length(fig(i).plots)
  % 		if any(bix==fig(i).plots(j).index), continue; end
      r 	= fig(i).plots(j).row;
      c 	= fig(i).plots(j).column;
      id 	= fig(i).plots(j).index;
      lbl = fig(i).plots(j).name;
      fprintf(fid,'%g %.6f %.6f %.6f %.6f %s\n',id,x(c),y(r),w(c),h(r),lbl);% id x y width height label
    end
  %   fprintf(fid,'%g %.6f %.6f %.6f %.6f %s\n',max(id)+1,max(x(1,:))+mean(diff(x(1,:))),max(y(:,1)),w,h,'COMNT');
    fprintf(fid,'%g %.6f %.6f %.6f %.6f %s\n',max([fig(i).plots(j).index])+1,(nc-1)/nc,max((nr-2)/nr,0),w(1),h(1),'COMNT');
    fclose(fid);
    clear w h x y k
    n = n + 1;
  end
elseif findstr(infile,'.asc')
  scalefactor = 1-padfactor;
  fig = ts_iEEG_readNSmontage(infile,sens);
  for i = 1:length(fig)
    if isempty(fig(i).plots), continue; end
    x   = fig(i).xpos   * scalefactor + (1-scalefactor)/2;
    y   = fig(i).ypos   * scalefactor + (1-scalefactor)/2;
    w   = fig(i).width  * scalefactor;
    h   = fig(i).height * scalefactor;
    
    layfile{n} = fullfile(filepath,sprintf('layout_page%g.lay',n));%fig(i).name)); % sprintf('%s/layout_%s.lay',filepath,fig(i).name);
    if exist(layfile{n},'file'), tempflag = 0; end
    fid = fopen(layfile{n},'wt');
    for j = 1:length(fig(i).plots)
      id 	= fig(i).plots(j).index;
      lbl = fig(i).plots(j).name;
      fprintf(fid,'%g %.6f %.6f %.6f %.6f %s\n',id,x(j),y(j),w(j),h(j),lbl);% id x y width height label
    end
    fprintf(fid,'%g %.6f %.6f %.6f %.6f %s\n',max([fig(i).plots(j).index])+1,max(x)+w(1),max(max(y)+h(1),0),w(1),h(1),'COMNT');
    fclose(fid);
    clear w h x y k
    n = n + 1;
  end
end


% % % generate ordered layout (code copied from ts_gui)
% layfile{end+1} = fullfile(filepath,sprintf('layout_ordered.lay'));
% if exist(layfile{end},'file'), tempflag = 0; end
% fid = fopen(layfile{end},'wt');
% 
%   nchan = length({sens.label});
%   ncol = ceil(sqrt(nchan))+1;
%   nrow = ceil(sqrt(nchan))+1;
%   k = 0;
%   for i=1:nrow
%     for j=1:ncol
%       k = k+1;
%       if k<=nchan
%         x = (j-1)/ncol;
%         y = (nrow-i-1)/nrow;
%         lay.pos(k,:) = [x y];
%         lay.width(k,1)  = 0.8 * 1/ncol;
%         lay.height(k,1) = 0.8 * 1/nrow;
%       end
%     end
%   end
%   lay.label = {sens.label};
% 
%   lay.label{end+1}  = 'SCALE';
%   lay.width(end+1)  = mean(lay.width);
%   lay.height(end+1) = mean(lay.height);
%   x = (ncol-2)/ncol;
%   y = 0/nrow;
%   lay.pos(end+1,:) = [x y];
% 
%   lay.label{end+1}  = 'COMNT';
%   lay.width(end+1)  = mean(lay.width);
%   lay.height(end+1) = mean(lay.height);
%   x = (ncol-1)/ncol;
%   y = 0/nrow;
%   lay.pos(end+1,:) = [x y];  
%   
%   for k = 1:length(lay.label)
%     fprintf(fid,'%g %f %f %f %f %s\n',k,lay.pos(k,1),lay.pos(k,2),lay.width(k),lay.height(k),lay.label{k});
%   end
%   
% fclose(fid);
