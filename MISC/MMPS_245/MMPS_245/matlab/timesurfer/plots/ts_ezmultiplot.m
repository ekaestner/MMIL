function [lay,err] = ts_ezmultiplot(data,varargin)
% Purpose:
%   Plot multiple channels of ERP/ERF or power for 1 condition
%
% Inputs: 
%   avg_data, epoch_data or timefreq_data; stat_data
% Other inputs: 
%   xlim: {[begin end],'all'} {default: all}
%   ylim: {[begin end],'maxmin','absmax','all'} {default: absmax}
%   zlim: {[begin end],'maxmin','absmax'} {default: absmax}
%   xparam
%   yparam
%   zparam
%
% Function call:
% To plot TFR: ts_xxx(timefreq_data,...)
% To plot ERP: ts_xxx(avg_data,...)
err   = 0;
parms = mmil_args2parms(varargin,...
						{'baseline','no',{'absolute','relative','relchange','zscore','no'},...
						 'baseline_start',-inf,[],...
             'baseline_end',0,[],...
             'xlim','maxmin',[],...
             'ylim',[],[],...
             'zlim','maxmin',[],...
             'xparam',[],[],...
             'yparam',[],[],...
             'zparam',[],[],...
             'layout',[],[],...
             'colorbar','no',{'yes','no'},...
             'showlabels','yes',{'yes','no'},...
             'zerolines','no',{'yes','no'},...
             'box','no',{'yes','no'},...
             'axes','yes',{'yes','no','all'},...
             'axis','tight',[],...
             'masknans','yes',{'yes','no'},...
             'comment',date,[],...
             'title','',[],...
             'figname',date,[],...
             'channels',[],[],...
             'chantype','all',[],...
             'conditions',[],[],...
             'events',[],[],...
             'mask',[],[],...
             'stat_data',[],[],...
             'grad',[],[],...
             'elec',[],[],...
             'fontsize',8,[],...
             'axisfontsize',4,[],...
             'newfig',0,{0,1},...
             'background','w',[],...
             'linewidth',.5,[],...
             'linestyle','-',[],...
             'marker','none',[],...
             'graphcolor','brgkywrgbkywrgbkywrgbkywbrgkywrgbkywrgbkywrgbkyw',[],...
             'statcolor','b',[],...
             'statmarker','*',[],...
             'transparency',0,[],... 0.1 % BQR 9/24/12
             'clip','yes',{'yes','no'},...
             'colormap','jet',[],...
             'cond_labels',[],[],...
             'footnotes','',[],...
             'autoscale',0,{0,1},...
             'alpha',.05,[],...
             'maxtickmarks',6,[],...
             'markersize',20,[],...
             'ticklength',.1,[],...
             'trials_flag',0,{1,0},...
             'vline',[],[],...
             'hline',[],[],...
             'vlinewidth',.5,[],...
             'hlinewidth',.5,[],...
             'vlinecolor','k',[],...
             'hlinecolor','k',[],...
             'zerolinetimes',0,[],...
             'zerolinewidth',.25,[],...
             'zerolinecolor','k',[],...
						},false);

          
% if ~isempty(parms.mask) && ndims(parms.mask)==2
%   parms.mask = parms.mask(~[data.sensor_info.badchan],:);
% end
parms.graphcolor = strrep(parms.graphcolor,parms.background,'');

data = ts_data_selection(data,varargin{:});
[datatype datafield dataparam] = ts_object_info(data);

mask_flag = ~isempty(parms.mask) || ~isempty(parms.stat_data);
temporary_files = {};
tfr_flag = 1;
npages   = 1;
% set default xparam, yparam, zparam, & a few flags
if (strcmp(datafield,'timefreq') || strcmp(datafield,'averages')) && any(strcmp(dataparam,'power'))
  % power
  if isempty(parms.xparam), parms.xparam = 'time';        end
  if isempty(parms.yparam), parms.yparam = 'frequencies'; end
  if isempty(parms.zparam), parms.zparam = 'power';       end
elseif strcmp(datafield,'timefreq') && any(strcmp(dataparam,'coh'))
  % coherence
  if isempty(parms.xparam), parms.xparam = 'time';        end
  if isempty(parms.yparam), parms.yparam = 'frequencies'; end
  if isempty(parms.zparam), parms.zparam = 'coh';         end
elseif strcmp(datafield,'timefreq') && any(strcmp(dataparam,'plv'))
  % phase-locking values
  if isempty(parms.xparam), parms.xparam = 'time';        end
  if isempty(parms.yparam), parms.yparam = 'frequencies'; end
  if isempty(parms.zparam), parms.zparam = 'plv';         end  
elseif strcmp(datafield,'timefreq') && any(strcmp(dataparam,'cmplx'))
  % convert complex spectra to power
  for k = 1:length(data.(datafield))
    data.(datafield)(k).power     = 2*abs(data.(datafield)(k).cmplx).^2;
  end
  data.(datafield)                = rmfield(data.(datafield),'cmplx');
  [datatype datafield dataparam]  = ts_object_info(data);
  if isempty(parms.xparam), parms.xparam = 'time';        end
  if isempty(parms.yparam), parms.yparam = 'frequencies'; end
  if isempty(parms.zparam), parms.zparam = 'power';       end  
elseif (strcmp(datafield,'cont') || strcmp(datafield,'averages')) && any(strcmp(dataparam,'data'))
  % average waveforms
  tfr_flag = 0;
  if isempty(parms.xparam),  parms.xparam = 'time';        end
  if isempty(parms.zparam),  parms.zparam = 'data';        end
  if ~isempty(parms.yparam), parms.zparam = parms.yparam;  end
  parms.yparam = [];
elseif strcmp(datafield,'epochs') && any(strcmp(dataparam,'data'))
  % convert epochs to average waveforms
  tfr_flag = 0;
  if ~parms.trials_flag,     data = ts_trials2avg(data);   end
  [datatype datafield dataparam]  = ts_object_info(data);
  if isempty(parms.xparam),  parms.xparam = 'time';        end
  if isempty(parms.zparam),  parms.zparam = 'data';        end  
  if ~isempty(parms.yparam), parms.zparam = parms.yparam;  end
  parms.yparam = [];  
else
  fprintf('Failed to determine the data type.\n');
  err = 1;
  return;
end

% Apply baseline correction:
if ~strcmp(parms.baseline, 'no')
  data = ts_blc(data,'baseline',parms.baseline,'baseline_start',parms.baseline_start,'baseline_end',parms.baseline_end);
%   data = ts_zscore(data,'baselinetype',parms.baseline);
end

ieeg_flag = (strcmpi(parms.chantype,'ieeg') || (ismember('eeg',unique({data.sensor_info.typestring})) ...
            && isempty(intersect({'grad1','grad2','mag','meg'},{data.sensor_info.typestring}))));

% Set up channel layout
if ~isempty(parms.layout) && isstruct(parms.layout)
  lay = parms.layout;
else
  if ~isempty(parms.layout) && ischar(parms.layout)   
    if any(findstr(parms.layout,'.lay'))              % use fieldtrip layout
      parms.layout = parms.layout;
    elseif any(findstr(parms.layout,'.lout'))         % use neuromag layout
      % convert neuromag layout to fieldtrip layout  
      [jnk1 layfile jnk2] = fileparts(parms.layout);
      layfile = [layfile '.lay'];
      unix(sprintf('cp %s %s',parms.layout,layfile));
      parms.layout = layfile;
      if ~exist(layfile,'file')
        temporary_files = {temporary_files{:} layfile};
      end
      clear layfile
    elseif any(findstr(parms.layout,'.asc'))          % use neuroscan layout
      % convert neuroscan montage to fieldtrip layout
    end
  elseif isempty(parms.layout) && ieeg_flag           % create ieeg layout
    parms.layout = write_layout_ieeg(pwd,data.sensor_info,10,10);
    for k = 1:length(parms.layout)
      if ~exist(parms.layout{k})
        temporary_files = {temporary_files{:} parms.layout{k}};
      end
    end
  elseif isempty(parms.layout) && ~isempty(intersect({'grad1','grad2','mag','meg'},{data.sensor_info.typestring}))
    try 
      parms.grad = getfield(ts_data2fieldtrip(data),'grad'); % create 'grad' structure for MEG
    catch
      parms.layout = 'ordered';                       % subplot all in order
    end  
  end
  if ~isempty(parms.layout) && iscell(parms.layout)
    cfg.layout = parms.layout{1};
    npages = length(parms.layout);
  elseif ~isempty(parms.layout)
    cfg.layout = parms.layout;
  elseif ~isempty(parms.grad)
    cfg.grad = parms.grad;
  elseif ~isempty(parms.elec)
    cfg.elec = parms.elec;
  else
    cfg.layout = 'ordered';
  end
  for p = 1:npages
    if p > 1
      cfg.layout = parms.layout{p};
    end
    if strcmp(cfg.layout,'ordered')
      tmpdat.label = {data.sensor_info.label};
      lay(p) = prepare_layout(cfg,tmpdat);
    else
      lay(p) = prepare_layout(cfg);
    end
  end
end

%% Set up z-limits
if tfr_flag
%   Get physical y-axis range:
  if isempty(parms.ylim) || ~isnumeric(parms.ylim)
    ymin = min([data.(datafield).(parms.yparam)]);
    ymax = max([data.(datafield).(parms.yparam)]);
  else
    ymin = parms.ylim(1);
    ymax = parms.ylim(2);
  end
end

% events
if isempty(parms.events) && isempty(parms.conditions)
  parms.conditions = 1:length(data.(datafield));
elseif isempty(parms.conditions)
  [jnk1 parms.conditions jnk2] = intersect([data.(datafield).event_code],parms.events);
end
if isempty(parms.cond_labels)
  for c = 1:length(parms.conditions)
    parms.cond_labels{c} = sprintf('Event %g',data.(datafield)(parms.conditions(c)).event_code);
  end
elseif length(parms.cond_labels) > length(parms.conditions)
  parms.cond_labels = parms.cond_labels(parms.conditions);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN LOOP OVER PAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:length(lay)
  if parms.newfig || length(lay) > 1
    screensize = get(0,'ScreenSize');
%     figure('Name',parms.figname,'NumberTitle','off','Color',parms.background,...
%            'Position',[1 screensize(4)/2 screensize(3)/2 screensize(4)/2]);
%     figwidth = min(screensize(3),screensize(4));
%     figure('Name',parms.figname,'NumberTitle','off','Color',parms.backgro
%     und,'Position',[1 1 figwidth figwidth]);
    figure('Name',parms.figname,'NumberTitle','on','Color',parms.background,'Position',[1 1 .9*screensize(3) .9*screensize(4)]);
  end

  % Select the channels in the data that match with the layout:
  [seldat, sellay] = match_str({data.sensor_info.label}, lay(p).label);
  if isempty(seldat)
    fprintf('%s: labels in data and labels in layout do not match\n',mfilename); 
    err = 1;
    return;
  end
  chanX       = lay(p).pos(sellay, 1); cx = max(chanX) - min(chanX);
  chanY       = lay(p).pos(sellay, 2); cy = max(chanY) - min(chanY);
  chanWidth   = lay(p).width(sellay);  
  chanHeight  = lay(p).height(sellay); ch = mean(chanHeight) / mean(chanWidth);
  chanLabels  = lay(p).label(sellay);

  if ~ieeg_flag
    if (any(chanX+chanWidth > 1) || any(chanY+chanHeight > 1))
      % translate to positive values
      if min(chanX) < 0, chanX = chanX - min(chanX); end
      if min(chanY) < 0, chanY = chanY - min(chanY); end
  %     if min(chanX) <= 0, chanX = chanX - min(chanX) + 0.05; end
  %     if min(chanY) <= 0, chanY = chanY - min(chanY) + 0.05; end

      % scale 0 to 1 and maintain aspect ratio
      if cy ~= 0
        chanY = chanY / cy;
      end
      if cx ~= 0
        chanX = chanX / cx;
        chanWidth  = chanWidth / cx;
      end
      chanHeight = chanWidth * ch;

      % if different chantypes overlap
      if ~strcmp(data.sensor_info(seldat(1)).typestring,data.sensor_info(seldat(2)).typestring)
        if abs(chanY(2)-chanY(1)) < chanHeight(1)
          chanHeight = .99*abs(chanY(2)-chanY(1))*ones(size(chanHeight));
          chanWidth  = chanHeight / ch;
        elseif abs(chanX(2)-chanX(1)) < chanWidth(1)
          chanWidth  = .99*abs(chanX(2)-chanX(1))*ones(size(chanWidth));
          chanHeight = chanWidth * ch;
        end
      end
    end
    % apply proportional scaling factor
    cf = .8;
    chanX = cf*chanX + (.9-cf)/2; % .9 instead of 1 b/c of legend
    chanY = cf*chanY + (.9-cf)/2;
    chanWidth  = chanWidth * cf;
    chanHeight = chanWidth * ch;
  else
%     % TEMPORARY CODE: remove after tests
%     chanX = (chanX - min(chanX))*cx + .05;
%     chanY = (chanY - min(chanY))*cy + .05;
%     chanWidth  = chanWidth * cx;
%     chanHeight = chanHeight * cy;
%     % end remove
  end
  
  % Get physical z-axis range
  if strcmp(parms.zlim,'maxmin') || strcmp(parms.zlim,'absmax')
    % get xlims over all conds to determine zlims
    if strcmp(parms.xlim,'maxmin')
      xmintmp = max(cellfun(@min,{data.(datafield).(parms.xparam)}));
      xmaxtmp = min(cellfun(@max,{data.(datafield).(parms.xparam)}));
    else
      xmintmp = parms.xlim(1);
      xmaxtmp = parms.xlim(2);
    end
    zmin=[]; zmax=[];
    for i = 1:length(parms.conditions)
      tix  = nearest(data.(datafield)(parms.conditions(i)).(parms.xparam),xmintmp):nearest(data.(datafield)(parms.conditions(i)).(parms.xparam),xmaxtmp);
      ztmp = [data.(datafield)(parms.conditions(i)).(parms.zparam)(seldat,tix)];
      zmin = min([zmin min([ztmp(:)])]);
      zmax = max([zmax max([ztmp(:)])]);
      clear ztmp
    end
    if strcmp(parms.zlim,'absmax')
      zmin = -max(abs(zmin),abs(zmax));
      zmax = max(abs(zmin),abs(zmax));
    end
  else
    zmin = parms.zlim(1);
    zmax = parms.zlim(2);
  end 
  if ~tfr_flag
    ymin = zmin;
    ymax = zmax;
  end    
  xscales = [inf*ones(length(seldat),1) -inf*ones(length(seldat),1)];
  yscales = [inf*ones(length(seldat),1) -inf*ones(length(seldat),1)];
%   xscales = zeros(length(seldat),2);
%   yscales = zeros(length(seldat),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BEGIN LOOP OVER CHANNELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  numbadchans    = 0;
  showscale_flag = 0;
  for k = 1:length(seldat)
%       if ~any(any(any(datavector(k,:))))
%         continue;
%       end
    plotpos = [chanX(k) chanY(k) chanWidth(k) chanHeight(k)];
    if all(isfinite(plotpos))
      subplot('Position',[chanX(k) chanY(k) chanWidth(k) chanHeight(k)])
      hold on
    else
      numbadchans = numbadchans + 1;
      if (k == (find((chanX==min(chanX)) & (chanY==min(chanY(find(chanX==min(chanX))))))))
        showscale_flag = 1;
      end
      continue;
    end
%%       
% 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % BEGIN LOOP OVER CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  for c = 1:length(parms.conditions)
    % Get physical x-axis range:
    if strcmp(parms.xlim,'maxmin')
      xmin = min(data.(datafield)(c).(parms.xparam));
      xmax = max(data.(datafield)(c).(parms.xparam));
    else
      xmin = parms.xlim(1);
      xmax = parms.xlim(2);
    end
    % Find corresponding x-axis bins:
    xidc = nearest(data.(datafield)(c).(parms.xparam),xmin):nearest(data.(datafield)(c).(parms.xparam),xmax);
    % Align physical x-axis range to the array bins:
    xmin = data.(datafield)(c).(parms.xparam)(xidc(1));
    xmax = data.(datafield)(c).(parms.xparam)(xidc(end));

    if tfr_flag
      yidc = nearest(data.(datafield)(c).(parms.yparam),ymin):nearest(data.(datafield)(c).(parms.yparam),ymax);
      % Align physical y-axis range to the array bins:
      ymin = data.(datafield)(c).(parms.yparam)(yidc(1));
      ymax = data.(datafield)(c).(parms.yparam)(yidc(end));
      datavector = data.(datafield)(c).(parms.zparam)(seldat,xidc,yidc);
      if ~isempty(parms.mask) && isequal(size(parms.mask),size(data.(datafield)(c).(parms.zparam)))
%         maskvector = parms.mask(seldat,yidc,xidc);
        maskvector = parms.mask(seldat,xidc,yidc);
      end
    else
      datavector = data.(datafield)(c).(parms.zparam)(seldat,xidc,:);
      if ~isempty(parms.mask) && isequal(size(parms.mask),size(data.(datafield)(c).(parms.zparam)))
        maskvector = parms.mask(seldat,xidc,:);
      end  
    end
    %% (copied from channel for-loop: remove all zeros or flat lines
      if ~any(any(any(datavector(k,:,:)))) || (~tfr_flag && all(all(datavector(k,:,1)==datavector(k,1,1))))
        axis off
        if (k == (find((chanX==min(chanX)) & (chanY==min(chanY(find(chanX==min(chanX))))))))
          showscale_flag = 1;
        end
        continue;
      end    
    %%
    
%%
      x = data.(datafield)(c).(parms.xparam)(xidc);
      
      if tfr_flag
        cdata = squeeze(datavector(k,:,:))';
        y     = data.(datafield)(c).(parms.yparam)(yidc);
        imagesc(x,y,cdata,[zmin zmax]);
      else
        y = squeeze(datavector(k,:,:));
        if strcmp(parms.clip,'yes')
          y(y > ymax) = ymax;
          y(y < ymin) = ymin;
        end
        if parms.trials_flag
          plot(x,y,'linewidth',parms.linewidth,'Clipping','off');
          plot(x,mean(y,2),'k-','linewidth',2,'Clipping','off')
        else
          plot(x,y,parms.graphcolor(c),'linestyle',parms.linestyle,'linewidth',parms.linewidth,'Clipping','off')
        end
      end
      
      % apply mask
      % single condition stats
      % comparison between conditions
      if mask_flag && ~tfr_flag && size(maskvector,1) >= k
        mask = squeeze(maskvector(k,:)) > 0;
        n = length(mask);b
        if mask(1) > 0, m0 = 1; else m0 = 0; end
        if mask(n) > 0, mf = 1; else mf = 0; end
        mask  = [m0 mask(2:end-1).*(((mask(1:n-2)==0) + (mask(3:n)==0)) > 0) mf];
        
        edges = find(mask);
        for kk = 1:length(edges)/2
          ix1 = x(edges(2*kk-1));
          ix2 = x(edges(2*kk));
          iy1 = ymin;
          iy2 = ymax; 
          % edited by BQR 9/24/12
          fill([ix1 ix1 ix2 ix2],[iy1 iy2 iy2 iy1],parms.statcolor,'EdgeColor',parms.statcolor,...
            'FaceAlpha',parms.transparency,'EdgeAlpha',parms.transparency);
%           fill([ix1 ix1 ix2 ix2],[iy1 iy2 iy2 iy1],parms.statcolor,'EdgeColor','k',...
%             'FaceAlpha',parms.transparency,'EdgeAlpha',0);
        end                      
      end
      

      xscales(k,:) = [min(xscales(k,1),min(x(:))) max(xscales(k,2),max(x(:)))];
      yscales(k,:) = [min(yscales(k,1),min(y(:))) max(yscales(k,2),max(y(:)))];
      
      % apply channel-specific plot options
      if c == length(parms.conditions)
        % Channel labels
        if strcmp(parms.showlabels,'yes')
          title(chanLabels{k},'FontSize',parms.fontsize,'VerticalAlignment','baseline','Clipping','off')
        end
        % Box
        if strcmp(parms.box,'no')
          box off
        else
          box on
          set(gca,'linewidth',parms.linewidth);
        end
        % Axis limits
        if parms.autoscale
          xlim([xscales(k,1) xscales(k,2)])
          ylim([yscales(k,1) yscales(k,2)])
          if tfr_flag
            caxis([min(cdata(:)) max(cdata(:))]);
          end
        else
          xlim([xmin xmax])
          ylim([ymin ymax])
          if tfr_flag
            caxis([zmin zmax]);
          end
        end
        if ~isempty(parms.vline), for l = 1:length(parms.vline), vline(parms.vline(l),'LineWidth',parms.vlinewidth,'Color',parms.vlinecolor); end; end
        if ~isempty(parms.hline)
          if tfr_flag
            parms.hline = intersect(parms.hline,y);
            for l = 1:length(parms.hline)
              yind = find(parms.hline(l) == y);
              hline(diff(ylim)*(yind/length(y)),'LineStyle','-','LineWidth',parms.hlinewidth,'Color',parms.hlinecolor);
            end
          else
            for l = 1:length(parms.hline), hline(parms.hline(l),'LineStyle','-','LineWidth',parms.hlinewidth,'Color',parms.hlinecolor); end
          end
        end
        % Axes
        if ~strcmp(parms.axes,'no') || strcmp(parms.zerolines,'yes')
          axis on
          % check the number of tick marks
          yticks = get(gca,'ytick');
          xticks = get(gca,'xtick');    
          if strcmp(parms.axes,'no')
            hline(0,'k');
            vline(0,'k');
          else
            if tfr_flag
              ntick = min(parms.maxtickmarks,length(yticks));
              ylims = get(gca,'ylim');
              ylocs = round(linspace(ylims(1),ylims(2),ntick));
              yinds = round(linspace(1,length(y),ntick));
              yvals = y(yinds);
              set(gca,'ytick',ylocs);
              set(gca,'yticklabel',yvals);
              scatter(xmin*ones(1,length(ylocs)),ylocs,'+k','linewidth',.25,'Clipping','on','sizedata',parms.markersize);          
            else
              if length(yticks) > parms.maxtickmarks, set(gca,'ytick',yticks(1:2:end)); end
              yticks = get(gca,'ytick');
              scatter(zeros(1,length(yticks)),yticks,'+k','linewidth',.25,'Clipping','on','sizedata',parms.markersize);
            end
            if length(xticks) > parms.maxtickmarks
              set(gca,'xtick',xticks(1:2:end)); 
            end 
          end
          set(gca,'FontSize',parms.axisfontsize);
%           plot(zeros(2),ylim,parms.linestyle,'color','black','linewidth',.25,'Clipping','on');
          for l = 1:length(parms.zerolinetimes)
            vline(parms.zerolinetimes(l),'LineWidth',parms.zerolinewidth,'Color',parms.zerolinecolor,'Clipping','on');
          end
          if strcmp(parms.axes,'yes') || strcmp(parms.axes,'all')
            xticks = get(gca,'xtick');
            if min(ylim) <= 0
%               plot(xlim,zeros(2),parms.linestyle,'color','black','linewidth',.25,'Clipping','on');
              plot(xlim,zeros(2),'-','color','black','linewidth',.25,'Clipping','on');
              xmark = 0;
            else
              xmark = min(ylim);
            end
            try scatter(xticks,xmark*ones(1,length(xticks)),'+k','linewidth',.25,'Clipping','on','sizedata',parms.markersize); end
          end
          if ~strcmp(parms.axes,'all') && showscale_flag==0 && any((k ~= (find((chanX==min(chanX)) & (chanY==min(chanY(find(chanX==min(chanX)))))))))
            if strcmpi(parms.box,'yes')
              set(gca,'ytick',[]);
              set(gca,'xtick',[]);
            else
              axis off
            end
          elseif showscale_flag == 1
            showscale_flag = 0;
          end
%           set(gca,'ticklength',[parms.ticklength 0]);
        else
          axis off
        end
      end
    end % end loop over conditions
    
   % edited by BQR 9/24/12
   % move masks to the back /lines to the front
   chldrn = get(gca,'Children');
   lnchldrn = findobj(gca,'Type','line');
   [~,iH] = sort(ismember(chldrn,lnchldrn),1,'descend');
   set(gca,'Children',chldrn(iH));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
  end % end loop over channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  if ~strcmpi(parms.comment,'none')
    % add condition to legend
    for c = 1:length(parms.conditions)
      annotation('textbox','Position',[.9 .85-(c-1)*.05 .1 .1],'VerticalAlignment','top',...
                 'HorizontalAlignment','left','FontSize',parms.fontsize,'Color',parms.graphcolor(c),'FitHeightToText','on','LineStyle','none',...
                 'String',sprintf('%s\nNum Trials: %g',parms.cond_labels{c},data.(datafield)(c).num_trials));
    end
    % add alpha value if supplied
    if mask_flag
      annotation('textbox','Position',[.9 .85-c*.05 .1 .1],'VerticalAlignment','top',...
                 'HorizontalAlignment','left','FontSize',parms.fontsize,'Color','k','FitHeightToText','on','LineStyle','none',...
                 'String',sprintf('\nSig: p < %g',parms.alpha));
      c = c + 1;
    end    
  end
  if numbadchans < length(seldat)
  %   % write comment:
    if ~isempty(parms.comment) && ~strcmpi(parms.comment,'none')
      comment = parms.comment;
      comment = sprintf('%0s\nxlim=[%.2g %.2g]', comment, xmin, xmax);
      if parms.autoscale
        comment = sprintf('%0s\nylim=autoscale', comment);
      elseif ymax >= 1E4 || round(ymax)~=ymax || round(ymin)~=ymin
        comment = sprintf('%0s\nylim=[%.2g %.2g]', comment, ymin, ymax);
      else
        comment = sprintf('%0s\nylim=[%g %g]', comment, ymin, ymax);
      end
      if tfr_flag
        comment = sprintf('%0s\nzlim=[%.2g %.2g]', comment, zmin, zmax);
      end
      comment = sprintf('%0s\n\n%s',comment,date);
      k = cellstrmatch('COMNT',lay(p).label);
      if ieeg_flag && ~isempty(k)
        annotation('textbox','Position',[.9 .85-c*.05 .1 .1],'VerticalAlignment','top',...
                   'HorizontalAlignment','left','FontSize',parms.fontsize,'Color','k','FitHeightToText','on','LineStyle','none',...
                   'String',comment);      
        c = c + 1;
  %       text(lay(p).pos(k,1), lay(p).pos(k,2), sprintf(comment), 'Fontsize', parms.fontsize);
      else
        annotation('textbox','Position',[.9 .85-c*.05 .1 .1],'VerticalAlignment','top',...
                   'HorizontalAlignment','left','FontSize',parms.fontsize,'Color','k','FitHeightToText','on','LineStyle','none',...
                   'String',comment);      
        c = c + 1;      
  %       comntx  = min(lay(p).pos(:,1));
  %       comnty  = min(lay(p).pos(:,2)); comnty = comnty + .15*(max(lay(p).pos(:,2))-comnty);
  %       text(comntx, comnty, sprintf(comment), 'fontsize', parms.fontsize);   
      end
    end

    % Title
    parms.title = strrep(parms.title,'_','\_');
    annotation('textbox','Position',[0 .9 1 .1],'VerticalAlignment','middle','HorizontalAlignment','center',...
               'Color','k','FontSize',10,'FontWeight','bold','FitHeightToText','on','LineStyle','none','String',parms.title);

    % Add footnotes (other annotations => for instance, filename)
    annotation('textbox','Position',[.9 .85-c*.05 1 .1],'VerticalAlignment','top','HorizontalAlignment','center',...
             'Color','k','FontSize',6,'FitHeightToText','on','LineStyle','none','String',parms.footnotes);

    % Colormap
    colormap(parms.colormap);

    % Colorbar
  end
  hold off

  % save figures
  % set axes
  
end % end loop over pages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove temporary files
for k = 1:length(temporary_files)
  unix(sprintf('rm %s',temporary_files{k}));
end

% SUBFUNCTIONS
function l = cellstrmatch(str,strlist)
l = [];
for k=1:length(strlist)
  if strcmp(char(str),char(strlist(k)))
    l = [l k];
  end
end
