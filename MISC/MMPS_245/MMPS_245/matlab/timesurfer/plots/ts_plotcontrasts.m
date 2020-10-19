function ts_plotcontrasts(data,stats,varargin)

parms = mmil_args2parms(varargin,...
						{ ...
             'channels',[],[],...
             'chanlabels',[],[],...
             'contrasts',[],[],...
						},false);

if ~iscell(data) , data  = {data} ; end
if ~iscell(stats), stats = {stats}; end

if numel(data)~=numel(stats) && ~isempty(stats{1})
  error('Stats must be supplied for all data sets or none.');
end
nsets = numel(data);

% check contrasts
if isempty(parms.contrasts)
  error('Must specify contrasts to plot.');
elseif ~iscell(parms.contrasts)
  parms.contrasts = {parms.contrasts};
end
events = unique([parms.contrasts{:}]);
ncontrasts = numel(parms.contrasts);
% check that all contrasts contain two event codes
% ...

% Check input data
% check that events in all contrasts are present in all data sets
% check that stats have been provided for all contrasts & data sets
% make sure data limits match between different data sets & stats

for k = 1:nsets
  data{k} = ts_data_selection(data{k},'events',events,'channels',parms.channels,'chanlabels',parms.chanlabels,'verbose',0);
  data{k} = ts_trials2avg(data{k},'verbose',0); % average over trials if given epochs
  [type,fld,parm] = ts_object_info(data{k});
  datatype{k}  = type;
  datafield{k} = fld;
  dataparam{k} = parm{1};
end
clear type fld parm
% make sure channels match between all data sets
% ...

nchan  = data{1}.num_sensors;
colors = 'brgkywrgbkywrgbkywrgbkywbrgkywrgbkywrgbkywrgbkyw';
alpha  = .05;
nrows  = ceil(ncontrasts/2);
ncols  = 2*nsets;

% Loop over channels (plot all contrasts for each in separate figure)
for ch = 1:nchan
  % select data
  for k = 1:nsets
    seldata{k} = ts_data_selection(data{k},'channels',ch,'verbose',0);
    for j = 1:ncontrasts
      selstats{k}(j) = ts_data_selection(stats{k}(j),'channels',ch,'verbose',0);
    end
  end
  % create figure
  figure('color','w');
  subhandles = btq_panels(nrows,ncols,'top');
  % plot all contrasts
  plot_contrasts; % nested function
  % save figures
  save_figures; % nested function
   
end

  % NESTED FUNCTIONS
  function plot_contrasts
    cnt   = 1;
    for j = 1:ncontrasts
      for k = 1:nsets
        % select data for this contrast and data set
        tmpdata = ts_data_selection(seldata{k},'events',parms.contrasts{j},'verbose',0); % data for one chan & two conds
        tmpstat = selstats{k}(j); % stats for one chan & two conds        
        tmpcond = [tmpdata.(datafield{k}).event_code];
        % plot waveforms (overlay conditions in this contrast)
        yi = ceil(cnt/(2*nsets));
        xi = mod(cnt,2*nsets); if xi==0, xi=2*nsets; end
        d = (yi-1)*ncols+xi;
        subplot(subhandles(d))      
%         subplot(nrows,ncols,cnt)
        for c = 1:2 % loop over events in this contrast
          t = tmpdata.(datafield{k})(c).time;
          y = tmpdata.(datafield{k})(c).data;
          err = tmpdata.(datafield{k})(c).stdev / sqrt(tmpdata.(datafield{c})(c).num_trials);
          btq_confplot_3andC(t,y,err,colors(k)); hold on
        end        
        % add between-condition statistics
        tix   = 1:length(tmpstat.stats(1).time);
        edges = shade_stats(squeeze(tmpstat.stats(1).prob(1,tix)<alpha),tmpstat.stats(1).time(tix));
        
        % add within-condition statistics
%         line(parms.blcwindow,[ylims(1) ylims(2)]*.95,'Color','g','LineWidth',2);
        for c = 1:2
          ind = find([tmpstat.stats.event_code]==tmpdata.(datafield{k})(c).event_code);
          tix = 1:length(tmpstat.stats(ind).time);
          edges = bl_line_stats(squeeze(tmpstat.stats(ind).prob(1,tix)<alpha),tmpstat.stats(ind).time(tix),colors(c),.95,3-(c-1));
        end % end loop over conditions in this contrast
        title(sprintf('contrast %g, set %g',j,k));
        legend(sprintf('event %g',tmpcond(1)),sprintf('event %g',tmpcond(2)));
        cnt = cnt + 1;
      end % end loop over data sets
    end % end loop over contrasts
  end

  function save_figures
    disp('This is when figures will be saved.');
    % ...
%     export_fig(sprintf('%s%s/html/images/%s%s.png',STUDYDIR,SUBJ,prefname,chan));
%     print('-depsc2','-tiff','-painters','-r300',sprintf('%s%s/html/images/%s%s.eps',STUDYDIR,SUBJ,prefname,chan));
  end

end
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edges=shade_stats(sigmask,sigtime)
edges=[];
if ~isempty(sigmask),
  ylims=get(gca,'YLim');
  xlims=get(gca,'XLim');
  onoffs=[sigmask 0]-[0 sigmask];
  ons=find(onoffs(1:end-1)==1);
  offs=find(onoffs(2:end)==-1);
  if sum(ons) > 0,
    for n = 1:length(ons),        
      fill([sigtime(ons(n)) sigtime(ons(n)) sigtime(offs(n)) sigtime(offs(n))], ...
          [ylims fliplr(ylims)],'m','EdgeColor','none','FaceAlpha',0.4);
     end
  end
  edges=[ons; offs];
end
end

function edges=bl_line_stats(sigmask,sigtime,color,offset,lw)
edges=[];
if ~isempty(sigmask),
  ylims=get(gca,'YLim');
  xlims=get(gca,'XLim');
  onoffs=[sigmask 0]-[0 sigmask];
  ons=find(onoffs(1:end-1)==1);
  offs=find(onoffs(2:end)==-1);
  if sum(ons) > 0,
    for n = 1:length(ons),
      line([sigtime(ons(n)) sigtime(offs(n))],[ylims(1) ylims(1)]*offset,'Color',color,'LineWidth',lw);
    end
  end
  edges=[ons; offs];
end
end

% function draw_axes()
% xticks = get(gca,'xtick');
% xticks=round(xticks.*1000)./1000;
% yticks = get(gca,'ytick');
% yticks_zero=find(abs(yticks)==min(abs(yticks)));
% yticks(yticks_zero)=[];
% ylims=get(gca,'YLim');
% xlims=get(gca,'XLim');
% line([0 0],repmat(ylims,2,1)','Color','k');
% line(repmat(xlims,2,1)',[0 0],'Color','k');
% scatter(zeros(1,length(yticks)),yticks,'+k','linewidth',.25,'Clipping','on','sizedata',20);
% scatter(xticks(2:end),zeros(1,length(xticks)-1),'+k','linewidth',.25,'Clipping','on','sizedata',25);
% text(xticks,zeros(1,length(xticks)),cellstr(num2str(xticks','%g')), ...
%     'FontSize',4,'HorizontalAlignment','Left','VerticalAlignment','Bottom');
% text(zeros(1,length(yticks)),yticks,cellstr(num2str(yticks','%g')), ...
%     'FontSize',4,'HorizontalAlignment','Right','VerticalAlignment','Bottom');
% end
% 
% function write_edges(fid,edges)
%     for nedge = 1:size(edges,2),
%         fprintf(fid,'%02i %g %g\n',nedge,edges(1,nedge),edges(2,nedge));
%     end
% end
% 
% function draw_legend(leg_str)
% leghand=legend(strrep(leg_str,'_','\_'),'FontSize',8,'Location','NorthWest','Color','none');
% legend(gca,'boxoff');
% legpos=get(leghand,'Position');
% olegpos=get(leghand,'OuterPosition');
% set(leghand,'Position',[olegpos(1) legpos(2:4)])
% legkids=findobj(get(leghand,'Children'));
% for lk=1:length(legkids),
%     switch lower(get(legkids(lk),'type'))
%         case 'line'
%             Xdata=get(legkids(lk),'XData');
%             if numel(get(legkids(lk),'XData'))==2,
%                 Xscale=Xdata(2)-Xdata(1);
%                 Xdata(2)=Xdata(1)+Xscale/10;
%                 set(legkids(lk),'XData',[0 Xdata(2)]);
%             end
%         case 'text'
%             tpos=get(legkids(lk),'Position');
%             lpos=get(legkids(lk-1),'XData');
%             set(legkids(lk),'Position',[lpos(2)*1.1 tpos(2) tpos(3)]);
%             set(legkids(lk),'HorizontalAlignment','left');
%     end
% end
% end
% 
% function zp_chnames=zeropad_chan_names(data_chan_names)
% zp_chnames=cell(size(data_chan_names));
% for cn=1:length(data_chan_names),
%   zp_chnames{cn}=sprintf('%s%02g',upper(data_chan_names{cn}(regexp(data_chan_names{cn},'[aA-zZ]'))),str2num(data_chan_names{cn}(regexp(data_chan_names{cn},'[0-9]'))));
% end
% end
% 
% function ylimits=calc_ylimits(avg_data,chans)
%     ch_inds=zeros(length(chans),1);
%     zp_chnames=zeropad_chan_names({avg_data.sensor_info(:).label});
%     for ch=1:length(chans),
%         ch_inds(ch)=strmatch(chans{ch},zp_chnames,'exact');
%     end
%     datchan_order=zeropad_chan_names({avg_data.sensor_info(:).label});
%     datm=horzcat(avg_data.averages(:).data);
%     daterr=[];
%     for ev=1:length(avg_data.averages),
%         tmperr=(avg_data.averages(ev).stdev./sqrt(avg_data.averages(ev).num_trials));
%         daterr=horzcat(daterr,tmperr);
%     end
%     dat=datm+abs(daterr);
%     datn=datm-abs(daterr);
%     ylimits=[min(datn,[],2) max(dat,[],2)];
%     sort_inds=zeros(length(chans),1);
%     for c=1:length(chans),
%         sort_inds(c)=strmatch(chans{c},datchan_order);
%     end
%     ylimits=ylimits(sort_inds,:);
% end  
%   
  