function ts_MEGtopoplot_callFT(param,data)

% Required inputs:
%   timefreq - timesurfer timefreq_data structure
%   param
%
% Optional inputs:
% 	param
% 	  zlim		[zmin zmax] (default: 'maxmin')
% 		nr			# subplot rows (when plotting a sequence of topoplots)
% 		nc 			# subplot columns (when plotting a sequence of topoplots)
%    style           = topoplot style (default = 'both')
%                       'straight' colormap only
%                       'contour' contour lines only
%                       'both' (default) both colormap and contour lines
%                       'fill' constant color between lines
%                       'blank' just head and electrodes
%     title         figure name (default: '');
%
% Note: param.nr and param.nc must be specified to plot a sequence of topoplots

% Created: 10-Sep-2008 by Jason Sherfey
% Revision: 05-Oct-2008 by Jason Sherfey
% Modifed to take timesurfer structure as input
% Revision: 11-Nov-2008 by Jason Sherfey
% Renamed; used to be called ts_MEGtopoplot

if isfield(data,'averages')
	datafield = 'averages';
	zparam = 'avg';
	freqs = 1;
    if ~isfield(param,'toilim'),param.toilim=[data.averages(1).time(1) data.averages(1).time(end)]; end;
elseif isfield(data,'timefreq')
	datafield = 'timefreq';
	zparam = 'powspctrm';
	freqs = data.timefreq(1).frequencies;
    if ~isfield(param,'toilim'),param.toilim=[data.timefreq(1).time(1) data.timefreq(1).time(end)]; end;
else
	error('Data must be avg_data or timefreq_data.');
end

% Defaults
if ~isfield(param,'showhighlight'), param.showhighlight = 0; 		 end
if ~isfield(param,'showmask'), 		param.showmask = 0; 
elseif ~isfield(param,'freqmask'), 	param.freqmask = 'matchany'; 	 end
if isfield(param,'masktype'),       param.freqmask = param.masktype; end
if ~isfield(param,'nr'), 			param.nr = 5;                    end
if ~isfield(param,'nc'), 			param.nc = 5;                    end
if ~isfield(param,'zlim'), 			param.zlim = 'maxmin';           end
if ~isfield(param,'chantype'), 		param.chantype = 'grad1';        end
if ~isfield(param,'foi'), 			param.foi = freqs;               end
if ~isfield(param,'title'),         param.title = '';                end
%if ~isfield(param,'condition'),     param.condition = 1;             end

% Convert to fieldtrip format

if strcmp(datafield,'timefreq')
	data = ts_data2fieldtrip(data, 'dimord', 'chan_freq_time', 'chantype',param.chantype);
else
	data = ts_data2fieldtrip(data,'chantype',param.chantype);
	data.avg = permute(data.avg,[1 3 2]);
end

% Plot all channels in data if channel_index is not specified
if ~isfield(param,'channel_index'), param.channel_index=1:size(data.(zparam),1); end

% Topoplot sequence parameters
toilim 	= param.toilim;
foi 	  = param.foi;
chidx   = param.channel_index; 
nr 			= param.nr;
nc 			= param.nc;

% Statistical masking parameters
showhighlight = param.showhighlight;
showmask      = param.showmask; 

dataplot 						= data;

try
    dataplot.grad.pnt 	= data.grad.pnt(chidx,:);
    dataplot.grad.label = data.grad.label(chidx);
end

% Set up statistical mask
if showmask || showhighlight
	if ~isfield(param,'mask'), error('Must supply mask when using showmask or showhighlight'); end
    % convert to fieldtrip dimensions
    try param.mask = permute(param.mask,[1 3 2]); end
	if ~isfield(param,'stat_time') && size(param.mask,3)~=length(data.time)
		error('Must specify stat_time when different from data.time');
	elseif ~isfield(param,'stat_time') && size(param.mask,3)==length(data.time)
		param.stat_time = data.time;
	end
	% get toi from stat_data
	t1 = find(param.stat_time(1)==data.time);
	t2 = find(param.stat_time(end)==data.time);
	dataplot.time 			= data.time(t1:t2);
	dataplot.(zparam) 	= data.(zparam)(:,:,t1:t2);
	clustermaskplot = param.mask(chidx,:,:);

    %if length(foi)>1
        % fix mismatch between mask freq indices and the data
        if max(foi)>size(param.mask,2) && length(foi)<=size(param.mask,2)
            % mask freqs start at min foi
            mask_fidx = 1:length(foi);
        elseif max(foi)>size(param.mask,2)
            error('Mask has too few frequencies for the frequencies of interest.');
        elseif ~isfield(param,'mask_fidx')
            %fprintf('Assuming mask freqs start at 1Hz and setting foi=1:%d\n',size(param.mask,2));
            %fprintf('If this is not true, specify mask indices with option "mask_fidx"\n');
            %foi = 1:size(param.mask,2);
            mask_fidx = foi;
        end
    %end
    %if length(foi)>1 && strcmpi(param.freqmask,'matchany')
    if strcmpi(param.freqmask,'matchany')    
		% average over frequency band (make entire band significant if any freq is)
		for ch=1:size(clustermaskplot,1)
			tt=find(sum(squeeze(clustermaskplot(ch,mask_fidx,:))==1));
			if ~isempty(tt), clustermaskplot(ch,mask_fidx,tt)=1; end
		end
	%elseif length(foi)>1 && strcmpi(param.freqmask,'matchall')
    elseif strcmpi(param.freqmask,'matchall')
        % average over frequency band (make entire band significant if all freqs are)
		for ch=1:size(clustermaskplot,1)
			maskvals = clustermaskplot(ch,mask_fidx,:);
			clustermaskplot(ch,mask_fidx,:) = zeros(length(ch),length(foi),size(clustermaskplot,3));
			if all(sum(squeeze(maskvals)'))
				% all freqs are significant at some times in this channel
				clear ttidx; sharedtimes = 1; sharedtidx = [];
				for f=1:length(mask_fidx)
					ttidx{f}=find(squeeze(maskvals(1,f,:)));
					if isempty(ttidx{f}), sharedtimes = 0;
					elseif sharedtimes && f==1, sharedtidx = ttidx{f};
					elseif sharedtimes, sharedtidx = intersect(sharedtidx,ttidx{f});
					end
				end
				if sharedtimes,clustermaskplot(ch,mask_fidx,sharedtidx) = 1; end;
			end
		end
	end
end
%mask_fidx = find(freqs==min(foi)):find(freqs==max(foi));
if ~exist('mask_fidx','var'), mask_fidx=1; end

% clip toi if stat interval is narrower than toi
if dataplot.time(1)   > toilim(1), toilim(1) = dataplot.time(1);   end;
if dataplot.time(end) < toilim(2), toilim(2) = dataplot.time(end); end;

% times that will be plotted
if nc*nr~=1
    di = floor((nearest(dataplot.time,toilim(2)) - nearest(dataplot.time,toilim(1))) / (nc*nr-1));
else
    di = 1;
end

tidx = [1:di:di*nr*nc] + nearest(dataplot.time,toilim(1));
if length(tidx) > nc*nr, tidx = tidx(1:nc*nr); end;
toi 	= [dataplot.time(tidx)]; 	% the times that will actually be plotted

cfg = [];
if isfield(param,'colormap'), cfg.colormap = param.colormap; end;
if isfield(param,'style'), 		cfg.style 	 = param.style;		 end;
if isfield(param,'electrodes'), cfg.electrodes = param.electrodes; end % {'on','off',default: 'on'}
cfg.comment = 'xlim';
cfg.commentpos = 'title';
%	cfg.showlabels = 'yes';

% set zlims
if ~isfield(param,'zscalefactor'), param.zscalefactor=1; end
if strcmp(param.zlim,'maxmin') || strcmp(param.zlim,'sym')
	fidx = nearest(dataplot.freq,foi(1)):nearest(dataplot.freq,foi(end));
%	zmin = min(min(min(min(dataplot.(zparam)(:,fidx,tidx))))); % channels already limited by ts_data2fieldtrip
%	zmax = max(max(max(max(dataplot.(zparam)(:,fidx,tidx)))));
	zmin = min(min(min(min(mean(dataplot.(zparam)(:,fidx,tidx),2))))); % channels already limited by ts_data2fieldtrip
	zmax = max(max(max(max(mean(dataplot.(zparam)(:,fidx,tidx),2))))); % averaged over frequency band
	if strcmp(param.zlim,'sym')
		% make limits symmetric around zero (useful for plotting power z-scores)
		zmax = max(abs(zmin),abs(zmax));	
		zmin = -zmax;
	end
	cfg.zlim = [zmin zmax]*param.zscalefactor;
else
	cfg.zlim = [param.zlim(1) param.zlim(2)]*param.zscalefactor;
end

figure('Name',param.title);
warning off;
for k=1:nr*nc
	subplot(nr,nc,k)
	cfg.xlim = [toi(k) toi(k)];
	cfg.ylim = [foi(1) foi(end)];
	if strcmp(datafield,'averages')
		cfg.ylim = cfg.zlim;
		dataplot.(zparam) = squeeze(dataplot.(zparam));
	end
  if showhighlight,	
		[chans jx jy]=find(clustermaskplot(:,mask_fidx,tidx(k))); 
		cfg.highlight=unique(chans); 
		if isfield(param,'hlmarkersize'),cfg.hlmarkersize=param.hlmarkersize; else cfg.hlmarkersize=6; end
		if isfield(param,'hllinewidth'),cfg.hllinewidth=param.hllinewidth; else cfg.hllinewidth=3; end
		if isfield(param,'hlmarker'),cfg.hlmarker=param.hlmarker; else cfg.hlmarker='o'; end
	end
  if showmask, cfg.mask = squeeze(clustermaskplot(:,mask_fidx(1),tidx(k))); end;

	ts_topoplotTFR(cfg,dataplot);      
    if k==1 && isfield(param,'colorbar'), colorbar; end;
end
warning on;


%filepath = '/home/jsherfey/matlab/development/MEG_TFstats/control_023/images';
%filename = 'proc_tf_epoch_data_cond2.3_effect_chan001-100.mat';
%fullpath = sprintf('%s/%s',filepath,filename); 
%print(gcf,'-djpeg',fullpath);
