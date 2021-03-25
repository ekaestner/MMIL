function varargout = ts_MEGtopoplot(data,varargin)

% Required inputs:
%   timefreq - timesurfer timefreq_data structure
%   parms
%
% Optional inputs:
% 	parms
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
% Note: parms.nr and parms.nc must be specified to plot a sequence of topoplots

% Created: 10-Sep-2008 by Jason Sherfey
%
% Revision: 05-Oct-2008 by Jason Sherfey
% Modifed to take timesurfer structure as input
% Revision: 11-Nov-2008 by Jason Sherfey
% This is now the wrapper around the main function ts_MEGtopoplot_callFT

parms = mmil_args2parms(varargin,{...
        'zlim','maxmin',[],...
        'nr',5,[],...
        'nc',7,[],...
        'style','both',[],...
        'datafile',[],[],...
        'zscore',0,{1,0},...
        'zscalefactor',1,[],...
        'chantype','grad1',[],...
        'foi',[],[],...
        },false);

if iscell(parms.datafile), parms.datafile = parms.datafile{1}; end      
if ~isempty(parms.datafile) && exist(parms.datafile,'file')
  load(parms.datafile);
end
if isempty(parms.foi), parms=rmfield(parms,'foi'); end
% datatype = ts_objecttype(data); 
% zparam   = 'power';
[object,datatype,dataparam] = ts_object_info(data);
zparam = dataparam{1};

if ~strcmp(datatype,'timefreq'),    error('%s: Does not support %s',mfilename,datatype); end
if ~isfield(parms,'get_zlim'),      parms.get_zlim      = 0;                end
if ~isfield(parms,'zlim'), 			parms.zlim          = 'maxmin';         end
if ~isfield(parms,'zscalefactor'),  parms.zscalefactor  = 1;                end
if ~isfield(parms,'chantype'),      parms.chantype      = 'grad1';          end
if ~isfield(parms,'toilim'),        parms.toilim    = [data.timefreq(1).time(1) data.timefreq(1).time(end)]; end;
if ~isfield(parms,'event'),         parms.event     = data.(datatype)(1).event_code;     end
if ~isfield(parms,'foi'),           parms.foi       = data.timefreq(1).frequencies; 
elseif length(parms.foi)==2,        parms.foi       = [parms.foi(1):parms.foi(2)]; end
if ~isfield(parms,'condition'),     [jnk1 parms.condition jnk2] ...
                                      = intersect([data.(datatype).event_code],parms.event); end

conds   = parms.condition;
chidx   = find(strcmp(parms.chantype,{data.sensor_info.typestring}));
fidx    = nearest(data.timefreq(1).frequencies,parms.foi(1)):...
          nearest(data.timefreq(1).frequencies,parms.foi(end));
tidx    = nearest(data.timefreq(1).time,parms.toilim(1)):...
          nearest(data.timefreq(1).time,parms.toilim(end));

if parms.zscore
  for c = 1:length(conds)
    data.(datatype)(conds(c)).(zparam) = ts_zscore(data,'cond',cond(c));
  end
end

if strcmp(parms.zlim,'maxmin') || strcmp(parms.zlim,'sym')
    zmin = []; 
    zmax = [];
    for i = 1:length(conds)
        zmin = [zmin min(min(min(mean(data.(datatype)(conds(i)).(zparam)(chidx,tidx,fidx),3))))];
        zmax = [zmax max(max(max(mean(data.(datatype)(conds(i)).(zparam)(chidx,tidx,fidx),3))))];
    end
    zmin = min(zmin);
    zmax = max(zmax);
%    parms.zlim = [min(zmin) max(zmax)];
    if strcmp(parms.zlim,'sym'), 
        zmax = max(abs(zmin),abs(zmax));
        zmin = -zmax;
    end
    parms.zlim = [zmin zmax];
end
parms.zlim = parms.zlim*parms.zscalefactor;
if parms.get_zlim && nargout==1, varargout{1} = parms.zlim; return; end

tmp_data = rmfield(data,datatype);
for i = 1:length(conds)
   tmp_data.(datatype) = data.(datatype)(conds(i));
   ts_MEGtopoplot_callFT(parms,tmp_data);
end

