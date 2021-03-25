function ts_multiplot(data,cfg)

% data: cell array of timesurfer data structures
% cfg
%   conditions - length = 1 or # data structures (default: 1)
%   print - save each figure
%   figpath - path to save figures
%   fignames - names of files
%   prefix - figures will be saved with name = prefix_#
%
% example: 
% cfg=[]; cfg.conditions=[1 2]; ts_multiplot(data,cfg) will save plots for
% both conditions
% cfg=[]; cfg.conditions=1; 

if ~exist('cfg','var'),       cfg=[];               end
if ~isfield(cfg,'print'),     cfg.print=0;          end
if ~isfield(cfg,'chantype'),  cfg.chantype='meg';   end
if isfield(cfg,'conditions'), c=cfg.conditions;  cfg=rmfield(cfg,'conditions'); else c=1; end

% set up the data
if ~iscell(data)
  tmpfield  = ts_objecttype(data);
  n = length(data.(tmpfield));
  if n==1
    data={data};
  else
    data_ = data; clear data;
    for k=1:n
      data{k} = ts_get_header(data_,'condition',k);
      data{k}.(tmpfield) = data_.(tmpfield)(k);
    end
    data = data(c);
    clear data_;
  end
end
n = length(data);

if length(c)==1, c=repmat(c,1,n); elseif length(c)~=n, c=ones(1,n); end
  
for s=1:length(data)

    % get data type info (epochs, averages, timefreq)
    datafield = ts_objecttype(data{s});

    % select one condition to plot
    c_ = c(s);

    % convert to fieldtrip
    dat = ts_data2fieldtrip(data{s},'chantype',cfg.chantype);

    % call the appropriate fieldtrip multiplot function
    switch datafield
      case 'timefreq'
        figure; multiplotTFR(cfg,dat);  clear dat;
      otherwise
        figure; multiplotER(cfg,dat); clear dat;
    end

    if cfg.print, printfig(cfg); end
end

function printfig(cfg)
% get filepath
if ~isfield(cfg,'figpath'), cfg.figpath=pwd; end
 
if ~exist(cfg.figpath,'dir')
  fprintf('%s: WARNING: creating %s...\n',mfilename,cfg.figpath);
  [s,m] = mkdir(fpath);
  if ~s
    error('failed to create %s:\n%s',cfg.figpath,m);
  end;
end;

% get filename
msk=''; if (isfield(cfg,'maskparameter') && isfield(data{s},cfg.maskparameter)) || isfield(data{s},'mask'), msk='_mask'; end
if ~isfield(cfg,'fignames') && ~isfield(cfg,'prefix'), figname = [sprintf('figure_%s_%d',datafield,s) msk];
elseif isfield(cfg,'prefix'), figname = [sprintf('%s_%s_%d',cfg.prefix,datafield,s) msk];
else figname = cfg.fignames{s};
end
keyboard
% save figure
print(gcf,'-djpeg',sprintf('%s/%s',figpath,figname));

end
end
