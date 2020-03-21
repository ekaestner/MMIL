%% ft_window_data
% 
% cfg.win   = window to use
% cfg.width = width of window
% 

function ft_dat = ft_window_data(cfg,ft_dat)
%% Create User-Specified Window
if ~isfield(cfg,'window')
    cfg.window = window(cfg.win_typ,cfg.width)/cfg.width;
end

%% Convolve Window With Data

for iTR = 1:numel(ft_dat.trial)
    for iCH = 1:size(ft_dat.trial{iTR},1)
        ft_dat.trial{iTR}(iCH,:) = conv(ft_dat.trial{iTR}(iCH,:),cfg.window,'same');
    end
end

%% Add Window information to .cfg
if isfield(cfg,'win_typ'); ft_dat.cfg.window_type = cfg.win_typ; end
if isfield(cfg,'width');   ft_dat.cfg.window_width = cfg.width; end
ft_dat.cfg.window = cfg.window;

end