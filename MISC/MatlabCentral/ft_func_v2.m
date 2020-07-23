%% Function to run ft functions but keep your cfg file constant
% Created by Erik Kaestner (ekaestne@ucsd.edu) 4/3/14

function ndat = ft_func_v2(func,cfg,dat,dat2)



if nargin == 2 % Typical input for loading data (function, cfg), however sometimes cfg will be data

    % Included in case need to work on only a ft data structure with no cfg
    if isfield(cfg,'trial'); ocfg = cfg.cfg; end
        
    ndat = feval(func,cfg);
    
    % Included in case need to work on only a ft data structure with no cfg
    if isfield(cfg,'trial'); ndat.cfg = catstruct(ocfg,ndat.cfg);
    elseif isfield(cfg,'task')
        ndat.cfg.task = cfg.task;
        ndat.cfg.subject = cfg.subject;
        ndat.cfg.orig_label = ndat.label;
    elseif isfield(cfg,'continuous')
        ndat.cfg.orig_label = ndat.label;
    end
    
elseif nargin == 3 % Typical input (function, cfg, and dataset) for fieldtrip
   
    if isfield(cfg,'empty') && strcmpi(cfg.empty,'yes') % Added for plotting functions that do not require an output
        feval(func,cfg,dat);
    else
        ocfg = dat.cfg;
        ndat = feval(func,cfg,dat);
        ndat.cfg = catstruct(ocfg,ndat.cfg);
    end
    
    % Added to keep track of hilbert range     
    if isfield(cfg,'hilbert') && strcmpi(cfg.hilbert,'yes')
        ndat.cfg.hilbertbp = [70 190]; % - EJK need to change to take actual filtering of hilbert
    end
    
    % Added in to modify ndat.cfg.orig_trl in case of readjusting the time
    if isfield(cfg,'toilim')      
        n_begin = find(dat.time{1}==ndat.time{1}(1))-1;
        n_end   = numel(dat.time{1}) - find(dat.time{1}==ndat.time{1}(end));
        
        ndat.cfg.orig_trl      = dat.cfg.orig_trl(:,1) + n_begin;
        ndat.cfg.orig_trl(:,2) = dat.cfg.orig_trl(:,2) - n_end;
    end
    
elseif nargin == 4 && isfield(cfg,'dat_sup')
    
    if isfield(cfg,'empty') && strcmpi(cfg.empty,'yes') % Added for plotting functions that do not require an output
        feval(func,cfg,dat1,dat2);
    else
        ocfg = dat2.cfg;
        ndat = feval(func,cfg,dat,dat2);
        ndat.cfg = catstruct(ocfg,ndat.cfg);
    end
    
    % Added to keep track of hilbert range
    if isfield(cfg,'hilbert') && strcmpi(cfg.hilbert,'yes')
        ndat.cfg.hilbertbp = [70 190]; % - EJK need to change to take actual filtering of hilbert
    end
    
    % Added in to modify ndat.cfg.orig_trl in case of readjusting the time
    if isfield(cfg,'toilim')
        n_begin = find(dat2.time{1}==ndat.time{1}(1))-1;
        n_end   = numel(dat2.time{1}) - find(dat2.time{1}==ndat.time{1}(end));
        
        ndat.cfg.orig_trl      = dat2.cfg.orig_trl(:,1) + n_begin;
        ndat.cfg.orig_trl(:,2) = dat2.cfg.orig_trl(:,2) - n_end;
    end
    
else % Added if needed to concatnate structures (must concatenate only two at a time) or supply >1 dataset (example: ft_erpac)
    % - EJK PTP: badtrl
    
    ocfg  = dat.cfg;
    ocfg2 = dat2.cfg;
    
    ndat = feval(func,cfg,dat,dat2);
    ndat.cfg = catstruct(ocfg,ocfg2,ndat.cfg);
    
    % Add in additional variables to preserve in concatenation here (applies only if channels, not trials)
    if isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'channels') && isfield(ndat.cfg,'orig_label') ; ndat.cfg.orig_label = unique([ocfg.orig_label; ocfg2.orig_label]); end
    if isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'channels') && isfield(ndat.cfg,'channel'); ndat.cfg.channel = [ocfg.channel; ocfg2.channel]; end
    if isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'channels') && isfield(ndat.cfg,'badchan') && isfield(ocfg,'badchan') && isfield(ocfg2,'badchan'); ndat.cfg.badchan = unique([ocfg.badchan; ocfg2.badchan]); end
    
    % Add in additional variables to preserve in concatenation here (applies only if appending trials, not channels)
    if isfield(ndat.cfg,'methapp') && isfield(dat2,'sampleinfo'); ndat.sampleinfo = [dat.sampleinfo; [dat2.sampleinfo(:,1)+dat.sampleinfo(end,1),dat2.sampleinfo(:,2)+dat.sampleinfo(end,2)]]; end
    if isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'trials') && isfield(ndat.cfg,'orig_trl');ndat.cfg.orig_trl = [ocfg.orig_trl; ocfg2.orig_trl];end
    if isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'trials') && isfield(ndat.cfg,'trl');ndat.cfg.trl = [ocfg.trl; ocfg2.trl];end
    if isfield(cfg,'cmb_cnt_dat') && isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'trials') && isfield(ndat.cfg,'orig_datafile');ndat.cfg.datafile = [ocfg.datafile {ocfg2.datafile}];end
    if isfield(cfg,'cmb_cnt_dat') && isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'trials') && isfield(ndat.cfg,'orig_dataset'); ndat.cfg.dataset = [ocfg.dataset {ocfg2.dataset}];end
    if isfield(ndat.cfg,'methapp') && strcmpi(cfg.methapp,'trials') && isfield(ndat.cfg,'stim_identity'); ndat.cfg.stim_identity = {ocfg.stim_identity{:} ocfg2.stim_identity{:}};end
        
    % Added for the CV/VC task to keep track of which consonants & vowels go with which trials when concatenating
    if isfield(cfg,'task') && isfield(cfg,'methapp') && strcmpi(cfg.methapp,'trials') && (strcmpi(ndat.cfg.task,'VC') || strcmpi(ndat.cfg.task,'CV'))
        ndat.cfg.con = [ocfg.con; ocfg2.con];
        ndat.cfg.vow = [ocfg.vow; ocfg2.vow];
    end
    
end

if ~isfield(cfg,'empty') || (isfield(cfg,'empty') && ~strcmpi(cfg.empty,'yes'))  
    ndat = mmil_update_data(ndat);
end
    
end