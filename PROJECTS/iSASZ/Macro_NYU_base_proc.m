clear; clc;
addpath /home/ekaestne/fieldtrip-20131031/;
ft_defaults

indir = '/space/mdeh3/9/halgdev/projects/mmilanguage/iSASZ/proc_data';
fls = subdir(fullfile(indir,'*proc.mat'));

for ifl = 1:numel(fls)
    
    % Make output folder
    subj = fls(ifl).name(strfind(fls(ifl).name,'NY'):strfind(fls(ifl).name,'NY')+4);
    if ~exist([indir '/' subj '/' 'base'],'dir')
        mkdir([indir '/' subj '/' 'base'])
    end
    dtype = strfind(fls(ifl).name,'hgp');
    if isempty(dtype)
        dtype = 'lfp';
    else
        dtype = 'hgp';
    end
    
    % Load Data
    dat = load(fls(ifl).name);
    dat_n = fieldnames(dat);
    dat = dat.(dat_n{1});
    
    % Preprocess Data
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [0.5 15];
    dat = ft_func(@ft_preprocessing,cfg,dat);
    
    cfg = [];
    cfg.toilim = [-0.35 1.2];
    cfg.trials = 'all';
    dat = ft_func(@ft_redefinetrial,cfg,dat);
        
    %% Calculate Baseline Differences
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {3 4};
    dat = ft_redefine_events(cfg,dat);
    
    cfg = []; % SA base stat check
    cfg.events = [4];
    cfg.basetime = [-0.35 0];
    stat_data_sa_base = ft_fdrstats(cfg,dat);
    save([indir '/' subj '/' 'stat_data_sa_base_' dtype '.mat'],'stat_data_sa_base')
    
    cfg = []; % SZ base stat check
    cfg.events = [3];
    cfg.basetime = [-0.35 0];
    stat_data_sz_base = ft_fdrstats(cfg,dat);
    save([indir '/' subj '/' 'stat_data_sz_base_' dtype '.mat'],'stat_data_sz_base')
    
    cfg = []; % SA base stat plot
    cfg.evs = [4];
    cfg.outpath = [indir '/' subj '/' 'base' '/' 'sa_base_' dtype];
    cfg.basetime = [-0.35 0];
    cfg.deviation = 1;
    if isempty(cell2mat(strfind(dat.label,'1>')))
        cfg.subchan = 'all';
    else
        ind_grid  = find(~cellfun(@isempty,strfind(dat.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(dat.label,'2>')));
        cfg.subchan = {ind_grid ind_strip};
    end
    cfg.stat_data = stat_data_sa_base;
    ft_macro_plot(cfg,dat)
    
    cfg = []; % SZ base stat plot
    cfg.evs = [3];
    cfg.outpath = [indir '/' subj '/' 'base' '/' 'sz_base_' dtype];
    cfg.basetime = [-0.35 0];
    cfg.deviation = 1;
    if isempty(cell2mat(strfind(dat.label,'1>')))
        cfg.subchan = 'all';
    else
        ind_grid  = find(~cellfun(@isempty,strfind(dat.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(dat.label,'2>')));
        cfg.subchan = {ind_grid ind_strip};
    end
    cfg.stat_data = stat_data_sz_base;
    ft_macro_plot(cfg,dat)
    
    cfg = [];
    cfg.return_events = 1;
    dat = ft_redefine_events(cfg,dat);
    
    keep fls indir
    
end