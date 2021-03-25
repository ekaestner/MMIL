clear; clc;
addpath /home/ekaestne/fieldtrip-20131031/;
ft_defaults

indir = '/space/mdeh3/9/halgdev/projects/mmilanguage/iSASZ/proc_data';
fls = subdir(fullfile(indir,'*proc.mat'));

for ifl = 1:numel(fls)
    
    % Make output folder
    subj = fls(ifl).name(strfind(fls(ifl).name,'NY'):strfind(fls(ifl).name,'NY')+4);
    if ~exist([indir '/' subj '/' 'manip'],'dir')
        mkdir([indir '/' subj '/' 'manip'])
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
    cfg.toilim = [-0.3 1.2];
    cfg.trials = 'all';
    dat = ft_func(@ft_redefinetrial,cfg,dat);
    
    %% Calculate Semantic Differences
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 3 5 7] [2 4 6 8] [11 13 15 17] [12 14 16 18]};
    cfg.new_events = {1 2 3 4};
    dat = ft_redefine_events(cfg,dat);
    
    cfg = []; % SA semantic stat check
    cfg.events = [3 4];
    stat_data_sa_semantic = ft_fdrstats(cfg,dat);
    save([indir '/' subj '/' 'stat_data_sa_semantic_' dtype '.mat'],'stat_data_sa_semantic')
    
    cfg = []; % SA semantic stat plot
    cfg.evs = [3 4];
    cfg.outpath = [indir '/' subj '/' 'manip' '/' 'sa_semantic_' dtype];
    cfg.deviation = 1;
    if isempty(cell2mat(strfind(dat.label,'1>')))
        cfg.subchan = 'all';
    else
        ind_grid  = find(~cellfun(@isempty,strfind(dat.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(dat.label,'2>')));
        cfg.subchan = {ind_grid ind_strip};
    end
    cfg.stat_data = stat_data_sa_semantic;
    cfg.col = {'g' 'r'};
    cfg.col_dev = {[0.3 0.7 0.3] [0.7 0.3 0.3]};
    ft_macro_plot(cfg,dat)
    
    cfg = []; % SZ semantic stat check
    cfg.events = [1 2];
    stat_data_sz_semantic = ft_fdrstats(cfg,dat);
    save([indir '/' subj '/' 'stat_data_sz_semantic_' dtype '.mat'],'stat_data_sz_semantic')
    
    cfg = []; % SZ semantic stat plot
    cfg.evs = [1 2];
    cfg.outpath = [indir '/' subj '/' 'manip' '/' 'sz_semantic_' dtype];
    cfg.deviation = 1;
    if isempty(cell2mat(strfind(dat.label,'1>')))
        cfg.subchan = 'all';
    else
        ind_grid  = find(~cellfun(@isempty,strfind(dat.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(dat.label,'2>')));
        cfg.subchan = {ind_grid ind_strip};
    end
    cfg.stat_data = stat_data_sz_semantic;
    cfg.col = {'g' 'r'};
   cfg.col_dev = {[0.3 0.7 0.3] [0.7 0.3 0.3]};
    ft_macro_plot(cfg,dat)
    
    cfg = [];
    cfg.return_events = 1;
    dat = ft_redefine_events(cfg,dat);
    
    %% Calculate Repeated Differences
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4] [5 6 7 8] [11 12 13 14] [15 16 17 18]};
    cfg.new_events = {1 2 3 4};
    dat = ft_redefine_events(cfg,dat);
    
    cfg = []; % SA repeated stat check
    cfg.events = [3 4];
    stat_data_sa_repeated = ft_fdrstats(cfg,dat);
    save([indir '/' subj '/' 'stat_data_sa_repeated_' dtype '.mat'],'stat_data_sa_repeated')
    
    cfg = []; % SA repeated stat plot
    cfg.evs = [3 4];
    cfg.outpath = [indir '/' subj '/' 'manip' '/' 'sa_repeated_' dtype];
    cfg.deviation = 1;
    if isempty(cell2mat(strfind(dat.label,'1>')))
        cfg.subchan = 'all';
    else
        ind_grid  = find(~cellfun(@isempty,strfind(dat.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(dat.label,'2>')));
        cfg.subchan = {ind_grid ind_strip};
    end
    cfg.stat_data = stat_data_sa_repeated;
    cfg.col = {'c' 'y'};
    cfg.col_dev = {[0.3 0.7 0.7] [0.7 0.7 0.3]};
    ft_macro_plot(cfg,dat)
    
    cfg = []; % SZ repeated stat check
    cfg.events = [1 2];
    stat_data_sz_repeated = ft_fdrstats(cfg,dat);
    save([indir '/' subj '/' 'stat_data_sz_repeated_' dtype '.mat'],'stat_data_sz_repeated')
    
    cfg = []; % SZ repeated stat plot
    cfg.evs = [1 2];
    cfg.outpath = [indir '/' subj '/' 'manip' '/' 'sz_repeated_' dtype];
    cfg.deviation = 1;
    if isempty(cell2mat(strfind(dat.label,'1>')))
        cfg.subchan = 'all';
    else
        ind_grid  = find(~cellfun(@isempty,strfind(dat.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(dat.label,'2>')));
        cfg.subchan = {ind_grid ind_strip};
    end
    cfg.stat_data = stat_data_sz_repeated;
    cfg.col = {'c' 'y'};
    cfg.col_dev = {[0.3 0.7 0.7] [0.7 0.7 0.3]};
    ft_macro_plot(cfg,dat)
    
    cfg = [];
    cfg.return_events = 1;
    dat = ft_redefine_events(cfg,dat);
    
    %% Calculate audvis Differences
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {3 4};
    dat = ft_redefine_events(cfg,dat);
    
    cfg = []; % SA audvis stat check
    cfg.events = [3 4];
    stat_data_audvis = ft_fdrstats(cfg,dat);
    save([indir '/' subj '/' 'stat_data_audvis_' dtype '.mat'],'stat_data_audvis')
    
    cfg = []; % SA audvis stat plot
    cfg.evs = [3 4];
    cfg.outpath = [indir '/' subj '/' 'manip' '/' 'sa_audvis_' dtype];
    cfg.deviation = 1;
    if isempty(cell2mat(strfind(dat.label,'1>')))
        cfg.subchan = 'all';
    else
        ind_grid  = find(~cellfun(@isempty,strfind(dat.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(dat.label,'2>')));
        cfg.subchan = {ind_grid ind_strip};
    end
    cfg.stat_data = stat_data_audvis;
    ft_macro_plot(cfg,dat)
    
    cfg = [];
    cfg.return_events = 1;
    dat = ft_redefine_events(cfg,dat);
    
end