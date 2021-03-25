%% Initial Analysis script for iSA/SZ project
% This script will perform basic pre-processing & look for electrodes with
% a response to either SA, SZ, or SA/SZ

isub = 8;

indir = indir;
subj = subj{isub};
infile = infile{isub};
outdir = outdir;

fcfg            = [];
fcfg.notch      = 1;
fcfg.preprocess = 1;
fcfg.init_plot  = 1;
fcfg.reref      = 1;
fcfg.rej        = 1;
fcfg.baseline   = 1;
fcfg.stats      = 1;
fcfg.plot       = 1;

keep isub indir subj infile outdir fcfg

function Macro_Init_proc(indir,subj,infile,outdir,fcfg)

if ~exist([outdir '/' subj],'dir')
    mkdir([outdir '/' subj])
end

load([indir '/' subj '/' infile]);

epoch_data.cfg.orig_label = epoch_data.label;

%% Preprocess Data
if fcfg.notch
    
    % Fix Noise Issues
    cont = 0;
    rchan = round(rand(1,8)*numel(epoch_data.label)); rtrl = round(rand(1,8)*numel(epoch_data.trial));
    for ir = 1:length(rchan);
        figure()
        periodogram(epoch_data.trial{rtrl(ir)}(rchan(ir),:),[],[],epoch_data.fsample)
    end
    while ~cont
        ln_ns = input('Vector of Line Noise []: ');
        bs_freq = input('Vectors of Bandstop {[] [] ...}: ');
        if ~isempty(ln_ns)
            cfg = [];
            cfg.dftfilter = 'yes';
            cfg.dftfreq = ln_ns;
            epoch_data = ft_func(@ft_preprocessing,cfg,epoch_data);
        end
        if ~isempty(bs_freq)
            for ibs = 1:length(bs_freq)
                cfg = [];
                cfg.bsfilter = 'yes';
                cfg.bsfreq = bs_freq{ibs};
                epoch_data = ft_func(@ft_preprocessing,cfg,epoch_data);
            end
        end
        for ir = 1:length(rchan);
            figure()
            periodogram(epoch_data.trial{rtrl(ir)}(rchan(ir),:),[],[],epoch_data.fsample)
        end
        cont = input('Continue (0=no 1=yes)?: '); % - EJK Just a place holder for now
    end
    epoch_data_notch = epoch_data;
    save([outdir '/' subj '/epoch_data_notch.mat'],'epoch_data_notch','-v7.3')
    close all
else
    load([outdir '/' subj '/epoch_data_notch.mat']);
end

if fcfg.preprocess
    % Filter LFP Data
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [0.5 30];
    epoch_data_lfp = ft_func(@ft_preprocessing,cfg,epoch_data_notch);
    
    epoch_data_lfp.sampleinfo(find(epoch_data_lfp.trialinfo>=10,1,'last')+1:end,1) = epoch_data_lfp.sampleinfo(find(epoch_data_lfp.trialinfo>=10,1,'last')+1:end,1) + epoch_data_lfp.sampleinfo(find(epoch_data_lfp.trialinfo>=10,1,'last'),1);
    epoch_data_lfp.sampleinfo(find(epoch_data_lfp.trialinfo>=10,1,'last')+1:end,2) = epoch_data_lfp.sampleinfo(find(epoch_data_lfp.trialinfo>=10,1,'last')+1:end,2) + epoch_data_lfp.sampleinfo(find(epoch_data_lfp.trialinfo>=10,1,'last'),1);
    epoch_data_lfp.cfg.orig_trl = epoch_data_lfp.sampleinfo; % Added because loading .eeg does not have a 'trl' field
    
    % Filter HGP Data
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [75 165];
    epoch_data_hgp = ft_func(@ft_preprocessing,cfg,epoch_data_notch);
    
    epoch_data_hgp.sampleinfo(find(epoch_data_hgp.trialinfo>=10,1,'last')+1:end,1) = epoch_data_hgp.sampleinfo(find(epoch_data_hgp.trialinfo>=10,1,'last')+1:end,1) + epoch_data_hgp.sampleinfo(find(epoch_data_hgp.trialinfo>=10,1,'last'),1);
    epoch_data_hgp.sampleinfo(find(epoch_data_hgp.trialinfo>=10,1,'last')+1:end,2) = epoch_data_hgp.sampleinfo(find(epoch_data_hgp.trialinfo>=10,1,'last')+1:end,2) + epoch_data_hgp.sampleinfo(find(epoch_data_hgp.trialinfo>=10,1,'last'),1);
    epoch_data_hgp.cfg.orig_trl = epoch_data_hgp.sampleinfo; % Added because loading .eeg does not have a 'trl' field
    
    cfg = [];
    cfg.hilbert = 'abs';
    epoch_data_hgp = ft_func(@ft_preprocessing,cfg,epoch_data_hgp);
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [0.5 30];
    epoch_data_hgp = ft_func(@ft_preprocessing,cfg,epoch_data_hgp);
    
    % Trim time to avoid edge artifacts
    cfg = [];
    cfg.toilim = [-0.5 1.5];
    cfg.trials = 'all';
    epoch_data_lfp = ft_func(@ft_redefinetrial,cfg,epoch_data_lfp);
    epoch_data_hgp = ft_func(@ft_redefinetrial,cfg,epoch_data_hgp);
    
    % Save initially processed data
    save([outdir '/' subj '/epoch_data_lfp.mat'],'epoch_data_lfp','-v7.3')
    save([outdir '/' subj '/epoch_data_hgp.mat'],'epoch_data_hgp','-v7.3')
else
    load([outdir '/' subj '/epoch_data_lfp.mat'])
    load([outdir '/' subj '/epoch_data_hgp.mat'])
end

%% Reject Trials
if fcfg.rej

    cfg           = [];
    cfg.ieegscale = 1;
    cfg.method    = 'trial';
    epoch_data_lfp = ft_func(@ft_rejectvisual,cfg,epoch_data_lfp);
    
    cfg = [];
    cfg.save = 1;
    cfg.saveloc = [outdir '/' subj '/reject_data_lfp'];
    ft_savereject(cfg,epoch_data_lfp);
    
    cfg = [];
    cfg.load = 1;
    cfg.saveloc = [outdir '/' subj '/reject_data_lfp.mat'];
    epoch_data_hgp = ft_savereject(cfg,epoch_data_hgp);
    
    cfg          = [];
    cfg.ieegscale = 1;
    cfg.method   = 'trial';
    epoch_data_hgp = ft_func(@ft_rejectvisual,cfg,epoch_data_hgp);
    
    cfg = [];
    cfg.save = 1;
    cfg.saveloc = [outdir '/' subj '/reject_data_hgp'];
    ft_savereject(cfg,epoch_data_hgp);
    
    epoch_data_lfp.cfg.orig_stim_identity = epoch_data_lfp.cfg.stim_identity;
    epoch_data_hgp.cfg.orig_stim_identity = epoch_data_hgp.cfg.stim_identity;
    cfg = [];
    cfg.special = {'stim_identity'};
    epoch_data_lfp = ft_special_reject(cfg,epoch_data_lfp);
    epoch_data_hgp = ft_special_reject(cfg,epoch_data_hgp); 
else
    cfg = [];
    cfg.load = 1;
    cfg.saveloc = [outdir '/' subj '/reject_data_lfp.mat'];
    epoch_data_lfp = ft_savereject(cfg,epoch_data_lfp);
    
    cfg = [];
    cfg.load = 1;
    cfg.saveloc = [outdir '/' subj '/reject_data_hgp.mat'];
    epoch_data_hgp = ft_savereject(cfg,epoch_data_hgp);
    
    epoch_data_lfp.cfg.orig_stim_identity = epoch_data_lfp.cfg.stim_identity;
    epoch_data_hgp.cfg.orig_stim_identity = epoch_data_hgp.cfg.stim_identity;
    cfg = [];
    cfg.special = {'stim_identity'};
    epoch_data_lfp = ft_special_reject(cfg,epoch_data_lfp);
    epoch_data_hgp = ft_special_reject(cfg,epoch_data_hgp);
end

%% Initial Plot
if fcfg.init_plot
%     rem_chan_ind = [];
%     cont = 0;
%     while ~cont
        
        if ~exist([outdir '/' subj '/' 'init'],'dir')
            mkdir([outdir '/' subj '/' 'init'])
        end
        
%         if ~isempty(rem_chan_ind)
%             rem_chan = strcat('-',num2str(rem_chan_ind));
%             cfg = []; cfg.channel = {'all',rem_chan};
%             epoch_data_lfp = ft_func(@ft_preprocessing,cfg,epoch_data_lfp);
%             epoch_data_hgp = ft_func(@ft_preprocessing,cfg,epoch_data_hgp);
%         end
        
        % LFP
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
        cfg.new_events = {1 2};
        epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
        
        ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'2>')));
        
        cfg = [];
        cfg.evs = [1 2];
        cfg.outpath = [outdir '/' subj '/' 'init' '/' 'epoch_data_lfp'];
        cfg.deviation = 1;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
        ft_macro_plot(cfg,epoch_data_lfp)
        
        cfg = [];
        cfg.return_events = 1;
        epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
        
        % HGP
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
        cfg.new_events = {1 2};
        epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
        
        ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'1>')));
        ind_strip = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'2>')));
        
        cfg = [];
        cfg.evs = [1 2];
        cfg.outpath = [outdir '/' subj '/' 'init' '/' 'epoch_data_hgp'];
        cfg.deviation = 1;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
        ft_macro_plot(cfg,epoch_data_hgp)
        
        cfg = [];
        cfg.return_events = 1;
        epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
        
%         rem_chan_ind = input('Channels to remove from data ({'' ''...}, or {} to continue)?: '); 
%         if isempty(rem_chan_ind); cont = 1; 
%         else
%             rem_chan_ind = cellfun(@find,cellfun(@(x) strcmpi(epoch_data.label,x),rem_chan_ind,'UniformOutput',0));
%         end
%     end
end

%% Rereference
if fcfg.reref
    
    % LFP Rereference
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'2>')));
    evs = [1];
    chans = {64:length(epoch_data_lfp.label)};
    
    for iev = 1:length(evs)
        for ich = 1:length(chans)
            cfg = [];
            cfg.ev_ind = evs(iev);
            cfg.chan_ind = chans{ich};
            epoch_data_lfp = ft_reref(cfg,epoch_data_lfp);
            
        end
    end
    
    cfg = [];
    cfg.return_events = 1;
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    % LFP Plot
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    cfg = [];
    cfg.evs = [1 2];
    cfg.outpath = [outdir '/' subj '/' 'init' '/' 'epoch_data_lfp_reref'];
    cfg.deviation = 1;
%     cfg.subchan = {ind_grid ind_strip};
    cfg.subchan = 'all';
    ft_macro_plot(cfg,epoch_data_lfp)
    
    cfg = [];
    cfg.return_events = 1;
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    % HGP Rereference
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
        
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'2>')));
    evs = [1];
%     chans = {ind_grid ind_strip};
    chans = {[1:32] [33:61] [62:118]};
    
    for iev = 1:length(evs)
        for ich = 1:length(chans)
            cfg = [];
            cfg.ev_ind = evs(iev);
            cfg.chan_ind = chans{ich};
            epoch_data_hgp = ft_reref(cfg,epoch_data_hgp);
        end
    end
    
    cfg = [];
    cfg.return_events = 1;
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
    
    % HGP - plot
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
    
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'2>')));
    
    cfg = [];
    cfg.evs = [1 2];
    cfg.outpath = [outdir '/' subj '/' 'init' '/' 'epoch_data_hgp_reref'];
    cfg.deviation = 1;
%     cfg.subchan = {ind_grid ind_strip};
    cfg.subchan = 'all';
    ft_macro_plot(cfg,epoch_data_hgp)
    
    cfg = [];
    cfg.return_events = 1;
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
    
end

%% Baseline Data?
if fcfg.baseline
    cfg = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = [-0.5 0];
    
    epoch_data_lfp   = ft_func(@ft_preprocessing,cfg,epoch_data_lfp);
    epoch_data_hgp   = ft_func(@ft_preprocessing,cfg,epoch_data_hgp);
end

%% Statistics
if fcfg.stats
    % SA LFP
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    cfg = [];
    cfg.basetime = [-0.5 0];
    cfg.events = [2];
    stat_data_sa_lfp_baseline = ft_fdrstats(cfg,epoch_data_lfp);
    
    save([outdir '/' subj '/' 'stat_data_sa_lfp_baseline'],'stat_data_sa_lfp_baseline')
    
    % SZ LFP
    cfg = [];
    cfg.basetime = [-0.5 0];
    cfg.events = [1];
    stat_data_sz_lfp_baseline = ft_fdrstats(cfg,epoch_data_lfp);
    
    save([outdir '/' subj '/' 'stat_data_sz_lfp_baseline'],'stat_data_sz_lfp_baseline')
    
    % SA/SZ LFP    
    cfg = [];
    cfg.events = [1 2];
    stat_data_sa_v_sz_lfp = ft_fdrstats(cfg,epoch_data_lfp);
    
    save([outdir '/' subj '/' 'stat_data_sa_v_sz_lfp'],'stat_data_sa_v_sz_lfp')
        
    % SA HGP
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
    
    cfg = [];
    cfg.basetime = [-0.5 0];
    cfg.events = [1];
    stat_data_sa_hgp_baseline = ft_fdrstats(cfg,epoch_data_hgp);
    
    save([outdir '/' subj '/' 'stat_data_sa_hgp_baseline'],'stat_data_sa_hgp_baseline')
       
    % SZ HGP    
    cfg = [];
    cfg.basetime = [-0.5 0];
    cfg.events = [2];
    stat_data_sz_hgp_baseline = ft_fdrstats(cfg,epoch_data_hgp);
    
    save([outdir '/' subj '/' 'stat_data_sz_hgp_baseline'],'stat_data_sz_hgp_baseline')
    
    % SA/SZ HGP   
    cfg = [];
    cfg.events = [1 2];
    stat_data_sa_v_sz_hgp = ft_fdrstats(cfg,epoch_data_hgp);
    
    save([outdir '/' subj '/' 'stat_data_sa_v_sz_hgp'],'stat_data_sa_v_sz_hgp')
    
    %     Return original events
    cfg = [];
    cfg.return_events = 1;
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    cfg = [];
    cfg.return_events = 1;
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);    
else
    load([outdir '/' subj '/' 'stat_data_sa_lfp_baseline']);
    load([outdir '/' subj '/' 'stat_data_sz_lfp_baseline']);
    load([outdir '/' subj '/' 'stat_data_sa_v_sz_lfp']);
    load([outdir '/' subj '/' 'stat_data_sa_hgp_baseline']);
    load([outdir '/' subj '/' 'stat_data_sz_hgp_baseline']);
    load([outdir '/' subj '/' 'stat_data_sa_v_sz_hgp']);
    
end

%% Plotting
if fcfg.plot
        
    if ~exist([outdir '/' subj '/' 'allevents'],'dir')
        mkdir([outdir '/' subj '/' 'allevents'])
    end
    
    % SA LFP Plot
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'2>')));
    
    cfg = [];
    cfg.evs = [2];
    cfg.basetime = [-0.5 0];
    cfg.stat_data = stat_data_sa_lfp_baseline;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
    cfg.deviation = 1;
    cfg.outpath = [outdir '/' subj '/' 'allevents' '/' 'epoch_data_sa_lfp'];
    ft_macro_plot(cfg,epoch_data_lfp)
    
    % SZ LFP Plot    
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'2>')));
    
    cfg = [];
    cfg.evs = [1];
    cfg.basetime = [-0.5 0];
    cfg.stat_data = stat_data_sz_lfp_baseline;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
    cfg.deviation = 1;
    cfg.outpath = [outdir '/' subj '/' 'allevents' '/' 'epoch_data_sz_lfp'];
    ft_macro_plot(cfg,epoch_data_lfp)
       
    % SA/SZ LFP Plot    
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_lfp.label,'2>')));
    
    cfg = [];
    cfg.evs = [1 2];
    cfg.stat_data = stat_data_sa_v_sz_lfp;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
    cfg.deviation = 1;
    cfg.outpath = [outdir '/' subj '/' 'allevents' '/' 'epoch_data_sa_v_sz_lfp'];
    ft_macro_plot(cfg,epoch_data_lfp)
    
    % SA HGP Plot    
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events = {[1 2 3 4 5 6 7 8] [11 12 13 14 15 16 17 18]};
    cfg.new_events = {1 2};
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);
    
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'2>')));
    
    cfg = [];
    cfg.evs = [2];
    cfg.basetime = [-0.5 0];
    cfg.stat_data = stat_data_sa_hgp_baseline;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
    cfg.deviation = 1;
    cfg.outpath = [outdir '/' subj '/' 'allevents' '/' 'epoch_data_sa_hgp'];
    ft_macro_plot(cfg,epoch_data_hgp)
    
    % SZ HGP Plot  
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'2>')));
    
    cfg = [];
    cfg.evs = [1];
    cfg.basetime = [-0.5 0];
    cfg.stat_data = stat_data_sz_hgp_baseline;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
    cfg.deviation = 1;
    cfg.outpath = [outdir '/' subj '/' 'allevents' '/' 'epoch_data_sz_hgp'];
    ft_macro_plot(cfg,epoch_data_hgp)
      
    % SA/SZ HGP Plot  
    ind_grid  = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'1>')));
    ind_strip = find(~cellfun(@isempty,strfind(epoch_data_hgp.label,'2>')));
    
    cfg = [];
    cfg.evs = [1 2];
    cfg.stat_data = stat_data_sa_v_sz_hgp;
%         cfg.subchan = {ind_grid ind_strip};
        cfg.subchan = 'all';
    cfg.deviation = 1;
    cfg.outpath = [outdir '/' subj '/' 'allevents' '/' 'epoch_data_sa_v_sz_hgp'];
    ft_macro_plot(cfg,epoch_data_hgp)
    
    %     Return original events
    cfg = [];
    cfg.return_events = 1;
    epoch_data_lfp = ft_redefine_events(cfg,epoch_data_lfp);
    
    cfg = [];
    cfg.return_events = 1;
    epoch_data_hgp = ft_redefine_events(cfg,epoch_data_hgp);    
    
end

epoch_data_lfp_preproc = epoch_data_lfp;
epoch_data_hgp_preproc = epoch_data_hgp;

save([outdir '/' subj '/epoch_data_lfp_preproc.mat'],'epoch_data_lfp_preproc','-v7.3')
save([outdir '/' subj '/epoch_data_hgp_preproc.mat'],'epoch_data_hgp_preproc','-v7.3')

end
