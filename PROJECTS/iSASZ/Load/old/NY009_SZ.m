%% SZ NY09 Analysis Script
clear; clc;

plt_spc = 1;

%% Data Paths
subj = 'NY009_SZ';

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

indir = '/space/mdeh3/7/halgdev/analysis/iEEG_NYU/NY09/NY09_SZ/eeg';
 
cln_dir = {'CLIN1'};

inpath_holder.(cln_dir{1}) = strsplit(ls([indir '/' '*' 'SZ' '*.eeg']),'.eeg');
inpath = [inpath_holder.(cln_dir{1})];

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
sem_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-'}),inpath,spl_end,'uni',0);

%% Load Visual & Auditory data
inpath = cellfun(@(x) [x '.eeg'],inpath,'uni',0);

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
cfg.continuous = 'no';
cfg.dataset    = inpath;
sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);

cfg      = [];
cfg.task = 'SZ';
sem_dat  = ft_func(@ft_sasz_addevent,cfg,sem_dat);

%% Remove Unimportant Channels
rmv_chn = find(~cellfun(@isempty,strfind(sem_dat.(sem_dat.data_name{1}).label,'DC')))';
rmv_chn = [rmv_chn find(~cellfun(@isempty,strfind(sem_dat.(sem_dat.data_name{1}).label,'EKG')))'];
rmv_chn = numel(sem_dat.(sem_dat.data_name{1}).label)-3:numel(sem_dat.(sem_dat.data_name{1}).label);

if ~isempty(rmv_chn)
    cfg = []; 
    cfg.channel = ['all',strcat('-',sem_dat.(sem_dat.data_name{1}).label(rmv_chn))'];
    sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
end

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' ];
    cfg.prefix = ['grp_men_' sem_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sem_dat);
end

%% Remove identified Bad Channels
rmv_chn = [63 70 117];

if ~isempty(rmv_chn)
    cfg = []; 
    cfg.channel = ['all',strcat('-',sem_dat.(sem_dat.data_name{1}).label(rmv_chn))'];
    sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
end

%% Remove identified Noise problems
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq   = {[56 64] [118 122] [178 182]};
cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
sem_dat      = ft_func(@ft_preprocessing,cfg,sem_dat);

%% Filter Data for LFP & Baseline
cfg            = [];
cfg.data_new   = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 20;
cfg.new_suffix = 'lfp';
sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);

%% Filter Data for HGP & Baseline
cfg            = [];
cfg.data_name  = [1:numel(inpath)];
cfg.hilbert    = 'amp';
cfg.freq_band  = {[70 170]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'hgp';
sem_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,sem_dat);

cfg           = [];
cfg.data_name = [numel(inpath)*2+1:numel(inpath)*3];
cfg.lpfilter  = 'yes';
cfg.lpfreq    = 20;
sem_dat       = ft_func(@ft_preprocessing,cfg,sem_dat);

%% Filter Data for Theta
cfg            = [];
cfg.data_name  = [1:numel(inpath)];
cfg.hilbert    = 'amp';
cfg.freq_band  = {[3 7]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'tht';
sem_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,sem_dat);

cfg           = [];
cfg.data_name = numel(inpath)*3+1:numel(inpath)*4;
cfg.lpfilter  = 'yes';
cfg.lpfreq    = 20;
sem_dat       = ft_func(@ft_preprocessing,cfg,sem_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1:numel(inpath);
cfg.rmfield   = 'yes';
sem_dat       = ft_func([],cfg,sem_dat);

%% Baseline Data & Remove Padding
cfg = [];
cfg.latency = [-0.8 1.6];
sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);

cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.4 0];
sem_dat            = ft_func(@ft_preprocessing,cfg,sem_dat);

%% Automatic rejection & Apply Rejections
cfg = [];
cfg.measures  = {'time' 'variance'};
cfg.thresh    = [0.98 0.98];
cfg.outdir    = [outpath '/' 'rejection' '/' subj '/' ];
cfg.prefix    = sem_dat.data_name;
cfg.specific  = {'prefix';[1 2 3]};
cfg.pad       = 0.4;
cfg.plot      = 1;
sem_dat       = ft_func(@auto_rej,cfg,sem_dat);

cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
sem_dat     = ft_func(@ft_apply_rejection,cfg,sem_dat);

%% Remove Padding
cfg = [];
cfg.latency = [-0.4 1.2];
sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);

%% Add in SNR measures
cfg = [];
cfg.return_events = 0;
cfg.old_events  = num2cell(1:8);
cfg.new_events  = num2cell(1:8);
cfg.crt_alt_eve = 'trialinfo';
sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);

cfg = [];
cfg.return_events = 0;
cfg.old_events  = {1:8};
cfg.new_events  = {101};
cfg.crt_alt_eve = 'ovr_all_evt';
sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);

cfg = [];
cfg.nse_tme     = [-0.4 0];
cfg.sig_tme     = [0 0.8];
cfg.alt_eve_loc = {'ovr_all_evt'};
cfg.eve         = {101};
cfg.snr_lab     = {'vis'};
sem_dat         = ft_func(@ft_SNR,cfg,sem_dat);

%% Label Upkeep
has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);

try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' subj]);
has_loc = 1;
catch
    has_loc = 0;
end
if has_loc
    chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' subj]);
    sem_dat.(sem_dat.data_name{2}).cfg.alt_lab.label = ft_correct_channel(chn_loc,sem_dat.(sem_dat.data_name{1}).label);
    sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{2}).cfg.alt_lab.label;
    sem_dat.(sem_dat.data_name{3}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{2}).cfg.alt_lab.label;
else
    sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{1}).label;
    sem_dat.(sem_dat.data_name{2}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{2}).label;
    sem_dat.(sem_dat.data_name{3}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{3}).label;
end

for iA = 1:numel(sem_dat.data_name)
   if has_aud && has_vis; sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.aud_vis_snr_lab = regexprep(strcat(sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label,{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(1),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(1,:)' * 100) / 100),{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(2),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(2,:)' * 100) / 100))',' +',' '); end
   if has_vis; sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.vis_snr_lab = regexprep(strcat(sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label,{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(1),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(1,:)' * 100) / 100))',' +',' '); end
   if has_aud; sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.aud_snr_lab = regexprep(strcat(sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label,{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(2),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(2,:)' * 100) / 100))',' +',' '); end
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'sem_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sem_dat);