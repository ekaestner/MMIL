%% iSASZ Analysis Script
clear; clc;

sbj_num = 5;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

indir       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/indir/' subj]);  
cln_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/cln_dir/' subj]); 
end_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/end_dir/' subj]); 
tsk         = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/tsk/' subj]);     
chn_nse_grp = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/cmn_nse/' subj]),'uni',0)'; 
rmv_chn     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/rmv_chn/' subj]),'uni',0)'; 
bse_frq     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/bse_frq/' subj]),'uni',0)';

% In Script Switches 
plt_spc = 1;

%% Data Paths
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

for iIN = 1:numel(indir)
    inpath_holder.(cln_dir{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
end
inpath = [];
for iIN = 1:numel(indir)
    inpath = [inpath inpath_holder.(cln_dir{iIN})];
end

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

for iTS = 1:numel(tsk)
    cfg        = [];
    cfg.task   = tsk{iTS};
    cfg.inpath = inpath{iTS};
    skp{iTS}     = mmil_check_eeg(cfg);
end

cfg          = [];
cfg.task     = [tsk(:)];
cfg.skp      = skp;
cfg.rep_fix  = 1;
cfg.specific = {'task' 'skp'; 1:numel(tsk) 1:numel(skp)};
sem_dat      = ft_func(@ft_sasz_addevent,cfg,sem_dat);

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' ];
    cfg.prefix = [sem_dat.data_name(1:numel(inpath))];
    cfg.specific  = {'prefix'; 1:numel(inpath)};
    ft_func(@ft_plot_spectrum,cfg,sem_dat);
end

%% Remove noise through removing common noise
if ~isempty(chn_nse_grp{1})
    cfg             = [];
    cfg.chn_nse_grp = chn_nse_grp;   
    sem_dat = ft_func(@ft_remove_common_noise,cfg,sem_dat);
end

%% Remove identified Unimportant & Bad Channels
% unimportant channels
if ~isempty(rmv_chn)
    cfg = [];
    cfg.channel = ['all',strcat('-',sem_dat.(sem_dat.data_name{1}).label(rmv_chn))'];
    sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
end

%% Remove identified Noise problems
if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sem_dat      = ft_func(@ft_preprocessing,cfg,sem_dat);
end

%% Combine SA & SZ trials
if numel(tsk) > 1
    cfg           = [];
    cfg.data_name = [1 ; 2];
    cfg.data_new  = 'yes';
    cfg.methapp   = 'trials';
    sem_dat = ft_func(@ft_appenddata,cfg,sem_dat);
    
    cfg = [];
    cfg.rmfield   = 'yes';
    cfg.data_name = [2];
    sem_dat = ft_func([],cfg,sem_dat);
    
    inpath(2) = [];
end

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

%% Filter Data for Control Beta
cfg            = [];
cfg.data_name  = [1:numel(inpath)];
cfg.hilbert    = 'amp';
cfg.freq_band  = {[23 27]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'tht';
sem_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,sem_dat);

cfg           = [];
cfg.data_name = numel(inpath)*4+1:numel(inpath)*5;
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
cfg.old_events  = num2cell(unique(sem_dat.(sem_dat.data_name{1}).trialinfo));
cfg.new_events  = num2cell(unique(sem_dat.(sem_dat.data_name{1}).trialinfo));
cfg.crt_alt_eve = 'trialinfo';
sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);

has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);

if has_vid && has_aud ; old_eve = {1:8 11:18}; new_eve = {101 102}; lab = {'vis' 'aud'};
    elseif has_vid; old_eve = {1:8}; new_eve = {101};   lab = {'vis'}; 
    elseif has_aud; old_eve = {11:18}; new_eve = {102}; lab = {'aud'}; end

cfg = [];
cfg.return_events = 0;
cfg.old_events  = old_eve;
cfg.new_events  = new_eve;
cfg.crt_alt_eve = 'ovr_all_evt';
sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);

cfg = [];
cfg.nse_tme     = [-0.4 0];
cfg.sig_tme     = [0 1];
cfg.alt_eve_loc = repmat({'ovr_all_evt'},1,numel(lab));
cfg.eve         = new_eve;
cfg.snr_lab     = lab;
sem_dat         = ft_func(@ft_SNR,cfg,sem_dat);

%% Label Upkeep
try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/'  'electrode_location' '/' subj]);
has_loc = 1;
catch
    has_loc = 0;
end
if has_loc
    chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
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
