%% iSASZ Analysis Script
clear; clc;

sbj_num = 20;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

fprintf('Starting on Subject %s \n\n\n',subj)

indir       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/indir/' subj]);  
cln_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/cln_dir/' subj]); 
end_dir     = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/end_dir/' subj]); 
tsk         = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/tsk/' subj]);     
rmv_chn     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/rmv_chn/' subj]),'uni',0)';
bse_frq     = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/bse_frq/' subj]),'uni',0)';
if any(~cellfun(@isempty,strfind(end_dir,'edf'))); trl_fun = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/trialfun/' subj]); end

chn_nse_grp = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/cmn_nse/' subj]);
tmp = strfind(chn_nse_grp,';');
if any(~isempty(tmp{:}))
    for i = 1:numel(chn_nse_grp)
        chn_nse_grp_dat_nme(i) = str2num(chn_nse_grp{i}(tmp{1}+1:end))';
        chn_nse_grp{i} = str2num(chn_nse_grp{i}(1:tmp{1}-1))'; 
    end
else
    chn_nse_grp = cellfun(@str2num,mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/cmn_nse/' subj]),'uni',0)';
end

% In Script Switches 
plt_spc = 1;

%% Data Paths
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

idn_fle = {};
for iIN = 1:numel(indir)
    if any(~cellfun(@isempty,strfind(end_dir,'eeg')));
        inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN}); if numel(inpath_holder.(tsk{iIN})) > 1; inpath_holder.(tsk{iIN})(2) = []; end;
    elseif any(~cellfun(@isempty,strfind(end_dir,'edf')));
        inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
        num_fle.(tsk{iIN}) = numel(inpath_holder.(tsk{iIN}));
        for i = 1:num_fle.(tsk{iIN}); idn_fle{end+1} = tsk{iIN}; end
    elseif any(~cellfun(@isempty,strfind(end_dir,'set')));
        inpath_holder.(tsk{iIN}) = strsplit(ls([indir{iIN} '/' '*' cln_dir{iIN} '*' end_dir{iIN}]),end_dir{iIN});
        num_fle.(tsk{iIN}) = numel(inpath_holder.(tsk{iIN}));
        for i = 1:num_fle.(tsk{iIN}); idn_fle{end+1} = tsk{iIN}; end
    end
end

inpath = [];
for iIN = 1:numel(indir)
    inpath = [inpath inpath_holder.(tsk{iIN})];
end

% Data_names and initialize data holder
spl_end           = strfind(inpath,'/');
sem_dat.data_name = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

trl.data_name   = cellfun(@(x,y) mmil_spec_char(x(y(end)+1:end),{'-' '.'}),inpath,spl_end,'uni',0);

%% Load Visual & Auditory data
if any(~cellfun(@isempty,strfind(end_dir,'eeg'))); inpath = cellfun(@(x) [x end_dir{1}],inpath,'uni',0); end
if any(~cellfun(@isempty,strfind(end_dir,'edf'))); inpath = cellfun(@(x) [x end_dir{1}],inpath,'uni',0); end
if any(~cellfun(@isempty,strfind(end_dir,'set'))); inpath = cellfun(@(x) [x end_dir{1}],inpath,'uni',0); end

cfg            = [];
cfg.specific   = {'dataset';1:numel(inpath)};
cfg.data_new   = 'yes';
if any(~cellfun(@isempty,strfind(end_dir,'set')) | ~cellfun(@isempty,strfind(end_dir,'edf'))); cfg.continuous = 'yes'; else cfg.continuous = 'no'; end
cfg.dataset    = inpath;
sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);

if any(~cellfun(@isempty,strfind(end_dir,'eeg')));
    
    for iTS = 1:numel(tsk)
        cfg        = [];
        cfg.task   = tsk{iTS};
        cfg.inpath = inpath{iTS};
        skp{iTS}     = mmil_check_eeg(cfg);
    end
    
    if strcmpi(subj,'NY017_SA_SZ'); cfg = []; cfg.trials = 1:538; sem_dat.(sem_dat.data_name{1}) = ft_preprocessing(cfg,sem_dat.(sem_dat.data_name{1})); end
    if strcmpi(subj,'NY019_SA_SZ'); cfg = []; cfg.trials = [1:338 340:400]; sem_dat.(sem_dat.data_name{1}) = ft_preprocessing(cfg,sem_dat.(sem_dat.data_name{1})); end
    
    cfg          = [];
    cfg.task     = [tsk(:)];
    cfg.skp      = skp;
    cfg.rep_fix  = 1;
    cfg.specific = {'task' 'skp'; 1:numel(tsk) 1:numel(skp)};
    sem_dat      = ft_func(@ft_sasz_addevent,cfg,sem_dat);
    
elseif any(~cellfun(@isempty,strfind(end_dir,'edf')));
    
    if ~exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'])
        cfg = [];
        cfg.trialfun = trl_fun{1};
        cfg.dataset  = inpath;
        cfg.specific = {'dataset' ; 1:numel(inpath)};
        cfg.minduration = 0.500;
        cfg.pre         = 1.5;
        cfg.pos         = 2.5;
        cfg.evt         = 1:8;
        trl = ft_func(@ft_definetrial,cfg,trl);
        
        trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end 
        
        for i = 1:numel(idn_fle); if strcmpi(idn_fle{i},'SA'); trl_hld{i}(:,4) = trl_hld{i}(:,4)+10; end; end
        
        tmp = input('Satisfied with defined trials?: '); clear tmp
        
        save(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
    else
        load(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
    end
    
    cfg = [];
    cfg.specific  = {'trl';1:numel(inpath)};
    cfg.trl       = trl_hld;
    sem_dat = ft_func(@ft_redefinetrial,cfg,sem_dat);
    
    % Combine as necessary
    for iFL = 1:numel(tsk)
        if num_fle.(tsk{iFL}) > 1
            
            idn = find(strcmpi(idn_fle,tsk{iFL}));
            
            cfg           = [];
            cfg.data_name = [repmat(idn(1),1,numel(idn)-1) ; idn(2):idn(end)];
            cfg.data_new  = 'yes';
            cfg.methapp   = 'trials';
            sem_dat = ft_func(@ft_appenddata,cfg,sem_dat);
            
            cfg = [];
            cfg.rmfield   = 'yes';
            cfg.data_name = idn(2):idn(end);
            sem_dat = ft_func([],cfg,sem_dat);
            
        end
    end
    
elseif any(~cellfun(@isempty,strfind(end_dir,'set')));
    
    if ~exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'])
        for iEP = 1:numel(inpath)
            cfg = [];
            cfg.indir    = inpath{iEP};
            cfg.tsk      = idn_fle{iEP};
            cfg.pre      = 0.95;
            cfg.pos      = 1.95;
            cfg.smp_rte  = sem_dat.(sem_dat.data_name{iEP}).fsample;
            trl_hld{iEP} = mmil_epoch_dat_set(cfg);
            
            for i = 1:numel(idn_fle); if strcmpi(idn_fle{i},'SA'); trl_hld{i}(:,4) = trl_hld{i}(:,4)+10; end; end
            
        end
        
        save(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
    else
        load(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/trialfun_output/' subj '.mat'],'trl_hld');
    end
    
    cfg = [];
    cfg.specific  = {'trl'; 1:numel(inpath)};
    cfg.trl       = trl_hld;
    sem_dat = ft_func(@ft_redefinetrial,cfg,sem_dat);
    
    % Combine as necessary
    for iFL = 1:numel(tsk)
        if num_fle.(tsk{iFL}) > 1
            
            idn = find(strcmpi(idn_fle,tsk{iFL}));
            
            cfg           = [];
            cfg.data_name = [repmat(idn(1),1,numel(idn)-1) ; idn(2):idn(end)];
            cfg.data_new  = 'yes';
            cfg.methapp   = 'trials';
            sem_dat = ft_func(@ft_appenddata,cfg,sem_dat);
            
            cfg = [];
            cfg.rmfield   = 'yes';
            cfg.data_name = idn(2):idn(end);
            sem_dat = ft_func([],cfg,sem_dat);
            
        end
    end
    
end

%% Misc
if strcmpi(subj,'NY008_SA_SZ')
    sem_dat.(sem_dat.data_name{2}).label = sem_dat.(sem_dat.data_name{1}).label;
elseif strcmpi(subj,'NY011_SA_SZ')
    sem_dat.(sem_dat.data_name{2}).label = sem_dat.(sem_dat.data_name{1}).label;
    
    cfg = [];
    cfg.latency = [-1.4 2.4];
    sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);
elseif strcmpi(subj,'NY068_SA')
    sem_dat.NY68_SAo.trialinfo = sem_dat.NY68_SAo.trialinfo';
    sem_dat.NY68_SAo.trial = sem_dat.NY68_SAo.trial';
elseif strcmpi(subj,'NY017_SA_SZ')
    sem_dat.NY17_SA.cfg.trl = [sem_dat.NY17_SA.sampleinfo sem_dat.NY17_SA.sampleinfo(:,2)];
end

%% Remove identified Unimportant & Bad Channels
% unimportant channels
if any(~cellfun(@isempty,rmv_chn))
    if numel(rmv_chn)==1
        cfg = [];
        cfg.channel = ['all',strcat('-',sem_dat.(sem_dat.data_name{1}).label(rmv_chn{1}))'];
        sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
    else
        for iRM = 1:numel(rmv_chn); tmp_rmv{iRM} = ['all',strcat('-',sem_dat.(sem_dat.data_name{iRM}).label(rmv_chn{iRM}))']; end
        cfg = [];
        cfg.channel = tmp_rmv;
        cfg.specific = {'channel' ; 1:numel(tmp_rmv)};
        sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
    end
    
end

%% Remove noise through removing common noise
if ~isvar('chn_nse_grp_dat_nme'); chn_nse_grp_dat_nme = repmat({1:numel(tsk)},1,numel(chn_nse_grp)); end
for i = 1:numel(chn_nse_grp)
    if ~isempty(chn_nse_grp{i})
        cfg             = [];
        cfg.data_name   = chn_nse_grp_dat_nme{1};
        cfg.chn_nse_grp = chn_nse_grp(i);
        sem_dat = ft_func(@ft_remove_common_noise,cfg,sem_dat);
    end
end

%% Remove identified Noise problems
if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sem_dat      = ft_func(@ft_preprocessing,cfg,sem_dat);
end

%% Examine Data for Noise
if plt_spc == 1;
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' ];
    cfg.prefix = [sem_dat.data_name(1:numel(tsk))];
    cfg.specific  = {'prefix'; 1:numel(tsk)};
    ft_func(@ft_plot_spectrum,cfg,sem_dat);
end

tmp = input([subj ' : ' 'Satisfied with spectrum?: ']); clear tmp

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
        
end

inpath(2:end) = [];

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
cfg.lpfreq    = 15;
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
cfg.new_suffix = 'bta';
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
cfg.specific  = {'prefix';1:numel(sem_dat.data_name)};
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

if has_vis && has_aud ; old_eve = {1:8 11:18}; new_eve = {101 102}; lab = {'vis' 'aud'};
    elseif has_vis; old_eve = {1:8}; new_eve = {101};   lab = {'vis'}; 
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

%% Keep track of meta-data
cfg = [];
cfg.sbj     = subj;
cfg.has_vis = has_vis;
cfg.has_aud = has_aud;
cfg.fle_loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';
mmil_sasz_log_event(cfg,sem_dat.(sem_dat.data_name{1}));

%% Label Upkeep
try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/'  'electrode_location' '/' subj]);
has_loc = 1;
catch
    has_loc = 0;
end
if has_loc
    chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
    sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label = ft_correct_channel(chn_loc,sem_dat.(sem_dat.data_name{1}).label);
    for iA = 2:numel(sem_dat.data_name); sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{1}).cfg.alt_lab.label; end
else
    for iA = 1:numel(sem_dat.data_name); sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{iA}).label; end
end

for iA = 1:numel(sem_dat.data_name)
   if has_aud && has_vis; sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.aud_vis_snr_lab = regexprep(strcat(sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label,{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(1),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(1,:)' * 100) / 100),{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(2),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(2,:)' * 100) / 100))',' +',' '); end
   if has_vis; sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.vis_snr_lab = regexprep(strcat(sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label,{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(1),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(1,:)' * 100) / 100))',' +',' '); end
   if has_aud; sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.aud_snr_lab = regexprep(strcat(sem_dat.(sem_dat.data_name{iA}).cfg.alt_lab.label,{' '},sem_dat.(sem_dat.data_name{iA}).cfg.snr_lab(has_vis+has_aud),{' '},num2str(round(sem_dat.(sem_dat.data_name{iA}).cfg.snr(has_vis+has_aud,:)' * 100) / 100))',' +',' '); end
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'sem_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sem_dat);
