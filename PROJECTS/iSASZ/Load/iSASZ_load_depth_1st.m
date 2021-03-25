function iSASZ_load_depth_1st(fcfg)

%% Initial Variables
subj = fcfg.sbj_nme;
outpath = fcfg.out_pth;

% In Script Switches
fprintf('Starting to Load Subject %s \n\n\n',subj)

%% Data Paths
indir       = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'indir'); %mmil_readtext([fcfg.clr_fld '/indir/' subj]);
cln_fld     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_fld');
cln_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cln_dir'); %mmil_readtext([fcfg.clr_fld '/cln_dir/' subj]);
end_dir     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir'); %mmil_readtext([fcfg.clr_fld '/end_dir/' subj]);
tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);

cfg = [];
cfg.indir   = indir;
cfg.cln_fld = cln_fld;
cfg.cln_dir = cln_dir;
cfg.end_dir = end_dir;
cfg.tsk     = tsk;
[inpath,sem_dat,trl] = mmil_find_files(cfg);

%% Load Visual & Auditory data
epc_tme      = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

for iL = 1:numel(sem_dat.data_name)
    
    if ~isempty(string_find(end_dir(iL),{'mat'}))
        ttt = load(inpath{iL});
        sem_dat.(sem_dat.data_name{iL}) = ttt.epoch_data;
        
        sem_dat.(sem_dat.data_name{iL}) = mmil_format_mat(sem_dat.(sem_dat.data_name{iL}));
        clear ttt
    elseif ~isempty(string_find(end_dir(iL),{'eeg'}))
        ttt = ts_load_data(inpath{iL},'accept_all_flag',1);
        sem_dat.(sem_dat.data_name{iL}) = ttt;
        
        sem_dat.(sem_dat.data_name{iL}) = mmil_format_mat(sem_dat.(sem_dat.data_name{iL}));
        clear ttt
    elseif ~isempty(string_find(end_dir(iL),{'set'}))
        ttt = ts_load_data(inpath{iL},'accept_all_flag',1);
        sem_dat.(sem_dat.data_name{iL}) = ttt;
        
        sem_dat.(sem_dat.data_name{iL}) = mmil_format_mat(sem_dat.(sem_dat.data_name{iL}));
        clear ttt
        
    elseif ~isempty(string_find(end_dir(iL),{'edf'}))
        
        cfg            = [];
        cfg.specific   = {'dataset';1:numel(inpath)};
        cfg.data_new   = 'yes';
        cfg.continuous = 'yes';
        cfg.dataset    = inpath;
        sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);
        
        if ~exist([fcfg.clr_fld '/trialfun_output/' subj '.mat'])
            
            trialfun       = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun'); %mmil_readtext([fcfg.clr_fld '/indir/' subj]);
            
            cfg = [];
            if numel(trialfun) == 1; cfg.trialfun = trialfun{1}; else cfg.trialfun = trialfun; end
            cfg.dataset  = inpath;
            if numel(trialfun) == 1; cfg.specific = {'dataset' ; 1:numel(inpath)}; else cfg.specific = {'dataset' 'trialfun' ; 1:numel(inpath) 1:numel(trialfun)}; end
            cfg.minduration = 0.500;
            cfg.pre         = epc_tme(1)-1;
            cfg.pos         = epc_tme(2)+1;
            cfg.evt         = 1:8;
            trl = ft_func(@ft_definetrial,cfg,trl);
            
            if any(strcmpi(subj,{'NY496_SA_SZ' 'NY503_SA_SZ' 'NY511_SA_SZ'}))
                trl = sasz_fix_trl(subj,trl);
            end
            
            trl_hld = cell(1,numel(inpath)); for iEP = 1:numel(inpath); trl_hld{iEP} = trl.(trl.data_name{iEP}).trl; end
            
            if any(strcmpi(subj,{'NY442_SA_SZ'})); trl_hld{5}(341:end,:) = []; end
            
            for i = 1:numel(tsk); if strcmpi(tsk{i},'SA'); trl_hld{i}(:,4) = trl_hld{i}(:,4)+10; end; end
            
            save([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
        else
            load([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
            if any(strcmpi(subj,{'NY511_SA_SZ'}));
                trl_hld{3}(145:end,:) = [];
            end
        end
        
        cfg = [];
        cfg.specific  = {'trl';1:numel(inpath)};
        cfg.trl       = trl_hld;
        sem_dat = ft_func(@ft_redefinetrial,cfg,sem_dat);
        
    end
    
end

%% Misc
for iL = 1:numel(sem_dat.data_name)
    if strcmpi(tsk{iL},'SA') && all(sem_dat.(sem_dat.data_name{iL}).trialinfo < 10)
        sem_dat.(sem_dat.data_name{iL}).trialinfo = sem_dat.(sem_dat.data_name{iL}).trialinfo + 10;
    end
end

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
    sem_dat.NY17_SZ.cfg.trl = [sem_dat.NY17_SZ.sampleinfo sem_dat.NY17_SZ.sampleinfo(:,2)];
elseif strcmpi(subj,'NY226_SA_SZ')
    cfg = [];
    cfg.latency = [-0.2 0.8];
    sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);
    
    sem_dat.(sem_dat.data_name{1}).fsample = upper(sem_dat.(sem_dat.data_name{1}).fsample);
    
    sem_dat.(sem_dat.data_name{1}).time = repmat({-0.2+0.002:0.002:0.8},1,numel(sem_dat.(sem_dat.data_name{1}).time));
    sem_dat.(sem_dat.data_name{2}).time = repmat({-0.2+0.002:0.002:0.8},1,numel(sem_dat.(sem_dat.data_name{2}).time));
    sem_dat.(sem_dat.data_name{3}).time = repmat({-0.2+0.002:0.002:0.8},1,numel(sem_dat.(sem_dat.data_name{3}).time));
    
end

if round(sem_dat.(sem_dat.data_name{1}).fsample) ~= sem_dat.(sem_dat.data_name{1}).fsample
    for iD = 1:numel(sem_dat.data_name);
        sem_dat.(sem_dat.data_name{iD}).fsample = round(sem_dat.(sem_dat.data_name{iD}).fsample);
        sem_dat.(sem_dat.data_name{iD}).time = repmat({sem_dat.(sem_dat.data_name{iD}).time{1}(1):1/sem_dat.(sem_dat.data_name{iD}).fsample:sem_dat.(sem_dat.data_name{iD}).time{1}(end)},1,numel(sem_dat.(sem_dat.data_name{iD}).time));
    end
end

%% First Combine
cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb'); cmb = cmb{1};

if isempty(cmb)
    inpath{:}
    disp(['Elaborate on which file(s) go together (cmb): ' ' '])
    pause();
    pause(5);
    cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb'); cmb = cmb{1};
end

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    sem_dat = mmil_combine_data2(cfg,sem_dat);
end

%% Remove identified Unimportant & Bad Channels
% identified channels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.rmv_chn = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rmv_chn'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';

sem_dat = mmil_remove_channels(cfg,sem_dat);

%% Combine Clinical/Day if necessary
cmb = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb_scd'); cmb = cmb{1};

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.sbj_nme = fcfg.sbj_nme;
    
    sem_dat = mmil_combine_data2(cfg,sem_dat);
end

%% Check Depth Existence
dpt_num = mmil_readtext([fcfg.clr_fld '/' 'ele_idn' '/' fcfg.sbj_nme '/' fcfg.sbj_nme '.csv']);

if sum(cell2mat(dpt_num(2:end,3))==2)~=0
    
    %% Fix Labels
    cfg = [];
    cfg.ovr_wrt = 0;
    cfg.sbj_nme = fcfg.sbj_nme;
    cfg.clr_fld = fcfg.clr_fld;
    cfg.chn_loc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'electrode_location');
    cfg.dpt = 1;
    sem_dat = mmil_fix_labels(cfg,sem_dat);
    
    %% Depth Bi-Polar
    cfg = [];
    sem_dat.(sem_dat.data_name{1}) = mmil_depth_bipolar(cfg,sem_dat.(sem_dat.data_name{1}));
    
    %% Remove identified noise problems
    bse_frq     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_frq'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';
    
    if ~isempty(bse_frq)
        cfg = [];
        cfg.bsfilter = 'yes';
        cfg.bsfreq   = bse_frq;
        cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
        sem_dat      = ft_func(@ft_preprocessing,cfg,sem_dat);
    end

    %% Filter Data for LFP
    cfg            = [];
    cfg.data_new   = 'yes';
    cfg.lpfilter   = 'yes';
    cfg.lpfreq     = 15;
    cfg.new_suffix = 'lfp';
    sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);
    
    %% Filter Data for HGP
    if ~strcmpi(fcfg.sbj_nme,'NY226_SA_SZ')
        cfg=[];
        cfg.data_name  = 1;
        cfg.data_new   = 'yes';
        cfg.new_suffix = 'hgp';
        cfg.foi    = [70 80 90 100 110 130 140 150 160 170];    %frequency of interest
        cfg.sf     = [repmat(10,1,numel(cfg.foi))]; %specific frequency
        cfg.width  = cfg.foi./10;
        cfg.gwidth = ones(size(cfg.foi))*pi; %wavelet
        cfg.keeptrials = 'yes';
        cfg.method = 'tfr';
        cfg.toi=sem_dat.(sem_dat.data_name{1}).time{1};
        cfg.keeptrials = 'yes';
        sem_dat        = ft_func(@mmil_hgp_freq_analysis,cfg,sem_dat);
    else
        cfg            = [];
        cfg.data_name  = 1;
        cfg.hilbert    = 'amp';
        cfg.freq_band  = {[70 170]};
        cfg.data_new   = 'yes';
        cfg.new_suffix = 'hgp_hlb';
        sem_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,sem_dat);
    end
    
    cfg             = [];
    cfg.data_name   = 3;
    HGP_smooth_msec = 0.025;
    w               = round(HGP_smooth_msec* sem_dat.(sem_dat.data_name{1}).fsample);
    gauss_x         = -w:1:w;
    gauss_y         = normpdf(gauss_x,0,round(0.016* sem_dat.(sem_dat.data_name{1}).fsample));
    cfg.window      = gauss_y/sum(gauss_y);
    sem_dat         = ft_func(@ft_window_data,cfg,sem_dat);
    
    %% Remove superfluous continuous data structures
    cfg = [];
    cfg.data_name = 1;
    cfg.rmfield   = 'yes';
    sem_dat       = ft_func([],cfg,sem_dat);
    
    %% Baseline Data
    bse_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_tme'); bse_tme = bse_tme{1};
    
    cfg = [];
    cfg.demean         = 'yes';
    cfg.baselinewindow = bse_tme;
    sem_dat            = ft_func(@ft_baseline,cfg,sem_dat);
    
    %% Automatic rejection & Apply Rejections
    epc_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
    rjt_plt = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rjt_plt'); rjt_plt = rjt_plt{1};
    
    cfg = [];
    cfg.measures  = {'time' 'time-1' 'variance'};
    cfg.thresh    = [0.98 0.98 0.98];
    cfg.outdir    = [outpath '/' 'rejection' '/' subj '/' ];
    cfg.prefix    = sem_dat.data_name;
    cfg.specific  = {'prefix';1:numel(sem_dat.data_name)};
    cfg.pad       = sem_dat.(sem_dat.data_name{1}).time{1}(end)-epc_tme(2);
    cfg.plot      = rjt_plt;
    sem_dat       = ft_func(@auto_rej,cfg,sem_dat);
    
    cfg = [];
    cfg.measure = 'all';
    cfg.apply   = 'ieeg';
    sem_dat     = ft_func(@ft_apply_rejection,cfg,sem_dat);
    
    %% Remove Padding
    epc_tme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
    
    cfg = [];
    cfg.latency = epc_tme;
    sem_dat     = ft_func(@ft_selectdata,cfg,sem_dat);
    
    %% Save Data
    cfg = [];
    cfg.str_nme  = 'sem_dat';
    cfg.save     = 'yes';
    cfg.filename =[outpath '/' subj '_overall_data_depth.mat'];
    ft_func([],cfg,sem_dat);
    
end