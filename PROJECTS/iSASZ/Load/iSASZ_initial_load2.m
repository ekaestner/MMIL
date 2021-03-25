function iSASZ_initial_load2(fcfg)

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
tot_tsk     = unique(tsk);

cfg = [];
cfg.indir   = indir;
cfg.cln_fld = cln_fld;
cfg.cln_dir = cln_dir;
cfg.end_dir = end_dir;
cfg.tsk     = tsk;
[inpath,sem_dat,trl] = mmil_find_files(cfg);

if numel(tsk) ~= numel(trl.data_name);
    for iT = 1:numel(tot_tsk)
        tsk(string_find(trl.data_name,{tot_tsk{iT}})) = repmat({tot_tsk{iT}},1,numel(string_find(trl.data_name,{tot_tsk{iT}})));
    end
end

%% Load Visual & Auditory data
epc_tme      = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

if ~any(~cellfun(@isempty,strfind(end_dir,'mat')));
    cfg            = [];
    cfg.specific   = {'dataset';1:numel(inpath)};
    cfg.data_new   = 'yes';
    if any(~cellfun(@isempty,strfind(end_dir,'set')) | ~cellfun(@isempty,strfind(end_dir,'edf'))); cfg.continuous = 'yes'; else cfg.continuous = 'no'; end
    cfg.dataset    = inpath;
    sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);
else
    for iL = 1:numel(sem_dat.data_name)
        ttt = load(inpath{iL});
        sem_dat.(sem_dat.data_name{iL}) = ttt.epoch_data;
        clear ttt
    end
end

%% Fix Events & Combine Data
if any(~cellfun(@isempty,strfind(end_dir,'eeg')));
    
    for iTS = 1:numel(tsk)
        cfg        = [];
        cfg.task   = tsk{iTS};
        cfg.inpath = inpath{iTS};
        skp{iTS}   = mmil_check_eeg(cfg);
    end
    
    if strcmpi(subj,'NY017_SA_SZ'); cfg = []; cfg.trials = 1:538; sem_dat.(sem_dat.data_name{1}) = ft_preprocessing(cfg,sem_dat.(sem_dat.data_name{1})); end
    if strcmpi(subj,'NY019_SA_SZ'); cfg = []; cfg.trials = [1:338 340:400]; sem_dat.(sem_dat.data_name{1}) = ft_preprocessing(cfg,sem_dat.(sem_dat.data_name{1})); end
    
    for iSK = 1:numel(skp); sem_dat.(sem_dat.data_name{iSK}).cfg.alt_eve.skp = skp{iSK}; end
    
    cfg          = [];
    cfg.clr_fld  = fcfg.clr_fld;
    cfg.end_dir  = end_dir;
    cfg.task     = [tsk(:)];
    cfg.skp      = skp;
    cfg.rep_fix  = 1;
    cfg.specific = {'task' 'skp'; 1:numel(tsk) 1:numel(skp)};
    sem_dat      = ft_func(@ft_sasz_addevent,cfg,sem_dat);
    
elseif any(~cellfun(@isempty,strfind(end_dir,'edf')));
    
    trialfun         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun');
    epc_tme          = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
    
    if ~exist([fcfg.clr_fld '/trialfun_output/' subj '.mat'])
        
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
    end
    
    cfg = [];
    cfg.specific  = {'trl';1:numel(inpath)};
    cfg.trl       = trl_hld;
    sem_dat = ft_func(@ft_redefinetrial,cfg,sem_dat);
    
    % Combine as necessary
    for iFL = 1:numel(tot_tsk)
        if sum(strcmpi(tsk,tot_tsk{iFL})) > 1
            
            if any(~cellfun(@isempty,strfind(sem_dat.data_name,'lin1')))
                
                for iCL = 1:2
                    
                    idn = find(strcmpi(tsk,tot_tsk{iFL}));
                    ttt.cl1 = find(~cellfun(@isempty,strfind(sem_dat.data_name,'lin1')));
                    ttt.cl2 = find(~cellfun(@isempty,strfind(sem_dat.data_name,'lin2')));
                    cln_idn = intersect(ttt.(['cl' num2str(iCL)]),idn);
                    
                    cfg           = [];
                    cfg.data_name = [repmat(cln_idn(1),1,numel(cln_idn)-1) ; cln_idn(2:end)];
                    cfg.data_new  = 'yes';
                    cfg.methapp   = 'trials';
                    sem_dat = ft_func(@ft_appenddata,cfg,sem_dat);
                    
                    cfg = [];
                    cfg.rmfield   = 'yes';
                    cfg.data_name =  cln_idn(2:end);
                    sem_dat = ft_func([],cfg,sem_dat);
                   
                    tsk(cln_idn(2:end)) = [];
                    
                end
                
            else
                
                idn = find(strcmpi(tsk,tot_tsk{iFL}));
                
                cfg = [];
                cfg.data_name = idn;
                cfg.channel   = ['all',strcat('-',intersect(sem_dat.(sem_dat.data_name{idn(1)}).label,sem_dat.(sem_dat.data_name{idn(2)}).label))'];
                sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
                
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
    
elseif any(~cellfun(@isempty,strfind(end_dir,'set')));
    
    if ~exist([fcfg.clr_fld '/trialfun_output/' subj '.mat'])
        
        trialfun         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun');
        epc_tme          = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
        
        for iEP = 1:numel(inpath)
            cfg = [];
            cfg.indir    = inpath{iEP};
            cfg.tsk      = tsk{iEP};
            cfg.pre      = 0.95;
            cfg.pos      = 1.95;
            cfg.smp_rte  = sem_dat.(sem_dat.data_name{iEP}).fsample;
            trl_hld{iEP} = mmil_epoch_dat_set(cfg);
        end
        
        for i = 1:numel(tsk); if strcmpi(tsk{i},'SA'); trl_hld{i}(:,4) = trl_hld{i}(:,4)+10; end; end
        
        save([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
    else
        load([fcfg.clr_fld '/trialfun_output/' subj '.mat'],'trl_hld');
    end
    
    cfg = [];
    cfg.specific  = {'trl'; 1:numel(inpath)};
    cfg.trl       = trl_hld;
    sem_dat = ft_func(@ft_redefinetrial,cfg,sem_dat);
    
    % Combine as necessary
    for iFL = 1:numel(tot_tsk)
        if sum(strcmpi(tsk,tot_tsk{iFL})) > 1
            
            if any(~cellfun(@isempty,strfind(sem_dat.data_name,'lin1')))
                
                for iCL = 1:2
                    
                    idn = find(strcmpi(tsk,tot_tsk{iFL}));
                    ttt.cl1 = find(~cellfun(@isempty,strfind(sem_dat.data_name,'lin1')));
                    ttt.cl2 = find(~cellfun(@isempty,strfind(sem_dat.data_name,'lin2')));
                    cln_idn = intersect(ttt.(['cl' num2str(iCL)]),idn);
                    
                    cfg           = [];
                    cfg.data_name = [repmat(cln_idn(1),1,numel(cln_idn)-1) ; cln_idn(2:end)];
                    cfg.data_new  = 'yes';
                    cfg.methapp   = 'trials';
                    sem_dat = ft_func(@ft_appenddata,cfg,sem_dat);
                    
                    cfg = [];
                    cfg.rmfield   = 'yes';
                    cfg.data_name =  cln_idn(2:end);
                    sem_dat = ft_func([],cfg,sem_dat);
                    
                    tsk(cln_idn(2:end)) = [];
                    
                end
                
            else
                
                idn = find(strcmpi(tsk,tot_tsk{iFL}));
                
                cfg = [];
                cfg.data_name = idn;
                cfg.channel   = ['all',strcat('-',intersect(sem_dat.(sem_dat.data_name{idn(1)}).label,sem_dat.(sem_dat.data_name{idn(2)}).label))'];
                sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
                
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
        
elseif any(~cellfun(@isempty,strfind(end_dir,'mat')));
    
    if strcmpi(subj,'NY226_SA_SZ')
        sem_dat = NY226_SA_SZ_fix(sem_dat);
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
% identified channels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.rmv_chn = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'rmv_chn'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';

sem_dat = mmil_remove_channels(cfg,sem_dat);

%% Examine Data for Noise - Initial
beg_plt_spc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'beg_plt_spc'); beg_plt_spc = beg_plt_spc{1};

if beg_plt_spc %&& ~exist([outpath '/' 'spectrum' '/' subj '/' 'Initial' '/' sem_dat.data_name{1} '_1.png'],'file');
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' 'Initial' ];
    cfg.prefix = sem_dat.data_name;
    cfg.specific  = {'prefix'; 1:numel(sem_dat.data_name)};
    ft_func(@ft_plot_spectrum,cfg,sem_dat);
end

%% Remove noise through removing common noise
chn_nse_grp = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
chn_nse_grp_nme = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nme');
nse_val         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'nse_val'); nse_val = nse_val{1};

% Make temporary ECOG/DEPTH Split
for iD = 1:numel(sem_dat.data_name); sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label = sem_dat.(sem_dat.data_name{iD}).label; end
scfg = [];
scfg.dat_nme = sem_dat.data_name;
scfg.clr_fld = fcfg.clr_fld;
scfg.sbj_nme = fcfg.sbj_nme;
scfg.alt_lab = 'label';
scfg.specific = {'dat_nme' ; 1:numel(sem_dat.data_name)};
ft_func(@mmil_create_depth2,scfg,sem_dat);

if nse_val
    
    if strcmpi(chn_nse_grp{1},'split'); chn_nse_grp = repmat(chn_nse_grp,1,numel(sem_dat.data_name)); end
    
    for iEJ = 1:numel(chn_nse_grp_nme);
        if strcmpi(chn_nse_grp{1},'split')
            for iD = 1:numel(sem_dat.data_name); pre_fix{iD}{1} = strcat(sem_dat.data_name{iD},'ecog'); pre_fix{iD}{2} = strcat(sem_dat.data_name{iD},'noisy_ecog'); pre_fix{iD}{3} = strcat(sem_dat.data_name{iD},'depth'); pre_fix{iD}{4} = strcat(sem_dat.data_name{iD},'noisy_depth'); end;
        else
            for iCN = 1:numel(chn_nse_grp_nme{iEJ});
                pre_fix{iEJ}{iCN} = strcat(sem_dat.data_name{iEJ},'_',chn_nse_grp_nme{iEJ}{iCN});
            end;
        end
    end
    
    cfg             = [];
    cfg.sbj_nme     = fcfg.sbj_nme;
    cfg.clr_fld     = fcfg.clr_fld;
    cfg.pre_fix     = pre_fix;
    cfg.dat_nme     = sem_dat.data_name;
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix     = pre_fix;
    cfg.specific    = {'pre_fix' 'chn_nse_grp' 'dat_nme' ; 1:numel(sem_dat.data_name) 1:numel(sem_dat.data_name) 1:numel(sem_dat.data_name)};
    cfg.out_dir     = [fcfg.out_pth '/' 'cmn_nse_rmv'];
    cfg.rmv_chn     = nse_val;
    sem_dat = ft_func(@ft_remove_common_noise2,cfg,sem_dat);   
    
    if nse_val == 1
        cfg.clr_fld = fcfg.clr_fld;
        cfg.sbj_nme = fcfg.sbj_nme;
        cfg.sub_fld_ind = num2cell(1:numel(sem_dat.data_name));
        cfg.specific    = {'pre_fix'  'sub_fld_ind' ; 1:numel(sem_dat.data_name) 1:numel(sem_dat.data_name)};
        sem_dat = ft_func(@mmil_update_cmn_nse2,cfg,sem_dat);
    end
    
end

%% Remove identified noise problems
bse_frq     = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_frq'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';

if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    sem_dat      = ft_func(@ft_preprocessing,cfg,sem_dat);
end

%% Examine Data for Noise
end_plt_spc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_plt_spc'); end_plt_spc = end_plt_spc{1};

if end_plt_spc == 1 %&& ~exist([outpath '/' 'spectrum' '/' subj '/' 'After' '/' sem_dat.data_name{1} '_1.png'],'file');
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [outpath '/' 'spectrum' '/' subj '/' 'After'];
    cfg.prefix = sem_dat.data_name;
    cfg.specific  = {'prefix'; 1:numel(sem_dat.data_name)};
    ft_func(@ft_plot_spectrum,cfg,sem_dat);
end

%% Combine if 2 clinical systems
if any(~cellfun(@isempty,strfind(sem_dat.data_name,'lin2')))
    
    for iFL = 1:numel(tot_tsk)
        
        idn = find(strcmpi(tsk,tot_tsk{iFL}));
        
        cfg = [];
        cfg.data_name = idn;
        cfg.channel   = ['all',strcat('-',intersect(sem_dat.(sem_dat.data_name{idn(1)}).label,sem_dat.(sem_dat.data_name{idn(2)}).label))'];
        sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
        
        cfg           = [];
        cfg.data_name = [repmat(idn(1),1,numel(idn)-1) ; idn(2):idn(end)];
        cfg.data_new  = 'yes';
        cfg.methapp   = 'trials';
        sem_dat = ft_func(@ft_appenddata,cfg,sem_dat);
        
        cfg = [];
        cfg.rmfield   = 'yes';
        cfg.data_name = idn(2):idn(end);
        sem_dat = ft_func([],cfg,sem_dat);
        
        tsk(idn(2:end)) = []; 
        
        sem_dat.(sem_dat.data_name{iFL}).cfg.alt_eve.trialinfo = sem_dat.(sem_dat.data_name{iFL}).trialinfo;
        
        cfg = [];
        cfg.tsk = tsk{iFL};
            
    end
end

for iSK = 1:numel(tsk)
    cfg = [];
    cfg.tsk = tsk{iSK};
    sem_dat.(sem_dat.data_name{iSK}).cfg.alt_eve.skp = mmil_missing_trials(cfg,sem_dat.(sem_dat.data_name{iSK}));
end

%% Fix Labels
cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.clr_fld = fcfg.clr_fld;
cfg.chn_loc = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'electrode_location');

sem_dat = mmil_fix_labels(cfg,sem_dat);

%% Filter Data for LFP
cfg            = [];
cfg.data_new   = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 15;
cfg.new_suffix = 'lfp';
sem_dat        = ft_func(@ft_preprocessing,cfg,sem_dat);

%% Filter Data for HGP
cfg=[];
cfg.data_name  = 1:numel(tsk);
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

cfg             = [];
cfg.data_name   = numel(tsk)*2+1:numel(tsk)*3;
HGP_smooth_msec = 0.025;
w               = round(HGP_smooth_msec* sem_dat.(sem_dat.data_name{1}).fsample);
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,round(0.016* sem_dat.(sem_dat.data_name{1}).fsample));
cfg.window      = gauss_y/sum(gauss_y);
sem_dat         = ft_func(@ft_window_data,cfg,sem_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1:numel(tsk);
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
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sem_dat);

end