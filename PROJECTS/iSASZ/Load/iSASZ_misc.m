%% Todo
% ignore?
if strcmpi(subj,'NY017_SA_SZ'); cfg = []; cfg.trials = 1:538; sem_dat.(sem_dat.data_name{1}) = ft_preprocessing(cfg,sem_dat.(sem_dat.data_name{1})); end
if strcmpi(subj,'NY019_SA_SZ'); cfg = []; cfg.trials = [1:338 340:400]; sem_dat.(sem_dat.data_name{1}) = ft_preprocessing(cfg,sem_dat.(sem_dat.data_name{1})); end
   
    



%% Old 
% Misc
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

% Load
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