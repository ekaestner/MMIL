function iSASZ_initial_statsv2(fcfg)

    sbj  = fcfg.sbj_nme;
    
    fprintf([fcfg.sbj_nme 'Starting stats work on %s \n'],sbj)
    
    infile = [fcfg.out_pth sbj '_overall_data.mat'];
    outpath = fcfg.out_pth;
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [outpath '/' sbj '_overall_data.mat'];
    sem_dat  = ft_func([],cfg);
    
%     for iD = 1:numel(sem_dat.data_name)
%         sem_dat.(sem_dat.data_name{iD}).cfg = rmfield(sem_dat.(sem_dat.data_name{iD}).cfg,'alt_stt');
%     end

%     if strcmpi(sbj,'NY226_SA_SZ')
%         for iD = 1:numel(sem_dat.data_name)
%             tmp_nme = sem_dat.data_name{iD};
%             tmp_fld = sem_dat.(tmp_nme);
%             sem_dat = rmfield(sem_dat,tmp_nme);
%             tmp_nme = tmp_nme([1:6 end-3:end]);
%             
%             sem_dat.data_name{iD} = tmp_nme;
%             sem_dat.(sem_dat.data_name{iD}) = tmp_fld;
%         end
%     end
    
    tic
    cfg = [];
    cfg.stt_fnc  = {'sasz_ovr_all_stt' 'sasz_rep_stt' 'sasz_sem_stt'};
    cfg.loc      = 'local';
    cfg.fld_nme  = sem_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    sem_dat      = ft_func(@mmil_cloud_stat,cfg,sem_dat);
    toc
    
    tic
    cfg     = [];
    cfg.stt_fnc  = {'sasz_ovr_all_stt'};
    cfg.fld_nme  = sem_dat.data_name;
    cfg.loc = 'local';
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    sem_dat = ft_func(@mmil_cloud_stat,cfg,sem_dat);
    toc
    
    %% Save Data & Stats
    cfg = [];
    cfg.str_nme  = 'sem_dat';
    cfg.save     = 'yes';
    cfg.sve_app  = 'app_all';
    cfg.filename =[outpath '/' sbj '_overall_data.mat'];
    ft_func([],cfg,sem_dat);
    
end
