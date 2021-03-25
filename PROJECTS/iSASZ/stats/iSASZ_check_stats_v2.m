function iSASZ_check_stats_v2(fcfg)
    
    subj  = fcfg.sbj_nme;
    
    fprintf('Starting work on %s \n',subj)
    
    infile = [fcfg.fle_out_pth subj '_overall_data.mat'];
    outpath = fcfg.clr_out_pth;
    fle_outpath = fcfg.fle_out_pth;
    
    cfg = [];
    cfg.load = 'yes';
    cfg.file = infile;
    sem_dat  = ft_func([],cfg);
    
    has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
    has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);
    
    %% Stat Mask
    if has_aud
        cfg     = [];
        cfg.stt     = {'aud_dif_rep_stt' 'aud_dif_sem_stt'};
        cfg.stt_msk = {'aud_ovr_all_stt' 'aud_ovr_all_stt'};
        sem_dat = ft_func(@ft_mask_stats,cfg,sem_dat);
    end
    if has_vis
        cfg     = [];
        cfg.stt     = {'vis_dif_rep_stt' 'vis_dif_sem_stt'};
        cfg.stt_msk = {'vis_ovr_all_stt' 'vis_ovr_all_stt'};
        sem_dat = ft_func(@ft_mask_stats,cfg,sem_dat);
    end
    
    %% Fix Stat Channels
    for iD = 1:numel(sem_dat.data_name)
        stt_fld_nme = fieldnames(sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt);
        for iST = 1:numel(stt_fld_nme)
            if numel(intersect(sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt.(stt_fld_nme{iST}).label,sem_dat.(sem_dat.data_name{iD}).label)) == numel(sem_dat.(sem_dat.data_name{iD}).label)
                sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt.(stt_fld_nme{iST}).old_lab = sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt.(stt_fld_nme{iST}).label;
                [~,~,stt_nme_ord] = intersect(sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt.(stt_fld_nme{iST}).label,sem_dat.(sem_dat.data_name{iD}).label);
                sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt.(stt_fld_nme{iST}).label = sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label(stt_nme_ord);
            else
                sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt.(stt_fld_nme{iST}).label = sem_dat.(sem_dat.data_name{iD}).cfg.alt_stt.(stt_fld_nme{iST}).old_lab;
            end
        end
    end
    
    
    %% Choose Channels based on hypotheses
    if fcfg.chn_chs == 1;
        
        if has_vis && has_aud
            
            cfg = [];
            cfg.ovr_wrt = fcfg.ovr_wrt;
            cfg.alt_stt = {'vis_ovr_all_stt' 'vis_ovr_all_stt' 'vis_ovr_all_stt' 'vis_dif_rep_stt_msk' 'vis_dif_rep_stt_msk' 'vis_dif_rep_stt_msk' ...
                           'aud_ovr_all_stt' 'aud_ovr_all_stt' 'aud_ovr_all_stt' 'aud_dif_rep_stt_msk' 'aud_dif_rep_stt_msk' 'aud_dif_rep_stt_msk'};
            cfg.alt_stt_col = {[0.7 0.7 0.7] [0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ...
                               [0.7 0.7 0.7] [0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow'))};
            cfg.cmp_stt = [1 2 3 4 5  6 ...
                           7 8 9 10 11 12];
            cfg.cmp_trl = {'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'rep_all_eve' ...
                           'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'rep_all_eve'};
            cfg.cmp     = {'101!999' '101!999' '111!112' '111!112' '121!122' '121!122'...
                           '102!999' '102!999' '113!114' '113!114' '123!124' '123!124'};
            cfg.cmp_nme = {'vis_ovr_all_pre' 'vis_ovr_all_lex' 'vis_ovr_all_pos' 'vis_dif_rep_pre' 'vis_dif_rep_lex' 'vis_dif_rep_pos' ...
                           'aud_ovr_all_pre' 'aud_ovr_all_lex' 'aud_ovr_all_pos' 'aud_dif_rep_pre' 'aud_dif_rep_lex' 'aud_dif_rep_pos'};
            cfg.tme_win = {[0.100 0.300] [0.300 0.600] [0.600 1.000] [0.100 0.300] [0.300 0.600] [0.600 1.000] ...
                           [0.100 0.300] [0.300 0.600] [0.600 1.000] [0.100 0.300] [0.300 0.600] [0.600 1.000]};
            cfg.clr_fld = outpath; %'/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn';
            cfg.sbj_nme = subj;
            cfg.typ      = sem_dat.data_name;
            cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
            sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
            
        elseif has_vis
            
            cfg = [];
            cfg.ovr_wrt = fcfg.ovr_wrt;
            cfg.alt_stt = {'vis_ovr_all_stt' 'vis_ovr_all_stt' 'vis_ovr_all_stt' 'vis_dif_rep_stt_msk' 'vis_dif_rep_stt_msk' 'vis_dif_rep_stt_msk'};
            cfg.alt_stt_col = {[0.7 0.7 0.7] [0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow'))};
            cfg.cmp_stt = [1 2 3 4 5 6];
            cfg.cmp_trl = {'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'rep_all_eve'};
            cfg.cmp     = {'101!999' '101!999' '101!999' '111!112' '111!112' '111!112'};
            cfg.cmp_nme = {'vis_ovr_all_pre' 'vis_ovr_all_lex' 'vis_ovr_all_pos' 'vis_dif_rep_pre' 'vis_dif_rep_lex' 'vis_dif_rep_pos'};
            cfg.tme_win = {[0.100 0.300] [0.300 0.600] [0.600 1.000] [0.100 0.300] [0.300 0.600] [0.600 1.000]};
            cfg.clr_fld = outpath; %'/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn';
            cfg.sbj_nme = subj;
            cfg.typ      = sem_dat.data_name;
            cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
            sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
            
        elseif has_aud
            
            cfg = [];
            cfg.ovr_wrt = fcfg.ovr_wrt;
            cfg.alt_stt = {'aud_ovr_all_stt' 'aud_ovr_all_stt' 'aud_ovr_all_stt' 'aud_dif_rep_stt_msk' 'aud_dif_rep_stt_msk' 'aud_dif_rep_stt_msk'};
            cfg.alt_stt_col = {[0.7 0.7 0.7] [0.7 0.7 0.7] [0.7 0.7 0.7] ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow')) ft_stt_col(rgb('dark yellow'))};
            cfg.cmp_stt = [1 2 3 4 5 6];
            cfg.cmp_trl = {'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt' 'rep_all_eve' 'rep_all_eve' 'rep_all_eve'};
            cfg.cmp     = {'102!999' '102!999' '102!999' '113!114' '113!114' '113!114'};
            cfg.cmp_nme = {'aud_ovr_all_pre' 'aud_ovr_all_lex' 'aud_dif_rep_pos' 'aud_dif_rep_pre' 'aud_dif_sem_lex' 'aud_dif_sem_pos'};
            cfg.tme_win = {[0.100 0.300] [0.300 0.600] [0.600 1.000] [0.100 0.300] [0.300 0.600] [0.600 1.000]};
            cfg.clr_fld = outpath; %'/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn';
            cfg.sbj_nme = subj;
            cfg.typ      = sem_dat.data_name;
            cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
            sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
            
        end
    end
    
    %% Save Data & Stats
    cfg = [];
    cfg.str_nme  = 'sem_dat';
    cfg.save     = 'yes';
    cfg.sve_app  = 'app_all';
    cfg.filename =[fle_outpath '/' subj '_overall_data.mat'];
    ft_func([],cfg,sem_dat);
    
end