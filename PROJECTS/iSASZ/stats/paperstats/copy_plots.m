

%% Move files around
if mve_plt == 1
    
    outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
    
    for iMV = 1:numel(sem_dat.data_name)
        
        if has_vis && has_aud
            cfg = [];
            cfg.plt_dim = [1 numel(sem_dat.(sem_dat.data_name{iMV}).label)];
            cfg.cmb     = 1;
            cfg.cmp_ind = [1 3 7 9];
            cfg.sig_fle = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}];
            [grp_chn,pre_fix,grp_stt,grp_stt_col] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iMV}));
        elseif has_vis
            cfg = [];
            cfg.plt_dim = [1 numel(sem_dat.(sem_dat.data_name{iMV}).label)];
            cfg.cmb     = 1;
            cfg.cmp_ind = [1 3];
            cfg.sig_fle = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}];
            [grp_chn,pre_fix,grp_stt,grp_stt_col] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iMV}));
        elseif has_aud
            cfg = [];
            cfg.plt_dim = [1 numel(sem_dat.(sem_dat.data_name{iMV}).label)];
            cfg.cmb     = 1;
            cfg.cmp_ind = [1 3];
            cfg.sig_fle = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/sigchn' '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iMV}];
            [grp_chn,pre_fix,grp_stt,grp_stt_col] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iMV}));
        end
        
        % Copyfiles over
        fld = [outpath '/' 'initial_plot' '/' 'ThreeStat' '/' subj '/' sem_dat.data_name{iMV}(end-2:end)];
        
        for iCP = 1:numel(pre_fix)
            if exist([fld '/' pre_fix{iCP}],'dir'); rmdir([fld '/' pre_fix{iCP}],'s'); mkdir([fld '/' pre_fix{iCP}]); else mkdir([fld '/' pre_fix{iCP}]); end
            for iFL = 1:numel(grp_chn{iCP});
                copyfile([fld '/' subj '__average_chan_' sem_dat.(sem_dat.data_name{iMV}).cfg.alt_lab.label{grp_chn{iCP}(iFL)} '.png'], ...
                    [fld '/' pre_fix{iCP} '/' subj '__average_chan_' sem_dat.(sem_dat.data_name{iMV}).cfg.alt_lab.label{grp_chn{iCP}(iFL)} '.png']);
            end
        end
        
    end
end