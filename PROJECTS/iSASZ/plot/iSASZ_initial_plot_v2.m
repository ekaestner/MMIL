function iSASZ_initial_plot_v2(fcfg)

    ovr_wrt = fcfg.ovr_wrt;

    sbj  = fcfg.sbj_nme;
    
    inpath = fcfg.in_pth;
    
    outpath = fcfg.out_pth;
    
    cln_dir = fcfg.clr_fld;
    
    try chn_loc = mmil_readtext([fcfg.clr_fld '/' sbj '_loc']);
    catch
        has_loc = 0;
    end
    if ~isvar('has_loc'); has_loc = 1; end
    
    cfg = [];
    cfg.load    = 'yes';
    cfg.file = [inpath '/' sbj '_overall_data.mat'];
    sem_dat     = ft_func([],cfg);
    
    if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo > 9); has_aud = 1; else has_aud = 0; end
    if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo < 9); has_vis = 1; else has_vis = 0; end
    
    fprintf('Starting work on %s \n',sbj)
    
    %% Plotting Events
    if has_vis && has_aud
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events  = {[2 4] [1 3] [5 6 7 8] [12 14] [11 13] [15 16 17 18]};
        cfg.new_events  = {131 132 133 134 135 136};
        cfg.crt_alt_eve = 'plt_all_eve';
        sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
    elseif has_vis
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events  = {[2 4] [1 3] [5 6 7 8]};
        cfg.new_events  = {131 132 133};
        cfg.crt_alt_eve = 'plt_all_eve';
        sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
    elseif has_aud
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events  = {[12 14] [11 13] [15 16 17 18]};
        cfg.new_events  = {134 135 136};
        cfg.crt_alt_eve = 'plt_all_eve';
        sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
    end
    
    %% Fix stat labels
    for iD = 1:numel(sem_dat.data_name)
        
%         sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.plt_lab = sem_dat.(sem_dat.data_name{iD}).label;
        sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.plt_lab = sem_dat.(sem_dat.data_name{iD}).cfg.alt_lab.label;
    
    end
    
    %% Plot Side-by-Side 3 Stat Plots
    if has_vis && has_aud
        
        cfg = [];
        
        cfg.type      = 'chan';
        cfg.dat       = {sem_dat.(sem_dat.data_name{2}) sem_dat.(sem_dat.data_name{1})};
        cfg.dat_loc   = [1 1 2 2 ; ...
            1 1 2 2 ];
        cfg.alt_eve   = { 'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt'; ...
            'rep_all_eve' 'rep_all_eve' 'rep_all_eve' 'rep_all_eve'};
        cfg.eve       = { [101]     [102]     [101]     [102]; ...
            [111 112] [113 114] [111 112] [113 114]};
        cfg.plt_lbl   = { 'Visual HGP' 'Auditory  HGP' 'Visual LFP' 'Auditory  LFP' ; ...
            ''           ''              ''           ''};
        cfg.ttl_num   = [1 1 1 1 ; ...
            0 0 0 0 ];
        cfg.plt_dim   = [2 4];
        
        cfg.lgd       = 0;
        cfg.std_err   = 1;
        cfg.lnstyle.col_ord = {rgb('reddish grey') rgb('bluish grey') rgb('green') rgb('red') rgb('violet') rgb('blue')};
        cfg.cnd_nme         = {'V Ovr' 'A Ovr' 'V Nov' 'V Rep'  'A Nov' 'A Rep'};
        cfg.leg_pos         = [1 5 2 3 6 7];
        cfg.alt_lbl         = {'plt_lab' 'plt_lab' 'plt_lab' 'plt_lab' ; ...
            'plt_lab' 'plt_lab' 'plt_lab' 'plt_lab' };
        
        cfg.stt_dat = { {'vis_ovr_all_stt'} {'aud_ovr_all_stt'} {'vis_ovr_all_stt'} {'aud_ovr_all_stt'} ; ...
            {'vis_dif_rep_stt_msk' 'vis_dif_rep_stt'} {'aud_dif_rep_stt_msk' 'aud_dif_rep_stt'} {'vis_dif_rep_stt_msk' 'vis_dif_rep_stt'} {'aud_dif_rep_stt_msk' 'aud_dif_rep_stt'} };
         cfg.stt_col = { {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}       ; ...
            {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]}};
       
        cfg.x_lim       = [-0.3 1.2];
        cfg.v_lne       = {0 0 0 0 ; ...
            0 0 0 0 };
        cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} ; ...
            {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} };
        cfg.v_lne       = { [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] ; ...
            [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] };
        cfg.v_lne_wdt   = { [3 1 1] [3 1 1] [3 1 1] [3 1 1] ; ...
            [3 1 1] [3 1 1] [3 1 1] [3 1 1] };
        cfg.axe_fnt_sze = [ 20 20 20 20; ...
            20 20 20 20];
        cfg.axe_lne_sze = [ 4 4 4 4 ; ...
            4 4 4 4 ];
        cfg.ttl_lne_sze = [50 50 50 50 ; ...
            50 50 50 50 ];
        cfg.y_axe       = [1 0 1 0; ...
            1 0 1 0 ];
        cfg.y_lnk = [1 1 2 2 ; ...
            1 1 2 2 ];
        cfg.x_axe = [0 0 0 0 ; ...
            1 1 1 1 ];
        
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = [outpath '/' 'paper_plot' '/' 'ThreeStat' '/' sbj '/' 'ovr'];
        cfg.prefix     = [sbj '_'];
        
        if strcmpi(sbj,'NY008_SA_SZ')
            cfg.y_lnk = [1 2 3 4 ; ...
                1 2 3 4 ];
            cfg.y_axe = [1 1 1 1; ...
                1 1 1 1];
        end
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
    elseif has_vis
        
        cfg = [];
        
        cfg.type      = 'chan';
        cfg.dat       = {sem_dat.(sem_dat.data_name{2}) sem_dat.(sem_dat.data_name{1})};
        cfg.dat_loc   = [1 1 2 2 ; ...
            1 1 2 2 ];
        cfg.alt_eve   = { 'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt'; ...
            'rep_all_eve' 'rep_all_eve' 'rep_all_eve' 'rep_all_eve'};
        cfg.eve       = { [101]     [102]     [101]     [102]; ...
            [111 112] [113 114] [111 112] [113 114]};
        cfg.plt_lbl   = { 'Visual HGP' 'Auditory  HGP' 'Visual LFP' 'Auditory  LFP' ; ...
            ''           ''              ''           ''};
        cfg.ttl_num   = [1 1 1 1 ; ...
            0 0 0 0 ];
        cfg.plt_dim   = [2 4];
        
        cfg.lgd       = 0;
        cfg.std_err   = 1;
        cfg.lnstyle.col_ord = {rgb('reddish grey') rgb('bluish grey') rgb('green') rgb('red') rgb('violet') rgb('blue')};
        cfg.cnd_nme         = {'V Ovr' 'A Ovr' 'V Nov' 'V Rep'  'A Nov' 'A Rep'};
        cfg.leg_pos         = [1 5 2 3 6 7];
        cfg.alt_lbl         = {'plt_lab' 'plt_lab' 'plt_lab' 'plt_lab' ; ...
            'plt_lab' 'plt_lab' 'plt_lab' 'plt_lab' };
        
        cfg.stt_dat = { {'vis_ovr_all_stt'} {} {'vis_ovr_all_stt'} {} ; ...
            {'vis_dif_rep_stt' 'vis_dif_rep_stt'} {} {'vis_dif_rep_stt' 'vis_dif_rep_stt'} {} };
        cfg.stt_col = { {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}       ; ...
            {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]}};
        
        cfg.x_lim       = [-0.3 1.2];
        cfg.v_lne       = {0 0 0 0 ; ...
            0 0 0 0 };
        cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} ; ...
            {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} };
        cfg.v_lne       = { [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] ; ...
            [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] };
        cfg.v_lne_wdt   = { [3 1 1] [3 1 1] [3 1 1] [3 1 1] ; ...
            [3 1 1] [3 1 1] [3 1 1] [3 1 1] };
        cfg.axe_fnt_sze = [ 20 20 20 20; ...
            20 20 20 20];
        cfg.axe_lne_sze = [ 4 4 4 4 ; ...
            4 4 4 4 ];
        cfg.ttl_lne_sze = [50 50 50 50 ; ...
            50 50 50 50 ];
        cfg.y_axe       = [1 0 1 0; ...
            1 0 1 0 ];
        cfg.y_lnk = [1 1 2 2 ; ...
            1 1 2 2 ];
        cfg.x_axe = [0 0 0 0 ; ...
            1 1 1 1 ];
        
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = [outpath '/' 'paper_plot' '/' 'ThreeStat' '/' sbj '/' 'ovr'];
        cfg.prefix     = [sbj '_'];
        
        if strcmpi(sbj,'NY008_SA_SZ')
            cfg.y_lnk = [1 2 3 4 ; ...
                1 2 3 4 ];
            cfg.y_axe = [1 1 1 1; ...
                1 1 1 1];
        end
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
    elseif has_aud
        
        cfg = [];
        
        cfg.type      = 'chan';
        cfg.dat       = {sem_dat.(sem_dat.data_name{2}) sem_dat.(sem_dat.data_name{1})};
        cfg.dat_loc   = [1 1 2 2 ; ...
            1 1 2 2 ];
        cfg.alt_eve   = { 'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt' 'ovr_all_evt'; ...
            'rep_all_eve' 'rep_all_eve' 'rep_all_eve' 'rep_all_eve'};
        cfg.eve       = { [101]     [102]     [101]     [102]; ...
            [111 112] [113 114] [111 112] [113 114]};
        cfg.plt_lbl   = { 'Visual HGP' 'Auditory  HGP' 'Visual LFP' 'Auditory  LFP' ; ...
            ''           ''              ''           ''};
        cfg.ttl_num   = [1 1 1 1 ; ...
            0 0 0 0 ];
        cfg.plt_dim   = [2 4];
        
        cfg.lgd       = 0;
        cfg.std_err   = 1;
        cfg.lnstyle.col_ord = {rgb('reddish grey') rgb('bluish grey') rgb('green') rgb('red') rgb('violet') rgb('blue')};
        cfg.cnd_nme         = {'V Ovr' 'A Ovr' 'V Nov' 'V Rep'  'A Nov' 'A Rep'};
        cfg.leg_pos         = [1 5 2 3 6 7];
        cfg.alt_lbl         = {'plt_lab' 'plt_lab' 'plt_lab' 'plt_lab' ; ...
            'plt_lab' 'plt_lab' 'plt_lab' 'plt_lab' };
        
        cfg.stt_dat = { {} {'aud_ovr_all_stt'} {} {'aud_ovr_all_stt'} ; ...
            {} {'aud_dif_rep_stt_msk' 'aud_dif_rep_stt'} {} {'aud_dif_rep_stt_msk' 'aud_dif_rep_stt'} };
         cfg.stt_col = { {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}        {[0.85 0.85 0.85]}       ; ...
            {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]} {ft_stt_col(rgb('dark yellow')) [0.9 0.9 0.9]}};
       
        cfg.x_lim       = [-0.3 1.2];
        cfg.v_lne       = {0 0 0 0 ; ...
            0 0 0 0 };
        cfg.v_lne_col   = { {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} ; ...
            {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} {rgb('red') rgb('black') rgb('black')} {rgb('blue') rgb('black') rgb('black')} };
        cfg.v_lne       = { [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] ; ...
            [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] [0.000 0.300 0.600] };
        cfg.v_lne_wdt   = { [3 1 1] [3 1 1] [3 1 1] [3 1 1] ; ...
            [3 1 1] [3 1 1] [3 1 1] [3 1 1] };
        cfg.axe_fnt_sze = [ 20 20 20 20; ...
            20 20 20 20];
        cfg.axe_lne_sze = [ 4 4 4 4 ; ...
            4 4 4 4 ];
        cfg.ttl_lne_sze = [50 50 50 50 ; ...
            50 50 50 50 ];
        cfg.y_axe       = [1 0 1 0; ...
            1 0 1 0 ];
        cfg.y_lnk = [1 1 2 2 ; ...
            1 1 2 2 ];
        cfg.x_axe = [0 0 0 0 ; ...
            1 1 1 1 ];
        
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = [outpath '/' 'paper_plot' '/' 'ThreeStat' '/' sbj '/' 'ovr'];
        cfg.prefix     = [sbj '_'];
        
        if strcmpi(sbj,'NY008_SA_SZ')
            cfg.y_lnk = [1 2 3 4 ; ...
                1 2 3 4 ];
            cfg.y_axe = [1 1 1 1; ...
                1 1 1 1];
        end
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end
