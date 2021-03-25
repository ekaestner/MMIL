%% Load each subjects data & sig_chn hgp
clear; clc;

clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
sbj_nme = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/' 'subjects']);

indir       = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';
abs_out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

for iS = [1:3 5:16]
    
    sbj = sbj_nme{iS};
    
    % Load data & sig_chn
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_hgp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    new_sig = mmil_readtext([clr_fld '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],['\t']);
    
    has_vis = any(strcmpi('visual_early',new_sig(1,:)));
    has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
    
    % Remove non-sig-chn
    if has_aud && has_vis
        cfg = [];
        cfg.channel = sem_dat.(sem_dat.data_name{1}).label(find(cell2mat(new_sig(2:end,2)) | cell2mat(new_sig(2:end,4))));
        sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
    elseif has_aud || has_vis
        cfg = [];
        cfg.channel = sem_dat.(sem_dat.data_name{1}).label(find(cell2mat(new_sig(2:end,2))));
        sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
    end
    
    if ~isempty(sem_dat.(sem_dat.data_name{1}).label)
        
        % Stats
        if has_aud && has_vis
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
        elseif has_vis
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
        elseif has_aud
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
        end
        
        % Plot Vis & Aud Channels side-by-side
        if has_aud && has_vis
            
            cfg = [];
            cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
            cfg.type      = 'chan';
            cfg.alt_eve   = {'rep_all_eve'};
            cfg.eve       = {[111 112] [113 114]};
            cfg.plt_lbl   = { 'Visual' 'Auditory'};
            cfg.plt_dim   = [1 2];
            cfg.lgd       = 0;
            cfg.std_err   = 1;
            cfg.lnstyle.col_ord = {rgb('green') rgb('red') rgb('yellow') rgb('blue')};
            cfg.cnd_nme         = {'V Nov' 'V Rep' 'A Nov'  'A Rep'};
            cfg.leg_pos = [1 2 5 6];
            cfg.stt_dat = {cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'vis_dif_rep_stt05' 'vis_dif_rep_stt10'}) ...
                cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'aud_dif_rep_stt05' 'aud_dif_rep_stt10'})};
            cfg.stt_col = {{ft_stt_col(rgb('dark red')) ft_stt_col(rgb('reddish grey'))} ...
                {ft_stt_col(rgb('dark blue')) ft_stt_col(rgb('greyish blue'))}};
            cfg.print      = 1;
            cfg.nofig      = 1;
            cfg.print_type = 'jpg';
            cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/vis_aud_cmp' '/' 'hgp'];
            cfg.prefix     = strcat(sbj,{'_'});
            mmil_ieeg_sensor_plot_v4(cfg)
            
        elseif has_vis
            
            cfg = [];
            cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
            cfg.type      = 'chan';
            cfg.alt_eve   = {'rep_all_eve'};
            cfg.eve       = {[111 112] [113 114]};
            cfg.plt_lbl   = { 'Visual' 'Auditory'};
            cfg.plt_dim   = [1 2];
            cfg.emp_plt   = [0 1];
            cfg.lgd       = 0;
            cfg.std_err   = 1;
            cfg.lnstyle.col_ord = {rgb('green') rgb('red') rgb('yellow') rgb('blue')};
            cfg.cnd_nme         = {'V Nov' 'V Rep' 'A Nov'  'A Rep'};
            cfg.leg_pos = [1 2];
            cfg.stt_dat = {cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'vis_dif_rep_stt05' 'vis_dif_rep_stt10'}),{}};
            cfg.stt_col = {{ft_stt_col(rgb('dark red')) ft_stt_col(rgb('reddish grey'))},{}};
            cfg.print      = 1;
            cfg.nofig      = 1;
            cfg.print_type = 'jpg';
            cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/vis_aud_cmp' '/' 'hgp'];
            cfg.prefix     = strcat(sbj,{'_'});
            mmil_ieeg_sensor_plot_v4(cfg)
            
        elseif has_aud
            
            cfg = [];
            cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
            cfg.type      = 'chan';
            cfg.alt_eve   = {'rep_all_eve'};
            cfg.eve       = {[111 112] [113 114]};
            cfg.plt_lbl   = {'Auditory' 'Visual'};
            cfg.plt_dim   = [1 2];
            cfg.emp_plt   = [1 0];
            cfg.lgd       = 0;
            cfg.std_err   = 1;
            cfg.lnstyle.col_ord = {rgb('green') rgb('red') rgb('yellow') rgb('blue')};
            cfg.cnd_nme         = {'V Nov' 'V Rep' 'A Nov'  'A Rep'};
            cfg.leg_pos = [1 2];
            cfg.stt_dat = {{},cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'aud_dif_rep_stt05' 'aud_dif_rep_stt10'})};
            cfg.stt_col = {{},{ft_stt_col(rgb('dark blue')) ft_stt_col(rgb('bluish grey'))}};
            cfg.print      = 1;
            cfg.nofig      = 1;
            cfg.print_type = 'jpg';
            cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/vis_aud_cmp' '/' 'hgp'];
            cfg.prefix     = strcat(sbj,{'_'});
            mmil_ieeg_sensor_plot_v4(cfg)
            
        end
    end
end

%% Load each subjects data & sig_chn lfp
clear; clc;

clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/clerical/';
sbj_nme = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/' '/' 'subjects']);

indir       = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';
abs_out_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract';

for iS = [1:3 5:16]%[1:3 5:16]
    
    sbj = sbj_nme{iS};
    
    % Load data & sig_chn
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [indir '/' sbj '_lfp_data.mat'];
    sem_dat  = ft_func([],cfg);
    
    new_sig = mmil_readtext([clr_fld '/'  'sig_chn' '/' sbj '/' sem_dat.data_name{1}],['\t']);
    
    has_vis = any(strcmpi('visual_early',new_sig(1,:)));
    has_aud = any(strcmpi('auditory_early',new_sig(1,:)));
    
    % Remove non-sig-chn
    if has_aud && has_vis
        cfg = [];
        cfg.channel = sem_dat.(sem_dat.data_name{1}).label(find(cell2mat(new_sig(2:end,2)) | cell2mat(new_sig(2:end,4))));
        sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
    elseif has_aud || has_vis
        cfg = [];
        cfg.channel = sem_dat.(sem_dat.data_name{1}).label(find(cell2mat(new_sig(2:end,2))));
        sem_dat = ft_func(@ft_preprocessing,cfg,sem_dat);
    end
    
    if ~isempty(sem_dat.(sem_dat.data_name{1}).label)
        
        % Stats
        if has_aud && has_vis
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
        elseif has_vis
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [111 112];
            cfg.add_stt  = 'vis_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
            
        elseif has_aud
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt05';
            cfg.alt_eve  = 'rep_all_eve';
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
            cfg = [];
            cfg.events   = [113 114];
            cfg.add_stt  = 'aud_dif_rep_stt10';
            cfg.alt_eve  = 'rep_all_eve';
            cfg.alpha    = 0.10;
            sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
            
        end
        
        % Plot Vis & Aud Channels side-by-side
        if has_aud && has_vis
            
            cfg = [];
            cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
            cfg.type      = 'chan';
            cfg.alt_eve   = {'rep_all_eve'};
            cfg.eve       = {[111 112] [113 114]};
            cfg.plt_lbl   = { 'Visual' 'Auditory'};
            cfg.plt_dim   = [1 2];
            cfg.lgd       = 0;
            cfg.std_err   = 1;
            cfg.lnstyle.col_ord = {rgb('green') rgb('red') rgb('yellow') rgb('blue')};
            cfg.cnd_nme         = {'V Nov' 'V Rep' 'A Nov'  'A Rep'};
            cfg.leg_pos = [1 2 5 6];
            cfg.stt_dat = {cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'vis_dif_rep_stt05' 'vis_dif_rep_stt10'}) ...
                cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'aud_dif_rep_stt05' 'aud_dif_rep_stt10'})};
            cfg.stt_col = {{ft_stt_col(rgb('dark red')) ft_stt_col(rgb('reddish grey'))} ...
                {ft_stt_col(rgb('dark blue')) ft_stt_col(rgb('greyish blue'))}};
            cfg.print      = 1;
            cfg.nofig      = 1;
            cfg.print_type = 'jpg';
            cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/vis_aud_cmp' '/' 'lfp'];
            cfg.prefix     = strcat(sbj,{'_'});
            mmil_ieeg_sensor_plot_v4(cfg)
            
        elseif has_vis
            
            cfg = [];
            cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
            cfg.type      = 'chan';
            cfg.alt_eve   = {'rep_all_eve'};
            cfg.eve       = {[111 112] [113 114]};
            cfg.plt_lbl   = { 'Visual' 'Auditory'};
            cfg.plt_dim   = [1 2];
            cfg.emp_plt   = [0 1];
            cfg.lgd       = 0;
            cfg.std_err   = 1;
            cfg.lnstyle.col_ord = {rgb('green') rgb('red') rgb('yellow') rgb('blue')};
            cfg.cnd_nme         = {'V Nov' 'V Rep' 'A Nov'  'A Rep'};
            cfg.leg_pos = [1 2];
            cfg.stt_dat = {cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'vis_dif_rep_stt05' 'vis_dif_rep_stt10'}),{}};
            cfg.stt_col = {{ft_stt_col(rgb('dark red')) ft_stt_col(rgb('reddish grey'))},{}};
            cfg.print      = 1;
            cfg.nofig      = 1;
            cfg.print_type = 'jpg';
            cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/vis_aud_cmp' '/' 'lfp'];
            cfg.prefix     = strcat(sbj,{'_'});
            mmil_ieeg_sensor_plot_v4(cfg)
            
        elseif has_aud
            
            cfg = [];
            cfg.dat       = {sem_dat.(sem_dat.data_name{1})};
            cfg.type      = 'chan';
            cfg.alt_eve   = {'rep_all_eve'};
            cfg.eve       = {[111 112] [113 114]};
            cfg.plt_lbl   = {'Auditory' 'Visual'};
            cfg.plt_dim   = [1 2];
            cfg.emp_plt   = [1 0];
            cfg.lgd       = 0;
            cfg.std_err   = 1;
            cfg.lnstyle.col_ord = {rgb('green') rgb('red') rgb('yellow') rgb('blue')};
            cfg.cnd_nme         = {'V Nov' 'V Rep' 'A Nov'  'A Rep'};
            cfg.leg_pos = [1 2];
            cfg.stt_dat = {{},cellfun(@(x) sem_dat.(sem_dat.data_name{1}).cfg.alt_stt.(x),{'aud_dif_rep_stt05' 'aud_dif_rep_stt10'})};
            cfg.stt_col = {{},{ft_stt_col(rgb('dark blue')) ft_stt_col(rgb('bluish grey'))}};
            cfg.print      = 1;
            cfg.nofig      = 1;
            cfg.print_type = 'jpg';
            cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_abstract/vis_aud_cmp' '/' 'lfp'];
            cfg.prefix     = strcat(sbj,{'_'});
            mmil_ieeg_sensor_plot_v4(cfg)
            
        end
    end
end

