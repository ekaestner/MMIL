clear; clc;

sbj_num = 1;
ovr_wrt = 1;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

cln_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';

try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical' '/' subj '_loc']);
catch
    has_loc = 0;
end
if ~isvar('has_loc'); has_loc = 1; end

cfg = [];
cfg.load    = 'yes';
cfg.file = [outpath '/' subj '_overall_data.mat'];
sem_dat     = ft_func([],cfg);

if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo > 9); has_aud = 1; else has_aud = 0; end
if any(sem_dat.(sem_dat.data_name{1}).cfg.alt_eve.trialinfo < 9); has_vis = 1; else has_vis = 0; end

fprintf('Starting work on %s \n',subj)

%% Plot Both Chosen & Not Chosen Channels
for iD = 1:numel(sem_dat.data_name) 
    
    % Create Plotting Events
    if has_vis && has_aud
        
        if ovr_wrt == 1;
            if exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Auditory_Visual' '/' sem_dat.data_name{iD}(end-2:end)]); rmdir(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Auditory_Visual' '/' sem_dat.data_name{iD}(end-2:end)],'s'); end
            if exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Visual' '/' sem_dat.data_name{iD}(end-2:end)]); rmdir(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Visual' '/' sem_dat.data_name{iD}(end-2:end)],'s'); end
            if exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Auditory' '/' sem_dat.data_name{iD}(end-2:end)]); rmdir(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Auditory' '/' sem_dat.data_name{iD}(end-2:end)],'s'); end
        end
        
        % Create Plotting Events
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events  = {[2 4] [1 3] [5 6 7 8] [12 14] [11 13] [15 16 17 18]};
        cfg.new_events  = {121 122 112 123 124 114};
        cfg.crt_alt_eve = 'plt_all_eve';
        sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
        
        % Create Auditory & Visual Plot
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{iD})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
            cfg.cmb = 1; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iD}));
        cfg.chn_grp   = grp_chn{iD};
        cfg.alt_eve         = {'rep_all_eve'};
        cfg.eve             = [111 112  ...
            113 114 ];
        cfg.lnstyle.col_ord = {rgb('green') rgb('red')  ...
            rgb('yellow') rgb('blue')};
        cfg.lnstyle.lin_ord = {'-' '-'  ...
            '-.' '-.'};
        cfg.cnd_nme         = {'V Nov' 'V Rep' ...
           'A Nov'  'A Rep'};
        cfg.leg_pos = [1 2 ...
            5 6];
        cfg.alt_lbl = 'aud_vis_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{iD};
        cfg.stt_col = grp_stt_col{iD};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Auditory_Visual' '/' sem_dat.data_name{iD}(end-2:end)];
        cfg.prefix     = strcat(sem_dat.data_name{iD},{'_'},pre_fix{iD});
        mmil_ieeg_sensor_plot_v4(cfg)
        
        % Create Visual Plot
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{iD})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             cfg.cmb = 1; cfg.cmp_ind = 1:3; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iD}));
        cfg.chn_grp   = grp_chn{iD};
        cfg.alt_eve         = {'plt_all_eve'};
        cfg.eve             = [112 121 122];
        cfg.lnstyle.col_ord = {rgb('red') rgb('dark green') rgb('light green')};
        cfg.cnd_nme         = {'V Rep' 'V Obj' 'V Ani'};
        cfg.alt_lbl = 'vis_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{iD};
        cfg.stt_col = grp_stt_col{iD};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' subj '/' 'Visual' '/' sem_dat.data_name{iD}(end-2:end)];
        cfg.prefix     = strcat(sem_dat.data_name{iD},{'_'},pre_fix{iD});
        mmil_ieeg_sensor_plot_v4(cfg)
        
        % Create Auditory Plot
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{iD})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             cfg.cmb = 1; cfg.cmp_ind = 4:6; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iD}));
        cfg.chn_grp   = grp_chn{iD};
        cfg.alt_eve         = {'plt_all_eve'};
        cfg.eve             = [114 123 124];
        cfg.lnstyle.col_ord = {rgb('blue') rgb('dark yellow') rgb('light yellow')};
        cfg.lnstyle.lin_ord = {'-.' '-.' '-.'};
        cfg.cnd_nme         = {'A Rep' 'A Obj' 'A Ani'};
        cfg.alt_lbl = 'aud_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{iD};
        cfg.stt_col = grp_stt_col{iD};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' subj '/' 'Auditory' '/' sem_dat.data_name{iD}(end-2:end)];
        cfg.prefix     = strcat(sem_dat.data_name{iD},{'_'},pre_fix{iD});
        mmil_ieeg_sensor_plot_v4(cfg)
        
    elseif has_vis
        
        if ovr_wrt == 1;
            if exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Visual' '/' sem_dat.data_name{iD}(end-2:end)]); rmdir(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Visual' '/' sem_dat.data_name{iD}(end-2:end)],'s'); end
        end
        
        % Create Plotting Events
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events  = {[2 4] [1 3] [5 6 7 8]};
        cfg.new_events  = {121 122 112};
        cfg.crt_alt_eve = 'plt_all_eve';
        sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
        
        % Create Visual Plot
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{iD})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             cfg.cmb = 1; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iD}));
        cfg.chn_grp   = grp_chn{iD};
        cfg.alt_eve         = {'plt_all_eve'};
        cfg.eve             = [112 121 122];
        cfg.lnstyle.col_ord = {rgb('red') rgb('dark green') rgb('light green')};
        cfg.cnd_nme         = {'V Rep' 'V Obj' 'V Ani'};
        cfg.alt_lbl = 'vis_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{iD};
        cfg.stt_col = grp_stt_col{iD};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' subj '/' 'Visual' '/' sem_dat.data_name{iD}(end-2:end)];
        cfg.prefix     = strcat(sem_dat.data_name{iD},{'_'},pre_fix{iD});
        cfg.print_type = 'jpg';
        mmil_ieeg_sensor_plot_v4(cfg)
        
    elseif has_aud
        
        if ovr_wrt == 1;
            if exist(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Auditory' '/' sem_dat.data_name{iD}(end-2:end)]); rmdir(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/'  subj '/' 'Auditory' '/' sem_dat.data_name{iD}(end-2:end)],'s'); end
        end
        
        % Create Plotting Events
        cfg = [];
        cfg.return_events = 0;
        cfg.old_events  = {[12 14] [11 13] [15 16 17 18]};
        cfg.new_events  = {123 124 114};
        cfg.crt_alt_eve = 'plt_all_eve';
        sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
        
        % Create Auditory Plot
        cfg = [];
        cfg.dat       = {sem_dat.(sem_dat.data_name{iD})};
        cfg.lgd       = 0;
        cfg.plt_dim   = [5 5];
             cfg.cmb = 1; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sem_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sem_dat.(sem_dat.data_name{iD}));
        cfg.chn_grp   = grp_chn{iD};
        cfg.alt_eve         = {'plt_all_eve'};
        cfg.eve             = [114 123 124];
        cfg.lnstyle.col_ord = {rgb('blue') rgb('dark yellow') rgb('light yellow')};
        cfg.lnstyle.lin_ord = {'-.' '-.' '-.'};
        cfg.cnd_nme         = {'A Rep' 'A Obj' 'A Ani'};
        cfg.alt_lbl = 'aud_snr_lab';
        cfg.plt_shf = 1;
        cfg.y_lim   = 'auto';
        cfg.stt_dat = grp_stt{iD};
        cfg.stt_col = grp_stt_col{iD};
        cfg.std_err = 1;
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'jpg';
        cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/initial_plot' '/' subj '/' 'Auditory' '/' sem_dat.data_name{iD}(end-2:end)];
        cfg.prefix     = strcat(sem_dat.data_name{iD},{'_'},pre_fix{iD});
        mmil_ieeg_sensor_plot_v4(cfg)
        
    end
       
end 
