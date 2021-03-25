clear; clc;

sbj_num = 17;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

fprintf('Starting work on %s \n',subj)

infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/' subj '_overall_data.mat'];
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' subj '_overall_data.mat'];
sem_dat  = ft_func([],cfg);

has_vis = any(sem_dat.(sem_dat.data_name{1}).trialinfo < 9);
has_aud = any(sem_dat.(sem_dat.data_name{1}).trialinfo > 9);

%% Baseline Stats
% Stats
if has_aud
    cfg = [];
    cfg.basetime = [-0.40 0];
    cfg.events   = [102];
    cfg.add_stt  = 'aud_ovr_all_stt';
    cfg.alt_eve  = 'ovr_all_evt';
    cfg.alpha = .01;
    sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
end
if has_vis
    cfg = [];
    cfg.basetime = [-0.40 0];
    cfg.events   = [101];
    cfg.add_stt  = 'vis_ovr_all_stt';
    cfg.alt_eve  = 'ovr_all_evt';
    cfg.alpha    = .01;
    sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
end

%% Repetition N4 stats
% Events
if has_vis && has_aud
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events  = {[1 2 3 4] [5 6 7 8] [11 12 13 14] [15 16 17 18]};
    cfg.new_events  = {111 112 113 114};
    cfg.crt_alt_eve = 'rep_all_eve';
    sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
elseif has_vis
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events  = {[1 2 3 4] [5 6 7 8]};
    cfg.new_events  = {111 112};
    cfg.crt_alt_eve = 'rep_all_eve';
    sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
elseif has_aud
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events  = {[11 12 13 14] [15 16 17 18]};
    cfg.new_events  = {113 114};
    cfg.crt_alt_eve = 'rep_all_eve';
    sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
end

% Stats
if has_aud
    cfg = [];
    cfg.events   = [113 114];
    cfg.add_stt  = 'aud_dif_rep_stt';
    cfg.alt_eve  = 'rep_all_eve';
    sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
end
if has_vis
    cfg = [];
    cfg.events   = [111 112];
    cfg.add_stt  = 'vis_dif_rep_stt';
    cfg.alt_eve  = 'rep_all_eve';
    sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
end

%% Semantic Category Stats
% Events
if has_vis && has_aud
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events  = {[2 4] [1 3] [12 14] [11 13]};
    cfg.new_events  = {121 122 123 124};
    cfg.crt_alt_eve = 'sem_all_eve';
    sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
elseif has_vis
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events  = {[2 4] [1 3]};
    cfg.new_events  = {121 122};
    cfg.crt_alt_eve = 'sem_all_eve';
    sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
elseif has_aud
    cfg = [];
    cfg.return_events = 0;
    cfg.old_events  = {[12 14] [11 13]};
    cfg.new_events  = {123 124};
    cfg.crt_alt_eve = 'sem_all_eve';
    sem_dat = ft_func(@ft_redefine_events,cfg,sem_dat);
end

% Stats
if has_aud
    cfg = [];
    cfg.events   = [123 124];
    cfg.add_stt  = 'aud_dif_sem_stt';
    cfg.alt_eve  = 'sem_all_eve';
    sem_dat      = ft_func(@ft_fdrstats,cfg,sem_dat);
end
if has_vis
    cfg = [];
    cfg.events   = [121 122];
    cfg.add_stt  = 'vis_dif_sem_stt';
    cfg.alt_eve  = 'sem_all_eve';
    sem_dat  = ft_func(@ft_fdrstats,cfg,sem_dat);
end

%% Choose Channels based on hypotheses
if has_vis && has_aud
    
    cfg = [];
    cfg.alt_stt = {'vis_ovr_all_stt' 'vis_dif_rep_stt' 'vis_dif_sem_stt' ...
        'aud_ovr_all_stt' 'aud_dif_rep_stt' 'aud_dif_sem_stt'};
    cfg.alt_stt_col = {[0.7 0.7 0.7] ft_stt_col(rgb('red')) ft_stt_col(rgb('green')) ...
        [0.7 0.7 0.7] ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))};
    cfg.cmp_stt = [1 2 3 4 5 6];
    cfg.cmp_trl = {'ovr_all_evt' 'rep_all_eve' 'sem_all_eve' ...
        'ovr_all_evt' 'rep_all_eve' 'sem_all_eve'};
    cfg.cmp     = {'101!999' '111!112' '121!122' ...
        '102!999' '113!114' '123!124'};
    cfg.cmp_nme = {'visual_overall_active' 'visual_repetition_differences' 'visual_semantic_differences' ...
        'auditory_overall_active' 'auditory_repetition_differences' 'auditory_semantic_differences'};
    cfg.tme_win = {[0.100 1.200] [0.100 1.200] [0.100 1.200] ...
        [0.100 1.200] [0.100 1.200] [0.100 1.200]};
    cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';
    cfg.sbj_nme = subj;
    cfg.typ      = sem_dat.data_name;
    cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
    sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
    
elseif has_vis
    
    cfg = [];
    cfg.alt_stt     = {'vis_ovr_all_stt' 'vis_dif_rep_stt' 'vis_dif_sem_stt'};
    cfg.alt_stt_col = {[0.7 0.7 0.7] ft_stt_col(rgb('red')) ft_stt_col(rgb('green'))};
    cfg.cmp_stt = [1 2 3];
    cfg.cmp_trl = {'ovr_all_evt' 'rep_all_eve' 'sem_all_eve'};
    cfg.cmp     = {'101!999' '111!112' '121!122'};
    cfg.cmp_nme = {'visual_overall_active' 'visual_repetition_differences' 'visual_semantic_differences'};
    cfg.tme_win = {[0.100 1.200] [0.100 1.200] [0.100 1.200]};
    cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';
    cfg.sbj_nme = subj;
    cfg.typ      = sem_dat.data_name;
    cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
    sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
    
elseif has_aud
    
    cfg = [];
    cfg.alt_stt     = {'aud_ovr_all_stt' 'aud_dif_rep_stt' 'aud_dif_sem_stt'};
    cfg.alt_stt_col = {[0.7 0.7 0.7] ft_stt_col(rgb('blue')) ft_stt_col(rgb('yellow'))};
    cfg.cmp_stt = [1 2 3];
    cfg.cmp_trl = {'ovr_all_evt' 'rep_all_eve' 'sem_all_eve'};
    cfg.cmp     = {'102!999' '113!114' '123!124'};
    cfg.cmp_nme = {'auditory_overall_active' 'auditory_repetition_differences' 'auditory_semantic_differences'};
    cfg.tme_win = {[0.100 1.200] [0.100 1.200] [0.100 1.200]};
    cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';
    cfg.sbj_nme = subj;
    cfg.typ      = sem_dat.data_name;
    cfg.specific = {'typ' ; 1:numel(sem_dat.data_name)};
    sem_dat = ft_func(@ft_choose_channel,cfg,sem_dat);
    
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'sem_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sem_dat);
