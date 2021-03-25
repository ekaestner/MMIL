clear; clc;

sbj_num = 1;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

fprintf('Starting work on %s \n',subj)

infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/' subj '_overall_data.mat'];
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' subj '_overall_data.mat'];
sll_dat  = ft_func([],cfg);

%% Baseline Stats
% stats
cfg = [];
cfg.events   = [101 102];
cfg.add_stt  = 'vis_nse_stt';
cfg.alt_eve  = 'vis_nse';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

cfg = [];
cfg.events   = [111 112];
cfg.add_stt  = 'aud_nse_stt';
cfg.alt_eve  = 'aud_nse';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

%% Mismatch Stats
% stats
cfg = [];
cfg.events   = [1 2];
cfg.add_stt  = 'vis_mtc_stt';
cfg.alt_eve  = 'trialinfo';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

cfg = [];
cfg.events   = [11 12];
cfg.add_stt  = 'aud_mtc_stt';
cfg.alt_eve  = 'trialinfo';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

%% All 4 stats
% stats
cfg = [];
cfg.events   = [1 2 3 4];
cfg.add_stt  = 'vis_ovr_stt';
cfg.alt_eve  = 'trialinfo';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

cfg = [];
cfg.events   = [11 12 13 14];
cfg.add_stt  = 'aud_ovr_stt';
cfg.alt_eve  = 'trialinfo';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

%% Choose Channels based on hypotheses
cfg = [];
cfg.alt_stt = {'vis_ovr_stt' 'vis_nse_stt' 'vis_mtc_stt' ...
               'aud_ovr_stt' 'aud_nse_stt' 'aud_mtc_stt'};
cfg.alt_stt_col = {[0.7 0.7 0.7] ft_stt_col(rgb('magenta')) ft_stt_col(rgb('green')) ...
                   [0.7 0.7 0.7] ft_stt_col(rgb('cyan')) ft_stt_col(rgb('green'))};
cfg.cmp_stt = [1 2 3 4 5 6];
cfg.cmp_trl = {'trialinfo' 'vis_nse' 'trialinfo' ...
               'trialinfo' 'aud_nse' 'trialinfo'};
cfg.cmp     = {'1!2'   '101!102' '121!122' ...
               '11!12' '111!112' '123!124'};
cfg.cmp_nme = {'Vis_All' 'Vis_Noise' 'Vis_Match' ...
               'Aud_All' 'Aud_Noise' 'Aud_Match'};
cfg.tme_win = {[0.100 1.400] [0.100 1.200] [0.100 1.200] ...
               [0.100 1.400] [0.100 1.200] [0.100 1.200]};
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical';
cfg.sbj_nme = subj;
cfg.typ      = sll_dat.data_name;
cfg.specific = {'typ' ; 1:numel(sll_dat.data_name)};
sll_dat = ft_func(@ft_choose_channel,cfg,sll_dat);

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'sll_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sll_dat);
