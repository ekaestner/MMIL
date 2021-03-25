clear; clc;

sbj_num = 1;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

fprintf('Starting work on %s \n',subj)

infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/' subj '_overall_data.mat'];
outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/';

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' subj '_overall_data.mat'];
sll_dat  = ft_func([],cfg);

%% Baseline Stats
% stats
cfg = [];
cfg.events   = [101 102];
cfg.add_stt  = 'vis_nse_stt';
cfg.alt_eve  = 'vis_tot_nse';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

cfg = [];
cfg.events   = [111 112];
cfg.add_stt  = 'aud_nse_stt';
cfg.alt_eve  = 'aud_tot_nse';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

%% Mismatch Stats
% stats
cfg = [];
cfg.events   = [1 2];
cfg.add_stt  = 'vis_mtc_stt';
cfg.alt_eve  = 'trialinfo';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

%% All 4 stats
% stats
cfg = [];
cfg.events   = [1 2 3 4];
cfg.add_stt  = 'vis_ovr_stt';
cfg.alt_eve  = 'trialinfo';
sll_dat      = ft_func(@ft_fdrstats,cfg,sll_dat);

%% Mask stats based on All 4 stats
cfg     = [];
cfg.stt     = {'vis_nse_stt' 'aud_nse_stt' 'vis_mtc_stt'};
cfg.stt_msk = {'vis_ovr_stt' 'vis_ovr_stt' 'vis_ovr_stt'};
sll_dat = ft_func(@ft_mask_stats,cfg,sll_dat);

%% Choose Channels based on hypotheses
cfg = [];
cfg.alt_stt = {'vis_ovr_stt' ...
    'vis_nse_stt_msk' 'vis_nse_stt_msk' ...
    'aud_nse_stt_msk' 'aud_nse_stt_msk' ...
    'vis_mtc_stt_msk' 'vis_mtc_stt_msk'};
cfg.alt_stt_col = {[0.7 0.7 0.7] ...
    ft_stt_col(rgb('magenta')) ft_stt_col(rgb('magenta')) ...
    ft_stt_col(rgb('cyan')) ft_stt_col(rgb('cyan')) ...
    ft_stt_col(rgb('green')) ft_stt_col(rgb('green'))};
cfg.cmp_stt = [1 ...
    2 3 ...
    4 5 ...
    6 7];
cfg.cmp_trl = {'trialinfo' ...
    'vis_nse' 'vis_nse' ...
    'aud_nse' 'aud_nse' ...
    'trialinfo' 'trialinfo'};
cfg.cmp     = {'1!2' ...
    '101!102' '101!102' ...
    '121!122' '121!122' ...
    '1!2' '1!2'};
cfg.cmp_nme = {'Vis_All' ...
    'Ely_Vis_Nse' 'Lte_Vis_Nse' ...
    'Ely_Aud_Nse' 'Lte_Aud_Nse'...
    'Ely_Vis_Mtc' 'Lte_Mtc'};
cfg.tme_win = {[0.100 1.800] ...
    [0.050 0.500] [0.500 1.800] ...
    [0.500 0.950] [0.950 1.800] ....
    [0.500 1.000] [1.000 1.800]};
cfg.chn_loc = 'repair_macro_ejk1st_meso_label';
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical';
cfg.sbj_nme = subj;
cfg.typ      = sll_dat.data_name;
cfg.specific = {'typ' ; 1:numel(sll_dat.data_name)};
cfg.ovr_wrt  = 1;
sll_dat = ft_func(@ft_choose_channel,cfg,sll_dat);

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'sll_dat';
cfg.save     = 'yes';
cfg.filename =[outpath '/' subj '_overall_data.mat'];
ft_func([],cfg,sll_dat);
