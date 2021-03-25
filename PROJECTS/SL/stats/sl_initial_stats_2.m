function sl_initial_stats_2(fcfg)
%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

infile = [fcfg.out_pth '/' sbj '_overall_data.mat'];
outpath = fcfg.out_pth;

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' sbj '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Run Stats
% Main events
cfg = [];
cfg.stt_fnc  = {'sll_ovr_all'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sll_vis_anv_stt' 'sll_vis_anv_stt_01'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sll_vis_nse'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sll_aud_nse'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sll_mtc_mss'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Vis events
cfg = [];
cfg.stt_fnc  = {'sl_vis_wrd'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% cfg = [];
% cfg.stt_fnc  = {'sl_vis_con'};
% cfg.loc      = 'local';
% cfg.fld_nme  = bcc_dat.data_name;
% cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
% bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sl_vis_big'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sl_vis_ngh'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Auditory events
cfg = [];
cfg.stt_fnc  = {'sl_aud_wrd'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% cfg = [];
% cfg.stt_fnc  = {'sl_aud_con'};
% cfg.loc      = 'local';
% cfg.fld_nme  = bcc_dat.data_name;
% cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
% bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sl_aud_big'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'sl_aud_ngh'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename =[outpath '/' sbj '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end
    