function FW_motor_stats(fcfg)

%% Load
cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.out_pth '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

%% Events
kep_trl = find(bcc_dat.(bcc_dat.data_name{1}).trialinfo==3 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==6 | bcc_dat.(bcc_dat.data_name{1}).trialinfo==7);

cfg = [];
cfg.trials = kep_trl;
bcc_dat    = ft_func(@ft_redefinetrial,cfg,bcc_dat);

for iD = 1:numel(bcc_dat.data_name); bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo = bcc_dat.(bcc_dat.data_name{iD}).trialinfo; end

%% Stats
cfg = [];
cfg.stt_fnc  = {'fw_vis_mtr_wrd_stt'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'fw_vis_mtr_trg_stt'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

%% Save
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename = [fcfg.out_pth '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);


end