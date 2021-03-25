function FW_add_ort_stt(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Events
for iD = 1:numel(bcc_dat.data_name)

    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh = nansum([bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.nwn_ngh bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.wrd_ngh]')';
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh( ...
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh>=5 & ...
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==3) = 401;
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh( ...
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh==0 & ...
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==3) = 402;
    
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh( ...
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh==0 & ...
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo==5) = 403;
   
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh( ...
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.tot_ngh<400) = nan;
    
end

%% Stats
% Make
cfg = [];
cfg.stt_fnc  = {'fw_vis_ort_tot'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'fw_vis_ort_nwd_stt'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'fw_vis_ort_wrd_stt'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end