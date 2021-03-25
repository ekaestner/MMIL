function FW_stats_word_char(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Setup Events
cfg = [];
cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl'};
cfg.new_events  = {[111 111 0 112 112] [121 121 0 122 122] [131 131 0 132 132] [141 141 0 142 142]};
cfg.crt_alt_eve = {'wrd_lng_med'       'wrd_frq_med'       'wrd_ngh_med'       'wrd_big_med'};
cfg.fld_nme     = {'wrd_lng'           'wrd_frq'           'wrd_ngh'           'wrd_big'};
bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);

eve_loc = { 'wrd_lng_med' 'wrd_big_med' 'wrd_frq_med' 'wrd_ngh_med' };
num_loc = { 'wrd_lng'     'wrd_big'     'wrd_frq'     'wrd_ngh' };

cmp = nchoosek(1:4,2);
cmp = [cmp ; fliplr(cmp)];
fix = size(cmp,1);

for iC = 1:fix
    
    eve_one     = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)})(~isnan(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.(eve_loc{cmp(iC,1)}))));
    
    cfg = [];
    cfg.eve_rul     = { 'eve_sql' };
    cfg.old_events  = { { eve_one(2) eve_one(3) } };
    cfg.use_alt_eve = { eve_loc{cmp(iC,1)} };
    cfg.fld_nme     = { num_loc{cmp(iC,2)} };
    cfg.crt_alt_eve = { [ eve_loc{cmp(iC,1)} '_crr_by_' num_loc{cmp(iC,2)}] };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
end

%% Orthographic Word Stats
% Length
cfg = [];
cfg.stt_fnc  = {'fw_wrd_lng'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'fw_wrd_lng_big' 'fw_wrd_lng_frq' 'fw_wrd_lng_ngh'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Bigram
cfg = [];
cfg.stt_fnc  = {'fw_wrd_big'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'fw_wrd_big_lng' 'fw_wrd_big_frq' 'fw_wrd_big_ngh'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Frequency
cfg = [];
cfg.stt_fnc  = {'fw_wrd_frq'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'fw_wrd_frq_lng' 'fw_wrd_frq_big' 'fw_wrd_frq_ngh'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

% Neighborhood
cfg = [];
cfg.stt_fnc  = {'fw_wrd_ngh'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

cfg = [];
cfg.stt_fnc  = {'fw_wrd_ngh_lng' 'fw_wrd_ngh_big' 'fw_wrd_ngh_frq'};
cfg.loc      = 'local';
cfg.fld_nme  = bcc_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end