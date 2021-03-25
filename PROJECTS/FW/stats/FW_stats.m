function FW_stats(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Stats
% Make
if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg,'alt_stt')
    vis_ovr_stt = 1;
    fw_ffn_stt = 1;
    vis_stm = 1;
    vis_stm_01  = 1;
    vis_old = 1;
    vis_wrd = 1;
    vis_ltr = 1;
else
    if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_stt,'vis_ovr_stt'); vis_ovr_stt = 1; else vis_ovr_stt = 0; end
    if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_stt,'fw_ffn_stt');  fw_ffn_stt = 1;  else fw_ffn_stt = 0;  end
    if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_stt,'vis_stm');     vis_stm = 1;     else vis_stm = 0; end
    if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_stt,'vis_stm_01');  vis_stm_01  = 1; else vis_stm_01  = 0; end
    if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_stt,'vis_old');     vis_old = 1;     else vis_old = 0; end
    if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_stt,'vis_wrd');     vis_wrd = 1;     else vis_wrd = 0; end
    if ~isfield(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_stt,'vis_ltr');     vis_ltr = 1;     else vis_ltr = 0; end
end

% Overall Baselines
if vis_ovr_stt
    cfg = [];
    cfg.stt_fnc  = {'fw_vis_ovr_stt'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

% Event Baselines
if fw_ffn_stt
    cfg = [];
    cfg.stt_fnc  = {'fw_ffn_stt' 'fw_ltr_stt' 'fw_wrd_stt' 'fw_rep_stt'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

% Main Stats
if vis_stm
    cfg = [];
    cfg.stt_fnc  = {'fw_vis_stm_stt'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

if vis_stm_01
    cfg = [];
    cfg.stt_fnc  = {'fw_vis_stm_stt_01'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

if vis_old
    cfg = [];
    cfg.stt_fnc  = {'fw_vis_old'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

if vis_wrd
    cfg = [];
    cfg.stt_fnc  = {'fw_vis_wrd'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

if vis_ltr
    cfg = [];
    cfg.stt_fnc  = {'fw_vis_ltr'};
    cfg.loc      = 'local';
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end