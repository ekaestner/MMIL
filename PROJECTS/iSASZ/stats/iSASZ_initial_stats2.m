function iSASZ_initial_stats2(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

infile  = [fcfg.fle_out_pth '/' 'epoch_data' '/' sbj '_overall_data.mat'];
outpath = [fcfg.fle_out_pth '/' 'epoch_data' '/'];

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' sbj '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

tsk         = mmil_load_subj_info([fcfg.clr_fld '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); %mmil_readtext([fcfg.clr_fld '/tsk/' subj]);     
tsk         = repmat(tsk,1,2);

cfg = [];
cfg.channel = bcc_dat.(bcc_dat.data_name{1}).label([78 91 49 21 22 47 55]);
bcc_dat     = ft_func(@ft_preprocessing,cfg,bcc_dat);

%% Run Stats
if any(strcmpi(tsk,'sz'))
    
    % Overall
    tic
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sz'));
    cfg.stt_fnc   = {'iSASZ_vis_ovr'}; % add in voice vs fst, voice vs vocode, voice vs tone, vocode vs tone
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);

    % Repetition
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sz'));
    cfg.stt_fnc   = {'iSASZ_vis_new_old' 'iSASZ_vis_new_old_sht'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantic
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sz'));
    cfg.stt_fnc   = {'iSASZ_vis_ani_obj' 'iSASZ_vis_cat_prm_new' 'iSASZ_vis_sem_mot' 'iSASZ_vis_sem_mnp'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Emotion
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sz'));
    cfg.stt_fnc   = {'iSASZ_vis_val' 'iSASZ_vis_ars'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Word Stats
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sz'));
    cfg.stt_fnc   = {'iSASZ_vis_FREQ_med' 'iSASZ_vis_LEN_med' 'iSASZ_vis_UN2_C_med'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    toc
    
end

if any(strcmpi(tsk,'sa'))
    
    % Overall
    tic
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sa'));
    cfg.stt_fnc   = {'iSASZ_aud_ovr'}; % add in voice vs fst, voice vs vocode, voice vs tone, vocode vs tone
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Repetition
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sa'));
    cfg.stt_fnc   = {'iSASZ_aud_new_old' 'iSASZ_aud_new_old_sht'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Semantic
    cfg.data_name = find(strcmpi(tsk,'sa'));
    cfg.stt_fnc   = {'iSASZ_aud_ani_obj' 'iSASZ_aud_cat_prm_new' 'iSASZ_aud_sem_mot' 'iSASZ_aud_sem_mnp'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Emotion
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sa'));
    cfg.stt_fnc   = {'iSASZ_aud_val' 'iSASZ_aud_ars'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Word Stats
    cfg = [];
    cfg.data_name = find(strcmpi(tsk,'sa'));
    cfg.stt_fnc   = {'iSASZ_aud_FREQ_med' 'iSASZ_aud_LEN_med' 'iSASZ_aud_UN2_C_med'};
    cfg.loc       = 'local';
    cfg.fld_nme   = bcc_dat.data_name(cfg.data_name);
    cfg.specific  = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat       = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    toc
    
end
 

%% Mask Stats
if any(strcmpi(tsk,'sz'))
    cfg     = [];
    cfg.stt     = { 'new_old_stt' 'new_old_sht_stt' 'ani_obj_stt' 'cat_prm_new_stt' 'FREQ_med_stt' 'LEN_med_stt' 'UN2_C_med_stt'};
    cfg.stt_msk = { 'ovr_stt'     'ovr_stt'         'ovr_stt'     'ovr_stt'         'ovr_stt'      'ovr_stt'     'ovr_stt'};
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
end

if any(strcmpi(tsk,'sa'))
    cfg     = [];
    cfg.stt     = { 'new_old_stt' 'new_old_sht_stt' 'ani_obj_stt' 'cat_prm_new_stt' 'FREQ_med_stt' 'LEN_med_stt' 'UN2_C_med_stt'};
    cfg.stt_msk = { 'ovr_stt'     'ovr_stt'         'ovr_stt'     'ovr_stt'         'ovr_stt'      'ovr_stt'     'ovr_stt'};
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
% cfg.sve_app  = 'app_all';
cfg.filename =[outpath '/' sbj '_overall_test_data.mat'];
ft_func([],cfg,bcc_dat);

end