function iSASZ_stats_word_char(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.sbj_dat_hld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Setup Events
if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo < 10)
    % Novel Setup
    cfg = [];
    cfg.eve_rul     = { 'new'         };
    cfg.old_events  = { { [1 2 3 4] } };
    cfg.new_events  = { [ 102 ]       };
    cfg.crt_alt_eve = { 'vis_new_bse' };
    cfg.use_alt_eve = { 'trialinfo'   };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    % Add in orthographic statistics
    cfg = [];
    cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl'};
    cfg.new_events  = {[131 131 0 132 132] [141 141 0 142 142] [151 151 0 152 152] [161 161 0 162 162]};
    cfg.crt_alt_eve = {'vis_lng_med'       'vis_big_med'       'vis_frq_med'       'vis_ngh_med'};
    cfg.fld_nme     = {'vis_lng'           'vis_big'           'vis_frq'           'vis_ngh'};
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    % Add in orthographic control statistics
    eve_loc = { 'vis_lng_med' 'vis_big_med' 'vis_frq_med' 'vis_ngh_med' };
    num_loc = { 'vis_lng'     'vis_big'     'vis_frq'     'vis_ngh' };
    
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
end

if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo > 10)
    % Novel Setup
    cfg = [];
    cfg.eve_rul     = { 'new'             };
    cfg.old_events  = { { [11 12 13 14] } };
    cfg.new_events  = { [ 202 ]           };
    cfg.crt_alt_eve = { 'aud_new_bse'     };
    cfg.use_alt_eve = { 'trialinfo'       };
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);

    % Add in phonetic statistics
    cfg = [];
    cfg.eve_rul     = {'med_spl'           'med_spl'           'med_spl'           'med_spl' };
    cfg.new_events  = {[231 231 0 232 232] [241 241 0 242 242] [251 251 0 252 252] [261 261 0 262 262]};
    cfg.crt_alt_eve = {'aud_lng_med'       'aud_big_med'       'aud_frq_med'       'aud_ngh_med'};
    cfg.fld_nme     = {'aud_lng'           'aud_big'           'aud_frq'           'aud_ngh'};
    bcc_dat = ft_func(@ft_redefine_events2,cfg,bcc_dat);
    
    % Add in phonetic control statistics
    eve_loc = { 'aud_lng_med' 'aud_big_med' 'aud_frq_med' 'aud_ngh_med' };
    num_loc = { 'aud_lng'     'aud_big'     'aud_frq'     'aud_ngh' };
    
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
    
end

%% Orthographic Word Stats
if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo < 10)
    
    % Novel Baseline
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_ovr_new'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Vis Length   
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_lng_big' 'sasz_vis_lng_frq' 'sasz_vis_lng_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Vis Bigram  
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_big_lng' 'sasz_vis_big_frq' 'sasz_vis_big_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Vis Frequency   
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_frq_lng' 'sasz_vis_frq_big' 'sasz_vis_frq_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Vis Neigborhood
    cfg = [];
    cfg.stt_fnc  = {'sasz_vis_ngh_lng' 'sasz_vis_ngh_big' 'sasz_vis_ngh_frq'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
end

%% Auditory Word Stats
if any(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo > 10)

    % Novel Baseline
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_ovr_new'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Aud Length   
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_lng_big' 'sasz_aud_lng_frq' 'sasz_aud_lng_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Aud Bigram
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_big_lng' 'sasz_aud_big_frq' 'sasz_aud_big_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Aud Frequency   
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_frq_lng' 'sasz_aud_frq_big' 'sasz_aud_frq_ngh'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
    % Aud Neigborhood    
    cfg = [];
    cfg.stt_fnc  = {'sasz_aud_ngh_lng' 'sasz_aud_ngh_big' 'sasz_aud_ngh_frq'};
    cfg.loc      = fcfg.loc;
    cfg.fld_nme  = bcc_dat.data_name;
    cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
    bcc_dat      = ft_func(@mmil_cloud_stat,cfg,bcc_dat);
    
end

%% Save Data
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.dat_fld '/' sbj '_overall_data.mat'];
ft_func([],cfg,bcc_dat);

end