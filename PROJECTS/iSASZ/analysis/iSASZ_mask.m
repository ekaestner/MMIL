function iSASZ_mask(fcfg)

%% Initial Load
fprintf([fcfg.sbj_nme ': Starting initial stats mask work on %s \n'],fcfg.sbj_nme)

infile = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
outpath = fcfg.dat_fld;

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

%% Remove mask
for iD = 1:numel(bcc_dat.data_name)
    msk_nme = fieldnames(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_stt);
    msk_nme = msk_nme(string_find(msk_nme,'_msk'));
    
    for iR = 1:numel(msk_nme)
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_stt = rmfield(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_stt,msk_nme{iR});
    end
end

eve_typ = unique(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo);

%% Mask Stats Visual
if any(eve_typ < 9)
    cfg     = [];
    cfg.stt     = { 'vis_new_old_stt' 'vis_ani_obj_stt' }; %
    cfg.stt_msk = { 'vis_ovr_stt'     'vis_ovr_stt'     }; %
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
end

%% Mask Stats Auditory
if any(eve_typ > 9)
    cfg     = [];
    cfg.stt     = { 'aud_new_old_stt' 'aud_ani_obj_stt' }; %
    cfg.stt_msk = { 'aud_ovr_stt'     'aud_ovr_stt' }; %
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
end

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename =[outpath '/' fcfg.sbj_nme '_overall_data'];
ft_func([],cfg,bcc_dat);

end