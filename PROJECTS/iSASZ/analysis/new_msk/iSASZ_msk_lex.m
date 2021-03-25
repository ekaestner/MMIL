function iSASZ_msk_lex(fcfg)

%% Initial Load
fprintf([fcfg.sbj_nme ': Starting iSASZ_msk_lex on %s \n'],fcfg.sbj_nme)

infile = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
outpath = fcfg.dat_fld;

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Mask Stats
eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

if any(eve_typ < 9)
    cfg     = [];
    cfg.stt     = { 'vis_frq_stt' }; %
    cfg.stt_msk = { 'vis_new_ovr_stt' }; %
    cfg.pst_fix = '_msk_ovr';
    bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);
end

if any(eve_typ > 9)
    cfg     = [];
    cfg.stt     = { 'aud_frq_stt' }; %
    cfg.stt_msk = { 'aud_new_ovr_stt' }; %
    cfg.pst_fix = '_msk_ovr';
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