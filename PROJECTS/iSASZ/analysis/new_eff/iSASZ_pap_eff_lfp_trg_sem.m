function iSASZ_pap_eff_lfp_trg_sem(fcfg)

fprintf([fcfg.sbj_nme ': Starting iSASZ_pap_eff_lfp_trg_sem on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name);
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label;
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Visual Responsive Electrodes
eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

if any(eve_typ < 9)
    
    cfg = [];
    
    cfg.ovr_wrt = 1;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.data_name = 1;
    
    cfg.stt_lab = 'stt_lab';
    
    cfg.alt_stt = { 'vis_sem_mot_mnp_stt_msk_ovr' };
    cfg.cmp_trl = { 'sem_mot_mnp' };
    cfg.cmp_stt = { [ 1 1 ] };
    cfg.cmp_nme = { { 'vis_trg_sem_500ms' 'vis_trg_sem_1200ms' } };
    cfg.cmp     = { { '901!902'       '901!902'        } };
    cfg.tme_win = { { [0.025 0.500]   [0.025 1.200]    } };
    cfg.stt_col = { { 'bright red' 'dark red' } };
    
    cfg.sbj_nme  = fcfg.sbj_nme;
    cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
    cfg.specific = {'typ' ; 1:numel(cfg.typ)};
    
    dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
    
    % Split effects into ecog/depth
    cfg = [];
    
    cfg.sbj_nme = fcfg.sbj_nme;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.dat_typ = {'lfp'};
    
    mmil_split_stat(cfg);
    
    % Split timing into ecog/depth
    cfg = [];
    
    cfg.sbj_nme = fcfg.sbj_nme;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.dat_typ = {'lfp'};
    
    mmil_split_timing(cfg);
    
end

%% Auditory Responsive Electrodes
eve_typ = unique(bcc_dat.(bcc_dat.data_name{1}).cfg.alt_eve.trialinfo);

if any(eve_typ > 9)
    
    cfg = [];
    cfg.ovr_wrt = 1;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.data_name = 1;
    
    cfg.stt_lab = 'stt_lab';
    
    cfg.alt_stt = { 'aud_sem_mot_mnp_stt_msk_ovr' };
    cfg.cmp_trl = { 'sem_mot_mnp' };
    cfg.cmp_stt = { [ 1 1 ] };
    cfg.cmp_nme = { { 'aud_trg_sem_500ms' 'aud_trg_sem_1200ms' } };
    cfg.cmp     = { { '911!912'       '911!912'        } };
    cfg.tme_win = { { [0.025 0.500]   [0.025 1.200]    } };
    cfg.stt_col = { { 'bright blue' 'dark blue' } };
    
    cfg.sbj_nme  = fcfg.sbj_nme;
    cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
    cfg.specific = {'typ' ; 1:numel(cfg.typ)};
    
    dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
    
    % Split effects into ecog/depth
    cfg = [];
    
    cfg.sbj_nme = fcfg.sbj_nme;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.dat_typ = {'lfp'};
    
    mmil_split_stat(cfg);
    
    % Split timing into ecog/depth
    cfg = [];
    
    cfg.sbj_nme = fcfg.sbj_nme;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.dat_typ = {'lfp'};
    
    mmil_split_timing(cfg);
    
end

end