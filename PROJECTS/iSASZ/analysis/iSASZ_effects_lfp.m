function iSASZ_effects_lfp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iD = 1:numel(bcc_dat.data_name)
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iD}).label;
end

eve_typ = unique(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_eve.trialinfo);

%% Visual Effects
if any(eve_typ < 9)
    
    cfg = [];
    cfg.ovr_wrt = 1;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.data_name = 1;
    
    cfg.stt_lab = 'stt_lab';
    cfg.chn_loc = 'stt_lab';
    
    cfg.alt_stt = { 'vis_ovr_stt' ...
                    'vis_ani_obj_stt_msk' ...
                    'vis_new_old_stt_msk' };
    cfg.cmp_trl = { 'vis_ovr' ...
                    'vis_ani_obj' ...
                    'vis_new_old' };
    cfg.cmp_stt = { [ 1 ] ...
                    [ 2 ] ...
                    [ 3 ] };
    cfg.cmp_nme = { { 'vis_ovr_800ms' } ; ...
                    { 'vis_ani_800ms' } ; ...
                    { 'vis_nov_800ms' } };
    cfg.cmp     = { { '101>999' } ...
                    { '121>122' } ...
                    { '111>112' } };
    cfg.tme_win = { { [0.025 0.800] } ...
                    { [0.025 0.800] } ...
                    { [0.025 0.800] } };
    cfg.stt_col = { { 'bright red'   } ...
                    { 'green' }  ...
                    { 'orange' }};
    
    cfg.sbj_nme  = fcfg.sbj_nme;
    cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
    if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'hgp'}; end
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

%% Auditory Effects
if any(eve_typ > 9)
    
    cfg = [];
    cfg.ovr_wrt = 1;
    cfg.clr_fld = fcfg.clr_fld;
    
    cfg.data_name = 1;
    
    cfg.stt_lab = 'stt_lab';
    cfg.chn_loc = 'stt_lab';
    
    cfg.alt_stt = { 'aud_ovr_stt' ...
                    'aud_ani_obj_stt_msk' ...
                    'aud_new_old_stt_msk' };
    cfg.cmp_trl = { 'aud_ovr' ...
                    'aud_ani_obj' ...
                    'aud_new_old' };
    cfg.cmp_stt = { [ 1 ] ...
                    [ 2 ] ...
                    [ 3 ] };
    cfg.cmp_nme = { { 'aud_ovr_800ms' } ; ...
                    { 'aud_ani_800ms' } ; ...
                    { 'aud_nov_800ms' } };
    cfg.cmp     = { { '201>999' } ...
                    { '221>122' } ...
                    { '211>112' } };
    cfg.tme_win = { { [0.025 0.800] } ...
                    { [0.025 0.800] } ...
                    { [0.025 0.800] } };
    cfg.stt_col = { { 'bright blue'   } ...
                    { 'green' }  ...
                    { 'orange' }};
    
    cfg.sbj_nme  = fcfg.sbj_nme;
    cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
    if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'hgp'}; end
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