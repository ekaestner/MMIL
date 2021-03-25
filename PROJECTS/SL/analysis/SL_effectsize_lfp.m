function SL_effectsize_lfp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Note Overall Effect
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 1;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'ovr_stt' };
cfg.cmp_trl = { 'ovr_all_eve' };
cfg.cmp_stt = { [ 1 ] };
cfg.cmp_nme = { { 'ovr_stt' } };  
cfg.cmp     = { { '101!999' } };
cfg.tme_win = { { [0.025 1.500] } };
cfg.stt_col = { { 'reddish grey' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
cfg.specific = {'typ' ; 1:numel(cfg.data_name)};

dat          = ft_func(@mmil_chn_eff_sze,cfg,bcc_dat);

typ     = cfg.typ;
alt_stt = cfg.alt_stt;
clear fle_chk

%% Linguistic Effects
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 1;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_nse_stt_msk' ...
                'aud_nse_stt_msk' };
cfg.cmp_trl = { 'vis_tot_nse' ...
                'aud_tot_nse' };
cfg.cmp_stt = { [ 1 ] ...
                [ 2 ] };
cfg.cmp_nme = { { 'ltr' } ...
                { 'vce'} };  
cfg.cmp     = { { '101!102' } ...
                { '111!112' } };
cfg.tme_win = { { [0.025 1.500] } ...
                { [0.025 1.500]  } };
cfg.stt_col = { { 'bright red'  } ...
                { 'bright blue' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
cfg.specific = {'typ' ; 1:numel(cfg.data_name)};

dat          = ft_func(@mmil_chn_eff_sze,cfg,bcc_dat);

typ     = cfg.typ;
alt_stt = cfg.alt_stt;
clear fle_chk

%% Mismatch/Match Effects
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 1;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_mtc_stt_msk' };
cfg.cmp_trl = { 'trialinfo' };
cfg.cmp_stt = { [ 1 ] };
cfg.cmp_nme = { { 'mtc' } };  
cfg.cmp     = { { '1!2' } };
cfg.tme_win = { { [0.025 1.450] } };
cfg.stt_col = { { 'green' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
cfg.specific = {'typ' ; 1:numel(cfg.data_name)};

dat          = ft_func(@mmil_chn_eff_sze,cfg,bcc_dat);

typ     = cfg.typ;
alt_stt = cfg.alt_stt;
clear fle_chk

end