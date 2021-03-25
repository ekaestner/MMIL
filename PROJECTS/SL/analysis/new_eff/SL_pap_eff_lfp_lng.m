function SL_pap_eff_lfp_lng(fcfg)

fprintf([fcfg.sbj_nme ': Starting SL_pap_eff_lfp_lng on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Note Lexical Decision Stats
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 1;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_nse_stt_msk_anv' ... 
                'aud_nse_stt_msk_anv' };
cfg.cmp_trl = { 'vis_tot_nse' ... 
                'aud_tot_nse' };
cfg.cmp_stt = { [1] ... 
                [2] };
cfg.cmp_nme = { { 'vis_lng_450ms'  } ; ... 
                { 'aud_lng_450ms'  } };  
cfg.cmp     = { { '101!102' } ... 
                { '111!112' } };
cfg.tme_win = { { [0.025 0.500] } ... 
                { [0.450 0.950] } };
cfg.stt_col = { { 'red'  } ... 
                { 'blue' }     };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(1),'uni',0);
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