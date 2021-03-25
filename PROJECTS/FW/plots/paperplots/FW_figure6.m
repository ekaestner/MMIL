cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/FW/clerical';
cfg.dat_fld = '/home/ekaestne/PROJECTS/OUTPUT/FW//epoch_data';
cfg.out_dir = 'figure6/version4';

cfg.stm_typ = {'Language' 'Motor' ''};

cfg.fst_ltr = 1;

cfg.tsk     = 'FW';
cfg.eff_typ = 'hgp';
cfg.eff_nme = { 'pap_wrd_600' };
cfg.eff_clm = { 3             };
cfg.eff_col = { 'bright red'  };
cfg.non_eff_col = {'purplish red'};

cfg.chn_plt = 0;

cfg.dat_plt = 2;

% cfg.alt_eve = 'trialinfo';
% cfg.eve     = [1 2 3 4];
% cfg.lnstyle.col_ord = {rgb('bright red') rgb('bright blue') rgb('reddish grey') rgb('bluish grey')} ;

% cfg.stt_lab = 'stt_lab';
% cfg.stt_dat = { 'vis_nse_stt_msk_anv'           'aud_nse_stt_msk_anv'          'vis_mtc_stt' };
% cfg.stt_col = { { ft_stt_col(rgb('bright red')) ft_stt_col(rgb('bright blue')) ft_stt_col(rgb('dark yellow')) } };
% cfg.stt_cmp = { { '0%3'                         '3%6'                          '6%9'   } };

% cfg.v_lne       = [0 0.450 0.900];
% cfg.v_lne_wdt   = [3 3 1];
% cfg.v_lne_col   = {rgb('red') rgb('blue') rgb('black')};

% cfg.xlm = [-0.2 1.3];

cfg.inc_reg = { 'Inferior Precentral' 'Middle Precentral' ...
                'Supramarginal' ...
                'Caudal STG' 'Middle STG' 'Rostral STG'...
                'Pars Opercularis' };

cfg.inc_non_stm = 1;
            
mmil_stimulation_plot(cfg)