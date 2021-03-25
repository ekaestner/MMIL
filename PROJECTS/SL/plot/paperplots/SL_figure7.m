% cfg = [];
% 
% cfg.chn_crr = 1;
% cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
% cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
% cfg.out_dir = 'figure7';
% 
% cfg.stm_typ = {'Language' 'Motor'};
% 
% cfg.fst_ltr = 1;
% 
% cfg.tsk     = 'SL';
% cfg.eff_typ = 'hgp';
% cfg.eff_nme = { 'pap_lng_950' 'pap_lng_950' 'pap_con_950'};
% cfg.eff_clm = { 1             2             2 };
% cfg.eff_col = { 'bright red'  'bright blue' 'cyan'};
% 
% cfg.chn_plt = 0;
% 
% cfg.dat_plt = 2;
% 
% cfg.alt_eve = 'trialinfo';
% cfg.eve     = [1 2 3 4];
% cfg.lnstyle.col_ord = {rgb('bright red') rgb('bright blue') rgb('reddish grey') rgb('bluish grey')} ;
% 
% cfg.stt_lab = 'stt_lab';
% cfg.stt_dat = { 'vis_nse_stt_msk_anv'           'aud_nse_stt_msk_anv'          'vis_mtc_stt' };
% cfg.stt_col = { { ft_stt_col(rgb('bright red')) ft_stt_col(rgb('bright blue')) ft_stt_col(rgb('dark yellow')) } };
% cfg.stt_cmp = { { '0%3'                         '3%6'                          '6%9'   } };
% 
% cfg.v_lne       = [0 0.450 0.900];
% cfg.v_lne_wdt   = [3 3 1];
% cfg.v_lne_col   = {rgb('red') rgb('blue') rgb('black')};
% 
% cfg.xlm = [-0.2 1.3];
% 
% cfg.rmv_sbj = '';
% 
% mmil_stimulation_plot(cfg)

%%
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
cfg.out_dir = 'figure7/version2';

cfg.stm_typ = {'Language' 'Motor' ''};

cfg.fst_ltr = 1;

cfg.tsk     = 'SL';
cfg.eff_typ = 'hgp';
cfg.eff_nme = { 'pap_lng_950' 'pap_lng_950' 'pap_con_950'};
cfg.eff_clm = { 1             2             2 };
cfg.eff_col = { 'bright red'  'bright blue' 'cyan'};

cfg.chn_plt = 0;

cfg.dat_plt = 2;

cfg.alt_eve = 'trialinfo';
cfg.eve     = [1 2 3 4];
cfg.lnstyle.col_ord = {rgb('bright red') rgb('bright blue') rgb('reddish grey') rgb('bluish grey')} ;

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { 'vis_nse_stt_msk_anv'           'aud_nse_stt_msk_anv'          'vis_mtc_stt' };
cfg.stt_col = { { ft_stt_col(rgb('bright red')) ft_stt_col(rgb('bright blue')) ft_stt_col(rgb('dark yellow')) } };
cfg.stt_cmp = { { '0%3'                         '3%6'                          '6%9'   } };

cfg.v_lne       = [0 0.450 0.900];
cfg.v_lne_wdt   = [3 3 1];
cfg.v_lne_col   = {rgb('red') rgb('blue') rgb('black')};

cfg.xlm = [-0.2 1.3];

cfg.rmv_sbj = '';

cfg.inc_reg = { 'Inferior Precentral' 'Middle Precentral' ...
                'Supramarginal' ...
                'Caudal STG' 'Middle STG' ...
                'Pars Opercularis' };

mmil_stimulation_plot(cfg)

%% For Presentation
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
cfg.out_dir = 'figure7/version_presentation';

cfg.stm_typ = {'Language' 'Motor' ''};

cfg.fst_ltr = 1;

cfg.tsk     = 'SL';
cfg.eff_typ = 'hgp';
cfg.eff_nme = { 'pap_lng_950' 'pap_lng_950' };
cfg.eff_clm = { 1             2             };
cfg.eff_col = { 'bright red'  'bright blue' };

cfg.chn_plt = 0;

cfg.dat_plt = 2;

cfg.alt_eve = 'trialinfo';
cfg.eve     = [1 2 3 4];
cfg.lnstyle.col_ord = {rgb('bright red') rgb('bright blue') rgb('reddish grey') rgb('bluish grey')} ;

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { 'vis_nse_stt_msk_anv'           'aud_nse_stt_msk_anv'          'vis_mtc_stt' };
cfg.stt_col = { { ft_stt_col(rgb('bright red')) ft_stt_col(rgb('bright blue')) ft_stt_col(rgb('dark yellow')) } };
cfg.stt_cmp = { { '0%3'                         '3%6'                          '6%9'   } };

cfg.v_lne       = [0 0.450 0.900];
cfg.v_lne_wdt   = [3 3 1];
cfg.v_lne_col   = {rgb('red') rgb('blue') rgb('black')};

cfg.xlm = [-0.2 1.3];

cfg.rmv_sbj = '';

cfg.inc_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral'...
                'Supramarginal' ...
                'Caudal STG' 'Middle STG' 'Middle MTG' ...
                'Pars Opercularis' };

mmil_stimulation_plot(cfg)

%
[pot_typ ; sel_ele]

intersect(sel_ele{3},sel_ele{4})

%







