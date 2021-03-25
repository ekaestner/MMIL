function iSASZ_Figure7

%% VISUAL OVERALL HGP
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/visualresponsive';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'hgp';
cfg.eff_nme = { 'pap_vis_act' 'pap_vis_rep' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'dark red'  'bright red' };

cfg.chn_plt = 0;

cfg.dat_plt = 2;

cfg.alt_eve = 'vis_new_old';
cfg.eve     = [111 112];
cfg.lnstyle.col_ord = { rgb('red')  rgb('bright red') };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { 'vis_new_ovr_stt'                 'vis_new_old_stt_msk_ovr'          };
cfg.stt_col = { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('red')) } };
cfg.stt_cmp = { { '0%3'                           '3%6'                          } };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

%% AUDITORY OVERALL HGP
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/auditoryresponsive';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'hgp';
cfg.eff_nme = { 'pap_aud_act' 'pap_aud_rep' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'dark blue' 'bright blue' };

cfg.chn_plt = 0;

cfg.dat_plt = 2;

cfg.alt_eve = 'aud_new_old';
cfg.eve     = [211 212];
cfg.lnstyle.col_ord = { rgb('blue')  rgb('bluish purple') };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { 'aud_new_ovr_stt'                'aud_new_old_stt_msk_ovr'          };
cfg.stt_col = { { ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('blue')) } };
cfg.stt_cmp = { { '0%3'                          '3%6'                          } };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('blue') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

%% VISUAL OVERALL LFP
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/visualresponsive_lfp';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'lfp';
cfg.eff_nme = { 'pap_vis_act' 'pap_vis_rep' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'bright red'  'orange' };

cfg.chn_plt = 1;

cfg.dat_plt = 2;

cfg.alt_eve = 'vis_new_old';
cfg.eve     = [111 112];
cfg.lnstyle.col_ord = { rgb('red')  rgb('orangish red') };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { 'vis_new_ovr_stt'                 'vis_new_old_stt_msk_ovr'          };
cfg.stt_col = { { ft_stt_col(rgb('reddish gray')) ft_stt_col(rgb('red')) } };
cfg.stt_cmp = { { '0%3'                           '3%6'                          } };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

%% AUDITORY OVERALL LFP
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/auditoryresponsive_lfp';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'lfp';
cfg.eff_nme = { 'pap_aud_act' 'pap_aud_rep' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'bright blue' 'bright purple' };

cfg.chn_plt = 1;

cfg.dat_plt = 2;

cfg.alt_eve = 'aud_new_old';
cfg.eve     = [211 212];
cfg.lnstyle.col_ord = { rgb('blue')  rgb('bluish purple') };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { 'aud_new_ovr_stt'                'aud_new_old_stt_msk_ovr'          };
cfg.stt_col = { { ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('blue')) } };
cfg.stt_cmp = { { '0%3'                          '3%6'                          } };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('blue') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

%% VIUSAL/AUDITORY Overall
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/bimodalresponsive';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'hgp';
cfg.eff_nme = { 'pap_vis_act' 'pap_aud_act' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'bright red'  'bright blue' };

cfg.chn_plt = 0;

cfg.dat_plt = 2;

cfg.dat_loc = [1 1 ; 0 0];
cfg.plt_dim = [2 2];

cfg.alt_eve = { 'vis_new_old' 'aud_new_old' ; '' '' };
cfg.eve     = { [111 112] [211 212] ; [] [] };
cfg.lnstyle.col_ord = { {  rgb('red')  rgb('reddish purple')} { rgb('blue')  rgb('bluish purple')} ; {} {} };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { { 'vis_new_ovr_stt'                'vis_new_old_stt_msk_ovr' } {'aud_new_ovr_stt'               'aud_new_old_stt_msk_ovr'} ; {} {} };
cfg.stt_col = { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('red')) } {ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('blue')) } ; {} {} };
cfg.stt_cmp = { { '0%3'                          '3%6'                   } { '0%3'                          '3%6' } ; {} {} };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('black') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

%% VIUSAL/AUDITORY Overall LFP
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/bimodalresponsive_lfp';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'lfp';
cfg.eff_nme = { 'pap_vis_act' 'pap_aud_rep' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'bright red'  'bright blue' };

cfg.chn_plt = 1;

cfg.dat_plt = 2;

cfg.dat_loc = [1 1 ; 0 0];
cfg.plt_dim = [2 2];

cfg.alt_eve = { 'vis_new_old' 'aud_new_old' ; '' '' };
cfg.eve     = { [111 112] [211 212] ; [] [] };
cfg.lnstyle.col_ord = { {  rgb('red')  rgb('reddish purple')} { rgb('blue')  rgb('bluish purple')} ; {} {} };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { { 'vis_new_ovr_stt'                'vis_new_old_stt_msk_ovr' } {'aud_new_ovr_stt'               'aud_new_old_stt_msk_ovr'} ; {} {} };
cfg.stt_col = { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('red')) } {ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('blue')) } ; {} {} };
cfg.stt_cmp = { { '0%3'                          '3%6'                   } { '0%3'                          '3%6' } ; {} {} };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('black') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

%% VIUSAL/AUDITORY Repetition
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/bimodalrepetition';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'hgp';
cfg.eff_nme = { 'pap_vis_rep' 'pap_aud_rep' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'bright orange'  'bright purple' };

cfg.chn_plt = 1;

cfg.dat_plt = 2;

cfg.dat_loc = [1 1 ; 0 0];
cfg.plt_dim = [2 2];

cfg.alt_eve = { 'vis_new_old' 'aud_new_old' ; '' '' };
cfg.eve     = { [111 112] [211 212] ; [] [] };
cfg.lnstyle.col_ord = { {  rgb('red')  rgb('reddish purple')} { rgb('blue')  rgb('bluish purple')} ; {} {} };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { { 'vis_new_ovr_stt'                'vis_new_old_stt_msk_ovr' } {'aud_new_ovr_stt'               'aud_new_old_stt_msk_ovr'} ; {} {} };
cfg.stt_col = { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('red')) } {ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('blue')) } ; {} {} };
cfg.stt_cmp = { { '0%3'                          '3%6'                   } { '0%3'                          '3%6' } ; {} {} };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('black') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

%% VIUSAL/AUDITORY Repetition LFP
cfg = [];

cfg.chn_crr = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.out_dir = '/figure7/bimodalresponsive_lfp';

cfg.stm_typ = {'Language' 'Motor'};

cfg.fst_ltr = 1;

cfg.tsk     = 'SA_SZ';
cfg.eff_typ = 'lfp';
cfg.eff_nme = { 'pap_vis_rep' 'pap_aud_rep' };
cfg.eff_clm = { 1             1 };
cfg.eff_col = { 'bright orange'  'bright purple' };

cfg.chn_plt = 1;

cfg.dat_plt = 2;

cfg.dat_loc = [1 1 ; 0 0];
cfg.plt_dim = [2 2];

cfg.alt_eve = { 'vis_new_old' 'aud_new_old' ; '' '' };
cfg.eve     = { [111 112] [211 212] ; [] [] };
cfg.lnstyle.col_ord = { {  rgb('red')  rgb('reddish purple')} { rgb('blue')  rgb('bluish purple')} ; {} {} };

cfg.stt_lab = 'stt_lab';
cfg.stt_dat = { { 'vis_new_ovr_stt'                'vis_new_old_stt_msk_ovr' } {'aud_new_ovr_stt'               'aud_new_old_stt_msk_ovr'} ; {} {} };
cfg.stt_col = { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('red')) } {ft_stt_col(rgb('bluish gray')) ft_stt_col(rgb('blue')) } ; {} {} };
cfg.stt_cmp = { { '0%3'                          '3%6'                   } { '0%3'                          '3%6' } ; {} {} };

cfg.v_lne       = [0          0.200        0.400 ];
cfg.v_lne_wdt   = [3          1            1 ];
cfg.v_lne_col   = {rgb('black') rgb('black') rgb('black') };

cfg.xlm = [-0.200 0.750];

cfg.rmv_sbj = '';

mmil_stimulation_plot(cfg)

end



















