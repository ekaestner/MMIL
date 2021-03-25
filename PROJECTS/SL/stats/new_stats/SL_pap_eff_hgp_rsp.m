function SL_pap_eff_hgp_rsp(fcfg)

fprintf([fcfg.sbj_nme ': Starting SL_pap_eff_hgp_lng on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
    
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_stt.pap_ovr_stt = bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_stt.ovr_stt;
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_stt.pap_anv_stt = bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_stt.vis_anv_stt;
end

%% Note Lexical Decision Stats
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 2;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'pap_ovr_stt' ... 
                'pap_anv_stt' };
cfg.cmp_trl = { 'ovr_all_eve' ... 
                'trialinfo' };
cfg.cmp_stt = { [ 1 1 1 1 ] ... 
                [ 2 2 2 2 ] };
cfg.cmp_nme = { { 'anv_500ms' 'anv_1000ms' 'anv_1500ms' 'anv_ovr' } ; ... 
                { 'ovr_500ms' 'ovr_1000ms' 'ovr_1500ms' 'anv_ovr' } };  
cfg.cmp     = { { '101>999'   '101>999'    '101>999'    '101>999' } ... 
                { '1!2'       '1!2'        '1!2'        '1!2' } };
cfg.tme_win = { { [0.025 0.500] [0.475 0.950] [0.975 1.500] [0.025 1.500] } ... 
                { [0.025 0.500] [0.475 0.950] [0.975 1.500] [0.025 1.500] } };
cfg.stt_col = { { 'red'         'blue'        'green'       'black' } ... 
                { 'red'         'blue'        'green'       'black' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'hgp'}; end
cfg.specific = {'typ' ; 1:numel(cfg.typ)};

dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);

% Split effects into ecog/depth
cfg = [];

cfg.sbj_nme = fcfg.sbj_nme;
cfg.clr_fld = fcfg.clr_fld;

cfg.dat_typ = {'hgp'};

mmil_split_stat(cfg);

% Split timing into ecog/depth
cfg = [];

cfg.sbj_nme = fcfg.sbj_nme;
cfg.clr_fld = fcfg.clr_fld;

cfg.dat_typ = {'hgp'};

mmil_split_timing(cfg);

end