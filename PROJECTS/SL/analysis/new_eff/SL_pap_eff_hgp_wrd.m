function SL_pap_eff_hgp_wrd(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%%
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 2;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_wrd_stt_msk_anv' ...
                'aud_wrd_stt_msk_anv' ...
                'vis_con_stt' };
cfg.cmp_trl = { 'vis_wrd' ...
                'aud_wrd' };
cfg.cmp_stt = { [ 1 1 ] ...
                [ 2 2 ] };
cfg.cmp_nme = { { 'vis_nwd' 'vis_wrd' } ...
                { 'aud_nwd' 'aud_wrd' } };  
cfg.cmp     = { { '101>102' '101<102' } ...
                { '201>202' '201<202'} };
cfg.tme_win = { { [0.025 0.450] [0.025 0.450] } ...
                { [0.475 0.900] [0.475 0.900] } };
cfg.stt_col = { { 'dark red'    'bright red'  } ...
                { 'dark blue'   'bright blue' } };
           
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