function SL_pap_eff_hgp_mtc(fcfg)

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

cfg.alt_stt = { 'vis_con_stt_msk_anv' ...
                'aud_con_stt_msk_anv' ...
                'vis_con_stt' ...
                'aud_con_stt' };
cfg.cmp_trl = { 'vis_con' ...
                'aud_con' ...
                'vis_con' ...
                'aud_con' };
cfg.cmp_stt = { [ 1 1 1 ] ...
                [ 2 2 2 ] ...
                [ 3 3 3 ] ...
                [ 4 4 4 ] };
cfg.cmp_nme = { { 'ely_vis_con_sel' 'lte_vis_con_sel' 'tot_vis_con_sel' } ...
                { 'ely_aud_con_sel' 'lte_aud_con_sel' 'tot_aud_con_sel' } ...
                { 'ely_vis_con_sel' 'lte_vis_con_sel' 'tot_vis_con_sel' } ...
                { 'ely_aud_con_sel' 'lte_aud_con_sel' 'tot_aud_con_sel' }};  
cfg.cmp     = { { '1!12'            '1!12'            '1!12' } ...
                { '1!12'            '1!12'            '1!12' } ...
                { '1!12'            '1!12'            '1!12' } ...
                { '1!12'            '1!12'            '1!12' } };
cfg.tme_win = { { [0.025 0.250] [0.250 0.450] [0.025 0.450] } ...
                { [0.475 0.700] [0.650 0.900] [0.475 0.900] } ...
                { [0.025 0.250] [0.250 0.450] [0.025 0.450] } ...
                { [0.475 0.700] [0.650 0.900] [0.475 0.900] } };
cfg.stt_col = { { 'red'         'red'         'red' } ...
                { 'blue'        'blue'        'blue' } ...
                { 'red'         'red'         'red' } ...
                { 'blue'        'blue'        'blue' } };
           
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