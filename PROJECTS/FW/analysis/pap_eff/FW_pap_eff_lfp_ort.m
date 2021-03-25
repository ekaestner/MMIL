function FW_pap_eff_hgp_ort(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end


cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 1;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_ort_tot' ...
                'vis_ort_nwd' ...
                'vis_ort_wrd' };
cfg.cmp_trl = { 'tot_ngh'   ...
                'tot_ngh' ...
                'tot_ngh' };
cfg.cmp_stt = { [ 1 ] ...
                [ 2 ] ...
                [ 3 ] };
cfg.cmp_nme = { { 'vis_ort_600ms'     } ...
                { 'vis_ort_nwd_600ms' } ...
                { 'vis_ort_wrd_600ms' } };  
cfg.cmp     = { { '401!403'     } ...
                { '402!403' } ...
                { '401!402' } };
cfg.tme_win = { { [0.025 0.600] } ...
                { [0.025 0.600] } ...
                { [0.025 0.600] } };
cfg.stt_col = { { 'bright red' } ...
                { 'bright red' } ...
                { 'bright red' } };
           
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