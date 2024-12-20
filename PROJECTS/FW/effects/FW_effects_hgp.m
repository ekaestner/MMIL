function FW_effects_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Responsive Electrodes
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 2;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_ovr_stt' ...
                'fw_ffn_stt' ...
                'fw_ltr_stt' ...
                'fw_wrd_stt' ...
                'fw_rep_stt' };                 
cfg.cmp_trl = { 'vis_ovr' ...
                'trialinfo' ...
                'trialinfo' ...
                'trialinfo' ...
                'trialinfo' };
cfg.cmp_stt = { [ 1 ] ...
                [ 2 ] ...
                [ 3 ] ...
                [ 4 ] ...
                [ 5 ] };
cfg.cmp_nme = { { 'ovr_600ms' }  ...
                { 'ffl_600ms' } ...
                { 'ltr_600ms' } ...
                { 'wrd_600ms' } ...
                { 'rep_600ms' } };
cfg.cmp     = { { '101>999' }  ...
                { '6>999' } ...
                { '5>999' } ...
                { '3>999' } ...
                { '4>999' } };
cfg.tme_win = { { [0.025 0.600] } ...
                { [0.025 0.600] } ...
                { [0.025 0.600] } ...
                { [0.025 0.600] } ...
                { [0.025 0.600] } };
cfg.stt_col = { { 'black' }  ...
                { 'black' } ...
                { 'black' } ...
                { 'black' } ...
                { 'black' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(2),'uni',0);
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


%% Note Overall Effect
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 2;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_stm' };
cfg.cmp_trl = { 'trialinfo' };
cfg.cmp_stt = { [ 1 ] };
cfg.cmp_nme = { { 'dff_600ms' } };  
cfg.cmp     = { { '1!4' } };
cfg.tme_win = { { [0.025 0.600] } };
cfg.stt_col = { { 'light red' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(2),'uni',0);
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

%% Note Overall Effect
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 2;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_stm_01' };
cfg.cmp_trl = { 'trialinfo' };
cfg.cmp_stt = { [ 1 ] };
cfg.cmp_nme = { { 'dff_01_600ms' } };  
cfg.cmp     = { { '1!4' } };
cfg.tme_win = { { [0.025 0.600] } };
cfg.stt_col = { { 'light red' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(2),'uni',0);
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

%% Note Lexical Decision Stats
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 2;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_ltr_msk' ... 
                'vis_wrd_msk' ...
                'vis_old_msk' };
cfg.cmp_trl = { 'trialinfo' ... 
                'trialinfo' ...
                'trialinfo' };
cfg.cmp_stt = { [1 1] ... 
                [2 2] ...
                [3 3] };
cfg.cmp_nme = { { 'vis_ltr_600ms' 'vis_ffn_600ms'  }  ; ... 
                { 'vis_wrd_600ms' 'vis_nwd_600ms'  } ; ...
                { 'vis_nov_600ms' 'vis_rep_600ms'  } };  
cfg.cmp     = { { '5>6' '5<6' }   ... 
                { '3>5' '3<5' } ...
                { '3>4' '3<4' } };
cfg.tme_win = { { [0.025 0.600] [0.025 0.600] } ... 
                { [0.025 0.600] [0.025 0.600] } ...
                { [0.025 0.600] [0.025 0.600] } };
cfg.stt_col = { { 'dark red'   'reddish grey' } ... 
                { 'bright red' 'dark red' }    ...
                { 'bright red' 'orange' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(2),'uni',0);
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