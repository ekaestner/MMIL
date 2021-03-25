function iSASZ_effects_hgp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iD = 1:numel(bcc_dat.data_name)
    bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iD}).label;
end

%% Note Overall Effect
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_ovr_stt' ... 
                'aud_ovr_stt' };
cfg.cmp_trl = { 'vis_ovr' ... 
                'aud_ovr' };
cfg.cmp_stt = { [ 1 ] ... 
                [ 2 ] };
cfg.cmp_nme = { { 'vis_ovr_800ms' }  ; ... 
                { 'aud_ovr_800ms' } };  
cfg.cmp     = { { '101>999' }   ... 
                { '201>999' } };
cfg.tme_win = { { [0.025 0.800] } ... 
                { [0.025 0.800] } };
cfg.stt_col = { { 'bright red'   } ... 
                { 'bright blue' }   };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name,'uni',0);
cfg.specific = {'typ' ; 1:numel(bcc_dat.data_name)};

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

%% Note Semantic Effect
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_ani_obj_stt_msk' ... 
                'aud_ani_obj_stt_msk' };
cfg.cmp_trl = { 'vis_ani_obj' ... 
                'aud_ani_obj' };
cfg.cmp_stt = { [ 1 1 ] ... 
                [ 2 2 ] };
cfg.cmp_nme = { { 'vis_ani_800ms' 'vis_obj_800ms'}  ; ... 
                { 'aud_ani_800ms' 'vis_obj_800ms'} };  
cfg.cmp     = { { '121>122' '121<122' }   ... 
                { '221>222' '221<222'} };
cfg.tme_win = { { [0.025 0.800] [0.025 0.800] } ... 
                { [0.025 0.800] [0.025 0.800] } };
cfg.stt_col = { { 'bright green' 'grey' } ... 
                { 'bright green' 'grey' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name,'uni',0);
cfg.specific = {'typ' ; 1:numel(bcc_dat.data_name)};

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

%% Note Repetition Effect
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_new_old_stt_msk' ... 
                'aud_new_old_stt_msk' };
cfg.cmp_trl = { 'vis_new_old' ... 
                'aud_new_old' };
cfg.cmp_stt = { [ 1 1 ] ... 
                [ 2 2 ] };
cfg.cmp_nme = { { 'vis_nov_800ms' 'vis_rep_800ms'}  ; ... 
                { 'aud_nov_800ms' 'vis_rep_800ms'} };  
cfg.cmp     = { { '111>112' '111<112' }   ... 
                { '211>212' '211<212'} };
cfg.tme_win = { { [0.025 0.800] [0.025 0.800] } ... 
                { [0.025 0.800] [0.025 0.800] } };
cfg.stt_col = { { 'bright red'  'orange' } ... 
                { 'bright blue' 'orange' } };
           
cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name,'uni',0);
cfg.specific = {'typ' ; 1:numel(bcc_dat.data_name)};

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