function SL_effects_lfp(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Note Overall Effect
% cfg = [];
% cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.data_name = 1;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.alt_stt = { 'ovr_stt' };
% cfg.cmp_trl = { 'ovr_all_eve' };
% cfg.cmp_stt = { [ 1 1 ] };
% cfg.cmp_nme = { { 'vis_erl_500ms' 'aud_lte_500ms'} };  
% cfg.cmp     = { { '101!999' '101!999'} };
% cfg.tme_win = { { [0.025 0.500] [0.475 0.950] } };
% cfg.stt_col = { { 'reddish grey' 'bluish grey'} };
%            
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'lfp'}; end
%     cfg.specific = {'typ' ; 1:numel(cfg.typ)};
% 
% dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
% 
% % Split effects into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_stat(cfg);
% 
% % Split timing into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_timing(cfg);

%% Note difference effects
% cfg = [];
% cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.data_name = 1;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.alt_stt = { 'vis_anv_stt' ...
%                 'vis_anv_stt_01' };
% cfg.cmp_trl = { 'trialinfo' ...
%                 'trialinfo' };
% cfg.cmp_stt = { [ 1 1 ] ...
%                 [ 2 2 ] };
% cfg.cmp_nme = { { 'anv_500ms' 'anv_500ms'} ...
%                 { 'anv_01_500ms' 'anv_01_500ms'} };  
% cfg.cmp     = { { '1!2' '1!2' } ...
%                 { '1!2' '1!2'} };
% cfg.tme_win = { { [0.025 0.500] [0.025 0.500] } ...
%                 { [0.475 0.950] [0.475 0.950] } };
% cfg.stt_col = { { 'orange' 'dark orange' } ...
%                 { 'orange' 'dark orange' } };
%            
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'hgp'}; end
%     cfg.specific = {'typ' ; 1:numel(cfg.typ)};
% 
% dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
% 
% % Split effects into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_stat(cfg);
% 
% % Split timing into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_timing(cfg);

%% Note difference effects
% cfg = [];
% cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.data_name = 1;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.alt_stt = { 'ovr_stt' };
% cfg.cmp_trl = { 'ovr_all_eve' };
% cfg.cmp_stt = { [ 1 1 ] };
% cfg.cmp_nme = { { 'vis_erl_500ms' 'aud_lte_500ms'} };  
% cfg.cmp     = { { '101!999' '101!999'} };
% cfg.tme_win = { { [0.025 0.500] [0.475 0.950] } };
% cfg.stt_col = { { 'reddish grey' 'bluish grey'} };
%            
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'lfp'}; end
%     cfg.specific = {'typ' ; 1:numel(cfg.typ)};
% 
% dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
% 
% % Split effects into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_stat(cfg);
% 
% % Split timing into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_timing(cfg);


%% Linguistic Effects
% cfg = [];
% cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.data_name = 1;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.alt_stt = { 'vis_nse_stt_msk' ...
%                 'aud_nse_stt_msk' };
% cfg.cmp_trl = { 'vis_tot_nse' ...
%                 'aud_tot_nse' };
% cfg.cmp_stt = { [ 1 ] ...
%                 [ 2 ] };
% cfg.cmp_nme = { { 'ltr_ffl_500ms' } ...
%                 { 'vce_nse_500ms' } };  
% cfg.cmp     = { { '101!102' } ...
%                 { '111!112' } };
% cfg.tme_win = { { [0.025 0.500] } ...
%                 { [0.475 0.950] } };
% cfg.stt_col = { { 'bright red'  } ...
%                 { 'bright blue' } };
%            
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'lfp'}; end
%     cfg.specific = {'typ' ; 1:numel(cfg.typ)};
% 
% dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
% 
% % Split effects into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_stat(cfg);
% 
% % Split timing into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_timing(cfg);

%% Mismatch/Match Effects
% cfg = [];
% cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.data_name = 1;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.alt_stt = { 'vis_mtc_stt_msk' };
% cfg.cmp_trl = { 'trialinfo' };
% cfg.cmp_stt = { [ 1 1 ] };
% cfg.cmp_nme = { { 'sns_mtc' 'lte_mtc' } };  
% cfg.cmp     = { { '1!2' '1!2' } };
% cfg.tme_win = { { [0.475 0.900] [0.900 1.450] } };
% cfg.stt_col = { { 'green' 'yellow' } };
%            
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'lfp'}; end
%     cfg.specific = {'typ' ; 1:numel(cfg.typ)};
% 
% dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
% 
% % Split effects into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_stat(cfg);
% 
% % Split timing into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_timing(cfg);

%% Word/NonWord Effects
% cfg = [];
% cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.data_name = 1;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.alt_stt = { 'vis_wrd_stt_msk' ...
%                 'aud_wrd_stt_msk' };
% cfg.cmp_trl = { 'vis_wrd' ...
%                 'aud_wrd' };
% cfg.cmp_stt = { [ 1 ] ...
%                 [ 2 ] };
% cfg.cmp_nme = { { 'vis_wrd_nwd_sel' } ...
%                 { 'aud_wrd_nwd_sel' } };  
% cfg.cmp     = { { '101!102' } ...
%                 { '201!202' } };
% cfg.tme_win = { { [0.025 0.500] } ...
%                 { [0.475 0.950] } };
% cfg.stt_col = { { 'bright red'  } ...
%                 { 'bright blue' } };
%            
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% if any(strcmpi(cfg.typ,'hlb')); cfg.typ(strcmpi(cfg.typ,'hlb')) = {'hgp'}; end
%     cfg.specific = {'typ' ; 1:numel(cfg.typ)};
% 
% dat          = ft_func(@mmil_chn_eff,cfg,bcc_dat);
% 
% % Split effects into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_stat(cfg);
% 
% % Split timing into ecog/depth
% cfg = [];
% 
% cfg.sbj_nme = fcfg.sbj_nme;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.dat_typ = {'lfp'};
% 
% mmil_split_timing(cfg);

%% Phonemes
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.data_name = 1;

cfg.stt_lab = 'stt_lab';

cfg.alt_stt = { 'vis_con_stt_msk' ...
                'aud_con_stt_msk' };
cfg.cmp_trl = { 'vis_con' ...
                'aud_con' };
cfg.cmp_stt = { [ 1 1 ] ...
                [ 2 2 ] };
cfg.cmp_nme = { { 'ely_vis_con_sel' 'lte_vis_con_sel' } ...
                { 'ely_aud_con_sel' 'lte_aud_con_sel' } };  
cfg.cmp     = { { '1!12' '1!12' } ...
                { '1!12' '1!12' } };
cfg.tme_win = { { [0.025 0.250] [0.250 0.450] } ...
                { [0.475 0.750] [0.250 0.450] } };
cfg.stt_col = { { 'red'  'red'  } ...
                { 'blue' 'blue' } };
           
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