function FW_effectsize_update(fcfg)

fprintf([fcfg.sbj_nme ': Starting initial effects work on %s \n'],fcfg.sbj_nme)

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Specific Effectsize
cfg = [];
cfg.ovr_wrt = 1;
cfg.clr_fld = fcfg.clr_fld;

cfg.stt_lab = 'stt_lab';

cfg.data_name = 2;

cfg.alt_stt = { 'vis_stm_01' };
cfg.cmp_trl = { 'trialinfo' };
cfg.cmp_stt = { [ 1 ] };
cfg.cmp_nme = { { 'specific_600ms' } };  
cfg.cmp     = { { '3!6' } };
cfg.tme_win = { { [0.025 1.000] } };
cfg.stt_col = { { 'dark red' } };

cfg.sbj_nme  = fcfg.sbj_nme;
cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
cfg.specific = {'typ' ; 1:numel(cfg.data_name)};

dat          = ft_func(@mmil_chn_eff_sze,cfg,bcc_dat);

typ     = cfg.typ;
alt_stt = cfg.alt_stt;
clear fle_chk

%% SELECTIVE EFFECTSIZE
% cfg = [];
%  cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.data_name = 2;
% 
% cfg.alt_stt = { 'vis_stm_01' ... 
%                 'vis_wrd_ffn_msk_01' ...
%                 'vis_wrd_msk_01' };
% cfg.cmp_trl = { 'trialinfo' ... 
%                 'trialinfo' ...
%                 'trialinfo' };
% cfg.cmp_stt = { [ 1 ] ... 
%                 [ 2 ] ...
%                 [ 3 ] };
% cfg.cmp_nme = { { 'stm_600ms' }  ; ... 
%                 { 'vis_ltr_600ms' } ; ...
%                 { 'vis_wrd_600ms' } };  
% cfg.cmp     = { { '3!6' } ... 
%                 { '3!6' } ...
%                 { '3!5' } };
% cfg.tme_win = { { [0.025 1.000] } ... 
%                 { [0.025 1.000] } ...
%                 { [0.025 1.000] } };
% cfg.stt_col = { { 'dark red' } ... 
%                 { 'bright red' }    ...
%                 { 'orange' } };
% 
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% cfg.specific = {'typ' ; 1:numel(cfg.data_name)};
% 
% dat          = ft_func(@mmil_chn_eff_sze,cfg,bcc_dat);
% 
% typ     = cfg.typ;
% alt_stt = cfg.alt_stt;
% clear fle_chk

%% Lexical
% cfg = [];
% cfg.ovr_wrt = 1;
% cfg.clr_fld = fcfg.clr_fld;
% 
% cfg.stt_lab = 'stt_lab';
% 
% cfg.data_name = 2;
% 
% cfg.alt_stt = { 'vis_old' ... 
%                 'vis_wrd_frq' };
% cfg.cmp_trl = { 'trialinfo' ... 
%                 'wrd_frq_med' };
% cfg.cmp_stt = { [ 1 ] ... 
%                 [ 2 ] };
% cfg.cmp_nme = { { 'new_old_600ms' }  ; ... 
%                 { 'lex_frq_600ms' } };  
% cfg.cmp     = { { '3!4' } ... 
%                 { '121!122' } };
% cfg.tme_win = { { [0.025 1.000] } ... 
%                 { [0.025 1.000] } };
% cfg.stt_col = { { 'dark red' }   ... 
%                 { 'bright red' } };
% 
% cfg.sbj_nme  = fcfg.sbj_nme;
% cfg.typ      = cellfun(@(x) x(end-2:end),bcc_dat.data_name(cfg.data_name),'uni',0);
% cfg.specific = {'typ' ; 1:numel(cfg.data_name)};
% 
% dat          = ft_func(@mmil_chn_eff_sze,cfg,bcc_dat);
% 
% typ     = cfg.typ;
% alt_stt = cfg.alt_stt;
% clear fle_chk

end