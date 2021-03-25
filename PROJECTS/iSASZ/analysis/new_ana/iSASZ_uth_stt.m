clear; clc;

% Load Data
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';

fcfg.sbj_nme = 'MG49_SA_SZ_utah';

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
ecg_dat  = ft_func([],cfg);

%% Average Response
% # of channels included
96 - numel( [2 65 33 67 35 22 48 71 44 75 76 51 57 89 94 95])

% time of significant differences
% uth_dat.(uth_dat.data_name{1}).cfg.alt_stt.aud_new_ovr_stt_avg.time{1}(find(uth_dat.(uth_dat.data_name{1}).cfg.alt_stt.aud_new_ovr_stt_avg.prob<0.001,5))
% uth_dat.(uth_dat.data_name{1}).cfg.alt_stt.vis_new_ovr_stt_avg.time{1}(find(uth_dat.(uth_dat.data_name{1}).cfg.alt_stt.vis_new_ovr_stt_avg.prob<0.001,1))

%% Individual Responses
% Total number of 4 electrodes
vis_ovr_eff = mmil_readtext([ fcfg.clr_fld '/' 'sig_chn' '/' 'MG49_SA_SZ_utah' '/' 'cmb' '/' 'effect' '/' 'vis_new_ovr_stt' ]);
aud_ovr_eff = mmil_readtext([ fcfg.clr_fld '/' 'sig_chn' '/' 'MG49_SA_SZ_utah' '/' 'cmb' '/' 'effect' '/' 'aud_new_ovr_stt' ]);

eff_typ{1} = find(cell2mat(vis_ovr_eff(2:end,3)));
eff_typ{2} = find(cell2mat(aud_ovr_eff(2:end,3)));
eff_typ{3} = find(cell2mat(vis_ovr_eff(2:end,3)) & cell2mat(aud_ovr_eff(2:end,3)));
eff_typ{4} = find(~cell2mat(vis_ovr_eff(2:end,3)) & ~cell2mat(aud_ovr_eff(2:end,3)));

%% Amplitude
% Max
max(vis_new_avg) / max(aud_new_avg)


% Correlation
chn_inc = unique([eff_typ{1}' eff_typ{2}' eff_typ{3}']);

aud_amp = cell2mat(mdl_ele(chn_inc,3));
vis_amp = cell2mat(mdl_ele(chn_inc,4));

[cor_val,cor_pvl] = corrcoef([aud_amp vis_amp]);

% [max_val,max_loc] = max(aud_amp);
% max_val / vis_amp(max_loc)
% 
% [max_val,max_loc] = max(vis_amp);
% max_val / aud_amp(max_loc)

max(aud_amp./vis_amp)

max(vis_amp./aud_amp)




