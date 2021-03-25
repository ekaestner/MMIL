function FW_mask(fcfg)

%% Initial Load
fprintf([fcfg.sbj_nme ': Starting initial stats mask work on %s \n'],fcfg.sbj_nme)

infile = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
outpath = fcfg.dat_fld;

cfg = [];
cfg.load = 'yes';
cfg.file = [outpath '/' fcfg.sbj_nme '_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

for iF = 1:numel(bcc_dat.data_name); 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.stt_lab = bcc_dat.(bcc_dat.data_name{iF}).label; 
    bcc_dat.(bcc_dat.data_name{iF}).cfg.alt_lab.label = bcc_dat.(bcc_dat.data_name{iF}).label;
end

%% Remove mask
for iD = 1:numel(bcc_dat.data_name)
    msk_nme = fieldnames(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_stt);
    msk_nme = msk_nme(string_find(msk_nme,'_msk'));
    
    for iR = 1:numel(msk_nme)
        bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_stt = rmfield(bcc_dat.(bcc_dat.data_name{iD}).cfg.alt_stt,msk_nme{iR});
    end    
end

%% Mask Stats
cfg     = [];
cfg.stt     = { 'vis_ltr' 'vis_wrd' 'vis_old' }; %   
cfg.stt_msk = { 'vis_stm' 'vis_stm' 'vis_stm' }; %
bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);

% cfg     = [];
% cfg.stt     = { 'vis_ltr'    'vis_wrd'    'vis_old' }; %   
% cfg.stt_msk = { 'vis_stm_01' 'vis_stm_01' 'vis_stm_01' }; %
% cfg.pst_fix = '_msk_01';
% bcc_dat = ft_func(@ft_mask_stats,cfg,bcc_dat);

%% Save Data & Stats
cfg = [];
cfg.str_nme  = 'bcc_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename =[outpath '/' fcfg.sbj_nme '_overall_data'];
ft_func([],cfg,bcc_dat);

end