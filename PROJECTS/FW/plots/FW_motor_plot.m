function FW_motor_plot(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial plots work on %s \n'],sbj)

infile = [fcfg.out_pth '/' sbj '_overall_data.mat'];

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.out_pth '/' sbj '_overall_data.mat'];
fwv_dat  = ft_func([],cfg);

for iD = 1:numel(fwv_dat.data_name); 
    fwv_dat.(fwv_dat.data_name{iD}).cfg.alt_lab.stt_lab = fwv_dat.(fwv_dat.data_name{iD}).label; 
end

%% Channel-by-Channel, 2x3 (VisNse/AudNse/Mtc x LFP/HGP)
for iD = 2;
    
    fprintf('Starting work on %s \n',fcfg.sbj_nme)
    
    cfg           = [];
    
    cfg.type      = 'chan';
    cfg.lgd       = 0;
    
    cfg.dat       = {fwv_dat.(fwv_dat.data_name{iD})};
    
    cfg.plt_dim   = [ 1 1 ];
    cfg.dat_loc   = [ 1 ];
    
    cfg.plt_lbl   = { 'Stimuli' };
    cfg.ttl_num   = zeros(1,1);
    
    cfg.alt_eve   = { 'trialinfo' };
    cfg.eve       = { [3 6 7] };
    
    cfg.std_err   = 1;
    cfg.lnstyle.col_ord = { { rgb('bright red') rgb('reddish grey') rgb('black') } };
    cfg.cnd_nme         = {};
    
    cfg.stt_dat = { { 'vis_mtr_wrd_stt' 'vis_mtr_trg_stt' } };
    cfg.stt_col = { { ft_stt_col(rgb('reddish grey')) ft_stt_col(rgb('black')) }  };
    cfg.stt_cmp = { { '3vs6' '3vs7' } };
    
    cfg.alt_lbl = 'label';
    cfg.stt_lab = 'stt_lab';
    
    cfg.y_lnk       = [ 1 ];
    cfg.x_lim       = [-0.15 0.75];
    cfg.v_lne       = { [0.000 0.200 0.400] };
    cfg.v_lne_col   = { {rgb('black') rgb('black') rgb('black')} };
    
    cfg.v_lne_wdt   = { [2 2 2] };
    cfg.axe_fnt_sze = [ 15 ];
    cfg.axe_lne_sze = [ 3 ];
    cfg.ttl_lne_sze = [ 36 ];
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data' '/' 'motor_plot' '/' fcfg.sbj_nme '/'];
    cfg.prefix     = [fcfg.sbj_nme '_'];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
end

end