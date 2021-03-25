function fw_plot(fcfg)

%% Initial Load
sbj  = fcfg.sbj_nme;

fprintf([fcfg.sbj_nme ': Starting initial plots work on %s \n'],sbj)

infile = [fcfg.sbj_dat_hld '/' sbj '_overall_data.mat'];

cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.sbj_dat_hld '/' sbj '_overall_data.mat'];
fwv_dat  = ft_func([],cfg);

for iD = 1:numel(fwv_dat.data_name); 
    fwv_dat.(fwv_dat.data_name{iD}).cfg.alt_lab.stt_lab = fwv_dat.(fwv_dat.data_name{iD}).label; 
end

%% Channel-by-Channel, 2x3 (VisNse/AudNse/Mtc x LFP/HGP)
for iD = 1:numel(fwv_dat.data_name);
    
    fprintf('Starting work on %s \n',fcfg.sbj_nme)
    
    cfg           = [];
    
    cfg.type      = 'chan';
    cfg.lgd       = 0;
    
    cfg.dat       = {fwv_dat.(fwv_dat.data_name{iD})};
    
    cfg.plt_dim   = [3 3];
    cfg.dat_loc   = [ 1 1 1 ; ...
                      1 1 1 ; ...
                      1 1 1 ];
    
    cfg.plt_lbl   = { 'Stimuli'    'Letters'     'Letter/FF'   ; ...
                      'Overall'    'Novel Words' 'Word/Letter' ; ...
                      'False Font' 'Repetition'  'Old / New'   };
    cfg.ttl_num   = zeros(4,4);
    
    cfg.alt_eve   = {'trialinfo' 'trialinfo' 'trialinfo' ; ...
                     'vis_ovr'   'trialinfo' 'trialinfo' ; ...
                     'trialinfo' 'trialinfo' 'trialinfo' };
    cfg.eve       = {[3 4 5 6] [5] [5 6] ; ...
                     [101]     [3] [3 5] ; ...
                     [6]       [4] [3 4] };
    
    cfg.std_err   = 1;
    cfg.lnstyle.col_ord = {{rgb('bright red') rgb('orange') rgb('red') rgb('grey')} {rgb('dark red')}   {rgb('dark red') rgb('grey')}       ; ...
                           {rgb('black')}                                           {rgb('bright red')} {rgb('bright red') rgb('dark red')} ; ...
                           {rgb('grey')}                                            {rgb('orange')}     {rgb('bright red') rgb('orange')}   };
    cfg.cnd_nme         = {};
    
    cfg.stt_dat = { {'vis_stm'}     {'fw_ltr_stt'} {'vis_ltr'} ; ...
                    {'vis_ovr_stt'} {'fw_wrd_stt'} {'vis_wrd'} ; ...
                    {'fw_ffn_stt'}  {'fw_rep_stt'} {'vis_old'}            };
    cfg.stt_col = { {ft_stt_col(rgb('red'))}   {ft_stt_col(rgb('dark red'))}   {ft_stt_col(rgb('dark red'))}   ; ...
                    {ft_stt_col(rgb('black'))} {ft_stt_col(rgb('bright red'))} {ft_stt_col(rgb('bright red'))} ; ...
                    {ft_stt_col(rgb('grey'))}  {ft_stt_col(rgb('orange'))}     {ft_stt_col(rgb('orange'))}     };
    
    cfg.alt_lbl = 'label';
    cfg.stt_lab = 'stt_lab';
    
    cfg.y_lnk       = [ 1 1 1 ; 1 1 1 ; 1 1 1 ];
    cfg.x_lim       = [-0.15 0.75];
    cfg.v_lne       = { [0.000 0.200 0.400] [0.000 0.200 0.400] [0.000 0.200 0.400] ; ...
                        [0.000 0.200 0.400] [0.000 0.200 0.400] [0.000 0.200 0.400] ; ...
                        [0.000 0.200 0.400] [0.000 0.200 0.400] [0.000 0.200 0.400] };
    cfg.v_lne_col   = {{rgb('black') rgb('black') rgb('black')} {rgb('black') rgb('black') rgb('black')} {rgb('black') rgb('black') rgb('black')} ; ...
                       {rgb('black') rgb('black') rgb('black')} {rgb('black') rgb('black') rgb('black')} {rgb('black') rgb('black') rgb('black')} ; ...
                       {rgb('black') rgb('black') rgb('black')} {rgb('black') rgb('black') rgb('black')} {rgb('black') rgb('black') rgb('black')} };
    
    cfg.v_lne_wdt   = {[2 2 2] [2 2 2] [2 2 2] ; ...
                       [2 2 2] [2 2 2] [2 2 2] ; ...
                       [2 2 2] [2 2 2] [2 2 2] };
    cfg.axe_fnt_sze = [15 15 15 ; ...
                       15 15 15 ; ...
                       15 15 15 ];
    cfg.axe_lne_sze = [3 3 3 ; ...
                       3 3 3 ; ...
                       3 3 3 ];
    cfg.ttl_lne_sze = [36 36 36 ; ...
                       36 36 36 ; ...
                       36 36 36 ];
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = [fcfg.sbj_dat_hld '/' 'initial_plot' '/' fcfg.sbj_nme '/' fwv_dat.data_name{iD}];
    cfg.prefix     = [fcfg.sbj_nme '_'];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
end

end