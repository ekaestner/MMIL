function ejk_stt_mst(cfg)

%% Load
sbj  = cfg.sbj_nme;

fprintf([cfg.sbj_nme ': Starting initial stats work on %s \n'],sbj)

infile = [cfg.out_pth '/' sbj '_overall_data.mat'];
outpath = cfg.out_pth;

fcfg = [];
fcfg.load = 'yes';
fcfg.file = [outpath '/' sbj '_overall_data.mat'];
stt_dat  = ft_func([],fcfg);

%% Stat Labels
for iD = 1:numel(stt_dat.data_name); 
    stt_dat.(stt_dat.data_name{iD}).cfg.alt_lab.stt_lab = stt_dat.(stt_dat.data_name{iD}).label; 
end

%% Stats
fcfg = [];
fcfg.stt_fnc  = {'mst_quk_lur_nov_stt'};
fcfg.chn      = [ 93 126 ];
fcfg.loc      = 'local';
fcfg.fld_nme  = stt_dat.data_name;
fcfg.specific = {'fld_nme' ; 1:numel(fcfg.fld_nme)};
stt_dat      = ft_func(@mmil_cloud_stat,fcfg,stt_dat);

fcfg = [];
fcfg.stt_fnc  = {'mst_quk_lur_rep_stt'};
fcfg.chn      = [ 93 126 ];
fcfg.loc      = 'local';
fcfg.fld_nme  = stt_dat.data_name;
fcfg.specific = {'fld_nme' ; 1:numel(fcfg.fld_nme)};
stt_dat      = ft_func(@mmil_cloud_stat,fcfg,stt_dat);

%% Plots
% LFP
fcfg = [];

fcfg.dat       = {stt_dat.(stt_dat.data_name{1})};
fcfg.alt_eve   = 'over_event';
fcfg.plt_dim   = [1 1];
fcfg.type      = 'chan';

eve_plt     = unique(stt_dat.(stt_dat.data_name{1}).cfg.alt_eve.(fcfg.alt_eve));
eve_plt(eve_plt==0) = [];

col_hld = distinguishable_colors(numel(eve_plt));
for iC = 1:size(col_hld,1); 
    lne_col_plt{iC} = col_hld(iC,:); 
    cnd_nme_plt{iC} = ['eve' num2str(iC)];
end
    
leg_pos_plt = 1:numel(cnd_nme_plt);

v_lne_plt     = [0];
v_lne_col_plt = {rgb('black')};
v_lne_wdt_plt = [1];

fcfg.eve       = eve_plt';

fcfg.alt_lbl = 'label';
fcfg.stt_lab = 'stt_lab';

fcfg.stt_dat = { {'mst_quk_lur_nov_stt' 'mst_quk_lur_rep_stt'} };
fcfg.stt_col = { {ft_stt_col(rgb('blue')) ft_stt_col(rgb('red'))} };
fcfg.stt_cmp = { { '0%4' '4%8' } };
  
fcfg.std_err         = 1;
fcfg.lnstyle.col_ord = lne_col_plt;
fcfg.cnd_nme         = cnd_nme_plt;
fcfg.leg_pos         = leg_pos_plt;

fcfg.x_lim       = [-0.300 2];
fcfg.y_lim       = 'maxmin';
fcfg.v_lne       = v_lne_plt;
fcfg.v_lne_col   = v_lne_col_plt;
fcfg.v_lne_wdt   = v_lne_wdt_plt;
fcfg.axe_fnt_sze = 10;
fcfg.axe_lne_sze = 1.5;
fcfg.ttl_lne_sze = 20;
fcfg.ttl_num     = 0;

fcfg.print      = 1;
fcfg.nofig      = 1;
fcfg.print_type = 'png';
fcfg.outdir     = [cfg.prj_dat_hld '/' 'epoch_data' '/' 'plot_quick_stats' '/' cfg.sbj_nme '/' 'channels' '/' stt_dat.data_name{1}];
fcfg.prefix     = [stt_dat.data_name{1}] ;

mmil_ieeg_sensor_plot_v5(fcfg)

%% Save Data & Stats
fcfg = [];
fcfg.str_nme  = 'stt_dat';
fcfg.save     = 'yes';
fcfg.sve_app  = 'app_all';
fcfg.filename =[outpath '/' sbj '_overall_data.mat'];
ft_func([],fcfg,stt_dat);


end