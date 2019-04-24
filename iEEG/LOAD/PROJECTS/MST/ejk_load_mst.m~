function ejk_load_mst(fcfg)

%% Initial Load
if ~exist(fcfg.out_pth); mkdir(fcfg.out_pth); end

%% Load
if ~exist([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/' fcfg.sbj_nme '_downsample.mat'])  
    end_dir     = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_dir');
    if strcmpi(end_dir,'.mat')
    [out_dat,trg_chn] = ejk_burke_raw_load(fcfg);
    end
else   
    load([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/' fcfg.sbj_nme '_downsample_micro.mat'])  
end

%% Trial Function
trl_fun     = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialfun'); %trl_fun = mmil_readtext([fcfg.clr_fld '/trialfun/' subj]); end
trialf_add  = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'trialf_add');
ignore      = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'ignore');
epc_tme     = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};
tsk         = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'tsk'); tsk = tsk{1};

if ~exist([fcfg.prj_dat_hld '/' 'clerical' '/' 'trialfun_output' '/' fcfg.sbj_nme '.mat'])
    
    cfg = [];
    cfg.sbj_nme  = fcfg.sbj_nme;
    cfg.ignore   = ignore;
    cfg.trialfun = trl_fun{1};
    cfg.prj_dat_hld  = fcfg.prj_dat_hld;
    cfg.trialfun_add = trialf_add;
    cfg.dataset      = {trg_chn};
    cfg.minduration  = 0.500;
    cfg.pre          = epc_tme(1)-1;
    cfg.pos          = epc_tme(2)+1;
    cfg.Fs          = out_dat.(out_dat.data_name{1}).fsample;
    cfg.time        = out_dat.(out_dat.data_name{1}).time{1};
    trl = ejk_trialfun(cfg);
        
    save([fcfg.prj_dat_hld '/' 'clerical' '/' 'trialfun_output' '/' fcfg.sbj_nme '.mat'],'trl');
    
else
    load([fcfg.prj_dat_hld '/' 'clerical' '/' 'trialfun_output' '/' fcfg.sbj_nme '.mat'],'trl');
end

cfg     = [];
cfg.specific  = {'trl';1:numel(trl)};
cfg.trl = trl;
ped_dat = ft_func(@ft_redefinetrial,cfg,out_dat);

clear out_dat

%% Adjust Events
ped_dat.(ped_dat.data_name{1}).cfg.alt_eve.trialinfo = ped_dat.(ped_dat.data_name{1}).trialinfo;

cfg = [];
cfg.return_events = 0;
cfg.old_events  = { [1:32 201:232] [101:132] [301:332] };
cfg.new_events  = {  1              2         3 };
cfg.crt_alt_eve = 'over_event';
cfg.use_alt_eve = 'trialinfo';
ped_dat = ft_func(@ejk_redefine_events,cfg,ped_dat);

%% Examine Data for Noise - Downsampled
beg_plt_spc = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'beg_plt_spc'); beg_plt_spc = beg_plt_spc{1};

if 1
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [fcfg.prj_dat_hld '/' 'epoch_data' '/' 'spectrum' '/' fcfg.sbj_nme '/' 'downsample_spectrum'];
    cfg.prefix = ['pedot_' fcfg.sbj_nme];
    ft_func(@ft_plot_spectrum,cfg,ped_dat);
end

%% Remove identified noise problems
bse_frq     = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'bse_frq'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';

if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    ped_dat      = ft_func(@ft_preprocessing,cfg,ped_dat);
end

%% Examine Data for Noise
end_plt_spc = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'end_plt_spc'); end_plt_spc = end_plt_spc{1};

if end_plt_spc
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [fcfg.prj_dat_hld '/' 'epoch_data' '/' 'spectrum' '/' fcfg.sbj_nme '/' 'After'];
    cfg.prefix = ped_dat.data_name;
    cfg.specific  = {'prefix'; 1:numel(ped_dat.data_name)};
    ft_func(@ft_plot_spectrum,cfg,ped_dat);
end

%% Filter Data
num_dta = 1:numel(ped_dat.data_name);

% Create LFP
cfg            = [];
cfg.data_name  = num_dta;
cfg.data_new   = 'yes';
cfg.lpfilter   = 'yes';
cfg.lpfreq     = 15;
cfg.new_suffix = 'lfp';
ped_dat        = ft_func(@ft_preprocessing,cfg,ped_dat);

% Create HGP Data
cfg            = [];
cfg.data_name  = num_dta;
cfg.hilbert    = 'amp';
cfg.freq_band  = {[70 160]};
cfg.data_new   = 'yes';
cfg.new_suffix = 'hgp_smt';
ped_dat        = ft_func(@ft_hilbert_freq_analsysis,cfg,ped_dat);

cfg             = [];
cfg.data_name   = (num_dta(end)*2)+1:num_dta(end)*3;
HGP_smooth_msec = 0.05;
w               = round(HGP_smooth_msec* ped_dat.(ped_dat.data_name{1}).fsample); 
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,w/2);
cfg.window      = gauss_y/sum(gauss_y);
ped_dat         = ft_func(@ft_window_data,cfg,ped_dat);

%% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = num_dta;
cfg.rmfield   = 'yes';
ped_dat       = ft_func([],cfg,ped_dat);

%% Baseline Data
cfg = [];
cfg.baselinewindow = [-0.300 0.000];
ped_dat            = ft_func(@ft_baseline,cfg,ped_dat);

%% Remove Artifacts
rjt_plt = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'rjt_plt'); rjt_plt = rjt_plt{1};%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);

cfg = [];
cfg.measures  = {'time' 'variance'};
cfg.thresh    = [0.95 0.95];
cfg.outdir    = [fcfg.out_pth '/' 'rejection_95' '/' fcfg.sbj_nme '/' ];
cfg.prefix    = ped_dat.data_name;
cfg.specific  = {'prefix';[1 2 3 4]};
cfg.pad       = 1;
cfg.plot      = 1;
ped_dat       = ft_func(@auto_rej,cfg,ped_dat);

cfg = [];
cfg.measure   = 'all';
cfg.apply     = 'ieeg';
ped_dat       = ft_func(@ft_apply_rejection,cfg,ped_dat);

%% Remove Padding
epc_tme = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'epc_tme'); epc_tme = epc_tme{1};

cfg = [];
cfg.latency = epc_tme;
ped_dat     = ft_func(@ft_selectdata,cfg,ped_dat);

%% Combine Data
cmb = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmb'); cmb = cmb{1};

if numel(cmb) > 1
    cfg = [];
    cfg.cmb = cmb;
    cfg.clr_fld = [fcfg.prj_dat_hld '/' 'clerical'];
    cfg.sbj_nme = fcfg.sbj_nme;
    ped_dat = mmil_combine_data2(cfg,ped_dat);
end

%% Preliminary Plots
ped_dat.(ped_dat.data_name{1}).cfg.alt_lab.label = ped_dat.(ped_dat.data_name{1}).label;
ped_dat.(ped_dat.data_name{2}).cfg.alt_lab.label = ped_dat.(ped_dat.data_name{2}).label;

% LFP
cfg = [];

cfg.dat       = {ped_dat.(ped_dat.data_name{1})};
cfg.alt_eve   = 'over_event';
cfg.plt_dim   = [1 1];
cfg.type      = 'chan';

eve_plt     = unique(ped_dat.(ped_dat.data_name{1}).cfg.alt_eve.(cfg.alt_eve));
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

cfg.eve       = eve_plt';

cfg.std_err         = 1;
cfg.lnstyle.col_ord = lne_col_plt;
cfg.cnd_nme         = cnd_nme_plt;
cfg.leg_pos         = leg_pos_plt;

cfg.x_lim       = [-0.300 2];
cfg.y_lim       = 'maxmin';
cfg.v_lne       = v_lne_plt;
cfg.v_lne_col   = v_lne_col_plt;
cfg.v_lne_wdt   = v_lne_wdt_plt;
cfg.axe_fnt_sze = 10;
cfg.axe_lne_sze = 1.5;
cfg.ttl_lne_sze = 20;
cfg.ttl_num     = 0;

cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'png';
cfg.outdir     = [fcfg.prj_dat_hld '/' 'epoch_data' '/' 'plot_explore' '/' fcfg.sbj_nme '/' 'channels' '/' ped_dat.data_name{1}];
cfg.prefix     = [ped_dat.data_name{1}] ;

mmil_ieeg_sensor_plot_v5(cfg)

% HGP
cfg.dat       = {ped_dat.(ped_dat.data_name{2})};
eve_plt     = unique(ped_dat.(ped_dat.data_name{2}).trialinfo);
cfg.outdir     = [fcfg.prj_dat_hld '/' 'epoch_data' '/' 'plot_explore' '/' fcfg.sbj_nme '/' 'channels' '/' ped_dat.data_name{2}];
cfg.prefix     = [ped_dat.data_name{2}] ;
mmil_ieeg_sensor_plot_v5(cfg)

%% Save
cfg = [];
cfg.str_nme  = 'ped_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.out_pth '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,ped_dat);

end