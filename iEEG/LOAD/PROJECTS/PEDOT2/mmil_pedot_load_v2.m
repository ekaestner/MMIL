function mmil_pedot_load_v2(fcfg)

%% Initial Load
if ~exist(fcfg.out_pth); mkdir(fcfg.out_pth); end

%% Load
if ~exist([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/' fcfg.sbj_nme '_downsample.mat'])  
    [out_dat,trg_chn] = mmil_pedot_raw_load(fcfg);    
else   
    load([fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/' fcfg.sbj_nme '_downsample.mat'])  
end

%% Electrode Layout - %%% EJK %%%
% [gdd_grd_lab,bdd_grd_lab,gdd_tal_lab,hld_dat] = mmil_pedot_chn(fcfg,hld_dat);

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
    cfg.dsfact       = out_dat.(out_dat.data_name{1}).dsfact;
    cfg.dataset      = trg_chn;
    cfg.minduration  = 0.500;
    cfg.pre          = epc_tme(1)-1;
    cfg.pos          = epc_tme(2)+1;
    cfg.Fs          = out_dat.(out_dat.data_name{1}).fsample*out_dat.(out_dat.data_name{1}).dsfact;
    cfg.time        = out_dat.(out_dat.data_name{1}).time{1};
    trl = mmil_pedot_trialfun(cfg);
        
    save([fcfg.prj_dat_hld '/' 'clerical' '/' 'trialfun_output' '/' fcfg.sbj_nme '.mat'],'trl');
    
else
    load([fcfg.prj_dat_hld '/' 'clerical' '/' 'trialfun_output' '/' fcfg.sbj_nme '.mat'],'trl');
end

cfg     = [];
cfg.specific  = {'trl';1:numel(trl)};
cfg.trl = trl;
ped_dat = ft_func(@ft_redefinetrial,cfg,out_dat);

clear out_dat

cfg = [];
cfg.str_nme  = 'ped_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.out_pth '/' fcfg.sbj_nme '_overall_data_epoched_only.mat'];
ft_func([],cfg,ped_dat);

%% Examine Data for Noise - Downsampled
beg_plt_spc = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'beg_plt_spc'); beg_plt_spc = beg_plt_spc{1};

if 0
    cfg = [];
    cfg.empty  = 'yes';
    cfg.outdir = [fcfg.prj_dat_hld '/' 'epoch_data' '/' 'spectrum' '/' fcfg.sbj_nme '/' 'downsample_spectrum'];
    cfg.prefix = ['pedot_' fcfg.sbj_nme];
    ft_func(@ft_plot_spectrum,cfg,ped_dat);
end

%% Remove noise through removing common noise
chn_nse_grp     = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
chn_nse_grp_nme = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'cmn_nme');
nse_val         = mmil_load_subj_info([fcfg.prj_dat_hld '/' 'clerical' '/' 'sbj_inf' '/' fcfg.sbj_nme],'nse_val'); nse_val = nse_val{1};

if nse_val
    
    if strcmpi(chn_nse_grp{1},'split') | strcmpi(chn_nse_grp{1},'pedot'); chn_nse_grp = repmat(chn_nse_grp,1,numel(ped_dat.data_name)); end
    
    for iEJ = 1:numel(chn_nse_grp_nme);
        if strcmpi(chn_nse_grp{1},'split')
            for iD = 1:numel(ped_dat.data_name); pre_fix{iD}{1} = strcat(ped_dat.data_name{iD},'ecog'); pre_fix{iD}{2} = strcat(ped_dat.data_name{iD},'noisy_ecog'); pre_fix{iD}{3} = strcat(ped_dat.data_name{iD},'depth'); pre_fix{iD}{4} = strcat(ped_dat.data_name{iD},'noisy_depth'); end;
        elseif strcmpi(chn_nse_grp{1},'pedot')
            for iD = 1:numel(ped_dat.data_name); pre_fix{iD}{1} = strcat(ped_dat.data_name{iD},'grid'); pre_fix{iD}{2} = strcat(ped_dat.data_name{iD},'noisy_grid'); end;
        else
            for iCN = 1:numel(chn_nse_grp_nme{iEJ});
                pre_fix{iEJ}{iCN} = strcat(ped_dat.data_name{iEJ},'_',chn_nse_grp_nme{iEJ}{iCN});
            end;
        end
    end
    
    cfg             = [];
    cfg.sbj_nme     = fcfg.sbj_nme;
    cfg.clr_fld     = fcfg.clr_fld;
    cfg.pre_fix     = pre_fix;
    cfg.dat_nme     = ped_dat.data_name;
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix     = pre_fix;
    cfg.specific    = {'pre_fix' 'chn_nse_grp' 'dat_nme' ; 1:numel(ped_dat.data_name) 1:numel(ped_dat.data_name) 1:numel(ped_dat.data_name)};
    cfg.out_dir     = [fcfg.out_pth '/' 'cmn_nse_rmv'];
    cfg.rmv_chn     = nse_val;
    ped_dat = ft_func(@ft_remove_common_noise2,cfg,ped_dat);
        
    if nse_val == 1
        cfg.clr_fld = fcfg.clr_fld;
        cfg.sbj_nme = fcfg.sbj_nme;
        cfg.sub_fld_ind = num2cell(1:numel(ped_dat.data_name));
        cfg.specific    = {'pre_fix'  'sub_fld_ind' ; 1:numel(ped_dat.data_name) 1:numel(ped_dat.data_name)};
        ped_dat = ft_func(@mmil_update_cmn_nse2,cfg,ped_dat);
    end
        
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

%% Save for Frequency/Machine Learning
% Save usable channels
% save([fcfg.sbj_dat_hld '/' 'raw_data' '/' 'cleaned' '/' fcfg.sbj_nme '_good_channels.mat'],'gdd_grd_lab','bdd_grd_lab','gdd_tal_lab','-v7.3');

% Save
% cfg = [];
% cfg.str_nme  = 'mch_dat';
% cfg.save     = 'yes';
% cfg.filename = [fcfg.prj_dat_hld '/' 'clerical' '/' 'raw_data' '/' 'cleaned' '/' fcfg.sbj_nme '_cleaned_data.mat'];
% ft_func([],cfg,ped_dat);

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
cfg.baselinewindow = [-0.200 0.000];
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

hld_dat = ped_dat;
ped_dat = hld_dat;

%% Preliminary Plots
ped_dat.(ped_dat.data_name{2}).cfg.alt_eve.trialinfo = ped_dat.(ped_dat.data_name{2}).trialinfo;
ped_dat.(ped_dat.data_name{2}).cfg.alt_lab.label = ped_dat.(ped_dat.data_name{2}).label;

cfg = [];

cfg.dat       = {ped_dat.(ped_dat.data_name{2})};
cfg.alt_eve   = 'trialinfo';
cfg.plt_dim   = [1 1];
cfg.type      = 'chan';

eve_plt     = 1:3; %unique(ped_dat.(ped_dat.data_name{2}).trialinfo);

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

cfg.x_lim       = [-0.300 1.000];
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
cfg.outdir     = [fcfg.prj_dat_hld '/' 'epoch_data' '/' 'plot_explore' '/' fcfg.sbj_nme '/' 'channels' '/' ped_dat.data_name{2}];
cfg.prefix     = [ped_dat.data_name{2}] ;

mmil_ieeg_sensor_plot_v5(cfg)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(unique(ped_dat.(ped_dat.data_name{2}).trialinfo)==11)
    
    ped_dat.(ped_dat.data_name{2}).cfg.alt_eve.trialinfo = ped_dat.(ped_dat.data_name{2}).trialinfo;
    ped_dat.(ped_dat.data_name{2}).cfg.alt_lab.label = ped_dat.(ped_dat.data_name{2}).label;
    
    cfg = [];
    
    cfg.dat       = {ped_dat.(ped_dat.data_name{2})};
    cfg.alt_eve   = 'trialinfo';
    cfg.plt_dim   = [1 1];
    cfg.type      = 'chan';
    
    eve_plt     = 11:14; %unique(ped_dat.(ped_dat.data_name{2}).trialinfo);
    
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
    
    cfg.x_lim       = [-0.300 1.000];
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
    cfg.outdir     = [fcfg.prj_dat_hld '/' 'epoch_data' '/' 'plot_explore' '/' fcfg.sbj_nme '/' 'channels' '/' ped_dat.data_name{2} '_broca'];
    cfg.prefix     = [ped_dat.data_name{2}] ;
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
end

%% Save
cfg = [];
cfg.str_nme  = 'ped_dat';
cfg.save     = 'yes';
cfg.filename = [fcfg.out_pth '/' fcfg.sbj_nme '_overall_data.mat'];
ft_func([],cfg,ped_dat);

end