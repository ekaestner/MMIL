clear; clc;

cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';

sbj_nme_hld = mmil_readtext([cfg.clr_fld 'subjects']);

%% Run Subjects
for sbj_num = 1:numel(sbj_nme_hld);
    
    %% Setting up Variables
    cfg.sbj_nme = sbj_nme_hld{sbj_num};
    cfg.fle_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst/';
       
    %% load the data
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst/';
    
    iSASZ_initial_load(cfg)
           
    %% run the initial stats
    cfg.out_pth =  '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst/';
    
    iSASZ_initial_statsv2(cfg)
    
    %% check initial stats
    cfg.chn_chs = 1;
    cfg.mve_plt = 1;
    cfg.ovr_wrt = 1;
    
    cfg.clr_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/';
    cfg.fle_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/';
    
    iSASZ_check_stats_v2(cfg)
    
    %% run the initial plots
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst/';
    cfg.ovr_wrt = 2;
    
    iSASZ_initial_plot_v2(cfg)
   
    %% Run Location Scripts    
    Fix Channel Names
   
    % Make ECOG/DEPTH Split
    cfg.alt_lab = 'label';
    mmil_create_depth(cfg)
    
    Split ECOG/Depth    
    cfg.plt_spl = [cfg.clr_out_pth '/' cfg.sbj_nme '/' 'ovr'];
    cfg.idn_spl = 1;
    mmil_split_plot(cfg);
    mmil_split_stat(cfg)
    
end

%% Fix problems
for sbj_num = [35 36];
   
% Load Subjects
infile = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/' sbj_nme_hld{sbj_num} '_overall_data.mat'];
cfg = [];
cfg.load = 'yes';
cfg.file = infile;
sem_dat  = ft_func([],cfg);

% Average Noise Subtraction
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';
chn_nse_grp = mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' sbj_nme_hld{sbj_num}],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
if ~isempty([chn_nse_grp{:}]); chn_nse_grp = cellfun(@str2num,mmil_load_subj_info([cfg.clr_fld '/' 'sbj_inf' '/' sbj_nme_hld{sbj_num}],'cmn_nse'),'uni',0)'; end

for i = 1:numel(chn_nse_grp)
    if ~isempty(chn_nse_grp{i})
        cfg             = [];
        cfg.chn_nse_grp = chn_nse_grp(i);
        sem_dat = ft_func(@ft_remove_common_noise,cfg,sem_dat);
    end
end

cfg = [];
cfg.str_nme  = 'sem_dat';
cfg.save     = 'yes';
cfg.sve_app  = 'app_all';
cfg.filename =['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst2' '/' sbj_nme_hld{sbj_num} '_overall_data.mat'];
ft_func([],cfg,sem_dat);

% Make New Stats
cfg.sbj_nme = sbj_nme_hld{sbj_num};
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';
cfg.out_pth =  '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst2/';
iSASZ_initial_statsv2(cfg)

% Make New Plots
cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst2/';
cfg.in_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master_tst2/';
cfg.ovr_wrt = 2;
iSASZ_initial_plot_v2(cfg)

end

%% Examine results
% Make stat table
typ = {'lfp' 'hgp'};

sig_in__pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/sig_chn';
plt_in__pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat';
loc_in__pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical';

tab_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/sig_chn';
tab_fle_nme = { 'overall'    ...
                'repetition' };          

tab_col = { {'vis_ovr_all_pre' 'vis_ovr_all_lex' 'vis_ovr_all_pos' 'aud_ovr_all_pre' 'aud_ovr_all_lex' 'aud_ovr_all_pos'} ...
            {'vis_dif_all_pre' 'vis_dif_all_lex' 'vis_dif_all_pos' 'aud_dif_all_pre' 'aud_dif_all_lex' 'aud_dif_all_pos'} };
 
ovr_lap_col = {[1 4] [2 5] [3 6]};
ovr_lap_lbl = { {'ovr_all_pre_lap' 'ovr_all_lex_lap' 'ovr_all_pos_lap'} ...
                {'dif_all_pre_lap' 'dif_all_lex_lap' 'dif_all_pos_lap'} };

plt_out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/sig_chn_plt';
plt_fle_nme = { 'overall'    ...
                'repetition' };

plt_grp     = {[1 4] [2 5] [3 6]};
plt_grp_nme = {'pre' 'lex' 'pos'};

for iT = 1:numel(tab_col) % Baseline & Repetition
    
    cfg = [];
    cfg.typ = typ; % LFP & HGP
    
    cfg.tab_col = tab_col{iT};
    
    cfg.ovr_lap_col = ovr_lap_col;
    cfg.ovr_lap_lbl = ovr_lap_lbl{iT};
    
    cfg.sig_in__pth = sig_in__pth;
    cfg.plt_in__pth = plt_in__pth;
    cfg.loc_in__pth = loc_in__pth;
    
    cfg.tab_out_pth = tab_out_pth;
    cfg.tab_fle_nme = tab_fle_nme{iT};
        
    cfg.plt_out_pth = plt_out_pth;
    cfg.plt_fle_nme = plt_fle_nme{iT};
    
    cfg.plt_grp     = plt_grp;
    cfg.plt_grp_nme = plt_grp_nme;
    
    cfg.cmb_plt_typ = 1;
    cfg.loc_tab     = 1;
    cfg.str_sbj     = 15;
    
    cfg.ecg = 1;
    cfg.dep = 1;
    
    iSASZ_stat_table(cfg);
    
end

% Move plots for visual inspection
cfg = [];
cfg.plt_sig_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/sig_chn';
cfg.plt_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/';
cfg.plt_out_fld = '';
cfg.sub_plt_fld = 'ovr';
mmil_move_sig_plot(cfg);

% Make individual location plots 
cfg = [];
cfg.idv_plt = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';
cfg.loc_hld = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/sig_chn' '/' 'overall_lfp_ecog_location_plot.mat'];
cfg.sve_loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/location_plot';
cfg.sve_pre = 'overall_lfp_ecog';
cfg.col_tab = {'vis_ovr_all_pre' rgb('light red') ; ...
               'vis_ovr_all_lex' rgb('red') ; ...
               'vis_ovr_all_pos' rgb('dark red') ; ...
               'aud_ovr_all_pre' rgb('light blue') ; ...
               'aud_ovr_all_lex' rgb('blue') ; ...
               'aud_ovr_all_pos' rgb('dark blue') ; ...
               'ovr_all_pre_lap' rgb('light purple') ; ...
               'ovr_all_lex_lap' rgb('purple') ; ...
               'ovr_all_pos_lap' rgb('dark purple')};
mmil_location_plot(cfg)

% Make group location plots
cfg = [];
cfg.grp_plt = 1;
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/clerical/';
cfg.loc_hld = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/sig_chn' '/' 'overall_lfp_ecog_location_plot.mat'];
cfg.sve_loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data_master/paper_plot/ThreeStat/location_plot';
cfg.sve_pre = 'overall_lfp_ecog';
cfg.col_tab = {'vis_ovr_all_pre' rgb('light red') ; ...
               'vis_ovr_all_lex' rgb('red') ; ...
               'vis_ovr_all_pos' rgb('dark red') ; ...
               'aud_ovr_all_pre' rgb('light blue') ; ...
               'aud_ovr_all_lex' rgb('blue') ; ...
               'aud_ovr_all_pos' rgb('dark blue') ; ...
               'ovr_all_pre_lap' rgb('light purple') ; ...
               'ovr_all_lex_lap' rgb('purple') ; ...
               'ovr_all_pos_lap' rgb('dark purple')};
mmil_location_plot(cfg)









