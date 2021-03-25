clear; clc;

prj_dat_hld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/';
clr_fld     = [prj_dat_hld '/clerical/'];

sbj_nme_hld = mmil_readtext([clr_fld 'subjects']);
sbj_clr_hld = mmil_readtext([clr_fld 'subjects_clerical']);
sbj_dat_hld = mmil_readtext([clr_fld 'subjects_data']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR THE PAPER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clinical Data
cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
cfg.inc_clm = {'age' 'onset' 'sex' 'handedness' 'wada'};
mmil_ovr_cln(cfg);

%% Overall Electrodes
sbj_nme = mmil_readtext([clr_fld '/' 'subjects' ]);
bim={}; vis={}; aud={};
for iS = 1:numel(sbj_nme)
    tsk = mmil_load_subj_info([clr_fld '/' 'sbj_inf' '/' sbj_nme{iS}],'tsk');
    if numel(tsk)>1;              bim = [bim(:)' sbj_nme(iS)]; end
    if any(strcmpi(tsk,'SZ'));    vis = [vis(:)' sbj_nme(iS)]; end
    if any(strcmpi(tsk,'SA'));    aud = [aud(:)' sbj_nme(iS)]; end
end

cfg = [];
cfg.loc = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
cfg.sbj = { bim vis aud };

cfg.nme = { 'bimodal' 'visual' 'auditory' };
mmil_ovr_ele_loc(cfg);

%% Figure2
iSASZ_figure2

%% Figure3
iSASZ_figure3

%% Figure4
iSASZ_figure4

%% Figure5
iSASZ_Figure5_v2

%% Figure6
iSASZ_Figure6_v2

%% Figure 7
iSASZ_Figure7

%% Figure 8


%% Figure 9
iSASZ_Figure9

%% Tables
iSASZ_eff_tbl

%% Stats
iSASZ_new_bse_stt
iSASZ_rep_stt_v2
iSASZ_n400_stt
iSASZ_Bi_VS_Uni
iSASZ_amplitude
% iSASZ_stm_stt
iSASZ_uth_stt
iSASZ_plv_stt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Connectivity
iSASZ_connectivity;

%% Redo originals
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'ani_obj_800_vis' 'ani_obj_800_aud' 'ani_obj_800_bim' ...
                'nov_rep_800_vis' 'nov_rep_800_aud' 'nov_rep_800_bim' ...
                'vis_lng_crr_600' 'vis_big_crr_600' 'vis_frq_crr_600' 'vis_ngh_crr_600' ...
                'aud_lng_crr_600' 'aud_big_crr_600' 'aud_frq_crr_600' 'aud_ngh_crr_600'};
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'ani_obj_800_vis' 'ani_obj_800_aud' 'ani_obj_800_bim' ...
                'nov_rep_800_vis' 'nov_rep_800_aud' 'nov_rep_800_bim' ...
                'vis_lng_crr_600' 'vis_big_crr_600' 'vis_frq_crr_600' 'vis_ngh_crr_600' ...
                'aud_lng_crr_600' 'aud_big_crr_600' 'aud_frq_crr_600' 'aud_ngh_crr_600' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);


%% Initial Channels
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
    iSASZ_pap_eff_hgp_act(cfg)

%     iSASZ_pap_eff_lfp_act(cfg)
    
    % Analysis
    iSASZ_pap_ana_hgp_act(cfg)    

%     iSASZ_pap_ana_lfp_act(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_act' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_act' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_act' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

%% Repetition Effects - HGP
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask
    iSASZ_msk_rep(cfg)    
    
    % Effects
    iSASZ_pap_eff_hgp_rep(cfg)

    % Analysis
    iSASZ_pap_ana_hgp_rep(cfg)    

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_rep' }; % 
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_rep' }; % 
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_rep'}; % 
cfg.typ     = { 'hgp' };
mmil_overall_analysis(cfg);

%% Repetition Effects - LFP
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
    iSASZ_pap_eff_lfp_rep(cfg)
    
    % Analysis
    iSASZ_pap_ana_lfp_rep(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_rep_nob' }; % 
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_rep_nob' }; % 
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_rep_nob'}; % 
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% N400
for sbj_num = 1:numel(sbj_nme_hld);
    
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
    iSASZ_pap_eff_lfp_n400(cfg)
    
    % Analysis
    iSASZ_pap_ana_lfp_n400(cfg)
    
end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% N400
for sbj_num = 1:numel(sbj_nme_hld);
    
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Effects
    iSASZ_pap_eff_lfp_n400(cfg)
    
    % Analysis
    iSASZ_pap_ana_lfp_n400(cfg)
    
end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_aud_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_bim_n400' };
cfg.typ     = { 'lfp' };
mmil_overall_analysis(cfg);

%% Utah Array
cfg = [];
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/Utah';
cfg.sbj_nme = {'MG49_SZ' 'MG49_SA'};
iSASZ_utah_array(cfg)

% Stats
cfg = [];
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
cfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
cfg.sbj_nme = 'MG49_SA_SZ_utah';
iSASZ_utah_array_stats(cfg)

%% Animal/Object
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask
    iSASZ_msk_sem(cfg)    
    
    % Effects
    iSASZ_pap_eff_hgp_sem(cfg)

    iSASZ_pap_eff_lfp_sem(cfg)
    
    % Analysis
    iSASZ_pap_ana_hgp_sem(cfg)    

    iSASZ_pap_ana_lfp_sem(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_sem' 'pap_aud_sem' 'pap_bim_sem' }; % 'pap_vis_sem' 'pap_aud_sem'
cfg.typ     = { 'hgp' 'lfp' };
mmil_overall_analysis(cfg);

%% Targetted Animal/Object
for sbj_num = 1:numel(sbj_nme_hld);
    
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Events, Stats, & Mask
    iSASZ_trg_sem_eve(cfg);
    
    % Effects
    iSASZ_pap_eff_hgp_trg_sem(cfg)

    iSASZ_pap_eff_lfp_trg_sem(cfg)
    
    % Analysis
    iSASZ_pap_ana_hgp_trg_sem(cfg)    

    iSASZ_pap_ana_lfp_trg_sem(cfg) 
        
end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_trg_sem' 'pap_aud_trg_sem' 'pap_bim_trg_sem' };
cfg.typ     = { 'hgp' 'lfp' };
mmil_overall_analysis(cfg);

%% Lexical
for sbj_num = 1:numel(sbj_nme_hld);
           
    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = [sbj_dat_hld{sbj_num} '/' 'epoch_data'];
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    % Mask
    iSASZ_msk_lex(cfg)    
    
    % Effects
    iSASZ_pap_eff_hgp_lex(cfg)

    iSASZ_pap_eff_lfp_lex(cfg)
    
    % Analysis
    iSASZ_pap_ana_hgp_lex(cfg)    

    iSASZ_pap_ana_lfp_lex(cfg) 

end

% total
cfg = [];
cfg.dat_fld     = sbj_dat_hld{1};
cfg.clr_fld     = sbj_clr_hld{1};
cfg.inc_sbj     = sbj_nme_hld;
cfg.ana_nme = { 'pap_vis_lex' 'pap_aud_lex' 'pap_bim_lex' };
cfg.typ     = { 'hgp' 'lfp' };
mmil_overall_analysis(cfg);

%% Depth Analysis
% Previously Analyzed
for sbj_num = 19:numel(sbj_nme_hld);

    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ';
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/depth';
    
    % Load
    iSASZ_load_depth(cfg)
    
    % Events
    cfg.sbj_dat_hld = cfg.out_pth;
    iSASZ_events_depth(cfg)
    
    % Stats
    cfg.sbj_dat_hld = cfg.out_pth;
    iSASZ_stats_depth(cfg)
    
    % Plot
    iSASZ_plot_depth(cfg)

end

% New Depth Only Analyzed
sbj_nme_hld = mmil_readtext([clr_fld '/' 'subjects_depth']);

for sbj_num = 3:numel(sbj_nme_hld);

    cfg = [];
    
    % Setting up Variables
    cfg.sbj_nme     = sbj_nme_hld{sbj_num};
    cfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ';
    cfg.clr_fld     = sbj_clr_hld{sbj_num};
    
    cfg.out_pth = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/depth';
    
    % Load
    iSASZ_load_depth_1st(cfg)
    
    % Events
    cfg.sbj_dat_hld = cfg.out_pth;
    iSASZ_events_depth(cfg)
    
    % Stats
    cfg.sbj_dat_hld = cfg.out_pth;
    iSASZ_stats_depth(cfg)
    
    % Plot
    iSASZ_plot_depth(cfg)

end

 