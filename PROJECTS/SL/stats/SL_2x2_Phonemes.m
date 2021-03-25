clear; clc;

sbj_num = 2;

%% Load
sbj_nme  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
sbj_nme  = sbj_nme{sbj_num};

cfg = [];
cfg.load    = 'yes';
cfg.file    = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/' '/' sbj_nme '_overall_data.mat'];
smt_rea_dat = ft_func([],cfg);

if strcmpi(sbj_nme,'NY439_SL'); smt_rea_dat.NY439_SL_Day3_Block1_1_Clin1_hgp.label = smt_rea_dat.NY439_SL_Day3_Block1_1_Clin1_hgp.cfg.alt_lab.repair_macro_ejk1st_meso_label; end

cfg = [];
cfg.data_name = [1 3];
cfg.rmfield   = 'yes';
smt_rea_dat   = ft_func([],cfg,smt_rea_dat);

%% for each stimuli, find the correct timing
if strcmpi(sbj_nme,'NY439_SL');
    tme_dat = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim/female_phoneme_timing.csv');
    
    con = {'N' 'B' 'S' 'W' 'M' 'T' 'H' 'F' 'L' 'Y' 'R' 'K'};
    vow = {'IH' 'UH' 'OH' 'AY' 'EE' 'OO' 'AH' 'UR' 'EH' 'AW' 'OY'};
else
    tme_dat = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/stim2/female_phoneme_timing.csv');
    
    con = {'N'  'T'  'S'  'W'  'M'  'D'  'H'  'F'  'L'  'Y'  'R'  'K'};
    vow = {'EE' 'UH' 'OH' 'AY' 'EE' 'OO' 'AH' 'UR' 'EH' 'AW' 'OY' 'IH'};
end

smt_off_set = zeros(numel(smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.alt_eve.a_con),1);
for iTM = 1:numel(smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.alt_eve.a_con)
    smt_off_set(iTM) = round(tme_dat{find(strcmpi(tme_dat(:,1),[con{smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.alt_eve.a_con(iTM)} vow{smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.alt_eve.a_vow(iTM)-100}])),3} * smt_rea_dat.(smt_rea_dat.data_name{1}).fsample);
end

%% change beggining and end %%%%%%%%%%%%%%%%%%%%%
smt_cur_tme            = smt_rea_dat.(smt_rea_dat.data_name{1}).sampleinfo;

% consonant trials
smt_con_new_trl = zeros(size(smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.trl,1),4);
for iTR = 1:size(smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.trl,1)
    smt_con_new_trl(iTR,1) = smt_cur_tme(iTR,1);
    smt_con_new_trl(iTR,2) = smt_cur_tme(iTR,2);
    smt_con_new_trl(iTR,3) = -486; %round(-(rea_dat.fsample*0.8)); round(0.15*rea_dat.fsample)+off_set(iTR)
    smt_con_new_trl(iTR,4) = smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.trl(iTR,4);
end
% vowel trials
smt_vow_new_trl = zeros(size(smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.trl,1),4);
for iTR = 1:size(smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.trl,1)
    smt_vow_new_trl(iTR,1) = smt_cur_tme(iTR,1);
    smt_vow_new_trl(iTR,2) = smt_cur_tme(iTR,2);
    smt_vow_new_trl(iTR,3) = -486 - smt_off_set(iTR); %round(-(rea_dat.fsample*0.8)); round(0.15*rea_dat.fsample)+off_set(iTR)
    smt_vow_new_trl(iTR,4) = smt_rea_dat.(smt_rea_dat.data_name{1}).cfg.trl(iTR,4);
end

%% re-epoch %%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.trl = smt_con_new_trl;
smt_con_alg_dat = ft_func(@ft_redefinetrial,cfg,smt_rea_dat);

cfg = [];
cfg.latency = [-0.2 0.6];
smt_con_alg_dat = ft_func(@ft_selectdata,cfg,smt_con_alg_dat);

cfg = [];
cfg.trl = smt_vow_new_trl;
smt_vow_alg_dat = ft_func(@ft_redefinetrial,cfg,smt_rea_dat);

cfg = [];
cfg.latency = [-0.2 0.6];
smt_vow_alg_dat = ft_func(@ft_selectdata,cfg,smt_vow_alg_dat);

cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = [-0.20 0];
smt_con_alg_dat            = ft_func(@ft_preprocessing,cfg,smt_con_alg_dat);
smt_vow_alg_dat            = ft_func(@ft_preprocessing,cfg,smt_vow_alg_dat);

smt_con_alg_dat = smt_con_alg_dat.NY439_SL_Day3_Block1_1_Clin1_hgp; %%%%%%%%%%%%%%%%%%%
smt_vow_alg_dat = smt_vow_alg_dat.NY439_SL_Day3_Block1_1_Clin1_hgp; %%%%%%%%%%%%%%%%%%%

smt_con_ovr_dat_hld = cat(3,smt_con_alg_dat.trial{:});
smt_vow_ovr_dat_hld = cat(3,smt_vow_alg_dat.trial{:});

%% Make pools of stimuli %%%%%%%%%%%%%%%%%%%%%
if strcmpi(sbj_nme,'NY439_SL');
    ttt.con1 = {'N' 'B' 'S' 'W'}; ttt.vow1 = {'IH' 'UH' 'OH' 'AY'}; cnt = 1; for iC = 1:numel(ttt.con1); for iV = 1:numel(ttt.vow1); ttt.stm1{cnt} = [ttt.con1{iC} ttt.vow1{iV}]; ttt.stm1_vow_idn(cnt) = iV; ttt.stm1_con_idn(cnt) = iC; cnt = cnt+1; end; end
    ttt.con2 = {'M' 'T' 'H' 'F'}; ttt.vow2 = {'EE' 'OO' 'AH' 'UR'}; cnt = 1; for iC = 1:numel(ttt.con2); for iV = 1:numel(ttt.vow2); ttt.stm2{cnt} = [ttt.con2{iC} ttt.vow2{iV}]; ttt.stm2_vow_idn(cnt) = iV; ttt.stm2_con_idn(cnt) = iC;cnt = cnt+1; end; end
    ttt.con3 = {'L' 'Y' 'R' 'K'}; ttt.vow3 = {'EH' 'AW' 'OY'}; cnt = 1; for iC = 1:numel(ttt.con3); for iV = 1:numel(ttt.vow3); ttt.stm3{cnt} = [ttt.con3{iC} ttt.vow3{iV}]; ttt.stm3_vow_idn(cnt) = iV; ttt.stm3_con_idn(cnt) = iC; cnt = cnt+1; end; end
else
    ttt.con1 = {'N'  'T'  'S'  'W'}; ttt.vow1 = {'EE' 'UH' 'OH' 'AY'}; cnt = 1; for iC = 1:numel(ttt.con1); for iV = 1:numel(ttt.vow1); ttt.stm1{cnt} = [ttt.con1{iC} ttt.vow1{iV}]; ttt.stm1_vow_idn(cnt) = iV; ttt.stm1_con_idn(cnt) = iC; cnt = cnt+1; end; end
    ttt.con2 = {'M'  'D'  'H'  'F'}; ttt.vow2 = {'EE' 'OO' 'AH' 'UR'}; cnt = 1; for iC = 1:numel(ttt.con2); for iV = 1:numel(ttt.vow2); ttt.stm2{cnt} = [ttt.con2{iC} ttt.vow2{iV}]; ttt.stm2_vow_idn(cnt) = iV; ttt.stm2_con_idn(cnt) = iC;cnt = cnt+1; end; end
    ttt.con3 = {'L'  'Y'  'R'  'K'}; ttt.vow3 = {'EH' 'AW' 'OY' 'IH'}; cnt = 1; for iC = 1:numel(ttt.con3); for iV = 1:numel(ttt.vow3); ttt.stm3{cnt} = [ttt.con3{iC} ttt.vow3{iV}]; ttt.stm3_vow_idn(cnt) = iV; ttt.stm3_con_idn(cnt) = iC; cnt = cnt+1; end; end 
end

int_num = [1  2  3  4  0  0  0  0  0  0  0  0;  ...
           5  6  7  8  0  0  0  0  0  0  0  0;  ...
           9  10 11 12 0  0  0  0  0  0  0  0;  ...
           13 14 15 16 0  0  0  0  0  0  0  0;  ...
           0  0  0  0  17 18 19 20 0  0  0  0;  ...
           0  0  0  0  21 22 23 24 0  0  0  0 ;  ...
           0  0  0  0  25 26 27 28 0  0  0  0;  ...
           0  0  0  0  29 30 31 32 0  0  0  0;  ...
           0  0  0  0  0  0  0  0  33 34 35 36; ...
           0  0  0  0  0  0  0  0  37 38 39 40; ...
           0  0  0  0  0  0  0  0  41 42 43 44; ...
           0  0  0  0  0  0  0  0  45 46 47 48]; 
       
for i = 1:numel(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_vow)
    smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_int(i) = int_num(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_con(i),smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_vow(i)-100)+200; smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_int(i) = int_num(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_con(i),smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_vow(i)-100)+200;
    smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_int(i) = int_num(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_con(i),smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_vow(i)-100)+200; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_int(i) = int_num(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_con(i),smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_vow(i)-100)+200;
end

smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_con_nse = smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_con; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_con_nse = smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_con;
smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_vow_nse = smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_vow; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_vow_nse = smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_vow;
smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_int_nse = smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_int; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_int_nse = smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_int;

smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_con_nse(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.trialinfo ~= 4) = 0; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_con_nse(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.trialinfo ~= 4) = 0;
smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_vow_nse(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.trialinfo ~= 4) = 0; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_vow_nse(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.trialinfo ~= 4) = 0;
smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_int_nse(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.trialinfo ~= 4) = 0; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_int_nse(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.trialinfo ~= 4) = 0;

smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_con(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.trialinfo == 4) = 0; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_con(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.trialinfo == 4) = 0;
smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_vow(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.trialinfo == 4) = 0; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_vow(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.trialinfo == 4) = 0;
smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.a_int(smt_con_alg_dat.(smt_con_alg_dat.data_name{1}).cfg.alt_eve.trialinfo == 4) = 0; smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.a_int(smt_vow_alg_dat.(smt_vow_alg_dat.data_name{1}).cfg.alt_eve.trialinfo == 4) = 0;

eve_for_plt = {[1 3 4] [5 6 7 8] [9 10 11 12]; ...
    [101 102 103 104] [105 106 107 108] [109 110 111 112]; ...
    201:216 217:232 233:248}';

plt_col     = distinguishable_colors(numel(unique(cat(2,eve_for_plt{[1 4 7 2 5 8]})))+4);
con_plt_col = plt_col([1:12],:);
vow_plt_col = plt_col(13:24,:);
[ttt2,ttt1] = find(int_num');
for j = 1:numel(ttt1)
    cmb_col(j,:) = sum([vow_plt_col(ttt2(j),:) ; con_plt_col(ttt1(j),:)]) / 2;
end







