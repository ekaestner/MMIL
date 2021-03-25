clear; clc;

sbj_num = 1;
ovr_wrt = 1;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/';

cln_dir = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical';

try chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical' '/' subj '_loc']);
catch
    has_loc = 0;
end
if ~isvar('has_loc'); has_loc = 1; end

cfg = [];
cfg.load    = 'yes';
cfg.file = [outpath '/' subj '_smoothed_hgp.mat'];
sll_dat     = ft_func([],cfg);

fprintf('Starting work on %s \n',subj)

%% Auditory 2-x-2 Anova
hgp_ind = 1;
anv_dat = sll_dat.(sll_dat.data_name{hgp_ind});

vow_ind = {1:4 5:8 9:11};
con_ind = {1:4 5:8 9:12};

tme_pnt = numel(anv_dat.time{1});

for iBL = 1:numel(vow_ind)
    for iVW = 1:numel(vow_ind{iBL})
        for iCN = 1:numel(con_ind{iBL})
            trl_ind{iBL}{iVW,iCN} = find(anv_dat.cfg.alt_eve.a_vow==vow_ind{iBL}(iVW) & anv_dat.cfg.alt_eve.a_con==con_ind{iBL}(iCN) & anv_dat.cfg.alt_eve.trialinfo~=4);
            vow_idn{iBL}{iVW,iCN} = repmat({num2str(iVW)},1,numel(trl_ind{iBL}{iVW,iCN}));
            con_idn{iBL}{iVW,iCN} = repmat({num2str(iCN)},1,numel(trl_ind{iBL}{iVW,iCN}));
        end
    end
end

ovr_dat_hld = cat(3,anv_dat.trial{:});

aud_vow_anv_stt = zeros(numel(vow_ind),numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label),tme_pnt);
aud_con_anv_stt = zeros(numel(vow_ind),numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label),tme_pnt);
aud_int_anv_stt = zeros(numel(vow_ind),numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label),tme_pnt);

for iCH = 1:numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label)
    for iBL = 1:numel(vow_ind)
                
        % Create dat
        p_hld = zeros(3,tme_pnt);
        dat_hld = squeeze(ovr_dat_hld(iCH,:,cat(1,trl_ind{iBL}{:})'));
        g_vow   = cat(2,vow_idn{iBL}{:});
        c_vow   = cat(2,con_idn{iBL}{:});
        
        % Run Tests
        for iTM = 1:tme_pnt
            p = anovan(dat_hld(iTM,:),{g_vow,c_vow},'model','interaction','display','off');
            aud_vow_anv_stt(iBL,iCH,iTM) = p(1);
            aud_con_anv_stt(iBL,iCH,iTM) = p(2);
            aud_int_anv_stt(iBL,iCH,iTM) = p(3);
        end
        
    end
end

anv_dat.cfg.alt_stt.aud_vow_blk1.prob = squeeze(aud_vow_anv_stt(1,:,:)); anv_dat.cfg.alt_stt.aud_vow_blk1.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_vow_blk1.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_vow_blk2.prob = squeeze(aud_vow_anv_stt(2,:,:)); anv_dat.cfg.alt_stt.aud_vow_blk2.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_vow_blk2.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_vow_blk3.prob = squeeze(aud_vow_anv_stt(3,:,:)); anv_dat.cfg.alt_stt.aud_vow_blk3.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_vow_blk3.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_con_blk1.prob = squeeze(aud_con_anv_stt(1,:,:)); anv_dat.cfg.alt_stt.aud_con_blk1.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_con_blk1.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_con_blk2.prob = squeeze(aud_con_anv_stt(2,:,:)); anv_dat.cfg.alt_stt.aud_con_blk2.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_con_blk2.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_con_blk3.prob = squeeze(aud_con_anv_stt(3,:,:)); anv_dat.cfg.alt_stt.aud_con_blk3.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_con_blk3.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_int_blk1.prob = squeeze(aud_int_anv_stt(1,:,:)); anv_dat.cfg.alt_stt.aud_int_blk1.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_int_blk1.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_int_blk2.prob = squeeze(aud_int_anv_stt(2,:,:)); anv_dat.cfg.alt_stt.aud_int_blk2.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_int_blk2.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.aud_int_blk3.prob = squeeze(aud_int_anv_stt(3,:,:)); anv_dat.cfg.alt_stt.aud_int_blk3.label = anv_dat.label; anv_dat.cfg.alt_stt.aud_int_blk3.time = anv_dat.time{1};

% Add FDR correction
ttt = fieldnames(anv_dat.cfg.alt_stt);
ttt = ttt(~cellfun(@isempty,regexpi(ttt,'aud.+blk')));
for iFD = 1:numel(ttt)
    anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']) = anv_dat.cfg.alt_stt.(ttt{iFD});
    for iCH = 1:size(anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob,1)
        anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob(iCH,:) = fdr_correction(anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob(iCH,:),0.05);
    end
end

%% Plot Auditory ANOVA
% Make interaction events
anv_dat.cfg.alt_eve.a_vow = anv_dat.cfg.alt_eve.a_vow+100;
int_num = [1  2  3  4  0  0  0  0  0  0  0;  ...
           5  6  7  8  0  0  0  0  0  0  0;  ...
           9  10 11 12 0  0  0  0  0  0  0;  ...
           13 14 15 16 0  0  0  0  0  0  0;  ...
           0  0  0  0  17 18 19 20 0  0  0;  ...
           0  0  0  0  21 22 23 24 0  0  0;  ...
           0  0  0  0  25 26 27 28 0  0  0;  ...
           0  0  0  0  29 30 31 32 0  0  0;  ...
           0  0  0  0  0  0  0  0  33 34 35; ...
           0  0  0  0  0  0  0  0  36 37 38; ...
           0  0  0  0  0  0  0  0  39 40 41; ...
           0  0  0  0  0  0  0  0  42 43 44];         
for i = 1:numel(anv_dat.cfg.alt_eve.a_vow)
    anv_dat.cfg.alt_eve.a_int(i) = int_num(anv_dat.cfg.alt_eve.a_con(i),anv_dat.cfg.alt_eve.a_vow(i)-100)+200;
end

anv_dat.cfg.alt_eve.a_con(anv_dat.cfg.alt_eve.trialinfo == 4) = 0;
anv_dat.cfg.alt_eve.a_vow(anv_dat.cfg.alt_eve.trialinfo == 4) = 0;
anv_dat.cfg.alt_eve.a_int(anv_dat.cfg.alt_eve.trialinfo == 4) = 0;

eve_for_plt = {[1 2 3 4] [5 6 7 8] [9 10 11 12]; ...
                 [101 102 103 104] [105 106 107 108] [109 110 111]; ...
                 201:216 217:232 233:244}';

% Colors
plt_col     = distinguishable_colors(numel(unique(cat(2,eve_for_plt{[1 4 7 2 5 8]}))));
con_plt_col = plt_col(1:12,:);
vow_plt_col = plt_col(13:23,:);
[ttt2,ttt1] = find(int_num');
for j = 1:numel(ttt1)
    cmb_col(j,:) = sum([vow_plt_col(ttt2(j),:) ; con_plt_col(ttt1(j),:)]) / 2;
end

% Interaction Plots
cfg           = [];
cfg.type      = 'chan';
cfg.plt_lbl   = {'Blk1 Con' 'Blk2 Con' 'Blk3 Con'; ...
                 'Blk1 Vow' 'Blk2 Vow' 'Blk3 Vow'; ...
                 'Blk1 Int' 'Blk2 Int' 'Blk3 Int'}';
cfg.dat       = {anv_dat};
cfg.dat_loc   = [1 1 1 ; ...
                 1 1 1 ; ...
                 1 1 1]';
cfg.lgd       = 0;
cfg.plt_dim   = [3 3];
cfg.alt_eve   = {'a_con' 'a_con' 'a_con'; ...
                 'a_vow' 'a_vow' 'a_vow'; ...
                 'a_int' 'a_int' 'a_int'}';
cfg.eve       = eve_for_plt;
cfg.xlim      = [-0.4 1]; %
cfg.y_lim      = 'auto';
cfg.lnstyle.col_ord = [num2cell(con_plt_col,2)' num2cell(vow_plt_col,2)' num2cell(cmb_col,2)'];
cfg.lnstyle.mrk_ord = repmat({''},1,numel(cfg.lnstyle.col_ord));
cfg.lnstyle.lin_ord = repmat({'-'},1,numel(cfg.lnstyle.col_ord));
cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
cfg.v_lne      = [0 0.450 0.900];
cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
% cfg.stt_dat    = {{anv_dat.cfg.alt_stt.aud_con_blk1 anv_dat.cfg.alt_stt.aud_con_blk1_fdr} {anv_dat.cfg.alt_stt.aud_con_blk2  anv_dat.cfg.alt_stt.aud_con_blk2_fdr} {anv_dat.cfg.alt_stt.aud_con_blk3  anv_dat.cfg.alt_stt.aud_con_blk3_fdr}; ...
%                   {anv_dat.cfg.alt_stt.aud_vow_blk1 anv_dat.cfg.alt_stt.aud_vow_blk1_fdr} {anv_dat.cfg.alt_stt.aud_vow_blk2  anv_dat.cfg.alt_stt.aud_vow_blk2_fdr} {anv_dat.cfg.alt_stt.aud_vow_blk3  anv_dat.cfg.alt_stt.aud_vow_blk3_fdr}; ...
%                   {anv_dat.cfg.alt_stt.aud_int_blk1 anv_dat.cfg.alt_stt.aud_int_blk1_fdr} {anv_dat.cfg.alt_stt.aud_int_blk2  anv_dat.cfg.alt_stt.aud_int_blk2_fdr} {anv_dat.cfg.alt_stt.aud_int_blk3  anv_dat.cfg.alt_stt.aud_int_blk3_fdr}}';
cfg.std_err    = 1;
cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'jpg';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'Aud2x2_smoothed' '/'];
cfg.prefix     = 'Aud_2x2_Phoneme';
mmil_ieeg_sensor_plot_v4(cfg)

save(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb' '/' subj '_2x2anova_smoothed.mat'],'anv_dat')

%% Visual 2-x-2 Anova
vow_ind = {1:4 5:8 9:11};
con_ind = {1:4 5:8 9:12};

tme_pnt = numel(anv_dat.time{1});

for iBL = 1:numel(vow_ind)
    for iVW = 1:numel(vow_ind{iBL})
        for iCN = 1:numel(con_ind{iBL})
            trl_ind{iBL}{iVW,iCN} = find(anv_dat.cfg.alt_eve.v_vow==vow_ind{iBL}(iVW) & anv_dat.cfg.alt_eve.v_con==con_ind{iBL}(iCN) & anv_dat.cfg.alt_eve.trialinfo~=3);
            vow_idn{iBL}{iVW,iCN} = repmat({num2str(iVW)},1,numel(trl_ind{iBL}{iVW,iCN}));
            con_idn{iBL}{iVW,iCN} = repmat({num2str(iCN)},1,numel(trl_ind{iBL}{iVW,iCN}));
        end
    end
end

ovr_dat_hld = cat(3,anv_dat.trial{:});

vis_vow_anv_stt = zeros(numel(vow_ind),numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label),tme_pnt);
vis_con_anv_stt = zeros(numel(vow_ind),numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label),tme_pnt);
vis_int_anv_stt = zeros(numel(vow_ind),numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label),tme_pnt);

for iCH = 1:numel(anv_dat.cfg.alt_lab.repair_macro_ejk1st_meso_label)
    for iBL = 1:numel(vow_ind)
                
        % Create dat
        p_hld = zeros(3,tme_pnt);
        dat_hld = squeeze(ovr_dat_hld(iCH,:,cat(1,trl_ind{iBL}{:})'));
        g_vow   = cat(2,vow_idn{iBL}{:});
        c_vow   = cat(2,con_idn{iBL}{:});
        
        % Run Tests
        for iTM = 1:tme_pnt
            p = anovan(dat_hld(iTM,:),{g_vow,c_vow},'model','interaction','display','off');
            vis_vow_anv_stt(iBL,iCH,iTM) = p(1);
            vis_con_anv_stt(iBL,iCH,iTM) = p(2);
            vis_int_anv_stt(iBL,iCH,iTM) = p(3);
        end
        
    end
end

anv_dat.cfg.alt_stt.vis_vow_blk1.prob = squeeze(vis_vow_anv_stt(1,:,:)); anv_dat.cfg.alt_stt.vis_vow_blk1.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_vow_blk1.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_vow_blk2.prob = squeeze(vis_vow_anv_stt(2,:,:)); anv_dat.cfg.alt_stt.vis_vow_blk2.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_vow_blk2.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_vow_blk3.prob = squeeze(vis_vow_anv_stt(3,:,:)); anv_dat.cfg.alt_stt.vis_vow_blk3.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_vow_blk3.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_con_blk1.prob = squeeze(vis_con_anv_stt(1,:,:)); anv_dat.cfg.alt_stt.vis_con_blk1.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_con_blk1.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_con_blk2.prob = squeeze(vis_con_anv_stt(2,:,:)); anv_dat.cfg.alt_stt.vis_con_blk2.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_con_blk2.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_con_blk3.prob = squeeze(vis_con_anv_stt(3,:,:)); anv_dat.cfg.alt_stt.vis_con_blk3.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_con_blk3.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_int_blk1.prob = squeeze(vis_int_anv_stt(1,:,:)); anv_dat.cfg.alt_stt.vis_int_blk1.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_int_blk1.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_int_blk2.prob = squeeze(vis_int_anv_stt(2,:,:)); anv_dat.cfg.alt_stt.vis_int_blk2.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_int_blk2.time = anv_dat.time{1};
anv_dat.cfg.alt_stt.vis_int_blk3.prob = squeeze(vis_int_anv_stt(3,:,:)); anv_dat.cfg.alt_stt.vis_int_blk3.label = anv_dat.label; anv_dat.cfg.alt_stt.vis_int_blk3.time = anv_dat.time{1};

% Add FDR correction
ttt = fieldnames(anv_dat.cfg.alt_stt);
ttt = ttt(~cellfun(@isempty,regexpi(ttt,'vis.+blk')));
for iFD = 1:numel(ttt)
    anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']) = anv_dat.cfg.alt_stt.(ttt{iFD});
    for iCH = 1:size(anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob,1)
        anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob(iCH,:) = fdr_correction(anv_dat.cfg.alt_stt.([ttt{iFD} '_fdr']).prob(iCH,:),0.05);
    end
end

%% Plot Visual ANOVA
% Make interaction events
anv_dat.cfg.alt_eve.v_vow = anv_dat.cfg.alt_eve.v_vow+100;
int_num = [1  2  3  4  0  0  0  0  0  0  0;  ...
           5  6  7  8  0  0  0  0  0  0  0;  ...
           9  10 11 12 0  0  0  0  0  0  0;  ...
           13 14 15 16 0  0  0  0  0  0  0;  ...
           0  0  0  0  17 18 19 20 0  0  0;  ...
           0  0  0  0  21 22 23 24 0  0  0;  ...
           0  0  0  0  25 26 27 28 0  0  0;  ...
           0  0  0  0  29 30 31 32 0  0  0;  ...
           0  0  0  0  0  0  0  0  33 34 35; ...
           0  0  0  0  0  0  0  0  36 37 38; ...
           0  0  0  0  0  0  0  0  39 40 41; ...
           0  0  0  0  0  0  0  0  42 43 44];        
for i = 1:numel(anv_dat.cfg.alt_eve.v_vow)
    anv_dat.cfg.alt_eve.v_int(i) = int_num(anv_dat.cfg.alt_eve.v_con(i),anv_dat.cfg.alt_eve.v_vow(i)-100)+200;
end

anv_dat.cfg.alt_eve.v_con(anv_dat.cfg.alt_eve.trialinfo == 3) = 0;
anv_dat.cfg.alt_eve.v_vow(anv_dat.cfg.alt_eve.trialinfo == 3) = 0;
anv_dat.cfg.alt_eve.v_int(anv_dat.cfg.alt_eve.trialinfo == 3) = 0;

eve_for_plt = {[1 2 3 4] [5 6 7 8] [9 10 11 12]; ...
                 [101 102 103 104] [105 106 107 108] [109 110 111]; ...
                 201:216 217:232 233:244}';

% Colors
plt_col     = distinguishable_colors(numel(unique(cat(2,eve_for_plt{[1 4 7 2 5 8]}))));
con_plt_col = plt_col(1:12,:);
vow_plt_col = plt_col(13:23,:);
[ttt2,ttt1] = find(int_num');
for j = 1:numel(ttt1)
    cmb_col(j,:) = sum([vow_plt_col(ttt2(j),:) ; con_plt_col(ttt1(j),:)]) / 2;
end

% Interaction Plots
cfg           = [];
cfg.type      = 'chan';
cfg.plt_lbl   = {'Blk1 Con' 'Blk2 Con' 'Blk3 Con'; ...
                 'Blk1 Vow' 'Blk2 Vow' 'Blk3 Vow'; ...
                 'Blk1 Int' 'Blk2 Int' 'Blk3 Int'}';
cfg.dat       = {anv_dat};
cfg.dat_loc   = [1 1 1 ; ...
                 1 1 1 ; ...
                 1 1 1]';
cfg.lgd       = 0;
cfg.plt_dim   = [3 3];
cfg.alt_eve   = {'v_con' 'v_con' 'v_con'; ...
                 'v_vow' 'v_vow' 'v_vow'; ...
                 'v_int' 'v_int' 'v_int'}';
cfg.eve       = eve_for_plt;
cfg.xlim      = [-0.4 1]; %
cfg.y_lim      = 'auto';
cfg.lnstyle.col_ord = [num2cell(con_plt_col,2)' num2cell(vow_plt_col,2)' num2cell(cmb_col,2)'];
cfg.lnstyle.mrk_ord = repmat({''},1,numel(cfg.lnstyle.col_ord));
cfg.lnstyle.lin_ord = repmat({'-'},1,numel(cfg.lnstyle.col_ord));
cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
cfg.v_lne      = [0 0.450 0.900];
cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
cfg.stt_dat    = {{anv_dat.cfg.alt_stt.vis_con_blk1 anv_dat.cfg.alt_stt.vis_con_blk1_fdr} {anv_dat.cfg.alt_stt.vis_con_blk2  anv_dat.cfg.alt_stt.vis_con_blk2_fdr} {anv_dat.cfg.alt_stt.vis_con_blk3  anv_dat.cfg.alt_stt.vis_con_blk3_fdr}; ...
                  {anv_dat.cfg.alt_stt.vis_vow_blk1 anv_dat.cfg.alt_stt.vis_vow_blk1_fdr} {anv_dat.cfg.alt_stt.vis_vow_blk2  anv_dat.cfg.alt_stt.vis_vow_blk2_fdr} {anv_dat.cfg.alt_stt.vis_vow_blk3  anv_dat.cfg.alt_stt.vis_vow_blk3_fdr}; ...
                  {anv_dat.cfg.alt_stt.vis_int_blk1 anv_dat.cfg.alt_stt.vis_int_blk1_fdr} {anv_dat.cfg.alt_stt.vis_int_blk2  anv_dat.cfg.alt_stt.vis_int_blk2_fdr} {anv_dat.cfg.alt_stt.vis_int_blk3  anv_dat.cfg.alt_stt.vis_int_blk3_fdr}}';
cfg.std_err    = 1;
cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'jpg';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'Vis2x2_smoothed' '/'];
cfg.prefix     = 'Vis_2x2_Phoneme';
mmil_ieeg_sensor_plot_v4(cfg)


anv_dat.cfg.alt_eve.v_con = sll_dat.(sll_dat.data_name{hgp_ind}).cfg.alt_eve.v_con;
anv_dat.cfg.alt_eve.v_vow = sll_dat.(sll_dat.data_name{hgp_ind}).cfg.alt_eve.v_vow;
anv_dat.cfg.alt_eve.a_con = sll_dat.(sll_dat.data_name{hgp_ind}).cfg.alt_eve.a_con;
anv_dat.cfg.alt_eve.a_vow = sll_dat.(sll_dat.data_name{hgp_ind}).cfg.alt_eve.a_vow;

save(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb' '/' subj '_2x2anova_smoothed.mat'],'anv_dat')

%% Noise plots
% Plot Auditory ANOVA
% Make interaction events
anv_dat.cfg.alt_eve.a_vow = anv_dat.cfg.alt_eve.a_vow+100;
int_num = [1  2  3  4  0  0  0  0  0  0  0;  ...
           5  6  7  8  0  0  0  0  0  0  0;  ...
           9  10 11 12 0  0  0  0  0  0  0;  ...
           13 14 15 16 0  0  0  0  0  0  0;  ...
           0  0  0  0  17 18 19 20 0  0  0;  ...
           0  0  0  0  21 22 23 24 0  0  0;  ...
           0  0  0  0  25 26 27 28 0  0  0;  ...
           0  0  0  0  29 30 31 32 0  0  0;  ...
           0  0  0  0  0  0  0  0  33 34 35; ...
           0  0  0  0  0  0  0  0  36 37 38; ...
           0  0  0  0  0  0  0  0  39 40 41; ...
           0  0  0  0  0  0  0  0  42 43 44];         
for i = 1:numel(anv_dat.cfg.alt_eve.a_vow)
    anv_dat.cfg.alt_eve.a_int(i) = int_num(anv_dat.cfg.alt_eve.a_con(i),anv_dat.cfg.alt_eve.a_vow(i)-100)+200;
end

anv_dat.cfg.alt_eve.a_con(anv_dat.cfg.alt_eve.trialinfo ~= 4) = 0;
anv_dat.cfg.alt_eve.a_vow(anv_dat.cfg.alt_eve.trialinfo ~= 4) = 0;
anv_dat.cfg.alt_eve.a_int(anv_dat.cfg.alt_eve.trialinfo ~= 4) = 0;

eve_for_plt = {[1 2 3 4] [5 6 7 8] [9 10 11 12]; ...
                 [101 102 103 104] [105 106 107 108] [109 110 111]; ...
                 201:216 217:232 233:244}';

% Colors
plt_col     = distinguishable_colors(numel(unique(cat(2,eve_for_plt{[1 4 7 2 5 8]}))));
con_plt_col = plt_col(1:12,:);
vow_plt_col = plt_col(13:23,:);
[ttt2,ttt1] = find(int_num');
for j = 1:numel(ttt1)
    cmb_col(j,:) = sum([vow_plt_col(ttt2(j),:) ; con_plt_col(ttt1(j),:)]) / 2;
end

% Interaction Plots
cfg           = [];
cfg.type      = 'chan';
cfg.plt_lbl   = {'Blk1 Con' 'Blk2 Con' 'Blk3 Con'; ...
                 'Blk1 Vow' 'Blk2 Vow' 'Blk3 Vow'; ...
                 'Blk1 Int' 'Blk2 Int' 'Blk3 Int'}';
cfg.dat       = {anv_dat};
cfg.dat_loc   = [1 1 1 ; ...
                 1 1 1 ; ...
                 1 1 1]';
cfg.lgd       = 0;
cfg.plt_dim   = [3 3];
cfg.alt_eve   = {'a_con' 'a_con' 'a_con'; ...
                 'a_vow' 'a_vow' 'a_vow'; ...
                 'a_int' 'a_int' 'a_int'}';
cfg.eve       = eve_for_plt;
cfg.xlim      = [-0.4 1]; %
cfg.y_lim      = 'auto';
cfg.lnstyle.col_ord = [num2cell(con_plt_col,2)' num2cell(vow_plt_col,2)' num2cell(cmb_col,2)'];
cfg.lnstyle.mrk_ord = repmat({''},1,numel(cfg.lnstyle.col_ord));
cfg.lnstyle.lin_ord = repmat({'-'},1,numel(cfg.lnstyle.col_ord));
cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
cfg.v_lne      = [0 0.450 0.900];
cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
% cfg.stt_dat    = {{anv_dat.cfg.alt_stt.aud_con_blk1 anv_dat.cfg.alt_stt.aud_con_blk1_fdr} {anv_dat.cfg.alt_stt.aud_con_blk2  anv_dat.cfg.alt_stt.aud_con_blk2_fdr} {anv_dat.cfg.alt_stt.aud_con_blk3  anv_dat.cfg.alt_stt.aud_con_blk3_fdr}; ...
%                   {anv_dat.cfg.alt_stt.aud_vow_blk1 anv_dat.cfg.alt_stt.aud_vow_blk1_fdr} {anv_dat.cfg.alt_stt.aud_vow_blk2  anv_dat.cfg.alt_stt.aud_vow_blk2_fdr} {anv_dat.cfg.alt_stt.aud_vow_blk3  anv_dat.cfg.alt_stt.aud_vow_blk3_fdr}; ...
%                   {anv_dat.cfg.alt_stt.aud_int_blk1 anv_dat.cfg.alt_stt.aud_int_blk1_fdr} {anv_dat.cfg.alt_stt.aud_int_blk2  anv_dat.cfg.alt_stt.aud_int_blk2_fdr} {anv_dat.cfg.alt_stt.aud_int_blk3  anv_dat.cfg.alt_stt.aud_int_blk3_fdr}}';
cfg.std_err    = 1;
cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'jpg';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'Aud2x2_smoothed' '/' 'noise'];
cfg.prefix     = 'Aud_2x2_Phoneme';
mmil_ieeg_sensor_plot_v4(cfg)

% Plot Visual ANOVA
% Make interaction events
anv_dat.cfg.alt_eve.v_vow = anv_dat.cfg.alt_eve.v_vow+100;
int_num = [1  2  3  4  0  0  0  0  0  0  0;  ...
           5  6  7  8  0  0  0  0  0  0  0;  ...
           9  10 11 12 0  0  0  0  0  0  0;  ...
           13 14 15 16 0  0  0  0  0  0  0;  ...
           0  0  0  0  17 18 19 20 0  0  0;  ...
           0  0  0  0  21 22 23 24 0  0  0;  ...
           0  0  0  0  25 26 27 28 0  0  0;  ...
           0  0  0  0  29 30 31 32 0  0  0;  ...
           0  0  0  0  0  0  0  0  33 34 35; ...
           0  0  0  0  0  0  0  0  36 37 38; ...
           0  0  0  0  0  0  0  0  39 40 41; ...
           0  0  0  0  0  0  0  0  42 43 44];        
for i = 1:numel(anv_dat.cfg.alt_eve.v_vow)
    anv_dat.cfg.alt_eve.v_int(i) = int_num(anv_dat.cfg.alt_eve.v_con(i),anv_dat.cfg.alt_eve.v_vow(i)-100)+200;
end

anv_dat.cfg.alt_eve.v_con(anv_dat.cfg.alt_eve.trialinfo ~= 3) = 0;
anv_dat.cfg.alt_eve.v_vow(anv_dat.cfg.alt_eve.trialinfo ~= 3) = 0;
anv_dat.cfg.alt_eve.v_int(anv_dat.cfg.alt_eve.trialinfo ~= 3) = 0;

eve_for_plt = {[1 2 3 4] [5 6 7 8] [9 10 11 12]; ...
                 [101 102 103 104] [105 106 107 108] [109 110 111]; ...
                 201:216 217:232 233:244}';

% Colors
plt_col     = distinguishable_colors(numel(unique(cat(2,eve_for_plt{[1 4 7 2 5 8]}))));
con_plt_col = plt_col(1:12,:);
vow_plt_col = plt_col(13:23,:);
[ttt2,ttt1] = find(int_num');
for j = 1:numel(ttt1)
    cmb_col(j,:) = sum([vow_plt_col(ttt2(j),:) ; con_plt_col(ttt1(j),:)]) / 2;
end

% Interaction Plots
cfg           = [];
cfg.type      = 'chan';
cfg.plt_lbl   = {'Blk1 Con' 'Blk2 Con' 'Blk3 Con'; ...
                 'Blk1 Vow' 'Blk2 Vow' 'Blk3 Vow'; ...
                 'Blk1 Int' 'Blk2 Int' 'Blk3 Int'}';
cfg.dat       = {anv_dat};
cfg.dat_loc   = [1 1 1 ; ...
                 1 1 1 ; ...
                 1 1 1]';
cfg.lgd       = 0;
cfg.plt_dim   = [3 3];
cfg.alt_eve   = {'v_con' 'v_con' 'v_con'; ...
                 'v_vow' 'v_vow' 'v_vow'; ...
                 'v_int' 'v_int' 'v_int'}';
cfg.eve       = eve_for_plt;
cfg.xlim      = [-0.4 1]; %
cfg.y_lim      = 'auto';
cfg.lnstyle.col_ord = [num2cell(con_plt_col,2)' num2cell(vow_plt_col,2)' num2cell(cmb_col,2)'];
cfg.lnstyle.mrk_ord = repmat({''},1,numel(cfg.lnstyle.col_ord));
cfg.lnstyle.lin_ord = repmat({'-'},1,numel(cfg.lnstyle.col_ord));
cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
cfg.v_lne      = [0 0.450 0.900];
cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
cfg.stt_dat    = {{anv_dat.cfg.alt_stt.vis_con_blk1 anv_dat.cfg.alt_stt.vis_con_blk1_fdr} {anv_dat.cfg.alt_stt.vis_con_blk2  anv_dat.cfg.alt_stt.vis_con_blk2_fdr} {anv_dat.cfg.alt_stt.vis_con_blk3  anv_dat.cfg.alt_stt.vis_con_blk3_fdr}; ...
                  {anv_dat.cfg.alt_stt.vis_vow_blk1 anv_dat.cfg.alt_stt.vis_vow_blk1_fdr} {anv_dat.cfg.alt_stt.vis_vow_blk2  anv_dat.cfg.alt_stt.vis_vow_blk2_fdr} {anv_dat.cfg.alt_stt.vis_vow_blk3  anv_dat.cfg.alt_stt.vis_vow_blk3_fdr}; ...
                  {anv_dat.cfg.alt_stt.vis_int_blk1 anv_dat.cfg.alt_stt.vis_int_blk1_fdr} {anv_dat.cfg.alt_stt.vis_int_blk2  anv_dat.cfg.alt_stt.vis_int_blk2_fdr} {anv_dat.cfg.alt_stt.vis_int_blk3  anv_dat.cfg.alt_stt.vis_int_blk3_fdr}}';
cfg.std_err    = 1;
cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'jpg';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'Vis2x2_smoothed' '/' 'noise'];
cfg.prefix     = 'Vis_2x2_Phoneme';
mmil_ieeg_sensor_plot_v4(cfg)

%% Mesgarani Approach

%% Initial Machine Learning Approach

%% Scratch
% figure()
% hold on;
% 
% cfg.iCG = 25;
% cfg.iPL = 7; %9;
% iEV     = [9 10 11 12];
% for iE = 1:4;
%     cfg.iEV = iEV(iE);
%     eve_ind = find(cell2mat(cfg.sty_tab(:,1)) == cfg.eve{cfg.iCG}{cfg.iPL}(cfg.iEV));
%     plot(cfg.dat_plt{cfg.iCG}{cfg.iPL}(cfg.iEV,:),'Color',cfg.sty_tab{eve_ind,4});
% end
% 
% figure()
% plot([1 2],[1 2],'Color',cmb_col(9,:),'LineWidth',10); hold on
% plot([1.1 2.1],[1 2],'Color',cmb_col(10,:),'LineWidth',10)
% plot([1.2 2.2],[1 2],'Color',cmb_col(11,:),'LineWidth',10)
% plot([1.3 2.3],[1 2],'Color',cmb_col(12,:),'LineWidth',10)
% 
% figure()
% plot([1 2],[1 2],'Color',cmb_col(3,:),'LineWidth',10); hold on
% plot([1.1 2.1],[1 2],'Color',cmb_col(7,:),'LineWidth',10)
% plot([1.2 2.2],[1 2],'Color',cmb_col(11,:),'LineWidth',10)
% plot([1.3 2.3],[1 2],'Color',cmb_col(15,:),'LineWidth',10)
