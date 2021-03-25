clear; clc;

sbj_num = 3;
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
cfg.file = [outpath '/' subj '_overall_data.mat'];
sll_dat     = ft_func([],cfg);

ttt = fieldnames(sll_dat.NY540_SL_Day3_1_1_Clin1_hgp.cfg.alt_stt);
ttt = ttt(8:end);
for iF = 1:numel(ttt)
    sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.(ttt{iF}).label = cellfun(@(x) x(5:end-4),sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.(ttt{iF}).label,'uni',0);
    jjj = cellfun(@(x) strfind(x,'_'),sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.(ttt{iF}).label,'uni',0);
    sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.(ttt{iF}).label = cellfun(@(x,y) x([1:y-1 y+1:end]),sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.(ttt{iF}).label,jjj,'uni',0);
end

fprintf('Starting work on %s \n',subj)

%% Initial Plots
% for iD = 1:numel(sll_dat.data_name)
%     %% Plot Both Chosen & Not Chosen 4-event Channels
%     % Plot Visual Events
%     cfg = [];
%     cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
%     cfg.lgd       = 0;
%     cfg.plt_dim   = [5 5];
%     
%     cfg.cmb = 1; 
%     cfg.cmp_ind = [2 4 6:7];
%     cfg.cmb_typ = {[6 7]};
%     cfg.cmb_lbl = {'Ely_Vis_Mtc'};    
%     cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; 
%     [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
%     
%     cfg.chn_grp   = grp_chn{iD};
%     cfg.alt_eve         = {'trialinfo'};
%     cfg.eve             = [1 2 3 4];
%     cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
%     cfg.lnstyle.lin_ord = {'-' '-' '-' '-'};
%     cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'  'AudNse'};
%     cfg.leg_pos = [1 2 5 6];
%     cfg.alt_lbl = 'repair_macro_ejk1st_meso_label';
%     cfg.plt_shf = 1;
%     cfg.y_lim   = 'auto';
%     cfg.stt_dat = grp_stt{iD};
%     cfg.stt_col = grp_stt_col{iD};
%     cfg.std_err = 1;
%     cfg.print      = 1;
%     cfg.nofig      = 1;
%     cfg.print_type = 'jpg';
%     cfg.outdir     = [outpath '/' 'initial_plot' '/' subj '/' sll_dat.data_name{iD}(end-2:end)];
%     cfg.prefix     = strcat(sll_dat.data_name{iD},{'_'},pre_fix{iD});
%     cfg.v_lne      = [0 0.450];
%     mmil_ieeg_sensor_plot_v4(cfg)
%     
% end
% 
% %% Early Visual!Noise Effects  
% % Data
% sig_chn_ind = 2;
% for iD = 1:numel(sll_dat.data_name)
%     cfg = [];
%     cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
%     cfg.lgd       = 0;
%     cfg.plt_dim   = [5 5];
%     
%     cfg.cmb = 1; 
%     cfg.cmp_ind = sig_chn_ind; 
%     cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; 
%     [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
%     
%     cfg.chn_grp   = grp_chn{iD};
%     cfg.alt_eve         = {'trialinfo'};
%     cfg.eve             = [1 2 3];
%     cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta')};
%     cfg.lnstyle.lin_ord = {'-' '-' '-'};
%     cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'};
%     cfg.leg_pos = [1 2 5];
%     cfg.alt_lbl = 'repair_macro_ejk1st_meso_label';
%     cfg.plt_shf = 1;
%     cfg.y_lim   = 'auto';
%     cfg.stt_dat = grp_stt{iD};
%     cfg.stt_col = grp_stt_col{iD};
%     cfg.std_err = 1;
%     cfg.print      = 1;
%     cfg.nofig      = 1;
%     cfg.print_type = 'jpg';
%     cfg.outdir     = [outpath '/' 'initial_plot' '/' subj '/' 'VisNse_' sll_dat.data_name{iD}(end-2:end)];
%     cfg.prefix     = strcat(sll_dat.data_name{iD},{'_VisNse_'},{'_'},pre_fix{iD});
%     cfg.v_lne      = [0 0.450 0.900];
%     cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
%     mmil_ieeg_sensor_plot_v4(cfg)
%     
% end
% 
% % Location
% for iD = 1:numel(sll_dat.data_name)
%     chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
%     
%     cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}];
%     sig_chn = mmil_readtext(cfg.sig_fle,[',']);
%     tmp = find(strcmpi(sig_chn(1,:),'stars'));
%     sig_chn = cell2mat(sig_chn(2:end,2:tmp-2));
% 
%     sig_col = 2;
%     
%     cfg = [];
%     
%     cfg.elec_text = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/electrode_location' '/' subj]);
%     cfg.pial_mat  = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/pial' '/' subj]);
%     cfg.hms       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/hemi' '/' subj]);
%     
%     cfg.sel_lbl   = {'VisualResponsive'};
%     cfg.hem       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/hemi/' subj]);
%     cfg.lab_pos   = 1;
%     cfg.sel_ele   = cellfun(@(x) sll_dat.(sll_dat.data_name{iD}).cfg.alt_lab.repair_macro_ejk1st_meso_label(x),cellfun(@find,num2cell(sig_chn(:,sig_chn_ind),1),'Uni',0),'Uni',0);
%     cfg.col       = {rgb('red')};
%     
%     cfg.sve_img   = 'jpg';
%     cfg.sve_loc   = [outpath '/' 'initial_plot' '/' subj '/' 'VisNse_' sll_dat.data_name{iD}(end-2:end)];
%     cfg.sve_pre   = [subj '_VisNse'];
%     
%     mmil_ieeg_sensor_location_plot_v1(cfg);
%     
% end
% 
% %% Early Auditory!Noise Effects
% % Data
% sig_chn_ind = 4;
% for iD = 1:numel(sll_dat.data_name)
%     cfg = [];
%     cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
%     cfg.lgd       = 0;
%     cfg.plt_dim   = [5 5];
%     
%     cfg.cmb = 1; 
%     cfg.cmp_ind = sig_chn_ind; 
%     cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; 
%     [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
%     
%     cfg.chn_grp   = grp_chn{iD};
%     cfg.alt_eve         = {'trialinfo'};
%     cfg.eve             = [1 2 4];
%     cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('cyan')};
%     cfg.lnstyle.lin_ord = {'-' '-' '-'};
%     cfg.cnd_nme         = {'Mtc' 'MssMtc' 'AudNse'};
%     cfg.leg_pos = [1 2 5];
%     cfg.alt_lbl = 'repair_macro_ejk1st_meso_label';
%     cfg.plt_shf = 1;
%     cfg.y_lim   = 'auto';
%     cfg.stt_dat = grp_stt{iD};
%     cfg.stt_col = grp_stt_col{iD};
%     cfg.std_err = 1;
%     cfg.print      = 1;
%     cfg.nofig      = 1;
%     cfg.print_type = 'jpg';
%     cfg.outdir     = [outpath '/' 'initial_plot' '/' subj '/' 'AudNse_' sll_dat.data_name{iD}(end-2:end)];
%     cfg.prefix     = strcat(sll_dat.data_name{iD},{'_AudNse_'},{'_'},pre_fix{iD});
%     cfg.v_lne      = [0 0.450 0.900];
%     cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
%     mmil_ieeg_sensor_plot_v4(cfg)
%     
% end
% 
% % Location
% for iD = 1:numel(sll_dat.data_name)
%     chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
%     
%     cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}];
%     sig_chn = mmil_readtext(cfg.sig_fle,[',']);
%     tmp = find(strcmpi(sig_chn(1,:),'stars'));
%     sig_chn = cell2mat(sig_chn(2:end,2:tmp-2));
% 
%     sig_col = 4;
%     
%     cfg = [];
%     
%     cfg.elec_text = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/electrode_location' '/' subj]);
%     cfg.pial_mat  = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/pial' '/' subj]);
%     cfg.hms       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/hemi' '/' subj]);
%     
%     cfg.sel_lbl   = {'VisualResponsive'};
%     cfg.hem       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/hemi/' subj]);
%     cfg.lab_pos   = 1;
%     cfg.sel_ele   = cellfun(@(x) sll_dat.(sll_dat.data_name{iD}).cfg.alt_lab.repair_macro_ejk1st_meso_label(x),cellfun(@find,num2cell(sig_chn(:,sig_chn_ind),1),'Uni',0),'Uni',0);
%     cfg.col       = {rgb('blue')};
%     
%     cfg.sve_img   = 'jpg';
%     cfg.sve_loc   = [outpath '/' 'initial_plot' '/' subj '/' 'AudNse_' sll_dat.data_name{iD}(end-2:end)];
%     cfg.sve_pre   = [subj '_AudNse'];
%     
%     mmil_ieeg_sensor_location_plot_v1(cfg);
%     
% end
% 
% %% Early & Late Match Effects
% % Data
% sig_chn_ind = [6 7];
% for iD = 1:numel(sll_dat.data_name)
%     cfg = [];
%     cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
%     cfg.lgd       = 0;
%     cfg.plt_dim   = [5 5];
%     
%     cfg.cmb = 1; 
%     cfg.cmp_ind = sig_chn_ind; 
%     cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; 
%     [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
%     
%     cfg.chn_grp   = grp_chn{iD};
%     cfg.alt_eve         = {'trialinfo'};
%     cfg.eve             = [1 2];
%     cfg.lnstyle.col_ord = {rgb('green') rgb('yellow')};
%     cfg.lnstyle.lin_ord = {'-' '-'};
%     cfg.cnd_nme         = {'Mtc' 'MssMtc'};
%     cfg.leg_pos = [1 2];
%     cfg.alt_lbl = 'repair_macro_ejk1st_meso_label';
%     cfg.plt_shf = 1;
%     cfg.y_lim   = 'auto';
%     cfg.stt_dat = grp_stt{iD};
%     cfg.stt_col = grp_stt_col{iD};
%     cfg.std_err = 1;
%     cfg.print      = 1;
%     cfg.nofig      = 1;
%     cfg.print_type = 'jpg';
%     cfg.outdir     = [outpath '/' 'initial_plot' '/' subj '/' 'Match_' sll_dat.data_name{iD}(end-2:end)];
%     cfg.prefix     = strcat(sll_dat.data_name{iD},{'_Match_'},{'_'},pre_fix{iD});
%     cfg.v_lne      = [0 0.450 0.900];
%     cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
%     mmil_ieeg_sensor_plot_v4(cfg)
%     
% end
% 
% % Location
% for iD = 1:numel(sll_dat.data_name)
%     chn_loc = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical' '/' 'electrode_location' '/' subj]);
%     
%     cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}];
%     sig_chn = mmil_readtext(cfg.sig_fle,[',']);
%     tmp = find(strcmpi(sig_chn(1,:),'stars'));
%     sig_chn = cell2mat(sig_chn(2:end,2:tmp-2));
% 
%     sig_col = 6:7;
%     
%     cfg = [];
%     
%     cfg.elec_text = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/electrode_location' '/' subj]);
%     cfg.pial_mat  = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/pial' '/' subj]);
%     cfg.hms       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/hemi' '/' subj]);
%     
%     cfg.sel_lbl   = {'Early Match' 'Late Match'};
%     cfg.hem       = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/hemi/' subj]);
%     cfg.lab_pos   = [1 2];
%     cfg.sel_ele   = cellfun(@(x) sll_dat.(sll_dat.data_name{iD}).cfg.alt_lab.repair_macro_ejk1st_meso_label(x),cellfun(@find,num2cell(sig_chn(:,sig_chn_ind),1),'Uni',0),'Uni',0);
%     cfg.col       = {rgb('light green') rgb('dark green')};
%     
%     cfg.sve_img   = 'jpg';
%     cfg.sve_loc   = [outpath '/' 'initial_plot' '/' subj '/' 'Match_' sll_dat.data_name{iD}(end-2:end)];
%     cfg.sve_pre   = [subj '_Match'];
%     
%     mmil_ieeg_sensor_location_plot_v1(cfg);
%     
% end

%% Channel-by-Channel, 3x3 (VisNse/AudNse/Mtc x LFP/HGP/Theta)
cfg           = [];
cfg.type      = 'chan';
cfg.plt_lbl   = {'LFP VisNse' 'HGP VisNse' 'THT VisNse'; ...
                 'LFP AudNse' 'HGP AudNse' 'THT AudNse'; ...
                 'LFP Match' 'HGP Match' 'THT Match'}';
cfg.dat       = {sll_dat.(sll_dat.data_name{1}) sll_dat.(sll_dat.data_name{2}) sll_dat.(sll_dat.data_name{3})};
cfg.dat_loc   = [1 1 1 ...
                 2 2 2 ...
                 3 3 3]';
cfg.lgd       = 0;
cfg.plt_dim   = [3 3];
cfg.alt_eve   = {'trialinfo' 'trialinfo' 'trialinfo' ...
                 'trialinfo' 'trialinfo' 'trialinfo'  ...
                 'trialinfo' 'trialinfo' 'trialinfo'}';
cfg.eve       = {[1 2 3] [1 2 3] [1 2 3]; ...
                 [1 2 4] [1 2 4] [1 2 4];...
                 [1 2]   [1 2]   [1 2]}';
cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
cfg.stt_dat = {{sll_dat.(sll_dat.data_name{1}).cfg.alt_stt.vis_nse_stt_msk} {sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.vis_nse_stt_msk} {sll_dat.(sll_dat.data_name{3}).cfg.alt_stt.vis_nse_stt_msk} ; ...
               {sll_dat.(sll_dat.data_name{1}).cfg.alt_stt.aud_nse_stt_msk} {sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.aud_nse_stt_msk} {sll_dat.(sll_dat.data_name{3}).cfg.alt_stt.aud_nse_stt_msk} ; ...
               {sll_dat.(sll_dat.data_name{1}).cfg.alt_stt.vis_mtc_stt_msk} {sll_dat.(sll_dat.data_name{2}).cfg.alt_stt.vis_mtc_stt_msk} {sll_dat.(sll_dat.data_name{3}).cfg.alt_stt.vis_mtc_stt_msk}}';
cfg.stt_col = {{ft_stt_col(rgb('magenta'))} {ft_stt_col(rgb('magenta'))} {ft_stt_col(rgb('magenta'))} ; ...
               {ft_stt_col(rgb('cyan'))}    {ft_stt_col(rgb('cyan'))}    {ft_stt_col(rgb('cyan'))} ; ...
               {ft_stt_col(rgb('green'))}   {ft_stt_col(rgb('green'))}   {ft_stt_col(rgb('green'))}}';
if strcmpi(subj,'NY439_SL'); cfg.alt_lbl = 'repair_macro_ejk1st_meso_label'; end;
cfg.y_lim      = 'auto';
cfg.std_err    = 1;
cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'jpg';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'Overall'];
cfg.prefix     = 'OverallPlots';
cfg.v_lne      = [0 0.450];
mmil_ieeg_sensor_plot_v4(cfg)

%% Phoneme Plots
% Make interaction events
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
for iD = 1:numel(sll_dat.data_name)
    for i = 1:numel(sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_vow)
        sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_int(i) = int_num(sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_con(i),sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_vow(i)-100)+200;
    end
end

sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_con_rel = sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_con;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_vow_rel = sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_vow;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_int_rel = sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_int;

sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_con_rel(sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.trialinfo == 4) = 0;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_vow_rel(sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.trialinfo == 4) = 0;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_int_rel(sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.trialinfo == 4) = 0;

eve_for_plt = {[1 2 3 4] [5 6 7 8] [9 10 11 12]; ...
                 [101 102 103 104] [105 106 107 108] [109 110 111 112]; ...
                 201:216 217:232 233:248};

% Colors
plt_col     = distinguishable_colors(numel(unique(cat(2,eve_for_plt{[1 4 7 2 5 8]}))));
con_plt_col = plt_col(1:12,:);
vow_plt_col = plt_col(13:24,:);
[ttt2,ttt1] = find(int_num');
for j = 1:numel(ttt1)
    cmb_col(j,:) = sum([vow_plt_col(ttt2(j),:) ; con_plt_col(ttt1(j),:)]) / 2;
end

% Interaction Plots
cfg           = [];
cfg.type      = 'chan';
cfg.plt_lbl   = {'Blk1 Con' 'Blk2 Con' 'Blk3 Con'; ...
                 'Blk1 Vow' 'Blk2 Vow' 'Blk3 Vow'; ...
                 'Blk1 Int' 'Blk2 Int' 'Blk3 Int'};
cfg.dat       = {sll_dat.(sll_dat.data_name{2})};
cfg.dat_loc   = [1 1 1 ; ...
                 1 1 1 ; ...
                 1 1 1]';
cfg.lgd       = 0;
cfg.plt_dim   = [3 3];
cfg.alt_eve   = {'a_con_rel' 'a_con_rel' 'a_con_rel'; ...
                 'a_vow_rel' 'a_vow_rel' 'a_vow_rel'; ...
                 'a_int_rel' 'a_int_rel' 'a_int_rel'};
cfg.eve       = eve_for_plt;
cfg.xlim      = [-0.4 1]; %
cfg.y_lim      = 'maxmin';
cfg.lnstyle.col_ord = [num2cell(con_plt_col,2)' num2cell(vow_plt_col,2)' num2cell(cmb_col,2)'];
cfg.lnstyle.mrk_ord = repmat({''},1,numel(cfg.lnstyle.col_ord));
cfg.lnstyle.lin_ord = repmat({'-'},1,numel(cfg.lnstyle.col_ord));
if strcmpi(subj,'NY439_SL'); cfg.alt_lbl = 'repair_macro_ejk1st_meso_label'; end;
cfg.v_lne      = [0 0.450 0.900];
cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
cfg.stt_dat    = {{'aud_con_blk1' 'aud_con_blk1_fdr'} {'aud_con_blk2'  'aud_con_blk2_fdr'} {'aud_con_blk3'  'aud_con_blk3_fdr'}; ...
                  {'aud_vow_blk1' 'aud_vow_blk1_fdr'} {'aud_vow_blk2'  'aud_vow_blk2_fdr'} {'aud_vow_blk3'  'aud_vow_blk3_fdr'}; ...
                  {'aud_int_blk1' 'aud_int_blk1_fdr'} {'aud_int_blk2'  'aud_int_blk2_fdr'} {'aud_int_blk3'  'aud_int_blk3_fdr'}};
cfg.std_err    = 1;
cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'jpg';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'Aud2x2_smoothed' '/'];
cfg.prefix     = 'Aud_2x2_Phoneme';
mmil_ieeg_sensor_plot_v5(cfg)

%% Noise plots
% Make interaction events
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
for iD = 1:numel(sll_dat.data_name)
    for i = 1:numel(sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_vow)
        sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_int(i) = int_num(sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_con(i),sll_dat.(sll_dat.data_name{iD}).cfg.alt_eve.a_vow(i)-100)+200;
    end
end

sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_con_nse = sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_con;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_vow_nse = sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_vow;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_int_nse = sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_int;

sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_con_nse(sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.trialinfo ~= 4) = 0;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_vow_nse(sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.trialinfo ~= 4) = 0;
sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.a_int_nse(sll_dat.(sll_dat.data_name{2}).cfg.alt_eve.trialinfo ~= 4) = 0;

eve_for_plt = {[1 2 3 4] [5 6 7 8] [9 10 11 12]; ...
                 [101 102 103 104] [105 106 107 108] [109 110 111 112]; ...
                 201:216 217:232 233:248};

% Colors
plt_col     = distinguishable_colors(numel(unique(cat(2,eve_for_plt{[1 4 7 2 5 8]}))));
con_plt_col = plt_col(1:12,:);
vow_plt_col = plt_col(13:24,:);
[ttt2,ttt1] = find(int_num');
for j = 1:numel(ttt1)
    cmb_col(j,:) = sum([vow_plt_col(ttt2(j),:) ; con_plt_col(ttt1(j),:)]) / 2;
end

% Interaction Plots
cfg           = [];
cfg.type      = 'chan';
cfg.plt_lbl   = {'Blk1 Con' 'Blk2 Con' 'Blk3 Con'; ...
                 'Blk1 Vow' 'Blk2 Vow' 'Blk3 Vow'; ...
                 'Blk1 Int' 'Blk2 Int' 'Blk3 Int'};
cfg.dat       = {sll_dat.(sll_dat.data_name{2})};
cfg.dat_loc   = [1 1 1 ; ...
                 1 1 1 ; ...
                 1 1 1];
cfg.lgd       = 0;
cfg.plt_dim   = [3 3];
cfg.alt_eve   = {'a_con_nse' 'a_con_nse' 'a_con_nse'; ...
                 'a_vow_nse' 'a_vow_nse' 'a_vow_nse'; ...
                 'a_int_nse' 'a_int_nse' 'a_int_nse'};
cfg.eve       = eve_for_plt;
cfg.xlim      = [-0.4 1]; %
cfg.y_lim      = 'maxmin';
cfg.lnstyle.col_ord = [num2cell(con_plt_col,2)' num2cell(vow_plt_col,2)' num2cell(cmb_col,2)'];
cfg.lnstyle.mrk_ord = repmat({''},1,numel(cfg.lnstyle.col_ord));
cfg.lnstyle.lin_ord = repmat({'-'},1,numel(cfg.lnstyle.col_ord));
if strcmpi(subj,'NY439_SL'); cfg.alt_lbl = 'repair_macro_ejk1st_meso_label'; end;
cfg.v_lne      = [0 0.450 0.900];
cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
% cfg.stt_dat    = {{'aud_con_blk1' 'aud_con_blk1_fdr'} {'aud_con_blk2'  'aud_con_blk2_fdr'} {'aud_con_blk3'  'aud_con_blk3_fdr'}; ...
%                   {'aud_vow_blk1' 'aud_vow_blk1_fdr'} {'aud_vow_blk2'  'aud_vow_blk2_fdr'} {'aud_vow_blk3'  'aud_vow_blk3_fdr'}; ...
%                   {'aud_int_blk1' 'aud_int_blk1_fdr'} {'aud_int_blk2'  'aud_int_blk2_fdr'} {'aud_int_blk3'  'aud_int_blk3_fdr'}};
cfg.std_err    = 1;
cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'jpg';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'Aud2x2_smoothed_Noise' '/'];
cfg.prefix     = 'Aud_2x2_Phoneme';
mmil_ieeg_sensor_plot_v5(cfg)









