clear; clc;

sbj_num = 1;
ovr_wrt = 1;

% Setting up Variables
subj  = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/clerical/subjects');
subj  = subj{sbj_num};

outpath = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/';

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

fprintf('Starting work on %s \n',subj)

%% Initial Plots
for iD = 1:numel(sll_dat.data_name)
    %% Plot Both Chosen & Not Chosen 4-event Channels
    % Plot Visual Events
    cfg = [];
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [5 5];
    cfg.cmb = 1; cfg.cmp_ind = 1; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
    cfg.chn_grp   = grp_chn{iD};
    cfg.alt_eve         = {'trialinfo'};
    cfg.eve             = [1 2 3 4];
    cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
    cfg.lnstyle.lin_ord = {'-' '-' '-' '-'};
    cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'  'AudNse'};
    cfg.leg_pos = [1 2 5 6];
    cfg.alt_lbl = 'label';
    cfg.plt_shf = 1;
    cfg.y_lim   = 'auto';
    cfg.stt_dat = grp_stt{iD};
    cfg.stt_col = grp_stt_col{iD};
    cfg.std_err = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/initial_plot' '/'  subj '/' 'VisOverall' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = strcat(sll_dat.data_name{iD},{'_'},pre_fix{iD});
    cfg.v_lne      = 0.450;
    mmil_ieeg_sensor_plot_v4(cfg)
    
    % Plot Auditory Events
    cfg = [];
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [5 5];
    cfg.cmb = 1; cfg.cmp_ind = 4; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
    cfg.chn_grp   = grp_chn{iD};
    cfg.alt_eve         = {'trialinfo'};
    cfg.eve             = [11 12 13 14];
    cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
    cfg.lnstyle.lin_ord = {'-' '-' '-' '-'};
    cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'  'AudNse'};
    cfg.leg_pos = [1 2 5 6];
    cfg.alt_lbl = 'label';
    cfg.plt_shf = 1;
    cfg.y_lim   = 'auto';
    cfg.stt_dat = grp_stt{iD};
    cfg.stt_col = grp_stt_col{iD};
    cfg.std_err = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/initial_plot' '/'  subj '/' 'AudOverall' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = strcat(sll_dat.data_name{iD},{'_'},pre_fix{iD});
    mmil_ieeg_sensor_plot_v4(cfg)
    
    %% Plot Both Chosen & Not Chosen mm/nse events Channels
    % Plot Visual Events
    cfg = [];
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [5 5];
    cfg.cmb = 1; cfg.cmp_ind = 2:3; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
    cfg.chn_grp   = grp_chn{iD};
    cfg.alt_eve         = {'trialinfo'};
    cfg.eve             = [1 2 3 4];
    cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
    cfg.lnstyle.lin_ord = {'-' '-' '-' '-'};
    cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'  'AudNse'};
    cfg.leg_pos = [1 2 5 6];
    cfg.alt_lbl = 'label';
    cfg.plt_shf = 1;
    cfg.y_lim   = 'auto';
    cfg.stt_dat = grp_stt{iD};
    cfg.stt_col = grp_stt_col{iD};
    cfg.std_err = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/initial_plot' '/'  subj '/' 'VisMatchMismatch' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = strcat(sll_dat.data_name{iD},{'_'},pre_fix{iD});
    mmil_ieeg_sensor_plot_v4(cfg)
    
    % Plot Auditory Events
    cfg = [];
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [5 5];
    cfg.cmb = 1; cfg.cmp_ind = 5:6; cfg.sig_fle = [cln_dir '/' 'sig_chn' '/' subj '/' sll_dat.data_name{iD}]; [grp_chn{iD},pre_fix{iD},grp_stt{iD},grp_stt_col{iD}] = plot_channel_setup(cfg,sll_dat.(sll_dat.data_name{iD}));
    cfg.chn_grp   = grp_chn{iD};
    cfg.alt_eve         = {'trialinfo'};
    cfg.eve             = [11 12 13 14];
    cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
    cfg.lnstyle.lin_ord = {'-' '-' '-' '-'};
    cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'  'AudNse'};
    cfg.leg_pos = [1 2 5 6];
    cfg.alt_lbl = 'label';
    cfg.plt_shf = 1;
    cfg.y_lim   = 'auto';
    cfg.stt_dat = grp_stt{iD};
    cfg.stt_col = grp_stt_col{iD};
    cfg.std_err = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/initial_plot' '/'  subj '/' 'AudMatchMismatch' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = strcat(sll_dat.data_name{iD},{'_'},pre_fix{iD});
    mmil_ieeg_sensor_plot_v4(cfg)
end

%% Initial Phoneme Plots
phn_lbl = {'N' 'B' 'S' 'W' 'M' 'T' 'H' 'F' 'L' 'Y' 'R' 'K' ...
           'IH' 'UH' 'OH' 'AY' 'EE' 'OO' 'AH' 'UR' 'EH' 'AW' 'OY' 'Activity' };
phn_num = [1:6 ; 7:12];
for iC = 1:6
    for iR = 1:2
        
        vis_lay_out{iR,iC}   = find(sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.v_con==phn_num(iR,iC));
        vis_lay_out{iR+2,iC} = find(sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.v_vow==phn_num(iR,iC));
        
        aud_lay_out{iR,iC}   = find(sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.a_con==phn_num(iR,iC));
        aud_lay_out{iR+2,iC} = find(sll_dat.(sll_dat.data_name{iST}).cfg.alt_eve.a_vow==phn_num(iR,iC));

    end
end

vis_lay_out{iR+2,iC} = [1:numel(sll_dat.(sll_dat.data_name{1}).trialinfo)]';
aud_lay_out{iR+2,iC} = [1:numel(sll_dat.(sll_dat.data_name{1}).trialinfo)]';

for iD = 1:numel(sll_dat.data_name)
    
    % Vis
    cfg           = [];
    cfg.type      = 'chan';
    cfg.plt_lbl   = phn_lbl;  
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [4 6];
    cfg.alt_eve   = {'vis_tot_nse'};
    cfg.eve    = [101 102];
    cfg.xlim      = [-0.4 0.6]; %
    cfg.ylink     = {[4 6]}; %
    cfg.lnstyle.col_ord = {rgb('red') rgb('magenta')};
    cfg.lnstyle.lin_ord = {'-' '-'};
    cfg.cnd_nme         = {'Vis' 'VisNse'};
    cfg.leg_pos = [1 2];
    cfg.alt_lbl = 'label';
    cfg.trl        = vis_lay_out;
    cfg.std_err    = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/initial_plot' '/'  subj '/' 'VisPhn' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = 'Vis_Grapheme';
    mmil_ieeg_sensor_plot_v4(cfg)
    
    % Aud
    cfg           = [];
    cfg.type      = 'chan';
    cfg.plt_lbl   = phn_lbl;
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [4 6];
    cfg.alt_eve   = {'aud_tot_nse'}; %
    cfg.eve       = [111 112]; %
    cfg.xlim      = [-0.4 0.6]; %
    cfg.ylink     = {[4 6]}; %
    cfg.lnstyle.col_ord = {rgb('blue') rgb('cyan')}; %
    cfg.lnstyle.lin_ord = {'-' '-'}; %
    cfg.cnd_nme         = {'Aud' 'AudNse'}; %
    cfg.leg_pos = [1 2]; %
    cfg.alt_lbl = 'label';
    cfg.trl        = aud_lay_out;
    cfg.std_err    = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data/initial_plot' '/'  subj '/' 'AudPhn' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = 'Aud_Phoneme';
    mmil_ieeg_sensor_plot_v4(cfg)
    
end











