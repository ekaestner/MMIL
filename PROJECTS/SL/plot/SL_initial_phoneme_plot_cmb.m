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

%% Initial Phoneme Plots
phn_lbl = {'N' 'B' 'S' 'W' 'M' 'T' 'H' 'F' 'L' 'Y' 'R' 'K' ...
    'IH' 'UH' 'OH' 'AY' 'EE' 'OO' 'AH' 'UR' 'EH' 'AW' 'OY' 'Activity' };
phn_num = [1:6 ; 7:12];
for iC = 1:6
    for iR = 1:2
        
        vis_lay_out{iR,iC}   = find(sll_dat.(sll_dat.data_name{1}).cfg.alt_eve.v_con==phn_num(iR,iC));
        vis_lay_out{iR+2,iC} = find(sll_dat.(sll_dat.data_name{1}).cfg.alt_eve.v_vow==phn_num(iR,iC));
        
        aud_lay_out{iR,iC}   = find(sll_dat.(sll_dat.data_name{1}).cfg.alt_eve.a_con==phn_num(iR,iC));
        aud_lay_out{iR+2,iC} = find(sll_dat.(sll_dat.data_name{1}).cfg.alt_eve.a_vow==phn_num(iR,iC));
        
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
    cfg.v_lne      = [0 0.450 0.900];
    cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
    cfg.std_err    = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'VisPhn' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = 'Vis_Grapheme';
    cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
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
    cfg.xlim      = [-0.4 1]; %
    cfg.ylink     = {[4 6]}; %
    cfg.lnstyle.col_ord = {rgb('blue') rgb('cyan')}; %
    cfg.lnstyle.lin_ord = {'-' '-'}; %
    cfg.cnd_nme         = {'Aud' 'AudNse'}; %
    cfg.leg_pos = [1 2]; %
    cfg.alt_lbl = 'repair_macro_ejk1st_meso_label';
    cfg.trl        = aud_lay_out;
    cfg.v_lne      = [0 0.450 0.900];
    cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
    cfg.std_err    = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'AudPhn' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = 'Aud_Phoneme';
    cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
    mmil_ieeg_sensor_plot_v4(cfg)
    
end

%% Plot with different events
for iD = 1:numel(sll_dat.data_name)
    
    % Vis
    cfg           = [];
    cfg.type      = 'chan';
    cfg.plt_lbl   = phn_lbl;  
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [4 6];
    cfg.alt_eve   = {'trialinfo'};
    cfg.eve       = [1 2 3 4];
    cfg.xlim      = [-0.4 1]; %
    cfg.ylink     = {[4 6]}; %
    cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
    cfg.lnstyle.lin_ord = {'-' '-' '-' '-'};
    cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'  'AudNse'};
    cfg.leg_pos = [1 2 5 6];
    cfg.alt_lbl = 'repair_macro_ejk1st_meso_label';
    cfg.trl        = vis_lay_out;
    cfg.v_lne      = [0 0.450 0.900];
    cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
    cfg.std_err    = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'VisPhnAllType' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = 'Vis_Grapheme';
    cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
    mmil_ieeg_sensor_plot_v4(cfg)
    
    % Aud
    cfg           = [];
    cfg.type      = 'chan';
    cfg.plt_lbl   = phn_lbl;
    cfg.dat       = {sll_dat.(sll_dat.data_name{iD})};
    cfg.lgd       = 0;
    cfg.plt_dim   = [4 6];
    cfg.alt_eve   = {'trialinfo'};
    cfg.eve       = [1 2 3 4];
    cfg.xlim      = [-0.4 1]; %
    cfg.ylink     = {[4 6]}; %
    cfg.lnstyle.col_ord = {rgb('green') rgb('yellow') rgb('magenta') rgb('cyan')};
    cfg.lnstyle.lin_ord = {'-' '-' '-' '-'};
    cfg.cnd_nme         = {'Mtc' 'MssMtc' 'VisNse'  'AudNse'};
    cfg.leg_pos = [1 2 5 6];
    cfg.alt_lbl = 'repair_macro_ejk1st_meso_label';
    cfg.trl        = aud_lay_out;
    cfg.v_lne      = [0 0.450 0.900];
    cfg.v_lne_col  = {rgb('red') rgb('blue') rgb('black')};
    cfg.std_err    = 1;
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'jpg';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data_cmb/initial_plot' '/'  subj '/' 'AudPhnAllType' '/' sll_dat.data_name{iD}(end-2:end)];
    cfg.prefix     = 'Aud_Phoneme';
    cfg.alt_lbl    = 'repair_macro_ejk1st_meso_label';
    mmil_ieeg_sensor_plot_v4(cfg)
    
end