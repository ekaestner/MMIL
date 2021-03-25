function iSASZ_Figure9

fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
fcfg.sbj_nme = 'MG49_SA_SZ_utah';

%% Neuroport Location
% Find Close electrodes
neu_prt_loc = mmil_readtext([fcfg.clr_fld '/' 'electrode_location_files' '/' 'MG049_SA_SZ' '/' 'output' '/' 'MG049_SA_SZ_lhs_ecog']);

neu_loc = cell2mat(neu_prt_loc(end,2:4));
for iR = 1:size(neu_prt_loc,1)-1 
    dst_hld(iR) = pdist([neu_loc ; cell2mat(neu_prt_loc(iR,2:4))]);
end
[~,dst_rnk] = sort(dst_hld);
neu_prt_loc(dst_rnk,1);

% Location Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'MG049_SA_SZ' '/'  'surf' '/' 'lh.pial']                [fcfg.clr_fld '/' 'electrode_location_files' '/' 'MG049_SA_SZ' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'MG049_SA_SZ' '/'  'output' '/' 'MG049_SA_SZ_lhs_ecog'] [fcfg.clr_fld '/' 'electrode_location_files' '/' 'MG049_SA_SZ' '/'  'output' '/' 'MG049_SA_SZ_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = {{'Neuroport'}}; % { { 'GR34' 'GR35' 'GR41' 'GR42' 'GR43' 'GR44' 'GR45' 'GR46' 'GR47' 'GR48' 'GR53' 'GR54' 'GR55' 'GR56' } {'Neuroport'}}; % 'GR27' 'GR28'
cfg.all_ele   = {'Neuroport'}; % {   'GR34' 'GR35' 'GR41' 'GR42' 'GR43' 'GR44' 'GR45' 'GR46' 'GR47' 'GR48' 'GR53' 'GR54' 'GR55' 'GR56'    'Neuroport'}; % 'GR27' 'GR28' 'GR33' 'GR34' 'GR35' 'GR36' 'GR41' 'GR42' 'GR43' 'GR44'

cfg.sel_lbl = {'Neuroport'};
cfg.col     = {rgb('green')};
cfg.nsl_col = {rgb('white')};

cfg.sep_str         = ',';

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'location_plot'];
cfg.sve_pre   = ['middle_pic_auditory'];
cfg.sep_str   = [','];

cfg.rad       = 2.25;

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

%% Additional Channel Plots
% Load Data
cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' 'MG049_SA_SZ' '_overall_data.mat'];
ecg_dat  = ft_func([],cfg);

ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.plt_eve(801:1200) = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.vis_new_bse(801:1200);
ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.plt_eve(1:801)    = ecg_dat.(ecg_dat.data_name{2}).cfg.alt_eve.aud_new_bse(1:801);

% Find Channels
ecg_stg = zeros(8,8);

ecg_stg(1,2) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR34'));
ecg_stg(1,3) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR35'));

ecg_stg(2,1) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR41'));
ecg_stg(2,2) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR42'));
ecg_stg(2,3) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR43'));
ecg_stg(2,4) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR44'));
ecg_stg(2,5) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR45'));
ecg_stg(2,6) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR46'));
ecg_stg(2,7) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR47'));
ecg_stg(2,8) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR48'));

ecg_stg(3,5) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR53'));
ecg_stg(3,6) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR54'));
ecg_stg(3,7) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR55'));
ecg_stg(3,8) = find(strcmpi(ecg_dat.(ecg_dat.data_name{1}).cfg.alt_lab.label,'GR56'));

% 
plt_lbl = cell(8,8);

plt_lbl{1,2} = 'GR34';
plt_lbl{1,3} = 'GR35';

plt_lbl{2,1} = 'GR41';
plt_lbl{2,2} = 'GR42';
plt_lbl{2,3} = 'GR43';
plt_lbl{2,4} = 'GR44';
plt_lbl{2,5} = 'GR45';
plt_lbl{2,6} = 'GR46';
plt_lbl{2,7} = 'GR47';
plt_lbl{2,8} = 'GR48';

plt_lbl{3,5} = 'GR53';
plt_lbl{3,6} = 'GR54';
plt_lbl{3,7} = 'GR55';
plt_lbl{3,8} = 'GR56';

% Make Plot
cfg = [];

cfg.dat = {ecg_dat.(ecg_dat.data_name{2})};

cfg.type      = 'macro';

cfg.chn_grp = { ecg_stg };
cfg.dat_loc = cfg.chn_grp{1}>0;
cfg.plt_dim = size(cfg.chn_grp{1});

cfg.eve       = { repmat({[102 202]},8,8) };

cfg.alt_eve   = { repmat({'plt_eve'},8,8) };

cfg.plt_lbl   = { plt_lbl };

cfg.y_lim     = [-10e5 70e5]; % [-1e7 20e7];

cfg.lgd       = 0;
cfg.std_err   = 1;

cfg.lnstyle.col_ord = repmat({{rgb('red')       rgb('blue')}},8,8);

cfg.stt_dat = repmat({{'vis_new_ovr_stt' 'aud_new_ovr_stt'}},8,8);

cfg.stt_col = repmat({{rgb('reddish grey') rgb('bluish grey') }},8,8);

cfg.stt_cmp = repmat({{'0%4' '4%8' }},8,8);

cfg.cnd_nme         = {};

cfg.alt_lbl         = 'label';

cfg.x_lim       = [-0.3 1.2];
cfg.v_lne       = 0;
cfg.v_lne_col   = repmat({rgb('black')},8,8);
cfg.v_lne_wdt   = 2;
cfg.axe_fnt_sze = [ 5 ];
cfg.axe_lne_sze = [ 0.5 ];
cfg.ttl_lne_sze = [ 1 ];

cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'png';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ' '/' 'clerical' '/' 'manuscript' '/' 'figure9' '/' 'chn_plt'];
cfg.prefix     = ['MG049_SASZ' '_' 'STG' '_' 'focus'];

mmil_ieeg_sensor_plot_v5(cfg)

%% AVERAGE HGP UTAH PLOT
% Load Data
cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' fcfg.sbj_nme '_overall_data.mat'];
uth_dat  = ft_func([],cfg);

uth_dat.(uth_dat.data_name{1}).cfg.alt_eve.plt_eve(801:1200) = uth_dat.(uth_dat.data_name{1}).cfg.alt_eve.vis_new_bse(801:1200);
uth_dat.(uth_dat.data_name{1}).cfg.alt_eve.plt_eve(1:801)    = uth_dat.(uth_dat.data_name{1}).cfg.alt_eve.aud_new_bse(1:801);
uth_dat.(uth_dat.data_name{1}).cfg.alt_eve.plt_eve(1:201)    = 0;

% Smooth Data
% cfg             = [];
% cfg.data_name   = 1;
% HGP_smooth_msec = 0.125;
% w               = round(HGP_smooth_msec* uth_dat.(uth_dat.data_name{1}).fsample);
% gauss_x         = -w:1:w;
% gauss_y         = normpdf(gauss_x,0,round(0.024* uth_dat.(uth_dat.data_name{1}).fsample));
% cfg.window      = gauss_y/sum(gauss_y);
% uth_dat         = ft_func(@ft_window_data,cfg,uth_dat);

% Stat Data


% Make Average Utah Plot
% Create Data
chn_num = eff_typ{3};

for iT = 1:numel(uth_dat.uth_cmb.trial)
    uth_dat_hld(iT,:) = nanmean(uth_dat.uth_cmb.trial{iT}(chn_num,:),1);
end

vis_avg = nanmean(uth_dat_hld(uth_dat.uth_cmb.cfg.alt_eve.vis_new_bse==102,:),1);
aud_avg = nanmean(uth_dat_hld(uth_dat.uth_cmb.cfg.alt_eve.aud_new_bse==202,:),1);

figure()
plot(uth_dat.uth_cmb.time{1},vis_avg,'r'); hold on;
plot(uth_dat.uth_cmb.time{1},aud_avg,'b');
vline(0,'k')
xlim([-0.150 1.2])

% Repetition
chn_num = unique([eff_typ{3}' eff_typ{1}']);
for iT = 1:numel(uth_dat.uth_cmb.trial)
    uth_dat_hld(iT,:) = nanmean(uth_dat.uth_cmb.trial{iT}(chn_num,:),1);
end
vis_old_avg = nanmean(uth_dat_hld(uth_dat.uth_cmb.cfg.alt_eve.vis_new_old==112,:),1);
vis_new_avg = nanmean(uth_dat_hld(uth_dat.uth_cmb.cfg.alt_eve.vis_new_old==111,:),1);

chn_num = unique([eff_typ{3}' eff_typ{2}']);
for iT = 1:numel(uth_dat.uth_cmb.trial)
    uth_dat_hld(iT,:) = nanmean(uth_dat.uth_cmb.trial{iT}(chn_num,:),1);
end
aud_old_avg = nanmean(uth_dat_hld(uth_dat.uth_cmb.cfg.alt_eve.aud_new_old==212,:),1);
aud_new_avg = nanmean(uth_dat_hld(uth_dat.uth_cmb.cfg.alt_eve.aud_new_old==211,:),1);

figure()
subplot(2,2,1)
plot(uth_dat.uth_cmb.time{1},vis_new_avg,'r'); hold on;
plot(uth_dat.uth_cmb.time{1},vis_old_avg,'Color',rgb('dark orange'));
vline(0,'k')
xlim([-0.150 1.2])
subplot(2,2,2)
plot(uth_dat.uth_cmb.time{1},aud_new_avg,'b'); hold on;
plot(uth_dat.uth_cmb.time{1},aud_old_avg,'Color',rgb('dark purple'));
vline(0,'k')
xlim([-0.150 1.2])

% Put together
for iT = 1:numel(uth_dat.(uth_dat.data_name{1}).trial)
    uth_dat.(uth_dat.data_name{1}).trial{iT}(97,:) = uth_dat_hld(iT,:);
end

uth_dat.(uth_dat.data_name{1}).label{97} = 'avg';
uth_dat.(uth_dat.data_name{1}).cfg.alt_lab.label{97} = 'avg';

cfg = [];
cfg.channel = ['all',strcat('-',uth_dat.(uth_dat.data_name{1}).label(1:96))'];
uth_dat     = ft_func(@ft_preprocessing,cfg,uth_dat);

uth_dat.(uth_dat.data_name{1}).label = [];
uth_dat.(uth_dat.data_name{1}).cfg.alt_lab.label = [];

uth_dat.(uth_dat.data_name{1}).label{1} = 'avg';
uth_dat.(uth_dat.data_name{1}).cfg.alt_lab.label{1} = 'avg';

% Stats
cfg = [];
cfg.stt_fnc  = {'sasz_vis_ovr_new' 'sasz_aud_ovr_new'};
cfg.loc      = 'local';
cfg.fld_nme  = uth_dat.data_name;
cfg.specific = {'fld_nme' ; 1:numel(cfg.fld_nme)};
uth_dat      = ft_func(@mmil_cloud_stat,cfg,uth_dat);

uth_dat.(uth_dat.data_name{1}).cfg.alt_stt.aud_new_ovr_stt_avg.prob(1:250) = 1;

% Make Plot
cfg = [];

cfg.dat = {uth_dat.(uth_dat.data_name{1})};

cfg.type      = 'macro';

cfg.chn_grp = { 1 };
cfg.dat_loc = 1;
cfg.plt_dim = [1 1];

cfg.eve       = {{[102 202]}};

cfg.alt_eve   = { {'plt_eve'} };

cfg.plt_lbl   = { '' };

cfg.y_lim     = 'maxmin'; %[-10e5 70e5]; % [-1e7 20e7];

cfg.lgd       = 0;
cfg.std_err   = 1;

cfg.lnstyle.col_ord = {{rgb('red')       rgb('blue')}};

cfg.stt_dat = {{'vis_new_ovr_stt_avg' 'aud_new_ovr_stt_avg'}};

cfg.stt_col = {{rgb('reddish grey') rgb('bluish grey') }};

cfg.stt_cmp = {{'0%4' '4%8' }};

cfg.cnd_nme         = {};

cfg.alt_lbl         = 'label';

cfg.x_lim       = [-0.150 1.0];
cfg.v_lne       = 0;
cfg.v_lne_col   = {rgb('black')};
cfg.v_lne_wdt   = 2;
cfg.axe_fnt_sze = [ 5 ];
cfg.axe_lne_sze = [ 0.5 ];
cfg.ttl_lne_sze = [ 1 ];

cfg.print      = 1;
cfg.nofig      = 1;
cfg.print_type = 'png';
cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ' '/' 'clerical' '/' 'manuscript' '/' 'figure9' '/' 'avg_chn_plt'];
cfg.prefix     = ['MG049_SASZ' '_' 'STG' '_' 'average'];

mmil_ieeg_sensor_plot_v5(cfg)

cfg.print_type = 'eps';
mmil_ieeg_sensor_plot_v5(cfg)

%% Find signififance
cfg = [];

% Setting up Variables
cfg.sbj_nme     = fcfg.sbj_nme;
cfg.dat_fld     = fcfg.dat_fld;
cfg.clr_fld     = fcfg.clr_fld;

iSASZ_utah_eff(cfg);



%% Create Dot Plot
% Visual & Auditory Overall
vis_ovr_eff = mmil_readtext([ fcfg.clr_fld '/' 'sig_chn' '/' 'MG49_SA_SZ_utah' '/' 'cmb' '/' 'effect' '/' 'vis_new_ovr_stt' ]);
aud_ovr_eff = mmil_readtext([ fcfg.clr_fld '/' 'sig_chn' '/' 'MG49_SA_SZ_utah' '/' 'cmb' '/' 'effect' '/' 'aud_new_ovr_stt' ]);

eff_typ{1} = find(cell2mat(vis_ovr_eff(2:end,3)));
eff_typ{2} = find(cell2mat(aud_ovr_eff(2:end,3)));
eff_typ{3} = find(cell2mat(vis_ovr_eff(2:end,3)) & cell2mat(aud_ovr_eff(2:end,3)));

eff_typ{1} = setxor(eff_typ{1},eff_typ{3});
eff_typ{2} = setxor(eff_typ{2},eff_typ{3});

eff_col{1} = 'red';
eff_col{2} = 'blue';
eff_col{3} = 'bright purple';

rmv_chn = [2 65 33 67 35 22 48 71 44 75 76 51 57 89 94 95];

figure('units','normalized','outerposition',[0 0 0.31 1])
for iR = 1:10
    for iC = 1:10
        
        pin_num = zeros(10,10);
        pin_num(iR,iC) = 1;
        pin_num = mmil_find_utaharray_pinnumber(pin_num);
        
        if (iR==1 && iC==1) || (iR==1 && iC==10) || (iR==10 && iC==1) || (iR==10 && iC==10) || any(pin_num==rmv_chn)
            
        elseif any(eff_typ{1} == pin_num)
            rectangle('Position',[iC-0.25 -iR-0.25 0.5 0.5],'Curvature',[1 1],'FaceColor',rgb(eff_col{1}));
        elseif any(eff_typ{2} == pin_num)
            rectangle('Position',[iC-0.25 -iR-0.25 0.5 0.5],'Curvature',[1 1],'FaceColor',rgb(eff_col{2}));
        elseif any(eff_typ{3} == pin_num)
            rectangle('Position',[iC-0.25 -iR-0.25 0.5 0.5],'Curvature',[1 1],'FaceColor',rgb(eff_col{3}));
        else
            rectangle('Position',[iC-0.25 -iR-0.25 0.5 0.5],'Curvature',[1 1],'FaceColor',rgb('black'));
        end
        
    end
end
axis('off')
tightfig()

cfg = [];
cfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'chn_plt' '/' 'sig_dis'];
cfg.prn_typ = 'png';
mmil_print_plot(cfg)

%% Create Trial Figures
pin_nme{1}  = 'OVR_TOP_RGH';
pin_nme{2}  = 'OVR_BOT_RGH';
pin_nme{3}  = 'OVR_MID_RGH';

pin_num{1}  = mmil_find_utaharray_pinnumber([ 0 0 0 0 1 1 1 1 1 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ]); pin_num{1} = repmat(pin_num{1}',5,1);
                                          
pin_num{2}  = mmil_find_utaharray_pinnumber([ 0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 1 1 1 1 1 0 ]); pin_num{2} = repmat(pin_num{2}',5,1);

pin_num{3}  = mmil_find_utaharray_pinnumber([ 0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 1 ; ...
                                              0 0 0 0 0 0 0 0 0 1 ; ...
                                              0 0 0 0 0 0 0 0 0 1 ; ...
                                              0 0 0 0 0 0 0 0 0 1 ; ...
                                              0 0 0 0 0 0 0 0 0 1 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ; ...
                                              0 0 0 0 0 0 0 0 0 0 ]); pin_num{3} = repmat(pin_num{3}',5,1);
ylm_lim{1} = [-3.5e8  15e8];
ylm_lim{2} = [-3.5e8  10e8];
ylm_lim{3} = [-3.5e8  10e8];                               

    

% SA/SZ - overall
for iP = 1:numel(pin_num)
           
    cfg = [];
    
    cfg.dat = {uth_dat.(uth_dat.data_name{1})};
    
    cfg.type      = 'macro';
    
    cfg.chn_grp = pin_num(iP);
    cfg.dat_loc = pin_num{iP}>0;
    cfg.plt_dim = size(pin_num{iP});
    
    cfg.eve       = {[repmat({[101 201]},1,size(pin_num{iP},2)) ; ...
                      repmat({[211 212]},1,size(pin_num{iP},2)) ; ...
                      repmat({[111 112]},1,size(pin_num{iP},2)) ; ...
                      repmat({[251 252]},1,size(pin_num{iP},2)) ; ...
                      repmat({[151 152]},1,size(pin_num{iP},2)) ]};
    
    cfg.alt_eve   = {[repmat({'ovr'},1,size(pin_num{iP},2)) ; ...
                      repmat({'aud_new_old'},1,size(pin_num{iP},2)) ; ...
                      repmat({'vis_new_old'},1,size(pin_num{iP},2)) ; ...
                      repmat({'aud_frq_med'},1,size(pin_num{iP},2)) ; ...
                      repmat({'vis_frq_med'},1,size(pin_num{iP},2)) ]};
    
    cfg.plt_lbl   = {repmat({''},5,size(pin_num{iP},2))};
    
    cfg.y_lim     = ylm_lim{iP};
    
    cfg.lgd       = 0;
    cfg.std_err   = 1;
    
    cfg.lnstyle.col_ord = {[repmat({{rgb('red')       rgb('blue')}},1,size(pin_num{iP},2)) ; ...
                            repmat({{rgb('blue')      rgb('dark orange')}},1,size(pin_num{iP},2)) ; ...
                            repmat({{rgb('red')       rgb('dark orange')}},1,size(pin_num{iP},2)) ; ...
                            repmat({{rgb('dark blue') rgb('bright blue')}},1,size(pin_num{iP},2)) ; ...
                            repmat({{rgb('dark red')  rgb('bright red')}},1,size(pin_num{iP},2)) ]};
    
    cfg.stt_dat = {[repmat({{'vis_new_ovr_stt' 'aud_new_ovr_stt'}},1,size(pin_num{iP},2)) ; ...
                    repmat({{'aud_new_old_stt'}},1,size(pin_num{iP},2)) ; ...
                    repmat({{'vis_new_old_stt'}},1,size(pin_num{iP},2)) ; ...
                    repmat({{'aud_frq_stt'}},1,size(pin_num{iP},2)) ; ...
                    repmat({{'vis_frq_stt'}},1,size(pin_num{iP},2)) ]};
    
    cfg.stt_col =  {[repmat({{rgb('reddish grey') rgb('bluish grey') }},1,size(pin_num{iP},2)) ; ...
                    repmat({{rgb('blue')      }},1,size(pin_num{iP},2)) ; ...
                    repmat({{rgb('red')       }},1,size(pin_num{iP},2)) ; ...
                    repmat({{rgb('dark blue') }},1,size(pin_num{iP},2)) ; ...
                    repmat({{rgb('dark red')  }},1,size(pin_num{iP},2)) ]};
    
    cfg.stt_cmp = {[repmat({{'0%4' '4%8' }},1,size(pin_num{iP},2)) ; ...
                    repmat({{'0%4'      }},1,size(pin_num{iP},2)) ; ...
                    repmat({{'0%4'      }},1,size(pin_num{iP},2)) ; ...
                    repmat({{'0%4' }},1,size(pin_num{iP},2)) ; ...
                    repmat({{'0%4'  }},1,size(pin_num{iP},2)) ]};
    
    cfg.cnd_nme         = {};
    
    cfg.alt_lbl         = 'label';
    
    cfg.x_lim       = [-0.150 1.0];
    cfg.v_lne       = 0;
    cfg.v_lne_col   = repmat({rgb('black')},5,size(pin_num{iP},2));
    cfg.v_lne_wdt   = 2;
    cfg.axe_fnt_sze = [ 5 ];
    cfg.axe_lne_sze = [ 0.5 ];
    cfg.ttl_lne_sze = [ 1 ];
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ' '/' 'clerical' '/' 'manuscript' '/' 'figure9' '/' 'utah_plot'];
    cfg.prefix     = ['MG049_SASZ' '_' pin_nme{iP}];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
end

% Chan Num
pin_num_chn = [pin_num{1}(1,2:5)  pin_num{2}(1,2:5) pin_num{3}(1,1:4)];
pin_nme = [repmat({'OVR_TOP_RGH'},1,4) repmat({'OVR_BOT_RGH'},1,4) repmat({'OVR_MID_RGH'},1,4)];
pin_ylm = [repmat({ylm_lim{1}},1,4) repmat({ylm_lim{2}},1,4) repmat({ylm_lim{3}},1,4)];
 
for iP = 1:numel(pin_num_chn)
           
    cfg = [];
    
    cfg.dat = {uth_dat.(uth_dat.data_name{1})};
    
    cfg.type      = 'macro';
    
    cfg.chn_grp = {pin_num_chn(iP)};
    cfg.dat_loc = 1;
    cfg.plt_dim = [1 1];
    
    cfg.eve       = { repmat({[101 201]},1,1) };
    
    cfg.alt_eve   = { repmat({'ovr'},1,1) };
    
    cfg.plt_lbl   = {repmat({''},5,1)};
    
    cfg.y_lim     = pin_ylm{iP}; %'maxmin';%[-0.5e8 2e8];
    
    cfg.lgd       = 0;
    cfg.std_err   = 1;
    
    cfg.lnstyle.col_ord = { repmat({{rgb('red')       rgb('blue')}},1,1) };
    
    cfg.stt_dat = { repmat({{'vis_new_ovr_stt' 'aud_new_ovr_stt'}},1,1) };
    
    cfg.stt_col =  { repmat({{rgb('reddish grey') rgb('bluish grey') }},1,1) };
    
    cfg.stt_cmp = { repmat({{'0%4' '4%8' }},1,1) };
    
    cfg.cnd_nme         = {};
    
    cfg.alt_lbl         = 'label';
    
    cfg.x_lim       = [-0.150 1.0];
    cfg.v_lne       = 0;
    cfg.v_lne_col   = repmat({rgb('black')},1,1);
    cfg.v_lne_wdt   = 2;
    cfg.axe_fnt_sze = [ 5 ];
    cfg.axe_lne_sze = [ 0.5 ];
    cfg.ttl_lne_sze = [ 1 ];
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ' '/' 'clerical' '/' 'manuscript' '/' 'figure9' '/' 'utah_plot'];
    cfg.prefix     = ['MG049_SASZ' '_' pin_nme{iP} '_' num2str(iP)];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
    cfg.print_type = 'eps';
    mmil_ieeg_sensor_plot_v5(cfg)
    
end

%% Amplitude
dta_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';

% Setup
fcfg = [];

fcfg.dta_fld = dta_fld;
fcfg.clr_fld = clr_fld;

fcfg.typ     = {'hgp'};
fcfg.ele_typ = {'ecog'};
fcfg.loc_typ = {'split'};

fcfg.alt_eve = { 'aud_new_bse' 'vis_new_bse'}; % 'aud_new_bse' % 'vis_new_bse' 'aud_new_bse'
fcfg.eve     = [ 202           102 ];     % 202           % 102           202
fcfg.eve_nme = { 'aud_act'     'vis_act' };     % 'aud_act'     % 

fcfg.dta_typ = 1;

fcfg.zsc_bse = [-0.150 0.000];
fcfg.zsc_epc = [0.050  1.000] ;

sbj_nme = 'MG49_SA_SZ_utah';

%
cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dta_fld '/' sbj_nme '_overall_data.mat'];
bcc_dat     = ft_func([],cfg);

mdl_ele(:,2) = bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_lab.label;

dta_hld = cat(3,bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).trial{:});

tme_hld = bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).time{1};

tme_bse_beg = dsearchn(tme_hld',fcfg.zsc_bse(1)); tme_bse_end = dsearchn(tme_hld',fcfg.zsc_bse(2));
tme_epc_beg = dsearchn(tme_hld',fcfg.zsc_epc(1)); tme_epc_end = dsearchn(tme_hld',fcfg.zsc_epc(2));

for iA = 1:numel(fcfg.alt_eve)
    
    if isfield(bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_eve,fcfg.alt_eve{iA})
        
        eve_hld = bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_eve.(fcfg.alt_eve{iA});
        eve_loc = find(eve_hld==fcfg.eve(iA));
        
        for iC = 1:numel(bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_lab.label)
            
            ele_num = find(strcmpi(mdl_ele(:,2),bcc_dat.(bcc_dat.data_name{fcfg.dta_typ}).cfg.alt_lab.label{iC}));
            
            if ~isempty(ele_num)
                
                dta_men = nanmean(dta_hld(iC,:,eve_loc),3);
                zsc_hld = ( dta_men - nanmean(dta_men(tme_bse_beg:tme_bse_end)) ) ./ std(dta_men(tme_bse_beg:tme_bse_end));
                
                mdl_ele{ele_num,iA+2} = max(zsc_hld(tme_epc_beg:tme_epc_end)); % max(nanmean(dta_zsc_hld(iC,tme_epc_beg:tme_epc_end,eve_loc),3));
                
            end
        end
    end
end

% Plot
chn_inc = unique([eff_typ{1}' eff_typ{2}' eff_typ{3}']);

aud_amp = cell2mat(mdl_ele(chn_inc,3));
vis_amp = cell2mat(mdl_ele(chn_inc,4));


figure('units','normalized','outerposition',[0 0 0.31 1])

scatter(vis_amp,aud_amp,102, ...
        'MarkerEdgeColor',rgb('black'), ...
        'MarkerFaceColor',rgb('light grey'), ...
        'LineWidth',1)
xlim([0 40])
ylim([0 40])

set(gca,'FontSize',22);
set(gca,'linewidth',2.5);
set(gca,'TickLength',[0 0]);
set(gcf,'color','w')

tightfig()

cfg = [];
cfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'amp_crr'];
cfg.prn_typ = 'eps';
mmil_print_plot(cfg)

close all

% Correlation
[cor_val,cor_pvl] = corrcoef([aud_amp vis_amp]);

end
