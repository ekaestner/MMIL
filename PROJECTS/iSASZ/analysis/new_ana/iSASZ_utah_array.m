function iSASZ_utah_array(fcfg)

%% Load Data
cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme{1} '_overall_data.mat'];
isz_dat     = ft_func([],cfg);

cfg = [];
cfg.load    = 'yes';
cfg.file    = [fcfg.dat_fld '/' fcfg.sbj_nme{2} '_overall_data.mat'];
isa_dat     = ft_func([],cfg);

uth_dat.data_name{1} = isz_dat.data_name{1};
uth_dat.data_name{2} = isa_dat.data_name{1};
uth_dat.(uth_dat.data_name{1}) = isz_dat.(isz_dat.data_name{1});
uth_dat.(uth_dat.data_name{2}) = isa_dat.(isa_dat.data_name{1});

% epoch data & Remove trigger chan
trl_hld{1} = ft_Utah_trialfunMG49SZ(0); trl_hld{1}(:,1) = trl_hld{1}(:,1)-500; trl_hld{1}(:,2) = trl_hld{1}(:,2)+500; trl_hld{1}(:,3) = trl_hld{1}(:,3)-500;
trl_hld{2} = ft_Utah_trialfunMG49SA(0); trl_hld{2}(:,1) = trl_hld{2}(:,1)-500; trl_hld{2}(:,2) = trl_hld{2}(:,2)+500; trl_hld{2}(:,3) = trl_hld{2}(:,3)-500;

cfg = [];
cfg.specific  = {'trl';1:numel(trl_hld)};
cfg.trl       = trl_hld;
uth_dat = ft_func(@ft_redefinetrial,cfg,uth_dat);

uth_dat.(uth_dat.data_name{1}).label = uth_dat.(uth_dat.data_name{1}).label';
uth_dat.(uth_dat.data_name{2}).label = uth_dat.(uth_dat.data_name{2}).label';

cfg = [];
cfg.sbj_nme = fcfg.sbj_nme;
cfg.rmv_chn = {[97 98]' [97 98]'}; %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/rmv_chn/' subj]),'uni',0)';
uth_dat     = mmil_remove_channels(cfg,uth_dat);

% Examine Data for Noise - Initial
% cfg = [];
% cfg.empty  = 'yes';
% cfg.outdir = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'spectrum' '/' 'MG49_Utah' '/' 'Initial' ];
% cfg.prefix = uth_dat.data_name;
% cfg.specific  = {'prefix'; 1:numel(uth_dat.data_name)};
% ft_func(@ft_plot_spectrum,cfg,uth_dat);

chn_nse_grp     = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'sbj_inf' '/' 'MG49_Utah'],'cmn_nse');%mmil_readtext([fcfg.clr_fld '/cmn_nse/' subj]);
chn_nse_grp_nme = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'sbj_inf' '/' 'MG49_Utah'],'cmn_nme');
nse_val         = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'sbj_inf' '/' 'MG49_Utah'],'nse_val'); nse_val = nse_val{1};

% Make temporary ECOG/DEPTH Split
for iD = 1:numel(uth_dat.data_name); uth_dat.(uth_dat.data_name{iD}).cfg.alt_lab.label = uth_dat.(uth_dat.data_name{iD}).label; end
scfg = [];
scfg.dat_nme = uth_dat.data_name;
scfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
scfg.sbj_nme = 'MG49_Utah';
scfg.alt_lab = 'label';
scfg.specific = {'dat_nme' ; 1:numel(uth_dat.data_name)};
scfg.add_nme  = '_tmp';
ft_func(@mmil_create_depth2,scfg,uth_dat);

if nse_val
    
    if strcmpi(chn_nse_grp{1},'split'); chn_nse_grp = repmat(chn_nse_grp,1,numel(uth_dat.data_name)); end
    
    for iEJ = 1:numel(chn_nse_grp_nme);
        if strcmpi(chn_nse_grp{1},'split')
            for iD = 1:numel(uth_dat.data_name); pre_fix{iD}{1} = strcat(uth_dat.data_name{iD},'ecog'); pre_fix{iD}{2} = strcat(uth_dat.data_name{iD},'noisy_ecog'); pre_fix{iD}{3} = strcat(uth_dat.data_name{iD},'depth'); pre_fix{iD}{4} = strcat(uth_dat.data_name{iD},'noisy_depth'); end;
        else
            for iCN = 1:numel(chn_nse_grp_nme{iEJ});
                pre_fix{iEJ}{iCN} = strcat(uth_dat.data_name{iEJ},'_',chn_nse_grp_nme{iEJ}{iCN});
            end;
        end
    end
    
    cfg             = [];
    cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
    cfg.sbj_nme = 'MG49_Utah';
    cfg.pre_fix     = pre_fix;
    cfg.dat_nme     = uth_dat.data_name;
    cfg.chn_nse_grp = chn_nse_grp;
    cfg.pre_fix     = pre_fix;
    cfg.specific    = {'pre_fix' 'chn_nse_grp' 'dat_nme' ; 1:numel(uth_dat.data_name) 1:numel(uth_dat.data_name) 1:numel(uth_dat.data_name)};
    cfg.out_dir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'cmn_nse_rmv'];
    cfg.rmv_chn     = nse_val;
    uth_dat = ft_func(@ft_remove_common_noise2,cfg,uth_dat);
    
    if nse_val == 1
        cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';
        cfg.sbj_nme = 'MG49_Utah';
        cfg.sub_fld_ind = num2cell(1:numel(uth_dat.data_name));
        cfg.specific    = {'pre_fix'  'sub_fld_ind' ; 1:numel(uth_dat.data_name) 1:numel(uth_dat.data_name)};
        uth_dat = ft_func(@mmil_update_cmn_nse2,cfg,uth_dat);
    end
    
end

% Remove identified noise problems
bse_frq     = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'sbj_inf' '/' 'MG49_Utah'],'bse_frq'); %cellfun(@str2num,mmil_readtext([fcfg.clr_fld '/bse_frq/' subj]),'uni',0)';

if ~isempty(bse_frq)
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = bse_frq;
    cfg.multi    = {'bsfreq'; 1:numel(cfg.bsfreq)};
    uth_dat      = ft_func(@ft_preprocessing,cfg,uth_dat);
end

% Filter Data for HGP - EJK
cfg=[];
cfg.data_name  = 1:2;
cfg.data_new   = 'yes';
cfg.new_suffix = 'hgp';
cfg.foi    = [70 80 90 100 110 130 140 150 160 170];    %frequency of interest
cfg.sf     = [repmat(10,1,numel(cfg.foi))]; %specific frequency
cfg.width  = cfg.foi./10;
cfg.gwidth = ones(size(cfg.foi))*pi; %wavelet
cfg.keeptrials = 'yes';
cfg.method = 'tfr';
cfg.toi=uth_dat.(uth_dat.data_name{1}).time{1};
cfg.keeptrials = 'yes';
uth_dat        = ft_func(@mmil_hgp_freq_analysis,cfg,uth_dat);

cfg             = [];
cfg.data_name   = 1:2;
HGP_smooth_msec = 0.125;
w               = round(HGP_smooth_msec* uth_dat.(uth_dat.data_name{1}).fsample);
gauss_x         = -w:1:w;
gauss_y         = normpdf(gauss_x,0,round(0.016* uth_dat.(uth_dat.data_name{1}).fsample));
cfg.window      = gauss_y/sum(gauss_y);
uth_dat         = ft_func(@ft_window_data,cfg,uth_dat);

% Remove superfluous continuous data structures
cfg = [];
cfg.data_name = 1:2;
cfg.rmfield   = 'yes';
uth_dat       = ft_func([],cfg,uth_dat);

% Baseline Data
bse_tme = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'sbj_inf' '/' 'MG49_Utah'],'bse_tme'); bse_tme = bse_tme{1};

cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = bse_tme;
uth_dat            = ft_func(@ft_baseline,cfg,uth_dat);

% Automatic rejection & Apply Rejections
epc_tme = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'sbj_inf' '/' 'MG49_Utah'],'epc_tme'); epc_tme = epc_tme{1};
rjt_plt = mmil_load_subj_info(['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical' '/' 'sbj_inf' '/' 'MG49_Utah'],'rjt_plt'); rjt_plt = rjt_plt{1};

cfg = [];
cfg.measures  = {'time' 'time-1' 'variance'};
cfg.thresh    = [0.9 0.9 0.9];
cfg.outdir    = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'rejection' '/' 'MG49_Utah' '/' ];
cfg.prefix    = uth_dat.data_name;
cfg.specific  = {'prefix';1:numel(uth_dat.data_name)};
cfg.pad       = uth_dat.(uth_dat.data_name{1}).time{1}(end)-epc_tme(2);
cfg.plot      = 0;
uth_dat       = ft_func(@auto_rej,cfg,uth_dat);

cfg = [];
cfg.measure = 'all';
cfg.apply   = 'ieeg';
uth_dat     = ft_func(@ft_apply_rejection,cfg,uth_dat);

%% Combine data
frs_num = 2;
scd_num = 1;

cfg           = [];
cfg.data_name = [frs_num ; scd_num];
cfg.data_new  = 'yes';
cfg.methapp   = 'trials';
dat = ft_func(@ft_appenddata,cfg,uth_dat);

uth_dat.data_name{3} = 'uth_cmb';
uth_dat.(uth_dat.data_name{3}) = dat.(dat.data_name{2});

clear dat

%% Remove channels
rmv_chn = [2 22 33 35 44 48 51 57 65 67 75 76 89 94 95];

for iD = 1:numel(uth_dat.data_name)
    %     for iRM = 1:numel(rmv_chn)
    for iT = 1:numel(uth_dat.(uth_dat.data_name{3}).trial)
        uth_dat.(uth_dat.data_name{iD}).trial{iT}(rmv_chn,:) = nan(size(uth_dat.(uth_dat.data_name{iD}).trial{iD}(rmv_chn,:)));
    end
    %     end
end

%% Events
cfg = [];
cfg.load    = 'yes';
cfg.file    = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'MG049_SA_SZ' '_overall_data.mat'];
eve_dat     = ft_func([],cfg);

eve_fld = fieldnames(eve_dat.(eve_dat.data_name{2}).cfg.alt_eve);
eve_dat.(eve_dat.data_name{2}).cfg.alt_eve.(eve_fld{numel(eve_fld)}) = eve_dat.(eve_dat.data_name{2}).cfg.alt_eve.(eve_fld{numel(eve_fld)})';
for iE = 1:numel(eve_fld)
    uth_dat.(uth_dat.data_name{1}).cfg.alt_eve.(eve_fld{iE}) = eve_dat.(eve_dat.data_name{2}).cfg.alt_eve.(eve_fld{iE})(801:end,:);
    uth_dat.(uth_dat.data_name{2}).cfg.alt_eve.(eve_fld{iE}) = eve_dat.(eve_dat.data_name{2}).cfg.alt_eve.(eve_fld{iE})(1:800,:);
    uth_dat.(uth_dat.data_name{3}).cfg.alt_eve.(eve_fld{iE}) = eve_dat.(eve_dat.data_name{2}).cfg.alt_eve.(eve_fld{iE});
end

uth_dat.(uth_dat.data_name{3}).cfg.alt_eve.ovr = [uth_dat.(uth_dat.data_name{3}).cfg.alt_eve.aud_ovr(1:800,:) ; uth_dat.(uth_dat.data_name{3}).cfg.alt_eve.vis_ovr(801:end,:)];

uth_dat.(uth_dat.data_name{1}).cfg.alt_lab.label = uth_dat.(uth_dat.data_name{1}).label;
uth_dat.(uth_dat.data_name{2}).cfg.alt_lab.label = uth_dat.(uth_dat.data_name{2}).label;
uth_dat.(uth_dat.data_name{3}).cfg.alt_lab.label = uth_dat.(uth_dat.data_name{3}).label;

%% Latency
cfg = [];
cfg.data_name = 1:2;
cfg.rmfield   = 'yes';
uth_dat       = ft_func([],cfg,uth_dat);

cfg = [];
cfg.latency = [-0.5 1.5];
uth_dat     = ft_func(@ft_selectdata,cfg,uth_dat);

%% Save
cfg = [];
cfg.str_nme  = 'uth_dat';
cfg.save     = 'yes';
cfg.filename = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'MG49_SA_SZ' '_utah_overall_data_v2.mat'];
ft_func([],cfg,uth_dat);

%% Channel Plot
% cfg = [];
% 
% cfg.type      = 'chan';
% 
% cfg.dat = {uth_dat.(uth_dat.data_name{1}) uth_dat.(uth_dat.data_name{2})};
% 
% cfg.plt_dim   = [4 4];
% cfg.dat_loc   = [1 1 2 2 ; 1 0 2 0 ; 1 1 2 2 ; 0 0 0 0];
% 
% cfg.alt_eve   = {'vis_ovr' 'vis_new_old' 'aud_ovr' 'aud_new_old' ; 'vis_ani_obj' '' 'aud_ani_obj' '' ; 'vis_big_med' 'vis_frq_med' 'aud_big_med' 'aud_frq_med' ; '' '' '' ''};
% cfg.eve       = {[101]     [111 112]     [201]     [211 212]     ; [121 122]     [] [221 222]     [] ; [141 142]     [151 152]     [241 242]     [251 252]     ; [] [] [] []};
% 
% cfg.plt_lbl   = { 'Visual Overall' 'Visual New/Old' 'Auditory Overall' 'Auditory New/Old' ; 'Visual Animal/Object' '' 'Auditory Animal/Object' '' ; 'Visual Bigram' 'Visual Frequency' 'Auditory BiPhoneme' 'Auditory Frequency' ; '' '' '' '' };
% 
% cfg.ttl_num   = zeros(4,4);
% 
% cfg.y_lim     = 'maxmin';
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = { {rgb('reddish grey')} {rgb('red') rgb('dark orange')} {rgb('bluish grey')} {rgb('blue') rgb('dark orange')} ; {rgb('red') rgb('grey')} {} {rgb('blue') rgb('grey')} {} ; {rgb('dark red') rgb('bright red')} {rgb('dark red') rgb('bright red')} {rgb('dark blue') rgb('bright blue')} {rgb('dark blue') rgb('bright blue')} ; {} {} {} {} };
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = { 0 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 0 };
% cfg.v_lne_col   = { {rgb('red')} {rgb('red')} {rgb('blue')} {rgb('blue')} ; {rgb('red')} {} {rgb('blue')} {} ; {rgb('red')} {rgb('red')} {rgb('blue')} {rgb('blue')} ; {} {} {} {} };
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [15 15 15 15 ; ...
%                    15 15 15 15 ; ...
%                    15 15 15 15 ; ...
%                    15 15 15 15 ];
% cfg.axe_lne_sze = [3 3 3 3 ; ...
%                    3 3 3 3 ; ...
%                    3 3 3 3 ; ...
%                    3 3 3 3 ];
% cfg.ttl_lne_sze = [36 36 36 36 ; ...
%                    36 36 36 36 ; ...
%                    36 36 36 36 ; ...
%                    36 36 36 36 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/'];
% cfg.prefix     = ['MG049_SA_SZ' '_'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)

%% Utah Plots
% SA/SZ - overall
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {uth_dat.(uth_dat.data_name{3})};
% 
% cfg.alt_eve   = 'ovr';
% cfg.eve       = [101 201];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('red') rgb('blue')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('black')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot'];
% cfg.prefix     = ['MG049_SASZ' '_'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)
% 
% % SA - New/Old
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {uth_dat.(uth_dat.data_name{2})};
% 
% cfg.alt_eve   = 'aud_new_old';
% cfg.eve       = [211 212];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('blue') rgb('dark orange')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('blue')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot'];
% cfg.prefix     = ['MG049_SA' '_'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)
% 
% % SA - Frequency
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {uth_dat.(uth_dat.data_name{2})};
% 
% cfg.alt_eve   = 'aud_frq_med';
% cfg.eve       = [251 252];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('dark blue') rgb('bright blue')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('blue')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot'];
% cfg.prefix     = ['MG049_SA' '_frequency'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)
% 
% % SZ - New/Old
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {uth_dat.(uth_dat.data_name{1})};
% 
% cfg.alt_eve   = 'vis_new_old';
% cfg.eve       = [111 112];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('red') rgb('dark orange')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('red')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot'];
% cfg.prefix     = ['MG049_SZ' '_'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)
% 
% % SZ - Frequency
% cfg = [];
% 
% cfg.type      = 'utah';
% 
% cfg.dat = {uth_dat.(uth_dat.data_name{1})};
% 
% cfg.alt_eve   = 'vis_frq_med';
% cfg.eve       = [151 152];
% 
% cfg.plt_lbl   = {repmat({''},10,10)};
% 
% cfg.y_lim     = [-0.5e9 2e9];
% 
% cfg.lgd       = 0;
% cfg.std_err   = 1;
% cfg.lnstyle.col_ord = {rgb('dark red') rgb('bright red')};
% cfg.cnd_nme         = {};
% 
% cfg.alt_lbl         = 'label';
% 
% cfg.x_lim       = [-0.3 1.5];
% cfg.v_lne       = 0;
% cfg.v_lne_col   = repmat({rgb('red')},10,10);
% cfg.v_lne_wdt   = 2;
% cfg.axe_fnt_sze = [ 5 ];
% cfg.axe_lne_sze = [ 0.5 ];
% cfg.ttl_lne_sze = [ 1 ];
% 
% cfg.print      = 1;
% cfg.nofig      = 1;
% cfg.print_type = 'png';
% cfg.outdir     = ['/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data' '/' 'utah_plot' '/' 'MG049' '/' 'utah_plot'];
% cfg.prefix     = ['MG049_SZ' '_frequency'];
% 
% mmil_ieeg_sensor_plot_v5(cfg)

end





