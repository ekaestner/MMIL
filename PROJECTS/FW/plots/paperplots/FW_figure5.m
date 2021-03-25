function FW_figure5

%% Set-up Places to look
fcfg = [];

fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/epoch_data';
fcfg.tsk     = 'FW';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = 'pap_wrd_600';
fcfg.eff_clm = 3;
fcfg.eff_col = {rgb('bright red')};
fcfg.eff_lbl = {'Words'};

imp_loc = { 'caudal-fusiform' 'middle-fusiform' 'lateraloccipital' 'inferior-precentral' 'middle-precentral' 'parsopercularis' 'caudal-STG' 'middle-STG' 'middle-STG' 'supramarginal' };

% Find relevant channels
eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme '/' 'subjects' '/' 'total' '/' fcfg.eff_nme '_plt']);
sel_ele = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(1)+2)))+1,2); 
sel_ele = sel_ele(string_find(sel_ele,'NY226'));
sel_ele = cellfun(@(x) x(10:end),sel_ele,'uni',0);

tmp1 = mmil_readtext([fcfg.clr_fld '/' 'electrode_location_files' '/' 'NY226_FW' '/' 'output' '/' 'NY226_FW' '_' 'split' '_' 'lhs' '_' 'ecog'],'[, ]');
ecg_loc{1} = tmp1(:,[1 5:end]);
for iR = 1:size(ecg_loc{1},1)
    if isempty(ecg_loc{1}{iR,2})
        ecg_loc{1}{iR,2} = 'unknown';
    else
        num_hld = cellfun(@(x) str2num(x(1:end-1)),ecg_loc{1}(iR,string_find(ecg_loc{1}(iR,~cellfun(@isempty,ecg_loc{1}(iR,:))),'%')),'uni',0);
        [~,max_hld] = max([num_hld{:}]);
        ecg_loc{1}{iR,2} = ecg_loc{1}{iR,max_hld+1+(1*max_hld)};
    end
end
ecg_loc{1} = ecg_loc{1}(:,1:2);

rmv_ind = [];
for iI = 1:numel(sel_ele)
    if ~any(ismember(ecg_loc{1}(strcmpi(ecg_loc{1}(:,1),sel_ele{iI}),2),imp_loc))
        rmv_ind = [rmv_ind iI];
    end
end
  
sel_ele(rmv_ind) = [];

%% Central Picture
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/']);

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'NY226_FW' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'NY226_FW' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'NY226_FW' '/'  'output' '/' 'NY226_FW_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'NY226_FW' '/'  'output' '/' 'NY226_FW_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = sel_ele;
cfg.all_ele   = sel_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_ele = {};

cfg.sep_str         = ',';

cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/'];
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';

mmil_ieeg_sensor_location_plot_v3(cfg);

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_ele;
mmil_loc_dot(pcfg)

%% PLV plots
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'plv_plots']);

% Load
cnn_dat = load([fcfg.clr_fld '/' 'sig_chn' '/' 'connectivity' '/' 'plv' '/' 'NY226_FW' '/' 'NY226' '_FW_sig_plv.mat']);

% Plot
for iP1 = 1:numel(sel_ele)
    for iP2 = iP1:numel(sel_ele)
        if sum(strcmpi(cnn_dat.plt(1).labelcmb,[sel_ele{iP1} '--' sel_ele{iP2}])) > 0
            
            figure('Visible','off')
            subplot(2,1,1)
            colormap(hot(200));
            surf(cnn_dat.plt(1).tme,cnn_dat.plt(1).frq,squeeze(double(cnn_dat.plt(1).plv(strcmpi(cnn_dat.plt(1).labelcmb,[sel_ele{iP1} '--' sel_ele{iP2}]),:,:))),'FaceColor','interp', 'EdgeColor','none', 'FaceLighting','phong');
            set(gca,'View',[0 90]); axis tight;
            set(gca,'YGrid','off','XGrid','off','clim',[0 0.4]);
            
            cfg = [];
            cfg.jpg = 1;
            cfg.prn_typ = 'png';
            cfg.fle_nme = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'plv_plots'  '/' 'NY226_FW_' sel_ele{iP1} '_' sel_ele{iP2}];
            mmil_print_plot(cfg)

            cfg.prn_typ = 'eps';
            mmil_print_plot(cfg)
            
            close all
            
        end
    end
end

%% Electrode Location Plots
mmil_chk_dir([fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'line_plots']);

% Load
cfg = [];
cfg.load = 'yes';
cfg.file = [fcfg.dat_fld '/' 'NY226' '_FW_overall_data.mat'];
bcc_dat  = ft_func([],cfg);

% Plot
for iT = 1:numel(sel_ele)
    
    cfg = [];
    
    cfg.type      = 'chan';
    cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,sel_ele{iT}))};
    cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
    cfg.dat_loc   = 1;
    
    cfg.plt_dim   = [1 1];
    
    cfg.y_lim     = 'maxmin';
    
    cfg.lgd       = 0;
    cfg.std_err   = 1;
    
    cfg.alt_eve = 'trialinfo';
    cfg.eve     = [3 5 6];
    cfg.lnstyle.col_ord = {rgb('bright red') rgb('purple') rgb('reddish grey')} ;
    
    cfg.stt_dat = { 'vis_stm_01'               'vis_wrd_msk_01'         'vis_wrd_ffn_msk_01'};
    cfg.stt_col = { { ft_stt_col(rgb('black')) ft_stt_col(rgb('black')) ft_stt_col(rgb('black')) } };
    cfg.stt_lab = 'stt_lab';
    cfg.stt_cmp = { { '0%5'                   '3vs5'          '3vs6'   } };
    
    cfg.v_lne       = [0 0.2 0.4];
    cfg.v_lne_wdt    = [3 1 1];
    cfg.v_lne_col   = {rgb('red') rgb('black') rgb('black')};
    
    cfg.x_lim = [-0.2 0.8];
    
    cfg.print      = 1;
    cfg.nofig      = 1;
    cfg.print_type = 'png';
    cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'line_plots' '/' ];
    cfg.prefix     = ['NY226_FW' '_' sel_ele{iT} '_'];
    
    mmil_ieeg_sensor_plot_v5(cfg)
    
    cfg.print_type = 'eps';
    mmil_ieeg_sensor_plot_v5(cfg)
    
end
    
    
    
    
    
    
    
    
    

end