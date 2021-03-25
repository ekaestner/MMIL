function SL_figure5

fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_lng_950' 'pap_con_950' 'pap_con_950' };
fcfg.eff_clm = [ 1 1 2 ]; % 2 
fcfg.eff_col = { rgb('red') rgb('dark magenta')-[0.125 0 0.085] rgb('cyan') }; % rgb('blue')
fcfg.eff_lbl = { 'Text Selective' 'False-Font Selective' 'Noise-Vocoded Selective' }; % 'Bigram Frequency'
fcfg.nsl_col = { rgb('purple') };

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = { 'pap_anv_1500' 'pap_anv_1500' 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1 2 3 ]; % 1

%% Middle Portion
% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

fcfg.sel_ele{1} = unique( [ intersect(fcfg.sel_ele{1},fcfg.sel_ele{2}) intersect(fcfg.sel_ele{1},fcfg.sel_ele{3}) ] );

% Overall Electrodes
for iEF = 1:numel(fcfg.ovr_clm)
    ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme{iEF} '_plt']);
    
    fcfg.all_ele{iEF} = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm(iEF)+2)))+1,2);
end
fcfg.all_ele = unique(cat(1,fcfg.all_ele{:}));

% Plot
cfg = [];

cfg.hms = {'lhs' 'rhs'};
cfg.hem = {'lhs' 'rhs'};

cfg.pial_mat  = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'lh.pial']               [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'surf' '/' 'rh.pial']};
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog']};

cfg.sml_vew = 1;

cfg.sel_ele   = fcfg.sel_ele;
cfg.all_ele   = fcfg.all_ele;

cfg.sel_lbl = fcfg.eff_lbl;
cfg.col     = fcfg.eff_col;
cfg.nsl_col = fcfg.nsl_col;

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/'];
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_ele;
mmil_loc_dot(pcfg)

%% Example Channels
% Example Channels
sbj = { 'NY439'           'NY598'           'NY439'           'NY540'   };
chn = { 'GR08'            'G28'             'LO08'            'G51'     };
loc = { 'lft_ifg'         'lft_pre'         'lft_stg'         'rgh_ifg' };
typ = [ 1                 1                 1                 1 ];
ylm = { [-6*10^7 16*10^7] [-3*10^7 8*10^7] [-3*10^7 8*10^7] [-6*10^7 16*10^7] };

tot_sbj = unique(sbj);

% Parameters
nme     = { 'controls'                                             };

alt_eve = { 'lng_tot_nse'                                          };
eve     = { [601                603           604]                 };
col_ord = { { rgb('light purple') rgb('reddish grey') rgb('bluish grey')} };
                     
stt_dat = { { 'vis_nse_stt'          'aud_nse_stt' }       };
stt_col = { { ft_stt_col(rgb('red')) ft_stt_col(rgb('cyan')) }     };
stt_cmp = { { '0%4'                  '4%8' }                       };

% Put together plot
for iS = 1:numel(tot_sbj)
    
    % Load
    cfg = [];
    cfg.load = 'yes';
    cfg.file = [fcfg.dat_fld '/' tot_sbj{iS} '_SL_overall_data.mat'];
    bcc_dat  = ft_func([],cfg);
        
    plt_loc = find(ismember(sbj,[tot_sbj{iS}]));
       
    for iP = 1:numel(plt_loc)
        
        cfg = [];
        
        cfg.type      = 'chan';
        cfg.chn_grp   = {find(strcmpi(bcc_dat.(bcc_dat.data_name{2}).cfg.alt_lab.label,chn{plt_loc(iP)}))};
        cfg.dat       = { bcc_dat.(bcc_dat.data_name{2}) };
        cfg.dat_loc   = 1;
        
        cfg.plt_dim   = [1 1];
        
        cfg.y_lim     = ylm{plt_loc(iP)}; % ylm{plt_loc(iP)};
        
        cfg.lgd       = 0;
        cfg.std_err   = 1;
        
        cfg.alt_eve = alt_eve{typ(plt_loc(iP))};
        cfg.eve     = eve{typ(plt_loc(iP))};
        cfg.lnstyle.col_ord = col_ord{typ(plt_loc(iP))};
        
        cfg.stt_dat = stt_dat(typ(plt_loc(iP)));
        cfg.stt_col = stt_col(typ(plt_loc(iP)));
        cfg.stt_lab = 'stt_lab';
        cfg.stt_cmp = stt_cmp(typ(plt_loc(iP)));
        
        cfg.v_lne       = {[0         0.450       0.900]};
        cfg.v_lne_wdt   = {[10        10          2.5]};
        cfg.v_lne_col   = {rgb('red') rgb('blue') rgb('black')};
        
        cfg.x_lim = [-0.1 1.2];
        
        cfg.print      = 1;
        cfg.nofig      = 1;
        cfg.print_type = 'png';
        cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'lines'];
        cfg.prefix     = [tot_sbj{iS} '_' 'control' '_' loc{plt_loc(iP)} '_'];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
        cfg.print_type = 'eps';
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end

% Line Legend
pcfg = [];
pcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/'];
pcfg.lne_col = {rgb('light purple') rgb('reddish grey') rgb('bluish grey')} ;
pcfg.lne_lbl = {'Language'          'False-Font'        'NoiseVocoding'};
mmil_lne_leg(pcfg)
