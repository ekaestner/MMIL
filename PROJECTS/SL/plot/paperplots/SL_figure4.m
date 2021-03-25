function SL_figure4

fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_mtc_1450'       'pap_lng_950'    'pap_lng_950'      'pap_con_950' }; %{ 'pap_mtc_1450' 'pap_mtc_1450' };
fcfg.eff_clm = [ 1                    1                2                  2 ];             %[ 1 2 ]; % 2 
fcfg.eff_col = { rgb('bright yellow') rgb('red')       rgb('blue')        rgb('cyan') };   %{ rgb('bright yellow') rgb('orangey yellow')-[0.35 0.2 0.08] }; % rgb('blue')
fcfg.eff_lbl = { 'StimulusMisMatch'   'VisualLanguage' 'AuditoryLanguage' 'AuditoryControl' }; %{ 'StimulusMisMatch' 'PostStimulusMisMatch' }; % 'Bigram Frequency'
fcfg.nsl_col = {rgb('purple')};

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = { 'pap_anv_1500' 'pap_anv_1500' 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1 2 3 ]; % 1

%% Middle Portion
% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

fcfg.sel_ele{2} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{2});
fcfg.sel_ele{3} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{3}); [~,~,rmv_ind] = intersect(fcfg.sel_ele{2},fcfg.sel_ele{3}); fcfg.sel_ele{3}(rmv_ind) = [];
fcfg.sel_ele{4} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{4}); [~,~,rmv_ind] = intersect(fcfg.sel_ele{2},fcfg.sel_ele{4}); fcfg.sel_ele{4}(rmv_ind) = [];

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
cfg.rad     = 3.25

cfg.sep_str         = ',';

cfg.inc_reg = { 'Lateral Occipital' 'Caudal Fusiform' 'Middle Fusiform' 'Middle ITG' 'Rostral ITG' ...
                'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' ...
                'Supramarginal' ...
                'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' 'Pars Triangularis' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'Rostral Middle Frontal' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/'];
cfg.sve_pre   = ['middle_pic'];
cfg.sep_str   = [','];

cfg.sve_img   = 'png';
mmil_ieeg_sensor_location_plot_v4(cfg);

% Put together balls for figure
pcfg = [];
pcfg.out_dir = cfg.sve_loc;
pcfg.sel_lbl = cfg.sel_lbl;
pcfg.col     = cfg.col;
pcfg.nsl_col = cfg.nsl_col;
mmil_loc_dot(pcfg)

%% Example Channels
% Example Channels
sbj = { 'NY439'          'NY439'          'NY590'          'NY598'          'NY590'           'NY540'           'NY439'           'NY537'          'NY439' };
chn = { 'LPTr09'         'LO08'           'G31'            'G14'            'G04'             'G51'             'G46'             'G23'            'RPT07' };
loc = { 'lhs_fus'        'lhs_stg'        'lhs_sup'        'lhs_pre'        'lhs_pop'         'rhs_pop'         'lhs_stg'         'rhs_stg'        'rhs_mtg' };
typ = [ 1                1                1                1                1                 1                 2                 2                2 ];
ylm = { [-1*10^8 3*10^8] [-1*10^7 6*10^7] [-1*10^7 6*10^7] [-1*10^8 3*10^8] [-3*10^6 18*10^6] [-3*10^7 18*10^7] [-3*10^7 18*10^7] [-1*10^8 3*10^8] [-1*10^7 6*10^7] };

tot_sbj = unique(sbj);

% Parameters
alt_eve = { 'trialinfo'                                               'trialinfo'                                             };
eve     = { [1 2 3]                                                   [1 2 4]                                                 };
col_ord = { { rgb('green') rgb('dark yellow') rgb('reddish grey') }   { rgb('green') rgb('dark yellow') rgb('bluish grey') } };

stt_dat = { { 'vis_nse_stt_msk_anv'  'vis_mtc_stt'}           { 'aud_nse_stt_msk_anv' 'vis_mtc_stt' }         };
stt_col = { { ft_stt_col(rgb('red')) ft_stt_col(rgb('dark yellow')) } { ft_stt_col(rgb('blue')) ft_stt_col(rgb('dark yellow')) }   };
stt_cmp = { { '0%4' '4%8' }                                           { '0%4' '4%8' }                                         };

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
        
        cfg.y_lim     =  ylm{plt_loc(iP)}; % 'maxmin'; % ylm{plt_loc(iP)};
        
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
        cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lines'];
        cfg.prefix     = [tot_sbj{iS} '_' 'match' '_' loc{plt_loc(iP)} '_'];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
        cfg.print_type = 'eps';
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end

% Line Legend
pcfg = [];
pcfg.out_dir = [ fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/'];
pcfg.lne_col = { rgb('green') rgb('dark yellow') rgb('reddish grey') rgb('bluish grey') } ;
pcfg.lne_lbl = { 'match'      'mismatch'         'false-font'        'noise-vocode' };
mmil_lne_leg(pcfg)
