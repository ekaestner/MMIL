function SL_figure3

fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_lng_950' 'pap_lng_950' };
fcfg.eff_clm = [ 1 2 ]; % 2
fcfg.eff_col = { rgb('red')-[0.3 0 0] rgb('bright blue')-[0 0.13 0.3] }; % rgb('blue')
fcfg.eff_lbl = { 'Text Selective' 'Voice Selective' }; % 'Bigram Frequency'
fcfg.nsl_col = { rgb('dark purple')};

fcfg.ovr_typ = 'hgp';
fcfg.ovr_nme = { 'pap_anv_1500' 'pap_anv_1500' 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1 2 3 ]; % 1

%% Middle Portion
% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

fcfg.sel_ele{3} = intersect(fcfg.sel_ele{1},fcfg.sel_ele{2});
fcfg.sel_ele{1} = setxor(fcfg.sel_ele{1},fcfg.sel_ele{3});
fcfg.sel_ele{2} = setxor(fcfg.sel_ele{2},fcfg.sel_ele{3});

fcfg.eff_col{3} = rgb('orange'); %rgb('bright violet')+[0.1 0.03 0];

fcfg.eff_lbl{3} = 'Bi-Modal Selective';

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
cfg.elec_text = {[fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_lhs_ecog_lat']      [fcfg.clr_fld '/' 'electrode_location_files' '/' 'total' '/'  'output' '/' 'total_rhs_ecog_lat']};

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
                'Pars Opercularis' 'Pars Triangularis' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/'];
cfg.sve_pre   = ['middle_pic_presentation_language'];
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
sbj = { 'NY591'            'NY439'            'NY439'            'NY439'            'NY523'            'NY439'            'NY523'             'NY591'            'NY591'            'NY590'             'NY540'          'NY540'            'NY537'};
chn = { 'PT01'             'GR12'             'GR35'             'G44'              'G28'              'GR53'             'GS55'              'G39'              'G43'              'G03'               'PT04'           'G51'              'G28' };
loc = { 'lft_fus'          'lft_pop'          'lft_pre'          'lft_stg'          'lft_pos'          'lft_stg'          'lft_pop'           'lft_sup'          'lft_stg'          'lft_pop'           'rgh_fus'        'rgh_pop'          'rgh_stg'};
typ = [ 1                  1                  1                  1                  2                  2                  2                   3                  3                  3                   1                1                  2 ];
ylm = { [-1*10^8 8*10^8]   [-0.5*10^8 4*10^8] [-1*10^8 8*10^8]   [-1*10^7 8*10^7]   [-0.5*10^8 4*10^8] [-0.5*10^8 4*10^8] [-0.5*10^8 4*10^8]  [-0.5*10^7 4*10^7] [-0.5*10^7 4*10^7] [-0.5*10^8 4*10^8]  [-1*10^8 8*10^8] [-0.5*10^8 4*10^8] [-0.5*10^8 4*10^8] };

tot_sbj = unique(sbj);

% Parameters
alt_eve = { 'lng_tot_nse'                                                  'lng_tot_nse'                                                  'lng_tot_nse' };
eve     = { [ 601           603                 604 ]                      [ 601           603                 604 ]                      [ 601           603                 604 ] };
col_ord = { { rgb('light purple') rgb('reddish grey') rgb('bluish grey') } { rgb('light purple') rgb('reddish grey') rgb('bluish grey') } { rgb('light purple') rgb('reddish grey') rgb('bluish grey') } };

stt_dat = { { 'vis_nse_stt_msk_anv'  }                               { 'aud_nse_stt_msk_anv'  }                               { 'vis_nse_stt_msk_anv'  'aud_nse_stt_msk_anv'} };
stt_col = { { ft_stt_col(rgb('bright red')) }                        { ft_stt_col(rgb('bright blue')) }                       { ft_stt_col(rgb('red')) ft_stt_col(rgb('bright blue')) }  };
stt_cmp = { { '0%4'                     }                            { '0%4'                     }                            { '0%4'                  '4%8' } };

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
        
        cfg.y_lim     = ylm{plt_loc(iP)};
        
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
        cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'lines'];
        cfg.prefix     = [tot_sbj{iS} '_' 'lexical_overlap' '_' loc{plt_loc(iP)} '_'];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
        cfg.print_type = 'eps';
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end

% Line Legend
pcfg = [];
pcfg.out_dir = [ fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/'];
pcfg.lne_col = { rgb('light purple') rgb('reddish grey') rgb('bluish grey')     } ;
pcfg.lne_lbl = { 'Language'    'False-Font'        'NoiseVocoded' };
mmil_lne_leg(pcfg)
