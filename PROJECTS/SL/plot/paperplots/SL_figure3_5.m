0function SL_figure3_5

fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_phn_950' 'pap_phn_950' };
fcfg.eff_clm = [ 1 2 ]; % 2 
fcfg.eff_col = { rgb('bright red') rgb('bright blue') }; 
fcfg.eff_lbl = { 'Letter-Selective' 'Phoneme-Selective' }; 
fcfg.nsl_col = { rgb('purple') };

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = { 'pap_anv_1500' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 4 ]; % 1

%% Middle Portion
% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end
fcfg.sel_ele{2}([20 63]) = [];

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
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/'];
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
sbj = { 'NY439'           'NY439'           'NY439'           'NY590'           'NY598'           'NY590'            'NY537'           'NY569' };
chn = { 'LPTr05'          'RPT02'           'GR55'            'G07'             'G14'             'G03'              'G28'             'G47' };
loc = { 'lhs_fus'         'rgh_fus'         'lhs_stg'         'lhs_pre'         'lhs_pre'         'lhs_opc'          'rgh_stg'         'rgh_pre' };
typ = [ 1                 1                 2                 1                 2                 1                  2                 1 ];
ylm = { [-2*10^8 20*10^8] [-1*10^8 10*10^8] [-1*10^8 10*10^8] [-1*10^8 10*10^8] [-1*10^8 10*10^8] [-0.2*10^8 2*10^8] [-1*10^8 10*10^8] [-0.2*10^8 2*10^8] };

tot_sbj = unique(sbj);

% Parameters
phn_col = distinguishable_colors(42);

for iC = 1:12
    phn_vis_col{iC} = phn_col(dsearchn(phn_col,rgb('red')),:);
    phn_col(dsearchn(phn_col,rgb('red')),:) = [];
end

for iC = 1:12
    phn_aud_col{iC} = phn_col(dsearchn(phn_col,rgb('blue')),:);
    phn_col(dsearchn(phn_col,rgb('blue')),:) = [];
end

nme     = { 'VisualLetter'                                         'AuditoryPhoneme'           };

alt_eve = { 'vis_con'                                              'aud_con'                   };
eve     = { [1:12]                                                 [1:12]                      };
col_ord = { phn_vis_col                                            phn_aud_col };

stt_dat = { { 'vis_con_stt' }                                      { 'aud_con_stt' }           };
stt_col = { { ft_stt_col(rgb('red')) }                             { ft_stt_col(rgb('blue')) } };
stt_cmp = { { '0%5' }                                              { '0%5' }                   };

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
        cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/' 'lines'];
        cfg.prefix     = [tot_sbj{iS} '_' 'lexical_overlap' '_' loc{plt_loc(iP)} '_'];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
        cfg.print_type = 'eps';
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end

% Line Legend
% pcfg = [];
% pcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/'];
% pcfg.lne_col = {rgb('') } ;
% pcfg.lne_lbl = {''      };
% mmil_lne_leg(pcfg)

%% Middle Portion Presentation
fcfg = [];
fcfg.chn_crr = 1;
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical';
fcfg.dat_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/epoch_data';
fcfg.tsk     = 'SL';

fcfg.eff_typ = 'hgp';
fcfg.eff_nme = { 'pap_phn_950' 'pap_phn_950' };
fcfg.eff_clm = [ 1 2 ]; % 2 
fcfg.eff_col = { rgb('bright red') rgb('bright blue') }; 
fcfg.eff_lbl = { 'Letter-Selective' 'Phoneme-Selective' }; 
fcfg.nsl_col = { rgb('purple') };

fcfg.ovr_typ = 'hgp'; 
fcfg.ovr_nme = { 'pap_phn_950' 'pap_phn_950' }; % 'pap_rsp_600'
fcfg.ovr_clm = [ 1             2 ]; % 1

% Selective Electrodes
for iEF = 1:numel(fcfg.eff_clm)
    eff_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.eff_typ '/' 'ecog' '/' 'split' '/' fcfg.eff_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.eff_nme{iEF} '_plt']);
    
    fcfg.sel_ele{iEF} = eff_txt(find(cell2mat(eff_txt(2:end,fcfg.eff_clm(iEF)+2)))+1,2);
end

fcfg.sel_ele{3} = fcfg.sel_ele{2}([20]);
fcfg.eff_col{3} = rgb('orange');
fcfg.eff_lbl{3} = 'bimodal';

fcfg.sel_ele{1}([4])           = [];
fcfg.sel_ele{2}([20 21 62 63]) = [];

% Overall Electrodes
% for iEF = 1:numel(fcfg.ovr_clm)
%     ovr_txt = mmil_readtext([fcfg.clr_fld '/' 'sig_chn' '/' fcfg.ovr_typ '/' 'ecog' '/' 'split' '/' fcfg.ovr_nme{iEF} '/' 'subjects' '/' 'total' '/' fcfg.ovr_nme{iEF} '_plt']);
%     
%     fcfg.all_ele{iEF} = ovr_txt(find(cell2mat(ovr_txt(2:end,fcfg.ovr_clm(iEF)+2)))+1,2);
% end
fcfg.all_ele = unique(cat(1,fcfg.sel_ele{:}));

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
                'Caudal STG' 'Middle STG' 'Rostral STG' ...
                'Pars Opercularis' };

cfg.sve_img   = 'eps';
cfg.sve_loc   = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/'];
cfg.sve_pre   = ['middle_pic_presentation'];
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

end