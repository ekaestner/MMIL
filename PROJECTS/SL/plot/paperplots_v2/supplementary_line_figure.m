fcfg.dat_fld = '';
fcfg.clr_fld = '';

%% Example Channels
% Example Channels
sbj = { 'NY439'           'NY439'           'NY439'           'NY590'           'NY598'           'NY590'            };
chn = { 'LPTr05'          'RPT02'           'GR55'            'G07'             'G14'             'G03'              };
loc = { 'lhs_fus'         'rgh_fus'         'lhs_stg'         'lhs_pre'         'lhs_pre'         'lhs_opc'          };
typ = [ 1                 1                 1                 1                 1                 1                  ];
ylm = { [-2*10^8 20*10^8] [-1*10^8 10*10^8] [-1*10^8 10*10^8] [-1*10^8 10*10^8] [-1*10^8 10*10^8] [-0.2*10^8 2*10^8] };

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