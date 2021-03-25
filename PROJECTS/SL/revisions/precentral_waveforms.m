%% Example Channels
% Example Channels
sbj = { 'NY439'            'NY590'            'NY591'            'NY598'            };
chn = { 'GR35'             'G07'              'G13'              'G06'              };
loc = { 'NY439'            'NY590'            'NY591'            'NY598'          };
typ = [ 1                  1                  1                  1                  ];
ylm = { [-1*10^8 8*10^8]   [-0.5*10^8 4*10^8] [-1*10^8 8*10^8]   [-1*10^7 8*10^7]   };

tot_sbj = unique(sbj);

% Parameters
alt_eve = { 'lng_tot_nse'                             };
eve     = { [ 601                 603                 ] };
col_ord = { { rgb('light purple') rgb('reddish grey') } };

stt_dat = { { 'vis_nse_stt_msk_anv'  }                };
stt_col = { { ft_stt_col(rgb('white')) }              };
stt_cmp = { { '0%4'                     }             };

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
        cfg.outdir     = [fcfg.clr_fld '/' 'manuscript' '/' 'revisions' '/'];
        cfg.prefix     = [tot_sbj{iS} '_' 'lexical_overlap' '_' loc{plt_loc(iP)} '_'];
        
        mmil_ieeg_sensor_plot_v5(cfg)
        
        cfg.print_type = 'eps';
        mmil_ieeg_sensor_plot_v5(cfg)
        
    end
    
end
