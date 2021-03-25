clear; clc;

%% Activation
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'    'pap_aud_act'   };
fcfg.ana_clm = [ 1                1                ];
fcfg.eff_clm = { [1]              [1]              };
fcfg.ana_col = { 'dark red'        'dark blue'     };
fcfg.ana_lbl = { 'TextResponsive' 'Auditory Responsive' };

fcfg.sub_tme = [ 0.000           0.000             ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    {'lateraloccipital'} ...
    {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
    {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
    {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
    {'supramarginal'} ...
    {'inferior-precentral' 'middle-precentral'} ...
    {'parstriangularis'} ...
    {'parsopercularis'} ...
    {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Supramarginal' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' ...
    'middle-middlefrontal' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'VisualHGP' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.ylm_exc = 0.060;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [2 1 ; 3 2 ]; % [ref_tbl tme_tab]
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Repetition HGP
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_rep'       'pap_aud_rep'        };
fcfg.ana_clm = [ 1                   1                    ];
fcfg.eff_clm = { [1]                 [1]                  };
fcfg.ana_col = { 'bright red'        'bright blue'        };
fcfg.ana_lbl = { 'VisualRepetition'  'AuditoryRepetition' };

fcfg.sub_tme = [ 0.000           0.000                    ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    {'lateraloccipital'} ...
    {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
    {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
    {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
    {'supramarginal'} ...
    {'inferior-precentral' 'middle-precentral'} ...
    {'parstriangularis'} ...
    {'parsopercularis'} ...
    {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Supramarginal' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' ...
    'middle-middlefrontal' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'AuditoryHGP' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [2 1 ; 3 2 ]; % [ref_tbl tme_tab]
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Activation N400
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'lfp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'    'pap_aud_act'   };
fcfg.ana_clm = [ 1                1                ];
fcfg.eff_clm = { [1]              [1]              };
fcfg.ana_col = { 'magenta'        'dark blue'     };
fcfg.ana_lbl = { 'TextResponsive' 'Auditory Responsive' };

fcfg.sub_tme = [ 0.000           0.000             ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    {'lateraloccipital'} ...
    {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
    {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
    {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
    {'supramarginal'} ...
    {'inferior-precentral' 'middle-precentral'} ...
    {'parstriangularis'} ...
    {'parsopercularis'} ...
    {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Supramarginal' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' ...
    'middle-middlefrontal' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'VisualHGP' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [2 1 ; 3 2 ]; % [ref_tbl tme_tab]
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Repetition N400
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_rep'       'pap_aud_rep'        };
fcfg.ana_clm = [ 1                   1                    ];
fcfg.eff_clm = { [1]                 [1]                  };
fcfg.ana_col = { 'bright red'        'bright blue'        };
fcfg.ana_lbl = { 'VisualRepetition'  'AuditoryRepetition' };

fcfg.sub_tme = [ 0.000           0.000                    ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    {'lateraloccipital'} ...
    {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
    {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
    {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
    {'supramarginal'} ...
    {'inferior-precentral' 'middle-precentral'} ...
    {'parstriangularis'} ...
    {'parsopercularis'} ...
    {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Supramarginal' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' ...
    'middle-middlefrontal' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'AuditoryHGP' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [2 1 ; 3 2 ]; % [ref_tbl tme_tab]
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Subtraction HGP
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_bim_act'         'pap_bim_rep'        };
fcfg.ana_clm = [ 1                     1                    ];
fcfg.eff_clm = { [1]                   [1]                  };
fcfg.ana_col = { 'dark purple'         'bright purple'        };
fcfg.ana_lbl = { 'Bi-Modal Activation' 'Bi-Modal Repetition' };

fcfg.sub_tme = [ 0.000                  0.000                ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    {'lateraloccipital'} ...
    {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
    {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
    {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
    {'supramarginal'} ...
    {'inferior-precentral' 'middle-precentral'} ...
    {'parstriangularis'} ...
    {'parsopercularis'} ...
    {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Supramarginal' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' ...
    'middle-middlefrontal' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'BiModalHGP' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.ylm_typ = 'maxmin';
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [2 1 ; 3 2 ]; % [ref_tbl tme_tab]
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
%         mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Time Table
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'    'pap_vis_rep'      'pap_aud_act'         'pap_aud_rep'        'pap_bim_act'         'pap_bim_rep'         'pap_vis_bim_act' 'pap_vis_bim_act' 'pap_vis_bim_act' 'pap_aud_bim_act' 'pap_aud_bim_act' 'pap_vis_bim_rep' 'pap_vis_bim_rep' 'pap_vis_bim_rep' 'pap_aud_bim_rep' 'pap_aud_bim_rep' };
fcfg.ana_clm = [ 1                1                  1                     1                    1                     1                     1                 2                 3                 1                 2                 1                 2                 3                 1                 2 ];
fcfg.eff_clm = { [1]              [1]                [1]                   [1]                  [1]                   [1]                   [1]               [2]               [3]               [1]               [2]               [1]               [2]               [3]               [1]               [2]};
fcfg.ana_col = { 'dark red'       'bright red'       'dark blue'           'bright blue'        'dark purple'         'bright purple'       'red'             'purple'          'black'           'blue'            'purple'          'red'             'purple'          'black'           'blue'            'purple' };
fcfg.ana_lbl = { 'TextResponsive' 'VisualRepetition' 'Auditory Responsive' 'AuditoryRepetition' 'Bi-Modal Activation' 'Bi-Modal Repetition' 'VisualUni'       'VisualBi'        'TotalUni'        'AuditoryUni'     'AuditoryBi'      'VisualUni'       'VisualBi'        'TotalUni'        'AuditoryUni'     'AuditoryBi'};

fcfg.sub_tme = [ 0.000           0.000               0.000                 0.000                0.000                 0.000                 0.000             0.000             0.000             0.000             0.000             0.000             0.000             0.000             0.000             0.000];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    {'lateraloccipital'} ...
    {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
    {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
    {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
    {'supramarginal'} ...
    {'inferior-precentral' 'middle-precentral'} ...
    {'parstriangularis'} ...
    {'parsopercularis'} ...
    {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Supramarginal' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' ...
    'middle-middlefrontal' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
%         mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        tme_hld = -2:0.020:2;
        for iR = 1:size(tme_tab{iH},1)
            for iC = 2:size(tme_tab{iH},2)
                if ~isempty(tme_tab{iH}{iR,iC})
                    for iTM = 1:numel(tme_tab{iH}{iR,iC})
                        tme_tab{iH}{iR,iC}(iTM) = tme_hld(dsearchn(tme_hld',tme_tab{iH}{iR,iC}(iTM)));
                    end
                end
            end
        end
        
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [1 1 ; 2 2 ; 3 3 ; 4 4 ; 5 5 ;    6 6]; % [ref_tbl tme_tab]
        cfg.ylm_typ = {''    ''    ''    ''    'maxmin' 'maxmin'};
        cfg.ylm_exc = [0.060 0.060 0.060 0.060 -999     -999];
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        cfg.ylm_typ = {''    ''    ''    ''    'maxmin' 'maxmin'};
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Time Table LFP
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'lfp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_rep_nob'   'pap_aud_rep_nob'    'pap_bim_rep_nob'     'pap_vis_bim_rep_nob' 'pap_vis_bim_rep_nob' 'pap_vis_bim_rep_nob' 'pap_aud_bim_rep_nob' 'pap_aud_bim_rep_nob'  };
fcfg.ana_clm = [ 1                   1                    1                     1                     2                     3                     1                     2 ];
fcfg.eff_clm = { [1]                 [1]                  [1]                   [1]                   [2]                   [3]                   [1]                   [2] };
fcfg.ana_col = { 'bright red'        'bright blue'        'bright purple'       'red'                 'purple'              'black'               'blue'                'purple'     };
fcfg.ana_lbl = { 'VisualRepetition'  'AuditoryRepetition' 'Bi-Modal Repetition' 'VisualUni'           'VisualBi'            'TotalUni'            'AuditoryUni'         'AuditoryBi' };

fcfg.sub_tme = [ 0.000               0.000                0.000                 0.000             0.000             0.000             0.000             0.000 ];

fcfg.inc_reg = { { 'lateraloccipital' 'caudal-fusiform' 'middle-fusiform' 'caudal-ITG' 'middle-ITG' 'rostral-ITG' } ...
                 { 'parsopercularis' 'parstriangularis' 'parsorbitalis' 'middle-middlefrontal' } ...
                 { 'inferior-postcentral' 'inferior-precentral' 'middle-precentral' 'caudal-STG' 'middle-STG' 'caudal-MTG' } };
fcfg.nme_reg = { 'Ventral' 'Frontal' 'Lateral' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' 'rostral-fusiform'} ...
    {'lateraloccipital'} ...
    {'caudal-ITG' 'middle-ITG' 'rostral-ITG'} ...
    {'caudal-MTG' 'middle-MTG' 'rostral-MTG'} ...
    {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
    {'supramarginal'} ...
    {'inferior-precentral' 'middle-precentral'} ...
    {'parstriangularis'} ...being 
    {'parsopercularis'} ...
    {'middle-middlefrontal'} };
fcfg.reg_nme = {'Fusiform' ...
    'Lateral Occipital' ...
    'ITG' ...
    'MTG' ...
    'STG' ...
    'Supramarginal' ...
    'Precentral' ...
    'Pars Triangularis' ...
    'Pars Operculatris' ...
    'middle-middlefrontal' };

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTableLFP' '/'];
mmil_chk_dir(fcfg.out_dir)

% Plot Construction
for iT = 1:numel(fcfg.typ)
    
    %% Load Data
    cfg = [];
    cfg.clr_fld = fcfg.clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sub_tme = fcfg.sub_tme;
    cfg.sbj_nme = 'total';
    [tme_tab,tme_plt] = mmil_tme_tab(cfg);
    
    save([fcfg.out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    %% Plot Figures
    for iH = 1:numel(fcfg.hms)
        for iP = 1%numel(fcfg.nme_reg)
            
            cfg = [];
            cfg.typ = fcfg.typ{iT};
            cfg.reg = fcfg.inc_reg{iP};
            cfg.nme = fcfg.nme_reg{iP};
            cfg.hms = fcfg.hms{iH};
            cfg.ana_col = fcfg.ana_col;
            cfg.ana_lbl = fcfg.ana_lbl;
            cfg.out_dir = fcfg.out_dir;
            %             mmil_tme_plt(cfg,tme_plt{iH})
            
        end
        
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
%         mmil_tme_plt(cfg,tme_plt{iH})
        
        %% Tables
        tme_hld = -2:0.020:2;
        for iR = 1:size(tme_tab{iH},1)
            for iC = 2:size(tme_tab{iH},2)
                if ~isempty(tme_tab{iH}{iR,iC})
                    for iTM = 1:numel(tme_tab{iH}{iR,iC})
                        tme_tab{iH}{iR,iC}(iTM) = tme_hld(dsearchn(tme_hld',tme_tab{iH}{iR,iC}(iTM)));
                    end
                end
            end
        end
        
        cfg.ref_tbl = [fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table2.csv'];
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '.csv'];
        cfg.ttl_col = 2;
        cfg.col_crr = [1 1 ; 2 2 ; 3 3 ; 4 4 ; 5 5 ;    6 6]; % [ref_tbl tme_tab]
        cfg.ylm_typ = {''    ''    ''    ''    'maxmin' 'maxmin'};
        cfg.ylm_exc = [0.060 0.060 0.060 0.060 -999     -999];
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
        cfg.out_tbl = [fcfg.out_dir '/' 'eff_tme_tbl_' fcfg.hms{iH} '_combined.csv'];
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.reg_nme = fcfg.reg_nme;
        cfg.ylm_typ = {''    ''    ''    ''    'maxmin' 'maxmin'};
        mmil_crt_tme_tab(cfg,tme_tab{iH})
        
    end
    
end

%% Table Construction
% lhs
hgp_eff = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/' 'eff_tme_tbl_' fcfg.hms{1} '_combined.csv']);
hgp_eff{2,3} = 'Text-Responsive HGP';
hgp_eff{2,4} = 'Text-Repetition HGP';
hgp_eff{2,5} = 'Voice-Responsive HGP';
hgp_eff{2,6} = 'Voice-Repetition HGP';
hgp_eff{2,7} = 'Bimodal-Responsive HGP';
hgp_eff{2,8} = 'Bimodal-Repetition HGP';

lfp_eff = mmil_readtext([fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTableLFP' '/' 'eff_tme_tbl_' fcfg.hms{1} '_combined.csv']);
lfp_eff{2,4} = 'Auditory Repetition';
lfp_eff{2,5} = '?';


tme_tbl = [hgp_eff lfp_eff];
ttl_col = 2;


col = {'' 'Text-Responsive HGP'    'Text-Repetition HGP'    'Text Repetition'  ...
       '' 'Voice-Responsive HGP'   'Voice-Repetition HGP'   'Auditory Repetition' }; % wonky LFP labeling

out_tbl = tme_tbl(:,1);
for iC = 1:numel(col)
    if ~strcmpi(col{iC},'')
        out_tbl(:,iC+1) = tme_tbl(:,string_find(tme_tbl(ttl_col,:),col(iC)));
    end
end

cell2csv([fcfg.clr_fld '/' 'manuscript' '/' 'tables' '/' 'table4.csv'],out_tbl);



















