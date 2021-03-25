clear; clc;

%% VISUAL
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'    'pap_vis_rep'    'pap_vis_lex' };
fcfg.ana_clm = [ 1                1                1    ];
fcfg.eff_clm = { [1]              [1]              [1]    };
fcfg.ana_col = { 'dark magenta'   'bright red'     'orangey yellow'  };
fcfg.ana_lbl = { 'TextResponsive' 'TextRepetition' 'TextLexical' };

fcfg.sub_tme = [ 0.000           0.000             0.000   ];

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

%% Auditory HGP
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_aud_act'        'pap_aud_rep'        'pap_aud_lex'};
fcfg.ana_clm = [ 1                    1                    1 ];
fcfg.eff_clm = { [1]                  [1]                  [1] };
fcfg.ana_col = { 'dark purple'        'bright blue'        'seafoam' };
fcfg.ana_lbl = { 'AuditoryResponsive' 'AuditoryRepetition' 'AuditoryLexical' };

fcfg.sub_tme = [ 0.000           0.000                     0.000 ];

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

%% VISUAL N400
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'lfp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'    'pap_vis_n400'   'pap_vis_lex' };
fcfg.ana_clm = [ 1                1                1    ];
fcfg.eff_clm = { [1]              [1]              [1]    };
fcfg.ana_col = { 'dark magenta'   'bright red'     'orangey yellow' };
fcfg.ana_lbl = { 'TextResponsive' 'TextRepetition' };

fcfg.sub_tme = [ 0.000           0.000             0.000   ];

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

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'VisualN400' '/'];
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

%% Auditory N400
% Setup
fcfg = [];
fcfg.dat_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data/';
fcfg.clr_fld     = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';

fcfg.typ     = { 'lfp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_aud_act'        'pap_aud_n400'       'pap_aud_n400' };
fcfg.ana_clm = [ 1                    1                    1   ];
fcfg.eff_clm = { [1]                  [1]                  [1] };
fcfg.ana_col = { 'dark purple'        'bright blue'        'seafoam' };
fcfg.ana_lbl = { 'AuditoryResponsive' 'AuditoryRepetition' 'AuditoryLexical' };

fcfg.sub_tme = [ 0.000           0.000                     0.000];

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

fcfg.out_dir = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'AuditoryN400' '/'];
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