clear; clc;

dta_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/epoch_data';
clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical';

%% Overall Plot
out_dir = [clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Modal_amplitude'];

% Setup
fcfg = [];

fcfg.sig_nme = 'pap_vis_act';
fcfg.sig_col = [ 1 ];

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'  'pap_aud_act' };
fcfg.ana_clm = [ 1              1 ];
fcfg.eff_clm = { [1]            [1] };
fcfg.ana_col = { 'dark red'     'blue' };
fcfg.ana_lbl = { 'VisualActive' 'AuditoryActive' };

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

% Get amp values 
for iT = 1:numel(fcfg.typ)

    cfg = [];
    cfg.clr_fld = clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sbj_nme = 'total';
    [amp_tab,amp_plt] = mmil_amp_tab(cfg);
       
    for iH = 1:numel(fcfg.hms)
        % Make plot
        cfg = [];
        cfg.typ = fcfg.typ{iT};
        cfg.hms = fcfg.hms{iH};
        cfg.ana_col = fcfg.ana_col;
        cfg.ana_lbl = fcfg.ana_lbl;
        cfg.out_dir = out_dir;
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.exc_amp = 100;
        mmil_amp_plt(cfg,amp_plt{iH})
    end
    
end

%% Bi-Modal Plot
out_dir = [clr_fld '/' 'manuscript' '/' 'figure6' '/' 'BiModal_amplitude'];

% Setup
fcfg = [];

fcfg.sig_nme = 'pap_vis_act';
fcfg.sig_col = [ 1 ];

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_bim_act'  };
fcfg.ana_clm = [ 1              ];
fcfg.eff_clm = { [1]            };
fcfg.ana_col = { 'dark purple'     };
fcfg.ana_lbl = { 'VisualActive' };

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

% Get amp values 
for iT = 1:numel(fcfg.typ)

    cfg = [];
    cfg.clr_fld = clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sbj_nme = 'total';
    [amp_tab,amp_plt] = mmil_amp_tab(cfg);
    
    for iH = 1:numel(fcfg.hms)
        % Make plot
        cfg = [];
        cfg.typ = fcfg.typ{iT};
        cfg.hms = fcfg.hms{iH};
        cfg.ana_col = fcfg.ana_col;
        cfg.ana_lbl = fcfg.ana_lbl;
        cfg.out_dir = out_dir;
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.exc_amp = [-100 100];
        cfg.ylm_typ = 'maxmin';
        mmil_amp_plt(cfg,amp_plt{iH})
    end
    
end

%% Overall Table
out_dir = [clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Table_amplitude'];

% Setup
fcfg = [];

fcfg.sig_nme = 'pap_vis_act';
fcfg.sig_col = [ 1 ];

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'  'pap_aud_act' };
fcfg.ana_clm = [ 1              1 ];
fcfg.eff_clm = { [1]            [1] };
fcfg.ana_col = { 'dark red'     'blue' };
fcfg.ana_lbl = { 'VisualActive' 'AuditoryActive' };

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

% Setup
fcfg = [];

fcfg.sig_nme = 'pap_vis_act';
fcfg.sig_col = [ 1 ];

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'pap_vis_act'  'pap_aud_act'  'pap_bim_act'    'pap_vis_bim_act' 'pap_vis_bim_act' 'pap_vis_bim_act' 'pap_aud_bim_act' 'pap_aud_bim_act'};
fcfg.ana_clm = [ 1              1              1                1                 2                 3                 1                 2  ];
fcfg.eff_clm = { [1]            [1]            [1]              [1]               [2]               [3]               [1]               [2] };
fcfg.ana_col = { 'blue'         'dark red'     'purple'         'red'             'purple'          'black'           'blue'            'purple'};
fcfg.ana_lbl = { 'VisualActive' 'VisualActive' 'AuditoryActive' 'VisualUni'       'VisualBi'        'TotalUni'        'AuditoryUni'     'AuditoryBi'};

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

% Get amp values 
for iT = 1:numel(fcfg.typ)

    cfg = [];
    cfg.clr_fld = clr_fld;
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.eff_clm = fcfg.eff_clm;
    cfg.sbj_nme = 'total';
    [amp_tab,amp_plt] = mmil_amp_tab(cfg);
    
    mmil_chk_dir(out_dir)
    tme_plt = amp_plt;
    save([out_dir '/' 'tme_pnt.mat'],'tme_plt')
    
    for iH = 1:numel(fcfg.hms)
        % Make plot
        cfg = [];
        cfg.typ = fcfg.typ{iT};
        cfg.hms = fcfg.hms{iH};
        cfg.ana_col = fcfg.ana_col;
        cfg.ana_lbl = fcfg.ana_lbl;
        cfg.out_dir = out_dir;
        cfg.reg_nme = fcfg.reg_nme;
        cfg.cmb_reg = fcfg.cmb_reg;
        cfg.exc_amp = [-100 100];
        cfg.ylm_typ = 'maxmin';
%         mmil_amp_plt(cfg,amp_plt{iH})
    end
    
end




























