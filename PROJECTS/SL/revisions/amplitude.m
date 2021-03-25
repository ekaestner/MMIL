clear; clc;

clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/SL/';

%% Overall Plot
out_dir = [clr_fld '/' 'manuscript' '/' 'figure6' '/' 'Modal_amplitude'];

% Setup
fcfg = [];

fcfg.sig_nme = { 'pap_anv_1500' 'pap_anv_1500'};
fcfg.sig_col = { 1              2 };

fcfg.typ     = { 'hgp'  };
fcfg.ele_typ = { 'ecog' };

fcfg.ana_nme = { 'vis_lng'      'aud_lng' };
fcfg.ana_clm = [ 1              1 ];
fcfg.eff_clm = { [1]            [1] };
fcfg.ana_col = { 'dark red'     'blue' };
fcfg.ana_lbl = { 'VisualActive' 'AuditoryActive' };

fcfg.hms     = {'lhs' 'rhs'};

fcfg.cmb_reg = { {'caudal-fusiform' 'middle-fusiform' } ...
                 {'caudal-STG' 'middle-STG' 'rostral-STG'} ...
                 {'supramarginal'} ...
                 {'inferior-precentral' 'middle-precentral'} ...
                 {'parsopercularis'} };
fcfg.reg_nme = {'Fusiform' ...
                'STG' ...
                'Supramarginal' ...
                'Precentral' ...
                'Pars Operculatris' };

% Get amp values 
for iT = 1:numel(fcfg.typ)

    cfg = [];
    cfg.clr_fld = clr_fld;
    cfg.sbj_nme = 'total';
    cfg.typ     = fcfg.typ{iT};
    cfg.ele_typ = fcfg.ele_typ{1};
    cfg.ana_nme = fcfg.ana_nme;
    cfg.ana_clm = fcfg.ana_clm;
    cfg.sig_nme = fcfg.sig_nme;
    cfg.sig_col = fcfg.sig_col;
    [amp_tab,amp_plt] = mmil_amp_tab_v2(cfg);
       
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