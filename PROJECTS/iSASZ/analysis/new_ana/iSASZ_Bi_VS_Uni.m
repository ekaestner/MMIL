clear; clc;

%% Responsive - Uni versus Bi - Proportion
clear; clc; 

reg.tot_reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.tot_reg.par_reg = { 'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tot_reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.tot_reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.tot_reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};
reg.tot_reg.ifg_reg = { 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

reg.cmb_reg.occ = { 'Lateral Occipital' };
reg.cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
reg.cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
reg.cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };
   
reg_nme = fieldnames(reg);

for iN = 1:numel(reg_nme)
    nme.(reg_nme{iN}) = fieldnames(reg.(reg_nme{iN}));
    
    reg.(reg_nme{iN}).tot_reg = cell(0);
    for iNM = 1:numel(nme.(reg_nme{iN}))
        reg.(reg_nme{iN}).tot_reg = [reg.(reg_nme{iN}).tot_reg reg.(reg_nme{iN}).(nme.(reg_nme{iN}){iNM})];
    end
    nme.(reg_nme{iN}){end+1} = 'tot_reg';
end

% Uni VS Bi WITHIN Proportion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_vis_bim_act' 'pap_vis_bim_act' 'pap_vis_bim_act' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3]             [7 8]             [12 13]};
fcfg.use_col = {[1 2] [3 2]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/' 'VisualBiModal'];
wth_lhs_act = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_act.cmb_reg.pvl_cnv(:,[1 3]) wth_lhs_act.cmb_reg.pvl(:,3)]

[wth_lhs_act.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_act.cmb_reg.pvl(:,2)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_aud_bim_act' 'pap_aud_bim_act' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8] };
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/' 'AuditoryBiModal'];
wth_lhs_act = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_act.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_act.cmb_reg.pvl(:,2)]

%% Responsive - Uni versus Bi - Timing
clc; clear

% Grainsize 
reg.tot_reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.tot_reg.par_reg = { 'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tot_reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.tot_reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.tot_reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};
reg.tot_reg.ifg_reg = { 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

reg.cmb_reg.occ = { 'Lateral Occipital' };
reg.cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
reg.cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
reg.cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };
   
reg_nme = fieldnames(reg);

for iN = 1:numel(reg_nme)
    nme.(reg_nme{iN}) = fieldnames(reg.(reg_nme{iN}));
    
    reg.(reg_nme{iN}).tot_reg = cell(0);
    for iNM = 1:numel(nme.(reg_nme{iN}))
        reg.(reg_nme{iN}).tot_reg = [reg.(reg_nme{iN}).tot_reg reg.(reg_nme{iN}).(nme.(reg_nme{iN}){iNM})];
    end
    nme.(reg_nme{iN}){end+1} = 'tot_reg';
end

% Bimodal Versus Uni-Modal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = { '' '' '' '' '' '' 'vis_uni' 'vis_bim' 'tot_uni' 'aud_uni' 'aud_bim' };
fcfg.use_col = { [7 8] [9 8] [10 11] };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/'];
wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_eff.cmb_reg.pvl

%% Repetition - Uni versus Bi - Proportion
% Uni VS Bi WITHIN Proportion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_vis_bim_rep' 'pap_vis_bim_rep' 'pap_vis_bim_rep' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3]             [7 8]             [12 13]};
fcfg.use_col = {[1 2] [3 2]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/' 'VisualBiModal'];
wth_lhs_rep = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 3]) wth_lhs_rep.cmb_reg.pvl(:,3)]

[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_rep.cmb_reg.pvl(:,2)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_aud_bim_rep' 'pap_aud_bim_rep' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8] };
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/' 'AuditoryBiModal'];
wth_lhs_rep = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_rep.cmb_reg.pvl(:,2)]

%% Repetition - Uni versus Bi - Timing
clc; clear

% Grainsize 
reg.tot_reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.tot_reg.par_reg = { 'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tot_reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.tot_reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.tot_reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};
reg.tot_reg.ifg_reg = { 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

reg.cmb_reg.occ = { 'Lateral Occipital' };
reg.cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
reg.cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
reg.cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };
   
reg_nme = fieldnames(reg);

for iN = 1:numel(reg_nme)
    nme.(reg_nme{iN}) = fieldnames(reg.(reg_nme{iN}));
    
    reg.(reg_nme{iN}).tot_reg = cell(0);
    for iNM = 1:numel(nme.(reg_nme{iN}))
        reg.(reg_nme{iN}).tot_reg = [reg.(reg_nme{iN}).tot_reg reg.(reg_nme{iN}).(nme.(reg_nme{iN}){iNM})];
    end
    nme.(reg_nme{iN}){end+1} = 'tot_reg';
end

% Bimodal Versus Uni-Modal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = { '' '' '' '' '' '' '' '' '' '' '' 'vis_uni' 'vis_bim' 'tot_uni' 'aud_uni' 'aud_bim' };
fcfg.use_col = { [12 13] [14 13] [15 16] };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/'];
wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_eff.cmb_reg.pvl

%% iN400 - Uni versus Bi - Proportion
clear; clc;

reg.tot_reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.tot_reg.par_reg = { 'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tot_reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.tot_reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.tot_reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};
reg.tot_reg.ifg_reg = { 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

reg.cmb_reg.occ = { 'Lateral Occipital' };
reg.cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
reg.cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
reg.cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };
   
reg_nme = fieldnames(reg);

for iN = 1:numel(reg_nme)
    nme.(reg_nme{iN}) = fieldnames(reg.(reg_nme{iN}));
    
    reg.(reg_nme{iN}).tot_reg = cell(0);
    for iNM = 1:numel(nme.(reg_nme{iN}))
        reg.(reg_nme{iN}).tot_reg = [reg.(reg_nme{iN}).tot_reg reg.(reg_nme{iN}).(nme.(reg_nme{iN}){iNM})];
    end
    nme.(reg_nme{iN}){end+1} = 'tot_reg';
end

% Within Region
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.eff_typ = 'lfp';
fcfg.loc_nme = {'pap_vis_bim_rep_nob' 'pap_vis_bim_rep_nob' 'pap_vis_bim_rep_nob' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3]             [7 8]             [12 13]};
fcfg.use_col = {[1 2] [3 2]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/' 'VisualBiModal'];
wth_lhs_rep_nob = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_rep_nob.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_rep_nob.cmb_reg.pvl(:,2)]
[wth_lhs_rep_nob.cmb_reg.pvl_cnv(:,[1 3]) wth_lhs_rep_nob.cmb_reg.pvl(:,3)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.eff_typ = 'lfp';
fcfg.loc_nme = {'pap_aud_bim_rep_nob' 'pap_aud_bim_rep_nob' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8] };
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/' 'AuditoryBiModal'];
wth_lhs_rep_nob = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_rep_nob.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_rep_nob.cmb_reg.pvl(:,2)]

%% iN400 - Uni versus Bi - Timing
clc; clear

% Grainsize 
reg.tot_reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.tot_reg.par_reg = { 'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tot_reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.tot_reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.tot_reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};
reg.tot_reg.ifg_reg = { 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

reg.cmb_reg.occ = { 'Lateral Occipital' };
reg.cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
reg.cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
reg.cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };
   
reg_nme = fieldnames(reg);

for iN = 1:numel(reg_nme)
    nme.(reg_nme{iN}) = fieldnames(reg.(reg_nme{iN}));
    
    reg.(reg_nme{iN}).tot_reg = cell(0);
    for iNM = 1:numel(nme.(reg_nme{iN}))
        reg.(reg_nme{iN}).tot_reg = [reg.(reg_nme{iN}).tot_reg reg.(reg_nme{iN}).(nme.(reg_nme{iN}){iNM})];
    end
    nme.(reg_nme{iN}){end+1} = 'tot_reg';
end

% Bimodal Versus Uni-Modal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = { '' '' '' 'vis_uni' 'vis_bim' 'tot_uni' 'aud_uni' 'aud_bim' };
fcfg.eff_typ = 'lfp';
fcfg.use_col = { [4 5] [6 5] [7 8] };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTableLFP' '/'];
wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_eff.cmb_reg.pvl

