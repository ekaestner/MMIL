clc; clear

%% Grainsize 
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

%% Distribution
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_vis_rep' 'pap_aud_rep' 'pap_bim_rep' };
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3]         [2 3]         [2 3]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
hms_sig_rep = mmil_fisher_exact_region_test(fcfg);

hms_sig_rep.cmb_reg.VisualRepetition
hms_sig_rep.cmb_reg.AuditoryRepetition
hms_sig_rep.cmb_reg.BimodalRepetition

% WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_vis_rep' 'pap_aud_rep' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3]         [2 3]         };
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
wth_lhs_rep = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_rep.cmb_reg.pvl(:,2)]
[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 3]) wth_lhs_rep.cmb_reg.pvl(:,3)]
[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 4]) wth_lhs_rep.cmb_reg.pvl(:,4)]

% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'wth_reg';
% fcfg.loc_nme = {'pap_vis_rep' 'pap_aud_rep' 'pap_bim_rep' };
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.col     = {[2 3]         [2 3]         [2 3]};
% fcfg.use_col = {[1 2] [1 3] [2 3]};
% fcfg.plt     = 1;
% fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
% wth_rhs_rep = mmil_fisher_exact_region_test(fcfg);
% 
% [wth_rhs_rep.cmb_reg.pvl_cnv(:,[1 2]) wth_rhs_rep.cmb_reg.pvl(:,2)]

% BETWEEN REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_vis_rep' 'pap_aud_rep' 'pap_bim_rep' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3]         [2 3]         [2 3]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
bet_lhs_rep = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_lhs_rep.cmb_reg.VisualRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.VisualRepetition.pvl{1,6}  bet_lhs_rep.cmb_reg.VisualRepetition.pvl(9,6) ; ...
            bet_lhs_rep.cmb_reg.VisualRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.VisualRepetition.pvl{1,7}  bet_lhs_rep.cmb_reg.VisualRepetition.pvl(9,7) ; ...
            bet_lhs_rep.cmb_reg.VisualRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.VisualRepetition.pvl{1,10} bet_lhs_rep.cmb_reg.VisualRepetition.pvl(9,10) ]
[bet_lhs_rep.cmb_reg.VisualRepetition.pvl([1 end],:) ; bet_lhs_rep.cmb_reg.VisualRepetition.pvl_cnv(end,:)]

spc_tst = [ bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl{1,6}  bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl(9,6) ; ...
            bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl{1,7}  bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl(9,7) ; ...
            bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl{1,10} bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl(9,10) ]
[bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl([1 end],:) ; bet_lhs_rep.cmb_reg.AuditoryRepetition.pvl_cnv(end,:)]

spc_tst = [ bet_lhs_rep.cmb_reg.BimodalRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.BimodalRepetition.pvl{1,6}  bet_lhs_rep.cmb_reg.BimodalRepetition.pvl(9,6) ; ...
            bet_lhs_rep.cmb_reg.BimodalRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.BimodalRepetition.pvl{1,7}  bet_lhs_rep.cmb_reg.BimodalRepetition.pvl(9,7) ; ...
            bet_lhs_rep.cmb_reg.BimodalRepetition.pvl{9,1} bet_lhs_rep.cmb_reg.BimodalRepetition.pvl{1,10} bet_lhs_rep.cmb_reg.BimodalRepetition.pvl(9,10) ]
[bet_lhs_rep.cmb_reg.BimodalRepetition.pvl([1 end],:) ; bet_lhs_rep.cmb_reg.BimodalRepetition.pvl_cnv(end,:)]

% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg';
% fcfg.loc_nme = {'pap_vis_rep' 'pap_aud_rep' 'pap_bim_rep' };
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.col     = {[2 3]         [2 3]         [2 3]};
% fcfg.plt     = 1;
% fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
% bet_rhs_rep = mmil_fisher_exact_region_test(fcfg);

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

[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 2]) wth_lhs_rep.cmb_reg.pvl(:,2)]
[wth_lhs_rep.cmb_reg.pvl_cnv(:,[1 3]) wth_lhs_rep.cmb_reg.pvl(:,3)]

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

% OVERLAP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pap_vis_rep = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_vis_rep/subjects/total/pap_vis_rep_plt');
pap_vis_rep = pap_vis_rep(2:end,:);

pap_aud_rep = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/sig_chn/hgp/ecog/split/pap_aud_rep/subjects/total/pap_aud_rep_plt');
pap_aud_rep = pap_aud_rep(2:end,:);

cfg = [];
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
cfg.sbj_nme = 'total';
cfg.loc_typ = 'split';
cfg.ele_typ = 'ecog';
cfg.ecg_hms =  {'lhs' 'rhs'};
tot_ele = mmil_location_calc(cfg);

cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [pap_vis_rep(2:end,2) pap_vis_rep(2:end,3) pap_aud_rep(2:end,3)];
cfg.eff_nme = {'Vis'                'Aud'                'Bim'};
vis_and_aud_tot = mmil_ovr_lap_ele_v3(cfg);

vis_and_aud_tot.lhs.Aud_in_Vis
vis_and_aud_tot.lhs.Aud_in_Vis_total

vis_and_aud_tot.lhs.Vis_in_Aud
vis_and_aud_tot.lhs.Vis_in_Aud_total

%% Timing
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

% HEMISPHERE
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = { 'vis_rep' 'aud_rep' 'bim_rep' };

fcfg.col     = [ 2          4         6 ];
cfg.ylm_typ  = { ''         ''        'maxmin' };

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/'];
hms_sig_eff  = mmil_ranksum_region_test(fcfg);

hms_sig_eff.cmb_reg

% WITHIN REGIONS
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = { '' 'vis_rep' '' 'aud_rep' '' 'bim_rep' };
fcfg.use_col = { [2 4] };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/'];
wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_eff.cmb_reg.pvl

% BETWEEN Regions
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = { 'vis_rep' 'aud_rep' 'bim_rep' };

fcfg.col     = [ 2          4         6 ];
fcfg.ylm_typ  = { ''   ''     ''  ''  'maxmin' 'maxmin' };

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/'];
bet_sig_lhs_eff = mmil_ranksum_region_test(fcfg);


spc_tst = [ bet_sig_lhs_eff.cmb_reg.vis_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.vis_rep.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.vis_rep.pvl(9,6) ; ...
            bet_sig_lhs_eff.cmb_reg.vis_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.vis_rep.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.vis_rep.pvl(9,8) ; ...
            bet_sig_lhs_eff.cmb_reg.vis_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.vis_rep.pvl{1,10} bet_sig_lhs_eff.cmb_reg.vis_rep.pvl(9,10) ]
          [bet_sig_lhs_eff.cmb_reg.vis_rep.pvl([1 end],:) ]

spc_tst = [ bet_sig_lhs_eff.cmb_reg.aud_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.aud_rep.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.aud_rep.pvl(9,6) ; ...
            bet_sig_lhs_eff.cmb_reg.aud_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.aud_rep.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.aud_rep.pvl(9,8) ; ...
            bet_sig_lhs_eff.cmb_reg.aud_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.aud_rep.pvl{1,10} bet_sig_lhs_eff.cmb_reg.aud_rep.pvl(9,10) ]
          [bet_sig_lhs_eff.cmb_reg.aud_rep.pvl([1 end],:) ]

% spc_tst = [ bet_sig_lhs_eff.cmb_reg.bim_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.bim_rep.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.bim_rep.pvl(9,6) ; ...
%             bet_sig_lhs_eff.cmb_reg.bim_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.bim_rep.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.bim_rep.pvl(9,8) ; ...
%             bet_sig_lhs_eff.cmb_reg.bim_rep.pvl{9,1} bet_sig_lhs_eff.cmb_reg.bim_rep.pvl{1,10} bet_sig_lhs_eff.cmb_reg.bim_rep.pvl(9,10) ]
%           [bet_sig_lhs_eff.cmb_reg.bim_rep.pvl([1 end],:) ]
%           
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/iSASZ/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg';
% fcfg.loc_nme = { 'vis_rep' 'aud_rep' 'bim_rep' };
% 
% fcfg.col     = [ 2          4         6 ];
% fcfg.ylm_typ  = { ''   ''     ''  ''  'maxmin' 'maxmin' };
% 
% fcfg.hms        = {'rhs' 'rhs'};
% fcfg.plt        = 1;
% fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'TimeTable' '/'];
% bet_sig_lhs_eff = mmil_ranksum_region_test(fcfg);





















