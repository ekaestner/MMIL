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

%% SPECIFIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specific Hemisphere
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'act'};

fcfg.col     = [1];

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'specific' '/'];
hms_sig_spc = mmil_ranksum_region_test(fcfg);

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'hem_lng';
% fcfg.loc_nme = {'act'       };
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_rsp_600' };
% fcfg.loc_col = {1             };
% fcfg.eff_typ = {'vis_stm_01'  };
% fcfg.loc_eve = {'!'           };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'specific' '/'];
% hms_sig_spc_lng = mmil_ranksum_region_test(fcfg);

% Specific Within Hemisphere
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'act'};
fcfg.col     = [1];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'specific' '/'];
bet_sig_lhs_spc = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_lhs_spc.cmb_reg.act.pvl{9,1} bet_sig_lhs_spc.cmb_reg.act.pvl{1,6} bet_sig_lhs_spc.cmb_reg.act.pvl(9,6) ; ...
            bet_sig_lhs_spc.cmb_reg.act.pvl{9,1} bet_sig_lhs_spc.cmb_reg.act.pvl{1,7} bet_sig_lhs_spc.cmb_reg.act.pvl(9,7) ; ...
            bet_sig_lhs_spc.cmb_reg.act.pvl{9,1} bet_sig_lhs_spc.cmb_reg.act.pvl{1,8} bet_sig_lhs_spc.cmb_reg.act.pvl(9,8) ]
bet_sig_lhs_spc.cmb_reg.act.pvl([1 end],:)

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'act'};
fcfg.col     = [1];

fcfg.hms        = {'rhs' 'rhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'specific' '/'];
bet_sig_rhs_spc = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_rhs_spc.cmb_reg.act.pvl{9,1} bet_sig_rhs_spc.cmb_reg.act.pvl{1,6} bet_sig_rhs_spc.cmb_reg.act.pvl(9,6) ; ...
            bet_sig_rhs_spc.cmb_reg.act.pvl{9,1} bet_sig_rhs_spc.cmb_reg.act.pvl{1,7} bet_sig_rhs_spc.cmb_reg.act.pvl(9,7) ; ...
            bet_sig_rhs_spc.cmb_reg.act.pvl{9,1} bet_sig_rhs_spc.cmb_reg.act.pvl{1,8} bet_sig_rhs_spc.cmb_reg.act.pvl(9,8) ]
bet_sig_rhs_spc.cmb_reg.act.pvl([1 end],:)

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg_lng';
% fcfg.loc_nme = {'act'};
% fcfg.col     = [1];
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_rsp_600' };
% fcfg.loc_col = {1             };
% fcfg.eff_typ = {'vis_stm_01'  };
% fcfg.loc_eve = {'!'           };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'lhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'specific' '/'];
% bet_sig_lhs_spc_lng = mmil_ranksum_region_test(fcfg);
% 
% 
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg_lng';
% fcfg.loc_nme = {'act'};
% fcfg.col     = [1];
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_rsp_600' };
% fcfg.loc_col = {1             };
% fcfg.eff_typ = {'vis_stm_01'  };
% fcfg.loc_eve = {'!'           };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'specific' '/'];
% bet_sig_rhs_spc_lng = mmil_ranksum_region_test(fcfg);

%% GRAINSIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAINSIZE COMPARING BETWEEN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'ltr' 'wrd' 'fls'};

fcfg.col     = [1 2 3];

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
hms_sig_grn = mmil_ranksum_region_test(fcfg);

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'hem_lng';
% fcfg.loc_nme = {'ltr'                'wrd'            'fls' };
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_wrd_600'        'pap_wrd_600'    'pap_con_600' };
% fcfg.loc_col = {[1 3]                3                1 };
% fcfg.eff_typ = {'vis_wrd_ffn_msk_01' 'vis_wrd_msk_01' 'vis_wrd_ffn_msk_01' };
% fcfg.loc_eve = {3                    3                6 };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
% hms_sig_lng_grn = mmil_ranksum_region_test(fcfg);

% GRAINSIZE COMPARING WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
fcfg.use_col = {[1 2] [3 2] [1 3]};

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
wth_sig_lhs_grn = mmil_ranksum_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
fcfg.use_col = {[1 2] [3 2] [1 3]};

fcfg.hms     = {'rhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
wth_sig_rhs_grn = mmil_ranksum_region_test(fcfg);

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'wth_reg_lng';
% fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
% fcfg.use_col = {[1 2] [3 2] [1 3]};
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_wrd_600'        'pap_wrd_600'    'pap_con_600' };
% fcfg.loc_col = {[1 3]                3                1 };
% fcfg.eff_typ = {'vis_wrd_ffn_msk_01' 'vis_wrd_msk_01' 'vis_wrd_ffn_msk_01' };
% fcfg.loc_eve = {3                    3                6 };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'lhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
% wth_sig_lhs_grn_lng = mmil_ranksum_region_test(fcfg);
% 
% 
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'wth_reg_lng';
% fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
% fcfg.use_col = {[1 2] [3 2] [1 3]};
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_wrd_600'        'pap_wrd_600'    'pap_con_600' };
% fcfg.loc_col = {[1 3]                3                1 };
% fcfg.eff_typ = {'vis_wrd_ffn_msk_01' 'vis_wrd_msk_01' 'vis_wrd_ffn_msk_01' };
% fcfg.loc_eve = {3                    3                6 };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
% wth_sig_rhs_grn_lng = mmil_ranksum_region_test(fcfg);

% GRAINSIZE COMPARING BETWEEN REGIONS
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
fcfg.col     = [1 2 3];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
bet_sig_lhs_grn = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_lhs_grn.cmb_reg.ltr.pvl{9,1} bet_sig_lhs_grn.cmb_reg.ltr.pvl{1,6} bet_sig_lhs_grn.cmb_reg.ltr.pvl(9,6) ; ...
            bet_sig_lhs_grn.cmb_reg.ltr.pvl{9,1} bet_sig_lhs_grn.cmb_reg.ltr.pvl{1,7} bet_sig_lhs_grn.cmb_reg.ltr.pvl(9,7) ; ...
            bet_sig_lhs_grn.cmb_reg.ltr.pvl{9,1} bet_sig_lhs_grn.cmb_reg.ltr.pvl{1,8} bet_sig_lhs_grn.cmb_reg.ltr.pvl(9,8) ]
bet_sig_lhs_grn.cmb_reg.ltr.pvl([1 end],:)

spc_tst = [ bet_sig_lhs_grn.cmb_reg.wrd.pvl{9,1} bet_sig_lhs_grn.cmb_reg.wrd.pvl{1,6} bet_sig_lhs_grn.cmb_reg.wrd.pvl(9,6) ; ...
            bet_sig_lhs_grn.cmb_reg.wrd.pvl{9,1} bet_sig_lhs_grn.cmb_reg.wrd.pvl{1,7} bet_sig_lhs_grn.cmb_reg.wrd.pvl(9,7) ; ...
            bet_sig_lhs_grn.cmb_reg.wrd.pvl{9,1} bet_sig_lhs_grn.cmb_reg.wrd.pvl{1,8} bet_sig_lhs_grn.cmb_reg.wrd.pvl(9,8) ]
bet_sig_lhs_grn.cmb_reg.wrd.pvl([1 end],:)

spc_tst = [ bet_sig_lhs_grn.cmb_reg.fls.pvl{9,1} bet_sig_lhs_grn.cmb_reg.fls.pvl{1,6} bet_sig_lhs_grn.cmb_reg.fls.pvl(9,6) ; ...
            bet_sig_lhs_grn.cmb_reg.fls.pvl{9,1} bet_sig_lhs_grn.cmb_reg.fls.pvl{1,7} bet_sig_lhs_grn.cmb_reg.fls.pvl(9,7) ; ...
            bet_sig_lhs_grn.cmb_reg.fls.pvl{9,1} bet_sig_lhs_grn.cmb_reg.fls.pvl{1,8} bet_sig_lhs_grn.cmb_reg.fls.pvl(9,8) ]
bet_sig_lhs_grn.cmb_reg.fls.pvl([1 end],:)

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
fcfg.col     = [1 2 3];

fcfg.hms        = {'rhs' 'rhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
bet_sig_rhs_grn = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_rhs_grn.cmb_reg.fls.pvl{9,1} bet_sig_rhs_grn.cmb_reg.fls.pvl{1,6} bet_sig_rhs_grn.cmb_reg.fls.pvl(9,6) ; ...
            bet_sig_rhs_grn.cmb_reg.fls.pvl{9,1} bet_sig_rhs_grn.cmb_reg.fls.pvl{1,7} bet_sig_rhs_grn.cmb_reg.fls.pvl(9,7) ; ...
            bet_sig_rhs_grn.cmb_reg.fls.pvl{9,1} bet_sig_rhs_grn.cmb_reg.fls.pvl{1,8} bet_sig_rhs_grn.cmb_reg.fls.pvl(9,8) ]
bet_sig_rhs_grn.cmb_reg.fls.pvl([1 end],:)

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg_lng';
% fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
% fcfg.col     = [1 2 3];
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_wrd_600'        'pap_wrd_600'    'pap_con_600' };
% fcfg.loc_col = {[1 3]                3                1 };
% fcfg.eff_typ = {'vis_wrd_ffn_msk_01' 'vis_wrd_msk_01' 'vis_wrd_ffn_msk_01' };
% fcfg.loc_eve = {3                    3                6 };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'lhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
% bet_sig_lhs_grn_lng = mmil_ranksum_region_test(fcfg);
% 
% 
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg_lng';
% fcfg.loc_nme = {'ltr' 'wrd' 'fls'};
% fcfg.col     = [1 2 3];
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_wrd_600'        'pap_wrd_600'    'pap_con_600' };
% fcfg.loc_col = {[1 3]                3                1 };
% fcfg.eff_typ = {'vis_wrd_ffn_msk_01' 'vis_wrd_msk_01' 'vis_wrd_ffn_msk_01' };
% fcfg.loc_eve = {3                    3                6 };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'];
% bet_sig_rhs_grn_lng = mmil_ranksum_region_test(fcfg);

%% REPETITION/FREQUENCY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPETITION/FREQUENCY COMPARING BETWEEN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'rep' 'frq'};
fcfg.col     = [1 2];

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
hms_sig_lex = mmil_ranksum_region_test(fcfg);

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'hem_lng';
% fcfg.loc_nme = {'rep'                'frq'            };
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_lex_600'        'pap_lex_600'    };
% fcfg.loc_col = {1                    3                };
% fcfg.eff_typ = {'vis_old'            'vis_wrd_frq'    };
% fcfg.loc_eve = {3                    121              };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
% hms_sig_lex_lng = mmil_ranksum_region_test(fcfg);

% REPETITION/FREQUENCY COMPARING WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'rep'                'frq'            };

fcfg.use_col = {[1 2]};

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
wth_sig_lhs_lex = mmil_ranksum_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'rep'                'frq'            };

fcfg.use_col = {[1 2]};

fcfg.hms     = {'rhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
wth_sig_rhs_lex = mmil_ranksum_region_test(fcfg);

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'wth_reg_lng';
% fcfg.loc_nme = {'rep'                'frq'            };
% fcfg.use_col = {[1 2]};
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_lex_600'        'pap_lex_600'    };
% fcfg.loc_col = {1                    3                };
% fcfg.eff_typ = {'vis_old'            'vis_wrd_frq'    };
% fcfg.loc_eve = {3                    121              };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'lhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
% wth_sig_lhs_lex_lng = mmil_ranksum_region_test(fcfg);
% 
% 
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'wth_reg_lng';
% fcfg.loc_nme = {'rep'                'frq'            };
% fcfg.use_col = {[1 2]};
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_lex_600'        'pap_lex_600'    };
% fcfg.loc_col = {1                    3                };
% fcfg.eff_typ = {'vis_old'            'vis_wrd_frq'    };
% fcfg.loc_eve = {3                    121              };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
% wth_sig_rhs_lex_lng = mmil_ranksum_region_test(fcfg);

% REPETITION/FREQUENCY COMPARING BETWEEN REGIONS
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'rep'                'frq'            };

fcfg.col     = [1 3];

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
bet_sig_lhs_lex = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_lhs_lex.cmb_reg.rep.pvl{9,1} bet_sig_lhs_lex.cmb_reg.rep.pvl{1,6} bet_sig_lhs_lex.cmb_reg.rep.pvl(9,6) ; ...
            bet_sig_lhs_lex.cmb_reg.rep.pvl{9,1} bet_sig_lhs_lex.cmb_reg.rep.pvl{1,7} bet_sig_lhs_lex.cmb_reg.rep.pvl(9,7) ; ...
            bet_sig_lhs_lex.cmb_reg.rep.pvl{9,1} bet_sig_lhs_lex.cmb_reg.rep.pvl{1,8} bet_sig_lhs_lex.cmb_reg.rep.pvl(9,8) ]
bet_sig_lhs_lex.cmb_reg.rep.pvl([1 end],:)

spc_tst = [ bet_sig_lhs_lex.cmb_reg.frq.pvl{9,1} bet_sig_lhs_lex.cmb_reg.frq.pvl{1,6} bet_sig_lhs_lex.cmb_reg.frq.pvl(9,6) ; ...
            bet_sig_lhs_lex.cmb_reg.frq.pvl{9,1} bet_sig_lhs_lex.cmb_reg.frq.pvl{1,7} bet_sig_lhs_lex.cmb_reg.frq.pvl(9,7) ; ...
            bet_sig_lhs_lex.cmb_reg.frq.pvl{9,1} bet_sig_lhs_lex.cmb_reg.frq.pvl{1,8} bet_sig_lhs_lex.cmb_reg.frq.pvl(9,8) ]
bet_sig_lhs_lex.cmb_reg.frq.pvl([1 end],:)

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'rep'                'frq'            };

fcfg.col     = [1 3];

fcfg.hms     = {'rhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
bet_sig_rhs_lex = mmil_ranksum_region_test(fcfg);

% effect length
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg_lng';
% fcfg.loc_nme = {'rep'                'frq'            };
% fcfg.col     = [1 2];
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_lex_600'        'pap_lex_600'    };
% fcfg.loc_col = {1                    3                };
% fcfg.eff_typ = {'vis_old'            'vis_wrd_frq'    };
% fcfg.loc_eve = {3                    121              };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'lhs' 'lhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
% bet_sig_lhs_lex_lng = mmil_ranksum_region_test(fcfg);
% 
% 
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg_lng';
% fcfg.loc_nme = {'rep'                'frq'            };
% fcfg.col     = [1 2];
% 
% fcfg.ele_typ = 'hgp';
% fcfg.loc_eff = {'pap_lex_600'        'pap_lex_600'    };
% fcfg.loc_col = {1                    3                };
% fcfg.eff_typ = {'vis_old'            'vis_wrd_frq'    };
% fcfg.loc_eve = {3                    121              };
% fcfg.tme_wdw = [0.100 0.750];
% 
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/'];
% bet_sig_rhs_lex_lng = mmil_ranksum_region_test(fcfg);

% Grainsize/Lexical Within Hemisphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'rep' 'frq' 'ltr' 'wrd' 'fls' };

fcfg.use_col = {[1 3]     [1 4]     [1 5]     [2 3]     [2 4]     [2 5]};

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = { [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/']   ...
                 [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'] };
wth_sig_lhs_lex_rep = mmil_ranksum_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'rep' 'frq' 'ltr' 'wrd' 'fls' };

fcfg.use_col = {[1 3] [1 4] [1 5] [2 3] [2 4] [2 5]};

fcfg.hms     = {'rhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = { [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'lexical' '/']   ...
                 [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'grainsize' '/'] };
wth_sig_rhs_lex_rep = mmil_ranksum_region_test(fcfg);

clear fcfg iN iNM nme reg reg_nme
save(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/manuscript/figure4/tme.mat'])



