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
reg.cmb_reg.itg = { 'Middle ITG' 'Rostral ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.pos = { 'Inferior Postcentral' 'Middle Postcentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
reg.cmb_reg.mid = { 'Middle Middle Frontal' 'caudalmiddlefrontal' };
   
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
% % Specific Hemisphere
% % onset
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'hem';
% fcfg.loc_nme = {'act'};
% 
% fcfg.col     = [1];
% 
% fcfg.hms     = {'lhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'specific' '/'];
% hms_sig_spc = mmil_ranksum_region_test(fcfg);
% 
% hms_sig_spc.cmb_reg.VisualSelective
% hms_sig_spc.cmb_reg.AuditorySelective
% 
% % Specific Within Hemisphere
% % onset
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg';
% fcfg.loc_nme = {'act'};
% fcfg.col     = [1];
% 
% fcfg.hms        = {'lhs' 'lhs'};
% fcfg.plt        = 1;
% fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'specific' '/'];
% bet_sig_lhs_spc = mmil_ranksum_region_test(fcfg);
% 
% spc_tst = [ bet_sig_lhs_spc.cmb_reg.act.pvl{9,1} bet_sig_lhs_spc.cmb_reg.act.pvl{1,6}  bet_sig_lhs_spc.cmb_reg.act.pvl(9,6) ; ...
%             bet_sig_lhs_spc.cmb_reg.act.pvl{9,1} bet_sig_lhs_spc.cmb_reg.act.pvl{1,7}  bet_sig_lhs_spc.cmb_reg.act.pvl(9,7) ; ...
%             bet_sig_lhs_spc.cmb_reg.act.pvl{9,1} bet_sig_lhs_spc.cmb_reg.act.pvl{1,10} bet_sig_lhs_spc.cmb_reg.act.pvl(9,10) ]
% [bet_sig_lhs_spc.cmb_reg.act.pvl([1 end],:)]
% 
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'bet_reg';
% fcfg.loc_nme = {'act'};
% fcfg.col     = [1];
% 
% fcfg.hms        = {'rhs' 'rhs'};
% fcfg.plt        = 1;
% fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'specific' '/'];
% bet_sig_rhs_spc = mmil_ranksum_region_test(fcfg);

%% LANGUAGE SELECTIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Text Hemisphere
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'txt' 'voc'};

fcfg.col     = [ 1 4 ];

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
hms_sig_eff  = mmil_ranksum_region_test(fcfg);

hms_sig_eff.cmb_reg

% Text WITHIN REGIONS
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'txt' 'fls' 'ltr' 'voc' 'phn' 'nse' 'mtc'};
fcfg.use_col = { [1 4] };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_eff.cmb_reg.pvl

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'txt' 'fls' 'ltr' 'voc' 'phn' 'nse' 'mtc'};
fcfg.use_col = { [1 4] };

fcfg.hms     = {'rhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
wth_sig_rhs_eff = mmil_ranksum_region_test(fcfg);

% Text BETWEEN Regions
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'txt' 'voc'};

fcfg.col     = [ 1 4 ];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_lhs_eff = mmil_ranksum_region_test(fcfg);


spc_tst = [ bet_sig_lhs_eff.cmb_reg.txt.pvl{9,1} bet_sig_lhs_eff.cmb_reg.txt.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.txt.pvl(9,6) ; ...
            bet_sig_lhs_eff.cmb_reg.txt.pvl{9,1} bet_sig_lhs_eff.cmb_reg.txt.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.txt.pvl(9,8) ; ...
            bet_sig_lhs_eff.cmb_reg.txt.pvl{9,1} bet_sig_lhs_eff.cmb_reg.txt.pvl{1,10} bet_sig_lhs_eff.cmb_reg.txt.pvl(9,10) ]
          [bet_sig_lhs_eff.cmb_reg.txt.pvl([1 end],:) ]

spc_tst = [ bet_sig_lhs_eff.cmb_reg.voc.pvl{9,1} bet_sig_lhs_eff.cmb_reg.voc.pvl{1,6}  bet_sig_lhs_eff.cmb_reg.voc.pvl(9,6) ; ...
            bet_sig_lhs_eff.cmb_reg.voc.pvl{9,1} bet_sig_lhs_eff.cmb_reg.voc.pvl{1,8}  bet_sig_lhs_eff.cmb_reg.voc.pvl(9,8) ; ...
            bet_sig_lhs_eff.cmb_reg.voc.pvl{9,1} bet_sig_lhs_eff.cmb_reg.voc.pvl{1,10} bet_sig_lhs_eff.cmb_reg.voc.pvl(9,10) ]
          [bet_sig_lhs_eff.cmb_reg.voc.pvl([1 end],:) ]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'txt' 'voc'};

fcfg.col     = [ 1 4 ];
fcfg.hms        = {'rhs' 'rhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_rhs_eff = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_rhs_eff.cmb_reg.txt.pvl{9,1} bet_sig_rhs_eff.cmb_reg.txt.pvl{1,6}  bet_sig_rhs_eff.cmb_reg.txt.pvl(9,6) ; ...
            bet_sig_rhs_eff.cmb_reg.txt.pvl{9,1} bet_sig_rhs_eff.cmb_reg.txt.pvl{1,8}  bet_sig_rhs_eff.cmb_reg.txt.pvl(9,8) ; ...
            bet_sig_rhs_eff.cmb_reg.txt.pvl{9,1} bet_sig_rhs_eff.cmb_reg.txt.pvl{1,10} bet_sig_rhs_eff.cmb_reg.txt.pvl(9,10) ]
          [bet_sig_rhs_eff.cmb_reg.txt.pvl([1 end],:) ]

spc_tst = [ bet_sig_rhs_eff.cmb_reg.voc.pvl{9,1} bet_sig_rhs_eff.cmb_reg.voc.pvl{1,6}  bet_sig_rhs_eff.cmb_reg.voc.pvl(9,6) ; ...
            bet_sig_rhs_eff.cmb_reg.voc.pvl{9,1} bet_sig_rhs_eff.cmb_reg.voc.pvl{1,8}  bet_sig_rhs_eff.cmb_reg.voc.pvl(9,8) ; ...
            bet_sig_rhs_eff.cmb_reg.voc.pvl{9,1} bet_sig_rhs_eff.cmb_reg.voc.pvl{1,10} bet_sig_rhs_eff.cmb_reg.voc.pvl(9,10) ]
          [bet_sig_rhs_eff.cmb_reg.voc.pvl([1 end],:) ]

%% PHONEME/LETTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Text Hemisphere
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = { 'ltr' 'phn' };

fcfg.col     = [1 2];

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'phoneme' '/'];
hms_sig_phn = mmil_ranksum_region_test(fcfg);

hms_sig_phn.cmb_reg

% Text BETWEEN Regions
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'ltr' 'phn'};
fcfg.col     = [1 2];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'phoneme' '/'];
bet_sig_lhs_phn = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_lhs_phn.cmb_reg.ltr.pvl{9,1} bet_sig_lhs_phn.cmb_reg.ltr.pvl{1,6}  bet_sig_lhs_phn.cmb_reg.ltr.pvl(9,6) ; ...
            bet_sig_lhs_phn.cmb_reg.ltr.pvl{9,1} bet_sig_lhs_phn.cmb_reg.ltr.pvl{1,8}  bet_sig_lhs_phn.cmb_reg.ltr.pvl(9,8) ; ...
            bet_sig_lhs_phn.cmb_reg.ltr.pvl{9,1} bet_sig_lhs_phn.cmb_reg.ltr.pvl{1,10} bet_sig_lhs_phn.cmb_reg.ltr.pvl(9,10) ]
          [bet_sig_lhs_phn.cmb_reg.ltr.pvl([1 end],:) ]

spc_tst = [ bet_sig_lhs_phn.cmb_reg.phn.pvl{9,1} bet_sig_lhs_phn.cmb_reg.phn.pvl{1,6}  bet_sig_lhs_phn.cmb_reg.phn.pvl(9,6) ; ...
            bet_sig_lhs_phn.cmb_reg.phn.pvl{9,1} bet_sig_lhs_phn.cmb_reg.phn.pvl{1,8}  bet_sig_lhs_phn.cmb_reg.phn.pvl(9,8) ; ...
            bet_sig_lhs_phn.cmb_reg.phn.pvl{9,1} bet_sig_lhs_phn.cmb_reg.phn.pvl{1,10} bet_sig_lhs_phn.cmb_reg.phn.pvl(9,10) ]
          [bet_sig_lhs_phn.cmb_reg.phn.pvl([1 end],:) ]
          
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'ltr' 'phn'};
fcfg.col     = [1 2];

fcfg.hms        = {'rhs' 'rhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'phoneme' '/'];
bet_sig_rhs_phn = mmil_ranksum_region_test(fcfg);

% Text BETWEEN CONDITIONS
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'ltr' 'phn'};
fcfg.use_col = {[1 2]};

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'phoneme' '/'];
wth_sig_lhs_phn = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_phn.cmb_reg.pvl

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'vis' 'fls' 'aud' 'nse' 'mtc'};
fcfg.use_col = {[1 2]};

fcfg.hms     = {'rhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'phoneme' '/'];
wth_sig_rhs_phn = mmil_ranksum_region_test(fcfg);

%% Control Selective
% Hemisphere
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'fls' 'nse'};

fcfg.col     = [ 3 6 ];

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
hms_sig_nse  = mmil_ranksum_region_test(fcfg);

hms_sig_nse.cmb_reg

% WITHIN REGIONS
% onset
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'wth_reg';
% fcfg.loc_nme = {'txt' 'fls' 'ltr' 'voc' 'phn' 'nse' 'mtc'};
% fcfg.use_col = { [3 6] };
% 
% fcfg.hms     = {'lhs' 'lhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
% wth_sig_lhs_eff = mmil_ranksum_region_test(fcfg);
% 
% wth_sig_lhs_eff.cmb_reg.pvl
% 
% fcfg = [];
% fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
% fcfg.nme     = nme;
% fcfg.reg     = reg;
% fcfg.typ     = 'wth_reg';
% fcfg.loc_nme = {'txt' 'fls' 'ltr' 'voc' 'phn' 'nse' 'mtc'};
% fcfg.use_col = { [3 6] };
% 
% fcfg.hms     = {'rhs' 'rhs'};
% fcfg.plt     = 1;
% fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
% wth_sig_rhs_nse = mmil_ranksum_region_test(fcfg);

% Text BETWEEN Regions
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'fls' 'nse'};

fcfg.col     = [ 3 6 ];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_lhs_nse = mmil_ranksum_region_test(fcfg);


spc_tst = [ bet_sig_lhs_nse.cmb_reg.fls.pvl{8,1} bet_sig_lhs_nse.cmb_reg.fls.pvl{1,6}  bet_sig_lhs_nse.cmb_reg.fls.pvl(8,6) ; ...
            bet_sig_lhs_nse.cmb_reg.fls.pvl{8,1} bet_sig_lhs_nse.cmb_reg.fls.pvl{1,7}  bet_sig_lhs_nse.cmb_reg.fls.pvl(8,7) ; ...
            bet_sig_lhs_nse.cmb_reg.fls.pvl{8,1} bet_sig_lhs_nse.cmb_reg.fls.pvl{1,10} bet_sig_lhs_nse.cmb_reg.fls.pvl(8,10) ]
          [bet_sig_lhs_nse.cmb_reg.fls.pvl([1 end],:) ]

spc_tst = [ bet_sig_lhs_nse.cmb_reg.nse.pvl{8,1} bet_sig_lhs_nse.cmb_reg.nse.pvl{1,6}  bet_sig_lhs_nse.cmb_reg.nse.pvl(8,6) ; ...
            bet_sig_lhs_nse.cmb_reg.nse.pvl{8,1} bet_sig_lhs_nse.cmb_reg.nse.pvl{1,7}  bet_sig_lhs_nse.cmb_reg.nse.pvl(8,7) ; ...
            bet_sig_lhs_nse.cmb_reg.nse.pvl{8,1} bet_sig_lhs_nse.cmb_reg.nse.pvl{1,10} bet_sig_lhs_nse.cmb_reg.nse.pvl(8,10) ]
          [bet_sig_lhs_nse.cmb_reg.nse.pvl([1 end],:) ]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'fls' 'nse'};

fcfg.col     = [ 3 6 ];
fcfg.hms        = {'rhs' 'rhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_rhs_nse = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,6}  bet_sig_rhs_nse.cmb_reg.txt.pvl(9,6) ; ...
            bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,8}  bet_sig_rhs_nse.cmb_reg.txt.pvl(9,8) ; ...
            bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,10} bet_sig_rhs_nse.cmb_reg.txt.pvl(9,10) ]
          [bet_sig_rhs_nse.cmb_reg.txt.pvl([1 end],:) ]

spc_tst = [ bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,6}  bet_sig_rhs_nse.cmb_reg.voc.pvl(9,6) ; ...
            bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,8}  bet_sig_rhs_nse.cmb_reg.voc.pvl(9,8) ; ...
            bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,10} bet_sig_rhs_nse.cmb_reg.voc.pvl(9,10) ]
            [bet_sig_rhs_nse.cmb_reg.voc.pvl([1 end],:) ]

%% MisMatch
% Hemisphere
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'mtc_sns'};

fcfg.col     = [ 7 ];

fcfg.hms     = {'lhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
hms_sig_phn  = mmil_ranksum_region_test(fcfg);

hms_sig_phn.cmb_reg

% BETWEEN Regions
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'mtc_sns'};

fcfg.col     = [ 7 ];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_lhs_phn = mmil_ranksum_region_test(fcfg);


spc_tst = [ bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl{8,1} bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl{1,6}  bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl(8,6) ; ...
            bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl{8,1} bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl{1,7}  bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl(8,7) ; ...
            bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl{8,1} bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl{1,10} bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl(8,10) ]
          [bet_sig_lhs_phn.cmb_reg.mtc_sns.pvl([1 end],:) ]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'fls' 'nse'};

fcfg.col     = [ 3 6 ];
fcfg.hms        = {'rhs' 'rhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_rhs_phn = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,6}  bet_sig_rhs_nse.cmb_reg.txt.pvl(9,6) ; ...
            bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,8}  bet_sig_rhs_nse.cmb_reg.txt.pvl(9,8) ; ...
            bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,10} bet_sig_rhs_nse.cmb_reg.txt.pvl(9,10) ]
          [bet_sig_rhs_nse.cmb_reg.txt.pvl([1 end],:) ]

spc_tst = [ bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,6}  bet_sig_rhs_nse.cmb_reg.voc.pvl(9,6) ; ...
            bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,8}  bet_sig_rhs_nse.cmb_reg.voc.pvl(9,8) ; ...
            bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,10} bet_sig_rhs_nse.cmb_reg.voc.pvl(9,10) ]
            [bet_sig_rhs_nse.cmb_reg.voc.pvl([1 end],:) ]
            

fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'txt' 'ltr' 'fls' 'voc' 'phn' 'nse' 'mtc'};
fcfg.use_col = { [2 5] };

fcfg.hms     = {'lhs' 'lhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
wth_sig_lhs_phn = mmil_ranksum_region_test(fcfg);

wth_sig_lhs_phn.cmb_reg.pvl

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'txt' 'ltr' 'fls' 'voc' 'phn' 'nse' 'mtc'};
fcfg.use_col = { [3 6] };

fcfg.hms     = {'rhs' 'rhs'};
fcfg.plt     = 1;
fcfg.out     = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
wth_sig_rhs_phn = mmil_ranksum_region_test(fcfg);

% Text BETWEEN Regions
% onset
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'ltr' 'phn'};

fcfg.col     = [ 2 5 ];

fcfg.hms        = {'lhs' 'lhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_lhs_phn = mmil_ranksum_region_test(fcfg);


spc_tst = [ bet_sig_lhs_phn.cmb_reg.ltr.pvl{8,1} bet_sig_lhs_phn.cmb_reg.ltr.pvl{1,6}  bet_sig_lhs_phn.cmb_reg.ltr.pvl(8,6) ; ...
            bet_sig_lhs_phn.cmb_reg.ltr.pvl{8,1} bet_sig_lhs_phn.cmb_reg.ltr.pvl{1,7}  bet_sig_lhs_phn.cmb_reg.ltr.pvl(8,7) ; ...
            bet_sig_lhs_phn.cmb_reg.ltr.pvl{8,1} bet_sig_lhs_phn.cmb_reg.ltr.pvl{1,10} bet_sig_lhs_phn.cmb_reg.ltr.pvl(8,10) ]
          [bet_sig_lhs_phn.cmb_reg.ltr.pvl([1 end],:) ]

spc_tst = [ bet_sig_lhs_phn.cmb_reg.phn.pvl{8,1} bet_sig_lhs_phn.cmb_reg.phn.pvl{1,6}  bet_sig_lhs_phn.cmb_reg.phn.pvl(8,6) ; ...
            bet_sig_lhs_phn.cmb_reg.phn.pvl{8,1} bet_sig_lhs_phn.cmb_reg.phn.pvl{1,7}  bet_sig_lhs_phn.cmb_reg.phn.pvl(8,7) ; ...
            bet_sig_lhs_phn.cmb_reg.phn.pvl{8,1} bet_sig_lhs_phn.cmb_reg.phn.pvl{1,10} bet_sig_lhs_phn.cmb_reg.phn.pvl(8,10) ]
          [bet_sig_lhs_phn.cmb_reg.phn.pvl([1 end],:) ]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'fls' 'nse'};

fcfg.col     = [ 3 6 ];
fcfg.hms        = {'rhs' 'rhs'};
fcfg.plt        = 1;
fcfg.out        = [fcfg.clr_fld '/' 'manuscript' '/' 'figure6' '/' 'AllEffect' '/'];
bet_sig_rhs_phn = mmil_ranksum_region_test(fcfg);

spc_tst = [ bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,6}  bet_sig_rhs_nse.cmb_reg.txt.pvl(9,6) ; ...
            bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,8}  bet_sig_rhs_nse.cmb_reg.txt.pvl(9,8) ; ...
            bet_sig_rhs_nse.cmb_reg.txt.pvl{9,1} bet_sig_rhs_nse.cmb_reg.txt.pvl{1,10} bet_sig_rhs_nse.cmb_reg.txt.pvl(9,10) ]
          [bet_sig_rhs_nse.cmb_reg.txt.pvl([1 end],:) ]

spc_tst = [ bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,6}  bet_sig_rhs_nse.cmb_reg.voc.pvl(9,6) ; ...
            bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,8}  bet_sig_rhs_nse.cmb_reg.voc.pvl(9,8) ; ...
            bet_sig_rhs_nse.cmb_reg.voc.pvl{9,1} bet_sig_rhs_nse.cmb_reg.voc.pvl{1,10} bet_sig_rhs_nse.cmb_reg.voc.pvl(9,10) ]
            [bet_sig_rhs_nse.cmb_reg.voc.pvl([1 end],:) ]