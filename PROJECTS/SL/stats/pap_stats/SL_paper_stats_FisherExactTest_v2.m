clc; clear;

%% FisherExactTest for Language-Sensitive Electrodes
% SETUP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg.reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.reg.par_reg = { 'Supramarginal' };
reg.reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

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

typ_nme = fieldnames(reg);

nme.(typ_nme{1}) = fieldnames(reg.(typ_nme{1}));
nme.(typ_nme{2}) = fieldnames(reg.(typ_nme{2}));

reg.(typ_nme{1}).tot_reg = cell(0);
for iN = 1:numel(nme.(typ_nme{1}))
    reg.(typ_nme{1}).tot_reg = [reg.(typ_nme{1}).tot_reg reg.(typ_nme{1}).(nme.(typ_nme{1}){iN})];
end
nme.(typ_nme{1}){end+1} = 'tot_reg';

reg.(typ_nme{2}).tot_reg = cell(0);
for iN = 1:numel(nme.(typ_nme{2}))
    reg.(typ_nme{2}).tot_reg = [reg.cmb_reg.tot_reg reg.(typ_nme{2}).(nme.(typ_nme{2}){iN})];
end
nme.(typ_nme{2}){end+1} = 'tot_reg';

%% SPECIFIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_anv_1500' 'pap_anv_1500' };
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
hms_sig_spc = mmil_fisher_exact_region_test(fcfg);

hms_sig_spc.cmb_reg.VisualSelective
hms_sig_spc.cmb_reg.AuditorySelective

% WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_anv_1500' 'pap_anv_1500' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
wth_lhs_spc = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_spc.cmb_reg.pvl_cnv wth_lhs_spc.cmb_reg.pvl(:,2)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_anv_1500' 'pap_anv_1500' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
wth_rhs_spc = mmil_fisher_exact_region_test(fcfg);

[wth_rhs_spc.cmb_reg.pvl_cnv wth_rhs_spc.cmb_reg.pvl(:,2)]

% BETWEEN REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_anv_1500' 'pap_anv_1500' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
bet_lhs_spc = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_lhs_spc.cmb_reg.VisualSelective.pvl{8,1} bet_lhs_spc.cmb_reg.VisualSelective.pvl{1,6}  bet_lhs_spc.cmb_reg.VisualSelective.pvl(8,6) ; ...
            bet_lhs_spc.cmb_reg.VisualSelective.pvl{8,1} bet_lhs_spc.cmb_reg.VisualSelective.pvl{1,7}  bet_lhs_spc.cmb_reg.VisualSelective.pvl(8,7) ; ...
            bet_lhs_spc.cmb_reg.VisualSelective.pvl{8,1} bet_lhs_spc.cmb_reg.VisualSelective.pvl{1,10} bet_lhs_spc.cmb_reg.VisualSelective.pvl(8,10) ]
[bet_lhs_spc.cmb_reg.VisualSelective.pvl([1 end],:) ; bet_lhs_spc.cmb_reg.VisualSelective.pvl_cnv(end,:)]

spc_tst = [ bet_lhs_spc.cmb_reg.AuditorySelective.pvl{8,1} bet_lhs_spc.cmb_reg.AuditorySelective.pvl{1,6}  bet_lhs_spc.cmb_reg.AuditorySelective.pvl(8,6) ; ...
            bet_lhs_spc.cmb_reg.AuditorySelective.pvl{8,1} bet_lhs_spc.cmb_reg.AuditorySelective.pvl{1,7}  bet_lhs_spc.cmb_reg.AuditorySelective.pvl(8,7) ; ...
            bet_lhs_spc.cmb_reg.AuditorySelective.pvl{8,1} bet_lhs_spc.cmb_reg.AuditorySelective.pvl{1,10} bet_lhs_spc.cmb_reg.AuditorySelective.pvl(8,10) ]
[bet_lhs_spc.cmb_reg.AuditorySelective.pvl([1 end],:) ; bet_lhs_spc.cmb_reg.AuditorySelective.pvl_cnv(end,:)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_anv_1500' 'pap_anv_1500' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
bet_rhs_spc = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_rhs_spc.cmb_reg.VisualSelective.pvl{8,1} bet_rhs_spc.cmb_reg.VisualSelective.pvl{1,6}  bet_rhs_spc.cmb_reg.VisualSelective.pvl(8,6) ; ...
            bet_rhs_spc.cmb_reg.VisualSelective.pvl{8,1} bet_rhs_spc.cmb_reg.VisualSelective.pvl{1,7}  bet_rhs_spc.cmb_reg.VisualSelective.pvl(8,7) ; ...
            bet_rhs_spc.cmb_reg.VisualSelective.pvl{8,1} bet_rhs_spc.cmb_reg.VisualSelective.pvl{1,10} bet_rhs_spc.cmb_reg.VisualSelective.pvl(8,10) ]
[bet_rhs_spc.cmb_reg.VisualSelective.pvl([1 end],:) ; bet_rhs_spc.cmb_reg.VisualSelective.pvl_cnv(end,:)]

spc_tst = [ bet_rhs_spc.cmb_reg.AuditorySelective.pvl{8,1} bet_rhs_spc.cmb_reg.AuditorySelective.pvl{1,6}  bet_rhs_spc.cmb_reg.AuditorySelective.pvl(8,6) ; ...
            bet_rhs_spc.cmb_reg.AuditorySelective.pvl{8,1} bet_rhs_spc.cmb_reg.AuditorySelective.pvl{1,7}  bet_rhs_spc.cmb_reg.AuditorySelective.pvl(8,7) ; ...
            bet_rhs_spc.cmb_reg.AuditorySelective.pvl{8,1} bet_rhs_spc.cmb_reg.AuditorySelective.pvl{1,10} bet_rhs_spc.cmb_reg.AuditorySelective.pvl(8,10) ]
[bet_rhs_spc.cmb_reg.AuditorySelective.pvl([1 end],:) ; bet_rhs_spc.cmb_reg.AuditorySelective.pvl_cnv(end,:)]

%% LANGUAGE-SPECIFIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_lng_950' 'pap_lng_950' };
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
hms_sig_lng = mmil_fisher_exact_region_test(fcfg);

hms_sig_lng.cmb_reg.VisualLanguage
hms_sig_lng.cmb_reg.AuditoryLanguage

% WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_lng_950' 'pap_lng_950' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
wth_lhs_lng = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_lng.cmb_reg.pvl_cnv wth_lhs_lng.cmb_reg.pvl(:,2)]
[wth_lhs_lng.cmb_reg.pvl_cnv([1 6 7 8 10],:) wth_lhs_lng.cmb_reg.pvl([1 6 7 8 10],2)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_lng_950' 'pap_lng_950' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
wth_rhs_lng = mmil_fisher_exact_region_test(fcfg);

[wth_rhs_lng.cmb_reg.pvl_cnv wth_rhs_lng.cmb_reg.pvl(:,2)]

% BETWEEN REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_lng_950' 'pap_lng_950' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
bet_lhs_lng = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_lhs_lng.cmb_reg.VisualLanguage.pvl{8,1} bet_lhs_lng.cmb_reg.VisualLanguage.pvl{1,6}  bet_lhs_lng.cmb_reg.VisualLanguage.pvl(8,6) ; ...
            bet_lhs_lng.cmb_reg.VisualLanguage.pvl{8,1} bet_lhs_lng.cmb_reg.VisualLanguage.pvl{1,7}  bet_lhs_lng.cmb_reg.VisualLanguage.pvl(8,7) ; ...
            bet_lhs_lng.cmb_reg.VisualLanguage.pvl{8,1} bet_lhs_lng.cmb_reg.VisualLanguage.pvl{1,10} bet_lhs_lng.cmb_reg.VisualLanguage.pvl(8,10) ]
[bet_lhs_lng.cmb_reg.VisualLanguage.pvl([1 end],:) ; bet_lhs_lng.cmb_reg.VisualLanguage.pvl_cnv(end,:)]

spc_tst = [ bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl{8,1} bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl{1,6}  bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl(8,6) ; ...
            bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl{8,1} bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl{1,7}  bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl(8,7) ; ...
            bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl{8,1} bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl{1,10} bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl(8,10) ]
[bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl([1 end],:) ; bet_lhs_lng.cmb_reg.AuditoryLanguage.pvl_cnv(end,:)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_lng_950' 'pap_lng_950' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
bet_rhs_lng = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_rhs_lng.cmb_reg.VisualSelective.pvl{8,1} bet_rhs_lng.cmb_reg.VisualSelective.pvl{1,6}  bet_rhs_lng.cmb_reg.VisualSelective.pvl(8,6) ; ...
            bet_rhs_lng.cmb_reg.VisualSelective.pvl{8,1} bet_rhs_lng.cmb_reg.VisualSelective.pvl{1,7}  bet_rhs_lng.cmb_reg.VisualSelective.pvl(8,7) ; ...
            bet_rhs_lng.cmb_reg.VisualSelective.pvl{8,1} bet_rhs_lng.cmb_reg.VisualSelective.pvl{1,10} bet_rhs_lng.cmb_reg.VisualSelective.pvl(8,10) ]
[bet_rhs_lng.cmb_reg.VisualSelective.pvl([1 end],:) ; bet_rhs_lng.cmb_reg.VisualSelective.pvl_cnv(end,:)]

spc_tst = [ bet_rhs_lng.cmb_reg.AuditorySelective.pvl{8,1} bet_rhs_lng.cmb_reg.AuditorySelective.pvl{1,6}  bet_rhs_lng.cmb_reg.AuditorySelective.pvl(8,6) ; ...
            bet_rhs_lng.cmb_reg.AuditorySelective.pvl{8,1} bet_rhs_lng.cmb_reg.AuditorySelective.pvl{1,7}  bet_rhs_lng.cmb_reg.AuditorySelective.pvl(8,7) ; ...
            bet_rhs_lng.cmb_reg.AuditorySelective.pvl{8,1} bet_rhs_lng.cmb_reg.AuditorySelective.pvl{1,10} bet_rhs_lng.cmb_reg.AuditorySelective.pvl(8,10) ]
[bet_rhs_lng.cmb_reg.AuditorySelective.pvl([1 end],:) ; bet_rhs_lng.cmb_reg.AuditorySelective.pvl_cnv(end,:)]

% BI-MODAL OVERLAP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lng_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_lng_950/subjects/total/pap_lng_950_plt');
lng_ovr_lap = lng_ovr_lap(2:end,:);

vis_ele = sum(cell2mat(lng_ovr_lap(:,3))); 
aud_ele = sum(cell2mat(lng_ovr_lap(:,4)));
bim_ele = sum(cell2mat(lng_ovr_lap(:,3)) & cell2mat(lng_ovr_lap(:,4)));

cfg = [];
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
cfg.sbj_nme = 'total';
cfg.loc_typ = 'split';
cfg.ele_typ = 'ecog';
cfg.ecg_hms =  {'lhs' 'rhs'};
tot_ele = mmil_location_calc(cfg);

% for iH = 1:numel(tot_ele)
%     
%     [~,vis_num] = intersect(tot_ele{iH}(:,1),lng_ovr_lap(logical(cell2mat(lng_ovr_lap(:,3))),2)); vis_num = tot_ele{iH}(vis_num,2);
%     [~,aud_num] = intersect(tot_ele{iH}(:,1),lng_ovr_lap(logical(cell2mat(lng_ovr_lap(:,4))),2)); aud_num = tot_ele{iH}(aud_num,2);
%     [~,bim_num] = intersect(tot_ele{iH}(:,1),lng_ovr_lap(cell2mat(lng_ovr_lap(:,3)) & cell2mat(lng_ovr_lap(:,4)),2)); bim_num = tot_ele{iH}(bim_num,2);
%     
%     vis_num = tabulate(vis_num); vis_num = vis_num(:,1:2); vis_num = mmil_order_table(vis_num);
%     aud_num = tabulate(aud_num); aud_num = aud_num(:,1:2); aud_num = mmil_order_table(aud_num);
%     bim_num = tabulate(bim_num); bim_num = bim_num(:,1:2); bim_num = mmil_order_table(bim_num);
%         
%     for iR = 1:numel(nme.cmb_reg)
%         
%         hld_dta{iR,1} = nme.cmb_reg{iR};
%         hld_dta{iR,2} = sum(cell2mat(vis_num(ismember(vis_num(:,1),reg.cmb_reg.(nme.cmb_reg{iR})),2)));
%         hld_dta{iR,3} = sum(cell2mat(aud_num(ismember(aud_num(:,1),reg.cmb_reg.(nme.cmb_reg{iR})),2)));
%         hld_dta{iR,4} = sum(cell2mat(bim_num(ismember(bim_num(:,1),reg.cmb_reg.(nme.cmb_reg{iR})),2)));
%         
%         hld_dta{iR,6} = myBinomTest(hld_dta{iR,4},hld_dta{iR,2},0.5,'lesser');                
%         hld_dta{iR,9}  = myBinomTest(hld_dta{iR,4},hld_dta{iR,3},0.5,'lesser');
%         hld_dta{iR,12} = myBinomTest(hld_dta{iR,4},hld_dta{iR,3}+hld_dta{iR,2}-hld_dta{iR,4},0.5,'lesser');
%         
%         if hld_dta{iR,4} > 0
%             hld_dta{iR,7} = myBinomTest(hld_dta{iR,4},hld_dta{iR,2},0.5,'greater');
%             hld_dta{iR,10} = myBinomTest(hld_dta{iR,4},hld_dta{iR,3},0.5,'greater');
%             hld_dta{iR,13} = myBinomTest(hld_dta{iR,4},hld_dta{iR,3}+hld_dta{iR,2}-hld_dta{iR,4},0.5,'greater');
%         else
%             hld_dta{iR,7} = 1;
%             hld_dta{iR,10} = 1;
%             hld_dta{iR,13} = 1;
%         end
%         
%     end
%     
%     dta.(cfg.ecg_hms{iH}) = hld_dta;
%     clear hld_dta
%     
% end

% 
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [lng_ovr_lap(2:end,2) lng_ovr_lap(2:end,3) lng_ovr_lap(2:end,4)];
cfg.eff_nme = {'Vis'                'Aud'                'Bim'};
vis_and_aud_tot = mmil_ovr_lap_ele_v3(cfg);

vis_and_aud_tot.lhs.Aud_in_Vis
vis_and_aud_tot.lhs.Aud_in_Vis_total

vis_and_aud_tot.lhs.Vis_in_Aud
vis_and_aud_tot.lhs.Vis_in_Aud_total

%% PHONEME/LETTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_phn_950' 'pap_phn_950' };
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] [7 8] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/' 'Fisher' '/'];
hms_sig_phn = mmil_fisher_exact_region_test(fcfg);

hms_sig_phn.cmb_reg.VisualLetter
hms_sig_phn.cmb_reg.AuditoryPhoneme

% WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_phn_950' 'pap_phn_950' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/' 'Fisher' '/'];
wth_lhs_phn = mmil_fisher_exact_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_phn_950' 'pap_phn_950'};
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/' 'Fisher' '/'];
wth_rhs_phn = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_phn.cmb_reg.pvl_cnv wth_lhs_phn.cmb_reg.pvl(:,2)]
[wth_rhs_phn.cmb_reg.pvl_cnv wth_rhs_phn.cmb_reg.pvl(:,2)]

% BETWEEN REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_phn_950' 'pap_phn_950' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3_5' '/' 'Fisher' '/'];
bet_lhs_phn = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_lhs_phn.cmb_reg.VisualLetter.pvl{8,1} bet_lhs_phn.cmb_reg.VisualLetter.pvl{1,6}  bet_lhs_phn.cmb_reg.VisualLetter.pvl(8,6) ; ...
            bet_lhs_phn.cmb_reg.VisualLetter.pvl{8,1} bet_lhs_phn.cmb_reg.VisualLetter.pvl{1,7}  bet_lhs_phn.cmb_reg.VisualLetter.pvl(8,7) ; ...
            bet_lhs_phn.cmb_reg.VisualLetter.pvl{8,1} bet_lhs_phn.cmb_reg.VisualLetter.pvl{1,10} bet_lhs_phn.cmb_reg.VisualLetter.pvl(8,10) ]
[bet_lhs_phn.cmb_reg.VisualLetter.pvl([1 end],:) ; bet_lhs_phn.cmb_reg.VisualLetter.pvl_cnv(end,:)]

spc_tst = [ bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl{8,1} bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl{1,6}  bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl(8,6) ; ...
            bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl{8,1} bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl{1,7}  bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl(8,7) ; ...
            bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl{8,1} bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl{1,10} bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl(8,10) ]
[bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl([1 end],:) ; bet_lhs_phn.cmb_reg.AuditoryPhoneme.pvl_cnv(end,:)]

%% MATCH/MISMATCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_mtc_1450' 'pap_mtc_1450' 'pap_mtc_1450'};
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] [7 8] [12 13]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'Fisher' '/'];
hms_sig_mtc = mmil_fisher_exact_region_test(fcfg);

hms_sig_mtc.cmb_reg.Mismatch_sensory
hms_sig_mtc.cmb_reg.Mismatch_total

% BETWEEN REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_mtc_1450' 'pap_mtc_1450' 'pap_mtc_1450'};
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8] [12 13]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'Fisher' '/'];
bet_lhs_mtc = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl{8,1} bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl{1,6}  bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl(8,6) ; ...
            bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl{8,1} bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl{1,7}  bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl(8,7) ; ...
            bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl{8,1} bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl{1,10} bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl(8,10) ]
[bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl([1 end],:) ; bet_lhs_mtc.cmb_reg.Mismatch_sensory.pvl_cnv(end,:)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_mtc_1450' 'pap_mtc_1450' 'pap_mtc_1450'};
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8] [12 13]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure4' '/' 'Fisher' '/'];
bet_rhs_mtc = mmil_fisher_exact_region_test(fcfg);

% MISMATCH OVERLAP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lng_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_lng_950/subjects/total/pap_lng_950_plt');
    con_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_con_950/subjects/total/pap_con_950_plt');
    phn_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_phn_950/subjects/total/pap_phn_950_plt');
    mtc_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_mtc_1450/subjects/total/pap_mtc_1450_plt');

cfg = [];
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
cfg.sbj_nme = 'total';
cfg.loc_typ = 'split';
cfg.ele_typ = 'ecog';
cfg.ecg_hms =  {'lhs' 'rhs'};
tot_ele = mmil_location_calc(cfg);

% Mtc & Vis
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [lng_ovr_lap(2:end,2) lng_ovr_lap(2:end,3) mtc_ovr_lap(2:end,3)];
cfg.eff_nme = {'Vis'                'Mtc'                'Bim'};
vis_and_mtc_tot = mmil_ovr_lap_ele_v3(cfg);

vis_and_mtc_tot.lhs.Mtc_in_Vis
vis_and_mtc_tot.lhs.Mtc_in_Vis_total

vis_and_mtc_tot.lhs.Vis_in_Mtc
vis_and_mtc_tot.lhs.Vis_in_Mtc_total

% Mtc & Fls
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [con_ovr_lap(2:end,2) con_ovr_lap(2:end,3) mtc_ovr_lap(2:end,3)];
cfg.eff_nme = {'Fls'                'Mtc'                'Bim'};
fls_and_mtc_tot = mmil_ovr_lap_ele_v3(cfg);

fls_and_mtc_tot.lhs.Mtc_in_Fls
fls_and_mtc_tot.lhs.Mtc_in_Fls_total

% Mtc & Aud
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [lng_ovr_lap(2:end,2) lng_ovr_lap(2:end,4) mtc_ovr_lap(2:end,3)];
cfg.eff_nme = {'Aud'                'Mtc'                'Bim'};
aud_and_mtc_tot = mmil_ovr_lap_ele_v3(cfg);

aud_and_mtc_tot.lhs.Mtc_in_Aud
aud_and_mtc_tot.lhs.Mtc_in_Aud_total

% Mtc & Nse
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [con_ovr_lap(2:end,2) con_ovr_lap(2:end,4) mtc_ovr_lap(2:end,3)];
cfg.eff_nme = {'Nse'                'Mtc'                'Bim'};
nse_and_mtc_tot = mmil_ovr_lap_ele_v3(cfg);

nse_and_mtc_tot.lhs.Mtc_in_Nse
nse_and_mtc_tot.lhs.Mtc_in_Nse_total

% Mtc & Ltr
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [phn_ovr_lap(2:end,2) phn_ovr_lap(2:end,3) mtc_ovr_lap(2:end,3)];
cfg.eff_nme = {'Ltr'                'Mtc'                'Bim'};
ltr_and_mtc_tot = mmil_ovr_lap_ele_v3(cfg);

ltr_and_mtc_tot.lhs.Mtc_in_Ltr
ltr_and_mtc_tot.lhs.Mtc_in_Ltr_total

% Mtc & Phn
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [phn_ovr_lap(2:end,2) phn_ovr_lap(2:end,4) mtc_ovr_lap(2:end,3)];
cfg.eff_nme = {'Phn'                'Mtc'                'Bim'};
phn_and_mtc_tot = mmil_ovr_lap_ele_v3(cfg);

phn_and_mtc_tot.lhs.Mtc_in_Phn
phn_and_mtc_tot.lhs.Mtc_in_Phn_total

%% CONTROL-SPECIFIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_con_950' 'pap_con_950' };
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'Fisher' '/'];
hms_sig_con = mmil_fisher_exact_region_test(fcfg);

% hms_sig_con.cmb_reg.VisualControl
hms_sig_con.cmb_reg.AuditoryControl

% WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_con_950' 'pap_con_950' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'Fisher' '/'];
wth_lhs_con = mmil_fisher_exact_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_con_950' 'pap_con_950' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.use_col = {[1 2] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'Fisher' '/'];
wth_rhs_con = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_con.cmb_reg.pvl_cnv wth_lhs_con.cmb_reg.pvl(:,2)]
[wth_rhs_con.cmb_reg.pvl_cnv wth_rhs_con.cmb_reg.pvl(:,2)]

% BETWEEN REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_con_950' 'pap_con_950' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'Fisher' '/'];
bet_lhs_con = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_lhs_con.cmb_reg.AuditoryControl.pvl{8,1} bet_lhs_con.cmb_reg.AuditoryControl.pvl{1,6}  bet_lhs_con.cmb_reg.AuditoryControl.pvl(8,6) ; ...
            bet_lhs_con.cmb_reg.AuditoryControl.pvl{8,1} bet_lhs_con.cmb_reg.AuditoryControl.pvl{1,7}  bet_lhs_con.cmb_reg.AuditoryControl.pvl(8,7) ; ...
            bet_lhs_con.cmb_reg.AuditoryControl.pvl{8,1} bet_lhs_con.cmb_reg.AuditoryControl.pvl{1,10} bet_lhs_con.cmb_reg.AuditoryControl.pvl(8,10) ]
[bet_lhs_con.cmb_reg.AuditoryControl.pvl([1 end],:) ; bet_lhs_con.cmb_reg.AuditoryControl.pvl_cnv(end,:)]

spc_tst = [ bet_lhs_con.cmb_reg.VisualControl.pvl{8,1} bet_lhs_con.cmb_reg.VisualControl.pvl{1,6}  bet_lhs_con.cmb_reg.VisualControl.pvl(8,6) ; ...
            bet_lhs_con.cmb_reg.VisualControl.pvl{8,1} bet_lhs_con.cmb_reg.VisualControl.pvl{1,7}  bet_lhs_con.cmb_reg.VisualControl.pvl(8,7) ; ...
            bet_lhs_con.cmb_reg.VisualControl.pvl{8,1} bet_lhs_con.cmb_reg.VisualControl.pvl{1,10} bet_lhs_con.cmb_reg.VisualControl.pvl(8,10) ]
[bet_lhs_con.cmb_reg.VisualControl.pvl([1 end],:) ; bet_lhs_con.cmb_reg.VisualControl.pvl_cnv(end,:)]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_con_950' 'pap_con_950' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure5' '/' 'Fisher' '/'];
bet_rhs_con = mmil_fisher_exact_region_test(fcfg);

% TEXT/CONTROL OVERLAP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
lng_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_lng_950/subjects/total/pap_lng_950_plt');
con_ovr_lap = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/sig_chn/hgp/ecog/split/pap_con_950/subjects/total/pap_con_950_plt');

vis_ele = sum(cell2mat(lng_ovr_lap(2:end,3))); 
nse_ele = sum(cell2mat(con_ovr_lap(2:end,4)));
bim_ele = sum(cell2mat(lng_ovr_lap(2:end,3)) & cell2mat(con_ovr_lap(2:end,4)));

cfg = [];
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
cfg.sbj_nme = 'total';
cfg.loc_typ = 'split';
cfg.ele_typ = 'ecog';
cfg.ecg_hms =  {'lhs' 'rhs'};
tot_ele = mmil_location_calc(cfg);

% Calculations
cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [lng_ovr_lap(2:end,2) lng_ovr_lap(2:end,3) con_ovr_lap(2:end,4)];
cfg.eff_nme = {'Vis'                'Nse'                'Bim'};
vis_and_nse = mmil_ovr_lap_ele_v3(cfg);

vis_and_nse.lhs.Nse_in_Vis
vis_and_nse.lhs.Nse_in_Vis_total

vis_and_nse.lhs.Vis_in_Nse
vis_and_nse.lhs.Vis_in_Nse_total

cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [lng_ovr_lap(2:end,2) lng_ovr_lap(2:end,4) con_ovr_lap(2:end,3)];
cfg.eff_nme = {'Fls'                'Aud'                'Bim'};
aud_and_fls = mmil_ovr_lap_ele_v3(cfg);

aud_and_fls.lhs.Aud_in_Fls_total
aud_and_fls.lhs.Fls_in_Aud_total

cfg = [];
cfg.nme    = nme;
cfg.reg    = reg;
cfg.reg_nme = 'cmb_reg';
cfg.tot_ele = tot_ele;
cfg.ecg_hms = {'lhs' 'rhs'};
cfg.eff_ele = [lng_ovr_lap(2:end,2) lng_ovr_lap(2:end,4) con_ovr_lap(2:end,4)];
cfg.eff_nme = {'Aud'                'Nse'                'Bim'};
aud_and_fls = mmil_ovr_lap_ele_v3(cfg);

aud_and_fls.lhs.Nse_in_Aud_total

%% WORD/PSEUDOWORD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950'};
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] [7 8] [12 13] [17 18]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'Fisher' '/'];
hms_sig_wrd = mmil_fisher_exact_region_test(fcfg);

% WITHIN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950'};
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8] [12 13] [17 18]};
fcfg.use_col = {[1 2] [3 4] [1 3] [2 4]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'Fisher' '/'];
wth_lhs_wrd = mmil_fisher_exact_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950'};
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8] [12 13] [17 18]};
fcfg.use_col = {[1 2] [3 4] [1 3] [2 4] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'Fisher' '/'];
wth_rhs_wrd = mmil_fisher_exact_region_test(fcfg);

[wth_lhs_phn.cmb_reg.pvl_cnv wth_lhs_phn.cmb_reg.pvl(:,2)]
[wth_rhs_phn.cmb_reg.pvl_cnv wth_rhs_phn.cmb_reg.pvl(:,2)]

% BETWEEN REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950'};
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [7 8] [12 13] [17 18]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'Fisher' '/'];
bet_lhs_wrd = mmil_fisher_exact_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/SL/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950' 'pap_wrd_950'};
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [7 8] [12 13] [17 18]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure9' '/' 'Fisher' '/'];
bet_rhs_wrd = mmil_fisher_exact_region_test(fcfg);



