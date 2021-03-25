clear; clc

%% SETUP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
reg.cmb_reg.pos = { 'Inferior Postcentral' 'Middle Postcentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
% reg.cmb_reg.orb = { 'Pars Orbitalis' };
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

%% OVERALL (ANOVA) COMPARISONS
% COMPARING BETWEEN HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_rsp_600' };
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = { [2 3] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
hms_sig = mmil_fisher_exact_region_test(fcfg);

% COMPARING BETWEEN REGIONS
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_rsp_600' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = { [2 3] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
bet_sig = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_sig.cmb_reg.Specific.pvl{9,1} bet_sig.cmb_reg.Specific.pvl{1,6} bet_sig.cmb_reg.Specific.pvl(9,6) ; ...
            bet_sig.cmb_reg.Specific.pvl{9,1} bet_sig.cmb_reg.Specific.pvl{1,7} bet_sig.cmb_reg.Specific.pvl(9,7) ; ...
            bet_sig.cmb_reg.Specific.pvl{9,1} bet_sig.cmb_reg.Specific.pvl{1,8} bet_sig.cmb_reg.Specific.pvl(9,8) ]

bet_sig.cmb_reg.Specific.pvl([1 end],:)
        
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_rsp_600' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = { [2 3] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
bet_sig_rhs = mmil_fisher_exact_region_test(fcfg);

%% GRAINSIZE COMPARISONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING BETWEEN HEMISPHERES 
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_wrd_600' 'pap_wrd_600' 'pap_con_600'};
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3]         [12 13]       [2 3]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
hms_sig = mmil_fisher_exact_region_test(fcfg);

% COMPARING BETWEEN REGIONS
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_wrd_600' 'pap_wrd_600' 'pap_con_600'};
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3]         [12 13]       [2 3]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
bet_sig = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_sig.cmb_reg.VisualOrthographic.pvl{9,1} bet_sig.cmb_reg.VisualOrthographic.pvl{1,6} bet_sig.cmb_reg.VisualOrthographic.pvl(9,6) ; ...
            bet_sig.cmb_reg.VisualOrthographic.pvl{9,1} bet_sig.cmb_reg.VisualOrthographic.pvl{1,7} bet_sig.cmb_reg.VisualOrthographic.pvl(9,7) ; ...
            bet_sig.cmb_reg.VisualOrthographic.pvl{9,1} bet_sig.cmb_reg.VisualOrthographic.pvl{1,8} bet_sig.cmb_reg.VisualOrthographic.pvl(9,8) ]
bet_sig.cmb_reg.VisualOrthographic.pvl([1 end],:)

spc_tst = [ bet_sig.cmb_reg.VisualTotal.pvl{9,1} bet_sig.cmb_reg.VisualTotal.pvl{1,6} bet_sig.cmb_reg.VisualTotal.pvl(9,6) ; ...
            bet_sig.cmb_reg.VisualTotal.pvl{9,1} bet_sig.cmb_reg.VisualTotal.pvl{1,7} bet_sig.cmb_reg.VisualTotal.pvl(9,7) ; ...
            bet_sig.cmb_reg.VisualTotal.pvl{9,1} bet_sig.cmb_reg.VisualTotal.pvl{1,8} bet_sig.cmb_reg.VisualTotal.pvl(9,8) ]
bet_sig.cmb_reg.VisualTotal.pvl([1 end],:)

spc_tst = [ bet_sig.cmb_reg.VisualSensory.pvl{9,1} bet_sig.cmb_reg.VisualSensory.pvl{1,6} bet_sig.cmb_reg.VisualSensory.pvl(9,6) ; ...
            bet_sig.cmb_reg.VisualSensory.pvl{9,1} bet_sig.cmb_reg.VisualSensory.pvl{1,7} bet_sig.cmb_reg.VisualSensory.pvl(9,7) ; ...
            bet_sig.cmb_reg.VisualSensory.pvl{9,1} bet_sig.cmb_reg.VisualSensory.pvl{1,8} bet_sig.cmb_reg.VisualSensory.pvl(9,8) ]
bet_sig.cmb_reg.VisualSensory.pvl([1 end],:)

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_wrd_600' 'pap_wrd_600' 'pap_con_600'};
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3]         [12 13]       [2 3]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
bet_sig_rhs = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_sig_rhs.cmb_reg.VisualSensory.pvl{9,1} bet_sig_rhs.cmb_reg.VisualSensory.pvl{1,6} bet_sig_rhs.cmb_reg.VisualSensory.pvl(9,6) ; ...
            bet_sig_rhs.cmb_reg.VisualSensory.pvl{9,1} bet_sig_rhs.cmb_reg.VisualSensory.pvl{1,7} bet_sig_rhs.cmb_reg.VisualSensory.pvl(9,7) ; ...
            bet_sig_rhs.cmb_reg.VisualSensory.pvl{9,1} bet_sig_rhs.cmb_reg.VisualSensory.pvl{1,8} bet_sig_rhs.cmb_reg.VisualSensory.pvl(9,8) ]
bet_sig_rhs.cmb_reg.VisualSensory.pvl([1 end],:)

% COMPARING WITHIN REGIONS
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_wrd_600' 'pap_wrd_600' 'pap_con_600'};
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3]         [12 13]       [2 3]};
fcfg.use_col = {[3 2] [1 2]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
wth_sig = mmil_fisher_exact_region_test(fcfg);

wth_sig.cmb_reg.pvl_cnv = wth_sig.cmb_reg.pvl;
for iR = 2:size(wth_sig.cmb_reg.pvl,1)
    for iC = 2:3
        if wth_sig.cmb_reg.pvl{iR,iC} < 0.05
            wth_sig.cmb_reg.pvl_cnv{iR,iC} = 9;
        elseif wth_sig.cmb_reg.pvl{iR,iC} < 0.10
            wth_sig.cmb_reg.pvl_cnv{iR,iC} = 3;
        elseif wth_sig.cmb_reg.pvl{iR,iC} < 0.20
             wth_sig.cmb_reg.pvl_cnv{iR,iC} = 1;
        else
             wth_sig.cmb_reg.pvl_cnv{iR,iC} = 0;
        end
    end
end

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'wth_reg';
fcfg.loc_nme = {'pap_wrd_600' 'pap_wrd_600' 'pap_con_600'};
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3]         [12 13]       [2 3]};
fcfg.use_col = {[3 2] [1 2]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure2' '/' 'Fisher' '/'];
wth_sig_rhs = mmil_fisher_exact_region_test(fcfg);

wth_sig_rhs.cmb_reg.pvl_cnv = wth_sig_rhs.cmb_reg.pvl;
for iR = 2:size(wth_sig_rhs.cmb_reg.pvl,1)
    for iC = 2:3
        if wth_sig_rhs.cmb_reg.pvl{iR,iC} < 0.05
            wth_sig_rhs.cmb_reg.pvl_cnv{iR,iC} = 9;
        elseif wth_sig_rhs.cmb_reg.pvl{iR,iC} < 0.10
            wth_sig_rhs.cmb_reg.pvl_cnv{iR,iC} = 3;
        elseif wth_sig_rhs.cmb_reg.pvl{iR,iC} < 0.20
             wth_sig_rhs.cmb_reg.pvl_cnv{iR,iC} = 1;
        else
             wth_sig_rhs.cmb_reg.pvl_cnv{iR,iC} = 0;
        end
    end
end

%% Lex Locations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_lex_600' 'pap_lex_600'};
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] [12 13]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
hms_sig = mmil_fisher_exact_region_test(fcfg);

% COMPARING BETWEEN REGIONS
fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_lex_600' 'pap_lex_600'};
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] [12 13]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
bet_sig = mmil_fisher_exact_region_test(fcfg);

spc_tst = [ bet_sig.cmb_reg.VisualFrq.pvl{9,1} bet_sig.cmb_reg.VisualFrq.pvl{1,6} bet_sig.cmb_reg.VisualFrq.pvl(9,6) ; ...
            bet_sig.cmb_reg.VisualFrq.pvl{9,1} bet_sig.cmb_reg.VisualFrq.pvl{1,7} bet_sig.cmb_reg.VisualFrq.pvl(9,7) ; ...
            bet_sig.cmb_reg.VisualFrq.pvl{9,1} bet_sig.cmb_reg.VisualFrq.pvl{1,8} bet_sig.cmb_reg.VisualFrq.pvl(9,8) ]

spc_tst = [ bet_sig.cmb_reg.VisualNovelPriming.pvl{9,1} bet_sig.cmb_reg.VisualNovelPriming.pvl{1,6} bet_sig.cmb_reg.VisualNovelPriming.pvl(9,6) ; ...
            bet_sig.cmb_reg.VisualNovelPriming.pvl{9,1} bet_sig.cmb_reg.VisualNovelPriming.pvl{1,7} bet_sig.cmb_reg.VisualNovelPriming.pvl(9,7) ; ...
            bet_sig.cmb_reg.VisualNovelPriming.pvl{9,1} bet_sig.cmb_reg.VisualNovelPriming.pvl{1,8} bet_sig.cmb_reg.VisualNovelPriming.pvl(9,8) ]

fcfg = [];
fcfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_lex_600' 'pap_lex_600'};
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] [12 13]};
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure3' '/' 'Fisher' '/'];
bet_sig_rhs = mmil_fisher_exact_region_test(fcfg);

%% Lex Overlap
frq_ovr = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/pap_lex_600/subjects/total' '/' 'pap_lex_600_plt']);
loc_ovr = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/electrode_location_files/total/output' '/' 'total_lhs_ecog']);
[~,frq_eff_lhs,~] = intersect(frq_ovr(:,2),loc_ovr(:,1));

lex_eff = cell2mat(frq_ovr(frq_eff_lhs,5)); % lex_eff = cell2mat(frq_ovr(2:end,5));
rep_eff = cell2mat(frq_ovr(frq_eff_lhs,3)); % rep_eff = cell2mat(frq_ovr(2:end,3));

sum(lex_eff & rep_eff) / sum(lex_eff);

%
pap_rep_act = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/pap_lex_600/subjects/total/pap_lex_600_plt');
pap_rep_act = pap_rep_act(2:end,:);

pap_frq_act = mmil_readtext('/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/pap_lex_600/subjects/total/pap_lex_600_plt');
pap_frq_act = pap_frq_act(2:end,:);

cfg = [];
cfg.clr_fld = '/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/';
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
cfg.eff_ele = [pap_rep_act(:,2)     pap_rep_act(:,3)     pap_frq_act(:,5)];
cfg.eff_nme = {'rep'                'frq'                'Bim'};
vis_and_aud_tot = mmil_ovr_lap_ele_v3(cfg);

% vis_and_aud_tot.lhs.frq_in_rep
vis_and_aud_tot.lhs.frq_in_rep_total

% vis_and_aud_tot.lhs.frq_in_rep
% vis_and_aud_tot.lhs.frq_in_rep_total

%% Lex Overlap with Grainsize
frq_ovr = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/pap_lex_600/subjects/total' '/' 'pap_lex_600_plt']);
wrd_ovr = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/pap_wrd_600/subjects/total' '/' 'pap_wrd_600_plt']);
con_ovr = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/pap_con_600/subjects/total' '/' 'pap_con_600_plt']);

% lhs
loc_lft = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/electrode_location_files/total/output' '/' 'total_lhs_ecog']);

[~,frq_eff_lhs,~] = intersect(frq_ovr(:,2),loc_lft(:,1));
[~,wrd_eff_lhs,~] = intersect(wrd_ovr(:,2),loc_lft(:,1));
[~,con_eff_lhs,~] = intersect(con_ovr(:,2),loc_lft(:,1));

nme_eff_lft = [frq_ovr(1,3) frq_ovr(1,5) wrd_ovr(1,3) wrd_ovr(1,5) con_ovr(1,3)];
lbl_eff_lft = [ frq_ovr(frq_eff_lhs,2) wrd_ovr(wrd_eff_lhs,2) con_ovr(con_eff_lhs,2)];
hld_eff_lft = cell2mat([ frq_ovr(frq_eff_lhs,3) frq_ovr(frq_eff_lhs,5) wrd_ovr(wrd_eff_lhs,3) wrd_ovr(wrd_eff_lhs,5) con_ovr(con_eff_lhs,3) ]); % lex_eff = cell2mat(frq_ovr(2:end,5));

% rhs
loc_rgh = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/electrode_location_files/total/output' '/' 'total_rhs_ecog']);

[~,frq_eff_rhs,~] = intersect(frq_ovr(:,2),loc_rgh(:,1));
[~,wrd_eff_rhs,~] = intersect(wrd_ovr(:,2),loc_rgh(:,1));
[~,con_eff_rhs,~] = intersect(con_ovr(:,2),loc_rgh(:,1));

nme_eff_rgh = [frq_ovr(1,3) frq_ovr(1,5) wrd_ovr(1,3) wrd_ovr(1,5) con_ovr(1,3)];
lbl_eff_rgh = [ frq_ovr(frq_eff_rhs,2) wrd_ovr(wrd_eff_rhs,2) con_ovr(con_eff_rhs,2)];
hld_eff_rgh = cell2mat([ frq_ovr(frq_eff_rhs,3) frq_ovr(frq_eff_rhs,5) wrd_ovr(wrd_eff_rhs,3) wrd_ovr(wrd_eff_rhs,5) con_ovr(con_eff_rhs,3) ]); % lex_eff = cell2mat(frq_ovr(2:end,5));

% Look at overlap with repetition & lexical
tot_nme{1} = 'Rep / Frq'; tot_pct(1,1) = round((sum(hld_eff_lft(:,2) & hld_eff_lft(:,1)) / sum(hld_eff_lft(:,2)))*100); tot_pct(1,2) = sum(hld_eff_lft(:,2) & hld_eff_lft(:,1)); tot_pct(1,3) =  sum(hld_eff_lft(:,2)); tot_pct(1,4) = sum(hld_eff_lft(:,1)); 
% Look at overlap of repetition with letterssum(hld_eff_lft(:,2) & hld_eff_lft(:,1)
tot_nme{2} = 'Rep / Ltr'; tot_pct(2,1) = round((sum(hld_eff_lft(:,3) & hld_eff_lft(:,1)) / sum(hld_eff_lft(:,1)))*100); tot_pct(2,2) = sum(hld_eff_lft(:,3) & hld_eff_lft(:,1)); tot_pct(2,3) =  sum(hld_eff_lft(:,1)); tot_pct(2,4) = sum(hld_eff_lft(:,3));
% Look at overlap of repetition with word
tot_nme{3} = 'Rep / Wrd'; tot_pct(3,1) = round((sum(hld_eff_lft(:,4) & hld_eff_lft(:,1)) / sum(hld_eff_lft(:,1)))*100); tot_pct(3,2) = sum(hld_eff_lft(:,4) & hld_eff_lft(:,1)); tot_pct(3,3) =  sum(hld_eff_lft(:,1)); tot_pct(3,4) = sum(hld_eff_lft(:,4));
% Look at overlap of repetition with false-font
tot_nme{4} = 'Rep / Flf'; tot_pct(4,1) = round((sum(hld_eff_lft(:,5) & hld_eff_lft(:,1)) / sum(hld_eff_lft(:,1)))*100); tot_pct(4,2) = sum(hld_eff_lft(:,5) & hld_eff_lft(:,1)); tot_pct(4,3) =  sum(hld_eff_lft(:,1)); tot_pct(4,4) = sum(hld_eff_lft(:,5));
% Look at overlap of lexical with letters
tot_nme{5} = 'Frq / Ltr'; tot_pct(5,1) = round((sum(hld_eff_lft(:,3) & hld_eff_lft(:,2)) / sum(hld_eff_lft(:,2)))*100); tot_pct(5,2) = sum(hld_eff_lft(:,3) & hld_eff_lft(:,2)); tot_pct(5,3) =  sum(hld_eff_lft(:,2)); tot_pct(5,4) = sum(hld_eff_lft(:,3));
% Look at overlap of lexical with word
tot_nme{6} = 'Frq / Wrd'; tot_pct(6,1) = round((sum(hld_eff_lft(:,4) & hld_eff_lft(:,2)) / sum(hld_eff_lft(:,2)))*100); tot_pct(6,2) = sum(hld_eff_lft(:,4) & hld_eff_lft(:,2)); tot_pct(6,3) =  sum(hld_eff_lft(:,2)); tot_pct(6,4) = sum(hld_eff_lft(:,4)); 
% Look at overlap of lexical with false-font
tot_nme{7} = 'Frq / Flf'; tot_pct(7,1) = round((sum(hld_eff_lft(:,5) & hld_eff_lft(:,2)) / sum(hld_eff_lft(:,2)))*100); tot_pct(7,2) = sum(hld_eff_lft(:,5) & hld_eff_lft(:,2)); tot_pct(7,3) =  sum(hld_eff_lft(:,2)); tot_pct(7,4) = sum(hld_eff_lft(:,5));
lft_out = [tot_nme'  num2cell(tot_pct)];

% Look at overlap with repetition & lexical
tot_nme{1} = 'Rep / Frq'; tot_pct(1,1) = round((sum(hld_eff_rgh(:,2) & hld_eff_rgh(:,1)) / sum(hld_eff_rgh(:,2)))*100); tot_pct(1,2) = sum(hld_eff_rgh(:,2) & hld_eff_rgh(:,1)); tot_pct(1,3) =  sum(hld_eff_rgh(:,2)); tot_pct(1,4) = sum(hld_eff_rgh(:,1)); 
% Look at overlap of repetition with letters
tot_nme{2} = 'Rep / Ltr'; tot_pct(2,1) = round((sum(hld_eff_rgh(:,3) & hld_eff_rgh(:,1)) / sum(hld_eff_rgh(:,1)))*100); tot_pct(2,2) = sum(hld_eff_rgh(:,3) & hld_eff_rgh(:,1)); tot_pct(2,3) =  sum(hld_eff_rgh(:,1)); tot_pct(2,4) = sum(hld_eff_rgh(:,3)); 
% Look at overlap of repetition with word
tot_nme{3} = 'Rep / Wrd'; tot_pct(3,1) = round((sum(hld_eff_rgh(:,4) & hld_eff_rgh(:,1)) / sum(hld_eff_rgh(:,1)))*100); tot_pct(3,2) = sum(hld_eff_rgh(:,4) & hld_eff_rgh(:,1)); tot_pct(3,3) =  sum(hld_eff_rgh(:,1)); tot_pct(3,4) = sum(hld_eff_rgh(:,4));
% Look at overlap of repetition with false-font
tot_nme{4} = 'Rep / Flf'; tot_pct(4,1) = round((sum(hld_eff_rgh(:,5) & hld_eff_rgh(:,1)) / sum(hld_eff_rgh(:,1)))*100); tot_pct(4,2) = sum(hld_eff_rgh(:,5) & hld_eff_rgh(:,1)); tot_pct(4,3) =  sum(hld_eff_rgh(:,1)); tot_pct(4,4) = sum(hld_eff_rgh(:,5));
% Look at overlap of lexical with letters
tot_nme{5} = 'Frq / Ltr'; tot_pct(5,1) = round((sum(hld_eff_rgh(:,3) & hld_eff_rgh(:,2)) / sum(hld_eff_rgh(:,2)))*100); tot_pct(5,2) = sum(hld_eff_rgh(:,3) & hld_eff_rgh(:,2)); tot_pct(5,3) =  sum(hld_eff_rgh(:,2)); tot_pct(5,4) = sum(hld_eff_rgh(:,3));
% Look at overlap of lexical with word 
tot_nme{6} = 'Frq / Wrd'; tot_pct(6,1) = round((sum(hld_eff_rgh(:,4) & hld_eff_rgh(:,2)) / sum(hld_eff_rgh(:,2)))*100); tot_pct(6,2) = sum(hld_eff_rgh(:,4) & hld_eff_rgh(:,2)); tot_pct(6,3) =  sum(hld_eff_rgh(:,2)); tot_pct(6,4) = sum(hld_eff_rgh(:,4));
% Look at overlap of lexical with false-font
tot_nme{7} = 'Frq / Flf'; tot_pct(7,1) = round((sum(hld_eff_rgh(:,5) & hld_eff_rgh(:,2)) / sum(hld_eff_rgh(:,2)))*100); tot_pct(7,2) = sum(hld_eff_rgh(:,5) & hld_eff_rgh(:,2)); tot_pct(7,3) =  sum(hld_eff_rgh(:,2)); tot_pct(7,4) = sum(hld_eff_rgh(:,5));
rgh_out = [tot_nme'  num2cell(tot_pct)];

% Compare Overall Overlap
myBinomTest(lft_out{1,3},lft_out{1,4},0.5,'greater')*2
myBinomTest(lft_out{1,3},lft_out{1,4},0.5,'lesser')*2

%
myBinomTest(lft_out{2,3},lft_out{2,4},0.5,'greater')*2
myBinomTest(lft_out{2,3},lft_out{2,4},0.5,'lesser')*2

myBinomTest(lft_out{3,3},lft_out{3,4},0.5,'greater')*2
myBinomTest(lft_out{3,3},lft_out{3,4},0.5,'lesser')*2

myBinomTest(lft_out{4,3},lft_out{4,4},0.5,'greater')*2
myBinomTest(lft_out{4,3},lft_out{4,4},0.5,'lesser')*2

%
myBinomTest(lft_out{5,3},lft_out{5,4},0.5,'greater')*2
myBinomTest(lft_out{5,3},lft_out{5,4},0.5,'lesser')*2

myBinomTest(lft_out{6,3},lft_out{6,4},0.5,'greater')*2
myBinomTest(lft_out{6,3},lft_out{6,4},0.5,'lesser')*2

%
myBinomTest(rgh_out{1,3},rgh_out{1,4},0.5,'greater')*2
myBinomTest(rgh_out{1,3},rgh_out{1,4},0.5,'lesser')*2

%
myBinomTest(rgh_out{2,3},rgh_out{2,4},0.5,'greater')*2
myBinomTest(rgh_out{2,3},rgh_out{2,4},0.5,'lesser')*2

myBinomTest(rgh_out{3,3},rgh_out{3,4},0.5,'greater')*2
myBinomTest(rgh_out{3,3},rgh_out{3,4},0.5,'lesser')*2

myBinomTest(rgh_out{4,3},rgh_out{4,4},0.5,'greater')*2
myBinomTest(rgh_out{4,3},rgh_out{4,4},0.5,'lesser')*2

% Compare Letter/Word/False-Font overlap 
ele_num_1st = lft_out{2,3}; tot_num_1st = lft_out{2,5}; ele_num_2nd = lft_out{3,3}; tot_num_2nd = lft_out{3,5};
yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
        ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
        ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
[sig,pvl] = FisherExactTest(yyy,zzz);

fprintf('%s : %f\n',lft_out{2,1},pvl)
crosstab(yyy,zzz)

ele_num_1st = lft_out{4,3}; tot_num_1st = lft_out{4,5}; ele_num_2nd = lft_out{3,3}; tot_num_2nd = lft_out{3,5};
yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
        ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
        ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
[sig,pvl] = FisherExactTest(yyy,zzz);

fprintf('%s : %f\n',lft_out{4,1},pvl)
crosstab(yyy,zzz)

ele_num_1st = lft_out{5,3}; tot_num_1st = lft_out{5,5}; ele_num_2nd = lft_out{6,3}; tot_num_2nd = lft_out{6,5};
yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
        ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
        ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
[sig,pvl] = FisherExactTest(yyy,zzz);

fprintf('%s : %f\n',lft_out{5,1},pvl)
crosstab(yyy,zzz)

ele_num_1st = lft_out{7,3}; tot_num_1st = lft_out{7,5}; ele_num_2nd = lft_out{6,3}; tot_num_2nd = lft_out{6,5};
yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
        ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
        ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
[sig,pvl] = FisherExactTest(yyy,zzz);

fprintf('%s : %f\n',lft_out{7,1},pvl)
crosstab(yyy,zzz)

