clc; clear;

%% FisherExactTest for Language-Sensitive Electrodes
% SETUP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg.tot_reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.tot_reg.par_reg = {  'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tot_reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.tot_reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.tot_reg.frn_reg = {'middle-middlefrontal' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'rostral-middlefrontal' 'caudal-middlefrontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};
reg.tot_reg.ifg_reg = { 'Pars Opercularis' 'Pars Triangularis'};

reg.cmb_reg.occ = { 'Lateral Occipital' };
reg.cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
reg.cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
reg.cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
% reg.cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
reg.cmb_reg.sup = { 'Supramarginal' };
reg.cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
reg.cmb_reg.pos = { 'Inferior Postcentral' 'Middle Postentral' };
reg.cmb_reg.opc = { 'Pars Opercularis' };
reg.cmb_reg.tri = { 'Pars Triangularis' };
reg.cmb_reg.mid = { 'middle-middlefrontal' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'rostral-middlefrontal' 'caudal-middlefrontal' };

reg.rst_reg.lat_occ = { 'Lateral Occipital' };
reg.rst_reg.fus     = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.rst_reg.inf_par = { 'Inferior Parietal' };
reg.rst_reg.sup     = { 'Supramarginal' };
reg.rst_reg.mtg     = { 'Caudal MTG' 'Middle MTG' };
reg.rst_reg.stg     = { 'Caudal STG' 'Middle STG' };
reg.rst_reg.prc     = { 'Inferior Precentral' 'Middle Precentral' };
reg.rst_reg.par_opc = { 'Pars Opercularis' };
reg.rst_reg.par_tri = { 'Pars Triangularis' };

reg.ana_reg.lat_occ = { 'Lateral Occipital' };
reg.ana_reg.fus     = { 'Caudal Fusiform' 'Middle Fusiform' };
reg.ana_reg.itg     = {'Caudal ITG' 'Middle ITG' };
reg.ana_reg.inf_par = { 'Inferior Parietal' };
reg.ana_reg.sup_par = { 'Superior Parietal' };
reg.ana_reg.sup     = { 'Supramarginal' };
reg.ana_reg.mtg     = { 'Caudal MTG' 'Middle MTG' };
reg.ana_reg.stg     = { 'Caudal STG' 'Middle STG' };
reg.ana_reg.prc     = { 'Inferior Precentral' 'Middle Precentral' };
reg.ana_reg.pos     = { 'Inferior Postcentral' 'Middle Postcentral' };
reg.ana_reg.par_opc = { 'Pars Opercularis' };
reg.ana_reg.par_tri = { 'Pars Triangularis' };
reg.ana_reg.par_orb = { 'Pars Orbitalis' };
reg.ana_reg.mid_frn = { 'middle-middlefrontal' 'Caudal Middle Frontal' 'Middle Middle Frontal' 'rostral-middlefrontal' 'caudal-middlefrontal' };

reg_nme = fieldnames(reg);

for iN = 1:numel(reg_nme)
    nme.(reg_nme{iN}) = fieldnames(reg.(reg_nme{iN}));
    
    reg.(reg_nme{iN}).tot_reg = cell(0);
    for iNM = 1:numel(nme.(reg_nme{iN}))
        reg.(reg_nme{iN}).tot_reg = [reg.(reg_nme{iN}).tot_reg reg.(reg_nme{iN}).(nme.(reg_nme{iN}){iNM})];
    end
    nme.(reg_nme{iN}){end+1} = 'tot_reg';
end

%% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcfg = [];
fcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'hem';
fcfg.loc_nme = {'pap_rsp_600' };
fcfg.hms     = {'lhs' 'rhs'};
fcfg.col     = {[2 3] };
fcfg.plt     = 0;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure1' '/' 'Fisher' '/'];
hms_sig = mmil_fisher_exact_region_test(fcfg);

%%
fcfg = [];
fcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_rsp_600' };
fcfg.hms     = {'lhs' 'lhs'};
fcfg.col     = {[2 3] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure1' '/' 'Fisher' '/'];
bet_sig = mmil_fisher_exact_region_test(fcfg);

fcfg = [];
fcfg.clr_fld = '/home/ekaestne/PROJECTS/OUTPUT/FW/clerical/';
fcfg.nme     = nme;
fcfg.reg     = reg;
fcfg.typ     = 'bet_reg';
fcfg.loc_nme = {'pap_rsp_600' };
fcfg.hms     = {'rhs' 'rhs'};
fcfg.col     = {[2 3] };
fcfg.plt     = 1;
fcfg.out_plt = [fcfg.clr_fld '/' 'manuscript' '/' 'figure1' '/' 'Fisher' '/'];
bet_sig_rhs = mmil_fisher_exact_region_test(fcfg);

%% COMPARING WITHIN HEMISPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tst_nme = 'pap_rsp_600';
hms     = {'lhs' 'lhs'};
col     = [2 3];

tst.lhs = mmil_readtext(['/home/ekaestne/PROJECTS/OUTPUT/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_lhs_table_plot']);
tst.lhs = tst.lhs(2:end-1,:);

tst.rhs = mmil_readtext(['/home/ekaestne/PROJECTS/OUTPUT/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_rhs_table_plot']);
tst.rhs = tst.rhs(2:end-1,:);

lef_sig = cell(1,1);
lef_pvl = cell(1,1);
lef_crs = cell(1,1);
ccc     = cell(0,0);
rrr     = cell(0,0);
crcrcrc = cell(0,0);
for iR = 1:numel(reg.cmb_reg.tot_reg)-1
    
    for iC = 1:numel(reg.cmb_reg.tot_reg)-1
        
        if iR ~= iC
            reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.cmb_reg.(reg.cmb_reg.tot_reg{iR}),'uni',0)));
            reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], reg.cmb_reg.(reg.cmb_reg.tot_reg{iC}),'uni',0)));
            
            ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
            for iN = 1:numel(reg_1st)
                ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col(1)};
                tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col(2)};
            end
            for iN = 1:numel(reg_2nd)
                ele_num_2nd = ele_num_2nd + tst.(hms{2}){reg_2nd(iN),col(1)};
                tot_num_2nd = tot_num_2nd + tst.(hms{2}){reg_2nd(iN),col(2)};
            end
                    
            yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
                    ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
            
            zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
                    ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
            
            [sig,pvl] = FisherExactTest(yyy,zzz);
            
            lef_sig{iR+1,1} = reg.cmb_reg.tot_reg{iR};
            lef_sig{1,iC+1} = reg.cmb_reg.tot_reg{iC};
            
            lef_pvl{iR+1,1} = cmb_nme{iR};
            lef_pvl{1,iC+1} = cmb_nme{iC};
            
            lef_crs{iR+1,1} = cmb_nme{iR};
            lef_crs{1,iC+1} = cmb_nme{iC};
            
            lef_sig{iR+1,iC+1} = sig;
            lef_pvl{iR+1,iC+1} = pvl;
            lef_crs{iR+1,iC+1} = crosstab(yyy,zzz);
            
            if ~sig && pvl>0.10 && ~any(strcmpi([cmb_nme{iR} cmb_nme{iC}],crcrcrc)) &&  ~any(strcmpi([cmb_nme{iC} cmb_nme{iR}],crcrcrc))
                rrr{end+1} = cmb_nme{iR};
                ccc{end+1} = cmb_nme{iC};
                crcrcrc{end+1} = [cmb_nme{iR} cmb_nme{iC}];
            end
            
        else
            lef_crs{iR+1,iC+1} = '-';
        end
        
    end
end

figure(); aaa = graph(rrr,ccc); plot(aaa);





