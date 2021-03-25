clc; clear;

%% FisherExactTest for Language-Sensitive Electrodes
% SETUP REGIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reg.oct_reg = { 'Caudal ITG' 'Middle ITG' 'Rostral ITG' 'Caudal Fusiform' 'Middle Fusiform' 'Rostral Fusiform' 'Lateral Occipital' };
reg.par_reg = { 'Supramarginal' 'Inferior Parietal' 'Superior Parietal' };
reg.tmp_reg = { 'Caudal MTG' 'Middle MTG' 'Rostral MTG' 'Caudal STG' 'Middle STG' 'Rostral STG' };
reg.rol_reg = { 'Inferior Precentral' 'Middle Precentral' 'Inferior Postcentral' 'Middle Postcentral' };
reg.frn_reg = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' 'Pars Opercularis' 'Pars Triangularis' 'Pars Orbitalis'};

cmb_reg.occ = { 'Lateral Occipital' };
cmb_reg.fus = { 'Caudal Fusiform' 'Middle Fusiform' };
cmb_reg.itg = { 'Caudal ITG' 'Middle ITG' };
cmb_reg.mtg = { 'Caudal MTG' 'Middle MTG' };
cmb_reg.stg = { 'Caudal STG' 'Middle STG' };
cmb_reg.par = { 'Inferior Parietal' 'Superior Parietal' };
cmb_reg.sup = { 'Supramarginal' };
cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
cmb_reg.pre = { 'Inferior Precentral' 'Middle Precentral' };
cmb_reg.opc = { 'Pars Opercularis' };
cmb_reg.tri = { 'Pars Triangularis' };
% cmb_reg.orb = { 'Pars Orbitalis' };
cmb_reg.mid = { 'rostral-middlefrontal' 'middle-middlefrontal' 'Caudal Middle Frontal' };

% LOAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nme = fieldnames(reg);
cmb_nme = fieldnames(cmb_reg);

reg.tot_reg = cell(0);
for iN = 1:numel(nme)
    reg.tot_reg = [reg.tot_reg reg.(nme{iN})];
end
nme{end+1} = 'tot_reg';

cmb_reg.tot_reg = cell(0);
for iN = 1:numel(cmb_nme)
    cmb_reg.tot_reg = [cmb_reg.tot_reg cmb_reg.(cmb_nme{iN})];
end
cmb_nme{end+1} = 'tot_reg';

% COMPARING HEMISPHERES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tst_nme = 'pap_rsp_600';
hms     = {'lhs' 'rhs'};
col     = [2 3];

tst.lhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_lhs_table_plot']);
tst.lhs = tst.lhs(2:end-1,:);

tst.rhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_rhs_table_plot']);
tst.rhs = tst.rhs(2:end-1,:);

reg_hld = cell(1,1);
for iR = 1:numel(nme)
    
    reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.(nme{iR}),'uni',0)));
    reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], reg.(nme{iR}),'uni',0)));
    
    ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
    for iN = 1:numel(reg_1st)
        ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col(1)};
        tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col(2)};
        
        ele_num_2nd = ele_num_2nd + tst.(hms{2}){reg_2nd(iN),col(1)};
        tot_num_2nd = tot_num_2nd + tst.(hms{2}){reg_2nd(iN),col(2)};
    end
    
    yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
            ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
  
    zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
            ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
    
    [sig,pvl] = FisherExactTest(yyy,zzz);
    
    reg_hld{iR,1} = nme{iR};
    reg_hld{iR,2} = sig;
    reg_hld{iR,3} = pvl;
    reg_hld{iR,4} = crosstab(yyy,zzz);
        
end
pct_hld = cat(3,reg_hld{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
set(gca,'XTickLabel',nme)

cmb_hld = cell(1,1);
for iR = 1:numel(cmb_nme)
    
    reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
    reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
    
    ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
    for iN = 1:numel(reg_1st)
        ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st(iN),col(1)};
        tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st(iN),col(2)};
        
        ele_num_2nd = ele_num_2nd + tst.(hms{2}){reg_2nd(iN),col(1)};
        tot_num_2nd = tot_num_2nd + tst.(hms{2}){reg_2nd(iN),col(2)};
    end
    
    yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
            ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
  
    zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
            ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
    
    [sig,pvl] = FisherExactTest(yyy,zzz);
    
    cmb_hld{iR,1} = cmb_nme{iR};
    cmb_hld{iR,2} = sig;
    cmb_hld{iR,3} = pvl;
    cmb_hld{iR,4} = crosstab(yyy,zzz);
        
end
pct_hld = cat(3,cmb_hld{:,4});
pct_bar = [squeeze(pct_hld(1,1,:))./(squeeze(pct_hld(1,2,:))+squeeze(pct_hld(1,1,:))) squeeze(pct_hld(2,1,:))./(squeeze(pct_hld(2,2,:))+squeeze(pct_hld(2,1,:)))];
figure()
bar(pct_bar)
set(gca,'XTickLabel',cmb_nme)

spf_hld = cell(1,1);
for iR = 1:numel(reg.tot_reg)
            
    reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], reg.tot_reg(iR),'uni',0)));
    reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], reg.tot_reg(iR),'uni',0)));
    
    ele_num_1st = 0; tot_num_1st = 0; ele_num_2nd = 0; tot_num_2nd = 0;
    ele_num_1st = ele_num_1st + tst.(hms{1}){reg_1st,col(1)};
    tot_num_1st = tot_num_1st + tst.(hms{1}){reg_1st,col(2)};
    
    ele_num_2nd = ele_num_2nd + tst.(hms{2}){reg_2nd,col(1)};
    tot_num_2nd = tot_num_2nd + tst.(hms{2}){reg_2nd,col(2)};

    yyy = [ ones(1,ele_num_1st)   ones(1,tot_num_1st-ele_num_1st)   ...
            ones(1,ele_num_2nd)*2 ones(1,tot_num_2nd-ele_num_2nd)*2 ];
  
    zzz = [ ones(1,ele_num_1st) ones(1,tot_num_1st-ele_num_1st)*2 ...
            ones(1,ele_num_2nd) ones(1,tot_num_2nd-ele_num_2nd)*2 ];
    
    [sig,pvl] = FisherExactTest(yyy,zzz);
    
    spf_hld{iR,1} = reg.tot_reg{iR};
    spf_hld{iR,2} = sig;
    spf_hld{iR,3} = pvl;
    spf_hld{iR,4} = crosstab(yyy,zzz);
            
end

% COMPARING WITHIN HEMISPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tst_nme = 'pap_rsp_600';
hms     = {'lhs' 'lhs'};
col     = [2 3];

tst.lhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_lhs_table_plot']);
tst.lhs = tst.lhs(2:end-1,:);

tst.rhs = mmil_readtext(['/space/mdeh4/1/halgdev/projects/mmilanguage/FW/clerical/sig_chn/hgp/ecog/split/' tst_nme '/total/total_' tst_nme '_rhs_table_plot']);
tst.rhs = tst.rhs(2:end-1,:);


lef_sig = cell(1,1);
lef_pvl = cell(1,1);
lef_crs = cell(1,1);
ccc     = cell(0,0);
rrr     = cell(0,0);
crcrcrc = cell(0,0);
for iR = 1:numel(cmb_nme)-1
    
    for iC = 1:numel(cmb_nme)-1
        
        if iR ~= iC
            reg_1st = find(ismember(tst.(hms{1})(:,1),cellfun(@(x) [hms{1} '_' x], cmb_reg.(cmb_nme{iR}),'uni',0)));
            reg_2nd = find(ismember(tst.(hms{2})(:,1),cellfun(@(x) [hms{2} '_' x], cmb_reg.(cmb_nme{iC}),'uni',0)));
            
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
            
            lef_sig{iR+1,1} = cmb_nme{iR};
            lef_sig{1,iC+1} = cmb_nme{iC};
            
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





